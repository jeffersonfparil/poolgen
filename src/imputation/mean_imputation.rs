use crate::base::*;
use ndarray::prelude::*;
use std::io;

impl GenotypesAndPhenotypes {
    pub fn mean_imputation(&mut self) -> io::Result<&mut Self> {
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (n_, l) = self.coverages.dim();
        let l_ = self.start_index_of_each_locus.len();
        assert_eq!(n, n_);
        assert_eq!(l, l_);
        for j in 0..l {
            // Use the indexes of each locus
            let idx_ini = self.start_index_of_each_locus[j] as usize;
            let idx_fin = if j < (l - 1) {
                self.start_index_of_each_locus[j + 1] as usize
            } else {
                p
            };
            let freqs: Array2<f64> = self
                .intercept_and_allele_frequencies
                .slice(s![.., idx_ini..idx_fin])
                .to_owned();
            let mean_freqs = mean_axis_ignore_nan(&freqs, 0).unwrap();
            let sum_freqs = mean_freqs.sum();
            // We need to correct for imputations resulting in a sum of allele frequencies greater or less than 1
            let mean_freqs = if sum_freqs != 1.0 {
                mean_freqs.map(|x| x / sum_freqs)
            } else {
                mean_freqs
            };
            // println!("mean_freqs={:?}", mean_freqs);
            for k in 0..freqs.ncols() {
                for i in 0..n {
                    if self.intercept_and_allele_frequencies[(i, idx_ini + k)].is_nan() {
                        self.intercept_and_allele_frequencies[(i, idx_ini + k)] = mean_freqs[k];
                    } else {
                        continue;
                    };
                }
            }
            // println!("self.start_index_of_each_locus={:?}", self.start_index_of_each_locus);
            // println!("n={}; p={}; n_={}; l={}; l_={}", n, p, n_, l, l_);
            // println!("self.start_index_of_each_locus.len()={:?}", self.start_index_of_each_locus.len());
            // println!("self.start_index_of_each_locus[l-2]={:?}", self.start_index_of_each_locus[l-2]);
            // println!("self.start_index_of_each_locus[l-1]={:?}", self.start_index_of_each_locus[l-1]);
            // println!("idx_ini={:?}", idx_ini);
            // println!("idx_fin={:?}", idx_fin);
            // println!("freqs={:?}", freqs);
            // println!("self.intercept_and_allele_frequencies.slice(s![.., idx_ini..idx_fin])={:?}", self.intercept_and_allele_frequencies.slice(s![.., idx_ini..idx_fin]));
        }
        Ok(self)
    }
}

// Impute using mean allele frequencies across pools
pub fn impute_mean(
    file_sync_phen: &FileSyncPhen,
    filter_stats: &FilterStats,
    min_depth_set_to_missing: &f64,
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {
    // All non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
    let keep_p_minus_1 = false;
    let mut genotypes_and_phenotypes = file_sync_phen
        .into_genotypes_and_phenotypes(filter_stats, keep_p_minus_1, n_threads)
        .unwrap();
    genotypes_and_phenotypes
        .set_missing_by_depth(min_depth_set_to_missing)
        .unwrap();
    genotypes_and_phenotypes.mean_imputation().unwrap();
    let out = genotypes_and_phenotypes
        .write_csv(filter_stats, keep_p_minus_1, out, n_threads)
        .unwrap();
    Ok(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_impute_mean() {
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "load".to_owned(),
        };
        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        println!(
            "Before simulating missing data:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let min_depth_set_to_missing = 5.0;
        let _ = frequencies_and_phenotypes
            .set_missing_by_depth(&min_depth_set_to_missing)
            .unwrap();
        println!(
            "Before imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let _ = frequencies_and_phenotypes.mean_imputation().unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );

        let outname = impute_mean(
            &file_sync_phen,
            &filter_stats,
            &min_depth_set_to_missing,
            &n_threads,
            &"test-impute_mean.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_mean.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)],
            Array1::from_vec(vec![0.3333333333333333, 0.2, 0.14285714285714285])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 1)],
            Array1::from_vec(vec![0.3333333333333333, 0.2, 0.14285714285714285])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 2)],
            Array1::from_vec(vec![0.6666666666666666, 0.8, 0.8571428571428571])
                .mean()
                .unwrap()
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 2)],
            Array1::from_vec(vec![0.6666666666666666, 0.8, 0.8571428571428571])
                .mean()
                .unwrap()
        );
        // assert_eq!(0, 1);
    }
}
