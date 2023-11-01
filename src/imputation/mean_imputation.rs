use crate::base::*;
use ndarray::prelude::*;
use std::io;

impl GenotypesAndPhenotypes {
    pub fn mean_imputation(&mut self) -> io::Result<&mut Self> {
        self.check().unwrap();
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (loci_idx, _loci_chr, _loci_pos) = self.count_loci().unwrap();
        let l = loci_idx.len() - 1;
        for j in 0..l {
            // Use the indexes of each locus
            let idx_ini = loci_idx[j] as usize;
            let idx_fin = if j < (l - 1) {
                loci_idx[j + 1] as usize
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
        }
        Ok(self)
    }
}

// Impute using mean allele frequencies across pools
pub fn impute_mean(
    file_sync_phen: &FileSyncPhen,
    filter_stats: &FilterStats,
    min_depth_set_to_missing: &f64,
    frac_top_missing_pools: &f64,
    frac_top_missing_loci: &f64,
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {
    // All non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
    let keep_p_minus_1 = false;
    let start = std::time::SystemTime::now();
    let mut genotypes_and_phenotypes = file_sync_phen
        .into_genotypes_and_phenotypes(filter_stats, keep_p_minus_1, n_threads)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Parsed the sync file into allele frequncies: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(min_depth_set_to_missing)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Set missing loci below the minimum depth: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_pools(frac_top_missing_pools)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest pools: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(frac_top_missing_loci)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest loci: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes.mean_imputation().unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Mean value imputation: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    // Remove 100% of the loci with missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&1.00)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Missing data removed, i.e. loci which cannot be imputed because of extreme sparsity: {} pools x {} loci | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        duration.as_secs()
    );
    // Output
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
            max_missingness_rate: 0.0,
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
            &0.2,
            &0.5,
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
        assert_eq!(0, 1);
    }
}
