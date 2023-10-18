use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use crate::base::*;

impl GenotypesAndPhenotypes {
    pub fn set_missing_by_depth(
        &mut self,
        min_depth_set_to_missing: &f64,
    ) -> io::Result<&mut Self> {
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (n_, l) = self.coverages.dim();
        let l_ = self.start_index_of_each_locus.len();
        // println!("n={}; p={}; n_={}; l={}; l_={}", n, p, n_, l, l_);
        assert_eq!(n, n_);
        assert_eq!(l, l_);
        for i in 0..n {
            for j in 0..l {
                if self.coverages[(i, j)] < *min_depth_set_to_missing {
                    // Use the indexes of the locus to set missing values to all alleles in the locus
                    let idx_ini = self.start_index_of_each_locus[j] as usize;
                    let idx_fin = if j < (l - 1) {
                        self.start_index_of_each_locus[j + 1] as usize
                    } else {
                        self.start_index_of_each_locus[l - 1] as usize
                    };
                    for k in idx_ini..idx_fin {
                        self.intercept_and_allele_frequencies[(i, k)] = f64::NAN;
                    }
                }
            }
        }
        Ok(self)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_set_missing() {
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
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        // NOTE: MAKE SURE THAT FilterStats is set to no filtering except by zero frequency alleles as in below:
        //       AND THAT WE KEEP ALL NON-ZERO ALLELES ACROSS ALL ENTRIES!
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        let filter_stats = FilterStats {
            remove_ns: false,
            min_quality: 1.0,
            min_coverage: 0,
            min_allele_frequency: 0.000001,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        let min_depth_set_to_missing = 5.0;
        let _ = frequencies_and_phenotypes
            .set_missing_by_depth(&min_depth_set_to_missing)
            .unwrap();
        // println!("frequencies_and_phenotypes={:?}", frequencies_and_phenotypes);
        println!(
            "frequencies_and_phenotypes.chromosome[0..5]={:?}",
            &frequencies_and_phenotypes.chromosome[0..5]
        );
        println!(
            "frequencies_and_phenotypes.position[0..5]={:?}",
            &frequencies_and_phenotypes.position[0..5]
        );
        println!(
            "frequencies_and_phenotypes.allele[0..5]={:?}",
            &frequencies_and_phenotypes.allele[0..5]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies={:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        println!(
            "frequencies_and_phenotypes.coverages={:?}",
            frequencies_and_phenotypes.coverages
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)].is_nan(),
            true
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 4)],
            1.0 / 21.0
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 13)],
            2.0 / 10.0
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 18)].is_nan(),
            true
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 18)].is_nan(),
            true
        );
        // assert_eq!(0, 1);
    }
}
