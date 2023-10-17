use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use crate::base::*;

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation_imputation(&mut self, min_depth_set_to_missing: &f64) -> io::Result<&mut Self> {
        let X = self.intercept_and_allele_frequencies.to_owned();
        let Y = self.coverages.to_owned();
        Ok(self)
    }    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_impute_aldknn() {
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
            .into_genotypes_and_phenotypes(&filter_stats, true, &n_threads)
            .unwrap();
        let min_depth_set_to_missing = 5.0;
        let _ = frequencies_and_phenotypes.adaptive_ld_knn_imputation_imputation(&min_depth_set_to_missing).unwrap();
        println!("frequencies_and_phenotypes={:?}", frequencies_and_phenotypes);
        // assert_eq!(0, 1);
    }
}
