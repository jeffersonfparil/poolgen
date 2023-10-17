use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use crate::base::*;

impl GenotypesAndPhenotypes {
    pub fn mean_imputation(&mut self) -> io::Result<&mut Self> {
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let (n_, l) = self.coverages.dim();
        let l_ = self.start_index_of_each_locus.len();
        assert_eq!(n, n_);
        assert_eq!(l, l_);
        for j in 0..l {
            // Use the indexes of the locus to set missing values to all alleles in the locus
            let idx_ini = self.start_index_of_each_locus[j] as usize;
            let idx_fin = if j < (l-1) {
                self.start_index_of_each_locus[j+1] as usize
            } else {
                self.start_index_of_each_locus[l-1] as usize
            };

            // We need to correct for imputations resulting in a sum of allele frequencies greater or less than 1
            let q = self.intercept_and_allele_frequencies.slice(s![.., idx_ini..idx_fin]);

            let idx_missing = self.intercept_and_allele_frequencies.column(j+1)
                .iter()
                .enumerate()
                .filter(|&(i, x)| x.is_nan())
                .map(|(i, x)| i)
                .collect::<Vec<usize>>();
            // println!("idx_missing={:?}", idx_missing);
            if idx_missing.len() > 0 {
                let x = self.intercept_and_allele_frequencies.column(j+1).to_owned();
                let mu = mean_axis_ignore_nan(&x, 0).unwrap();
                println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
                println!("x={:?}", x);
                println!("idx_missing={:?}", idx_missing);
                println!("mu={}", mu);
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
        let min_depth_set_to_missing = 5.0;
        let _ = frequencies_and_phenotypes.set_missing_by_depth(&min_depth_set_to_missing).unwrap();
        println!("frequencies_and_phenotypes.intercept_and_allele_frequencies={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies);
        let _ = frequencies_and_phenotypes.mean_imputation().unwrap();
        // println!("frequencies_and_phenotypes.intercept_and_allele_frequencies={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies);
        assert_eq!(0, 1);
    }
}
