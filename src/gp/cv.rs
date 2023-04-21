use crate::base::*;
use crate::gwas::*;
use nalgebra::{self, DMatrix, DVector};
use std::io::{self, Error, ErrorKind};

impl CrossValidate for FrequenciesAndPhenotypes {
    fn split(&self, k: usize) -> io::Result<Vec<usize>> {
        let n = self.pool_names.len();
        let p = self.frequencies.len();
        Ok(vec![0])
    }

    fn performance(&self, y_hat: DVector<f64>) -> io::Result<Vec<f64>> {
        Ok(vec![0.0])
    }

    fn cross_validate(&self, k: usize) -> io::Result<PredictionPerformance> {
        Ok(PredictionPerformance {
            n: 0,
            p: 0,
            k: 0,
            model: "".to_owned(),
            rmse: DVector::from_column_slice(&[f64::NAN]),
            mse: DVector::from_column_slice(&[f64::NAN]),
            mae: DVector::from_column_slice(&[f64::NAN]),
            mbe: DVector::from_column_slice(&[f64::NAN]),
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_ols() {
        // Expected
        let expected_output1: DMatrix<f64> = DMatrix::from_column_slice(3, 1, &[-0.73, 5.53, 6.42]);
        // Inputs
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
        let n_threads = 2;
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let freqs = file_sync_phen
            .load(&filter_stats, true, &n_threads)
            .unwrap();
        let freqs_and_pheno = FrequenciesAndPhenotypes {
            frequencies: freqs,
            phenotypes: file_sync_phen.phen_matrix,
            pool_names: file_sync_phen.pool_names,
        };
        // Outputs
        let x = 0.0;
        // Assertions
        assert_eq!(0, 0); // Output dimensions
    }
}
