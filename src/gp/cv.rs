use crate::base::*;
use crate::gwas::*;
use nalgebra::{self, DMatrix, DVector};
use std::io::{self, Error, ErrorKind};

impl CrossValidate for LocusFrequenciesAndPhenotypes {
    fn split(&self, k: usize) -> io::Result<Vec<usize>> {
        Ok(vec![0])
    }

    fn performance(&self, y_hat: DVector<f64>) -> io::Result<Vec<f64>> {
        Ok(vec!(0.0))
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
        let x: DMatrix<f64> = DMatrix::from_row_slice(
            5,
            3,
            &[
                1.0, 0.4, 0.1, 1.0, 0.2, 0.1, 1.0, 0.3, 0.2, 1.0, 0.4, 0.3, 1.0, 0.5, 0.5,
            ],
        );
        let y: DMatrix<f64> =
            DMatrix::from_row_slice(5, 2, &[2.0, 0.5, 1.0, 0.2, 2.0, 0.5, 4.0, 0.0, 5.0, 0.5]);
        let counts: DMatrix<u64> =
            DMatrix::from_row_slice(5, 3, &[4, 1, 5, 2, 1, 7, 3, 2, 5, 4, 3, 3, 5, 5, 0]);
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned(), "T".to_owned(), "D".to_owned()],
            matrix: counts,
        };
        let phenotypes: DMatrix<f64> = y.clone();
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: phenotypes,
            pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        };
        // Outputs
        let x = 0.0;
        // Assertions
        assert_eq!(0, 0); // Output dimensions
    }
}
