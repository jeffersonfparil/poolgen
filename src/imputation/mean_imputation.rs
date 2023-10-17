use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use crate::base::*;

pub fn mean_imputation(genotypes_and_phenotypes: &mut GenotypesAndPhenotypes) -> io::Result<&mut GenotypesAndPhenotypes> {
    Ok(genotypes_and_phenotypes)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_fisher() {
        // Expected
        let expected_output1: f64 = 2.0791812460476247;
        let expected_output2: f64 = 0.24705882352941286;
        let expected_output3 =
            "Chromosome1,12345,TC,0.24705882352941286,0.6073529411764731\n".to_owned();
        // Inputs
        let x: f64 = 5.0;
        let counts_f64: Array2<f64> =
            Array2::from_shape_vec((3, 2), vec![0.0, 3.0, 1.0, 5.0, 2.0, 6.0]).unwrap();
        let counts_u64: Array2<u64> =
            Array2::from_shape_vec((3, 2), vec![0, 3, 1, 5, 2, 6]).unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let mut locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["T".to_owned(), "C".to_owned()],
            matrix: counts_u64,
        };
        // Outputs
        let log10factorial = factorial_log10(x).unwrap();
        let hypergeom_pval = hypergeom_ratio(&counts_f64, &19.959563872703743).unwrap();
        let fisher_line = fisher(&mut locus_counts, &filter_stats).unwrap();
        // Assertions
        assert_eq!(expected_output1, log10factorial);
        assert_eq!(expected_output2, hypergeom_pval);
        assert_eq!(expected_output3, fisher_line);
    }
}
