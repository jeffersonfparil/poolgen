use std::io::{self, Error, ErrorKind};

use crate::base::*;

impl Parse<FilterStats> for FilterStats {
    fn lparse(&self) -> io::Result<Box<FilterStats>> {
        let remove_ns = self.remove_ns.clone();
        let remove_monoallelic = self.remove_monoallelic.clone();
        let keep_lowercase_reference = self.keep_lowercase_reference.clone();
        let max_base_error_rate = if self.max_base_error_rate <= 1.0 && self.max_base_error_rate >= 0.0 {
            self.max_base_error_rate.clone()
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid range. max_base_error_rate must be between 0.0 and 1.0",
            ));
        };
        let min_coverage_depth = self.min_coverage_depth;
        let min_coverage_breadth = if self.min_coverage_breadth <= 1.0 && self.max_base_error_rate >= 0.0 {
            self.min_coverage_breadth.clone()
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid range. min_coverage_breadth must be between 0.0 and 1.0",
            ));
        };
        let min_allele_frequency = if self.min_allele_frequency <= 1.0 && self.min_allele_frequency >= 0.0 {
            self.min_allele_frequency.clone()
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid range. min_allele_frequency must be between 0.0 and 1.0",
            ));
        };
        let max_missingness_rate = if self.max_missingness_rate <= 1.0 && self.max_missingness_rate >= 0.0 {
            self.max_missingness_rate.clone()
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid range. max_missingness_rate must be between 0.0 and 1.0",
            ));
        };
        let pool_sizes = self.pool_sizes.clone();
        return Ok(Box::new(FilterStats {
            remove_ns,
            remove_monoallelic,
            keep_lowercase_reference,
            max_base_error_rate,
            min_coverage_depth,
            min_coverage_breadth,
            min_allele_frequency,
            max_missingness_rate,
            pool_sizes,
        }));
    }
}

