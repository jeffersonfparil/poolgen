use crate::base::*;
use ndarray::prelude::*;
use std::io;
use statrs::distribution::{MultivariateNormal, Continuous};

// Simulate genotypes with some linkage disequilibrium

fn simulate_genotypes(n: u64, p: u64, n_chr: u64, max_bp: u64, r2_50_perc_bp: u64) -> io::Result<Array2<f64>> {
    // Define approximately equally sized chromosomes and coordinates of the characterised loci or markers or SNPs
    let mut chromosome_sizes = vec![p/n_chr; n_chr];
    let p_ = chromosome_sizes.iter().fold(0, |sum, &x| sum + x);
    if p_ < p {
        chromosome_sizes[n_chr-1] += p - p_;
    }

    let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();
    return(Array2::from_shape_elem((n, p), f64::NAN))
}