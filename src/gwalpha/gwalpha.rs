use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix};
use crate::base::*;

use statrs::distribution::{StudentsT, ContinuousCDF};

fn neg_loglik_beta() -> io::Result<f64> {
    Ok(0.0)
}

pub fn gwalpha(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
                                                            .locus_counts
                                                            .filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None
    };
    // Extract the genotype and phenotypes
    let x_matrix = locus_frequencies.matrix.clone();
    let bins = locus_counts_and_phenotypes.phenotypes.column(0).clone();
    let q = locus_counts_and_phenotypes.phenotypes.column(1).clone();
    let sig = locus_counts_and_phenotypes.phenotypes[(0,2)];
    let min = locus_counts_and_phenotypes.phenotypes[(1,2)];
    let max = locus_counts_and_phenotypes.phenotypes[(2,2)];
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, p) =  x_matrix.shape();
    let m = bins.len();
    if n != m {
        return None
    }
    // Iterate across alleles
    for j in 0..p {
        let freqs_a = x_matrix.column(j).clone(); // allele frequencies per pool
        let p_a = (freqs_a * bins.transpose())[(0,0)]; // mean allele frequency across pools
        // Quantiles per pool (for least squares estimation)
        let mut q_prime: DMatrix<f64> = DMatrix::from_element(n+1, 1, 0.0);
        q_prime[0] = 0.0;
        for i in 0..n {
            q_prime[i+1] = (q[i] - min) / (max-min);
        }
        // Bins (sums up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut bins_a = DMatrix::from_element(n, 1, 0.0);
        let mut bins_b = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n {
            bins_a[i] = (0.0 - freqs_a[i]) * bins[i] / (0.0 - p_a);
            bins_b[i] = (1.0 - freqs_a[i]) * bins[i] / (1.0 - p_a);
        }
        // Percentiles (cummulative bins summing up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut percs_a = bins_a.clone();
        let mut percs_b = bins_b.clone();
        for i in 0..n {
            percs_a[i] = bins_a.view((0,0), (i, 0)).sum();
            percs_b[i] = bins_b.view((0,0), (i, 0)).sum();
        }
        // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
        let mut percs_a0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        let mut percs_b0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n-1 {
            percs_a0[i+1] = percs_a[i];
            percs_b0[i+1] = percs_b[i];
        }
    }

    println!("sig={}", sig);

    let out = vec!["x", "y", "z"].join(",").replace("\n,", "\n");
    Some(out)
}
