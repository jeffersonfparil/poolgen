use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use nalgebra::{self, DMatrix, DVector};
use std::sync::{Arc, Mutex};
use std::fs::{File, OpenOptions};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::base::*;
use crate::io::sync::AlleleCountsOrFrequencies;
use crate::io::phen::{Phenotypes, load_phen};

use statrs::distribution::{Normal, ContinuousCDF};

pub fn gwalpha_base(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
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
    let X = locus_frequencies.matrix.clone();
    let bins = locus_counts_and_phenotypes.phenotypes.column(0).clone();
    let q_prime = locus_counts_and_phenotypes.phenotypes.column(1).clone();
    let sig = locus_counts_and_phenotypes.phenotypes[(0,2)];
    let min = locus_counts_and_phenotypes.phenotypes[(1,2)];
    let max = locus_counts_and_phenotypes.phenotypes[(2,2)];
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, p) =  X.shape();
    let m = bins.len();
    if n != m {
        return None
    }
    // Iterate across alleles
    for j in 0..p {
        let freqA = X.column(j).clone();
        let pA = (freqA * bins.transpose())[(0,0)];
        let mut binsA = DMatrix::from_element(n, 1, 0.0);
        let mut binsB = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n {
            binsA[i] = (0.0 - freqA[i]) * bins[i] / (0.0 - pA);
            binsB[i] = (1.0 - freqA[i]) * bins[i] / (1.0 - pA);
        }
        let mut percA = binsA.clone();
        let mut percB = binsB.clone();
        for i in 0..n {
            percA[i] = percA.view((0,0), (i, 0)).sum();
            percB[i] = percB.view((0,0), (i, 0)).sum();
        }
    }

    let out = vec!["x", "y", "z"].join(",").replace("\n,", "\n");
    Some(out)
}
