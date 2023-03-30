use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DVector};
use crate::base::*;

use statrs::distribution::{StudentsT, ContinuousCDF};

fn pearsons_correlation(x: &DVector<f64>, y: &DVector<f64>) -> io::Result<(f64, f64)> {
    let n = x.len();
    // println!("x={:?}; y={:?}", x, y);
    if n != y.len() {
        return Err(Error::new(ErrorKind::Other, "Input vectors are not the same size."));
    }
    let mu_x = x.mean();
    let mu_y = y.mean();
    let x_less_mu_x = x.map(|x| x-mu_x);
    let y_less_mu_y = y.map(|y| y-mu_y);
    let x_less_mu_x_squared = x_less_mu_x.map(|x| x.powf(2.0));
    let y_less_mu_y_squared = y_less_mu_y.map(|y| y.powf(2.0));
    // println!("x_less_mu_x={:?}", x_less_mu_x);
    // println!("y_less_mu_y={:?}", y_less_mu_y);
    let numerator = x_less_mu_x.component_mul(&y_less_mu_y).sum();
    let denominator = x_less_mu_x_squared.sum().sqrt() * y_less_mu_y_squared.sum().sqrt();
    let r_tmp = numerator / denominator;
    let r = match r_tmp.is_nan() {
        true => 0.0,
        false => r_tmp,
    };
    // println!("numeratorr={:?}; demonitatorr={:?}; r={:?}", numerator, denominator, r);
    let sigma_r = ((1.0 - r.powf(2.0)) / (n as f64 - 2.0)).sqrt();
    let t = r / sigma_r;
    let d = StudentsT::new(0.0, 1.0, n as f64 - 1.0).unwrap();
    let pval = 2.00 * ( 1.00 - d.cdf(t.abs()));
    Ok((r, pval))
}

pub fn correlation(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
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
    let mut x_matrix = locus_frequencies.matrix.clone();
    let y_matrix = locus_counts_and_phenotypes.phenotypes.clone();
    // Keep p-1 alleles for conciseness
    let (_, p) =  x_matrix.shape();
    if p >= 2 {
        x_matrix = x_matrix.clone().remove_columns(p-1, 1);
    }
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, p) =  x_matrix.shape();
    let (m, k) = y_matrix.shape();
    if n != m {
        return None
    }
    // Iterate across alleles
    let (mut corr, mut pval): (f64, f64);
    let first_2_col = vec![locus_frequencies.chromosome, locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    for i in 0..p {
        let x = DVector::from(x_matrix.column(i));
        for j in 0..k {
            line.append(&mut first_2_col.clone());
            line.push(locus_frequencies.alleles_vector[i].clone());
            line.push(x.mean().to_string());
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            let y  = DVector::from(y_matrix.column(j));
            (corr, pval) = pearsons_correlation(&x, &y).unwrap();
            line.push(corr.to_string());
            line.push(pval.to_string() + "\n");
        }
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}
