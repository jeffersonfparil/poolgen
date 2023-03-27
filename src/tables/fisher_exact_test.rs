use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix};

use crate::base::*;

fn factorial_log10(x: f64) -> io::Result<f64> {
    if x > 34.0 {
        return Err(Error::new(ErrorKind::Other, "\u{1F494} Input is far too big \u{1F6AB}"))
    }
    let mut out: f64 = 0.0;
    for i in 1..x as usize {
        out = out + f64::log10(i as f64);
    }
    Ok(out)
}

fn hypergeom_ratio(counts: &DMatrix<f64>, log_prod_fac_marginal_sums: &f64) -> io::Result<f64> {
    // Log-Product of counts
    let mut prod_fac_sums = 1 as f64;
    for i in counts.iter() {
        prod_fac_sums = prod_fac_sums + factorial_log10(*i).unwrap();
    }
    prod_fac_sums = prod_fac_sums + factorial_log10(counts.sum()).unwrap();
    // Calculate the p-value
    let p = f64::powf(10.0, log_prod_fac_marginal_sums - prod_fac_sums);
    Ok(p)
}

pub fn fisher(locus_counts: &mut LocusCounts, filter_stats: &FilterStats) -> Option<String> {
    let locus_counts = match locus_counts.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None
    };
    let (n, p) = locus_frequencies.matrix.shape();
    let frequencies: DMatrix<f64> = locus_frequencies.matrix.clone();
    let mut counts: DMatrix<f64> = DMatrix::from_element(n, p, 0.0);
    // Find the minimum frequency to get the maximum natural number restricted by f64, i.e. n=34 if n! > f64::MAX
    let mut s = (1.00 / frequencies.max()).ceil() as usize;
    s = match s < 34  {
        true => s,
        false => 34,
    };
    // Populate the counts matrix
    for i in 0..n {
        for j in 0..p {
            counts[(i, j)] = (s as f64 * frequencies[(i, j)]).ceil() as f64;
        }
    }
    // println!("COUNTS: {:?}", C);
    // Log-Product of the marginal sums (where C.row_sum() correspond to the the column marginal sums and vice versa)
    let row_sums = counts.column_sum().clone_owned();
    let col_sums = counts.row_sum().clone_owned();
    let mut log_prod_fac_marginal_sums = 1 as f64;
    for r in row_sums.iter() {
        log_prod_fac_marginal_sums = log_prod_fac_marginal_sums + factorial_log10(*r).unwrap();
    }
    for c in col_sums.iter() {
        log_prod_fac_marginal_sums = log_prod_fac_marginal_sums + factorial_log10(*c).unwrap();
    }
    // Define the observed hypergeometric ratio, i.e. p of the observed data
    let p_observed = hypergeom_ratio(&counts, &log_prod_fac_marginal_sums).unwrap();
    // Iterate across all possible combinations of counts with the same marginal sums
    let mut p_extremes: f64 = 0.0;
    let mut max: f64;
    for max_i in 0..n {
        for max_j in 0..p {
            for i in 0..n {
                for j in 0..p {
                    max = vec![(row_sums[(i, 0)] - counts.index((i, 0..j)).sum().ceil()) as usize,
                               (col_sums[(0, j)] - counts.index((0..i, j)).sum().ceil()) as usize]
                               .into_iter().min().unwrap() as f64;
                    if (i==(n-1)) | (j==(p-1)) {
                        counts[(i,j)] = max;
                    } else if (i < max_i) | (j < max_j) {
                        counts[(i,j)] = 0.0;
                    } else {
                        counts[(i,j)] = max;
                    }
                }
            }
            let mut j: usize;
            let mut i: usize;
            for inv_j in 0..p {
                for inv_i in 0..n {
                    j = p - (inv_j+1);
                    i = n - (inv_i+1);
                    max = vec![(row_sums[(i, 0)] - counts.index((i, 0..p)).sum().ceil()) as usize,
                               (col_sums[(0, j)] - counts.index((0..n, j)).sum().ceil()) as usize]
                               .into_iter().min().unwrap() as f64;
                    if max > 0.0 {
                        counts[(i,j)] = max;
                    }
                }
            }
            // Make sure we kept the marginal sums constant
            assert!(row_sums == counts.column_sum().clone_owned());
            assert!(col_sums == counts.row_sum().clone_owned());
            // Calculate hypergeometric ratio
            p_extremes += hypergeom_ratio(&counts, &log_prod_fac_marginal_sums).unwrap();
        }
    }
    let out = vec![locus_frequencies.chromosome.clone(),
                           locus_frequencies.position.clone().to_string(),
                           locus_frequencies.alleles_vector.clone().join(""),
                           p_observed.to_string(),
                           (p_observed + p_extremes).to_string()]
                      .join(",") + "\n";

    Some(out)
}
