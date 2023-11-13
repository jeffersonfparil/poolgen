use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use crate::base::*;

fn factorial_log10(x: f64) -> io::Result<f64> {
    if x > 34.0 {
        return Err(Error::new(
            ErrorKind::Other,
            "\u{1F494} Input is far too big \u{1F6AB}",
        ));
    }
    let mut out: f64 = 0.0;
    for i in 2..(x + 1.0) as usize {
        out = out + f64::log10(i as f64);
    }
    Ok(out)
}

fn hypergeom_ratio(counts: &Array2<f64>, log_prod_fac_marginal_sums: &f64) -> io::Result<f64> {
    // Log-Product of counts
    let mut prod_fac_sums = 0.0;
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
        Err(_) => return None,
    };
    // Restrict so that the sum is less than or equal to 34, i.e. at n>34 : n! > f64::MAX
    let n = locus_counts.matrix.nrows();
    let p = locus_counts.matrix.ncols();
    let mut counts: Array2<f64> = Array2::from_shape_vec(
        (n, p),
        locus_counts
            .matrix
            .clone()
            .into_iter()
            .map(|x| x as f64)
            .collect::<Vec<f64>>(),
    )
    .unwrap();
    let total = counts.sum();
    if total > 34.0 {
        let coef: f64 = 34.0 / total;
        for i in 0..n {
            for j in 0..p {
                counts[(i, j)] = (counts[(i, j)] * coef).floor();
            }
        }
    }
    // Log-Product of the marginal sums (where C.row_sum() correspond to the the column marginal sums and vice versa)
    let row_sums: Array1<f64> = counts.sum_axis(Axis(1));
    let col_sums: Array1<f64> = counts.sum_axis(Axis(0));
    let mut log_prod_fac_marginal_sums = 0.0;
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
                    max = vec![
                        (row_sums[i] - counts.slice(s![i, 0..j]).sum()) as usize,
                        (col_sums[j] - counts.slice(s![0..i, j]).sum()) as usize,
                    ]
                    .into_iter()
                    .min()
                    .unwrap() as f64;
                    if (i == (n - 1)) | (j == (p - 1)) {
                        counts[(i, j)] = max;
                    } else if (i < max_i) | (j < max_j) {
                        counts[(i, j)] = 0.0;
                    } else {
                        counts[(i, j)] = max;
                    }
                }
            }
            let mut j: usize;
            let mut i: usize;
            for inv_j in 0..p {
                for inv_i in 0..n {
                    j = p - (inv_j + 1);
                    i = n - (inv_i + 1);
                    max = vec![
                        (row_sums[i] - counts.slice(s![i, 0..p]).sum()) as usize,
                        (col_sums[j] - counts.slice(s![0..n, j]).sum()) as usize,
                    ]
                    .into_iter()
                    .min()
                    .unwrap() as f64;
                    if max > 0.0 {
                        counts[(i, j)] = max;
                    }
                }
            }
            // Make sure we kept the marginal sums constant
            assert!(row_sums == counts.sum_axis(Axis(1)));
            assert!(col_sums == counts.sum_axis(Axis(0)));
            // Calculate hypergeometric ratio
            p_extremes += hypergeom_ratio(&counts, &log_prod_fac_marginal_sums).unwrap();
        }
    }
    let out = vec![
        locus_counts.chromosome.clone(),
        locus_counts.position.clone().to_string(),
        locus_counts.alleles_vector.clone().join(""),
        p_observed.to_string(),
        (p_observed + p_extremes).to_string(),
    ]
    .join(",")
        + "\n";

    Some(out)
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
            max_base_error_rate: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2, 0.2, 0.2],
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
