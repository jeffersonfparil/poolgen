use crate::base::*;
use ndarray::prelude::*;
use statrs::distribution::{ChiSquared, ContinuousCDF};

pub fn chisq(locus_counts: &mut LocusCounts, filter_stats: &FilterStats) -> Option<String> {
    let locus_counts = match locus_counts.filter(filter_stats) {
        Ok(Some(x)) => x,
        _ => return None,
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Marginal sums and total count
    let n = locus_frequencies.matrix.nrows();
    let p = locus_frequencies.matrix.ncols();
    let t = (n * p) as f64;
    let total = locus_frequencies.matrix.sum() as f64;
    let row_sums = locus_frequencies.matrix.sum_axis(Axis(1));
    let col_sums = locus_frequencies.matrix.sum_axis(Axis(0));
    // Caculate the chi-square test statistic
    let mut chi2: f64 = 0.0;
    let mut observed: f64;
    let mut expected: f64;
    for i in 0..n {
        for j in 0..p {
            observed = locus_frequencies.matrix[(i, j)] as f64;
            expected = (row_sums[i] * col_sums[j]) as f64 / total;
            chi2 += f64::powf(observed - expected, 2.0) / expected;
        }
    }
    // Calculate the one-tailed significance, i.e. how far away are we from the expected mean of the chi-squared distribution?
    let d = ChiSquared::new(t - 1.0).unwrap();
    // i.e. using the upper tail cdf of the chi2 distribution
    let pval = 1.00 - d.cdf(chi2);
    // Output line
    let out = vec![
        locus_frequencies.chromosome.clone(),
        locus_frequencies.position.to_string(),
        locus_frequencies.alleles_vector.join(""),
        parse_f64_roundup_and_own(chi2, 6),
        pval.to_string(),
    ]
    .join(",")
        + "\n";
    Some(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_chisq() {
        // Expected
        let expected_line = "Chromosome1,12345,AT,4,0.7797774084757156\n".to_owned(); // where df=7, then the pvalue is calculated as the lower tail because if chi2 < df
                                                                                      // Inputs
        let mut locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A", "T"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            matrix: Array2::from_shape_vec((4, 2), vec![0, 20, 20, 0, 0, 20, 20, 0]).unwrap(),
        };
        let filter_stats = FilterStats {
            remove_ns: true,
            keep_lowercase_reference: false,
            max_base_error_rate: 0.01,
            min_coverage_depth: 1,
            min_coverage_breadth: 1.0,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2],
        };
        // Outputs
        let out = chisq(&mut locus_counts, &filter_stats).unwrap();
        // Assertions
        assert_eq!(expected_line, out);
    }
}
