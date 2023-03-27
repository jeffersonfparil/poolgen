use crate::base::*;
use statrs::distribution::{ChiSquared, ContinuousCDF};

pub fn chisq(locus_counts: &mut LocusCounts, filter_stats: &FilterStats) -> Option<String> {
    let locus_counts = match locus_counts.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    };
    // Marginal sums and total count
    let (n, p) = locus_counts.matrix.shape();
    let t = (n * p) as f64;
    let total = locus_counts.matrix.sum() as f64;
    let row_sums = locus_counts.matrix.column_sum();
    let col_sums = locus_counts.matrix.row_sum();
    // Caculate the chi-square test statistic
    let mut denominator: f64 = 0.0;
    let mut observed: f64;
    let mut expected: f64;
    for i in 0..n {
        for j in 0..p {
            observed = locus_counts.matrix[(i,j)] as f64;
            expected = (row_sums[i] * col_sums[j]) as f64 / total;
            denominator += f64::powf( observed - expected, 2.0);
        }
    }
    let chi2 = denominator / t;
    // Calculate the significance, i.e. how far away are we from the expected mean of the chi-squared distribution?
    let d = ChiSquared::new(t - 1.0).unwrap();
    let pval = d.cdf(chi2);
    // Output line
    let out = vec![locus_counts.chromosome.clone(),
                           locus_counts.position.to_string(),
                           locus_counts.alleles_vector.join(""),
                           chi2.to_string(),
                           pval.to_string()]
                      .join(",") + "\n";

    Some(out)
}
