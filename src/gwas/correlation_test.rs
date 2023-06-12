use crate::base::*;
use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use statrs::distribution::{ContinuousCDF, StudentsT};

pub fn pearsons_correlation(
    x: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
    y: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
) -> io::Result<(f64, f64)> {
    let n = x.len();
    if n != y.len() {
        return Err(Error::new(
            ErrorKind::Other,
            "Input vectors are not the same size.",
        ));
    }
    let mu_x = x.mean().unwrap();
    let mu_y = y.mean().unwrap();
    let x_less_mu_x = x.map(|x| x - mu_x);
    let y_less_mu_y = y.map(|y| y - mu_y);
    let x_less_mu_x_squared = x_less_mu_x.map(|x| x.powf(2.0));
    let y_less_mu_y_squared = y_less_mu_y.map(|y| y.powf(2.0));
    let numerator = (x_less_mu_x * y_less_mu_y).sum();
    let denominator = x_less_mu_x_squared.sum().sqrt() * y_less_mu_y_squared.sum().sqrt();
    let r_tmp = numerator / denominator;
    let r = match r_tmp.is_nan() {
        true => 0.0,
        false => r_tmp,
    };
    let sigma_r_denominator = (1.0 - r.powf(2.0)) / (n as f64 - 2.0);
    if sigma_r_denominator <= 0.0 {
        // Essentially no variance in r2, hence very significant
        return Ok((r, f64::EPSILON));
    }
    let sigma_r = sigma_r_denominator.sqrt();
    let t = r / sigma_r;
    let d = StudentsT::new(0.0, 1.0, n as f64 - 2.0).unwrap();
    let pval = 2.00 * (1.00 - d.cdf(t.abs()));
    Ok((r, pval))
}

pub fn correlation(
    locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes,
    filter_stats: &FilterStats,
) -> Option<String> {
    // Check struct
    locus_counts_and_phenotypes.check().unwrap();
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
        .locus_counts
        .filter(filter_stats)
    {
        Ok(x) => x.clone(),
        Err(_) => return None,
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Extract the genotype and phenotypes
    let mut x_matrix = locus_frequencies.matrix.clone();
    let y_matrix = locus_counts_and_phenotypes.phenotypes.clone();
    // Keep p-1 alleles for conciseness
    let p = x_matrix.ncols();
    if p >= 2 {
        x_matrix.remove_index(Axis(1), p - 1);
    }
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let n = x_matrix.nrows();
    let p = x_matrix.ncols();
    let m = y_matrix.nrows();
    let k = y_matrix.ncols();
    if n != m {
        return None;
    }
    // Iterate across alleles
    let (mut corr, mut pval): (f64, f64);
    let first_2_col = vec![
        locus_frequencies.chromosome,
        locus_frequencies.position.to_string(),
    ];
    let mut line: Vec<String> = vec![];
    for i in 0..p {
        let x: ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>> = x_matrix.column(i);
        for j in 0..k {
            line.append(&mut first_2_col.clone());
            line.push(locus_frequencies.alleles_vector[i].clone());
            line.push(x.mean().unwrap().to_string());
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            let y: ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>> = y_matrix.column(j);
            (corr, pval) = pearsons_correlation(&x, &y).unwrap();
            line.push(parse_f64_roundup_and_own(corr, 6));
            line.push(pval.to_string() + "\n");
        }
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_correlation() {
        // Expected
        let expected_output1: f64 = 0.3849001794597505;
        let expected_output2: f64 = 0.5223146158470686;
        let expected_output3: String =
            "Chromosome1,12345,A,0.3,Pheno_0,0.3849,0.5223146158470686\n".to_owned();
        // Inputs
        let x_ = Array1::from_vec(vec![0.1, 0.2, 0.3, 0.4, 0.5]);
        let x = x_.view();
        let y_ = Array1::from_vec(vec![2.0, 1.0, 1.0, 5.0, 2.0]);
        let y = y_.view();
        let counts: Array2<u64> =
            Array2::from_shape_vec((5, 2), vec![1, 9, 2, 8, 3, 7, 4, 6, 5, 5]).unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned(), "T".to_owned()],
            matrix: counts,
        };
        let phenotypes: Array2<f64> = y_.clone().into_shape((y_.len(), 1)).unwrap();
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: phenotypes,
            pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        };
        // Outputs
        let (corr, pval) = pearsons_correlation(&x, &y).unwrap();
        let correlation_line =
            correlation(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        // Assertions
        assert_eq!(expected_output1, corr);
        assert_eq!(expected_output2, pval);
        assert_eq!(expected_output3, correlation_line);
    }
}
