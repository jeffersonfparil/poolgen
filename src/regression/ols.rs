use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix, DVector};
use crate::base::*;

use statrs::distribution::{StudentsT, ContinuousCDF};

fn ols(x_matrix: &DMatrix<f64>, y_matrix: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, DMatrix<f64>)> {
    let (n, p) = x_matrix.shape();
    let (n_, k) = y_matrix.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let mut b: DVector<f64>;
    let mut c_matrix: DMatrix<f64>;
    let mut beta = DMatrix::from_element(p, k, 0.0);
    let mut var_beta = DMatrix::from_element(p, k, 0.0);
    for j in 0..k {
        let y_vector = y_matrix.column(j);
        if n < p {
            let xt_matrix = x_matrix.transpose();
            let inv_xxt_matrix = match (x_matrix * &xt_matrix).try_inverse() {
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xxt_matrix.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"))
            }
            b = &xt_matrix * &inv_xxt_matrix * y_vector;
            c_matrix = &xt_matrix * &inv_xxt_matrix * &inv_xxt_matrix * x_matrix;
        } else {
            let xt_matrix = x_matrix.transpose();
            let inv_xtx_matrix = match (&xt_matrix * x_matrix).try_inverse(){
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xtx_matrix.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"))
            }
            b = &inv_xtx_matrix * &xt_matrix * y_vector;
            c_matrix = inv_xtx_matrix;
        }
        let e = y_vector - (x_matrix * &b);
        let se = (&e.transpose() * &e).sum() / (n as f64 - p as f64);
        let vb = se * &c_matrix;
        for i in 0..p {
            beta[(i,j)] = b[i];
            var_beta[(i,j)] = vb[(i, i)];
        }
        // println!("#################################");
        // println!("b={:?}", b);
        // println!("c_matrix={:?}", c_matrix);
        // println!("e={:?}", e);
        // println!("se={:?}", se);
        // println!("vb={:?}", vb);
        // println!("beta={:?}", beta);
        // println!("var_beta={:?}", var_beta);
    }
    Ok((beta, var_beta))
}

pub fn ols_iterate(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
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
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, mut p) =  x_matrix.shape();
    let (m, k) = y_matrix.shape();
    if n != m {
        return None
    }
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if p >= 2 {
        x_matrix = x_matrix.clone().remove_columns(p-1, 1);
        p -= 1;
    }
    x_matrix = x_matrix.clone().insert_column(0, 1.0);
    p += 1;
    // OLS and compute the p-values associated with each estimate
    let (beta, var_beta) = match ols(&x_matrix, &y_matrix) {
        Ok(x) => x,
        Err(_) => return None,
    };
    let d = StudentsT::new(0.0, 1.0, p as f64 - 1.0).unwrap();
    let mut t: f64;
    let mut pval = DMatrix::from_element(p, k, 0.0);
    for i in 0..p {
        for j in 0..k {
            t = beta[(i, j)] / var_beta[(i, j)];
            if t.is_infinite() {
                pval[(i,j)] = 0.0
            } else if t.is_nan() {
                pval[(i,j)] = 1.0
            } else {
                pval[(i, j)] = 2.00 * (1.00 - d.cdf(t.abs()));
            }
        }
    }

    // Iterate across alleles
    let first_2_col = vec![locus_frequencies.chromosome, locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    for i in 1..p {
        // excluding the intercept
        for j in 0..k {
            line.append(&mut first_2_col.clone());
            line.push(locus_frequencies.alleles_vector[i-1].clone());
            line.push(x_matrix.column(i).mean().to_string());
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            line.push(beta[(i,j)].to_string());
            line.push(pval[(i,j)].to_string() + "\n");
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
    fn test_ols() {
        // Expected
        let expected_output1: DMatrix<f64> =  DMatrix::from_column_slice(3, 1, &[-3.00,  0.92, 10.26]);
        let expected_output2: DMatrix<f64> =  DMatrix::from_column_slice(3, 1, &[12.74, 31.37, 65.35]);
        // Inputs
        let x: DMatrix<f64> = DMatrix::from_column_slice(5, 3, &[1.0, 1.0, 1.0, 1.0, 1.0,
                                                                                    0.1, 0.2, 0.3, 0.4, 0.5,
                                                                                    0.5, 0.3, 0.5, 0.6, 0.5]);
        let y: DMatrix<f64> = DMatrix::from_column_slice(5, 2, &[2.0, 1.0, 1.0, 5.0, 2.0,
                                                                                    0.1, 0.2, 0.1, 0.2, 0.0]);
        // let counts: DMatrix<u64> = DMatrix::from_row_slice(5, 2, &[1, 9,
        //                                                                               2, 8, 
        //                                                                               3, 7,
        //                                                                               4, 6,
        //                                                                               5, 5]);
        // let filter_stats = FilterStats{remove_ns: true, min_quality: 0.005, min_coverage: 1, min_allele_frequency: 0.005, pool_sizes: vec![0.2,0.2,0.2,0.2,0.2]};
        // let locus_counts = LocusCounts{chromosome: "Chromosome1".to_owned(), position: 12345, alleles_vector: vec!["A".to_owned(), "T".to_owned()], matrix:counts};
        // let phenotypes: DMatrix<f64> = y.clone();
        // let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes{locus_counts: locus_counts, phenotypes: phenotypes, pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"].into_iter().map(|x| x.to_owned()).collect::<Vec<String>>()};
        // Outputs
        let (beta, var_beta) = ols(&x, &y).unwrap();
        // Assertions
        assert_eq!((expected_output1[(0,0)]*100.0).round(), (beta[(0,0)]*100.0).round());
        assert_eq!((expected_output1[(1,0)]*100.0).round(), (beta[(1,0)]*100.0).round());
        assert_eq!((expected_output1[(2,0)]*100.0).round(), (beta[(2,0)]*100.0).round());

        assert_eq!((expected_output2[(0,0)]*100.0).round(), (var_beta[(0,0)]*100.0).round());
        assert_eq!((expected_output2[(1,0)]*100.0).round(), (var_beta[(1,0)]*100.0).round());
        assert_eq!((expected_output2[(2,0)]*100.0).round(), (var_beta[(2,0)]*100.0).round());
    }
}
