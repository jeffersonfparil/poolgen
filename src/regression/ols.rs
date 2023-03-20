use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix, DVector};
use crate::base::*;
use crate::io::sync::AlleleCountsOrFrequencies;
use crate::io::phen::{Phenotypes, load_phen};

use statrs::distribution::{StudentsT, ContinuousCDF};

fn ols(X: &DMatrix<f64>, Y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, DMatrix<f64>)> {
    let (n, p) = X.shape();
    let (n_, k) = Y.shape();
    // println!("X={:?}; Y={:?}", X, Y);
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let mut b: DVector<f64>;
    let mut C: DMatrix<f64>;
    let mut beta = DMatrix::from_element(p, k, 0.0);
    let mut var_beta = DMatrix::from_element(p, k, 0.0);
    for j in 0..k {
        let y = Y.column(j);
        if n < p {
            let Xt = X.transpose();
            let inv_XXt = match (X * &Xt).try_inverse() {
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible X")),
            };
            if inv_XXt.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible X"))
            }
            b = &Xt * &inv_XXt * y;
            C = &Xt * &inv_XXt * &inv_XXt * X;
        } else {
            let Xt = X.transpose();
            let inv_XtX = match (&Xt * X).try_inverse(){
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible X")),
            };
            if inv_XtX.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible X"))
            }
            b = &inv_XtX * &Xt * y;
            C = inv_XtX;
        }
        let e = y - (X * &b);
        let se = (&e.transpose() * &e).sum() / (n as f64 - p as f64);
        let vb = se * &C;
        for i in 0..p {
            beta[(i,j)] = b[i];
            var_beta[(i,j)] = vb[(i, i)];
        }
        // println!("#################################");
        // println!("b={:?}", b);
        // println!("C={:?}", C);
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
    let mut X = locus_frequencies.matrix.clone();
    let Y = locus_counts_and_phenotypes.phenotypes.clone();
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, mut p) =  X.shape();
    let (m, k) = Y.shape();
    if n != m {
        return None
    }
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if p >= 2 {
        X = X.clone().remove_columns(p-1, 1);
        p -= 1;
    }
    X = X.clone().insert_column(0, 1.0);
    p += 1;
    // OLS and compute the p-values associated with each estimate
    let (beta, var_beta) = match ols(&X, &Y) {
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
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            line.push(beta[(i,j)].to_string());
            line.push(pval[(i,j)].to_string() + "\n");
        }
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}
