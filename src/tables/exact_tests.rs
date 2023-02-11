use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DVector, DMatrix, DMatrixView, U3, U4};

use crate::io::sync::Sync;
use crate::io::sync::AlleleCountsOrFrequencies;

fn factorial_log10(x: f64) -> io::Result<f64> {
    let mut out: f64 = 0.0;
    for i in 1..x as usize {
        out = out + f64::log10(i as f64);
    }
    Ok(out)
}

fn hypergeom_ratio(C: &DMatrix<f64>, log_prod_fac_marginal_sums: &f64) -> io::Result<f64> {
    // Log-Product of counts
    let mut prod_fac_sums = 1 as f64;
    for i in C.iter() {
        prod_fac_sums = prod_fac_sums + factorial_log10(*i).unwrap();
    }
    prod_fac_sums = prod_fac_sums + factorial_log10(C.sum()).unwrap();
    // Calculate the p-value
    let p = f64::powf(10.0, log_prod_fac_marginal_sums - prod_fac_sums);
    println!("PROD_FAC_MAR_SUMS: {:?}", log_prod_fac_marginal_sums);
    println!("PROD_FAC_SUMS: {:?}", prod_fac_sums);
    println!("OUT: {:?}", p);
    Ok(p)
}

fn fisher_base(X: &DMatrix<f64>) -> io::Result<f64> {
    // Instatiate the counts matrix
    let (n, m) = X.shape();
    let mut C = DMatrix::from_element(n, m, 0 as f64);
    // Find the minimum frequency to get the maximum natural number restricted by f64, i.e. n=34 if n! > f64::MAX
    let mut x = X.iter()
                                            .filter(|a| *a > &0.0)
                                            .into_iter()
                                            .map(|a| a.to_owned())
                                            .collect::<Vec<f64>>();
    x.sort_by(|a, b| a.partial_cmp(&b).unwrap());
    let mut s = (1.0 / x[0]).ceil() as usize;
    s = match s < 34  {
        true => s,
        false => 34,
    };
    // Populate the counts matrix
    for i in 0..n {
        for j in 0..m {
            C[(i, j)] = (s as f64 * X[(i, j)]).ceil() as f64;
        }
    }
    println!("COUNTS: {:?}", C);
    // Log-Product of the marginal sums (where C.row_sum() correspond to the the column marginal sums and vice versa)
    let row_sums = C.column_sum().clone_owned();
    let col_sums = C.row_sum().clone_owned();
    let mut log_prod_fac_marginal_sums = 1 as f64;
    for r in row_sums.iter() {
        log_prod_fac_marginal_sums = log_prod_fac_marginal_sums + factorial_log10(*r).unwrap();
    }
    for c in col_sums.iter() {
        log_prod_fac_marginal_sums = log_prod_fac_marginal_sums + factorial_log10(*c).unwrap();
    }
    // Define the observed hypergeometric ratio, i.e. p of the observed data
    let p_observed = hypergeom_ratio(&C, &log_prod_fac_marginal_sums).unwrap();
    // Iterate across all possible combinations of counts with the same marginal sums
    println!("C shape: {:?}", C.shape());
    println!("row_sums shape: {:?}", row_sums.shape());
    println!("col_sums shape: {:?}", col_sums.shape());
    println!("n: {:?}", n);
    println!("m: {:?}", m);

    let mut m_: usize;
    let mut p: f64;
    for max_i in 0..n
        for max_j in 0..m {
            for i in 0..n {
                for j in 0..m {
                    let max: usize = vec![(row_sums[(i, 0)] - C.index((i, 0..j)).sum().ceil()) as usize,
                                          (col_sums[(0, j)] - C.index((0..i, j)).sum().ceil()) as usize]
                                          .into_iter().min().unwrap();
                    println!("i: {:?}", i);
                    println!("j: {:?}", j);
                    println!("m_: {:?}", max);
                    if (max_i != i) & (max_j != j) {
                        if i < max_i {
                            C[(i,j)] = 0;
                        } else {
                            C[(i,j)] = max;
                        }
                    }
                }
            }
            // Make sure we kept the marginal sums constant
            assert!(row_sums == C.column_sum().clone_owned());
            asset!(col_sums == C.row_sum().clone_owned());
            // Calculate hypergeometric ratio
            p = hypergeom_ratio(&C, &log_prod_fac_marginal_sums).unwrap();
            println!("New C: {:?}", C);
            println!("row_sums: {:?}", row_sums);
            println!("col_sums: {:?}", col_sums);
            println!("p: {:?}", p);
        }
    }

    Ok(p_observed)
}

pub fn fisher(vec_acf: &mut Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
    vec_acf.counts_to_frequencies().unwrap();
    for acf in vec_acf {
        if acf.matrix.shape().1 > 1 {
            let X = acf.matrix.clone();
            println!("X: {:?}", X);
            println!("FISHER_BASE: {:?}", fisher_base(&X));
            println!("NALLELES: {:?}", X.shape().1);
            println!("ROW 0: {:?}", X.row(0).clone_owned());
        }

    }


    Ok(0)
}

pub fn barnard(vec_acf: &mut Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
    println!("CONVERT TO MATRIX");
    println!("X0: {:?}", vec_acf[0]);
    println!("CONVERT TO COUNTS TO FREQS");
    vec_acf.counts_to_frequencies().unwrap();
    println!("X1: {:?}", vec_acf[0]);
    let (c, p, a, x) = vec_acf.convert_to_matrix(false).unwrap();
    let x_tmp = x.view((0,0), (2,2));
    println!("MATRIX: {:?}", x);
    println!("MATRIX: {:?}", x);
    println!("CHROM: {:?}", &c[0..10]);
    println!("POS: {:?}", &p[0..10]);
    println!("ALLELE: {:?}", &a[0..10]);

    Ok(0)
}

pub fn boschloo(vec_acf: &mut Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
    println!("CONVERT TO MATRIX");
    let (c, p, a, x) = vec_acf.convert_to_matrix(true).unwrap();
    println!("X: {:?}", x);

    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[0], p[0], a[0]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[1], p[1], a[1]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[2], p[2], a[2]);
    
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[2000], p[2000], a[2000]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[2001], p[2001], a[2001]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[2002], p[2002], a[2002]);

    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[5000], p[5000], a[5000]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[5001], p[5001], a[5001]);
    // println!("chr: {:?}; pos: {:?}; allele: {:?}", c[5002], p[5002], a[5002]);

    // println!("CONVERT TO COUNTS TO FREQS");
    // println!("X0: {:?}", vec_acf[0]);
    // vec_acf.counts_to_frequencies();
    // println!("X1: {:?}", vec_acf[0]);

    // println!("FILTER BY MAF");
    // println!("X0: {:?}", vec_acf[0]);
    // println!("X0: n=: {:?}", vec_acf.len());
    // vec_acf.filter(0.01);
    // println!("X1: {:?}", vec_acf[0]);
    // println!("X1: n=: {:?}", vec_acf.len());

    Ok(0)
}