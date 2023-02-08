use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DVector, DMatrix, DMatrixView, U3, U4};

use crate::io::sync::Sync;
use crate::io::sync::AlleleCountsOrFrequencies;

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

fn factorial(x: u128) -> io::Result<u128> {
    match x {
        0 => Ok(1),
        1 => Ok(1),
        _ => Ok(factorial(x - 1).unwrap() * x)
    }
}

fn fisher_base(X: &DMatrix<f64>) -> io::Result<f64> {
    let (n, m) = X.shape();
    let mut C = DMatrix::from_element(n, m, 0 as u128);
    
    // Find the minimum frequency to get the maximum natural number restricted by u128, i.e. n=34 if n! > u128::MAX
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
            C[(i, j)] = (s as f64 * X[(i, j)]).ceil() as u128;
        }
    }
    println!("COUNTS: {:?}", C);

    let row_sums = X.row_sum();
    let col_sums = X.column_sum();
    

    Ok(0.0)
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