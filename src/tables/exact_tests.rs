use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DVector, DMatrix, DMatrixView, U3, U4};

use crate::io::sync::Sync;
use crate::io::sync::AlleleCountsOrFrequencies;

pub fn fisher(vec_acf: &mut Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
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

pub fn barnard(vec_acf: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
    println!("CONVERT TO MATRIX");
    let (c, p, a, x) = vec_acf.convert_to_matrix(true).unwrap();
    println!("X: {:?}", x);
    Ok(0)
}

pub fn boschloo(vec_acf: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>) -> io::Result<i32> {
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