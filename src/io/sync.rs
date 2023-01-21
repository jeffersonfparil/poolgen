use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader};
use std::time::{SystemTime, UNIX_EPOCH};

use nalgebra::{DVector, DMatrix};

pub fn load(fname: &String) -> io::Result<i64> {
    // let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    // let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
    //                  .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
    //                  .join(".");
    // let out = bname + "-" + &time.to_string() + "." + file_format;

    // Test linear algebra
    let n: usize = 10;
    let p: usize = 5;
    let row_slice: Vec<f64> = (0..(n as i32 * p as i32)).map(|x| x as f64).collect();
    let matrix = DMatrix::from_row_slice(n, p, &row_slice);
    let vector = DVector::from_column_slice(&vec![1.0, 0.43, 0.3, 0.3456, 0.134]);
    let b = &matrix * &vector;
    println!("matrix variance: {:?}", matrix.variance());
    println!("vector variance: {:?}", vector.variance());
    println!("matrix multiplication product: {:?}", b);
    Ok(0)
}