use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DVector, DMatrix, DMatrixView, U3, U4};
use std::sync::{Arc, Mutex};
use std::fs::File;
use std::io::{prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::sync::Sync;
use crate::io::sync::load;
use crate::io::sync::AlleleCountsOrFrequencies;

#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct SummaryStatistics {
    coordinate: String,
    chromosome: String,
    position: u64,
    alleles: String,
    statistic: f64,
}

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
    // println!("PROD_FAC_MAR_SUMS: {:?}", log_prod_fac_marginal_sums);
    // println!("PROD_FAC_SUMS: {:?}", prod_fac_sums);
    // println!("OUT: {:?}", p);
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
    // Return error if we have no coverage
    if x.len() == 0 {
        return Err(Error::new(ErrorKind::Other, "No coverage"))
    }
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
    // println!("COUNTS: {:?}", C);
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
    // println!("C shape: {:?}", C.shape());
    // println!("row_sums: {:?}", row_sums);
    // println!("col_sums: {:?}", col_sums);
    // println!("n: {:?}", n);
    // println!("m: {:?}", m);
    // println!("###################################################");

    let mut p_extremes: f64 = 0.0;
    let mut max: f64;
    for max_i in 0..n {
        for max_j in 0..m {
            for i in 0..n {
                for j in 0..m {
                    max = vec![(row_sums[(i, 0)] - C.index((i, 0..j)).sum().ceil()) as usize,
                               (col_sums[(0, j)] - C.index((0..i, j)).sum().ceil()) as usize]
                               .into_iter().min().unwrap() as f64;
                    if (i==(n-1)) | (j==(m-1)) {
                        C[(i,j)] = max;
                    } else if (i < max_i) | (j < max_j) {
                        C[(i,j)] = 0.0;
                    } else {
                        C[(i,j)] = max;
                    }
                    // println!("i={:?}; max_i={:?}", i, max_i);
                    // println!("j={:?}; max_j={:?}", j, max_j);
                    // println!("max={:?}; C[(i,j)]={:?}", max, C[(i,j)]);
                }
            }
            // println!("NEW C: {:?}", C);
            let mut j: usize;
            let mut i: usize;
            for inv_j in 0..m {
                for inv_i in 0..n {
                    j = m - (inv_j+1);
                    i = n - (inv_i+1);
                    max = vec![(row_sums[(i, 0)] - C.index((i, 0..m)).sum().ceil()) as usize,
                               (col_sums[(0, j)] - C.index((0..n, j)).sum().ceil()) as usize]
                               .into_iter().min().unwrap() as f64;
                    if max > 0.0 {
                        C[(i,j)] = max;
                    }
                    // println!("i={:?}; max_i={:?}", i, max_i);
                    // println!("j={:?}; max_j={:?}", j, max_j);
                    // println!("max={:?}; C[(i,j)]={:?}", max, C[(i,j)]);
                }
            }
            // Make sure we kept the marginal sums constant
            // println!("NEW C: {:?}", C);
            // println!("row_sum: {:?}", C.column_sum().clone_owned());
            // println!("col_sum: {:?}", C.row_sum().clone_owned());
            assert!(row_sums == C.column_sum().clone_owned());
            assert!(col_sums == C.row_sum().clone_owned());
            // Calculate hypergeometric ratio
            p_extremes += hypergeom_ratio(&C, &log_prod_fac_marginal_sums).unwrap();
            // println!("New C: {:?}", C);
            // println!("row_sums: {:?}", row_sums);
            // println!("col_sums: {:?}", col_sums);
            // println!("p_out: {:?}", p_extremes);
        }
    }
    Ok(p_observed + p_extremes)
}

pub fn fisher(fname: &String, out: &String, n_threads: &u64) -> io::Result<String> {
    let mut vec_acf = load(fname, n_threads).unwrap();
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                        .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                        .join(".");
        out = bname + "-Fisher_exact_test-" + &time.to_string() + ".csv";
    }

    vec_acf.counts_to_frequencies().unwrap();
    let n = vec_acf.len();
    // println!("VEC_ACF: {:?}; len={:?}", vec_acf, n);
    // for acf in vec_acf {
    //     if acf.matrix.shape().1 > 1 {
    //         let X = acf.matrix.clone();
    //         println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
    //         // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
    //         // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
    //         println!("chr-pos: {:?}-{:?}", acf.chromosome, acf.position);
    //         println!("X: {:?}", X);
    //         println!("FISHER_BASE: {:?}", fisher_base(&X));
    //         println!("NALLELES: {:?}", X.shape().1);
    //         println!("ROW 0: {:?}", X.row(0).clone_owned());
    //     }

    // }
    // Instantiate thread object for parallel execution
    let mut thread_objects = Vec::new();
    // Vector holding all returns from read_chunk()
    let mut thread_ouputs: Arc<Mutex<Vec<SummaryStatistics>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..n {
        // Clone read_chunk parameters
        let idx = vec_acf[i].coordinate.clone();
        let chr = vec_acf[i].chromosome.clone();
        let pos = vec_acf[i].position.clone();
        let ale = vec_acf[i].alleles_vector.clone().join("");
        let X = vec_acf[i].matrix.clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            // let vec_out_per_thread = 0.01;
            // println!("OUT: {:?}", fisher_base(&X).unwrap());
            // let vec_out_per_thread = fisher_base(&X).unwrap();
            let vec_out_per_thread = match fisher_base(&X) {
                Ok(x) => x,
                Err(_) => 1.0,
            };
            // println!("OUT: {:?}", vec_out_per_thread);
            thread_ouputs_clone.lock().unwrap().push(SummaryStatistics{coordinate: idx,
                                                                       chromosome: chr,
                                                                       position: pos,
                                                                       alleles: ale,
                                                                       statistic: vec_out_per_thread});
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }

    // Trying to sort the output
    let mut p: Vec<SummaryStatistics> = Vec::new();
    for i in thread_ouputs.lock().unwrap().iter() {
        p.push(i.clone());
    }
    p.sort_by(|a, b| a.coordinate.partial_cmp(&b.coordinate).unwrap());

// Instatiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &out;
    let mut file_out = File::create(&out).expect(&error_writing_file);
    file_out.write_all("chr,pos,alleles,Fisher_exact_test_pval\n".to_owned().as_bytes()).unwrap();
    for v in p.into_iter() {
        let mut line = vec![v.chromosome,
                                    v.position.to_string(),
                                    v.alleles,
                                    v.statistic.to_string()]
                                .join(",") + "\n";
        file_out.write_all(line.to_owned().as_bytes()).unwrap();
    }
    Ok(out)
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