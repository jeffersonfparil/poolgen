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

fn chisq_base(X: &DMatrix<f64>) -> io::Result<f64> {
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
    Ok(0.0)
}

pub fn chisq(fname: &String, out: &String, n_threads: &u64) -> io::Result<String> {
    let mut vec_acf = load(fname, n_threads).unwrap();
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                        .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                        .join(".");
        out = bname + "-Chisq_exact_test-" + &time.to_string() + ".csv";
    }

    vec_acf.counts_to_frequencies().unwrap();
    let n = vec_acf.len();
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
            let vec_out_per_thread = match chisq_base(&X) {
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
    file_out.write_all("chr,pos,alleles,Chisq_test_pval\n".to_owned().as_bytes()).unwrap();
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
