use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use nalgebra::{self, DMatrix};
use std::sync::{Arc, Mutex};
use std::fs::{File, OpenOptions};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::sync::{Sync, sync_analyser_and_writer_single_thread};
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
    // println!("PROD_FAC_MAR_SUMS: {:?}", log_prod_fac_marginal_sums);
    // println!("PROD_FAC_SUMS: {:?}", prod_fac_sums);
    // println!("OUT: {:?}", p);
    Ok(p)
}

pub fn fisher_base(acf: &mut AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>) -> Option<String> {
    acf.counts_to_frequencies().unwrap();
    let idx = acf.coordinate.clone();
    let chr = acf.chromosome.clone();
    let pos = acf.position.clone();
    let ale = acf.alleles_vector.clone().join("");
    let X = acf.matrix.clone();

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
        return None
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
    let out = vec![chr,
                           pos.to_string(),
                           ale,
                           (p_observed + p_extremes).to_string()]
                      .join(",") + "\n";

    Some(out)
}

pub fn fisher(fname: &String, out: &String, n_threads: &u64) -> io::Result<String> {
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                        .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                        .join(".");
        out = bname + "-Fisher_exact_test-" + &time.to_string() + ".csv";
    }
    // Instatiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &out;
    let mut file_out = OpenOptions::new().create_new(true)
                                               .write(true)
                                               .append(false)
                                               .open(&out)
                                               .expect(&error_writing_file);
    let chunks = crate::io::find_file_splits(fname, n_threads).unwrap();
    let n_digits = chunks[*n_threads as usize].to_string().len();
    println!("Chunks: {:?}", chunks);

    // Determine the format of the input file as well as the number of pools
    let file = File::open(fname).unwrap();
    let mut reader = BufReader::new(file);
    let mut format: String = "sync".to_owned();
    let mut caught_1_line = false;
    let mut n_pools: usize = 0;
    while caught_1_line == false {
        let mut line = String::new();
        let _ = reader.read_line(&mut line).unwrap();
        if line.as_bytes()[0] == 35 as u8 {
            continue;
        } else {
            let vec_line = line.split("\t").collect::<Vec<&str>>();
            let allele_column = vec_line[3]
                                                .split(":")
                                                .collect::<Vec<&str>>()
                                                .into_iter().map(|a| a.to_owned())
                                                .collect::<Vec<String>>();
            // println!("{:?}", allele_column[0]);
            // println!("{:?}", allele_column[0].split("|").collect::<Vec<&str>>()[1]);
            format = match allele_column[0].parse::<i64>() {
                Ok(_) => "sync".to_owned(),
                Err(_) => match allele_column[0].split("|").collect::<Vec<&str>>()[1].parse::<f64>() {
                    Ok(_) => "syncx".to_owned(),
                    Err(_) => return Err(Error::new(ErrorKind::Other, "Please check the format of the input file: ".to_owned() + fname + ". Please use a sync or syncx file.")),
                },
            };
            if format == "sync" {
                n_pools = vec_line.len() - 3;
            } else {
                let vec_line_parse = vec_line[3].split("|").collect::<Vec<&str>>();
                let freqs = vec_line_parse[1]
                                    .split(":")
                                    .collect::<Vec<&str>>()
                                    .into_iter()
                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                    .collect::<Vec<f64>>();
                n_pools = freqs.len();
            }
            // println!("VEC_LINE: {:?}", vec_line);
            // println!("N_POOLS: {:?}", n_pools);
            caught_1_line = true;
        }
    }

    // Instantiate thread object for parallel execution
    let mut thread_objects = Vec::new();
    // Vector holding all returns from read_chunk()
    let mut thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..(*n_threads as usize) {
        // Clone read_chunk parameters
        let fname_clone = fname.clone();
        let format_clone = format.clone();
        let start = chunks[i].clone();
        let end = chunks[i+1].clone();
        let n_pools_clone = n_pools.clone();
        let n_digits_clone = n_digits.clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let vec_out_per_thread = sync_analyser_and_writer_single_thread(&fname_clone, &format_clone, &n_pools_clone, &start, &end, &n_digits_clone, fisher_base).unwrap();
            thread_ouputs_clone.lock().unwrap().push(vec_out_per_thread);
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }
    // Extract output filenames from each thread into a vector and sort them
    let mut fnames_out: Vec<String> = Vec::new();
    for f in thread_ouputs.lock().unwrap().iter() {
        fnames_out.push(f.to_owned());
    }
    fnames_out.sort();
    // println!("{:?}", fnames_out);
    // Add header
    let header = "#chr,pos,alleles,pvalue\n".to_owned();
    file_out.write_all(header.as_bytes()).unwrap();
    // Iterate across output files from each thread, and concatenate non-empty files
    for f in fnames_out {
        let mut file: File = File::open(&f).unwrap();
        if file.metadata().unwrap().len() == 0 {
        } else {
            io::copy(&mut file, &mut file_out).unwrap();
            // println!("{:?}", f);
        }
        // Clean-up: remove temporary output files from each thread
        let error_deleting_file = "Unable to remove file: ".to_owned() + &f;
        std::fs::remove_file(f).expect(&error_deleting_file);
    }
    
    Ok(out)
}
