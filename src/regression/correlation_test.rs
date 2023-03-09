use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use nalgebra::{self, DMatrix, DVector};
use std::sync::{Arc, Mutex};
use std::fs::{File, OpenOptions};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::sync::{Sync, sync_and_pheno_analyser_and_writer_single_thread};
use crate::io::sync::AlleleCountsOrFrequencies;
use crate::io::phen::{Phenotypes, load_phen};

use statrs::distribution::{StudentsT, ContinuousCDF};

fn pearsons_correlation(x: &DVector<f64>, y: &DVector<f64>) -> io::Result<(f64, f64)> {
    let n = x.len();
    // println!("x={:?}; y={:?}", x, y);
    if n != y.len() {
        return Err(Error::new(ErrorKind::Other, "Input vectors are not the same size."));
    }
    let mu_x = x.mean();
    let mu_y = y.mean();
    let x_less_mu_x = x.map(|x| x-mu_x);
    let y_less_mu_y = y.map(|y| y-mu_y);
    let x_less_mu_x_squared = x_less_mu_x.map(|x| x.powf(2.0));
    let y_less_mu_y_squared = y_less_mu_y.map(|y| y.powf(2.0));
    // println!("x_less_mu_x={:?}", x_less_mu_x);
    // println!("y_less_mu_y={:?}", y_less_mu_y);
    let numerator = x_less_mu_x.component_mul(&y_less_mu_y).sum();
    let denominator = x_less_mu_x_squared.sum().sqrt() * y_less_mu_y_squared.sum().sqrt();
    let r_tmp = numerator / denominator;
    let r = match r_tmp.is_nan() {
        true => 0.0,
        false => r_tmp,
    };
    // println!("numeratorr={:?}; demonitatorr={:?}; r={:?}", numerator, denominator, r);
    let sigma_r = ((1.0 - r.powf(2.0)) / (n as f64 - 2.0)).sqrt();
    let t = r / sigma_r;
    let d = StudentsT::new(0.0, 1.0, n as f64 - 1.0).unwrap();
    let pval = 1.00 - d.cdf(t.abs());
    Ok((r, pval))
}

pub fn correlation_base(acf: &mut AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>, phen: &Phenotypes<f64, nalgebra::Dyn, nalgebra::Dyn>) -> Option<String> {
    acf.counts_to_frequencies().unwrap();
    let idx = acf.coordinate.clone();
    let chr = acf.chromosome.clone();
    let pos = acf.position.clone();
    let ale = acf.alleles_vector.clone();
    let X = acf.matrix.clone();
    let nam = phen.name.clone();
    let Y = phen.phen.clone();
    // println!("ACF={:?}", acf);
    // println!("Y={:?}", Y);
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, p) =  X.shape();
    let (m, k) = Y.shape();
    if n != m {
        return None
    }

    // Iterate across alleles
    let (mut corr, mut pval): (f64, f64);
    let first_2_col = vec![chr, pos.to_string()];
    // println!("first_2_col: {:?}", first_2_col);
    // println!("ale: {:?}", ale);
    // println!("p: {:?}", p);
    let mut line: Vec<String> = vec![];
    for i in 0..p {
        let x = DVector::from(X.column(i));
        for j in 0..k {
            line.append(&mut first_2_col.clone());
            line.push(ale[i].clone());
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            let y  = DVector::from(Y.column(j));
            // println!("LINE: {:?}", line);
            // println!("x={:?}; y={:?}", x, y);
            (corr, pval) = pearsons_correlation(&x, &y).unwrap();
            line.push(corr.to_string());
            line.push(pval.to_string() + "\n");
        }
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

pub fn correlation(fname: &String, phen_fname: &String, delim: &String, header: &bool, name_col: &usize, phen_col: &Vec<usize>, out: &String, n_threads: &u64) -> io::Result<String> {
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                        .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                        .join(".");
        out = bname + "-Pearsons_correlation_test-" + &time.to_string() + ".csv";
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

    // Load the phnoetypes
    let phen = load_phen(phen_fname, delim, header, name_col, phen_col).unwrap();
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
        let phen_clone = phen.clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let vec_out_per_thread = sync_and_pheno_analyser_and_writer_single_thread(&fname_clone, &format_clone, &n_pools_clone, &start, &end, &n_digits_clone, &phen_clone, correlation_base).unwrap();
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
    let header = "#chr,pos,allele,Pearsons_correlation,pvalue\n".to_owned();
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
