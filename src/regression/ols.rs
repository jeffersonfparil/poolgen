use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use nalgebra::{self, DMatrix, DVector};
use std::sync::{Arc, Mutex};
use std::fs::{File, OpenOptions};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::io::sync::{Sync, sync_and_pheno_analyser_and_writer_single_thread};
use crate::io::sync::AlleleCountsOrFrequencies;
use crate::io::phen::{Phenotypes, load_phen};

use statrs::distribution::{Normal, ContinuousCDF};

fn ols(X: &DMatrix<f64>, Y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, DMatrix<f64>)> {
    let (n, p) = X.shape();
    let (n_, m) = Y.shape();
    // println!("x={:?}; y={:?}", x, y);
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let beta: DMatrix<f64>;
    let var_cov: DMatrix<f64>;
    if n < p {
        let Xt = X.transpose();
        let inv_XXt = (X * &Xt).try_inverse().unwrap();
        beta = &Xt * &inv_XXt * Y;
        let error = Y - (X * &beta);
        let se = (&error.transpose() * &error) / (n as f64 - p as f64);
        var_cov = &Xt * &inv_XXt * &inv_XXt * X  * &se;
    } else {
        let Xt = X.transpose();
        let inv_XtX = (&Xt * X).try_inverse().unwrap();
        beta = &inv_XtX * &Xt * Y;
        let error = Y - (X * &beta);
        let se = (&error.transpose() * &error) / (n as f64 - p as f64);
        var_cov = &inv_XtX * &se;
    }
    Ok((beta, var_cov))
}

pub fn ols_iterate_base(acf: &mut AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>, phen: &Phenotypes<f64, nalgebra::Dyn, nalgebra::Dyn>) -> Option<String> {
    acf.counts_to_frequencies().unwrap();
    let idx = acf.coordinate.clone();
    let chr = acf.chromosome.clone();
    let pos = acf.position.clone();
    let ale = acf.alleles_vector.clone();
    let mut X = acf.matrix.clone();
    let nam = phen.name.clone();
    let Y = phen.phen.clone();
    // println!("ACF={:?}", acf);
    // println!("Y={:?}", Y);
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

    let (beta, var_cov) = ols(&X, &Y).unwrap();
    
    // Iterate across alleles
    let first_3_col = vec![chr, pos.to_string(), ale.join("")];
    let mut line: Vec<String> = vec![];
    for i in 0..k {
        line.append(&mut first_3_col.clone());
        line.push("Pheno_".to_string() + &(i.to_string())[..]);
        line.push(beta.column(i).iter().map(|x| x.to_string()).collect::<Vec<String>>().join(","));
        // line.push(pval.to_string() + "\n");
    }
    let out = line.join(",").replace("\n,", "\n");


    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

pub fn ols_iterate(fname: &String, phen_fname: &String, delim: &String, header: &bool, name_col: &usize, phen_col: &Vec<usize>, out: &String, n_threads: &u64) -> io::Result<String> {
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
            let vec_out_per_thread = sync_and_pheno_analyser_and_writer_single_thread(&fname_clone, &format_clone, &n_pools_clone, &start, &end, &n_digits_clone, &phen_clone, ols_iterate_base).unwrap();
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
