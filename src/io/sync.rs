use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader, BufWriter};
use std::time::{SystemTime, UNIX_EPOCH};
use std::io::{Error, ErrorKind};

use nalgebra::{DVector, DMatrix};

// File splitting for thread allocation for parallele computation
fn find_start_of_next_line(fname: &String, pos: u64) -> u64 {
    let mut out = pos.clone();
    if out > 0 {
        let mut file = File::open(fname).unwrap();
        let _ = file.seek(SeekFrom::Start(out));
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let _ = reader.read_line(&mut line).unwrap();
        out = reader.seek(SeekFrom::Current(0)).unwrap();
    }
    return out
}

fn find_file_splits(fname: &String, n_threads: &u64) -> Vec<u64> {
    let mut file = File::open(fname).unwrap();
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.seek(SeekFrom::Current(0)).unwrap();
    let mut out = (0..end).step_by((end/n_threads) as usize).collect::<Vec<u64>>();
    out.push(end);
    // println!("{:?}", end);
    // println!("{:?}", out);
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
    }
    out.dedup();
    // println!("{:?}", out);
    return out
}


fn syn2sync2_per_chunk (fname: &String, start: u64, end: u64, n_digits: usize, min_qual: &f64, min_cov: &u64) -> io::Result<String> {
    // Add leading zeroes to the start-of-the-chunk index so we can propoerly sort the output files after parallele processing    
    let mut start_string = start.to_string();
    for i in 0..(n_digits - start_string.len()) {
        start_string = "0".to_string() + &start_string;
    }
    // Add leading zeroes to the end-of-the-chunk index so we can propoerly sort the output files after parallele processing
    let mut end_string = end.to_string();
    for i in 0..(n_digits - end_string.len()) {
        end_string = "0".to_string() + &end_string;
    }
    // Output file name for the current chunk
    let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + ".tmp";
    let out = fname_out.clone();
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_out;
    let error_writing_line = "Unable to write line into file: ".to_owned() + &fname_out;
    // println!("{}", fname_out);
    let file_out = File::create(fname_out).expect(&error_writing_file);
    let mut file_out = BufWriter::new(file_out);
    // Input file chunk
    let file = File::open(fname).unwrap();

    let file = File::open(fname).unwrap();
    let mut reader = BufReader::new(file);
    // Navigate to the start of the chunk
    let mut i: u64 = start;
    reader.seek(SeekFrom::Start(start)).unwrap();
    // Read and parse until the end of the chunk
    while i < end {
        // Instantiate the line
        let mut line = String::new();
        // Read the line which automatically movesthe cursor position to the next line
        let _ = reader.read_line(&mut line).unwrap();
        // Find the new cursor position
        i = reader.seek(SeekFrom::Current(0)).unwrap();
        // Remove trailing newline character in Unix-like (\n) and Windows (\r)
        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }
        // println!("i: {} | {:?}", i, line);
        // Parse the sync line
        let vec_line = line.split("\t").collect::<Vec<&str>>();
        let chr = vec_line[0].to_owned();
        let pos = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check format of the file: ".to_owned() + &fname + " as the position is not a valid integer (i.e. u64) at the line whose first 20 characters are: " + &line[0..20] + ".")),
        };
        let ref_allele = match vec_line[2] {
            "A" => "A".to_owned(),
            "T" => "T".to_owned(),
            "C" => "C".to_owned(),
            "G" => "G".to_owned(),
            _ => return Err(Error::new(ErrorKind::Other, "Please check format of the file: ".to_owned() + &fname + " as the reference allele is neither A, T, C, nor G, at the line whose first 20 characters are: " + &line[0..20] + ".")),
        };
     
        // let mut p = parse(&line).expect(&("Input file error, i.e. '".to_owned() + fname + &"' at line with the first 20 characters as: ".to_owned() + &line[0..20] + &".".to_owned()));
    }


    Ok(out)
}


pub fn sync2syncx(fname: &String, n_threads: &u64) -> io::Result<String> {
    let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                     .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                     .join(".");
    let out = bname + "-" + &time.to_string() + ".syncx";

    let chunks = find_file_splits(fname, n_threads);


    
    Ok(out)
}


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
