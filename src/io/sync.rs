use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader, BufWriter};
use std::time::{SystemTime, UNIX_EPOCH};
use std::io::{Error, ErrorKind};
use std::sync::{Arc, Mutex};
use nalgebra::{DVector, DMatrix};

// Struct for ordering the allele columns by variance in allele frequencies across pools for the syncx format
#[derive(Debug, PartialEq, PartialOrd)]
struct SyncxAlleleFreqs {
    var_x_ave: f64,
    freqs: String,
}

fn sync2syncx_per_chunk (fname: &String, start: &u64, end: &u64, n_digits: &usize, min_cov: &u64) -> io::Result<String> {
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
    let mut i: u64 = *start;
    reader.seek(SeekFrom::Start(*start)).unwrap();
    // Read and parse until the end of the chunk
    'lines: while i < *end {
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
        // Ignore commented-out lines (i.e. '#' => 35)
        if line.as_bytes()[0] == 35 as u8 {
            continue
        }
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
        // Read the allele counts and convert them into frequencies excluding the Ns
        let n_pools = vec_line.len();
        let mut a: Vec<f64> = Vec::new();
        let mut t: Vec<f64> = Vec::new();
        let mut c: Vec<f64> = Vec::new();
        let mut g: Vec<f64> = Vec::new();
        let mut d: Vec<f64> = Vec::new();
        for i in 3..n_pools {
            let counts = vec_line[i]
                                                    .split(":")
                                                    .collect::<Vec<&str>>()
                                                    .into_iter()
                                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                                    .collect::<Vec<f64>>();
            let sum = counts[0] + counts[1] + counts[2] + counts[3] + counts[5]; // Exclude Ns
            if sum >= *min_cov as f64 {
                a.push(counts[0] / sum);
                t.push(counts[1] / sum);
                c.push(counts[2] / sum);
                g.push(counts[3] / sum);
                d.push(counts[5] / sum);
            } else {
                continue 'lines;
            }
        }
        // println!("{:?}", a.len());
        // Include only the polymorphic alleles
        let frequencies = vec![a, t, c, g, d];
        // Temporarily store the allele frequencies and note of their variances so we can sort them later by the variance x mean
        let mut x_tmp: Vec<SyncxAlleleFreqs> = Vec::new();
        for i in 0..frequencies.len() {
            let freqs = frequencies[i].clone();
            let var = DVector::from_column_slice(&freqs).variance();
            let ave = DVector::from_column_slice(&freqs).mean();
            // println!("{:?}", var);
            if var > 0.0 {
                let a = match i {
                    0 => "a".to_owned(),
                    1 => "t".to_owned(),
                    2 => "c".to_owned(),
                    3 => "g".to_owned(),
                    4 => "n".to_owned(),
                    _ => "d".to_owned(),
                };
                let mut column = vec![a];
                column.push(freqs.iter().map(|y| y.to_string()).collect::<Vec<String>>().join(":"));
                x_tmp.push(SyncxAlleleFreqs{var_x_ave: var*ave, freqs: column.join("|")});
            } else {
                continue;
            }
        }
        // Sort the allele columns by decreasing variance such that the most polymorphic allele across pools is in the first column
        x_tmp.sort_by(|x, y| y.var_x_ave.partial_cmp(&x.var_x_ave).unwrap());
        // Instantiate the output syncx line
        let mut x = vec![chr, pos.to_string(), ref_allele];
        // Append the polymorphic and sorted allele frequencies
        for xi in x_tmp.iter() {
            x.push(xi.freqs.to_owned());
        }
        // Write the syncx line
        let data = x.join("\t") + "\n";
        if x.len() > 3 {
            file_out.write_all(data.as_bytes()).expect(&error_writing_line);
        } else {
            continue;
        }
    }
    Ok(out)
}

pub fn sync2syncx(fname: &String, min_cov: &u64, n_threads: &u64) -> io::Result<String> {
    let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                     .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                     .join(".");
    let out = bname + "-" + &time.to_string() + ".syncx";

    let chunks = crate::io::find_file_splits(fname, n_threads);
    let n_digits = chunks[*n_threads as usize].to_string().len();
    println!("Chunks: {:?}", chunks);

    // Instantiate thread object for parallel execution
    let mut thread_objects = Vec::new();
    // Vector holding all returns from read_chunk()
    let mut thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..(*n_threads as usize) {
        // Clone read_chunk parameters
        let fname_clone = fname.clone();
        let start = chunks[i].clone();
        let end = chunks[i+1].clone();
        let n_digits_clone = n_digits.clone();
        let min_cov_clone = min_cov.clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let fname_out_per_thread = sync2syncx_per_chunk(&fname_clone, &start, &end, &n_digits_clone, &min_cov_clone).unwrap();
            thread_ouputs_clone.lock().unwrap().push(fname_out_per_thread);
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }
    // Instatiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &out;
    let mut file_out = File::create(&out).expect(&error_writing_file);
    // But first, extract header lines, i.e. starting with '#' => 35 ascii u8 code
    let file_in = File::open(fname).unwrap();
    let mut file_in = BufReader::new(file_in);
    let mut header_end = false;
    while header_end == false {
        // Instantiate the line
        let mut line = String::new();
        // Read the line which automatically movesthe cursor position to the next line
        let _ = file_in.read_line(&mut line).unwrap();
        if line.as_bytes()[0] == 35 as u8 {
            file_out.write_all(line.to_owned().as_bytes()).unwrap();
        } else {
            header_end = true;
        }
    }

    // Extract output filenames from each thread into a vector and sort them
    let mut fnames_out: Vec<String> = Vec::new();
    for f in thread_ouputs.lock().unwrap().iter() {
        fnames_out.push(f.to_owned());
    }
    fnames_out.sort();
    println!("{:?}", fnames_out);
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
