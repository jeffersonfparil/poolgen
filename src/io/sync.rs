use std;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, SeekFrom, BufReader, BufWriter};
use std::time::{SystemTime, UNIX_EPOCH};
use std::io::{Error, ErrorKind};
use std::sync::{Arc, Mutex};
use nalgebra::{self, DVector, DMatrix};

// Struct for ordering the allele columns by variance in allele frequencies across pools for the syncx format
#[derive(Debug, PartialEq, PartialOrd)]
struct SyncxAlleleFreqs {
    var_x_ave: f64,
    freqs: String,
}

#[derive(Debug, Clone)]
pub struct AlleleCountsOrFrequencies <T, R: nalgebra::Dim, C: nalgebra::Dim> {
    pub coordinate: String,
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: nalgebra::Matrix<T, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<T, R, C>>, // n pools x p alleles
}

pub trait Sync {
    fn counts_to_frequencies(&mut self) -> io::Result<&mut Self>;
    fn filter(&mut self, maf: f64) -> io::Result<&mut Self>;
    fn convert_to_matrix(&self, keep_n_minus_1: bool) -> io::Result<(Vec<String>,
                                                                    Vec<u64>,
                                                                    Vec<String>,
                                                                    nalgebra::Matrix<f64,
                                                                                     nalgebra::Dyn,
                                                                                     nalgebra::Dyn,
                                                                                     nalgebra::VecStorage<f64,
                                                                                                          nalgebra::Dyn,
                                                                                                          nalgebra::Dyn
                                                                                                         >
                                                                                    >
                                                                    )>;
}

impl Sync for Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>> {
    fn counts_to_frequencies(&mut self) -> io::Result<&mut Self> {
        for a in self.into_iter() {
            let n = a.matrix.nrows();
            let p = a.matrix.ncols();
            for i in 0..n {
                let s = a.matrix.row(i).sum();
                for j in 0..p {
                    a.matrix[(i, j)] = a.matrix[(i, j)] / s;
                    if f64::is_nan(a.matrix[(i, j)]) {
                        a.matrix[(i, j)] = 0.0;
                    }
                }
            }
        }
        Ok(self)
    }

    fn filter(&mut self, maf: f64) -> io::Result<&mut Self> {
        let mut i: usize = 0;
        let mut n = self.len();
        while i < n-1 {
            i += 1;
            let a = self[i].clone();
            'per_allele: for j in 0..a.alleles_vector.len() {
                // println!("i={:?}; n={:?}", i, n);
                let _ = match (a.matrix.column(j).mean() < maf) | (a.matrix.column(j).mean() > (1.0-maf)) {
                    true => {self.remove(i); i-=1; n-=1; break 'per_allele},
                    false => 0,
                };
            }
        }
        Ok(self)
    }
    
    fn convert_to_matrix(&self, keep_n_minus_1: bool) -> io::Result<(Vec<String>,Vec<u64>,Vec<String>,nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>>)> {
        // Extract allele counts or freqs matrix
        let mut chromosomes: Vec<String> = Vec::new();
        let mut positions: Vec<u64> = Vec::new();
        let mut alleles: Vec<String> = Vec::new();
        let mut X: nalgebra::Matrix<f64, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<f64, nalgebra::Dyn, nalgebra::Dyn>> = self[0].matrix.clone();
        for i in 1..self.len() {
            let v = &self[i];
            let mut p = v.alleles_vector.len();
            if keep_n_minus_1 & (p>1){
                p -= 1;
            }
            for j in 0..p {
                let a = v.alleles_vector[j].clone();
                chromosomes.push(v.chromosome.to_owned());
                positions.push(v.position);
                alleles.push(a);
                X = X.insert_column(i, 0.0);
                X.set_column(i, &v.matrix.column(j));
            }
        }
        // println!("X: {:?}", X);
        // println!("X0: {:?}", X.column(0));
        // println!("X1: {:?}", X.column(100));
        Ok((chromosomes, positions, alleles, X))
    }
}

fn parse_sync_line (reader: &mut BufReader<File>) -> Option<Vec<String>> {
    // Instatiate the line
    let mut line = String::new();
    // Read the line which automatically movesthe cursor position to the next line
    let _ = reader.read_line(&mut line).unwrap();
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
        return None
    }
    // Parse the sync line
    let vec_line = line.split("\t").collect::<Vec<&str>>().into_iter().map(|x| x.to_owned()).collect::<Vec<String>>();
    Some(vec_line)
}

fn sync2syncx_per_chunk (fname: &String, start: &u64, end: &u64, n_digits: &usize, min_cov: &u64) -> io::Result<String> {
    // Add leading zeroes to the start-of-the-chunk index so we can propoerly sort the output files after parallele processing    
    let mut start_string = start.to_string();
    for _i in 0..(n_digits - start_string.len()) {
        start_string = "0".to_owned() + &start_string;
    }
    // Add leading zeroes to the end-of-the-chunk index so we can propoerly sort the output files after parallele processing
    let mut end_string = end.to_string();
    for _i in 0..(n_digits - end_string.len()) {
        end_string = "0".to_owned() + &end_string;
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
    let mut reader = BufReader::new(file);
    // Navigate to the start of the chunk
    let mut i: u64 = *start;
    reader.seek(SeekFrom::Start(*start)).unwrap();
    // Read and parse until the end of the chunk
    // Instantiate the variables to be reused in each loop
    let mut chr: String;
    let mut pos: u64;
    let mut vec_line: Vec<String>;
    let mut ref_allele: String;
    'lines: while i < *end {
        vec_line = match parse_sync_line(&mut reader) {
            None => continue 'lines,
            Some(x) => x,
        };
        chr = vec_line[0].to_owned();
        pos = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check format of the file: ".to_owned() + &fname + " as the position is not a valid integer (i.e. u64) at the line whose first 20 characters are: " + &vec_line.join("\t")[0..20] + ".")),
        };
        ref_allele = vec_line[2].clone();
        // Read the allele counts and convert them into frequencies excluding the Ns
        let mut coverage: f64;
        let mut a: Vec<f64> = Vec::new();
        let mut t: Vec<f64> = Vec::new();
        let mut c: Vec<f64> = Vec::new();
        let mut g: Vec<f64> = Vec::new();
        let mut d: Vec<f64> = Vec::new();
        for j in 3..vec_line.len() {
            let counts = vec_line[j]
                                                    .split(":")
                                                    .collect::<Vec<&str>>()
                                                    .into_iter()
                                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &vec_line.join("\t")[0..20] + ".")))
                                                    .collect::<Vec<f64>>();
            coverage = counts[0] + counts[1] + counts[2] + counts[3] + counts[5]; // Exclude Ns
            if coverage >= *min_cov as f64 {
                a.push(counts[0] / coverage);
                t.push(counts[1] / coverage);
                c.push(counts[2] / coverage);
                g.push(counts[3] / coverage);
                // Exclude Ns
                d.push(counts[5] / coverage);
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
        // Find the new cursor position
        i = reader.seek(SeekFrom::Current(0)).unwrap();
    }
    Ok(out)
}

pub fn sync2syncx(fname: &String, out: &String, min_cov: &u64, n_threads: &u64) -> io::Result<String> {
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                        .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                        .join(".");
        out = bname + "-" + &time.to_string() + ".syncx";
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
            let vec_out_per_thread = sync2syncx_per_chunk(&fname_clone, &start, &end, &n_digits_clone, &min_cov_clone).unwrap();
            thread_ouputs_clone.lock().unwrap().push(vec_out_per_thread);
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }
    // Extract header lines, i.e. starting with '#' => 35 ascii u8 code
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
    // println!("{:?}", fnames_out);
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

fn load_per_chunk(fname: &String, format: &String, n_pools: &usize, start: &u64, end: &u64) -> io::Result<Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>> {
    let file = File::open(fname).unwrap();
    let mut reader = BufReader::new(file);
    reader.seek(SeekFrom::Start(*start)).unwrap();
    let mut i = *start;
     // Instantiate vector of allele counts across loci
    let mut vec_allele_out: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>> = Vec::new();
    while i < *end {
        let mut line = String::new();
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
        // Split and extract genome coordinate
        let vec_line = line.split("\t").collect::<Vec<&str>>();
        let chr = vec_line[0].to_owned();
        let pos = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check format of the file: ".to_owned() + &fname + " as the position is not a valid integer (i.e. u64) at the line whose first 20 characters are: " + &line[0..20] + ".")),
        };
        // Load counts or frequencies into a vector
        let p = vec_line.len() - (3+0);
        let mut vec_counts_or_freqs: Vec<f64> = Vec::new();
        let mut alleles: Vec<String> = Vec::new();
        if format == &"sync".to_owned() {
            // Extract and merge alleles IDs and counts per pool into a single vector
            let vec_alleles: Vec<String> = vec!["A", "T", "C", "G", "D"].into_iter().map(|x| x.to_owned()).collect::<Vec<String>>(); // Exclude Ns as they are filtered out or ambiguous reads
            for j in 3..(3+p) {
                alleles.push(vec_alleles[j-3].to_owned());
                let counts = vec_line[j]
                                    .split(":")
                                    .collect::<Vec<&str>>()
                                    .into_iter()
                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                    .collect::<Vec<f64>>();
                for k in 0..counts.len() {
                    if k == 4 {
                        continue;
                    }
                    vec_counts_or_freqs.push(counts[k]);
                }
            }
        } else if format == &"syncx".to_owned() {
            // Extract and merge alleles IDs and frequencies per pool into a single vector
            for j in 3..(3+p) {
                // Extract allele
                let vec_line_parse = vec_line[j].split("|").collect::<Vec<&str>>();
                alleles.push(vec_line_parse[0].to_owned());
                let freqs = vec_line_parse[1]
                                    .split(":")
                                    .collect::<Vec<&str>>()
                                    .into_iter()
                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                    .collect::<Vec<f64>>();
                // println!("line: {:?}", freqs);
                for k in 0..freqs.len() {
                    vec_counts_or_freqs.push(freqs[k]);
                }
                // println!("vec_freqs: {:?}", vec_counts_or_freqs);
                // println!("n_pools: {:?}", n_pools);
                // println!("p: {:?}", p);
            }
        } else {
            return Err(Error::new(ErrorKind::Other, "Unrecognised format: ".to_owned() + format + ". Please check the input file: " + fname + " at the line whose first 20 characters are: " + &line[0..20]));
        }
        // Reshape the vector into a matrix of counts where each row is a pool and each column is an allele excluding Ns
        let mut mat_counts_or_freqs = DMatrix::from_row_slice(*n_pools, p, &vec_counts_or_freqs); // generates the matrix row-wise
        if format == &"syncx".to_owned() {
            mat_counts_or_freqs = DMatrix::from_column_slice(*n_pools, p, &vec_counts_or_freqs); // generates the matrix column-wise, i.e. n_pools x 5 alleles
        }
        // Remove non-polymorphic alleles [counts: axn * ones: nx1 = sum: ax1]
        let var = mat_counts_or_freqs.row_variance();
        // println!("#################################");
        // println!("LINE: {:?}", line);
        // println!("p: {:?}", p);
        // println!("MAT: {:?}", mat_counts_or_freqs);
        // println!("VAR: {:?}", var);
        let mut vec_alleles: Vec<String> = Vec::new();
        let mut idx_counter = 0;
        for j in 0..p {
            if var[j] == 0.0 {
                mat_counts_or_freqs = mat_counts_or_freqs.clone().remove_columns(j-idx_counter, 1);
                idx_counter += 1;
            } else {
                vec_alleles.push(alleles[j].clone());
            }
        }
        // println!("Counts term: {:?}", counts);            
        // Put everything together into the allele counts sruct
        let coordinate = chr.to_owned() + "-" + &pos.to_string();
        // println!("COOR: {:?}", coordinate);
        vec_allele_out.push(AlleleCountsOrFrequencies {coordinate: coordinate, chromosome: chr, position: pos, alleles_vector: vec_alleles, matrix: mat_counts_or_freqs});
    }
    // println!("OUT: {:?}", vec_allele_out[vec_allele_out.len()-1]);
    // let x = vec_allele_out[vec_allele_out.len()-1].matrix.clone() * DVector::from_element(3, 1.0);
    // println!("x: {:?}", x);
    Ok(vec_allele_out)
}

pub fn load(fname: &String, n_threads: &u64) -> io::Result<Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>> {

    let chunks = crate::io::find_file_splits(fname, n_threads).unwrap();
    // println!("Chunks: {:?}", chunks);

    // Determine the format of the input file
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
    let mut thread_ouputs: Arc<Mutex<Vec<Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..(*n_threads as usize) {
        // Clone read_chunk parameters
        let fname_clone = fname.clone();
        let format_clone = format.clone();
        let start = chunks[i].clone();
        let end = chunks[i+1].clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let vec_out_per_thread = load_per_chunk(&fname_clone, &format_clone, &n_pools, &start, &end).unwrap();
            thread_ouputs_clone.lock().unwrap().push(vec_out_per_thread);
        });
        thread_objects.push(thread);
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }

    // Trying to sort the output
    let mut vec_acf: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>> = Vec::new();
    for v in thread_ouputs.lock().unwrap().iter() {
        // println!("counter");
        // let q: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>> = t;
        for e in v.into_iter() {
            vec_acf.push(e.to_owned());
        }
    }
    vec_acf.sort_by(|a, b| a.coordinate.partial_cmp(&b.coordinate).unwrap());

    // println!("CONVERT TO MATRIX");
    // let (c, p, a, x) = vec_acf.convert_to_matrix(true).unwrap();
    // println!("X: {:?}", x);

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

    Ok(vec_acf)
}

pub fn loader_with_fun<F>(fname: &String, format: &String, n_pools: &usize, start: &u64, end: &u64, function: F) -> io::Result<Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>> 
where
    F: Fn(&Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>>),
 {
    let file = File::open(fname).unwrap();
    let mut reader = BufReader::new(file);
    reader.seek(SeekFrom::Start(*start)).unwrap();
    let mut i = *start;
     // Instantiate vector of allele counts across loci
    let mut vec_allele_out: Vec<AlleleCountsOrFrequencies<f64, nalgebra::Dyn, nalgebra::Dyn>> = Vec::new();
    while i < *end {
        let mut line = String::new();
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
        // Split and extract genome coordinate
        let vec_line = line.split("\t").collect::<Vec<&str>>();
        let chr = vec_line[0].to_owned();
        let pos = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check format of the file: ".to_owned() + &fname + " as the position is not a valid integer (i.e. u64) at the line whose first 20 characters are: " + &line[0..20] + ".")),
        };
        // Load counts or frequencies into a vector
        let p = vec_line.len() - (3+0);
        let mut vec_counts_or_freqs: Vec<f64> = Vec::new();
        let mut alleles: Vec<String> = Vec::new();
        if format == &"sync".to_owned() {
            // Extract and merge alleles IDs and counts per pool into a single vector
            let vec_alleles: Vec<String> = vec!["A", "T", "C", "G", "D"].into_iter().map(|x| x.to_owned()).collect::<Vec<String>>(); // Exclude Ns as they are filtered out or ambiguous reads
            for j in 3..(3+p) {
                alleles.push(vec_alleles[j-3].to_owned());
                let counts = vec_line[j]
                                    .split(":")
                                    .collect::<Vec<&str>>()
                                    .into_iter()
                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                    .collect::<Vec<f64>>();
                for k in 0..counts.len() {
                    if k == 4 {
                        continue;
                    }
                    vec_counts_or_freqs.push(counts[k]);
                }
            }
        } else if format == &"syncx".to_owned() {
            // Extract and merge alleles IDs and frequencies per pool into a single vector
            for j in 3..(3+p) {
                // Extract allele
                let vec_line_parse = vec_line[j].split("|").collect::<Vec<&str>>();
                alleles.push(vec_line_parse[0].to_owned());
                let freqs = vec_line_parse[1]
                                    .split(":")
                                    .collect::<Vec<&str>>()
                                    .into_iter()
                                    .map(|x| x.to_string().parse::<f64>().expect(&("Please check format of the file: ".to_owned() + &fname + " as the allele counts are not numbers (i.e. f64), at the line whose first 20 characters are: " + &line[0..20] + ".")))
                                    .collect::<Vec<f64>>();
                // println!("line: {:?}", freqs);
                for k in 0..freqs.len() {
                    vec_counts_or_freqs.push(freqs[k]);
                }
                // println!("vec_freqs: {:?}", vec_counts_or_freqs);
                // println!("n_pools: {:?}", n_pools);
                // println!("p: {:?}", p);
            }
        } else {
            return Err(Error::new(ErrorKind::Other, "Unrecognised format: ".to_owned() + format + ". Please check the input file: " + fname + " at the line whose first 20 characters are: " + &line[0..20]));
        }
        // Reshape the vector into a matrix of counts where each row is a pool and each column is an allele excluding Ns
        let mut mat_counts_or_freqs = DMatrix::from_row_slice(*n_pools, p, &vec_counts_or_freqs); // generates the matrix row-wise
        if format == &"syncx".to_owned() {
            mat_counts_or_freqs = DMatrix::from_column_slice(*n_pools, p, &vec_counts_or_freqs); // generates the matrix column-wise, i.e. n_pools x 5 alleles
        }
        // Remove non-polymorphic alleles [counts: axn * ones: nx1 = sum: ax1]
        let var = mat_counts_or_freqs.row_variance();
        // println!("#################################");
        // println!("LINE: {:?}", line);
        // println!("p: {:?}", p);
        // println!("MAT: {:?}", mat_counts_or_freqs);
        // println!("VAR: {:?}", var);
        let mut vec_alleles: Vec<String> = Vec::new();
        let mut idx_counter = 0;
        for j in 0..p {
            if var[j] == 0.0 {
                mat_counts_or_freqs = mat_counts_or_freqs.clone().remove_columns(j-idx_counter, 1);
                idx_counter += 1;
            } else {
                vec_alleles.push(alleles[j].clone());
            }
        }
        // println!("Counts term: {:?}", counts);            
        // Put everything together into the allele counts sruct
        let coordinate = chr.to_owned() + "-" + &pos.to_string();
        // println!("COOR: {:?}", coordinate);
        vec_allele_out.push(AlleleCountsOrFrequencies {coordinate: coordinate, chromosome: chr, position: pos, alleles_vector: vec_alleles, matrix: mat_counts_or_freqs});
        let out = function(&vec_allele_out);
    }
    // println!("OUT: {:?}", vec_allele_out[vec_allele_out.len()-1]);
    // let x = vec_allele_out[vec_allele_out.len()-1].matrix.clone() * DVector::from_element(3, 1.0);
    // println!("x: {:?}", x);
    Ok(vec_allele_out)
}
