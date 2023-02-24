use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader, BufWriter};
use std::io::{Error, ErrorKind};
use std::str;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};
use nalgebra::DVector;

// Struct for a locus from a pileup line
#[derive(Debug)]
pub struct PileupLine {
    chromosome: String,             // chromosome or scaffold name
    position: u64,                  // position in number of bases
    reference_allele: char,         // reference allele
    coverages: Vec<u64>,            // number of times the locus was covered
    read_codes: Vec<Vec<u8>>,       // utf8 read codes corresponding to 'A', 'T', 'C', or 'G' (252 other alleles can be accommodated)
    read_qualities: Vec<Vec<u8>>,   // utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

#[derive(Debug)]
pub struct AlleleCounts {
    chromosome: String, // chromosome or scaffold name
    position: u64,      // position in number of bases
    a: Vec<u64>,        // allele A counts
    t: Vec<u64>,        // allele T counts
    c: Vec<u64>,        // allele C counts
    g: Vec<u64>,        // allele G counts
    n: Vec<u64>,        // allele N or ambiguous allele counts
    d: Vec<u64>,        // allele DEL or deletion counts counts
}

#[derive(Debug)]
pub struct AlleleFrequencies {
    chromosome: String, // chromosome or scaffold name
    position: u64,      // position in number of bases
    a: Vec<f64>,        // allele A counts
    t: Vec<f64>,        // allele T counts
    c: Vec<f64>,        // allele C counts
    g: Vec<f64>,        // allele G counts
    n: Vec<f64>,        // allele N or ambiguous allele counts
    d: Vec<f64>,        // allele DEL or deletion counts counts
}

// Struct for tracking the insertions and deletions prepended with +/i which we want to skip as they refer to the next base position and not the current locus
#[derive(Debug)]
struct IndelMarker {
    indel: bool,    // insertion or deletion, i.e. +/- codes
    count: usize,   // size of the indel, i.e. how many bases
    left: usize,    // how many indel bases are left to be removed (initialised with the maximum possible value of 4294967295)
}

// Struct for ordering the allele columns by variance in allele frequencies across pools for the syncx format
#[derive(Debug, PartialEq, PartialOrd)]
struct SyncxAlleleFreqs {
    var_x_ave: f64,
    freqs: String,
}

// Parse each line
fn parse(line: &String) -> io::Result<Box<PileupLine>> {
    let raw_locus_data: Box<Vec<&str>> = Box::new(line.split("\t").collect());
    // Chromosome or scaffold name
    let chromosome: String = raw_locus_data[0].to_owned();
    // Position or locus coordinate in the genome assembly
    let position = match raw_locus_data[1].parse::<u64>() {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "Please check the format of the input pileup file as position is not a valid integer (i.e. u64).".to_owned())),
    };
    // Allele in the reference genome assembly
    let reference_allele  = match raw_locus_data[2].to_owned().parse::<char>() {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "Please check the format of the input pileup file as the reference allele is not a valid nucleotide base (i.e. not a valid single character).".to_owned())),
    };
    // List of the number of times the locus was read in each pool
    let mut coverages: Vec<u64> = Vec::new();
    for i in (3..raw_locus_data.len()).step_by(3) {
        let cov = match raw_locus_data[i].to_owned().parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check the format of the input pileup file as coverage field/s is/are not valid integer/s (i.e. u64) at pool: ".to_owned() + &(i/3).to_string() + ".")),
        };
        coverages.push(cov);
    }
    // List of alleles that were read in each pool
    let mut read_codes: Vec<Vec<u8>> = Vec::new();
    for i in (4..raw_locus_data.len()).step_by(3) {
        // Parse if the current pool was read at least once for the current locus (i.e. line)
        if coverages[((i-1)/3)-1] > 0 {
            let raw_read_codes = raw_locus_data[i].as_bytes().to_owned();
            let mut alleles: Vec<u8> = Vec::new();
            let mut indel_marker: IndelMarker = IndelMarker{indel: false, count: 0, left: 4294967295}; // Instantiate the indel marker with no indel (i.e. indel==false and count==0) and the maximum number of indels left (i.e. left==4294967295 which is the maximum value for usize type of left)
            'per_pool: for j in 0..raw_read_codes.len() {
                let code = raw_read_codes[j];
                // Are we dealing with indels?
                if indel_marker.indel {
                    // Find the firs digit of the number of indel/s
                    if (indel_marker.count == 0) & (indel_marker.left == 4294967295) {
                        indel_marker.count = match str::from_utf8(&[code]).unwrap().parse::<usize>() {
                            Ok(x) => x,
                            Err(_) => return Err(Error::new(ErrorKind::Other, "Check the codes for insertions and deletion, i.e. they must be integers after '+' and '-' at pool: ".to_owned() + &((i-1)/3).to_string() + ".")),
                        };
                        continue 'per_pool;
                    }
                    // Find the next digit of the number of indels, if we have more than 9 indels
                    if (indel_marker.count > 0) & (indel_marker.left == 4294967295) {
                        let x = str::from_utf8(&[code]).unwrap().parse::<i64>();
                        let y = match x {
                            Ok(z) => z,
                            Err(_) => -1,
                        };
                        if y >= 0 {
                            let z = indel_marker.count.to_string() + "0";
                            indel_marker.count = z.parse::<usize>().unwrap() + (y as usize);
                            // println!("{}", z.parse::<usize>().unwrap() + y);
                            continue 'per_pool;
                        } else {
                            indel_marker.left = indel_marker.count - 1;
                            continue 'per_pool;
                        }
                        // println!("{:?}", y);
                    }
                    // Remove indels
                    if indel_marker.left > 0 {
                        indel_marker.left -= 1;
                        continue 'per_pool;
                    } else {
                        // Reset after removing all indels and proceed to testing
                        indel_marker.indel = false;
                        indel_marker.count = 0;
                        indel_marker.left = 4294967295;
                    }
                }
                // Is the current code another indel?
                // Remove insertion and deletion codes (unicodes: '+' == 43 and '-' == 45)
                if (code == 43) | (code == 45) {
                    indel_marker.indel = true;
                    indel_marker.count = 0;
                    continue 'per_pool;
                }
                // Is the current code a read start/end marker?
                // Remove read start and end codes (unicodes: '^' == 94 and '$' == 36)
                if (code == 94) | (code == 36){
                    if code == 94 {
                        // If we have a start marker then we use the IndelMarker struct to get rig of the next character which encodes for the mapping quality of the read
                        indel_marker.indel = true;
                        indel_marker.count = 1;
                        indel_marker.left = 1;
                    }
                    continue 'per_pool;
                }
                // Reference allele (unicodes: ',' == 44 and '.' == 46) and 
                // non-reference alleles (unicodes: 'A'/'a' == 65/97,
                //                                  'T'/'t' == 84/116,
                //                                  'C'/'c' == 67/99,
                //                                  'G'/'g' == 71/103, and
                //                                  '*' == 42 codes for deletion on the current locus which is different from the \-[0-9]+[ACGTNacgtn]+ which encodes for indels in the next position/s)
if (code == 44) | (code == 46) {
  reference_allele.to_string().as_bytes()[0]; // whatever the reference allele is
}
let a: u8 = match code {
                    65 => 65,   // A
                    97 => 65,   // a -> A
                    84 => 84,   // T
                    116 => 84,  // t -> T
                    67 => 67,   // C
                    99 => 67,   // c -> C
                    71 => 71,   // G
                    103 => 71,  // g -> G
                    42 => 68,   // * -> D
                    _ => 78,    // N
                };
                alleles.push(a);
            }
            read_codes.push(alleles);
        } else {
            read_codes.push(Vec::new());
        }
    }
    // println!("{:?}", read_codes);
    // List of the qualities of the read alleles in each pool
    let mut read_qualities: Vec<Vec<u8>> = Vec::new();
    for i in (5..raw_locus_data.len()).step_by(3) {
        if coverages[((i-1)/3)-1] > 0 {
                let qualities = raw_locus_data[i].as_bytes().to_owned();
                read_qualities.push(qualities);
            } else {
                read_qualities.push(Vec::new());
            }
    }
    // println!("{:?}", read_qualities);
    // Output PileupLine struct
    let out = Box::new(PileupLine {
        chromosome: chromosome,
        position: position,
        reference_allele: reference_allele,
        coverages: coverages,
        read_codes: read_codes,
        read_qualities: read_qualities,
    });
    // Sanity check to see if the coverage, number of alleles and quality codes match per pool
    let mut c: u64;
    let mut a: u64;
    let mut q: u64;
    for i in 0..out.coverages.len() {
        c = out.coverages[i];
        a = out.read_codes[i].len() as u64;
        q = out.read_qualities[i].len() as u64;
        if (c != a) | (c != q) | (a != q) {
            return Err(Error::new(ErrorKind::Other, "Please check the format of the input pileup file as the coverages, number of read alleles and read qualities do not match at pool: ".to_owned() + &(i+1).to_string() + "."));
        }
    }
    return Ok(out);
}

impl PileupLine {
    fn mean_quality(&self) -> io::Result<f64> {
        let mut s: f64 = 0.0;
        let mut n: f64 = 0.0;
        for q in &self.read_qualities {
            for x in q.iter().map(|&x| f64::from(x)).collect::<Vec<f64>>().iter() {
                // println!("{:?}", x);
                if *x < 33.0 {
                    return Err(Error::new(ErrorKind::Other, "Phred score out of bounds."));
                } else {
                    s += x - 33.0; // Assumes PHRED 33 (i.e. !=10^(-0/10) to I=10^(-40/10))
                    n += 1.0;
                }
            }
        }
        let out = f64::powf(10.0, -s/(10.0*n));
        // println!("{:?}", out);
        Ok(out)
    }

    fn filter(&mut self, min_coverage: &u64, min_quality: &f64, remove_ns: &bool) -> io::Result<&mut Self> {
        // Convert low quality bases into Ns
        let n = &self.read_qualities.len();
        for i in 0..*n {
            let pool = &self.read_qualities[i];
            let mut j: i64 = 0;
            while j < self.read_codes[i].len() as i64 {
                if pool[j as usize] < 33 {
                    return Err(Error::new(ErrorKind::Other, "Phred score out of bounds."));
                } else {
                    let q = f64::powf(10.0, -(pool[j as usize] as f64 - 33.0) / 10.0); 
                    if q > *min_quality {
                        if *remove_ns {
                            self.read_codes[i].remove(j as usize);
                            j -= 1;
                        } else {
                            self.read_codes[i][j as usize] = 78; // convert to N
                        }
                        self.coverages[i] = self.coverages[i] - 1; // remove the coverage for ambiguous alleles
                    }
                }
                j += 1;
            }
        }
        // All the pools needs be have been covered at least min_coverage times
        for c in &self.coverages {
            // println!("{:?}", c);
            if c < min_coverage {
                return Err(Error::new(ErrorKind::Other, "Filtered out."));
            }
        }
        Ok(self)
    }

    fn reads_to_counts(&self) -> io::Result<Box<AlleleCounts>> {
        let mut out = Box::new(AlleleCounts {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            a: Vec::new(),
            t: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            n: Vec::new(),
            d: Vec::new()});
        for pool in &self.read_codes {
            let mut a_per_pool: u64 = 0;
            let mut t_per_pool: u64 = 0;
            let mut c_per_pool: u64 = 0;
            let mut g_per_pool: u64 = 0;
            let mut n_per_pool: u64 = 0;
            let mut d_per_pool: u64 = 0;    
            for allele in pool {
                match allele {
                    65 => a_per_pool += 1,
                    84 => t_per_pool += 1,
                    67 => c_per_pool += 1,
                    71 => g_per_pool += 1,
                    68 => d_per_pool += 1,
                    _  => n_per_pool += 1,
                };
            }
            out.a.push(a_per_pool);
            out.t.push(t_per_pool);
            out.c.push(c_per_pool);
            out.g.push(g_per_pool);
            out.n.push(n_per_pool);
            out.d.push(d_per_pool);
        }
        Ok(out)
    }
    
    fn reads_to_frequencies(&self, remove_n: bool) -> io::Result<Box<AlleleFrequencies>> {
        let mut out = Box::new(AlleleFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            a: Vec::new(),
            t: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            n: Vec::new(),
            d: Vec::new()});
        let counts = self.reads_to_counts().unwrap();
        let n = counts.a.len();
        for i in 0..n {
            let sum = match remove_n {
                true => (counts.a[i] + counts.t[i] + counts.c[i] + counts.g[i] + counts.d[i]) as f64,
                false => (counts.a[i] + counts.t[i] + counts.c[i] + counts.g[i] + counts.n[i] + counts.d[i]) as f64,
            };
            out.a.push((counts.a[i] as f64) / sum);
            out.t.push((counts.t[i] as f64) / sum);
            out.c.push((counts.c[i] as f64) / sum); 
            out.g.push((counts.g[i] as f64) / sum);
            if remove_n {
                out.n.push((0.000000000 as f64) / sum);
            } else {
                out.n.push((counts.n[i] as f64) / sum);
            } 
            out.d.push((counts.d[i] as f64) / sum); 
        }
        Ok(out)
    }

}

fn pileup2sync_chunk(fname: &String, start: &u64, end: &u64, n_digits: &usize, min_qual: &f64, min_cov: &u64, remove_ns: &bool, fmt: &String) -> io::Result<String>{
    // Add leading zeros in front of the start file position so that we can sort the output files per chuck or thread properly
    let mut start_string = start.to_string();
    for  i in 0..(n_digits - start_string.len()) {
        start_string = "0".to_owned() + &start_string;
    }
    // Add leading zeros in front of the end file position so that we can sort the output files per chuck or thread properly
    let mut end_string = end.to_string();
    for  i in 0..(n_digits - end_string.len()) {
        end_string = "0".to_owned() + &end_string;
    }
    // Output temp file for the chunk    
    let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + "." + fmt + ".tmp";
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
    while i < *end {
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
        // Parse the pileup line
        let mut p = parse(&line).expect(&("Input file error, i.e. '".to_owned() + fname + "' at line with the first 20 characters as: " + &line[0..20] + "."));
        // println!("i: {} | Raw line: {:?}", i, p);
        // Filter
        match p.filter(min_cov, min_qual, remove_ns) {
            Ok(x) => x,
            Err(_) => continue,
        };
        // println!("i: {} | Filtered: {:?}", i, p);
        // println!("i: {} | Counts: {:?}", i, r);
        // println!("i: {} | Frequencies: {:?}", i, f);
        // println!("#########################################");
        // Instantiate the output line
        let mut x = vec![p.chromosome.clone(), p.position.to_string(), p.reference_allele.to_string()];
        if fmt == "sync" {
            // Convert to a vector of counts (Note: outputs an "empty" struct if the locus has been filtered out by minimum coverage)
            let r = p.reads_to_counts().unwrap();
            // Sync canonical popoolation2 format
            for i in 0..r.a.len() {
                let column = vec![r.a[i].to_string(),
                                                r.t[i].to_string(),
                                                r.c[i].to_string(),
                                                r.g[i].to_string(),
                                                r.n[i].to_string(),
                                                r.d[i].to_string()];
                x.push(column.join(":"));
            }
        } else if fmt == "syncx" {
            // Convert to a vector of frequencies (Note: outputs an "empty" struct if the locus has been filtered out by minimum coverage; also do not count Ns hence allele frequencies sum up to one across alleles A, T, C, G, and DEL)
            let f = p.reads_to_frequencies(true).unwrap();
            // Include only the polymorphic alleles
            let frequencies = vec![f.a, f.t, f.c, f.g, f.n, f.d];
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
            // if x_tmp.len() > 2 {
            //     println!("{:?}", x_tmp);
            // }
            // Append the polymorphic and sorted allele frequencies
            for xi in x_tmp.iter() {
                x.push(xi.freqs.to_owned());
            }
        } else {
            return Err(Error::new(ErrorKind::Other, "Format: ".to_owned() + &fmt + " not reconised. Please use: 'sync', or 'syncx'."));
        }
        // Write the line
        let data = x.join("\t") + "\n";
        if x.len() > 3 {
            file_out.write_all(data.as_bytes()).expect(&error_writing_line);
        } else {
            continue;
        }
    }
    Ok(out)
}

// Convert pileup into sync or syncx
pub fn pileup2sync(fname: &String, out: &String, pool_names: &String, min_qual: &f64, min_cov: &u64, remove_ns: &bool, file_format: &String, n_threads: &u64) -> io::Result<String> {
    // Output filename
    let mut out = out.to_owned();
    if out == "".to_owned() {
        let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
        let bname = fname.split(".").collect::<Vec<&str>>().into_iter().map(|a| a.to_owned())
                                    .collect::<Vec<String>>().into_iter().rev().collect::<Vec<String>>()[1..].to_owned()
                                    .into_iter().rev().collect::<Vec<String>>().join(".");
        out = bname.to_owned() + "-" + &time.to_string() + "." + &(file_format.to_owned());
    }
    // Pool names
    let mut names: Vec<String> = Vec::new();
    let file_names = match File::open(pool_names) {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "Pool names file not found which is required for converting pileup into sync/syncx. Please make sure to include a valid file for the --pool_names parameter.")),
    };
    let reader_names = BufReader::new(file_names);
    for line in reader_names.lines() {
        names.push(line.unwrap().to_owned());
    }
    let names = names.join("\t");

    // // Find the positions whereto split the file into n_threads pieces
    let chunks = crate::io::find_file_splits(fname, n_threads).unwrap();
    let n_digits = chunks[*n_threads as usize].to_string().len();
    println!("Chunks: {:?}", chunks);
    // Tuple arguments of pileup2sync_chunks
    // Instantiate thread object for parallel execution
    let mut thread_objects = Vec::new();
    // Vector holding all returns from pileup2sync_chunk()
    let mut thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    // Making four separate threads calling the `search_for_word` function
    for i in 0..(*n_threads as usize) {
        // Clone pileup2sync_chunk parameters
        let fname_clone = fname.clone();
        let start = chunks[i].clone();
        let end = chunks[i+1].clone();
        let n_digits_clone = n_digits.clone();
        let min_qual_clone = min_qual.clone();
        let min_cov_clone = min_cov.clone();
        let remove_ns = remove_ns.clone();
        let file_format_clone = file_format.clone();
        let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
        let thread = std::thread::spawn(move || {
            let fname_out_per_thread = pileup2sync_chunk(&fname_clone, &start, &end, &n_digits_clone, &min_qual_clone, &min_cov_clone, &remove_ns, &file_format_clone).unwrap();
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
    let _ = match &file_format[..] {
        "sync" => file_out.write_all(("#chr\tpos\tref\t".to_owned() + &names + "\n").as_bytes()),
        "syncx" => file_out.write_all(("#chr\tpos\tref\t".to_owned() + &names + "\n").as_bytes()),
        _ => return Err(Error::new(ErrorKind::Other, "Format: ".to_owned() + &file_format + " not reconised. Please use: 'sync', or 'syncx'.")),
    };
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
