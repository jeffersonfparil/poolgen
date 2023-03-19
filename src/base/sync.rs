use std::io::{self, prelude::*, Error, ErrorKind, BufReader, BufWriter, SeekFrom};
use std::fs::{File, OpenOptions};
use std::str;
use nalgebra::DMatrix;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};
use crate::base::*;

impl Parse<LocusCounts> for String {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<LocusCounts>> {
        // Remove trailing newline character in Unix-like (\n) and Windows (\r)
        let mut line = self.clone();
        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }
        // Ignore commented-out lines (i.e. '#' => 35)
        if line.as_bytes()[0] == 35 as u8 {
            return Err(Error::new(ErrorKind::Other, "Commented out line"))
        }
        // Parse the sync line
        let vec_line = line.split("\t")
                                        .collect::<Vec<&str>>()
                                        .into_iter()
                                        .map(|x| x.to_owned())
                                        .collect::<Vec<String>>();
        let n: usize = vec_line.len() - 3;
        let p: usize = 6;
        let chromosome = vec_line[0].to_owned();
        let position = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check format of the file: position is not and integer.")),
        };
        let alleles_vector = vec!["A", "T", "C", "G", "N", "D"].into_iter()
                                                                       .map(|x| x.to_owned())
                                                                       .collect::<Vec<String>>();
        // Read the allele counts into the counts matrix
        let mut matrix: DMatrix<u64> = DMatrix::from_element(n, p, 0);
        let mut counts: Vec<u64>;
        for i in 0..n {
            counts = vec_line[i+3].split(":")
                                .map(|x| x.to_string().parse::<u64>().expect("Please check the input sync file as the allele counts are not valid integers."))
                                .collect::<Vec<u64>>();
            for j in 0..p {
                matrix[(i, j)] = counts[j];
            }
        }
        return Ok(Box::new(LocusCounts { chromosome: chromosome,
                                         position: position,
                                         alleles_vector: alleles_vector,
                                         matrix: matrix }))
    }   
}

impl Filter for LocusCounts {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let out = self.clone();
        Ok(Box::new(out))
    }
    
    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        let (n, p) = self.matrix.shape();
        let row_sums = self.matrix.column_sum(); // summation across the columns which means sum of all elements per row
        let mut matrix: DMatrix<f64> = DMatrix::from_element(n, p, 0.0 as f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = self.matrix[(i, j)] as f64 / row_sums[i] as f64;
            }
        }
        Ok(Box::new(LocusFrequencies { chromosome: self.chromosome.clone(),
                                       position: self.position.clone(),
                                       alleles_vector: self.alleles_vector.clone(),
                                       matrix: matrix }))
    }
    
    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        // Remove Ns
        if filter_stats.remove_ns {
            let matrix: DMatrix<u64>;
            let i = match self.alleles_vector.iter().position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned())) {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                matrix = self.matrix.clone().remove_column(i as usize);
                self.matrix = matrix;
            }
        }
        // Filter by minimum coverage
        let mean_coverage = self.matrix.column_sum(); // summation across the columns which means sum of all elements per row
        if mean_coverage.min() < filter_stats.min_coverage {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter by minimum allele frequency
        let allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Cannot convert locus counts to locus frequencies.")),
        };
        let mean_allele_freqs = allele_frequencies.matrix.row_mean(); // average across the rows which means sum of all elements per column
        if (1.00 - mean_allele_freqs.max()) < filter_stats.min_allele_frequency {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Remove zero columns
        let mut p: usize = self.matrix.ncols();
        let mut i: usize = 0;
        let mut j: usize = 0;
        let mut matrix: DMatrix<u64> = self.matrix.clone();
        while j < p {
            if mean_allele_freqs[i] == 0.0 {
                self.alleles_vector.remove(j);
                matrix = matrix.remove_column(j);
                p -= 1
            } else {
                j += 1
            }
            i += 1;
        }
        self.matrix = matrix;
        Ok(self)
    }
}

impl Filter for LocusFrequencies {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let (n, p) = self.matrix.shape();
        let mut matrix: DMatrix<u64> = DMatrix::from_element(n, p, 0);
        let mut max_n: f64;
        for i in 0..n {
            max_n = 1.00 / self.matrix.row(i).min();
            for j in 0..p {
                matrix[(i, j)] = (max_n * self.matrix[(i, j)]).round() as u64;
            }
        }
        Ok(Box::new(LocusCounts { chromosome: self.chromosome.clone(),
                                  position: self.position.clone(),
                                  alleles_vector: self.alleles_vector.clone(),
                                  matrix: matrix }))
    }
    
    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        // Recompute the frequencies using frequencies when the number of colulmns or one or more alleles have been filtered out/removed
        let (n, p) = self.matrix.shape();
        let row_sums = self.matrix.column_sum(); // summation across the columns which means sum of all elements per row
        let mut matrix: DMatrix<f64> = DMatrix::from_element(n, p, 0.0 as f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = self.matrix[(i, j)] / row_sums[i];
            }
        }
        Ok(Box::new(LocusFrequencies { chromosome: self.chromosome.clone(),
                                       position: self.position.clone(),
                                       alleles_vector: self.alleles_vector.clone(),
                                       matrix: matrix }))
    }
    
    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        // Also, cannot filter by minimum coverage as that data is lost from counts to frequencies conversion
        // Remove Ns
        if filter_stats.remove_ns {
            let matrix: DMatrix<f64>;
            let i = match self.alleles_vector.iter().position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned())) {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                matrix = self.matrix.clone().remove_column(i as usize);
                self.matrix = matrix;
            }
        }
        // Recompute frequencies after removing Ns
        let recomputed_self = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "T_T Cannot convert locus counts to locus frequencies.")),
        };
        self.alleles_vector = recomputed_self.alleles_vector;
        self.matrix = recomputed_self.matrix;
        // Filter by minimum allele frequency
        let mean_allele_freqs = self.matrix.row_mean(); // average across the rows which means sum of all elements per column
        if (1.00 - mean_allele_freqs.max()) < filter_stats.min_allele_frequency {
            return Err(Error::new(ErrorKind::Other, "-_- Filtered out."));
        }
        // Remove zero columns
        let mut p: usize = self.matrix.ncols();
        let mut i: usize = 0;
        let mut j: usize = 0;
        let mut matrix: DMatrix<f64> = self.matrix.clone();
        while j < p {
            if mean_allele_freqs[i] == 0.0 {
                self.alleles_vector.remove(j);
                matrix = matrix.remove_column(j);
                p -= 1
            } else {
                j += 1
            }
            i += 1;
        }
        self.matrix = matrix;
        // Recompute frequencies after removing zero columns
        let recomputed_self = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "T_T Cannot convert locus counts to locus frequencies.")),
        };
        self.alleles_vector = recomputed_self.alleles_vector;
        self.matrix = recomputed_self.matrix;
        Ok(self)
    }
}

impl ChunkyReadAnalyseWrite<LocusCounts, fn(&mut LocusCounts, &FilterStats) -> Option<String>> for FileSync {
    fn per_chunk(&self, start: &u64, end: &u64, outname_ndigits: &usize, filter_stats: &FilterStats, function: fn(&mut LocusCounts, &FilterStats) -> Option<String>) -> io::Result<String> {
        let fname = self.filename.clone();
        // Add leading zeros in front of the start file position so that we can sort the output files per chuck or thread properly
        let mut start_string = start.to_string();
        for _ in 0..(outname_ndigits - start_string.len()) {
            start_string = "0".to_owned() + &start_string;
        }
        // Add leading zeros in front of the end file position so that we can sort the output files per chuck or thread properly
        let mut end_string = end.to_string();
        for _ in 0..(outname_ndigits - end_string.len()) {
            end_string = "0".to_owned() + &end_string;
        }
        // Output temp file for the chunk    
        let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + ".fisher_test.tmp";
        let out = fname_out.clone();
        let error_writing_file = "T_T Unable to create file: ".to_owned() + &fname_out;
        let error_writing_line = "T_T Unable to write line into file: ".to_owned() + &fname_out;
        let file_out = File::create(fname_out).expect(&error_writing_file);
        let mut file_out = BufWriter::new(file_out);
        // Input file chunk
        let file = File::open(fname.clone()).unwrap();
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
            // Parse the pileup line
            let mut locus_counts: Box<LocusCounts> = match line.lparse() {
                Ok(x) => x,
                Err(x) => match x.kind() {
                    ErrorKind::Other => continue,
                    _ => return Err(Error::new(ErrorKind::Other, "T_T Input sync file error, i.e. '".to_owned() + &fname + "' at line with the first 20 characters as: " + &line[0..20] + "."))
                }
            };
            // Write the line
            let _ = match function(&mut locus_counts, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }
    
    fn read_analyse_write(&self, filter_stats: &FilterStats, out: &String, n_threads: &u64, function: fn(&mut LocusCounts, &FilterStats) -> Option<String>) -> io::Result<String> {
        // Unpack pileup and pool names filenames
        let fname = self.filename.clone();
        let test = self.test.clone();
        // Output filename
        let mut out = out.to_owned();
        if out == "".to_owned() {
            let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
            let bname = fname.split(".").collect::<Vec<&str>>().into_iter().map(|a| a.to_owned())
                                        .collect::<Vec<String>>().into_iter().rev().collect::<Vec<String>>()[1..].to_owned()
                                        .into_iter().rev().collect::<Vec<String>>().join(".");
            out = bname.to_owned() + "-" + &time.to_string() + "-" + &test + ".csv";
        }
        // Instatiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        // let mut file_out = File::create(&out).expect(&error_writing_file);
        let mut file_out = OpenOptions::new().create_new(true)
                                                .write(true)
                                                .append(false)
                                                .open(&out)
                                                .expect(&error_writing_file);
        // // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let outname_ndigits = chunks[*n_threads as usize].to_string().len();
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let mut thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
        // Making four separate threads calling the `search_for_word` function
        for i in 0..(*n_threads as usize) {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i].clone();
            let end = chunks[i+1].clone();
            let outname_ndigits = outname_ndigits.clone();
            let filter_stats = filter_stats.clone();
            let mut thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let fname_out_per_thread = self_clone.per_chunk(&start, &end, &outname_ndigits, &filter_stats, function).unwrap();
                thread_ouputs_clone.lock().unwrap().push(fname_out_per_thread);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            let _ = thread.join().expect("Unknown thread error occured.");
        }
        // Write out
        file_out.write_all(("#chr,pos,alleles,pvalue\n").as_bytes()).unwrap();
        // Extract output filenames from each thread into a vector and sort them
        let mut fnames_out: Vec<String> = Vec::new();
        for f in thread_ouputs.lock().unwrap().iter() {
            fnames_out.push(f.to_owned());
        }
        fnames_out.sort();
        // Iterate across output files from each thread, and concatenate non-empty files
        for f in fnames_out {
            let mut file: File = File::open(&f).unwrap();
            if file.metadata().unwrap().len() > 0 {
                io::copy(&mut file, &mut file_out).unwrap();
            }
            // Clean-up: remove temporary output files from each thread
            let error_deleting_file = "Unable to remove file: ".to_owned() + &f;
            std::fs::remove_file(f).expect(&error_deleting_file);
        }
        Ok(out)
    }
}
