use std::io::{self, prelude::*, Error, ErrorKind, BufReader, BufWriter, SeekFrom};
use std::fs::{File, OpenOptions};
use std::str;
use nalgebra::DMatrix;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};
use crate::base::*;

impl Parse<PileupLine> for String {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<PileupLine>> {
        let raw_locus_data: Box<Vec<&str>> = Box::new(self.split("\t").collect());
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
                                continue 'per_pool;
                            } else {
                                indel_marker.left = indel_marker.count - 1;
                                continue 'per_pool;
                            }
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
                    let a: u8 = if (code == 44) | (code == 46) {
                        // whatever the reference allele is
                        reference_allele.to_string().as_bytes()[0]
                    } else {
                        match code {
                            65  => 65,  // A
                            97  => 65,  // a -> A
                            84  => 84,  // T
                            116 => 84,  // t -> T
                            67  => 67,  // C
                            99  => 67,  // c -> C
                            71  => 71,  // G
                            103 => 71,  // g -> G
                            42  => 68,  // * -> D
                            _   => 78,  // N
                        }
                    };
                    alleles.push(a);
                }
                read_codes.push(alleles);
            } else {
                read_codes.push(Vec::new());
            }
        }
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
        return Ok(out)
    }   
}

impl Filter for PileupLine {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let n: usize = self.coverages.len();
        let p: usize = 6;
        let mut matrix: DMatrix<u64> = DMatrix::from_element(n, p, 0 as u64);
        let mut counts: Vec<Vec<u64>> = vec![Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()];
        let alleles_vector = vec!["A", "T", "C", "G", "D", "N"].into_iter().map(|x| x.to_owned()).collect::<Vec<String>>();
        for pool in &self.read_codes {
            let mut counts_per_pool: Vec<u64> = vec![0, 0, 0, 0, 0, 0];
            for allele in pool {
                match allele {
                    65 => counts_per_pool[0] += 1, // A
                    84 => counts_per_pool[1] += 1, // T
                    67 => counts_per_pool[2] += 1, // C
                    71 => counts_per_pool[3] += 1, // G
                    68 => counts_per_pool[4] += 1, // D
                    _  => counts_per_pool[5] += 1, // N
                };
            }
            for j in 0..p {
                counts[j].push(counts_per_pool[j]);
            }
        }
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = counts[j][i];
            }
        }
        Ok(Box::new(LocusCounts{ chromosome: self.chromosome.clone(),
                                 position: self.position.clone(),
                                 alleles_vector: alleles_vector,
                                 matrix: matrix }))
    }
    
    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        let locus_counts = self.to_counts().unwrap();
        let (n, p) = locus_counts.matrix.shape();
        let row_sums = locus_counts.matrix.column_sum(); // summation across the columns which means sum of all elements per row
        let mut matrix: DMatrix<f64> = DMatrix::from_element(n, p, 0.0 as f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = locus_counts.matrix[(i, j)] as f64 / row_sums[i] as f64;
            }
        }
        Ok(Box::new(LocusFrequencies { chromosome: locus_counts.chromosome,
                                       position: locus_counts.position,
                                       alleles_vector: locus_counts.alleles_vector,
                                       matrix: matrix }))
    }
    
    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Convert low quality bases into Ns
        let n = self.read_qualities.len();
        for i in 0..n {
            let mut j: usize = 0;
            while j < self.read_codes[i].len() {
                if self.read_qualities[i][j] < 33 {
                    return Err(Error::new(ErrorKind::Other, "Phred score out of bounds."));
                } else {
                    let q = f64::powf(10.0, -(self.read_qualities[i][j] as f64 - 33.0) / 10.0); 
                    if q > filter_stats.min_quality {
                        self.read_codes[i][j] = 78; // convert to N
                    }
                    if filter_stats.remove_ns & (self.read_codes[i][j] == 78) {
                        self.read_codes[i].remove(j);
                        self.read_qualities[i].remove(j);
                        self.coverages[i] = self.coverages[i] - 1; // remove the coverage for ambiguous alleles
                    } else {
                        j += 1;
                    }
                }
            }
        }
        // All the pools needs be have been covered at least min_coverage times
        for c in &self.coverages {
            if c < &filter_stats.min_coverage {
                return Err(Error::new(ErrorKind::Other, "Filtered out."));
            }
        }
        // Filter by minimum allele frequency,
        // i.e. let maximum allele frequency be p: if 1-p < min_allele_frequency then remove
        let allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Cannot convert pileup line into allele frequencies.")),
        };
        let mean_allele_freqs = allele_frequencies.matrix.row_mean(); // summation across the rows which means sum of all elements per column
        if (1.00 - mean_allele_freqs.max()) < filter_stats.min_allele_frequency {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        Ok(self)
    }
}

pub fn pileup_to_sync(pileup_line: &mut PileupLine, filter_stats: &FilterStats) -> Option<String> {
    // let mut pileup_line: Box<PileupLine> = line.lparse().unwrap();
    // Filter
    match pileup_line.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Convert to counts
    let locus_counts = pileup_line.to_counts().unwrap();
    let (n, _p) = locus_counts.matrix.shape();
    // Instantiate the output line
    let mut x = vec![pileup_line.chromosome.clone(),
                                    pileup_line.position.to_string(),
                                    pileup_line.reference_allele.to_string()];
    for i in 0..n {
        x.push(locus_counts.matrix.row(i)
                                    .into_iter()
                                    .map(|w| w.to_string()).collect::<Vec<String>>()
                                    .join(":"));
    }
    let out = x.join("\t") + "\n";
    Some(out)
}

impl ChunkyReadAnalyseWrite<PileupLine, fn(&mut PileupLine, &FilterStats) -> Option<String>>
for FilePileup {
    fn per_chunk(&self, start: &u64, end: &u64, outname_ndigits: &usize, filter_stats: &FilterStats, function: fn(&mut PileupLine, &FilterStats) -> Option<String>) -> io::Result<String> {
        let fname = self.filename.clone();
        // Add leading zeros in front of the start file position so that we can sort the output files per chuck or thread properly
        let mut start_string = start.to_string();
        for  _i in 0..(outname_ndigits - start_string.len()) {
            start_string = "0".to_owned() + &start_string;
        }
        // Add leading zeros in front of the end file position so that we can sort the output files per chuck or thread properly
        let mut end_string = end.to_string();
        for  _i in 0..(outname_ndigits - end_string.len()) {
            end_string = "0".to_owned() + &end_string;
        }
        // Output temp file for the chunk    
        let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + ".sync.tmp";
        let out = fname_out.clone();
        let error_writing_file = "Unable to create file: ".to_owned() + &fname_out;
        let error_writing_line = "Unable to write line into file: ".to_owned() + &fname_out;
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
            let mut pileup_line: Box<PileupLine> = line.lparse()
                                                .expect(&("Input file error, i.e. '".to_owned() + &fname + "' at line with the first 20 characters as: " + &line[0..20] + "."));

            // Write the line
            let _ = match function(&mut pileup_line, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }
    
    fn read_analyse_write(&self, filter_stats: &FilterStats, out: &String, n_threads: &u64, function: fn(&mut PileupLine, &FilterStats) -> Option<String>) -> io::Result<String> {
        // Unpack pileup and pool names filenames
        let fname = self.filename.clone();
        let pool_names = self.pool_names.clone();
        // Output filename
        let mut out = out.to_owned();
        if out == "".to_owned() {
            let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
            let bname = fname.split(".").collect::<Vec<&str>>().into_iter().map(|a| a.to_owned())
                                        .collect::<Vec<String>>().into_iter().rev().collect::<Vec<String>>()[1..].to_owned()
                                        .into_iter().rev().collect::<Vec<String>>().join(".");
            out = bname.to_owned() + "-" + &time.to_string() + ".sync";
        }
        // Instatiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        // let mut file_out = File::create(&out).expect(&error_writing_file);
        let mut file_out = OpenOptions::new().create_new(true)
                                                .write(true)
                                                .append(false)
                                                .open(&out)
                                                .expect(&error_writing_file);
        // Pool names
        let mut names: Vec<String> = Vec::new();
        let file_names = match File::open(pool_names) {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Pool names file not found which is required for converting pileup into sync/syncx. Please make sure to include a valid file for the --pool-names parameter.")),
        };
        let reader_names = BufReader::new(file_names);
        for line in reader_names.lines() {
            names.push(line.unwrap().to_owned());
        }
        let names = names.join("\t");
        // // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let outname_ndigits = chunks[*n_threads as usize].to_string().len();
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
        // Making four separate threads calling the `search_for_word` function
        for i in 0..(*n_threads as usize) {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i].clone();
            let end = chunks[i+1].clone();
            let outname_ndigits = outname_ndigits.clone();
            let filter_stats = filter_stats.clone();
            let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
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
        file_out.write_all(("#chr\tpos\tref\t".to_owned() + &names + "\n").as_bytes()).unwrap();
        // Extract output filenames from each thread into a vector and sort them
        let mut fnames_out: Vec<String> = Vec::new();
        for f in thread_ouputs.lock().unwrap().iter() {
            fnames_out.push(f.to_owned());
        }
        println!("fnames_out={:?}", fnames_out);
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