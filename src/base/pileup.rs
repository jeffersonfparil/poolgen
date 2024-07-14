//! Pileup data (i.e. reference genome-aligned DNA sequences) processing (`Parse` and `Filter` traits and `pileup_to_sync` format conversion function) and parallel I/O (`ChunkyReadAnalyseWrite` trait)

use crate::base::*;
use ndarray::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, BufReader, BufWriter, Error, ErrorKind, SeekFrom};
use std::str;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

impl Parse<PileupLine> for String {
    /// Parse a line of pileup into a `PileupLine` struct corresponding to a locus, i.e. representing the allele counts of one or more pools in a single locus
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
            if coverages[((i - 1) / 3) - 1] > 0 {
                let raw_read_codes = raw_locus_data[i].as_bytes().to_owned();
                let mut alleles: Vec<u8> = Vec::new();
                let mut indel_marker: IndelMarker = IndelMarker {
                    indel: false,
                    count: 0,
                    left: 4294967295,
                }; // Instantiate the indel marker with no indel (i.e. indel==false and count==0) and the maximum number of indels left (i.e. left==4294967295 which is the maximum value for usize type of left)
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
                    if (code == 94) | (code == 36) {
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
                            65 => 65,  // A
                            97 => 65,  // a -> A
                            84 => 84,  // T
                            116 => 84, // t -> T
                            67 => 67,  // C
                            99 => 67,  // c -> C
                            71 => 71,  // G
                            103 => 71, // g -> G
                            42 => 68,  // * -> D
                            _ => 78,   // N
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
            if coverages[((i - 1) / 3) - 1] > 0 {
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
        return Ok(out);
    }
}

impl Filter for PileupLine {
    /// Parse the `PileupLine` into `AlleleCounts`
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let n: usize = self.coverages.len();
        let p: usize = 6;
        let mut matrix: Array2<u64> = Array2::from_elem((n, p), 0);
        let mut counts: Vec<Vec<u64>> = vec![
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ];
        let alleles_vector = vec!["A", "T", "C", "G", "D", "N"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        for pool in &self.read_codes {
            let mut counts_per_pool: Vec<u64> = vec![0, 0, 0, 0, 0, 0];
            for allele in pool {
                match allele {
                    65 => counts_per_pool[0] += 1, // A
                    84 => counts_per_pool[1] += 1, // T
                    67 => counts_per_pool[2] += 1, // C
                    71 => counts_per_pool[3] += 1, // G
                    68 => counts_per_pool[4] += 1, // D
                    _ => counts_per_pool[5] += 1,  // N
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
        Ok(Box::new(LocusCounts {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            alleles_vector: alleles_vector,
            matrix: matrix,
        }))
    }

    /// Parse `PileupLine` into `AlleleFrequencies`
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        let locus_counts = self.to_counts().unwrap();
        let n = locus_counts.matrix.nrows();
        let p = locus_counts.matrix.ncols();
        let row_sums = locus_counts.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
        let mut matrix: Array2<f64> = Array2::from_elem((n, p), 0.0);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = locus_counts.matrix[(i, j)] as f64 / row_sums[i] as f64;
            }
        }
        Ok(Box::new(LocusFrequencies {
            chromosome: locus_counts.chromosome,
            position: locus_counts.position,
            alleles_vector: locus_counts.alleles_vector,
            matrix: matrix,
        }))
    }

    /// Filter `PileupLine` by:
    /// - removing the entire locus if the locus is fixed, i.e. only 1 allele was found or retained after filterings
    /// Note that we are not removing alleles per locus if they fail the minimum allele frequency threshold, only if all alleles fail this threshold, i.e. when the locus is close to being fixed
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<Option<&mut Self>> {
        // First, make sure we have the correct the correct number of expected pools as the phenotype file
        // TODO: Make the error pop-out because it's currently being consumed as None in the only function calling it below.
        if self.coverages.len() != filter_stats.pool_sizes.len() {
            return Err(Error::new(ErrorKind::Other, "The number of pools in the pileup file does not correspond to the number of pools in the phenotype file."));
        }
        // Convert low quality bases into Ns
        let n = self.read_qualities.len();
        for i in 0..n {
            let mut j: usize = 0;
            while j < self.read_codes[i].len() {
                if self.read_qualities[i][j] < 33 {
                    return Err(Error::new(ErrorKind::Other, "Phred score out of bounds."));
                } else {
                    let q = f64::powf(10.0, -(self.read_qualities[i][j] as f64 - 33.0) / 10.0);
                    if q > filter_stats.max_base_error_rate {
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

        // Coverage depth and breadth requirement
        let min_coverage_breadth = (filter_stats.min_coverage_breadth * filter_stats.pool_sizes.len() as f64).ceil() as u32;
        let pools_covered = self.coverages.iter()
            .filter(|&&c| c >= filter_stats.min_coverage_depth)
            .take(min_coverage_breadth as usize)
            .count();

        if pools_covered != min_coverage_breadth as usize {
            return Ok(None)
        }

        // Remove monoallelic loci (each loci must have coverage of at least 2 alleles)
        if filter_stats.remove_monoallelic {
            let mut unique_alleles = Vec::new();

            'outer: for pool in &self.read_codes {
                for &read in pool {
                    if !unique_alleles.contains(&read) {
                        unique_alleles.push(read);
                        if unique_alleles.len() == 2 {
                            break 'outer;
                        }
                    }
                }
            }
            if unique_alleles.len() < 2 {
                return Ok(None)
            }
        }

        if filter_stats.keep_lowercase_reference {
            for pool in &mut self.read_codes {
                for read in pool.iter_mut() {
                    *read = match *read {
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
                }
            }
        }

        // Filter by minimum allele frequency,
        //// First convert the counts per pool into frequencies
        let allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert pileup line into allele frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let mut m = allele_frequencies.matrix.ncols();
        let mut q: f64;
        let mut j: usize = 1;
        // let mut matrix_new: Array2<f64>;
        // let mut alleles_vector_new: Vec<String>;
        while j < m {
            q = 0.0;
            for i in 0..n {
                q += allele_frequencies.matrix[(i, j)] * filter_stats.pool_sizes[i];
                // We've made sure the pool_sizes sum up to one in phen.rs
            }
            if (q < filter_stats.min_allele_frequency)
                | (q > (1.00 - filter_stats.min_allele_frequency))
            {
                // allele_frequencies.matrix = allele_frequencies.matrix.remove_column(j);
                m -= 1;
            } else {
                // if matrix_new.len() == 0 {
                //     matrix_new = allele_frequencies.matrix.slice(s![..,j..j]).to_owned();
                //     alleles_vector_new = vec![allele_frequencies.alleles_vector[j]];
                // } else {
                //     matrix_new = concatenate![Axis(0),  matrix_new, allele_frequencies.matrix.slice(s![..,j..j]).to_owned()];
                //     alleles_vector_new.push(allele_frequencies.alleles_vector[j]);
                // }
                j += 1;
            }
        }
        // Filter the whole locus depending on whether or not we have retained at least 2 alleles
        if m < 2 {
            return Ok(None)
        }
        Ok(Some(self))
    }
}

/// Convert `PileupLine` into a string representing a line or locus in a `sync` file.
pub fn pileup_to_sync(pileup_line: &mut PileupLine, filter_stats: &FilterStats) -> Option<String> {
    // let mut pileup_line: Box<PileupLine> = line.lparse().unwrap();
    // Filter
    match pileup_line.filter(filter_stats) {
        Ok(Some(x)) => x,
        _ => return None,
    };
    // Convert to counts
    let locus_counts = match pileup_line.to_counts() {
        Ok(x) => x,
        Err(_) => return None,
    };
    let n = locus_counts.matrix.nrows();
    // Instantiate the output line
    let mut x = vec![
        pileup_line.chromosome.clone(),
        pileup_line.position.to_string(),
        pileup_line.reference_allele.to_string(),
    ];
    for i in 0..n {
        x.push(
            locus_counts
                .matrix
                .row(i)
                .into_iter()
                .map(|w| w.to_string())
                .collect::<Vec<String>>()
                .join(":"),
        );
    }
    let out = x.join("\t") + "\n";
    Some(out)
}

impl ChunkyReadAnalyseWrite<PileupLine, fn(&mut PileupLine, &FilterStats) -> Option<String>>
    for FilePileup
{
    /// Load a `FilePileup`, process/analyse `PileupLine`s using some `function`, and output a file
    /// This function performs the I/O and processing/analyses for a chunk of the input file, such that each chunk is allocated to a dedicated thread for computational efficiency
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: fn(&mut PileupLine, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        let fname = self.filename.clone();
        // Add leading zeros in front of the start file position so that we can sort the output files per chuck or thread properly
        let mut start_string = start.to_string();
        for _i in 0..(outname_ndigits - start_string.len()) {
            start_string = "0".to_owned() + &start_string;
        }
        // Add leading zeros in front of the end file position so that we can sort the output files per chuck or thread properly
        let mut end_string = end.to_string();
        for _i in 0..(outname_ndigits - end_string.len()) {
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
            let mut pileup_line: Box<PileupLine> = line.lparse().expect(
                &("Input file error, i.e. '".to_owned()
                    + &fname
                    + "' at line with the first 20 characters as: "
                    + &line[0..20]
                    + "."),
            );

            // Write the line
            let _ = match function(&mut pileup_line, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }

    /// I/O and processing/analyses of a `FilePileup` using multiple threads for computational efficiency
    /// Each thread is allocated a chunk of the file, and the number of threads is set by the user with a default of 2.
    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &String,
        n_threads: &usize,
        function: fn(&mut PileupLine, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        // Unpack pileup and pool names filenames
        let fname = self.filename.clone();
        // Output filename
        let mut out = out.to_owned();
        if out == "".to_owned() {
            let time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs_f64();
            let bname = fname
                .split(".")
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|a| a.to_owned())
                .collect::<Vec<String>>()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()[1..]
                .to_owned()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()
                .join(".");
            out = bname.to_owned() + "-" + &time.to_string() + ".sync";
        }
        // Instatiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        // let mut file_out = File::create(&out).expect(&error_writing_file);
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&out)
            .expect(&error_writing_file);
        // Pool names
        let names = self.pool_names.join("\t");
        // // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let outname_ndigits = chunks[*n_threads].to_string().len();
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
                                                                                       // Making four separate threads calling the `search_for_word` function
        for i in 0..*n_threads {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i].clone();
            let end = chunks[i + 1].clone();
            let outname_ndigits = outname_ndigits.clone();
            let filter_stats = filter_stats.clone();
            let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let fname_out_per_thread = self_clone
                    .per_chunk(&start, &end, &outname_ndigits, &filter_stats, function)
                    .unwrap();
                thread_ouputs_clone
                    .lock()
                    .unwrap()
                    .push(fname_out_per_thread);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            let _ = thread.join().expect("Unknown thread error occured.");
        }
        file_out
            .write_all(("#chr\tpos\tref\t".to_owned() + &names + "\n").as_bytes())
            .unwrap();
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_pileup_methods() {
        let line = "Chromosome1\t456527\tC\t4\t....+1c\tJJJJ\t3\t.T.-3atg\tJJJ\t7\t.*.T..T\tJFJFJFJ\t5\tT....\tJJJJJ\t7\t...T...\tJJJJ<7J".to_owned();
        let pileup_line: PileupLine = *(line.lparse().unwrap());
        let counts = *(pileup_line.to_counts().unwrap());
        let frequencies = *(pileup_line.to_frequencies().unwrap());
        let filter_stats = FilterStats {
            remove_ns: true,
            remove_monoallelic: false,
            keep_lowercase_reference: false,
            max_base_error_rate: 0.005,
            min_coverage_depth: 1,
            min_coverage_breadth: 1.0,
            min_allele_frequency: 0.0,
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let mut filtered_pileup = pileup_line.clone();
        filtered_pileup.filter(&filter_stats).unwrap();
        // Assertions
        assert_eq!(
            pileup_line,
            PileupLine {
                chromosome: "Chromosome1".to_owned(),
                position: 456527,
                reference_allele: "C".to_owned().parse::<char>().unwrap(),
                coverages: vec![4, 3, 7, 5, 7],
                read_codes: vec![
                    vec![67, 67, 67, 67],
                    vec![67, 84, 67],
                    vec![67, 68, 67, 84, 67, 67, 84],
                    vec![84, 67, 67, 67, 67],
                    vec![67, 67, 67, 84, 67, 67, 67],
                ],
                read_qualities: vec![
                    vec![74, 74, 74, 74],
                    vec![74, 74, 74],
                    vec![74, 70, 74, 70, 74, 70, 74],
                    vec![74, 74, 74, 74, 74],
                    vec![74, 74, 74, 74, 60, 55, 74],
                ],
            }
        );
        assert_eq!(
            filtered_pileup,
            PileupLine {
                chromosome: "Chromosome1".to_owned(),
                position: 456527,
                reference_allele: "C".to_owned().parse::<char>().unwrap(),
                coverages: vec![4, 3, 7, 5, 6],
                read_codes: vec![
                    vec![67, 67, 67, 67],
                    vec![67, 84, 67],
                    vec![67, 68, 67, 84, 67, 67, 84],
                    vec![84, 67, 67, 67, 67],
                    vec![67, 67, 67, 84, 67, 67],
                ],
                read_qualities: vec![
                    vec![74, 74, 74, 74],
                    vec![74, 74, 74],
                    vec![74, 70, 74, 70, 74, 70, 74],
                    vec![74, 74, 74, 74, 74],
                    vec![74, 74, 74, 74, 60, 74],
                ],
            }
        );
        let counts_matrix: Array2<u64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                0, 0, 4, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 2, 4, 0, 1, 0, 0, 1, 4, 0, 0, 0, 0, 1, 6, 0,
                0, 0,
            ],
        )
        .unwrap();
        let mut frequencies_matrix: Array2<f64> =
            Array2::from_elem((counts_matrix.nrows(), counts_matrix.ncols()), 0.0);
        let row_sums: Vec<f64> = counts_matrix
            .sum_axis(Axis(1))
            .into_iter()
            .map(|x| x as f64)
            .collect::<Vec<f64>>();
        for i in 0..counts_matrix.nrows() {
            for j in 0..counts_matrix.ncols() {
                frequencies_matrix[(i, j)] = counts_matrix[(i, j)] as f64 / row_sums[i];
            }
        }
        assert_eq!(
            counts,
            LocusCounts {
                chromosome: "Chromosome1".to_owned(),
                position: 456527,
                alleles_vector: vec!["A", "T", "C", "G", "D", "N"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
                matrix: counts_matrix,
            }
        );
        assert_eq!(
            frequencies,
            LocusFrequencies {
                chromosome: "Chromosome1".to_owned(),
                position: 456527,
                alleles_vector: vec!["A", "T", "C", "G", "D", "N"]
                    .into_iter()
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>(),
                matrix: frequencies_matrix,
            }
        );
    }
}
