//! Vcf data (i.e. variant calling format cannonical in human genetics) processing (`Parse` and `Filter` traits and `vcf_to_sync` format conversion function) and parallel I/O (`ChunkyReadAnalyseWrite` trait)
use crate::base::*;
use ndarray::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, BufReader, BufWriter, Error, ErrorKind, SeekFrom};
use std::str;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

impl Parse<VcfLine> for String {
    /// Parse a line of vcf into a `VcfLine` struct corresponding to a locus, i.e. representing the allele counts of one or more pools in a single locus
    fn lparse(&self) -> io::Result<Box<VcfLine>> {
        let raw_locus_data: Box<Vec<&str>> = Box::new(self.split("\t").collect());
        // Chromosome or scaffold name
        let chromosome: String = raw_locus_data[0].to_owned();
        // Position or locus coordinate in the genome assembly
        let position = match raw_locus_data[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Please check the format of the input vcf file as position is not a valid integer (i.e. u64).".to_owned())),
        };
        // Allele in the reference genome assembly
        let reference_allele = match raw_locus_data[3].to_owned().parse::<char>() {
            Ok(x) => x,
            Err(_) => 'D', // Set multi-character reference alleles as an indel (We cannot distinguish indels of different types at the moment, i.e. both AATGTG and AATCG will be classified as the same indel)
        };
        // Alternative alleles
        let alternative_alleles = raw_locus_data[4]
            .split(",")
            .map(|x| match x.parse::<char>() {
                Ok(x) => x,
                Err(_) => 'D', // Set multi-character reference alleles as an indel (We cannot distinguish indels of different types at the moment, i.e. both AATGTG and AATCG will be classified as the same indel)
            })
            .collect::<Vec<char>>();
        // Allele depths per pool
        let mut allele_depths = vec![];
        // First find the index of each genotype column (delimited by ":") referes to the allele depth, i.e. field "AD"
        let idx = raw_locus_data[8]
            .split(":")
            .enumerate()
            .filter(|(_idx, x)| x == &"AD")
            .map(|(idx, _x)| idx)
            .collect::<Vec<usize>>();
        let idx = if idx.len() == 1 {
            idx[0]
        } else {
            return Err(Error::new(ErrorKind::Other, "Please check the format of the input vcf file as the allele depths (AD attribute) were not generated.".to_owned()));
        };
        // Iterate across pools to extract the allele depths
        for i in 9..raw_locus_data.len() {
            allele_depths.push(
                raw_locus_data[i].split(":").collect::<Vec<&str>>()[idx]
                    .split(",")
                    .collect::<Vec<&str>>()
                    .iter()
                    .map(|&x| x.parse::<u64>().unwrap())
                    .collect::<Vec<u64>>(),
            );
        }
        // Output VcfLine struct
        let out = Box::new(VcfLine {
            chromosome: chromosome,
            position: position,
            reference_allele: reference_allele,
            alternative_alleles: alternative_alleles,
            allele_depths: allele_depths,
        });
        return Ok(out);
    }
}

impl Filter for VcfLine {
    /// Parse the `VcfLine` into `AlleleCounts`
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let mut alleles_vector = vec![self.reference_allele.to_string()];
        for a in &self.alternative_alleles {
            alleles_vector.push(a.to_string());
        }
        let n = self.allele_depths.len();
        let m = self.allele_depths[0].len(); // All pools should have the same number of alleles in the current locus
        let mut matrix: Array2<u64> = Array2::from_elem((n, m), 0);
        for i in 0..n {
            for j in 0..m {
                matrix[(i, j)] = self.allele_depths[i][j];
            }
        }
        Ok(Box::new(LocusCounts {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            alleles_vector: alleles_vector,
            matrix: matrix,
        }))
    }

    /// Parse `VcfLine` into `AlleleFrequencies`
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

    /// Filter `VcfLine` by:
    /// - removing the entire locus if the locus is fixed, i.e. only 1 allele was found or retained after filterings
    /// Note that we are not removing alleles per locus if they fail the minimum allele frequency threshold, only if all alleles fail this threshold, i.e. when the locus is close to being fixed
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // All the pools needs be have been covered at least min_coverage times
        for i in 0..self.allele_depths.len() {
            let c = &self.allele_depths[i].iter().fold(0, |sum, x| sum + x);
            if c < &filter_stats.min_coverage {
                return Err(Error::new(ErrorKind::Other, "Filtered out."));
            }
        }
        // Filter by minimum allele frequency,
        //// First convert the counts per pool into frequencies
        let allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert vcf line into allele frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let mut m = allele_frequencies.matrix.ncols();
        let mut q: f64;
        let mut j: usize = 1;
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
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        Ok(self)
    }
}

/// Convert `VcfLine` into a string representing a line or locus in a `sync` file.
pub fn vcf_to_sync(vcf_line: &mut VcfLine, filter_stats: &FilterStats) -> Option<String> {
    // Filter
    match vcf_line.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Convert to counts
    let locus_counts = match vcf_line.to_counts() {
        Ok(x) => x,
        Err(_) => return None,
    };
    let n = locus_counts.matrix.nrows();

    // Setup the A:T:C:G:DEL:N allele counts per pool (where DEL refers to an indel where we are limited to a single class, i.e. multiple different indels at a locus are currently all just classified as if they were a single type of variant class)
    let expected_alleles = vec!["A", "T", "C", "G", "D", "N"]
        .iter()
        .map(|&x| x.to_owned())
        .collect::<Vec<String>>();
    let mut atcgdn_matrix = Array2::from_elem((n, 6), 0);
    for i in 0..n {
        for j in 0..6 {
            for k in 0..locus_counts.alleles_vector.len() {
                if locus_counts.alleles_vector[k] == expected_alleles[j] {
                    atcgdn_matrix[(i, j)] = locus_counts.matrix[(i, k)];
                    break;
                }
            }
        }
    }

    // Instantiate the output line
    let mut x = vec![
        vcf_line.chromosome.clone(),
        vcf_line.position.to_string(),
        vcf_line.reference_allele.to_string(),
    ];
    for i in 0..n {
        x.push(
            atcgdn_matrix
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

impl ChunkyReadAnalyseWrite<VcfLine, fn(&mut VcfLine, &FilterStats) -> Option<String>> for FileVcf {
    /// Load a `FileVcf`, process/analyse `VcfLine`s using some `function`, and output a file
    /// This function performs the I/O and processing/analyses for a chunk of the input file, such that each chunk is allocated to a dedicated thread for computational efficiency
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: fn(&mut VcfLine, &FilterStats) -> Option<String>,
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
            // Ignore commented-out lines (i.e. '#' => 35)
            if line.as_bytes()[0] == 35 as u8 {
                continue;
            }
            // Parse the vcf line
            let mut vcf_line: Box<VcfLine> = line.lparse().expect(
                &("Input file error, i.e. '".to_owned()
                    + &fname
                    + "' at line with the first 20 characters as: "
                    + &line[0..20]
                    + "."),
            );

            // Write the line
            let _ = match function(&mut vcf_line, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }

    /// I/O and processing/analyses of a `FileVcf` using multiple threads for computational efficiency
    /// Each thread is allocated a chunk of the file, and the number of threads is set by the user with a default of 2.
    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &String,
        n_threads: &usize,
        function: fn(&mut VcfLine, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        // Unpack vcf and pool names filenames
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
        // Find the pool names from the last header line of the vcf file
        let file = File::open(fname.clone()).unwrap();
        let reader = BufReader::new(file);
        let mut pool_names: Vec<String> = vec![];
        for line in reader.lines() {
            let line = line.unwrap();
            if &line[0..6] == "#CHROM" {
                pool_names = line
                    .split("\t")
                    .enumerate()
                    .filter(|&(idx, _x)| idx > 8)
                    .map(|(_idx, x)| x.to_owned())
                    .collect();
            }
        }
        let names = if pool_names.len() > 0 {
            pool_names.join("\t")
        } else {
            return Err(Error::new(ErrorKind::Other, "Pool names not found, please check the header line of the vcf file. Make sure `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT...` header exists.".to_owned()));
        };
        // Find the positions where to split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let outname_ndigits = chunks[*n_threads].to_string().len();
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of vcf2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from vcf2sync_chunk()
        let thread_ouputs: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        // Mutated within each thread worker
        // Making four separate threads calling the `search_for_word` function
        for i in 0..*n_threads {
            // Clone vcf2sync_chunk parameters
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
    fn test_vcf_methods() {
        let line = "chrA\t323723\t.\tC\tA,CCT\t327.237\t.\tDP=469\tGT:PL:DP:AD\t0/0:0,24,239:7:7,0,1\t0/0:0,15,194:5:5,0,2\t0/1:42,0,61:6:4,2,1\t0/1:23,0,71:9:7,1,0\t0/1:75,0,132:14:6,5,2\t0/0:0,37,202:7:7,0,1\t0/0:0,15,235:5:5,0,1\t0/1:59,0,237:9:5,4,0\t0/0:0,27,123:9:9,0,0\t0/1:59,0,123:13:23,3,10".to_owned();
        let vcf_line: VcfLine = *(line.lparse().unwrap());
        println!("vcf_line={:?}", vcf_line);
        assert_eq!(
            vcf_line,
            VcfLine {
                chromosome: "chrA".to_owned(),
                position: 323723,
                reference_allele: 'C',
                alternative_alleles: vec!['A', 'D'],
                allele_depths: vec![
                    vec![7, 0, 1],
                    vec![5, 0, 2],
                    vec![4, 2, 1],
                    vec![7, 1, 0],
                    vec![6, 5, 2],
                    vec![7, 0, 1],
                    vec![5, 0, 1],
                    vec![5, 4, 0],
                    vec![9, 0, 0],
                    vec![23, 3, 10],
                ]
            }
        );
        let counts = *(vcf_line.to_counts().unwrap());
        println!("counts={:?}", counts);
        assert_eq!(
            counts,
            LocusCounts {
                chromosome: "chrA".to_owned(),
                position: 323723,
                alleles_vector: vec!["C".to_owned(), "A".to_owned(), "D".to_owned()],
                matrix: Array2::from_shape_vec(
                    (10, 3),
                    vec![
                        7, 0, 1, 5, 0, 2, 4, 2, 1, 7, 1, 0, 6, 5, 2, 7, 0, 1, 5, 0, 1, 5, 4, 0, 9,
                        0, 0, 23, 3, 10,
                    ]
                )
                .unwrap()
            }
        );
        let frequencies = *(vcf_line.to_frequencies().unwrap());
        println!("frequencies={:?}", frequencies);
        assert_eq!(
            frequencies,
            LocusFrequencies {
                chromosome: "chrA".to_owned(),
                position: 323723,
                alleles_vector: vec!["C".to_owned(), "A".to_owned(), "D".to_owned()],
                matrix: Array2::from_shape_vec(
                    (10, 3),
                    vec![
                        7. / 8.,
                        0. / 8.,
                        1. / 8.,
                        5. / 7.,
                        0. / 7.,
                        2. / 7.,
                        4. / 7.,
                        2. / 7.,
                        1. / 7.,
                        7. / 8.,
                        1. / 8.,
                        0. / 8.,
                        6. / 13.,
                        5. / 13.,
                        2. / 13.,
                        7. / 8.,
                        0. / 8.,
                        1. / 8.,
                        5. / 6.,
                        0. / 6.,
                        1. / 6.,
                        5. / 9.,
                        4. / 9.,
                        0. / 9.,
                        9. / 9.,
                        0. / 9.,
                        0. / 9.,
                        23. / 36.,
                        3. / 36.,
                        10. / 36.,
                    ]
                )
                .unwrap()
            }
        );
        let filter_stats_1 = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.0,
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2; 10],
        };
        let filter_stats_2 = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 10,
            min_allele_frequency: 0.0,
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2; 10],
        };
        let mut filtered_vcf_1 = vcf_line.clone();
        let mut filtered_vcf_2 = vcf_line.clone();
        println!(
            "filtered_vcf_1={:?}",
            filtered_vcf_1.filter(&filter_stats_1)
        );
        println!(
            "filtered_vcf_2={:?}",
            filtered_vcf_2.filter(&filter_stats_2)
        );
        filtered_vcf_1.filter(&filter_stats_1).unwrap();
        let err = match filtered_vcf_2.filter(&filter_stats_2) {
            Ok(_) => 0,
            Err(_) => 1,
        };
        assert_eq!(filtered_vcf_1, vcf_line);
        assert_eq!(err, 1);

        let mut vcf_line_2_sync = vcf_line.clone();
        let sync_line = vcf_to_sync(&mut vcf_line_2_sync, &filter_stats_1).unwrap();
        println!("sync_line={:?}", sync_line);
        assert_eq!(sync_line, "chrA\t323723\tC\t0:0:7:0:1:0\t0:0:5:0:2:0\t2:0:4:0:1:0\t1:0:7:0:0:0\t5:0:6:0:2:0\t0:0:7:0:1:0\t0:0:5:0:1:0\t4:0:5:0:0:0\t0:0:9:0:0:0\t3:0:23:0:10:0\n".to_owned());

        let file_vcf = FileVcf {
            filename: "./tests/test.vcf".to_owned(),
        };
        let output = file_vcf
            .read_analyse_write(
                &filter_stats_1,
                &"./tests/test-vcf2sync.sync".to_owned(),
                &2,
                vcf_to_sync,
            )
            .unwrap();
        println!("output={:?}", output);
        // assert_eq!(0, 1);
    }
}
