use crate::base::*;
use ndarray::prelude::*;
use std::fs::{File, OpenOptions};
use std::io::{self, prelude::*, BufReader, BufWriter, Error, ErrorKind, SeekFrom};
use std::str;
use std::sync::{Arc, Mutex};
use std::time::{SystemTime, UNIX_EPOCH};

impl CheckStruct for LocusCounts {
    fn check(&self) -> io::Result<()> {
        let (_n, p) = self.matrix.dim();
        let a = self.alleles_vector.len();
        match p == a {
            true => return Ok(()),
            false => {
                let locus = vec![self.chromosome.clone(), self.position.to_string()].join("-");
                return Err(Error::new(ErrorKind::Other, locus + " ::: LocusCounts error: the counts matrix does not have the same the number columns as the number of alleles."));
            }
        }
    }
}

impl CheckStruct for LocusFrequencies {
    fn check(&self) -> io::Result<()> {
        let (_n, p) = self.matrix.dim();
        let a = self.alleles_vector.len();
        match p == a {
            true => return Ok(()),
            false => {
                let locus = vec![self.chromosome.clone(), self.position.to_string()].join("-");
                return Err(Error::new(ErrorKind::Other, locus + " ::: LocusFrequencies error: the frequencies matrix does not have the same the number columns as the number of alleles."));
            }
        }
    }
}

impl CheckStruct for LocusCountsAndPhenotypes {
    fn check(&self) -> io::Result<()> {
        self.locus_counts.check().unwrap();
        let (n, _p) = self.locus_counts.matrix.dim();
        let (n_, _k) = self.phenotypes.dim();
        let n__ = self.pool_names.len();
        match (n == n_) & (n_ == n__) {
            true => return Ok(()),
            false => {
                let locus = vec![
                    self.locus_counts.chromosome.clone(),
                    self.locus_counts.position.to_string(),
                ]
                .join("-");
                return Err(Error::new(ErrorKind::Other, locus + " ::: LocusCountsAndPhenotypes error: the number of pools are inconsistent in the locus counts, and/or phenotypes matrix and/or pool names."));
            }
        }
    }
}

impl CheckStruct for GenotypesAndPhenotypes {
    fn check(&self) -> io::Result<()> {
        let p = self.chromosome.len();
        let p_ = self.position.len();
        let (n, p__) = self.intercept_and_allele_frequencies.dim();
        let (n_, _k) = self.phenotypes.dim();
        let n__ = self.pool_names.len();
        let (n___, l) = self.coverages.dim();
        match (p == p_) & (p_ == p__) & (n == n_) & (n_ == n__) & (n__ == n___) & (l <= p) {
            true => return Ok(()),
            false => return Err(Error::new(ErrorKind::Other, "GenotypesAndPhenotypes: there are at least 1 mismatch in the number of pools and loci."))
        }
    }
}

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
            return Err(Error::new(ErrorKind::Other, "Commented out line"));
        }
        // Parse the sync line
        let vec_line = line
            .split("\t")
            .collect::<Vec<&str>>()
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        let n: usize = vec_line.len() - 3;
        let p: usize = 6;
        let chromosome = vec_line[0].to_owned();
        let position = match vec_line[1].parse::<u64>() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Please check format of the file: position is not and integer.",
                ))
            }
        };
        let alleles_vector = vec!["A", "T", "C", "G", "N", "D"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        // Read the allele counts into the counts matrix
        let mut matrix: Array2<u64> = Array2::from_elem((n, p), 0);
        let mut counts: Vec<u64>;
        for i in 0..n {
            counts = vec_line[i+3].split(":")
                                .map(|x| x.to_string().parse::<u64>().expect("Please check the input sync file as the allele counts are not valid integers."))
                                .collect::<Vec<u64>>();
            for j in 0..p {
                matrix[(i, j)] = counts[j];
            }
        }
        return Ok(Box::new(LocusCounts {
            chromosome: chromosome,
            position: position,
            alleles_vector: alleles_vector,
            matrix: matrix, // n pools x 6 alleles
        }));
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
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let row_sums = self.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
                                                      // Make sure all pools have been convered
        if row_sums
            .iter()
            .fold(row_sums[0], |min, &x| if x < min { x } else { min })
            == 0
        {
            return Err(Error::new(
                ErrorKind::Other,
                "At least one pool have no coverage.",
            ));
        }
        let mut matrix: Array2<f64> = Array2::from_elem((n, p), 0.0 as f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = self.matrix[(i, j)] as f64 / row_sums[i] as f64;
            }
        }
        Ok(Box::new(LocusFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            alleles_vector: self.alleles_vector.clone(),
            matrix: matrix,
        }))
    }

    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        self.check().unwrap(); // preliminary check of the structure format
                               // Remove Ns
        if filter_stats.remove_ns {
            let i = match self
                .alleles_vector
                .iter()
                .position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned()))
            {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                self.matrix.remove_index(Axis(1), i as usize);
            }
        }
        // println!("self={:?}", self);
        // Filter by minimum coverage
        let mean_coverage = self.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
        let min_mean_coverage =
            mean_coverage
                .iter()
                .fold(mean_coverage[0], |min, &x| if x < min { x } else { min });
        if min_mean_coverage < filter_stats.min_coverage {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        // Filter by minimum allele frequency
        // Before anything else, we clone matrix of allele counts
        let mut matrix = self.matrix.clone();
        //// First convert allele counts into frequencies
        let mut allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let mut p = allele_frequencies.matrix.ncols();
        let mut q: f64;
        let mut j: usize = 0;
        while j < p {
            q = 0.0;
            for i in 0..n {
                q += allele_frequencies.matrix[(i, j)] * filter_stats.pool_sizes[i];
                // We've made sure the pool_sizes sum up to one in phen.rs
            }
            if (q < filter_stats.min_allele_frequency)
                | (q > (1.00 - filter_stats.min_allele_frequency))
            {
                allele_frequencies.matrix.remove_index(Axis(1), j);
                matrix.remove_index(Axis(1), j);
                self.alleles_vector.remove(j);
                p -= 1;
            } else {
                j += 1;
            }
        }
        if p < 2 {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        self.matrix = matrix;
        Ok(self)
    }
}

impl Filter for LocusFrequencies {
    // PileupLine to AlleleCounts
    fn to_counts(&self) -> io::Result<Box<LocusCounts>> {
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let mut matrix: Array2<u64> = Array2::from_elem((n, p), 0);
        let mut max_n: f64;
        for i in 0..n {
            let row = self.matrix.row(i);
            let min = row
                .iter()
                .filter(|&&x| x != 0.0)
                .fold(row[0], |min, &x| if x < min { x } else { min });
            if min == 0.0 {
                return Err(Error::new(
                    ErrorKind::Other,
                    "At least one of the pools have no coverage.",
                ));
            }
            max_n = 1.00 / min;
            for j in 0..p {
                matrix[(i, j)] = (max_n * self.matrix[(i, j)]).round() as u64;
            }
        }
        Ok(Box::new(LocusCounts {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            alleles_vector: self.alleles_vector.clone(),
            matrix: matrix,
        }))
    }

    // PileupLine to AlleleFrequencies
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>> {
        // Recompute the frequencies using frequencies when the number of colulmns or one or more alleles have been filtered out/removed
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let row_sums = self.matrix.sum_axis(Axis(1)); // summation across the columns which means sum of all elements per row
                                                      // Make sure all pools have been convered
        if row_sums
            .iter()
            .fold(row_sums[0], |min, &x| if x < min { x } else { min })
            == 0.0
        {
            return Err(Error::new(
                ErrorKind::Other,
                "At least one pool did not have coverage.",
            ));
        }
        let mut matrix: Array2<f64> = Array2::from_elem((n, p), 0.0 as f64);
        for i in 0..n {
            for j in 0..p {
                matrix[(i, j)] = self.matrix[(i, j)] / row_sums[i];
            }
        }
        Ok(Box::new(LocusFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            alleles_vector: self.alleles_vector.clone(),
            matrix: matrix,
        }))
    }

    // Filter PileupLine by minimum coverage, minimum quality
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self> {
        // Cannot filter by base qualities as this information is lost and we are assuming this has been performed during pileup to sync conversion
        // Also, cannot filter by minimum coverage as that data is lost from counts to frequencies conversion
        self.check().unwrap(); // preliminary check of the structure format
                               // Remove Ns
        if filter_stats.remove_ns {
            let i = match self
                .alleles_vector
                .iter()
                .position(|x| (x == &"N".to_owned()) | (x == &"n".to_owned()))
            {
                Some(x) => x as i32,
                None => -1,
            };
            if i != -1 {
                self.alleles_vector.remove(i as usize);
                self.matrix.remove_index(Axis(1), i as usize);
            }
        }
        // println!("self={:?}", self);
        // Recompute frequencies after removing Ns
        let recomputed_self = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "T_T Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        self.alleles_vector = recomputed_self.alleles_vector;
        self.matrix = recomputed_self.matrix;
        // Filter by minimum allele frequency
        // Before anything else, we clone matrix of allele frequencies
        let mut matrix = self.matrix.clone();
        //// First convert allele counts into frequencies
        let mut allele_frequencies = match self.to_frequencies() {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "Cannot convert locus counts to locus frequencies.",
                ))
            }
        };
        //// Next account for pool sizes to get the proper minmum allele frequency across all pools
        let n = allele_frequencies.matrix.nrows();
        let p_ = allele_frequencies.matrix.ncols();
        let mut p = p_;
        let mut q: f64;
        let mut j: usize = 0;
        while j < p {
            q = 0.0;
            for i in 0..n {
                q += allele_frequencies.matrix[(i, j)] * filter_stats.pool_sizes[i];
                // We've made sure the pool_sizes sum up to one in phen.rs
            }
            if (q < filter_stats.min_allele_frequency)
                | (q > (1.00 - filter_stats.min_allele_frequency))
            {
                allele_frequencies.matrix.remove_index(Axis(1), j);
                matrix.remove_index(Axis(1), j);
                self.alleles_vector.remove(j);
                p -= 1;
            } else {
                j += 1;
            }
        }
        if p < 2 {
            return Err(Error::new(ErrorKind::Other, "Filtered out."));
        }
        self.matrix = matrix;
        Ok(self)
    }
}

impl Sort for LocusFrequencies {
    fn sort_by_allele_freq(&mut self, decreasing: bool) -> io::Result<&mut Self> {
        let n = self.matrix.nrows();
        let p = self.matrix.ncols();
        let mut sorted_matrix: Array2<f64> = Array2::from_elem((n, p), f64::NAN);
        let mut sorted_alleles_vector: Vec<String> = vec![];
        let mut idx = (0..self.matrix.ncols()).collect::<Vec<usize>>();
        let column_sums = self.matrix.sum_axis(Axis(0));
        // println!("self={:?}", self);
        if decreasing {
            idx.sort_by(|&a, &b| column_sums[b].partial_cmp(&column_sums[a]).unwrap());
        } else {
            idx.sort_by(|&a, &b| column_sums[a].partial_cmp(&column_sums[b]).unwrap());
        }
        for i in 0..p {
            sorted_matrix
                .column_mut(i)
                .assign(&self.matrix.column(idx[i]));
            sorted_alleles_vector.push(self.alleles_vector[idx[i]].clone());
        }
        self.matrix = sorted_matrix;
        self.alleles_vector = sorted_alleles_vector;
        Ok(self)
    }
}

impl ChunkyReadAnalyseWrite<LocusCounts, fn(&mut LocusCounts, &FilterStats) -> Option<String>>
    for FileSync
{
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: fn(&mut LocusCounts, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
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
        let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + ".tmp";
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
                    _ => {
                        return Err(Error::new(
                            ErrorKind::Other,
                            "T_T Input sync file error, i.e. '".to_owned()
                                + &fname
                                + "' at line with the first 20 characters as: "
                                + &line[0..20]
                                + ".",
                        ))
                    }
                },
            };
            // Write the line
            let _ = match function(&mut locus_counts, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }

    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &String,
        n_threads: &usize,
        function: fn(&mut LocusCounts, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        // Unpack pileup and pool names filenames
        let fname = self.filename.clone();
        let test = self.test.clone();
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
            out = bname.to_owned() + "-" + &time.to_string() + "-" + &test + ".csv";
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
        // // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        let n_chunks = chunks.len();
        let n_threads = if n_chunks < *n_threads {
            &n_chunks
        } else {
            n_threads
        };
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
        // Write out
        file_out
            .write_all(("#chr,pos,alleles,statistic,pvalue\n").as_bytes())
            .unwrap();
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

impl
    ChunkyReadAnalyseWrite<
        LocusCountsAndPhenotypes,
        fn(&mut LocusCountsAndPhenotypes, &FilterStats) -> Option<String>,
    > for FileSyncPhen
{
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: fn(&mut LocusCountsAndPhenotypes, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        let fname = self.filename_sync.clone();
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
        let fname_out = fname.to_owned() + "-" + &start_string + "-" + &end_string + ".tmp";
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
            let locus_counts: Box<LocusCounts> = match line.lparse() {
                Ok(x) => x,
                Err(x) => match x.kind() {
                    ErrorKind::Other => continue,
                    _ => {
                        return Err(Error::new(
                            ErrorKind::Other,
                            "T_T Input sync file error, i.e. '".to_owned()
                                + &fname
                                + "' at line with the first 20 characters as: "
                                + &line[0..20]
                                + ".",
                        ))
                    }
                },
            };
            let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
                locus_counts: *locus_counts.clone(),
                pool_names: self.pool_names.clone(),
                phenotypes: self.phen_matrix.clone(),
            };
            // Write the line
            let _ = match function(&mut locus_counts_and_phenotypes, filter_stats) {
                Some(x) => file_out.write_all(x.as_bytes()).expect(&error_writing_line),
                None => continue,
            };
        }
        Ok(out)
    }

    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &String,
        n_threads: &usize,
        function: fn(&mut LocusCountsAndPhenotypes, &FilterStats) -> Option<String>,
    ) -> io::Result<String> {
        // Unpack pileup and pool names filenames
        let fname = self.filename_sync.clone();
        let test = self.test.clone();
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
            out = bname.to_owned() + "-" + &time.to_string() + "-" + &test + ".csv";
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
        // Write out
        file_out
            .write_all(("#chr,pos,alleles,freq,pheno,statistic,pvalue\n").as_bytes())
            .unwrap();
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

impl LoadAll for FileSyncPhen {
    fn per_chunk_load(
        &self,
        start: &u64,
        end: &u64,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        // Input syn file
        let fname = self.filename_sync.clone();

        // Prepare output vectors
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let mut cnts: Vec<LocusCounts> = Vec::new();
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
            let mut locus_counts: LocusCounts = match line.lparse() {
                Ok(x) => *x,
                Err(x) => match x.kind() {
                    ErrorKind::Other => continue,
                    _ => {
                        return Err(Error::new(
                            ErrorKind::Other,
                            "T_T Input sync file error, i.e. '".to_owned()
                                + &fname
                                + "' at line with the first 20 characters as: "
                                + &line[0..20]
                                + ".",
                        ))
                    }
                },
            };
            match locus_counts.filter(filter_stats) {
                Ok(x) => x,
                Err(_) => continue,
            };
            let mut locus_frequencies = *locus_counts.to_frequencies().unwrap();
            // let mut locus_frequencies = *locus_counts.to_frequencies().unwrap();
            // // println!("locus_frequencies={:?}", locus_frequencies);
            // match locus_frequencies.filter(filter_stats) {
            //     Ok(x) => x,
            //     Err(_) => continue,
            // };
            // Remove minimum allele
            if keep_p_minus_1 {
                locus_frequencies.sort_by_allele_freq(true).unwrap();
                locus_frequencies.matrix.remove_index(Axis(1), 0);
                locus_frequencies.alleles_vector.remove(0);
            }
            freq.push(locus_frequencies);
            cnts.push(locus_counts);
        }
        Ok((freq, cnts))
    }

    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)> {
        let fname = self.filename_sync.clone();
        // Find the positions whereto split the file into n_threads pieces
        let chunks = find_file_splits(&fname, n_threads).unwrap();
        println!("Chunks: {:?}", chunks);
        // Tuple arguments of pileup2sync_chunks
        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs_freq: Arc<Mutex<Vec<LocusFrequencies>>> =
            Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
        let thread_ouputs_cnts: Arc<Mutex<Vec<LocusCounts>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
                                                                                                 // Making four separate threads calling the `search_for_word` function
        for i in 0..*n_threads {
            // Clone pileup2sync_chunk parameters
            let self_clone = self.clone();
            let start = chunks[i].clone();
            let end = chunks[i + 1].clone();
            let filter_stats = filter_stats.clone();
            let thread_ouputs_freq_clone = thread_ouputs_freq.clone(); // Mutated within the current thread worker
            let thread_ouputs_cnts_clone = thread_ouputs_cnts.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let (mut freq, mut cnts) = self_clone
                    .per_chunk_load(&start, &end, &filter_stats, keep_p_minus_1)
                    .unwrap();
                thread_ouputs_freq_clone.lock().unwrap().append(&mut freq);
                thread_ouputs_cnts_clone.lock().unwrap().append(&mut cnts);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            let _ = thread.join().expect("Unknown thread error occured.");
        }
        // Extract output filenames from each thread into a vector and sort them
        let mut freq: Vec<LocusFrequencies> = Vec::new();
        let mut cnts: Vec<LocusCounts> = Vec::new();
        for x in thread_ouputs_freq.lock().unwrap().iter() {
            freq.push(x.clone());
        }
        for x in thread_ouputs_cnts.lock().unwrap().iter() {
            cnts.push(x.clone());
        }
        freq.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });
        cnts.sort_by(|a, b| {
            a.chromosome
                .cmp(&b.chromosome)
                .then(a.position.cmp(&b.position))
        });

        Ok((freq, cnts))
    }

    fn write_csv(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        out: &String,
        n_threads: &usize,
    ) -> io::Result<String> {
        // Output filename
        let out = if *out == "".to_owned() {
            let time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_secs_f64();
            let bname = self
                .filename_sync
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
            bname.to_owned() + "-" + &time.to_string() + "-allele_frequencies.csv"
        } else {
            out.clone()
        };
        // Instatiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &out;
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&out)
            .expect(&error_writing_file);
        file_out
            .write_all(
                ("#chr,pos,allele,".to_owned() + &self.pool_names.join(",") + "\n").as_bytes(),
            )
            .unwrap();
        // Load the full sync file in parallel and sort
        let (freqs, _cnts) = self.load(filter_stats, keep_p_minus_1, n_threads).unwrap();
        for f in freqs.iter() {
            for i in 0..f.alleles_vector.len() {
                let freqs_per_pool = f
                    .matrix
                    .column(i)
                    .iter()
                    .map(|x| parse_f64_roundup_and_own(*x, 6))
                    .collect::<Vec<String>>()
                    .join(",");
                let line = vec![
                    f.chromosome.to_owned(),
                    f.position.to_string(),
                    f.alleles_vector[i].to_owned(),
                    freqs_per_pool,
                ]
                .join(",")
                    + "\n";
                file_out.write_all(line.as_bytes()).unwrap();
            }
        }

        Ok(out)
    }

    fn into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes> {
        let (freqs, cnts) = self.load(filter_stats, keep_p_minus_1, n_threads).unwrap();
        let n = self.pool_names.len();
        let m = freqs.len(); // total number of loci
                             // Find the total number of alleles across all loci
        let mut p = 1; // start with the intercept
        for f in freqs.iter() {
            p += f.matrix.ncols();
        }
        // println!("p={}", p);
        let mut chromosome: Vec<String> = Vec::with_capacity(p);
        chromosome.push("intercept".to_owned());
        let mut position: Vec<u64> = Vec::with_capacity(p);
        position.push(0);
        let mut allele: Vec<String> = Vec::with_capacity(p);
        allele.push("intercept".to_owned());
        let mut coverages: Array2<f64> = Array2::from_elem((n, m), f64::NAN);
        let mut l: usize = 0; // locus index
        let mut mat: Array2<f64> = Array2::from_elem((n, p), 1.0);
        let mut j: usize = 1; // SNP index across loci, start after the intercept
        assert_eq!(
            freqs.len(),
            cnts.len(),
            "Frequencies and counts not the same length."
        );
        for i in 0..freqs.len() {
            // Allele frequencies
            let f = &freqs[i];
            for j_ in 0..f.matrix.ncols() {
                chromosome.push(f.chromosome.clone());
                position.push(f.position);
                allele.push(f.alleles_vector[j_].clone());
                for i in 0..f.matrix.nrows() {
                    mat[(i, j)] = f.matrix[(i, j_)];
                }
                j += 1; // next allele
            }
            // Coverages
            let c = &cnts[i];
            let cov: Array1<f64> = c.matrix.sum_axis(Axis(1)).map(|&x| x as f64);
            for l_ in 0..cov.len() {
                coverages[(l_, l)] = cov[l_];
            }
            l += 1; // next locus
        }
        // println!("mat={:?}", mat.slice(s![0..5, 0..4]));
        // println!("chromosome[0]={:?}", chromosome[0]);
        // println!("chromosome[1]={:?}", chromosome[1]);
        // println!("chromosome[2]={:?}", chromosome[2]);
        // println!("chromosome[3]={:?}", chromosome[3]);
        // println!("position[0]={:?}", position[0]);
        // println!("position[1]={:?}", position[1]);
        // println!("position[2]={:?}", position[2]);
        // println!("position[3]={:?}", position[3]);
        Ok(GenotypesAndPhenotypes {
            chromosome: chromosome,
            position: position,
            allele: allele,
            intercept_and_allele_frequencies: mat,
            phenotypes: self.phen_matrix.clone(),
            pool_names: self.pool_names.clone(),
            coverages: coverages,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use rand::prelude::Distribution;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_sync() {
        // (01) CheckStruct: LocusCounts
        let mut locus_counts = LocusCounts {
            chromosome: "chr1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned()], // Missing an allele in alleles_vector
            matrix: Array2::from_shape_vec((5, 2), vec![10, 90, 20, 80, 50, 50, 80, 20, 90, 10])
                .unwrap(),
        };
        // Missing an allele in alleles_vector
        assert_eq!(true, locus_counts.check().is_err());
        // Append the missing allele
        locus_counts.alleles_vector.push("D".to_owned());
        // No error
        assert_eq!((), locus_counts.check().unwrap());
        // (02) CheckStruct: LocusFrequencies
        let mut locus_frequncies = LocusFrequencies {
            chromosome: "chr1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned()], // Missing an allele in alleles_vector
            matrix: Array2::from_shape_vec(
                (5, 2),
                vec![0.1, 0.9, 0.2, 0.8, 0.5, 0.5, 0.8, 0.2, 0.9, 0.1],
            )
            .unwrap(),
        };
        // Missing an allele in alleles_vector
        assert_eq!(true, locus_frequncies.check().is_err());
        // Append the missing allele
        locus_frequncies.alleles_vector.push("D".to_owned());
        // No error
        assert_eq!((), locus_frequncies.check().unwrap());
        // (03) CheckStruct: LocusCountsAndPhenotypes
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: Array2::from_shape_vec((4, 1), vec![0.1, 0.2, 0.3, 0.4]).unwrap(),
            pool_names: (1..5)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(),
        };
        // Missing phenotype, and pool name
        assert_eq!(true, locus_counts_and_phenotypes.check().is_err());
        // Replace phenotype array with the complete data
        locus_counts_and_phenotypes.phenotypes =
            Array2::from_shape_vec((5, 1), vec![0.1, 0.2, 0.3, 0.4, 0.5]).unwrap();
        // Missing pool name
        assert_eq!(true, locus_counts_and_phenotypes.check().is_err());
        // Append missing pool name
        locus_counts_and_phenotypes
            .pool_names
            .push("Pop5".to_owned());
        // No error
        assert_eq!((), locus_counts_and_phenotypes.check().unwrap());
        // (04) CheckStruct: GenotypesAndPhenotypes
        let mut rng = rand::thread_rng();
        let dist_unif = statrs::distribution::Uniform::new(0.0, 1.0).unwrap();
        let dist_gaus = statrs::distribution::Normal::new(0.0, 1.0).unwrap();
        let mut genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: vec![
                "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chrX", "chrX", "chrX",
            ]
            .iter()
            .map(|&x| x.to_owned())
            .collect(),
            position: vec![123, 456, 789, 10005, 10010, 10015, 701, 702],
            allele: vec!["A", "T", "C", "A", "T", "C", "A", "T"]
                .iter()
                .map(|&x| x.to_owned())
                .collect(),
            intercept_and_allele_frequencies: Array2::from_shape_vec(
                (5, 9),
                dist_unif.sample_iter(rng.clone()).take(5 * 9).collect(),
            )
            .unwrap(),
            phenotypes: Array2::from_shape_vec(
                (4, 2),
                dist_gaus.sample_iter(rng.clone()).take(4 * 2).collect(),
            )
            .unwrap(),
            pool_names: (1..5)
                .map(|x| "Pop".to_owned() + &x.to_string()[..])
                .collect(),
            coverages: Array2::from_shape_vec(
                (4, 2),
                std::iter::repeat(100.0).take(4 * 2).collect(),
            )
            .unwrap(),
        };

        // Expected output
        let counts_matrix: Array2<u64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                1, 0, 999, 0, 4, 0, 0, 1, 2, 0, 0, 0, 0, 2, 4, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 1, 6,
                0, 0, 0,
            ],
        )
        .unwrap();
        let mut frequencies_matrix: Array2<f64> =
            Array2::zeros((counts_matrix.nrows(), counts_matrix.ncols()));
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
        // Notice the difference in the default order of alleles between PileupLine-derived counts and frequences and Sync-derived counts and frequencies, i.e. "D" before "N" in the former.
        let expected_output1 = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 456527,
            alleles_vector: vec!["A", "T", "C", "G", "N", "D"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            matrix: counts_matrix,
        };
        let expected_output2 = LocusFrequencies {
            chromosome: "Chromosome1".to_owned(),
            position: 456527,
            alleles_vector: vec!["A", "T", "C", "G", "N", "D"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            matrix: frequencies_matrix,
        };
        let mut expected_output3 = expected_output1.clone();
        expected_output3.alleles_vector = vec!["T", "C"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        expected_output3.matrix.remove_index(Axis(1), 0);
        expected_output3.matrix.remove_index(Axis(1), 2);
        expected_output3.matrix.remove_index(Axis(1), 2);
        expected_output3.matrix.remove_index(Axis(1), 2);
        let mut expected_output4 = expected_output2.clone();
        expected_output4.alleles_vector = vec!["T", "C"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        expected_output4.matrix.remove_index(Axis(1), 0);
        expected_output4.matrix.remove_index(Axis(1), 2);
        expected_output4.matrix.remove_index(Axis(1), 2);
        expected_output4.matrix.remove_index(Axis(1), 2);
        let expected_output4 = *(expected_output4.to_frequencies().unwrap());
        let mut expected_output5 = expected_output4.clone();
        expected_output5.alleles_vector = vec!["C", "T"]
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<String>>();
        let mut new_freqs = expected_output5
            .matrix
            .column(1)
            .into_iter()
            .map(|x| x.to_owned())
            .collect::<Vec<f64>>();
        for i in 0..5 {
            new_freqs.push(expected_output5.matrix[(i, 0)]);
        }
        expected_output5.matrix = Array2::from_shape_vec((2, 5), new_freqs)
            .unwrap()
            .reversed_axes();
        let expected_output6 = LocusFrequencies {
            chromosome: "Chromosome1".to_owned(),
            position: 456527,
            alleles_vector: ["T"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            matrix: Array2::from_shape_vec(
                (1, 5),
                vec![
                    0.0,
                    0.3333333333333333,
                    0.3333333333333333,
                    0.2,
                    0.14285714285714285,
                ],
            )
            .unwrap()
            .reversed_axes(),
        };
        let expected_output7 = LocusFrequencies {
            chromosome: "Chromosome1".to_owned(),
            position: 1041321,
            alleles_vector: ["G"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            matrix: Array2::from_shape_vec(
                (1, 5),
                vec![
                    0.047619047619047616,
                    0.0,
                    0.0625,
                    0.043478260869565216,
                    0.037037037037037035,
                ],
            )
            .unwrap()
            .reversed_axes(),
        };
        // Inputs
        let line = "Chromosome1\t456527\tC\t1:0:999:0:4:0\t0:1:2:0:0:0\t0:2:4:0:0:0\t0:1:4:0:0:0\t0:1:6:0:0:0".to_owned();
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "load".to_owned(),
        };
        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        // Outputs
        let counts: LocusCounts = *(line.lparse().unwrap());
        let frequencies = *(counts.to_frequencies().unwrap());
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let mut filtered_counts = counts.clone();
        filtered_counts.filter(&filter_stats).unwrap();
        let filtered_frequencies = *(frequencies
            .clone()
            .filter(&filter_stats)
            .unwrap()
            .to_frequencies()
            .unwrap());
        let mut sorted_filtered_frequencies = filtered_frequencies.clone();
        sorted_filtered_frequencies
            .sort_by_allele_freq(true)
            .unwrap();
        let n_threads = 2;
        let (loaded_freqs, loaded_counts) = file_sync_phen
            .load(&filter_stats, true, &n_threads)
            .unwrap();
        let frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, true, &n_threads)
            .unwrap();
        // println!("loaded_freqs={:?}", loaded_freqs);
        // println!("len(loaded_freqs)={:?}", loaded_freqs.len());
        // Assertions
        assert_eq!(expected_output1, counts);
        assert_eq!(expected_output2, frequencies);
        assert_eq!(expected_output3, filtered_counts);
        assert_eq!(expected_output4, filtered_frequencies);
        assert_eq!(expected_output5, sorted_filtered_frequencies);
        assert_eq!(expected_output6, loaded_freqs[0]);
        assert_eq!(
            expected_output6.chromosome,
            frequencies_and_phenotypes.chromosome[1]
        );
        assert_eq!(
            expected_output6.position,
            frequencies_and_phenotypes.position[1]
        );
        println!(
            "frequencies_and_phenotypes={:?}",
            frequencies_and_phenotypes
        );
        assert_eq!(
            expected_output6.matrix[(0, 0)],
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)]
        );
        // TO DO: IMPROVE UNIT TESTS!!!! 2023/06/09
        // assert_eq!(
        //     expected_output6.matrix[(1, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(2, 0)]
        // );
        // assert_eq!(
        //     expected_output6.matrix[(2, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 0)]
        // );
        // assert_eq!(
        //     expected_output6.matrix[(3, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(4, 0)]
        // );
        // assert_eq!(
        //     expected_output6.matrix[(4, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(5, 0)]
        // );
        // assert_eq!(
        //     expected_output7.matrix[(0, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)]
        // );
        // assert_eq!(
        //     expected_output7.matrix[(1, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 1)]
        // );
        // assert_eq!(
        //     expected_output7.matrix[(2, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(2, 1)]
        // );
        // assert_eq!(
        //     expected_output7.matrix[(3, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 1)]
        // );
        // assert_eq!(
        //     expected_output7.matrix[(4, 0)],
        //     frequencies_and_phenotypes.intercept_and_allele_frequencies[(4, 1)]
        // );
    }
}
