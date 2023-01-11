use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::result::Result;
// use std::error::Error;
use std::io::{Error, ErrorKind};
use std::str;

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

// Parse each line
fn parse(line: &String) -> Result<Box<PileupLine>, String> {
    let raw_locus_data: Box<Vec<&str>> = Box::new(line.split("\t").collect());
    // Chromosome or scaffold name
    let chromosome: String = raw_locus_data[0].to_string();
    // Position or locus coordinate in the genome assembly
    let position = match raw_locus_data[1].to_string().parse::<u64>() {
        Ok(x) => x,
        Err(_) => return Err("Please check the format of the input pileup file as position is not a valid integer (i.e. u64).".to_string()),
    };
    // Allele in the reference genome assembly
    let reference_allele  = match raw_locus_data[2].to_string().parse::<char>() {
        Ok(x) => x,
        Err(_) => return Err("Please check the format of the input pileup file as the reference allele is not a valid nucleotide base (i.e. not a valid single character).".to_string()),
    };
    // List of the number of times the locus was read in each pool
    let mut coverages: Vec<u64> = Vec::new();
    for i in (3..raw_locus_data.len()).step_by(3) {
        let cov = match raw_locus_data[i].to_string().parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err("Please check the format of the input pileup file as coverage field/s is/are not valid integer/s (i.e. u64) at pool: ".to_owned() + &(i/3).to_string() + &".".to_owned()),
        };
        coverages.push(cov);
    }
    // List of alleles that were read in each pool
    let mut read_codes: Vec<Vec<u8>> = Vec::new();
    for i in (4..raw_locus_data.len()).step_by(3) {
        // Parse if the current pool was read at least once for the current locus (i.e. line)
        if coverages[((i-1)/3)-1] > 0 {
            let error_message = "Check the codes for insertions and deletion, i.e. they must be integers after '+' and '-' at pool: ".to_owned() + &((i-1)/3).to_string() + &".".to_owned();
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
                            Err(_) => return Err(error_message),
                        };
                        continue 'per_pool;
                    }
                    // Find the next digit of the number of indels, if we have more than 9 indels
                    if (indel_marker.count > 0) & (indel_marker.left == 4294967295) {
                        let x = str::from_utf8(&[code]).unwrap().parse::<usize>();
                        let y = match x {
                            Ok(z) => z,
                            Err(_) => 0,
                        };
                        if y > 0 {
                            let z = indel_marker.count.to_string() + "0";
                            indel_marker.count = z.parse::<usize>().unwrap() + y;
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
                    continue 'per_pool;
                }
                // Reference allele (unicodes: ',' == 44 and '.' == 46) and 
                // non-reference alleles (unicodes: 'A'/'a' == 65/97,
                //                                  'T'/'t' == 84/116,
                //                                  'C'/'c' == 67/99,
                //                                  'G'/'g' == 71/103, and
                //                                  '*' == 42 codes for deletion on the current locus which is different from the \-[0-9]+[ACGTNacgtn]+ which encodes for indels in the next position/s)
                let a: u8 = match code {
                    44 => reference_allele.to_string().as_bytes()[0], // whatever the reference allele is
                    46 => reference_allele.to_string().as_bytes()[0], // whatever the reference allele is
                    65 => 65,   // A
                    97 => 65,   // a -> A
                    84 => 84,   // T
                    116 => 84,  // t -> T
                    67 => 67,   // C
                    99 => 67,   // c -> c
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
            return Err("Please check the format of the input pileup file as the coverages, number of read alleles and read qualities do not match at pool: ".to_owned() + &(i+1).to_string() + &".".to_owned());
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
                s += x - 33.0; // Assumes PHRED 33 (i.e. !=10^(-0/10) to I=10^(-40/10))
                n += 1.0;
            }
        }
        let out = f64::powf(10.0, -s/(10.0*n));
        // println!("{:?}", out);
        Ok(out)
    }

    fn filter(&mut self, min_quality: &f64) -> io::Result<&mut Self> {
        let n = &self.read_qualities.len();
        for i in 0..*n {
            let pool = &self.read_qualities[i];
            let m = pool.len();
            for j in 0..m {
                let q = f64::powf(10.0, -(pool[j] as f64 - 33.0) / 10.0); 
                if q > *min_quality {
                    self.read_codes[i][j] = 78; // convert to N
                }
            }
        }
        Ok(self)
    }

    fn reads_to_counts(&self, min_coverage: &u64) -> io::Result<Box<AlleleCounts>> {
        for c in &self.coverages {
            // println!("{:?}", c);
            if c < min_coverage {
                return Err(Error::new(ErrorKind::Other, "Filtered out."));
            }
        }
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
            let mut A: u64 = 0;
            let mut T: u64 = 0;
            let mut C: u64 = 0;
            let mut G: u64 = 0;
            let mut N: u64 = 0;
            let mut D: u64 = 0;    
            for allele in pool {
                match allele {
                    65 => A += 1,
                    84 => T += 1,
                    67 => C += 1,
                    71 => G += 1,
                    68 => D += 1,
                    _ => N += 1,
                };
            }
            out.a.push(A);
            out.t.push(T);
            out.c.push(C);
            out.g.push(G);
            out.n.push(N);
            out.d.push(D);
        }
        Ok(out)
    }
    
    fn reads_to_frequencies(&self, min_coverage: &u64) -> io::Result<Box<AlleleFrequencies>> {
        let mut out = Box::new(AlleleFrequencies {
            chromosome: self.chromosome.clone(),
            position: self.position.clone(),
            a: Vec::new(),
            t: Vec::new(),
            c: Vec::new(),
            g: Vec::new(),
            n: Vec::new(),
            d: Vec::new()});
        let counts = match self.reads_to_counts(min_coverage) {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Filtered out.")),
        };
        let n = counts.a.len();
        for i in 0..n {
            let sum = (counts.a[i] + counts.t[i] + counts.c[i] + counts.g[i] + counts.n[i] + counts.d[i]) as f64;
            out.a.push((counts.a[i] as f64) / sum);
            out.t.push((counts.t[i] as f64) / sum);
            out.c.push((counts.c[i] as f64) / sum); 
            out.g.push((counts.g[i] as f64) / sum); 
            out.n.push((counts.n[i] as f64) / sum); 
            out.d.push((counts.d[i] as f64) / sum); 
        }
        Ok(out)
    }

}

// Read pileup file
pub fn read(fname: &str, min_qual: &f64, min_cov: &u64) -> io::Result<Vec<Box<PileupLine>>> {
    // let fname: &str = "/home/jeffersonfparil/Documents/poolgen/tests/test.pileup";
    let file = File::open(fname).expect(&("Input file: '".to_owned() + fname + &"' not found.".to_owned()));
    let reader:BufReader<File> = BufReader::new(file);
    let mut i: i64 = 0;
    let mut out: Vec<Box<PileupLine>> = Vec::new();
    for line in reader.lines() {
        i += 1;
        let mut p = parse(&line.unwrap()).expect(&("Input file error, i.e. '".to_owned() + fname + &"' at line: ".to_owned() + &i.to_string() + &".".to_owned()));
        let q = p.mean_quality().unwrap();
        // println!("{}: Before: {:?}", i, p);
        p.filter(min_qual).unwrap();
        let r = match p.reads_to_counts(min_cov) {
            Ok(x) => x,
            Err(_) => Box::new(AlleleCounts { chromosome: "".to_owned(), position: 0, a: vec![0], t: vec![0], c: vec![0], g: vec![0], n: vec![0], d: vec![0] }),
        };
        let f = match p.reads_to_frequencies(min_cov) {
            Ok(x) => x,
            Err(_) => Box::new(AlleleFrequencies { chromosome: "".to_owned(), position: 0, a: vec![0.0], t: vec![0.0], c: vec![0.0], g: vec![0.0], n: vec![0.0], d: vec![0.0] }),
        };
        // println!("{:?}", p);
        // println!("{:?}", r.position);
        if r.position != 0 {
            println!("{}: Counts: {:?}", i, r);
        }
        // println!("{}: After: {:?}", i, p);
        // println!("idx: {}| read: {:?} | qualities: {:?}---{:?}", i, &r, &q, min_qual);
        if &q <= min_qual {
            out.push(p);
        }
    }
    Ok(out)
}
