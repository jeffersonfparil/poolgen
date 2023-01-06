use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use std::str;

struct PileupLine {
    chromosome: String,             // chromosome or scaffold name
    position: u64,                  // position in number of bases
    reference_allele: char,         // reference allele
    coverages: Vec<u64>,            // number of times the locus was covered
    read_codes: Vec<Vec<u8>>,       // utf8 read codes corresponding to 'A', 'T', 'C', or 'G' (252 other alleles can be accommodated)
    read_qualities: Vec<Vec<u8>>,   // utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

struct IndelMarker {
    indel: bool,    // insertion or deletion, i.e. +/- codes
    count: usize,   // size of the indel, i.e. how many bases
    left: usize,    // how many indel bases are left to be removed (initialised with the maximum possible value of 4294967295)
}

fn parse(line: &String) -> io::Result<PileupLine> {
    let raw_locus_data: Vec<&str>= line.split("\t").collect();
    // Chromosome or scaffold name
    let chromosome: String = raw_locus_data[0].to_string();
    // Position or locus coordinate in the genome assembly
    let position: u64 = raw_locus_data[1].to_string().parse::<u64>().unwrap();
    // Allele in the reference genome assembly
    let reference_allele: char = raw_locus_data[2].to_string().parse::<char>().unwrap();
    // List of the number of times the locus was read in each pool
    let mut coverages: Vec<u64> = vec![];
    for i in (3..raw_locus_data.len()).step_by(3) {
        coverages.push(raw_locus_data[i].to_string().parse::<u64>().unwrap());
    }
    // List of alleles that were read in each pool
    let mut read_codes: Vec<Vec<u8>> = vec![];
    for i in (4..raw_locus_data.len()).step_by(3) {
        if coverages[((i-1)/3)-1] > 0 {
                let raw_read_codes = raw_locus_data[i].as_bytes().to_owned();
                let mut alleles: Vec<u8> = vec![];
                let mut indel_marker: IndelMarker = IndelMarker{indel: false, count: 0, left: 4294967295};
                'per_pool: for j in 0..raw_read_codes.len() {
                    let code = raw_read_codes[j];
                    // Are we dealing with indels?
                    if indel_marker.indel {
                        // Find the firs digit of the number of indel/s
                        if (indel_marker.count == 0) & (indel_marker.left == 4294967295) {
                            indel_marker.count = str::from_utf8(&[code]).unwrap().parse::<usize>().unwrap();
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
                read_codes.push(vec![]);
            }
    }
    // println!("{:?}", read_codes);
    // List of the quailities of the read alleles in each pool
    let mut read_qualities: Vec<Vec<u8>> = vec![];
    for i in (5..raw_locus_data.len()).step_by(3) {
        if coverages[((i-1)/3)-1] > 0 {
                let qualities = raw_locus_data[i].as_bytes().to_owned();
                read_qualities.push(qualities);
            } else {
                read_qualities.push(vec![]);
            }
    }
    // Output PileupLine struct
    let out = PileupLine {
        chromosome: chromosome,
        position: position,
        reference_allele: reference_allele,
        coverages: coverages,
        read_codes: read_codes,
        read_qualities: read_qualities,
    };
    Ok(out)
}

pub fn read(fname: &str) -> io::Result<i32> {
    // let fname: &str = "/home/jeffersonfparil/Documents/poolgen/tests/test.pileup";
    let error_message = "File: ".to_owned() + fname + &" not found.".to_owned();
    let file = File::open(fname).expect(&error_message);
    let reader:BufReader<File> = BufReader::new(file);
    let mut i: i64 = 0;
    for line in reader.lines() {
        i += 1;
        let x = line.unwrap();
        let p = parse(&x).unwrap();
        if i < 5 {
            println!("{:?}", x);
            println!("chr:{:?}-pos:{:?}-ref:{:?}-covs:{:?}-reads:{:?}-quals:{:?}", p.chromosome, p.position, p.reference_allele, p.coverages, p.read_codes, p.read_qualities);
        }
    }
    Ok(0)
}
