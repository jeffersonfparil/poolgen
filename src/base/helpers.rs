use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader};
use std::io::{Error, ErrorKind};

// File splitting for thread allocation for parallele computation
fn find_start_of_next_line(fname: &String, pos: u64) -> u64 {
    let mut out = pos.clone();
    if out > 0 {
        let mut file = File::open(fname).unwrap();
        let _ = file.seek(SeekFrom::Start(out));
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let _ = reader.read_line(&mut line).unwrap();
        out = reader.seek(SeekFrom::Current(0)).unwrap();
    }
    return out
}

pub fn find_file_splits(fname: &String, n_threads: &u64) -> io::Result<Vec<u64>> {
    let mut file = match File::open(fname) {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "The input file: ".to_owned() + fname + " does not exist. Please make sure you are entering the correct filename and/or the correct path.")),
    };
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.seek(SeekFrom::Current(0)).unwrap();
    let mut out = (0..end).step_by((end/n_threads) as usize).collect::<Vec<u64>>();
    out.push(end);
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
    }
    out.dedup();
    return Ok(out)
}

pub fn parse_f64_roundup_and_own(x: f64, n_digits: usize) -> String {
    let s = x.to_string();
    if s.len() < n_digits {
        return s
    }
    s[0..n_digits].parse::<f64>().unwrap().to_string()
}

pub fn bound_parameters_with_logit(params: &Vec<f64>, lower_limit: f64, upper_limit: f64) -> Vec<f64> {
    // Map parameters with a logistic regression to bound them between 0 and PARAMETER_UPPER_LIMIT
    params.into_iter()
        .map(|x| lower_limit + 
                          ( (upper_limit-lower_limit) / 
                            (1.00 + (-x).exp()) 
                          ) 
            )
        .collect::<Vec<f64>>()
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_file_splits() {
        let expected_output1: Vec<u64> = vec![0, 3563430, 7125955];
        let expected_output2: String = "0.42".to_owned();
        // Inputs
        let fname: &String = &"./tests/test.pileup".to_owned();
        let n_threads: &u64 = &2;
        let number: f64 = 0.420000012435;
        // Output
        let splits = find_file_splits(fname, n_threads).unwrap();
        let string_f64 = parse_f64_roundup_and_own(number, 4);
        // Assertion
        assert_eq!(expected_output1, splits);
        assert_eq!(expected_output2, string_f64);
    }
}