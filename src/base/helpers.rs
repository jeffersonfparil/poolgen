use argmin::solver::neldermead::NelderMead;
use nalgebra::DVector;
use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, SeekFrom};
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
    return out;
}

pub fn find_file_splits(fname: &String, n_threads: &usize) -> io::Result<Vec<u64>> {
    let mut file = match File::open(fname) {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "The input file: ".to_owned() + fname + " does not exist. Please make sure you are entering the correct filename and/or the correct path.")),
    };
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.seek(SeekFrom::Current(0)).unwrap();
    let mut out = (0..end).step_by((end as usize) / n_threads).collect::<Vec<u64>>();
    out.push(end);
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
    }
    out.dedup();
    return Ok(out);
}

pub fn sensible_round(x: f64, n_digits: usize) -> f64 {
    let factor = ("1e".to_owned() + &n_digits.to_string())
        .parse::<f64>()
        .unwrap();
    (x * factor).round() / factor
}

pub fn parse_f64_roundup_and_own(x: f64, n_digits: usize) -> String {
    let s = x.to_string();
    if s.len() < n_digits {
        return s;
    }
    sensible_round(x, n_digits).to_string()
}

pub fn bound_parameters_with_logit(
    params: &Vec<f64>,
    lower_limit: f64,
    upper_limit: f64,
) -> Vec<f64> {
    // Map parameters with a logistic regression to bound them between 0 and PARAMETER_UPPER_LIMIT
    params
        .into_iter()
        .map(|x| lower_limit + ((upper_limit - lower_limit) / (1.00 + (-x).exp())))
        .collect::<Vec<f64>>()
}

pub fn prepare_solver_neldermead(p: f64, h: f64) -> NelderMead<Vec<f64>, f64> {
    let mut init_param: Vec<Vec<f64>> = Vec::new();
    for i in 0..(p as usize + 1) {
        init_param.push(vec![]);
        for j in 0..p as usize {
            if i == j {
                // init_param[i].push(1.5 * h)
                init_param[i].push(h + 0.5)
            } else {
                init_param[i].push(h)
            }
        }
    }
    NelderMead::new(init_param)
}

pub fn histogram(x: Vec<f64>, nbins: usize) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let max = x.iter().max_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let min = x.iter().min_by(|x, y| x.partial_cmp(y).unwrap()).unwrap();
    let bin_size = (max - min) / (nbins as f64);
    let mut bins_start: Vec<f64> = vec![];
    let mut bins_end: Vec<f64> = vec![];
    let mut counts: Vec<usize> = vec![];
    for i in 0..nbins {
        bins_start.push(sensible_round(min + ((i + 0) as f64 * bin_size), 7));
        bins_end.push(sensible_round(min + ((i + 1) as f64 * bin_size), 7));
        counts.push(
            x.iter()
                .map(|xi| {
                    if (*xi >= bins_start[i]) & (*xi < bins_end[i]) {
                        1 as usize
                    } else {
                        0 as usize
                    }
                })
                .sum(),
        );
    }
    counts[nbins - 1] += 1; // add the x.max()
    (bins_start, bins_end, counts)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_helpers() {
        let expected_output1: Vec<u64> = vec![0, 3563430, 7125955];
        let expected_output2: String = "0.42".to_owned();
        let expected_output3 = (
            vec![0.0, 0.2, 0.4, 0.6, 0.8],
            vec![0.2, 0.4, 0.6, 0.8, 1.0],
            vec![1, 2, 3, 4, 5],
        );
        // Inputs
        let fname: &String = &"./tests/test.pileup".to_owned();
        let n_threads = 2;
        let number: f64 = 0.420000012435;
        let betas = vec![
            0.0, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.75, 0.77, 0.8, 0.85, 0.86, 0.89, 1.0,
        ];
        // Output
        let splits = find_file_splits(fname, &n_threads).unwrap();
        let string_f64 = parse_f64_roundup_and_own(number, 4);
        let binning = histogram(betas, 5);
        // Assertion
        assert_eq!(expected_output1, splits);
        assert_eq!(expected_output2, string_f64);
        assert_eq!(expected_output3, binning);
    }
}
