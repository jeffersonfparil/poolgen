use std;
use std::fs::File;
use std::io::{prelude::*, SeekFrom, BufReader};

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

pub fn find_file_splits(fname: &String, n_threads: &u64) -> Vec<u64> {
    let mut file = File::open(fname).unwrap();
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.seek(SeekFrom::Current(0)).unwrap();
    let mut out = (0..end).step_by((end/n_threads) as usize).collect::<Vec<u64>>();
    out.push(end);
    // println!("{:?}", end);
    // println!("{:?}", out);
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
    }
    out.dedup();
    // println!("{:?}", out);
    return out
}

