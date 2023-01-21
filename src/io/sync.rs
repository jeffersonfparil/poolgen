use std;
use std::fs::File;
use std::io::{self, prelude::*, SeekFrom, BufReader};
use std::time::{SystemTime, UNIX_EPOCH};

pub fn load(fname: &String) -> io::Result<f64> {
    let time = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs_f64();
    let bname = fname.split(".").into_iter().map(|a| a.to_owned()).collect::<Vec<String>>()
                     .into_iter().rev().collect::<Vec<String>>()[1..].to_owned().into_iter().rev().collect::<Vec<String>>()
                     .join(".");
    Ok(time)
}