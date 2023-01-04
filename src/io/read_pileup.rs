mod read_pileup;

use std::fs::File;
use std::io::{self, prelude::*, BufReader};

pub fn read() {
    fname = "/data-weedomics-2/poolgen/tests/test.pileup";
    let file: File = File::open(fname)?;
    let reader:BufReader<File> = BufReader::new(file);
    let x = "dfgfg";

}