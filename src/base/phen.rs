use std::io::{self, prelude::*, Error, ErrorKind, BufReader};
use std::fs::File;
use std::{str, vec};
use nalgebra::DMatrix;

use crate::base::*;

impl Parse<FileSyncPhen> for (FileSync, FilePhen) {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<FileSyncPhen>> {
        let filename_sync = self.0.filename.clone();
        let test = self.0.test.clone();
        let filename_phen = self.1.filename.clone();
        let phen_delim = self.1.phen_delim.clone();
        let phen_name_col = self.1.phen_name_col.clone();
        let phen_value_col = self.1.phen_value_col.clone();
        
        let k = phen_value_col.len();
        let mut pool_names: Vec<String> = vec![];
        let mut phen_vec: Vec<f64> = vec![];
        
        let file = File::open(filename_phen).unwrap();
        let reader = BufReader::new(file);
        for l in reader.lines() {
            let mut line = l.unwrap();
            // Remove trailing newline character in Unix-like (\n) and Windows (\r)
            if line.ends_with('\n') {
                line.pop();
                if line.ends_with('\r') {
                    line.pop();
                }
            }
            // Ignore commented-out lines (i.e. '#' => 35)
            if line.as_bytes()[0] == 35 as u8 {
                continue
            }
            // Parse the sync line
            let vec_line = line.split(&phen_delim[..])
                                            .collect::<Vec<&str>>()
                                            .into_iter()
                                            .map(|x| x.to_owned())
                                            .collect::<Vec<String>>();
            pool_names.push(vec_line[phen_name_col].clone());
            for j in 0..k {
                phen_vec.push(vec_line[phen_value_col[j]].parse::<f64>()
                                                         .expect("T_T Error parsing the phenotype file. The trait values specified cannot be casted into float64."));
            }
        }
        // Reshape the vector of phenotypes into an nxk matrix
        let n = phen_vec.len() / k;
        let mut phen_matrix: DMatrix<f64> = DMatrix::from_element(n, k, 0.0);
        for i in 0..n {
            for j in 0..k {
                phen_matrix[(i, j)] = phen_vec[(i*k) + j]
            }
        }
        Ok(Box::new(FileSyncPhen{ filename_sync: filename_sync,
                                  pool_names: pool_names,
                                  phen_matrix: phen_matrix,
                                  test: test }))
        
    }   
}