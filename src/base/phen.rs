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
        let mut pool_names: Vec<String> = vec![];
        if self.1.format == "default".to_string() {
            ////////////////////
            // Default format //
            ////////////////////
            let filename_phen = self.1.filename.clone();
            let phen_delim = self.1.phen_delim.clone();
            let phen_name_col = self.1.phen_name_col.clone();
            let phen_value_col = self.1.phen_value_col.clone();
            let k = phen_value_col.len();
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
            let mut phen_matrix = DMatrix::from_element(n, k, 0.0);
            for i in 0..n {
                for j in 0..k {
                    phen_matrix[(i, j)] = phen_vec[(i*k) + j]
                }
            }
            return Ok(Box::new(FileSyncPhen{ filename_sync: filename_sync,
                pool_names: pool_names,
                phen_matrix: phen_matrix,
                test: test }))
        } else if self.1.format == "gwalpha_fmt".to_string() {
            ////////////////////
            // GWAlpha format //
            ////////////////////
            let filename_phen = self.1.filename.clone();
            let file = File::open(filename_phen).unwrap();
            let reader = BufReader::new(file);
            let mut all_lines: Vec<String> = vec![];
            for line in reader.lines() {
                all_lines.push(line.expect("T_T Phenotype file in GWAlpha format is missing some lines, e.g. Pheno_name, sig, MIN, MAX, perc and/or q."));
            }
            let name = all_lines[0].split("=").collect::<Vec<&str>>()[1].replace(";", "").trim().to_string();
            let sig = all_lines[1].split("=").collect::<Vec<&str>>()[1].replace(";", "").trim().parse::<f64>().expect("T_T Error parsing the standard deviation of the trait as f64 in the GWAlpha formatted phenotype file.");
            let min = all_lines[2].split("=").collect::<Vec<&str>>()[1].replace(";", "").trim().parse::<f64>().expect("T_T Error parsing the minimum value of the trait as f64 in the GWAlpha formatted phenotype file.");
            let max = all_lines[3].split("=").collect::<Vec<&str>>()[1].replace(";", "").trim().parse::<f64>().expect("T_T Error parsing the maximum value of the trait as f64 in the GWAlpha formatted phenotype file.");
            let perc = all_lines[4].split("=").collect::<Vec<&str>>()[1].replace(";", "").replace("[", "").replace("]", "").trim().to_string()
                                             .split(",").collect::<Vec<&str>>().into_iter().map(|x| x.trim().to_string())
                                             .map(|x| x.parse::<f64>().expect("T_T Error parsing the pool percentiles as f64 in the GWAlpha formatted phenotype file."))
                                             .collect::<Vec<f64>>();
            let q = all_lines[5].split("=").collect::<Vec<&str>>()[1].replace(";", "").replace("[", "").replace("]", "").trim().to_string()
                                          .split(",").collect::<Vec<&str>>().into_iter().map(|x| x.trim().to_string())
                                          .map(|x| x.parse::<f64>().expect("T_T Error parsing the pool quantiles as f64 in the GWAlpha formatted phenotype file."))
                                          .collect::<Vec<f64>>();
            // For ML
            let mut _perc0 = perc.clone(); _perc0.push(1.0);
            let perc0: DMatrix<f64> = DMatrix::from_vec(_perc0.len(), 1, _perc0);
            let mut _perc1 = vec![0.0]; _perc1.append(&mut perc.clone());
            let perc1: DMatrix<f64> = DMatrix::from_vec(_perc1.len(), 1, _perc1);
            let bins = perc0 - perc1;
            let mut n = bins.len() + 1;
            if n < 3 {
                n = 3;
            }
            // For LS
            let _q: DMatrix<f64> = DMatrix::from_vec(q.len(), 1, q.clone());
            let _q_min: DMatrix<f64> = DMatrix::from_element(q.len(), 1, min);
            let mut q_prime: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
            for i in 0..n {
                q_prime[i] = (_q[i] - _q_min[i]) / (max-min);
            }
            let mut phen_matrix = DMatrix::from_element(n, 3, f64::NEG_INFINITY);
            for i in 0..bins.nrows() {
                phen_matrix[(i, 0)] = bins[i];
                phen_matrix[(i, 1)] = q_prime[i];
            }
            phen_matrix[(0, 2)] = sig;
            phen_matrix[(1, 2)] = min;
            phen_matrix[(2, 2)] = max;
            return Ok(Box::new(FileSyncPhen{ filename_sync: filename_sync,
                                             pool_names: pool_names,
                                             phen_matrix: phen_matrix,
                                             test: test }))
        } else {
            return Err(Error::new(ErrorKind::Other, "Invalid phenotype format. PLease select: 'default' or 'gwalpha_fmt'"))
        }
    }   
}
