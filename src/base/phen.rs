use nalgebra::DMatrix;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, Error, ErrorKind};
use std::{str, vec};

use crate::base::*;

impl Parse<Phen> for FilePhen {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<Phen>> {
        let mut pool_names: Vec<String> = vec![];
        let mut pool_sizes: Vec<f64> = vec![];
        if self.format == "default".to_string() {
            ////////////////////
            // Default format //
            ////////////////////
            let filename_phen = self.filename.clone();
            let delim = self.delim.clone();
            let names_column_id = self.names_column_id.clone();
            let sizes_column_id = self.sizes_column_id.clone();
            let trait_values_column_ids = self.trait_values_column_ids.clone();
            let k = trait_values_column_ids.len();
            let mut phen_vec: Vec<f64> = vec![];
            let file = File::open(filename_phen.clone())
                .expect(&("Input phenotype file not found: ".to_owned() + &filename_phen[..]));
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
                    continue;
                }
                // Parse the sync line
                let vec_line = line
                    .split(&delim[..])
                    .collect::<Vec<&str>>()
                    .into_iter()
                    .map(|x| x.trim()) // remove pesky leading and trailing whitespace
                    .map(|x| x.to_owned())
                    .collect::<Vec<String>>();
                pool_names.push(vec_line[names_column_id].clone());
                pool_sizes.push(
                    vec_line[sizes_column_id]
                        .parse::<f64>()
                        .expect("T_T Pool sizes column is not a valid number."),
                );
                for j in 0..k {
                    phen_vec.push(vec_line[trait_values_column_ids[j]].parse::<f64>()
                                                            .expect("T_T Error parsing the phenotype file. The trait values specified cannot be casted into float64."));
                }
            }
            // Make sure that the pool sizes sum up to 1.00
            let s: f64 = pool_sizes.iter().sum();
            let pool_sizes = pool_sizes.into_iter().map(|x| x / s).collect::<Vec<f64>>();
            // Reshape the vector of phenotypes into an nxk matrix
            let n = phen_vec.len() / k;
            let mut phen_matrix = DMatrix::from_element(n, k, 0.0);
            for i in 0..n {
                for j in 0..k {
                    phen_matrix[(i, j)] = phen_vec[(i * k) + j]
                }
            }
            return Ok(Box::new(Phen {
                pool_names: pool_names,
                pool_sizes: pool_sizes,
                phen_matrix: phen_matrix,
            }));
        } else if self.format == "gwalpha_fmt".to_string() {
            ////////////////////
            // GWAlpha format //
            ////////////////////
            let filename_phen = self.filename.clone();
            let file = File::open(filename_phen.clone())
                .expect(&("Input phenotype file not found: ".to_owned() + &filename_phen[..]));
            let reader = BufReader::new(file);
            let mut all_lines: Vec<String> = vec![];
            for line in reader.lines() {
                all_lines.push(line.expect("T_T Phenotype file in GWAlpha format is missing some lines, e.g. Pheno_name, sig, MIN, MAX, perc and/or q."));
            }
            let _name = all_lines[0].split("=").collect::<Vec<&str>>()[1]
                .replace(";", "")
                .trim()
                .to_string();
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
            let mut _perc0 = perc.clone();
            _perc0.push(1.0);
            let perc0: DMatrix<f64> = DMatrix::from_vec(_perc0.len(), 1, _perc0);
            let mut _perc1 = vec![0.0];
            _perc1.append(&mut perc.clone());
            let perc1: DMatrix<f64> = DMatrix::from_vec(_perc1.len(), 1, _perc1);
            let bins = perc0 - perc1;
            let mut n = bins.len();
            if n < 3 {
                n = 3;
            }
            // For LS
            let mut q_prime: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
            for i in 0..q.len() {
                q_prime[i + 1] = (q[i] - min) / (max - min);
            }
            let mut phen_matrix = DMatrix::from_element(n, 3, f64::NEG_INFINITY);
            for i in 0..bins.nrows() {
                phen_matrix[(i, 0)] = bins[i];
                phen_matrix[(i, 1)] = q_prime[i];
            }
            phen_matrix[(0, 2)] = sig;
            phen_matrix[(1, 2)] = min;
            phen_matrix[(2, 2)] = max;
            return Ok(Box::new(Phen {
                pool_names: pool_names,
                pool_sizes: bins.iter().copied().collect::<Vec<f64>>(), // pool sizes, i.e. bins in the GWAlpha.py format is found in phen_matrix
                phen_matrix: phen_matrix,
            }));
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid phenotype format. PLease select: 'default' or 'gwalpha_fmt'",
            ));
        }
    }
}

impl Parse<FileSyncPhen> for (FileSync, FilePhen) {
    // Parse a line of pileup into PileupLine struct
    fn lparse(&self) -> io::Result<Box<FileSyncPhen>> {
        let filename_sync = self.0.filename.clone();
        let test = self.0.test.clone();
        if self.1.format == "default".to_string() {
            ////////////////////
            // Default format //
            ////////////////////
            let phen = self.1.lparse().unwrap();
            return Ok(Box::new(FileSyncPhen {
                filename_sync: filename_sync,
                pool_names: phen.pool_names,
                pool_sizes: phen.pool_sizes,
                phen_matrix: phen.phen_matrix,
                test: test,
            }));
        } else if self.1.format == "gwalpha_fmt".to_string() {
            ////////////////////
            // GWAlpha format //
            ////////////////////
            let phen = self.1.lparse().unwrap();
            return Ok(Box::new(FileSyncPhen {
                filename_sync: filename_sync,
                pool_names: phen.pool_names,
                pool_sizes: vec![f64::NAN], // pool sizes, i.e. bins in the GWAlpha.py format is found in phen_matrix
                phen_matrix: phen.phen_matrix,
                test: test,
            }));
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Invalid phenotype format. PLease select: 'default' or 'gwalpha_fmt'",
            ));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_phen() {
        let expected_output = FileSyncPhen {
            filename_sync: "./tests/test.sync".to_owned(),
            pool_names: vec!["G1", "G2", "G3", "G4", "G5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
            phen_matrix: DMatrix::from_column_slice(
                5,
                2,
                &[0.1, 0.3, 0.5, 0.7, 0.9, 83.2, 75.3, 49.8, 23.9, 12.0],
            ),
            test: "".to_owned(),
        };
        // Inputs
        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "".to_owned(),
        };
        // Output
        let output = *(file_sync, file_phen).lparse().unwrap();
        // Assertion
        assert_eq!(expected_output, output);
    }
}
