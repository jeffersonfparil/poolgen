use crate::base::*;
use crate::gp::*;
use nalgebra::{self, DMatrix, DVector};
use nalgebra::base::storage::StorageMut;
use std::io::{self, Error, ErrorKind};

impl CrossValidate<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<DMatrix<f64>>> for FrequenciesAndPhenotypes {
    fn k_split(&self, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
        let (n, _) = self.frequencies.shape();
        if k >= n {
            return Err(Error::new(ErrorKind::Other, "The number of splits, i.e. k needs to be less than the number of pools, n. We are aiming for fold sizes of 10 or greater."))
        }
        let mut s = (n as f64 / k as f64).floor() as usize;
        while s < 10 {
            if n < 20 {
                println!("Warning: number of pools is less than 20, so we're using k=2.");
                k = 2;
                s = (n as f64 / k as f64).floor() as usize;
                break;
            }
            k -= 1;
            s = (n as f64 / k as f64).floor() as usize;    
        }
        let mut g = (0..k).flat_map(|x| std::iter::repeat(x).take(s)).collect::<Vec<usize>>();
        if n-s > 0 {
            for i in 0..(n-s) {
                g.push(k);
            }
        }
        let mut rng = rand::thread_rng();
        let shuffle = rand::seq::index::sample(&mut rng, n, n)
            .into_iter()
            .map(|x| x as usize)
            .collect::<Vec<usize>>();
        let mut out: Vec<usize> = Vec::new();
        for i in 0..n {
            out.push(g[shuffle[i]]);
        }
        Ok((out, k, s))
    }

    fn performance(&self, y_true: &DMatrix<f64>, y_pred: &DMatrix<f64>) -> io::Result<Vec<f64>> {
        let mbe = (y_true - y_pred).mean();
        let mae = (y_true - y_pred).norm();
        let mse = (y_true - y_pred).norm_squared();
        let rmse = mse.sqrt();
        Ok(vec![mbe, mae, mse, rmse])
    }

    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        function: fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<DMatrix<f64>>,
    ) -> io::Result<PredictionPerformance> {
        let (n, p) = self.frequencies.shape();
        let (groupings, k, s) = self.k_split(k).unwrap();
        let mut x_matrix_training = DMatrix::from_element(n-s, p+1, 1.00); // including intercept in the first column
        let mut x_matrix_validation = DMatrix::from_element(s, p+1, 1.00); // including intercept in the first column
        let mut y_matrix_training = DMatrix::from_element(n-s, 1, f64::NAN);
        let mut y_matrix_validation = DMatrix::from_element(s, 1, f64::NAN);
        let mut mbe = DVector::from_element(r*k, f64::NAN);
        let mut mae = DVector::from_element(r*k, f64::NAN);
        let mut mse = DVector::from_element(r*k, f64::NAN);
        let mut rmse = DVector::from_element(r*k, f64::NAN);
        for _rep in 0..r {
            for fold in 0..k {
                for j in 0..p {
                    let mut idx_training: usize = 0;
                    let mut idx_validation: usize = 0;
                    for i in 0..n {
                        if fold == groupings[i] {
                            x_matrix_validation[(idx_validation,j+1)] = self.frequencies[(i,j)];
                            if j == 0 {
                                y_matrix_validation[(idx_validation,j)] = self.phenotypes[(i,j)];
                            }
                            idx_validation += 1;
                        } else {
                            x_matrix_training[(idx_training,j+1)] = self.frequencies[(i,j)];
                            if j == 0 {
                                y_matrix_training[(idx_training,j)] = self.phenotypes[(i,j)];
                            }
                            idx_training += 1;
                        }
                    }
                }
                let b_hat = function(&x_matrix_training, &y_matrix_training).unwrap();
                // println!("b_hat={:?}", b_hat);
                let y_pred = &x_matrix_validation * b_hat;
                let metrics = self.performance(&y_matrix_validation, &y_pred).unwrap();
                mbe[fold]  = metrics[0];
                mae[fold]  = metrics[1];
                mse[fold]  = metrics[2];
                rmse[fold]  = metrics[3];
            }
        }
        Ok(PredictionPerformance {
            n: n,
            p: p,
            k: k,
            r: r,
            model: stringify!(function).to_owned(),
            mbe: mbe,
            mae: mae,
            mse: mse,
            rmse: rmse,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::{Distribution, Uniform};
    #[test]
    fn test_ols() {
        // Expected
        let expected_output1: DMatrix<f64> = DMatrix::from_column_slice(3, 1, &[-0.73, 5.53, 6.42]);
        // Inputs
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "load".to_owned(),
        };
        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        let n_threads = 2;
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let freqs = file_sync_phen
            .load(&filter_stats, true, &n_threads)
            .unwrap();

        // Outputs
        let n = 100;
        let p = 1_000;
        let mut rng = rand::thread_rng();
        let dist_unif = Uniform::new(0.0, 1.0);
        let f: Vec<f64> = dist_unif.sample_iter(&mut rng).take(n*p).collect::<Vec<f64>>();
        let y: Vec<f64> = dist_unif.sample_iter(&mut rng).take(n*2).collect::<Vec<f64>>();
        let frequencies_and_phenotypes = FrequenciesAndPhenotypes {
            chromosome: vec!["".to_owned()],
            position: vec![0],
            frequencies: DMatrix::from_column_slice(100, 1000, &f),
            phenotypes: DMatrix::from_column_slice(100, 2, &y),
            pool_names: vec!["".to_owned()],
        };
        let (a, k, s) = frequencies_and_phenotypes.k_split(10).unwrap();
        let prediction_performance = frequencies_and_phenotypes.cross_validate(10, 1, ols).unwrap();
        println!("prediction_performance={:?}", prediction_performance);
        // Assertions
        assert_eq!(0, 0); // Output dimensions
    }
}
