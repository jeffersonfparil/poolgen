use crate::base::*;
use crate::gwas::*;
use nalgebra::{self, DMatrix, DVector};
use std::io::{self, Error, ErrorKind};

impl EstmateAndPredict<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>>
    for (&DMatrix<f64>, &DMatrix<f64>)
{
    fn estimate_effects(
        &self,
        function: fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>,
    ) -> io::Result<(DMatrix<f64>, String)>
    where
        fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<DMatrix<f64>>:
            Fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<DMatrix<f64>>,
    {
        // b_hat
        function(self.0, self.1)
    }

    fn predict_phenotypes(&self) -> io::Result<DMatrix<f64>> {
        // y_hat
        Ok(self.0 * self.1)
    }
}

impl CrossValidate<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>>
    for FrequenciesAndPhenotypes
{
    fn k_split(&self, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
        let (n, _) = self.intercept_and_frequencies.shape();
        if (k >= n) | (n <= 2) {
            return Err(Error::new(ErrorKind::Other, "The number of splits, i.e. k, needs to be less than the number of pools, n, and n > 2. We are aiming for fold sizes of 10 or greater."));
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
        let mut g = (0..k)
            .flat_map(|x| std::iter::repeat(x).take(s))
            .collect::<Vec<usize>>();
        if n - s > 0 {
            for _i in 0..(n - s) {
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
        let n = y_true.len();
        let min = y_true.min();
        let max = y_true.max();
        let (cor, _pval) = pearsons_correlation(
            &DVector::from_iterator(n, y_true.column(0).into_iter().map(|x| *x)),
            &DVector::from_iterator(n, y_pred.column(0).into_iter().map(|x| *x)),
        )
        .unwrap();
        let mbe = (y_true - y_pred).mean() / (max - min);
        let mae = (y_true - y_pred).norm() / (max - min);
        let mse = (y_true - y_pred).norm_squared() / f64::powf(max - min, 2.0);
        let rmse = mse.sqrt() / (max - min);
        Ok(vec![cor, mbe, mae, mse, rmse])
    }

    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        functions: Vec<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>>,
    ) -> io::Result<PredictionPerformance> {
        let (n, p) = self.intercept_and_frequencies.shape();
        let (groupings, k, s) = self.k_split(k).unwrap();
        let mut x_matrix_training = DMatrix::from_element(n - s, p, f64::NAN);
        let mut x_matrix_validation = DMatrix::from_element(s, p, f64::NAN);
        let mut y_matrix_training = DMatrix::from_element(n - s, 1, f64::NAN);
        let mut y_matrix_validation = DMatrix::from_element(s, 1, f64::NAN);
        let mut models: Vec<String> = vec![];
        let mut cor: Vec<f64> = vec![];
        let mut mbe: Vec<f64> = vec![];
        let mut mae: Vec<f64> = vec![];
        let mut mse: Vec<f64> = vec![];
        let mut rmse: Vec<f64> = vec![];
        for _rep in 0..r {
            for fold in 0..k {
                for j in 0..p {
                    let mut idx_training: usize = 0;
                    let mut idx_validation: usize = 0;
                    for i in 0..n {
                        if fold == groupings[i] {
                            x_matrix_validation[(idx_validation, j)] =
                                self.intercept_and_frequencies[(i, j)];
                            if j == 0 {
                                y_matrix_validation[(idx_validation, j)] = self.phenotypes[(i, j)];
                            }
                            idx_validation += 1;
                        } else {
                            x_matrix_training[(idx_training, j)] =
                                self.intercept_and_frequencies[(i, j)];
                            if j == 0 {
                                y_matrix_training[(idx_training, j)] = self.phenotypes[(i, j)];
                            }
                            idx_training += 1;
                        }
                    }
                }
                for f in functions.iter() {
                    let (b_hat, model_name) = (&x_matrix_training, &y_matrix_training)
                        .estimate_effects(*f)
                        .unwrap();
                    // println!("b_hat={:?}", b_hat);
                    let y_pred: DMatrix<f64> =
                        (&x_matrix_validation, &b_hat).predict_phenotypes().unwrap();
                    let metrics = self.performance(&y_matrix_validation, &y_pred).unwrap();
                    models.push(model_name);
                    cor.push(metrics[0]);
                    mbe.push(metrics[1]);
                    mae.push(metrics[2]);
                    mse.push(metrics[3]);
                    rmse.push(metrics[4]);
                }
            }
        }
        Ok(PredictionPerformance {
            n: n,
            p: p,
            k: k,
            r: r,
            models: models,
            cor: cor,
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
    use crate::gp::*;
    use rand::prelude::*;
    use statrs;
    // use statrs::statistics::Distribution;
    #[test]
    fn test_cv() {
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
        let q = 1;
        let mut rng = rand::thread_rng();
        let dist_unif = statrs::distribution::Uniform::new(0.0, 1.0).unwrap();
        let dist_gaus = statrs::distribution::Normal::new(0.0, 0.01).unwrap();
        // Simulate allele frequencies
        let mut f: DMatrix<f64> = DMatrix::from_column_slice(
            n,
            p,
            &dist_unif
                .sample_iter(&mut rng)
                .take(n * p)
                .collect::<Vec<f64>>(),
        );
        // Simulate effects
        let mut b: DMatrix<f64> = DMatrix::from_element(p, 1, 0.0);
        let idx_b: Vec<usize> = dist_unif
            .sample_iter(&mut rng)
            .take(q)
            .map(|x| (x * p as f64).floor() as usize)
            .collect::<Vec<usize>>();
        for i in idx_b.into_iter() {
            b[(i, 0)] = 1.00;
        }
        // Insert intercept
        f = f.insert_column(0, 1.0);
        b = b.insert_row(0, 0.0);
        // Simulate phenotype
        let e: DMatrix<f64> = DMatrix::from_column_slice(
            n,
            1,
            &dist_gaus
                .sample_iter(&mut rng)
                .take(n)
                .collect::<Vec<f64>>(),
        );
        let y = (&f * &b) + e;
        // let y = DMatrix::from_column_slice(n, 1, &dist_gaus.sample_iter(&mut rng).take(n).collect::<Vec<f64>>());
        let frequencies_and_phenotypes = FrequenciesAndPhenotypes {
            chromosome: vec!["".to_owned()],
            position: vec![0],
            intercept_and_frequencies: f,
            phenotypes: y,
            pool_names: vec!["".to_owned()],
        };
        println!(
            "frequencies_and_phenotypes.intercept_and_frequencies[(0, 1)]={:?}",
            frequencies_and_phenotypes.intercept_and_frequencies[(0, 1)]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_frequencies[(1, 2)]={:?}",
            frequencies_and_phenotypes.intercept_and_frequencies[(1, 2)]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_frequencies[(2, 3)]={:?}",
            frequencies_and_phenotypes.intercept_and_frequencies[(2, 3)]
        );
        let (a, k, s) = frequencies_and_phenotypes.k_split(10).unwrap();
        let (k, r) = (10, 1);
        let models: Vec<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>> =
            // vec![ols_test, ols, ols2];
            vec![ols_new, ols, penalise_lasso_like, penalise_ridge_like];
        let m = models.len();
        let prediction_performance = frequencies_and_phenotypes
            .cross_validate(k, r, models)
            .unwrap();
        // println!("prediction_performance={:?}", prediction_performance);
        let cor = DMatrix::from_row_slice(k * r, m, &prediction_performance.cor);
        let rmse = DMatrix::from_row_slice(k * r, m, &prediction_performance.rmse);
        println!(
            "prediction_performance.models={:?}",
            prediction_performance.models
        );
        println!("cor.row_mean()={:?}", cor.row_mean());
        println!("rmse.row_mean()={:?}", rmse.row_mean());
        // Assertions
        assert_eq!(0, 0); // Output dimensions
    }
}
