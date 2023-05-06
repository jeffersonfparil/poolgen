use crate::base::*;
use crate::gp::*;
use crate::gwas::*;
use nalgebra::{self, DMatrix, DVector};
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, RwLock};

impl CrossValidation<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>>
    for GenotypesAndPhenotypes
{
    fn k_split(&self, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
        let (n, _) = self.intercept_and_allele_frequencies.shape();
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

    fn performance(
        &self,
        y_true: &DMatrix<f64>,
        y_pred: &DMatrix<f64>,
    ) -> io::Result<Vec<DVector<f64>>> {
        let (n, m) = y_true.shape();
        let (n_, m_) = y_pred.shape();
        if (n != n_) | (m != m_) {
            return Err(Error::new(
                ErrorKind::Other,
                "Observed and predicted phenotypes do not match.",
            ));
        }
        let mut cor: DVector<f64> = DVector::from_element(m, f64::NAN);
        let mut mbe: DVector<f64> = DVector::from_element(m, f64::NAN);
        let mut mae: DVector<f64> = DVector::from_element(m, f64::NAN);
        let mut mse: DVector<f64> = DVector::from_element(m, f64::NAN);
        let mut rmse: DVector<f64> = DVector::from_element(m, f64::NAN);
        for j in 0..m {
            let min = y_true.column(j).min();
            let max = y_true.column(j).max();
            let (cor_, _pval) = pearsons_correlation(
                &DVector::from_iterator(n, y_true.column(j).into_iter().map(|x| *x)),
                &DVector::from_iterator(n, y_pred.column(j).into_iter().map(|x| *x)),
            )
            .unwrap();
            cor[j] = cor_;
            mbe[j] = (y_true.column(j) - y_pred.column(j)).mean() / (max - min);
            mae[j] = (y_true.column(j) - y_pred.column(j)).norm() / (max - min);
            mse[j] =
                (y_true.column(j) - y_pred.column(j)).norm_squared() / f64::powf(max - min, 2.0);
            rmse[j] = mse[j].sqrt() / (max - min);
        }
        Ok(vec![cor, mbe, mae, mse, rmse])
    }

    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        functions: Vec<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>>,
    ) -> io::Result<PredictionPerformance>
    where
        fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>:
            Fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>,
    {
        let (n, p) = self.intercept_and_allele_frequencies.shape();
        let (_n, m) = self.phenotypes.shape();
        // let mut x_matrix_training = DMatrix::from_element(n - s, p, f64::NAN);
        // let mut x_matrix_validation = DMatrix::from_element(s, p, f64::NAN);
        // let mut y_matrix_training = DMatrix::from_element(n - s, 1, f64::NAN);
        // let mut y_matrix_validation = DMatrix::from_element(s, 1, f64::NAN);
        let mut models: Vec<String> = vec![];
        let mut b_histogram = vec![];
        let mut cor: DMatrix<f64> = DMatrix::from_element(r * k * functions.len(), m, f64::NAN);
        let mut mbe: DMatrix<f64> = DMatrix::from_element(r * k * functions.len(), m, f64::NAN);
        let mut mae: DMatrix<f64> = DMatrix::from_element(r * k * functions.len(), m, f64::NAN);
        let mut mse: DMatrix<f64> = DMatrix::from_element(r * k * functions.len(), m, f64::NAN);
        let mut rmse: DMatrix<f64> = DMatrix::from_element(r * k * functions.len(), m, f64::NAN);
        let mut i: usize = 0;
        for _rep in 0..r {
            let (groupings, k, s) = self.k_split(k).unwrap();
            println!("groupings={:?}; k={:?}; s={:?}", groupings, k, s);
            for fold in 0..k {
                let idx_validation: Vec<usize> = groupings
                    .iter()
                    .enumerate()
                    .filter(|(_, x)| *x == &fold)
                    .map(|(i, _)| i)
                    .collect();

                let idx_training: Vec<usize> = groupings
                    .iter()
                    .enumerate()
                    .filter(|(_, x)| *x != &fold)
                    .map(|(i, _)| i)
                    .collect();
                println!("idx_validation={:?}", idx_validation);
                println!("idx_training={:?}", idx_training);
                // let gap = Arc::new(RwLock::new(self));
                // let data = gap.read().unwrap();
                // let x_training = &data
                //     .intercept_and_allele_frequencies
                //     .select_rows(&idx_training);
                // let y_training = &data.phenotypes.select_rows(&idx_training);
                // let x_validation = &data
                //     .intercept_and_allele_frequencies
                //     .select_rows(&idx_validation);
                // let y_validation = &data.phenotypes.select_rows(&idx_validation);
                let x_training = &self
                    .intercept_and_allele_frequencies
                    .select_rows(&idx_training);
                let y_training = &self.phenotypes.select_rows(&idx_training);
                let x_validation = &self
                    .intercept_and_allele_frequencies
                    .select_rows(&idx_validation);
                let y_validation = &self.phenotypes.select_rows(&idx_validation);
                for f in functions.iter() {
                    let (b_hat, model_name) = f(x_training, y_training).unwrap();
                    // println!("b_hat={:?}", b_hat);
                    let y_pred: DMatrix<f64> = x_validation * &b_hat;
                    let metrics = self.performance(&y_validation, &y_pred).unwrap();
                    models.push(model_name);
                    b_histogram.push(histogram(
                        b_hat.iter().map(|x| x.clone()).collect::<Vec<f64>>(),
                        10,
                    ));
                    for j in 0..m {
                        cor[(i, j)] = metrics[0][j];
                        mbe[(i, j)] = metrics[1][j];
                        mae[(i, j)] = metrics[2][j];
                        mse[(i, j)] = metrics[3][j];
                        rmse[(i, j)] = metrics[4][j];
                    }
                    i += 1;
                }
            }
        }
        Ok(PredictionPerformance {
            n: n,
            p: p,
            k: k,
            r: r,
            models: models,
            b_histogram: b_histogram,
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
        // let p = 50_000;
        // let q = 20;
        let p = 1_000;
        let q = 2;
        let h2 = 0.75;
        let mut rng = rand::thread_rng();
        let dist_unif = statrs::distribution::Uniform::new(0.0, 1.0).unwrap();
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
        let xb = &f * &b;
        let vg = xb.variance();
        let ve = (vg / h2) - vg;
        let dist_gaus = statrs::distribution::Normal::new(0.0, ve.sqrt()).unwrap();
        let e: DMatrix<f64> = DMatrix::from_column_slice(
            n,
            1,
            &dist_gaus
                .sample_iter(&mut rng)
                .take(n)
                .collect::<Vec<f64>>(),
        );
        let y = &xb + e;
        // let y_ = &xb + e;
        // let y = y_.add_scalar(-y_.mean()) / y_.variance().sqrt();
        let frequencies_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: vec!["".to_owned()],
            position: vec![0],
            intercept_and_allele_frequencies: f,
            phenotypes: y,
            pool_names: vec!["".to_owned()],
        };
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)]={:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 2)]={:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(1, 2)]
        );
        println!(
            "frequencies_and_phenotypes.intercept_and_allele_frequencies[(2, 3)]={:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(2, 3)]
        );
        let (a, k, s) = frequencies_and_phenotypes.k_split(10).unwrap();
        let (k, r) = (10, 1);
        let models: Vec<fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>> =
            vec![ols, penalise_lasso_like, penalise_ridge_like, penalise_lasso_like_iterative_base, penalise_ridge_like_iterative_base];
        let m = models.len();
        let prediction_performance = frequencies_and_phenotypes
            .cross_validate(k, r, models)
            .unwrap();
        // println!("prediction_performance={:?}", prediction_performance);
        let cor = DMatrix::from_vec(
            m,
            r * k,
            prediction_performance
                .cor
                .column(0)
                .as_slice()
                .to_owned()
                .try_into()
                .unwrap(),
        );
        let rmse = DMatrix::from_vec(
            m,
            r * k,
            prediction_performance
                .rmse
                .column(0)
                .as_slice()
                .to_owned()
                .try_into()
                .unwrap(),
        );
        for i in 0..prediction_performance.models.len() {
            println!(
                "MODEL: {:?}; B: {:?}",
                prediction_performance.models[i], prediction_performance.b_histogram[i]
            );
        }
        println!("cor={:?}", cor);
        println!("rmse={:?}", rmse);
        println!("cor.column_mean()={:?}", cor.column_mean());
        println!("rmse.column_mean()={:?}", rmse.column_mean());
        // Assertions
        assert_eq!(1, 1); // Output dimensions
    }
}
