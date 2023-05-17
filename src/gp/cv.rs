use crate::base::*;
use crate::gwas::*;
use ndarray::{prelude::*, stack};
use std::io::{self, Error, ErrorKind};
use std::ops::Sub;

// use ndarray::{prelude::*, ArrayView, ArrayView2, arr2};

impl
    CrossValidation<
        fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
    > for GenotypesAndPhenotypes
{
    fn k_split(&self, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
        let n = self.intercept_and_allele_frequencies.nrows();
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
        y_true: &Array2<f64>,
        y_pred: &Array2<f64>,
    ) -> io::Result<Vec<Array1<f64>>> {
        let n = y_true.nrows();
        let m = y_true.ncols();
        let n_ = y_pred.nrows();
        let m_ = y_pred.ncols();
        if (n != n_) | (m != m_) {
            return Err(Error::new(
                ErrorKind::Other,
                "Observed and predicted phenotypes do not match.",
            ));
        }
        let mut cor: Array1<f64> = Array1::from_elem(m, f64::NAN);
        let mut mbe: Array1<f64> = Array1::from_elem(m, f64::NAN);
        let mut mae: Array1<f64> = Array1::from_elem(m, f64::NAN);
        let mut mse: Array1<f64> = Array1::from_elem(m, f64::NAN);
        let mut rmse: Array1<f64> = Array1::from_elem(m, f64::NAN);
        for j in 0..m {
            let min = y_true
                .column(j)
                .iter()
                .fold(y_true[(0, j)], |min, &x| if x < min { x } else { min });
            let max = y_true
                .column(j)
                .iter()
                .fold(y_true[(0, j)], |max, &x| if x > max { x } else { max });
            let (cor_, _pval) = pearsons_correlation(&y_true.column(j), &y_pred.column(j)).unwrap();
            cor[j] = cor_;
            mbe[j] = y_true.column(j).sub(&y_pred.column(j)).mean().unwrap() / (max - min);
            mae[j] = y_true
                .column(j)
                .sub(&y_pred.column(j))
                .iter()
                .map(|&x| x.abs())
                .fold(0.0, |sum, x| sum + x)
                / (max - min);
            mse[j] = y_true
                .column(j)
                .sub(&y_pred.column(j))
                .iter()
                .map(|&x| x.powf(2.0))
                .fold(0.0, |sum, x| sum + x)
                / (max - min).powf(2.0);
            rmse[j] = mse[j].sqrt() / (max - min);
        }
        Ok(vec![cor, mbe, mae, mse, rmse])
    }

    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        functions: Vec<
            fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
        >,
    ) -> io::Result<PredictionPerformance>
    where
        fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>:
            Fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
    {
        let n = self.intercept_and_allele_frequencies.nrows();
        let p = self.intercept_and_allele_frequencies.ncols();
        let m = self.phenotypes.ncols();
        // let mut x_matrix_training = Array2::from_element(n - s, p, f64::NAN);
        // let mut x_matrix_validation = Array2::from_element(s, p, f64::NAN);
        // let mut y_matrix_training = Array2::from_element(n - s, 1, f64::NAN);
        // let mut y_matrix_validation = Array2::from_element(s, 1, f64::NAN);
        let l = functions.len();
        let mut models: Vec<String> = vec![];
        let mut b_histogram = vec![];
        let mut cor: Array2<f64> = Array2::from_elem((r * k * l, m), f64::NAN);
        let mut mbe: Array2<f64> = Array2::from_elem((r * k * l, m), f64::NAN);
        let mut mae: Array2<f64> = Array2::from_elem((r * k * l, m), f64::NAN);
        let mut mse: Array2<f64> = Array2::from_elem((r * k * l, m), f64::NAN);
        let mut rmse: Array2<f64> = Array2::from_elem((r * k * l, m), f64::NAN);
        let idx_cols: Vec<usize> = (0..p).collect::<Vec<usize>>();
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

                let y_validation: Array2<f64> = stack(
                    Axis(0),
                    idx_validation
                        .iter()
                        .map(|&x| self.phenotypes.slice(s![x, ..]))
                        .collect::<Vec<_>>()
                        .as_slice(),
                )
                .unwrap();

                for f in functions.iter() {
                    let (b_hat, model_name) = f(
                        &self.intercept_and_allele_frequencies,
                        &self.phenotypes,
                        &idx_training,
                    )
                    .unwrap();
                    // println!("b_hat={:?}", b_hat);
                    let y_pred: Array2<f64> = multiply_views_xx(
                        &self.intercept_and_allele_frequencies,
                        &b_hat,
                        &idx_validation,
                        &idx_cols,
                        &idx_cols,
                        &(0..b_hat.ncols()).collect::<Vec<usize>>(),
                    )
                    .unwrap();
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
        let expected_output1: Array2<f64> =
            Array2::from_shape_vec((3, 1), vec![-0.73, 5.53, 6.42]).unwrap();
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
        let _freqs = file_sync_phen
            .load(&filter_stats, true, &n_threads)
            .unwrap();

        // Outputs
        let n = 100;
        // let p = 50_000;
        // let q = 10;
        let p = 1_000;
        let q = 3;
        let h2 = 0.75;
        let (k, r) = (10, 5);
        let mut rng = rand::thread_rng();
        let dist_unif = statrs::distribution::Uniform::new(0.0, 1.0).unwrap();
        // Simulate allele frequencies
        let mut f: Array2<f64> = Array2::ones((n, p + 1));
        for i in 0..n {
            for j in 1..(p + 1) {
                f[(i, j)] = dist_unif.sample(&mut rng);
            }
        }
        // Simulate effects
        let mut b: Array2<f64> = Array2::zeros((p + 1, 1));
        let idx_b: Vec<usize> = dist_unif
            .sample_iter(&mut rng)
            .take(q)
            .map(|x| (x * p as f64).floor() as usize)
            .collect::<Vec<usize>>();
        for i in idx_b.into_iter() {
            b[(i, 0)] = 1.00;
        }
        // Simulate phenotype
        let xb = multiply_views_xx(
            &f,
            &b,
            &(0..n).collect::<Vec<usize>>(),
            &(0..(p + 1)).collect::<Vec<usize>>(),
            &(0..(p + 1)).collect::<Vec<usize>>(),
            &vec![0 as usize],
        )
        .unwrap();
        let vg = xb.var_axis(Axis(0), 0.0)[0];
        let ve = (vg / h2) - vg;
        let dist_gaus = statrs::distribution::Normal::new(0.0, ve.sqrt()).unwrap();
        let e: Array2<f64> = Array2::from_shape_vec(
            (n, 1),
            dist_gaus
                .sample_iter(&mut rng)
                .take(n)
                .collect::<Vec<f64>>(),
        )
        .unwrap();
        let y = &xb + e;
        let frequencies_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: vec!["".to_owned()],
            position: vec![0],
            allele: vec!["".to_owned()],
            intercept_and_allele_frequencies: f.clone(),
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
        let (_a, _k, _s) = frequencies_and_phenotypes.k_split(10).unwrap();
        let models: Vec<
            fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
        > = vec![
            ols,
            penalise_lasso_like,
            penalise_ridge_like,
            penalise_lasso_like_with_iterative_proxy_norms,
            penalise_ridge_like_with_iterative_proxy_norms,
        ];
        let m = models.len();
        let prediction_performance = frequencies_and_phenotypes
            .cross_validate(k, r, models)
            .unwrap();
        // println!("prediction_performance={:?}", prediction_performance);
        let cor = prediction_performance.cor.into_shape((k * r, m)).unwrap();
        let rmse = prediction_performance.rmse.into_shape((k * r, m)).unwrap();
        for i in 0..prediction_performance.models.len() {
            println!(
                "MODEL: {:?}; B: {:?}",
                prediction_performance.models[i], prediction_performance.b_histogram[i]
            );
        }
        println!("cor={:?}", cor);
        println!("rmse={:?}", rmse);
        println!("cor.column_mean()={:?}", cor.mean_axis(Axis(0)));
        println!("rmse.column_mean()={:?}", rmse.mean_axis(Axis(0)));
        println!("n={:?}", f.nrows());
        println!("p={:?}", f.ncols());
        println!("q={:?}", q);
        // Assertions
        assert_eq!(0, 0); // Output dimensions
    }
}

// // In R testing the relationships between RMSE, etc...
// mbe = c()
// mae = c()
// mse = c()
// rmse = c()
// n = 1000
// for (i in 1:n) {
//     x = runif(0,1,n=10000)
//     y = runif(0,1,n=10000)
//     mbe = c(mbe, mean(x-y))
//     mae = c(mae, mean(abs(x-y)))
//     mse = c(mse, mean((x-y)^2))
//     rmse = c(rmse, sqrt(mean((x-y)^2)))
// }
// df = data.frame(mbe, mae, mse, rmse)
// plot(df)
