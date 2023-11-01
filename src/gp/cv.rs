use crate::base::*;
use crate::gwas::*;

use ndarray::{prelude::*, stack};
use std::fs::OpenOptions;
use std::io::{self, prelude::*, Error, ErrorKind};
use std::ops::Sub;
use std::time::{SystemTime, UNIX_EPOCH};

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
            let _min = y_true
                .column(j)
                .iter()
                .fold(y_true[(0, j)], |min, &x| if x < min { x } else { min });
            let _max = y_true
                .column(j)
                .iter()
                .fold(y_true[(0, j)], |max, &x| if x > max { x } else { max });
            let (cor_, _pval) = pearsons_correlation(
                &y_true.column(j),
                &y_pred.column(j),
                &"sensible_corr".to_owned(),
            )
            .unwrap();
            cor[j] = cor_;
            mbe[j] = y_true.column(j).sub(&y_pred.column(j)).mean().unwrap();
            mae[j] = y_true
                .column(j)
                .sub(&y_pred.column(j))
                .iter()
                .map(|&x| x.abs())
                .fold(0.0, |sum, x| sum + x);
            mse[j] = y_true
                .column(j)
                .sub(&y_pred.column(j))
                .iter()
                .map(|&x| x.powf(2.0))
                .fold(0.0, |sum, x| sum + x);
            rmse[j] = mse[j].sqrt();
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
        // Check struct
        self.check().unwrap();
        // Prepare metric arrays
        let n = self.intercept_and_allele_frequencies.nrows();
        let p = self.intercept_and_allele_frequencies.ncols();
        let m = self.phenotypes.ncols();
        let l = functions.len();
        let mut models: Vec<String> = vec![];
        let mut cor: Array4<f64> = Array4::from_elem((r, k, l, m), f64::NAN);
        let mut mbe: Array4<f64> = Array4::from_elem((r, k, l, m), f64::NAN);
        let mut mae: Array4<f64> = Array4::from_elem((r, k, l, m), f64::NAN);
        let mut mse: Array4<f64> = Array4::from_elem((r, k, l, m), f64::NAN);
        let mut rmse: Array4<f64> = Array4::from_elem((r, k, l, m), f64::NAN);

        let mut y_validation_and_predicted: Array4<f64> =
            Array4::from_elem((r, l, n, m + m), f64::NAN);

        let idx_cols: Vec<usize> = (0..p).collect::<Vec<usize>>();
        for rep in 0..r {
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
                // println!("idx_validation={:?}", idx_validation);
                // println!("idx_training={:?}", idx_training);

                let y_validation: Array2<f64> = stack(
                    Axis(0),
                    idx_validation
                        .iter()
                        .map(|&x| self.phenotypes.slice(s![x, ..]))
                        .collect::<Vec<_>>()
                        .as_slice(),
                )
                .unwrap();

                // for f in functions.iter() {
                for model in 0..l {
                    let (b_hat, model_name) = functions[model](
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
                    // Save model name for the first rep and first fold only for brevity
                    if (rep == 0) & (fold == 0) {
                        models.push(model_name.clone());
                    }
                    // Save expected and predicted phenotypes (reps x n x models x traits+traits, i.e. prediced phenotypes from axis(4) field 0 to m-1 and expected phenotypes afterwards)
                    for i_ in 0..idx_validation.len() {
                        for j_ in 0..(m + m) {
                            if j_ >= m {
                                y_validation_and_predicted[(rep, model, idx_validation[i_], j_)] =
                                    y_validation[(i_, (j_ - m))];
                            } else {
                                y_validation_and_predicted[(rep, model, idx_validation[i_], j_)] =
                                    y_pred[(i_, j_)];
                            }
                        }
                    }
                    // Extract prediction performance metrics
                    let metrics = self.performance(&y_validation, &y_pred).unwrap();
                    for phe in 0..m {
                        cor[(rep, fold, model, phe)] = metrics[0][phe];
                        mbe[(rep, fold, model, phe)] = metrics[1][phe];
                        mae[(rep, fold, model, phe)] = metrics[2][phe];
                        mse[(rep, fold, model, phe)] = metrics[3][phe];
                        rmse[(rep, fold, model, phe)] = metrics[4][phe];
                    }
                }
            }
        }
        Ok(PredictionPerformance {
            n: n,
            p: p,
            k: k,
            r: r,
            models: models,
            y_validation_and_predicted: y_validation_and_predicted,
            cor: cor,
            mbe: mbe,
            mae: mae,
            mse: mse,
            rmse: rmse,
        })
    }

    fn tabulate_predict_and_output(
        &self,
        prediction_performance: &PredictionPerformance,
        functions: Vec<
            fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
        >,
        fname_input: &String,
        fname_output: &String,
    ) -> io::Result<(String, String, Vec<String>)>
    where
        fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>:
            Fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
    {
        // Prepare output basename
        let mut fname_output = fname_output.to_owned();
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        // Write tabulated performance output
        if fname_output == "".to_owned() {
            let bname = fname_input
                .split(".")
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|a| a.to_owned())
                .collect::<Vec<String>>()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()[1..]
                .to_owned()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()
                .join(".");
            fname_output = bname.to_owned() + "-cross_validation-" + &time.to_string() + ".csv";
        }
        // Instatiate output file
        let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&fname_output)
            .expect(&error_writing_file);
        file_out
            .write_all(("#rep,fold,model,phenotype,pearsons_correlation,mean_bias_error,mean_absolute_error,mean_square_error,root_mean_square_error\n").as_bytes())
            .unwrap();
        // Tabulate and write-out
        let (r, k, l, m) = prediction_performance.cor.dim();
        for rep in 0..r {
            for fold in 0..k {
                for idx_model in 0..l {
                    for phen in 0..m {
                        let line = vec![
                            rep.to_string(),
                            fold.to_string(),
                            prediction_performance.models[idx_model].clone(),
                            phen.to_string(),
                            prediction_performance.cor[(rep, fold, idx_model, phen)].to_string(),
                            prediction_performance.mbe[(rep, fold, idx_model, phen)].to_string(),
                            prediction_performance.mae[(rep, fold, idx_model, phen)].to_string(),
                            prediction_performance.mse[(rep, fold, idx_model, phen)].to_string(),
                            prediction_performance.rmse[(rep, fold, idx_model, phen)].to_string(),
                        ]
                        .join(",")
                            + "\n";
                        file_out.write_all(line.as_bytes()).unwrap();
                    }
                }
            }
        }
        // Write expected and predicted phenotypes across replications, models, pools, and traits+traits (i.e. prediced phenotypes from axis(4) field 0 to m-1 and expected phenotypes afterwards)
        let (r, l, n, m_twice) = prediction_performance.y_validation_and_predicted.dim();
        let m = m_twice / 2; // total number of traits
        let bname = fname_output
            .split(".")
            .collect::<Vec<&str>>()
            .into_iter()
            .map(|a| a.to_owned())
            .collect::<Vec<String>>()
            .into_iter()
            .rev()
            .collect::<Vec<String>>()[1..]
            .to_owned()
            .into_iter()
            .rev()
            .collect::<Vec<String>>()
            .join(".");
        let predicted_and_expected_fname =
            bname.to_owned() + "-expected_and_predicted_phenotypes.csv";
        let error_writing_file =
            "Unable to create file: ".to_owned() + &predicted_and_expected_fname;
        let mut file_out = OpenOptions::new()
            .create_new(true)
            .write(true)
            .append(false)
            .open(&predicted_and_expected_fname)
            .expect(&error_writing_file);
        let header = vec![
            "#rep,model,pool".to_owned(),
            (0..m)
                .map(|x| "predicted_trait_".to_owned() + &x.to_string()[..])
                .collect::<Vec<String>>()
                .join(","),
            (0..m)
                .map(|x| "expected_trait_".to_owned() + &x.to_string()[..])
                .collect::<Vec<String>>()
                .join(","),
        ]
        .join(",")
            + "\n";
        file_out.write_all(header.as_bytes()).unwrap();
        for rep in 0..r {
            for idx_model in 0..l {
                for pool in 0..n {
                    let line = vec![
                        rep.to_string(),
                        prediction_performance.models[idx_model].clone(),
                        self.pool_names[pool].clone(),
                        prediction_performance
                            .y_validation_and_predicted
                            .slice(s![rep, idx_model, pool, ..])
                            .iter()
                            .map(|&x| x.to_string())
                            .collect::<Vec<String>>()
                            .join(","),
                    ]
                    .join(",")
                        + "\n";
                    file_out.write_all(line.as_bytes()).unwrap();
                }
            }
        }
        // Generate the predictors for all the models tested and write-out
        let (n, p) = self.intercept_and_allele_frequencies.dim();
        let idx_all = (0..n).collect::<Vec<usize>>();
        let mut model_fit_fnames: Vec<String> = vec![];
        for f in functions.iter() {
            // Fit
            let (b_hat, model_name) = f(
                &self.intercept_and_allele_frequencies,
                &self.phenotypes,
                &idx_all,
            )
            .unwrap();
            let bname = fname_output
                .split(".")
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|a| a.to_owned())
                .collect::<Vec<String>>()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()[1..]
                .to_owned()
                .into_iter()
                .rev()
                .collect::<Vec<String>>()
                .join(".");
            let model_fit_fname = bname.to_owned() + "-genomic_predictors-" + &model_name + ".csv";
            model_fit_fnames.push(model_fit_fname.clone());
            let error_writing_file = "Unable to create file: ".to_owned() + &model_fit_fname;
            let mut file_out = OpenOptions::new()
                .create_new(true)
                .write(true)
                .append(false)
                .open(&model_fit_fname)
                .expect(&error_writing_file);
            file_out
                .write_all(("#chromosome,position,allele,phenotype,predictor\n").as_bytes())
                .unwrap();
            for i in 0..p {
                for j in 0..m {
                    let line = vec![
                        self.chromosome[i].clone(),
                        self.position[i].to_string(),
                        self.allele[i].clone(),
                        j.to_string(),
                        b_hat[(i, j)].to_string(),
                    ]
                    .join(",")
                        + "\n";
                    file_out.write_all(line.as_bytes()).unwrap();
                }
            }
        }
        Ok((fname_output, predicted_and_expected_fname, model_fit_fnames))
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
        let _expected_output1: Array2<f64> =
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
            max_missingness_rate: 0.0,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let _freqs = file_sync_phen
            .load(&filter_stats, true, &n_threads)
            .unwrap();

        // Outputs
        let n = 100;
        // let p = 50_000;
        // let q = 50;
        let p = 1_000;
        let q = 2;
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
        let mut b: Array2<f64> = Array2::zeros((p + 1, 2));
        let idx_b: Vec<usize> = dist_unif
            .sample_iter(&mut rng)
            .take(q)
            .map(|x| (x * p as f64).floor() as usize)
            .collect::<Vec<usize>>();
        for i in idx_b.into_iter() {
            b[(i, 0)] = 1.00;
        }
        // Some genome-wide polygenic trait. We are actually getting better prediction with lasso that OLS with this trait given the less complex trait, i.e. q=2 trait above. Probably because of some covariation between traits? It's helping with the prediction accuracy of the other? Or coding errors?!
        for i in 0..p {
            b[(i, 1)] = 1.0;
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
            (n, 2),
            dist_gaus
                .sample_iter(&mut rng)
                .take(2 * n)
                .collect::<Vec<f64>>(),
        )
        .unwrap();
        let y = &xb + e;
        let frequencies_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: std::iter::repeat("dummy_chr".to_owned())
                .take(p + 1)
                .collect(),
            position: (0..p + 1).map(|x| x as u64).collect(),
            allele: std::iter::repeat("A".to_owned()).take(p + 1).collect(),
            intercept_and_allele_frequencies: f.clone(),
            phenotypes: y,
            pool_names: (0..n)
                .map(|x| "pool-".to_owned() + &x.to_string()[..])
                .collect(),
            coverages: Array2::from_elem((n, p), 100.0),
        };
        let (_a, _k, _s) = frequencies_and_phenotypes.k_split(10).unwrap();
        let models: Vec<
            fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
        > = vec![
            ols,
            penalise_lasso_like,
            penalise_ridge_like,
            // penalise_lasso_like_with_iterative_proxy_norms,
            // penalise_ridge_like_with_iterative_proxy_norms,
            penalise_glmnet,
        ];
        let _m = models.len();
        let prediction_performance = frequencies_and_phenotypes
            .cross_validate(k, r, models.clone())
            .unwrap();
        let mean_cor = prediction_performance
            .cor
            .mean_axis(Axis(0))
            .unwrap()
            .mean_axis(Axis(0))
            .unwrap();
        let mean_rmse = prediction_performance
            .rmse
            .mean_axis(Axis(0))
            .unwrap()
            .mean_axis(Axis(0))
            .unwrap();
        println!("cor.column_mean()={:?}", mean_cor);
        println!("rmse.column_mean()={:?}", mean_rmse);
        println!("n={:?}", f.nrows());
        println!("p={:?}", f.ncols());
        println!("q={:?}", q);
        // Assertions
        // assert_eq!(0, 1);
        assert_eq!(mean_cor[(3, 0)].round(), 1.0);
        assert_eq!(mean_cor[(3, 1)].round(), 1.0);

        let (tabulated, pred_v_expe, predictor_files) = frequencies_and_phenotypes
            .tabulate_predict_and_output(
                &prediction_performance,
                models,
                &"test-cv-input.tmp".to_owned(),
                &"".to_owned(),
            )
            .unwrap();
        println!("tabulated={}", tabulated);
        println!("pred_v_expe={}", pred_v_expe);
        println!("predictor_files={:?}", predictor_files);
    }
}
