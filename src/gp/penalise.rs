use crate::base::*;
use crate::gp::*;
use crate::gwas::*;
use ndarray::{prelude::*, Zip};
use rand::prelude::*;
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, Mutex};

#[function_name::named]
pub fn penalise_lasso_like(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, lambdas) = penalised_lambda_path_with_k_fold_cross_validation(
        x,
        y,
        row_idx,
        false,
        "Lasso-like".to_owned(),
        0.1,
        10,
    )
    .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("-"),
    ))
}

#[function_name::named]
pub fn penalise_ridge_like(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, lambdas) = penalised_lambda_path_with_k_fold_cross_validation(
        x,
        y,
        row_idx,
        false,
        "Ridge-like".to_owned(),
        0.1,
        10,
    )
    .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("-"),
    ))
}

fn expand_and_contract(b_hat: &Array2<f64>, norm: String, lambda: f64) -> io::Result<Array2<f64>> {
    // Clone b_hat
    let mut b_hat: Array2<f64> = b_hat.clone();
    let (p, k) = (b_hat.nrows(), b_hat.ncols());
    for j in 0..k {
        //Exclude the intercept from penalisation
        let mut intercept = b_hat[(0, j)];
        // Norm 1 or norm 2 (exclude the intercept)
        let normed: Array1<f64> = if norm == "Lasso-like".to_owned() {
            b_hat.column(j).slice(s![1..p]).map(|&x| x.abs())
            // b_hat.rows(1, p - 1).map(|x| x.abs())
        } else if norm == "Ridge-like".to_owned() {
            b_hat.column(j).slice(s![1..p]).map(|&x| x.powf(2.0))
            // b_hat.rows(1, p - 1).map(|x| x.powf(2.0))
        } else {
            return Err(Error::new(
                ErrorKind::Other,
                "Please enter: 'Lasso-like' or 'Ridge-like' norms.",
            ));
        };
        // Find estimates that will be penalised
        let normed_max = normed
            .iter()
            .fold(normed[0], |max, &x| if x > max { x } else { max });
        let normed_scaled: Array1<f64> = &normed / normed_max;
        let idx_penalised = normed_scaled
            .iter()
            .enumerate()
            .filter(|(_, &value)| value < lambda)
            .map(|(index, _)| index)
            .collect::<Vec<usize>>();
        let idx_depenalised = normed_scaled
            .iter()
            .enumerate()
            .filter(|(_, &value)| value >= lambda)
            .map(|(index, _)| index)
            .collect::<Vec<usize>>();

        // Penalise: contract
        let mut subtracted_penalised = 0.0;
        let mut added_penalised = 0.0;
        for i in idx_penalised.into_iter() {
            if b_hat[(i + 1, j)] >= 0.0 {
                b_hat[(i + 1, j)] -= normed[i];
                subtracted_penalised += normed[i];
            } else {
                b_hat[(i + 1, j)] += normed[i];
                added_penalised += normed[i];
            }
        }
        // Find total depenalised (expanded) values
        let mut subtracted_depenalised = 0.0;
        let mut added_depenalised = 0.0;
        for i in idx_depenalised.clone().into_iter() {
            if b_hat[(i + 1, j)] >= 0.0 {
                subtracted_depenalised += normed[i];
            } else {
                added_depenalised += normed[i];
            }
        }
        // Account for the absence of available slots to transfer the contracted effects into
        if (subtracted_penalised > 0.0) & (subtracted_depenalised == 0.0) {
            added_penalised -= subtracted_penalised;
            subtracted_penalised = 0.0;
        } else if (added_penalised > 0.0) & (added_depenalised == 0.0) {
            subtracted_penalised -= added_penalised;
            added_penalised = 0.0;
        }
        if (subtracted_penalised < 0.0)
            | ((subtracted_depenalised == 0.0) & (added_depenalised == 0.0))
        {
            intercept += subtracted_penalised;
            intercept -= added_penalised;
            subtracted_penalised = 0.0;
            added_penalised = 0.0;
        }
        // Depenalise: expand
        for i in idx_depenalised.into_iter() {
            if b_hat[(i + 1, j)] >= 0.0 {
                b_hat[(i + 1, j)] += subtracted_penalised * (normed[i] / subtracted_depenalised);
            } else {
                b_hat[(i + 1, j)] -= added_penalised * (normed[i] / added_depenalised);
            }
        }
        // Insert the unpenalised intercept
        b_hat[(0, j)] = intercept;
    }
    Ok(b_hat)
}

fn error_index(
    b_hat: &Array2<f64>,
    x: &Array2<f64>,
    y: &Array2<f64>,
    idx_validation: &Vec<usize>,
) -> io::Result<Vec<f64>> {
    let (n, p) = (idx_validation.len(), x.ncols());
    let k = y.ncols();
    let (p_, k_) = (b_hat.nrows(), b_hat.ncols());
    if p != p_ {
        return Err(Error::new(
            ErrorKind::Other,
            "The X matrix is incompatible with b_hat.",
        ));
    }
    if k != k_ {
        return Err(Error::new(
            ErrorKind::Other,
            "The y matrix/vector is incompatible with b_hat.",
        ));
    }
    let idx_b_hat = &(0..p).collect();
    let mut error_index: Vec<f64> = Vec::with_capacity(k);
    for j in 0..k {
        // let y_pred: Array2<f64> = x * b_hat;
        let y_true: Array2<f64> = Array2::from_shape_vec(
            (n, 1),
            idx_validation
                .iter()
                .map(|&i| y[(i, j)])
                .collect::<Vec<f64>>(),
        )
        .unwrap();
        let y_pred: Array2<f64> =
            multiply_views_xx(x, b_hat, idx_validation, idx_b_hat, idx_b_hat, &vec![j]).unwrap();
        let min = y_true
            .iter()
            .fold(y_true[(0, j)], |min, &x| if x < min { x } else { min });
        let max = y_true
            .iter()
            .fold(y_true[(0, j)], |max, &x| if x > max { x } else { max });
        let (cor, _pval) = pearsons_correlation(&y_true.column(j), &y_pred.column(j)).unwrap();
        // let mbe = (y_true - &y_pred).mean() / (max - min);vec![0.0]
        let mae = (&y_true - &y_pred)
            .iter()
            .fold(0.0, |norm, &x| norm + x.abs())
            / (max - min);
        let mse = (&y_true - &y_pred)
            .iter()
            .fold(0.0, |norm, &x| norm + x.powf(2.0))
            / (max - min).powf(2.0);
        let rmse = mse.sqrt() / (max - min);
        error_index.push(((1.0 - cor.abs()) + mae + mse + rmse) / 4.0);
    }
    Ok(error_index)
}

fn k_split(row_idx: &Vec<usize>, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
    let n = row_idx.len();
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
    let shuffle = row_idx.iter().map(|&x| x).choose_multiple(&mut rng, n);
    let mut out: Vec<usize> = Vec::new();
    for i in 0..n {
        out.push(g[shuffle[i]]);
    }
    Ok((out, k, s))
}

fn penalised_lambda_path_with_k_fold_cross_validation(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
    iterative: bool,
    norm: String,
    lambda_step_size: f64,
    r: usize,
) -> io::Result<(Array2<f64>, Vec<f64>)> {
    let (n, p) = (row_idx.len(), x.ncols());
    let k = y.ncols();
    let max_usize: usize = (1.0 / lambda_step_size).round() as usize;
    let lambda_path: Array1<f64> = (0..max_usize)
        .into_iter()
        .map(|x| (x as f64) / (max_usize as f64))
        .collect();
    let (_, nfolds, s) = k_split(row_idx, 10).unwrap();
    let mut performances: Array2<f64> =
        Array2::from_elem((r * nfolds * lambda_path.len(), k), f64::NAN);
    for rep in 0..r {
        let (groupings, _, _) = k_split(row_idx, 10).unwrap();
        for fold in 0..nfolds {
            let idx_validation: Vec<usize> = groupings
                .iter()
                .enumerate()
                .filter(|(_, x)| *x == &fold)
                .map(|(i, _)| row_idx[i])
                .collect();
            let idx_training: Vec<usize> = groupings
                .iter()
                .enumerate()
                .filter(|(_, x)| *x != &fold)
                .map(|(i, _)| row_idx[i])
                .collect();
            // let (b_hat, _) = if iterative {
            //     ols_iterative_for_penalisation(&x, &y, &idx_training).unwrap()
            // } else {
            //     ols(&x, &y, &idx_training).unwrap()
            // };
            let (b_hat, _) = ols(&x, &y, &idx_training).unwrap();
            let mut errors: Array1<Vec<f64>> = Array1::from_elem(lambda_path.len(), vec![]);
            Zip::from(&mut errors)
                .and(&lambda_path)
                .par_for_each(|err, &lambda| {
                    let b_hat_new: Array2<f64> =
                        expand_and_contract(&b_hat, norm.clone(), lambda).unwrap();
                    *err = error_index(&b_hat_new, x, y, &idx_validation).unwrap();
                });

            let start = (((rep * nfolds) + fold) * lambda_path.len()) + 0;
            let end = (((rep * nfolds) + fold) * lambda_path.len()) + lambda_path.len();
            let idx_performances = (start..end).collect::<Vec<usize>>();

            for i in 0..lambda_path.len() {
                let i_ = idx_performances[i];
                for j in 0..k {
                    performances[(i_, j)] = errors[i][j];
                }
            }
        }
    }

    // let error_indices_across_folds_x_lambdas: Array2<f64> =
    //     Array2::from_row_slice(r * k, lambda_path.len(), &performances);
    // Array2::from_row_slice(k, lambda_path.len(), &(0..(k*lambda_path.len())).into_iter().map(|x| x as f64).collect::<Vec<f64>>());
    // Find best lambda and estimate effects on the full dataset
    let (b_hat, _) = ols(x, y, row_idx).unwrap();
    let mut b_hat_penalised = b_hat.clone();
    let mut lambdas = vec![];
    for j in 0..k {
        let perf = performances
            .column(j)
            .into_shape((r * nfolds, lambda_path.len()))
            .unwrap();
        let mean_error: Array1<f64> = perf.mean_axis(Axis(0)).unwrap();
        let min_error = mean_error
            .iter()
            .fold(mean_error[0], |min, &x| if x < min { x } else { min });
        let idx = mean_error.iter().position(|&x| x == min_error).unwrap();
        lambdas.push(lambda_path[idx]);
        println!("mean_error={:?}", mean_error);
        println!("lambdas={:?}", lambdas);
        let b_hat_penalised_2d: Array2<f64> =
            expand_and_contract(&b_hat, norm.clone(), lambdas[j]).unwrap();
        for i in 0..p {
            b_hat_penalised[(i, j)] = b_hat_penalised_2d[(i, j)];
        }
    }
    Ok((b_hat_penalised, lambdas))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::concatenate;
    #[test]
    fn test_penalised() {
        let b: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![5.0, -0.4, 0.0, 1.0, -0.1, 1.0, 0.0]).unwrap();
        let new_b: Array2<f64> = expand_and_contract(&b, "Lasso-like".to_owned(), 0.5).unwrap();
        println!("new_b={:?}", new_b);
        let expected_output1: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![4.5, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0]).unwrap();
        assert_eq!(expected_output1, new_b);

        // let n = 100;
        // // let p = 50_000;
        // // let q = 20;
        // let p = 1_000;
        // let q = 2;
        // let h2 = 0.75;
        // let mut rng = rand::thread_rng();
        // let dist_unif = statrs::distribution::Uniform::new(0.0, 1.0).unwrap();
        // // Simulate allele frequencies
        // let mut f: Array2<f64> = Array2::ones((n, p + 1));
        // for i in 0..n {
        //     for j in 1..(p + 1) {
        //         f[(i, j)] = dist_unif.sample(&mut rng);
        //     }
        // }
        // // Simulate effects
        // let mut b: Array2<f64> = Array2::zeros((p + 1, 1));
        // let idx_b: Vec<usize> = dist_unif
        //     .sample_iter(&mut rng)
        //     .take(q)
        //     .map(|x| (x * p as f64).floor() as usize)
        //     .collect::<Vec<usize>>();
        // for i in idx_b.into_iter() {
        //     b[(i, 0)] = 1.00;
        // }
        // // Simulate phenotype
        // let xb = multiply_views_xx(
        //     &f,
        //     &b,
        //     &(0..n).collect::<Vec<usize>>(),
        //     &(0..(p + 1)).collect::<Vec<usize>>(),
        //     &(0..(p + 1)).collect::<Vec<usize>>(),
        //     &vec![0 as usize],
        // )
        // .unwrap();
        // let vg = xb.var_axis(Axis(0), 0.0)[0];
        // let ve = (vg / h2) - vg;
        // let dist_gaus = statrs::distribution::Normal::new(0.0, ve.sqrt()).unwrap();
        // let e: Array2<f64> = Array2::from_shape_vec(
        //     (n, 1),
        //     dist_gaus
        //         .sample_iter(&mut rng)
        //         .take(n)
        //         .collect::<Vec<f64>>(),
        // )
        // .unwrap();
        // let y = &xb + e;

        // let idx_train = (0..90).collect::<Vec<usize>>();
        // let (b_ols, _) = ols(&f, &y, &idx_train).unwrap();
        // println!("b_ols.slice(s![0..10])={:?}", b_ols.slice(s![0..10, ..]));
        // let (b_ridge_like, _) = penalise_ridge_like(&f, &y, &idx_train).unwrap();
        // println!(
        //     "b_ridge_like.slice(s![0..10])={:?}",
        //     b_ridge_like.slice(s![0..10, ..])
        // );
        // let (b_lasso_like, _) = penalise_lasso_like(&f, &y, &idx_train).unwrap();
        // println!(
        //     "b_lasso_like.slice(s![0..10])={:?}",
        //     b_lasso_like.slice(s![0..10, ..])
        // );

        // assert_eq!(0, 1);
    }
}
