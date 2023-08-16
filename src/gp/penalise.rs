use crate::base::*;
use crate::gp::*;
use crate::gwas::*;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::*;
use rand::prelude::*;

use std::io::{self, Error, ErrorKind};

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

fn soft_thresholding(rho_j: f64, z_j: f64, lambda: f64) -> io::Result<f64> {
    if rho_j < -lambda {
        return Ok((rho_j + lambda) / z_j);
    } else if rho_j > lambda {
        return Ok((rho_j - lambda) / z_j);
    } else {
        return Ok(0.0);
    }
}

fn coordinate_descent(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
    lambda: f64,
    convergence_threshold: f64,
    max_iterations: usize,
) -> io::Result<Array2<f64>> {
    let (_n, p) = x.dim();
    let n = row_idx.len();
    let mut y_true = Array2::from_elem((n, 1), f64::NAN);
    for i in 0..n {
        y_true[(i, 0)] = y[(row_idx[i], 0)];
    }
    // let mut beta = Array2::from_elem((p, 1), 0.0);
    // beta[(0, 0)] = y_true.mean().unwrap();
    let col_idx = (0..p).collect::<Vec<usize>>();
    let row_new_idx = (0..row_idx.len()).collect::<Vec<usize>>();
    let mut beta = multiply_views_xx(
        &multiply_views_xtx(
            x,
            &multiply_views_xxt(x, x, row_idx, &col_idx, row_idx, &col_idx)
                .unwrap()
                .pinv()
                .unwrap(),
            row_idx,
            &col_idx,
            &row_new_idx,
            &row_new_idx,
        )
        .unwrap(),
        y,
        &col_idx,
        &row_new_idx,
        row_idx,
        &vec![0],
    )
    .unwrap();
    for iter in 0..max_iterations {
        let mut change_in_beta = 0.0;
        for j in 0..p {
            if (iter > 5) & (beta[(j, 0)].abs() <= 1e-9) {
                beta[(j, 0)] = 0.0;
                continue;
            }
            let col_idx = (0..p).filter(|&i| i != j).collect::<Vec<usize>>();
            let beta_idx = (0..p).filter(|&i| i != j).collect::<Vec<usize>>();
            let yhat_notj =
                multiply_views_xx(x, &beta, row_idx, &col_idx, &beta_idx, &vec![0]).unwrap();
            let error = &y_true - &yhat_notj;
            let rho_j = multiply_views_xtx(
                x,
                &error,
                row_idx,
                &vec![j],
                &(0..error.len()).collect::<Vec<usize>>(),
                &vec![0],
            )
            .unwrap()[(0, 0)];
            let z_j =
                multiply_views_xtx(x, x, row_idx, &vec![j], row_idx, &vec![j]).unwrap()[(0, 0)];
            let new_b_j = soft_thresholding(rho_j, z_j, lambda).unwrap();
            change_in_beta += (new_b_j - beta[(j, 0)]).abs();
            beta[(j, 0)] = new_b_j;
        }
        if change_in_beta <= convergence_threshold {
            break;
        }
    }
    Ok(beta)
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

#[function_name::named]
pub fn penalise_lasso_like(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, alphas, lambdas) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, row_idx, 1.00, false, 0.1, 10)
            .unwrap();
    // println!("##############################");
    // let (train, valid) = linfa_datasets::diabetes().split_with_ratio(0.90);
    // println!("valid={:?}", valid);
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-alphas_"
            + &alphas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_")
            + "-lambdas_"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_"),
    ))
}

#[function_name::named]
pub fn penalise_ridge_like(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, alphas, lambdas) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, row_idx, 0.00, false, 0.1, 10)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-alphas_"
            + &alphas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_")
            + "-lambdas_"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_"),
    ))
}

#[function_name::named]
pub fn penalise_glmnet(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, alphas, lambdas) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, row_idx, -0.1, false, 0.1, 10)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-alphas_"
            + &alphas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_")
            + "-lambdas_"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_"),
    ))
}

#[function_name::named]
pub fn penalise_lasso_like_with_iterative_proxy_norms(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, alphas, lambdas) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, row_idx, 1.00, true, 0.1, 10)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-alphas_"
            + &alphas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_")
            + "-lambdas_"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_"),
    ))
}

#[function_name::named]
pub fn penalise_ridge_like_with_iterative_proxy_norms(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let (b_hat, alphas, lambdas) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, row_idx, 1.00, true, 0.1, 10)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned()
            + "-alphas_"
            + &alphas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_")
            + "-lambdas_"
            + &lambdas
                .iter()
                .map(|&x| x.to_string())
                .collect::<Vec<String>>()
                .join("_"),
    ))
}

fn expand_and_contract(
    b_hat: &Array2<f64>,
    b_hat_proxy: &Array2<f64>,
    alpha: f64,
    lambda: f64,
) -> io::Result<Array2<f64>> {
    // Clone b_hat
    let mut b_hat: Array2<f64> = b_hat.clone();
    let (p, k) = (b_hat.nrows(), b_hat.ncols());
    for j in 0..k {
        //Exclude the intercept from penalisation
        let intercept = b_hat[(0, j)];
        // Norm 1 or norm 2 (exclude the intercept)
        let normed1: Array1<f64> = b_hat.column(j).slice(s![1..p]).map(|&x| x.abs());
        let normed2 = b_hat.column(j).slice(s![1..p]).map(|&x| x.powf(2.0));
        let normed = ((1.00 - alpha) * normed2 / 1.00) + (alpha * normed1);

        // Proxy norm 1 or norm 2 (exclude the intercept) for finding the loci that need to be penalised
        let normed1_proxy: Array1<f64> = b_hat_proxy.column(j).slice(s![1..p]).map(|&x| x.abs());
        let normed2_proxy = b_hat_proxy.column(j).slice(s![1..p]).map(|&x| x.powf(2.0));
        let normed_proxy = ((1.00 - alpha) * normed2_proxy / 1.00) + (alpha * normed1_proxy);

        // // Partition the norms between actual and proxy betas
        // let normed_proxy = normed
        //     .iter()
        //     .zip(&normed_proxy)
        //     .map(|(&x, &y)| (x * (1.00 - 0.5)) + (y * (0.00 - 0.5)))
        //     .collect::<Array1<f64>>();

        // Find estimates that will be penalised using the proxy b_hat norms
        let normed_proxy_max =
            normed_proxy
                .iter()
                .fold(normed_proxy[0], |max, &x| if x > max { x } else { max });
        let normed_proxy_scaled: Array1<f64> = &normed_proxy / normed_proxy_max;
        let idx_penalised = normed_proxy_scaled
            .iter()
            .enumerate()
            .filter(|(_, &value)| value < lambda)
            .map(|(index, _)| index)
            .collect::<Vec<usize>>();
        let idx_depenalised = normed_proxy_scaled
            .iter()
            .enumerate()
            .filter(|(_, &value)| value >= lambda)
            .map(|(index, _)| index)
            .collect::<Vec<usize>>();

        // Penalise: contract using the non-proxy b_hat norms
        let mut subtracted_penalised = 0.0;
        let mut added_penalised = 0.0;
        for i in idx_penalised.into_iter() {
            if b_hat[(i + 1, j)] >= 0.0 {
                if (b_hat[(i + 1, j)] - normed[i]) < 0.0 {
                    subtracted_penalised += b_hat[(i + 1, j)];
                    b_hat[(i + 1, j)] = 0.0;
                } else {
                    subtracted_penalised += normed[i];
                    b_hat[(i + 1, j)] -= normed[i];
                };
            } else {
                if (b_hat[(i + 1, j)] + normed[i]) > 0.0 {
                    added_penalised += b_hat[(i + 1, j)].abs();
                    b_hat[(i + 1, j)] = 0.0;
                } else {
                    added_penalised += normed[i];
                    b_hat[(i + 1, j)] += normed[i];
                }
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
        // if (subtracted_penalised < 0.0)
        //     | ((subtracted_depenalised == 0.0) & (added_depenalised == 0.0))
        // {
        //     intercept += subtracted_penalised;
        //     intercept -= added_penalised;
        //     subtracted_penalised = 0.0;
        //     added_penalised = 0.0;
        // }
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
        let y_true_j: Array2<f64> = Array2::from_shape_vec(
            (n, 1),
            idx_validation
                .iter()
                .map(|&i| y[(i, j)])
                .collect::<Vec<f64>>(),
        )
        .unwrap();
        let b_hat_j: Array2<f64> =
            Array2::from_shape_vec((p, 1), b_hat.column(j).to_owned().to_vec()).unwrap();
        let y_pred_j: Array2<f64> =
            multiply_views_xx(x, &b_hat_j, idx_validation, idx_b_hat, idx_b_hat, &vec![0]).unwrap();
        let min = y_true_j
            .iter()
            .fold(y_true_j[(0, 0)], |min, &x| if x < min { x } else { min });
        let max = y_true_j
            .iter()
            .fold(y_true_j[(0, 0)], |max, &x| if x > max { x } else { max });
        let (cor, _pval) = pearsons_correlation(&y_true_j.column(0), &y_pred_j.column(0)).unwrap();
        // let mbe = (y_true_j - &y_pred_j).mean() / (max - min);vec![0.0]
        let mae = (&y_true_j - &y_pred_j)
            .iter()
            .fold(0.0, |norm, &x| norm + x.abs())
            / (max - min);
        let mse = (&y_true_j - &y_pred_j)
            .iter()
            .fold(0.0, |norm, &x| norm + x.powf(2.0))
            / (max - min).powf(2.0);
        let rmse = mse.sqrt() / (max - min);
        error_index.push(((1.0 - cor.abs()) + mae + mse + rmse) / 4.0);
        // error_index.push(((1.0 - cor.abs()) + mae + mse) / 3.0);
        // error_index.push(((1.0 - cor.abs()) + rmse) / 2.0);
        // error_index.push(1.0 - cor.abs());
        // error_index.push(rmse);
        // error_index.push(mae);
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
    alpha: f64,
    iterative: bool,
    lambda_step_size: f64,
    r: usize,
) -> io::Result<(Array2<f64>, Vec<f64>, Vec<f64>)> {
    let (_n, p) = (row_idx.len(), x.ncols());
    let k = y.ncols();
    let max_usize: usize = (1.0 / lambda_step_size).round() as usize;
    let parameters_path: Array1<f64> = (0..(max_usize + 1))
        .into_iter()
        .map(|x| (x as f64) / (max_usize as f64))
        .collect();
    let l = parameters_path.len();

    let (alpha_path, a): (Array2<f64>, usize) = if alpha >= 0.0 {
        // ridge or lasso optimise for lambda only
        (
            Array2::from_shape_vec((1, l), std::iter::repeat(alpha).take(l).collect()).unwrap(),
            1,
        )
    } else {
        // glmnet optimise for both alpha and lambda
        (
            Array2::from_shape_vec(
                (l, l),
                parameters_path
                    .clone()
                    .iter()
                    .flat_map(|&x| std::iter::repeat(x).take(l))
                    .collect(),
            )
            .unwrap(),
            l,
        )
    };
    let lambda_path: Array2<f64> = Array2::from_shape_vec(
        (a, l),
        std::iter::repeat(parameters_path.clone())
            .take(a)
            .flat_map(|x| x)
            .collect(),
    )
    .unwrap();

    let (_, nfolds, _s) = k_split(row_idx, 10).unwrap();
    let mut performances: Array5<f64> = Array5::from_elem((r, nfolds, a, l, k), f64::NAN);
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
            let (b_hat, _) = ols(&x, &y, &idx_training).unwrap();
            let mut errors: Array2<Vec<f64>> = Array2::from_elem((a, l), vec![]);
            let mut b_hats: Array2<Array2<f64>> =
                Array2::from_elem((a, l), Array2::from_elem((1, 1), f64::NAN));
            if iterative == false {
                Zip::from(&mut errors)
                    .and(&mut b_hats)
                    .and(&alpha_path)
                    .and(&lambda_path)
                    .par_for_each(|err, b, &alfa, &lambda| {
                        let b_hat_new: Array2<f64> =
                            expand_and_contract(&b_hat, &b_hat, alfa, lambda).unwrap();
                        *err = error_index(&b_hat_new, x, y, &idx_validation).unwrap();
                        *b = b_hat_new;
                    });
            } else {
                let (b_hat_proxy, _) =
                    ols_iterative_with_kinship_pca_covariate(x, y, row_idx).unwrap();
                Zip::from(&mut errors)
                    .and(&mut b_hats)
                    .and(&alpha_path)
                    .and(&lambda_path)
                    .par_for_each(|err, b, &alfa, &lambda| {
                        let b_hat_new: Array2<f64> =
                            expand_and_contract(&b_hat, &b_hat_proxy, alfa, lambda).unwrap();
                        *err = error_index(&b_hat_new, x, y, &idx_validation).unwrap();
                        *b = b_hat_new;
                    });
            }

            // Append performances, i.e. error index: f(1-cor, rmse, mae, etc...)
            for i0 in 0..a {
                for i1 in 0..l {
                    for j in 0..k {
                        performances[(rep, fold, i0, i1, j)] = errors[(i0, i1)][j];
                    }
                }
            }
        }
    }

    // Find best alpha, lambda and beta on the full dataset
    // let mean_error_across_reps_and_folds: Array3<f64> = performances
    //     .mean_axis(Axis(0))
    //     .unwrap()
    //     .mean_axis(Axis(0))
    //     .unwrap();
    let (b_hat, _) = ols(x, y, row_idx).unwrap();
    let mut b_hat_penalised = b_hat.clone();
    let mut alphas = vec![];
    let mut lambdas = vec![];
    for j in 0..k {
        ///////////////////////////////////
        // TODO: Account for overfit cross-validation folds, i.e. filter them out, or just use mode of the lambda and alphas?
        let mut alpha_path_counts: Array1<usize> = Array1::from_elem(l, 0);
        let mut lambda_path_counts: Array1<usize> = Array1::from_elem(l, 0);
        for rep in 0..r {
            let mean_error_per_rep_across_folds: Array2<f64> = performances
                .slice(s![rep, .., .., .., j])
                .mean_axis(Axis(0))
                .unwrap();
            let min_error = mean_error_per_rep_across_folds.fold(
                mean_error_per_rep_across_folds[(0, 0)],
                |min, &x| {
                    if x < min {
                        x
                    } else {
                        min
                    }
                },
            );
            // println!("min_error={:?}", min_error);
            // println!("mean_error_per_rep_across_folds={:?}", mean_error_per_rep_across_folds);
            // println!("lambda_path_counts={:?}", lambda_path_counts);
            let ((idx_0, idx_1), _) = mean_error_per_rep_across_folds
                .indexed_iter()
                .find(|((_i, _j), &x)| x == min_error)
                .unwrap();
            // println!("lambda_path[(idx_0, idx_1)]={:?}", lambda_path[(idx_0, idx_1)]);
            for a in 0..l {
                if alpha_path[(idx_0, idx_1)] == parameters_path[a] {
                    alpha_path_counts[a] += 1;
                }
                if lambda_path[(idx_0, idx_1)] == parameters_path[a] {
                    lambda_path_counts[a] += 1;
                }
            }
        }
        // Find the mode alpha and lambda
        // println!("lambda_path_counts={:?}", lambda_path_counts);
        let alpha_max_count = alpha_path_counts.fold(0, |max, &x| if x > max { x } else { max });
        let (alpha_idx, _) = alpha_path_counts
            .indexed_iter()
            .find(|(_a, &x)| x == alpha_max_count)
            .unwrap();
        let lambda_max_count = lambda_path_counts.fold(0, |max, &x| if x > max { x } else { max });
        let (lambda_idx, _) = lambda_path_counts
            .indexed_iter()
            .find(|(_a, &x)| x == lambda_max_count)
            .unwrap();
        alphas.push(parameters_path[alpha_idx]);
        lambdas.push(parameters_path[lambda_idx]);
        ///////////////////////////////////

        // let min_error = mean_error_across_reps_and_folds
        //     .index_axis(Axis(2), j)
        //     .iter()
        //     .fold(mean_error_across_reps_and_folds[(0, 0, j)], |min, &x| {
        //         if x < min {
        //             x
        //         } else {
        //             min
        //         }
        //     });

        // let ((idx_0, idx_1), _) = mean_error_across_reps_and_folds
        //     .index_axis(Axis(2), j)
        //     .indexed_iter()
        //     .find(|((_i, _j), &x)| x == min_error)
        //     .unwrap();

        // alphas.push(alpha_path[(idx_0, idx_1)]);
        // lambdas.push(lambda_path[(idx_0, idx_1)]);

        // Note: Lasso/Ridge and glmnet have different best paramaters even though for example lasso seems to get the best performance while glmnet failed to get the same result even though it should given it searches other alphas including alpha=1.00 in lasso.
        //       This is because of the stochasticity per fold, i.e. glmnet might get an alpha  not equal to 1 which results in better performance.
        // println!("min_error={}; alpha={:?}; lambda={:?}", min_error, alphas, lambdas);
        let b_hat_penalised_2d: Array2<f64> = if iterative == false {
            expand_and_contract(&b_hat, &b_hat, alphas[j], lambdas[j]).unwrap()
        } else {
            let (b_hat_proxy, _) = ols_iterative_with_kinship_pca_covariate(x, y, row_idx).unwrap();
            expand_and_contract(&b_hat, &b_hat_proxy, alphas[j], lambdas[j]).unwrap()
        };

        for i in 0..p {
            b_hat_penalised[(i, j)] = b_hat_penalised_2d[(i, j)];
        }
    }
    // println!("#########################################");
    // println!("alphas={:?}", alphas);
    // println!("lambdas={:?}", lambdas);
    Ok((b_hat_penalised, alphas, lambdas))
    // Ok((b_hat_proxy, vec![0.0], vec![0.0]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::concatenate;
    use rand::distributions::{Bernoulli, Distribution};

    #[test]
    fn test_penalised() {
        let n: usize = 100 as usize;
        let p: usize = 10_000 as usize + 1;
        let intercept: Array2<f64> = Array2::ones((50, 1));
        let d = Bernoulli::new(0.5).unwrap();
        let frequencies = d.sample(&mut rand::thread_rng()) as u64;
        let mut x = Array2::from_shape_fn((n, p), |(i, j)| d.sample(&mut rand::thread_rng()) as u64 as f64);
        for i in 0..n {
            x[(i, 0)] = 1.00
        }
        let mut b: Array2<f64> = Array2::from_elem((1000, 1), 0.0);
        b[(1, 0)] = 1.0;
        b[(10, 0)] = 1.0;
        b[(100, 0)] = 1.0;
        b[(500, 0)] = 1.0;
        let y: Array2<f64> = multiply_views_xx(
            &x,
            &b,
            &(0..50).collect::<Vec<usize>>(),
            &(0..1000).collect::<Vec<usize>>(),
            &(0..1000).collect::<Vec<usize>>(),
            &vec![0],
        )
        .unwrap();
        let row_idx: Vec<usize> = (0..50).collect();

        let (b_hat, name) = penalise_lasso_like(&x, &y, &row_idx).unwrap();
        // assert_eq!(0.009667742247346768, b_hat[(0,0)]);

        let b: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![5.0, -0.4, 0.0, 1.0, -0.1, 1.0, 0.0]).unwrap();
        let new_b: Array2<f64> = expand_and_contract(&b, &b, 1.00, 0.5).unwrap();
        let c: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![5.0, 0.4, 0.0, -1.0, 0.1, -1.0, 0.0]).unwrap();
        let new_c: Array2<f64> = expand_and_contract(&c, &c, 1.00, 0.5).unwrap();
        let expected_output1: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![5.0, 0.0, 0.0, 0.75, 0.0, 0.75, 0.0]).unwrap();
        let expected_output2: Array2<f64> =
            Array2::from_shape_vec((7, 1), vec![5.0, 0.0, 0.0, -0.75, 0.0, -0.75, 0.0]).unwrap();
        assert_eq!(expected_output1, new_b);
        assert_eq!(expected_output2, new_c);

        let b_hat = coordinate_descent(
            &x,
            &y,
            &(0..50).collect::<Vec<usize>>(),
            10.0,
            0.1,
            100 as usize,
        )
        .unwrap();

        println!("b_hat={:?}", b_hat);
        println!("n_non_zero={:?}", b_hat.fold(0, |sum, x| if x.abs() > 1e-7 {sum + 1}else{sum + 0}));
        println!("b_hat[(1, 0)]={:?}", b_hat[(1, 0)]);
        println!("b_hat[(10, 0)]={:?}", b_hat[(10, 0)]);
        println!("b_hat[(100, 0)]={:?}", b_hat[(100, 0)]);
        println!("b_ha500[(500, 0)]={:?}", b_hat[(500, 0)]);

        // assert_eq!(0, 1);
    }
}
