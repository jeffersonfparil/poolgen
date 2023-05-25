use crate::base::*;
use crate::gp::*;
use crate::gwas::*;
use ndarray::{prelude::*, Zip};
use rand::prelude::*;
use statrs::statistics::Statistics;
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, Mutex};

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
        let mut intercept = b_hat[(0, j)];
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
        let rmse = mse.sqrt(); // / (max - min);
        error_index.push(((1.0 - cor.abs()) + mae + mse + rmse) / 4.0);
        // error_index.push(((1.0 - cor.abs()) + mae + mse) / 3.0);
        // error_index.push(((1.0 - cor.abs()) + rmse) / 2.0);
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
    let (n, p) = (row_idx.len(), x.ncols());
    let k = y.ncols();
    let max_usize: usize = (1.0 / lambda_step_size).round() as usize;
    let lambda_path: Array1<f64> = (0..(max_usize + 1))
        .into_iter()
        .map(|x| (x as f64) / (max_usize as f64))
        .collect();
    let l = lambda_path.len();

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
                lambda_path
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
        std::iter::repeat(lambda_path)
            .take(a)
            .flat_map(|x| x)
            .collect(),
    )
    .unwrap();

    let (_, nfolds, s) = k_split(row_idx, 10).unwrap();
    let mut performances: Array5<f64> = Array5::from_elem((r, nfolds, a, l, k), f64::NAN);
    let mut effects: Array5<Array1<f64>> =
        Array5::from_elem((r, nfolds, a, l, k), Array1::from_elem(1, f64::NAN));
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

            for i0 in 0..a {
                for i1 in 0..l {
                    for j in 0..k {
                        performances[(rep, fold, i0, i1, j)] = errors[(i0, i1)][j];
                        effects[(rep, fold, i0, i1, j)] = b_hats[(i0, i1)].column(j).to_owned();
                        // reps x folds x alpha x lambda x traits
                    }
                }
            }
        }
    }

    // // Estimate a mean beta accounting for their associated errors
    // let min_error = performances.iter().fold(performances[(0,0,0,0,0)], |min, &x| if x<min{x}else{min});
    // let max_error = performances.iter().fold(performances[(0,0,0,0,0)], |max, &x| if x>max{x}else{max});
    // let one_less_scaled_error = performances.iter().map(|&x| 1.00 - (x-min_error)/(max_error-min_error)).collect::<Vec<f64>>();
    // let sum_perf = one_less_scaled_error.iter().fold(0.0, |sum, &x| sum + x);

    let mut errors = performances.iter().map(|&x| x).collect::<Vec<f64>>();
    errors.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let e = errors.len();
    let idx = (0.01 * e as f64).ceil() as usize;
    let max_error_threshold = errors[idx];
    let mut b_hat_proxy: Array2<f64> = Array2::from_elem((p, k), 0.0);
    // let mut perf = 0.0;
    // // let nall = performances.len();
    // // let (mut b, _) = ols(x, y, row_idx).unwrap();
    for rep in 0..r {
        for fold in 0..nfolds {
            for alfa in 0..a {
                for lam in 0..l {
                    for phe in 0..k {
                        if performances[(rep, fold, alfa, lam, phe)] <= max_error_threshold {
                            for i in 0..p {
                                b_hat_proxy[(i, phe)] += &effects[(rep, fold, alfa, lam, phe)][i]
                                    * (1.00 / (idx as f64 + 1.00));
                            }
                        }
                        // let error = (&performances[(rep, fold, alfa, lam, phe)] - min_error)/(max_error-min_error);
                        // let weight = (1.00-error) / sum_perf;
                        // let b_new: Array1<f64> = b.column(phe).to_owned() + (weight * &effects[(rep, fold, alfa, lam, phe)]);
                        // if (1.00-error) > perf {
                        //     perf = 1.00 - error;
                        //
                        //
                        //
                        // }
                    }
                }
            }
        }
    }
    // println!("sum_perf={}, min_error={:?}; max_error={:?}; b={:?}", sum_perf, min_error, max_error, &b);

    // Find best lambda and estimate effects on the full dataset
    let mean_error_across_reps_and_folds: Array3<f64> = performances
        .mean_axis(Axis(0))
        .unwrap()
        .mean_axis(Axis(0))
        .unwrap();
    // let mean_error_lambdas: Array2<f64> =
    //     mean_error_across_reps_and_folds.mean_axis(Axis(0)).unwrap();
    // let mean_error_alphas: Array2<f64> =
    //     mean_error_across_reps_and_folds.mean_axis(Axis(1)).unwrap();
    let (b_hat, _) = ols(x, y, row_idx).unwrap();
    let mut b_hat_penalised = b_hat.clone();
    let mut alphas = vec![];
    let mut lambdas = vec![];
    for j in 0..k {
        let min_error = mean_error_across_reps_and_folds
            .index_axis(Axis(2), j)
            .iter()
            .fold(mean_error_across_reps_and_folds[(0, 0, j)], |min, &x| {
                if x < min {
                    x
                } else {
                    min
                }
            });

        let ((idx_0, idx_1), _) = mean_error_across_reps_and_folds
            .index_axis(Axis(2), j)
            .indexed_iter()
            .find(|((_i, _j), &x)| x == min_error)
            .unwrap();

        alphas.push(alpha_path[(idx_0, idx_1)]);
        lambdas.push(lambda_path[(idx_0, idx_1)]);

        let b_hat_penalised_2d: Array2<f64> =
            // expand_and_contract(&b_hat, &b_hat, alphas[j], lambdas[j]).unwrap();
            expand_and_contract(&b_hat, &b_hat_proxy, alphas[j], lambdas[j]).unwrap();
        for i in 0..p {
            b_hat_penalised[(i, j)] = b_hat_penalised_2d[(i, j)];
        }
    }
    Ok((b_hat_penalised, alphas, lambdas))
    // Ok((b_hat_proxy, vec![0.0], vec![0.0]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::concatenate;
    #[test]
    fn test_penalised() {
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
    }
}
