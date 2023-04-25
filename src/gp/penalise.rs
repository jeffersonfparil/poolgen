use crate::base::*;
use crate::gp::*;
use crate::gwas::*;
use nalgebra::{self, DMatrix, DVector};
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, Mutex};

#[function_name::named]
pub fn penalise_lasso_like(
    x: &DMatrix<f64>,
    y: &DMatrix<f64>,
) -> io::Result<(DMatrix<f64>, String)> {
    let (b_hat, lambda) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, "Lasso-like".to_owned(), 0.1)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned() + "-" + &lambda.to_string(),
    ))
}

#[function_name::named]
pub fn penalise_ridge_like(
    x: &DMatrix<f64>,
    y: &DMatrix<f64>,
) -> io::Result<(DMatrix<f64>, String)> {
    let (b_hat, lambda) =
        penalised_lambda_path_with_k_fold_cross_validation(x, y, "Ridge-like".to_owned(), 0.1)
            .unwrap();
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((
        b_hat,
        function_name!().to_owned() + "-" + &lambda.to_string(),
    ))
}

fn expand_and_contract(
    b_hat: &DMatrix<f64>,
    norm: String,
    lambda: f64,
) -> io::Result<DMatrix<f64>> {
    // Clone b_hat
    let mut b_hat: DMatrix<f64> = b_hat.clone();
    // Exclude the intercept from penalisation
    let mut intercept = b_hat[(0, 0)];
    let (p, k) = b_hat.shape();
    // Norm 1 or norm 2
    let normed: DMatrix<f64> = if norm == "Lasso-like".to_owned() {
        b_hat.rows(1, p - 1).map(|x| x.abs())
    } else if norm == "Ridge-like".to_owned() {
        b_hat.rows(1, p - 1).map(|x| x.powf(2.0))
    } else {
        return Err(Error::new(
            ErrorKind::Other,
            "Please enter: 'Lasso-like' or 'Ridge-like' norms.",
        ));
    };
    // Find estimates that will be penalised
    let normed_max = normed.max();
    let normed_scaled: DMatrix<f64> = &normed / normed_max;
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
        if b_hat[i + 1] >= 0.0 {
            b_hat[i + 1] -= normed[(i, 0)];
            subtracted_penalised += normed[(i, 0)];
        } else {
            b_hat[i + 1] += normed[(i, 0)];
            added_penalised += normed[(i, 0)];
        }
    }
    // Find total depenalised values
    let mut subtracted_depenalised = 0.0;
    let mut added_depenalised = 0.0;
    for i in idx_depenalised.clone().into_iter() {
        if b_hat[i + 1] >= 0.0 {
            subtracted_depenalised += normed[(i, 0)];
        } else {
            added_depenalised += normed[(i, 0)];
        }
    }
    // Account for the absence of available slots to transfer the contracted effects into
    if (subtracted_penalised > 0.0) & (subtracted_depenalised == 0.0) {
        added_penalised += subtracted_penalised;
    } else if (added_penalised > 0.0) & (added_depenalised == 0.0) {
        subtracted_penalised -= added_penalised;
    }
    if (subtracted_penalised < 0.0) | ((subtracted_depenalised == 0.0) & (added_depenalised == 0.0))
    {
        intercept += subtracted_penalised;
        intercept += added_penalised;
    }
    // Depenalise: expand
    for i in idx_depenalised.into_iter() {
        if b_hat[i + 1] >= 0.0 {
            b_hat[i + 1] += subtracted_penalised * (normed[(i, 0)] / subtracted_depenalised);
        } else {
            b_hat[i + 1] -= added_penalised * (normed[(i, 0)] / added_depenalised);
        }
    }
    // Insert the unpenalised intercept
    b_hat[(0, 0)] = intercept;
    Ok(b_hat)
}

fn error_index(b_hat: &DMatrix<f64>, x: &DMatrix<f64>, y_true: &DMatrix<f64>) -> io::Result<f64> {
    let (n, p) = x.shape();
    if p != b_hat.nrows() {
        return Err(Error::new(
            ErrorKind::Other,
            "The X matrix is incompatible with b_hat.",
        ));
    }
    if n != y_true.nrows() {
        return Err(Error::new(
            ErrorKind::Other,
            "The X matrix is incompatible with y.",
        ));
    }
    if y_true.ncols() != b_hat.ncols() {
        return Err(Error::new(
            ErrorKind::Other,
            "The y matrix/vector is incompatible with b_hat.",
        ));
    }
    let y_pred: DMatrix<f64> = x * b_hat;
    // Assumes y is a column-vector //TODO: Make it matrix-compatible
    let n = y_true.len();
    let min = y_true.min();
    let max = y_true.max();
    let (cor, _pval) = pearsons_correlation(
        &DVector::from_iterator(n, y_true.column(0).into_iter().map(|x| *x)),
        &DVector::from_iterator(n, y_pred.column(0).into_iter().map(|x| *x)),
    )
    .unwrap();
    // let mbe = (y_true - &y_pred).mean() / (max - min);vec![0.0]
    let mae = (y_true - &y_pred).norm() / (max - min);
    let mse = (y_true - &y_pred).norm_squared() / f64::powf(max - min, 2.0);
    let rmse = mse.sqrt() / (max - min);
    let error_index = ((1.0 - cor) + mae + mse + rmse) / 4.0;
    // let error_index = rmse;
    Ok(error_index)
}

fn k_split(x: &DMatrix<f64>, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
    let (n, _) = x.shape();
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

fn penalised_lambda_path_with_k_fold_cross_validation(
    x: &DMatrix<f64>,
    y: &DMatrix<f64>,
    norm: String,
    lambda_step_size: f64,
) -> io::Result<(DMatrix<f64>, f64)> {
    let (n, p) = x.shape();
    if n != y.nrows() {
        return Err(Error::new(
            ErrorKind::Other,
            "The X matrix is incompatible with y.",
        ));
    }

    let (groupings, k, s) = k_split(x, 10).unwrap();
    let max_usize: usize = (1.0 / lambda_step_size).round() as usize;
    let lambda_path: Vec<f64> = (0..max_usize)
        .into_iter()
        .map(|x| (x as f64) / (max_usize as f64))
        .collect::<Vec<f64>>();

    let mut x_matrix_training = DMatrix::from_element(n - s, p + 1, 1.00); // including intercept in the first column
    let mut x_matrix_validation = DMatrix::from_element(s, p + 1, 1.00); // including intercept in the first column
    let mut y_matrix_training = DMatrix::from_element(n - s, 1, f64::NAN);
    let mut y_matrix_validation = DMatrix::from_element(s, 1, f64::NAN);
    let lambda_path = lambda_path.clone();
    let mut performances: Vec<f64> = vec![];
    for fold in 0..k {
        for j in 0..p {
            let mut idx_training: usize = 0;
            let mut idx_validation: usize = 0;
            for i in 0..n {
                if fold == groupings[i] {
                    x_matrix_validation[(idx_validation, j + 1)] = x[(i, j)]; // Including the intercept, that's why we use j+1
                    if j == 0 {
                        y_matrix_validation[(idx_validation, j)] = y[(i, j)];
                    }
                    idx_validation += 1;
                } else {
                    x_matrix_training[(idx_training, j + 1)] = x[(i, j)]; // Including the intercept, that's why we use j+1
                    if j == 0 {
                        y_matrix_training[(idx_training, j)] = y[(i, j)];
                    }
                    idx_training += 1;
                }
            }
        }
        let (b_hat, _model_name) = ols(&x_matrix_training, &y_matrix_training).unwrap();

        // Instantiate thread object for parallel execution
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs: Arc<Mutex<Vec<LambdaError>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker

        for lambda in lambda_path.clone().into_iter() {
            // let b_hat_new: DMatrix<f64> =
            //     expand_and_contract(&b_hat, norm.clone(), lambda).unwrap();
            // performances
            //     .push(error_index(&b_hat_new, &x_matrix_validation, &y_matrix_validation).unwrap());
            let b_hat_clone: DMatrix<f64> = b_hat.clone();
            let norm_clone: String = norm.clone();
            let x_matrix_validation_clone: DMatrix<f64> = x_matrix_validation.clone();
            let y_matrix_validation_clone: DMatrix<f64> = y_matrix_validation.clone();
            let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
            let thread = std::thread::spawn(move || {
                let b_hat_new: DMatrix<f64> =
                    expand_and_contract(&b_hat_clone, norm_clone, lambda).unwrap();
                let error = error_index(
                    &b_hat_new,
                    &x_matrix_validation_clone,
                    &y_matrix_validation_clone,
                )
                .unwrap();
                let lambda_error = LambdaError {
                    lambda: lambda,
                    error: error,
                };
                thread_ouputs_clone.lock().unwrap().push(lambda_error);
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            let _ = thread.join().expect("Unknown thread error occured.");
        }
        // Extract Vec<LambdaError>
        let mut lambda_errors = vec![];
        for l_e in thread_ouputs.lock().unwrap().iter() {
            lambda_errors.push(l_e.clone());
        }
        // Sort by lambda prior to pushing into the performance vector
        // println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
        // println!("BEFORE: lambda_errors={:?}", lambda_errors);
        lambda_errors.sort_by(|a, b| a.lambda.partial_cmp(&b.lambda).unwrap());
        // println!("AFTER: lambda_errors={:?}", lambda_errors);
        for l_e in lambda_errors.into_iter() {
            performances.push(l_e.error);
        }
    }

    let error_indices_across_folds_x_lambdas: DMatrix<f64> =
        DMatrix::from_row_slice(k, lambda_path.len(), &performances);
    // DMatrix::from_row_slice(k, lambda_path.len(), &(0..(k*lambda_path.len())).into_iter().map(|x| x as f64).collect::<Vec<f64>>());
    // Find best lambda and estimate effects on the full dataset
    let mean_error = error_indices_across_folds_x_lambdas.row_mean();
    let idx = mean_error
        .iter()
        .position(|x| *x == mean_error.min())
        .unwrap();
    let lambda = lambda_path[idx];
    // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
    // println!("model={:?}; mean_error={:?}; lambda_path={:?}; lambda={:?}", norm, mean_error, lambda_path, lambda);
    let (b_hat, _model_name) = ols(x, y).unwrap();
    let b_hat_penalised = expand_and_contract(&b_hat, norm, lambda).unwrap();
    Ok((b_hat_penalised, lambda))
}
