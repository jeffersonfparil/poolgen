// use crate::base::*;
// use crate::gp::*;
// use crate::gwas::*;
// use ndarray::{prelude::*, Zip};
// use statrs::statistics::Statistics;
// use std::io::{self, Error, ErrorKind};
// use std::sync::{Arc, Mutex};

// #[function_name::named]
// pub fn penalise_lasso_like(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
// ) -> io::Result<(Array2<f64>, String)> {
//     let (b_hat, lambda) = penalised_lambda_path_with_k_fold_cross_validation(
//         x,
//         y,
//         row_idx,
//         false,
//         "Lasso-like".to_owned(),
//         0.1,
//         10,
//     )
//     .unwrap();
//     // println!("##############################");
//     // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
//     Ok((
//         b_hat,
//         function_name!().to_owned() + "-" + &lambda.to_string(),
//     ))
// }

// #[function_name::named]
// pub fn penalise_ridge_like(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
// ) -> io::Result<(Array2<f64>, String)> {
//     let (b_hat, lambda) = penalised_lambda_path_with_k_fold_cross_validation(
//         x,
//         y,
//         row_idx,
//         false,
//         "Ridge-like".to_owned(),
//         0.1,
//         10,
//     )
//     .unwrap();
//     // println!("##############################");
//     // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
//     Ok((
//         b_hat,
//         function_name!().to_owned() + "-" + &lambda.to_string(),
//     ))
// }

// #[function_name::named]
// fn ols_iterative_for_penalisation(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
// ) -> io::Result<(Array2<f64>, String)> {
//     let (n, p) = (row_idx.len(), x.ncols());
//     let m = y.ncols();
//     if x.column(0).sum() < n as f64 {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "Please add the intercept in the X matrix.",
//         ));
//     }
//     let mut b_hat = Array2::from_elem((p, m), f64::NAN); // No need to fill the intercept effect
//     let mut x_sub = Array2::from_elem((n, 2), 1.0);
//     let mut y_sub = Array2::from_elem((n, 1), f64::NAN);
//     for k in 0..m {
//         b_hat[(0, k)] = y.mean_axis(Axis(0)).unwrap()[k]; // set the intercept
//         for i in 0..n {
//             y_sub[(i,0)] = y[(i,k)]; // fill the single column trait array
//         }
//         for j in 1..p {
//             for i in 0..n {
//                 x_sub[(i,1)] = x[(i,j)]; // fill the 2-column genotype array for the locus of interest with an intercept in front
//             }
//             let (b, _) = ols(&x_sub, &y_sub, row_idx).unwrap();
//             b_hat[(j, k)] = b[(1,0)]; // locus effect
//         }
//     }
//     Ok((b_hat, function_name!().to_owned()))
// }

// #[function_name::named]
// pub fn penalise_lasso_like_iterative_base(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
// ) -> io::Result<(Array2<f64>, String)> {
//     let (b_hat, lambda) = penalised_lambda_path_with_k_fold_cross_validation(
//         x,
//         y,
//         row_idx,
//         true,
//         "Lasso-like".to_owned(),
//         0.1,
//         10,
//     )
//     .unwrap();
//     // println!("##############################");
//     // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
//     Ok((
//         b_hat,
//         function_name!().to_owned() + "-" + &lambda.to_string(),
//     ))
// }

// #[function_name::named]
// pub fn penalise_ridge_like_iterative_base(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
// ) -> io::Result<(Array2<f64>, String)> {
//     let (b_hat, lambda) = penalised_lambda_path_with_k_fold_cross_validation(
//         x,
//         y,
//         row_idx,
//         true,
//         "Ridge-like".to_owned(),
//         0.1,
//         10,
//     )
//     .unwrap();
//     // println!("##############################");
//     // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
//     Ok((
//         b_hat,
//         function_name!().to_owned() + "-" + &lambda.to_string(),
//     ))
// }

// fn expand_and_contract(
//     b_hat: &Array2<f64>,
//     norm: String,
//     lambda: f64,
// ) -> io::Result<Array2<f64>> {
//     // Clone b_hat
//     let mut b_hat: Array2<f64> = b_hat.clone();
//     // Exclude the intercept from penalisation
//     let mut intercept = b_hat[(0, 0)];
//     let (p, k) = (b_hat.nrows(), b_hat.ncols());
//     for j in 0..k {
//         // Norm 1 or norm 2
//         let normed: Array1<f64> = if norm == "Lasso-like".to_owned() {
//             b_hat.column(j).map(|&x| x.abs())
//             // b_hat.rows(1, p - 1).map(|x| x.abs())
//         } else if norm == "Ridge-like".to_owned() {
//             b_hat.column(j).map(|&x| x.powf(2.0))
//             // b_hat.rows(1, p - 1).map(|x| x.powf(2.0))
//         } else {
//             return Err(Error::new(
//                 ErrorKind::Other,
//                 "Please enter: 'Lasso-like' or 'Ridge-like' norms.",
//             ));
//         };
//         // Find estimates that will be penalised
//         let normed_max = normed.max();
//         let normed_scaled: Array1<f64> = &normed / normed_max;
//         let idx_penalised = normed_scaled
//             .iter()
//             .enumerate()
//             .filter(|(_, &value)| value < lambda)
//             .map(|(index, _)| index)
//             .collect::<Vec<usize>>();
//         let idx_depenalised = normed_scaled
//             .iter()
//             .enumerate()
//             .filter(|(_, &value)| value >= lambda)
//             .map(|(index, _)| index)
//             .collect::<Vec<usize>>();

//         // Penalise: contract
//         let mut subtracted_penalised = 0.0;
//         let mut added_penalised = 0.0;
//         for i in idx_penalised.into_iter() {
//             if b_hat[(i+1,j)] >= 0.0 {
//                 b_hat[(i+1,j)] -= normed[i];
//                 subtracted_penalised += normed[i];
//             } else {
//                 b_hat[(i+1,j)] += normed[i];
//                 added_penalised += normed[i];
//             }
//         }
//         // Find total depenalised (expanded) values
//         let mut subtracted_depenalised = 0.0;
//         let mut added_depenalised = 0.0;
//         for i in idx_depenalised.clone().into_iter() {
//             if b_hat[(i+1,j)] >= 0.0 {
//                 subtracted_depenalised += normed[i];
//             } else {
//                 added_depenalised += normed[i];
//             }
//         }
//         // Account for the absence of available slots to transfer the contracted effects into
//         if (subtracted_penalised > 0.0) & (subtracted_depenalised == 0.0) {
//             added_penalised -= subtracted_penalised;
//             subtracted_penalised = 0.0;
//         } else if (added_penalised > 0.0) & (added_depenalised == 0.0) {
//             subtracted_penalised -= added_penalised;
//             added_penalised = 0.0;
//         }
//         if (subtracted_penalised < 0.0) | ((subtracted_depenalised == 0.0) & (added_depenalised == 0.0))
//         {
//             intercept += subtracted_penalised;
//             intercept -= added_penalised;
//             subtracted_penalised = 0.0;
//             added_penalised = 0.0;
//         }
//         // Depenalise: expand
//         for i in idx_depenalised.into_iter() {
//             if b_hat[(i+1,j)] >= 0.0 {
//                 b_hat[(i+1,j)] += subtracted_penalised * (normed[i] / subtracted_depenalised);
//             } else {
//                 b_hat[(i+1,j)] -= added_penalised * (normed[i] / added_depenalised);
//             }
//         }
//         // Insert the unpenalised intercept
//         b_hat[(0, k)] = intercept;
//     }
//     Ok(b_hat)
// }

// fn error_index(b_hat: &Array2<f64>, x: &Array2<f64>, y_true: &Array2<f64>) -> io::Result<Vec<f64>> {
//     let (n, p) = (x.nrows(), x.ncols());
//     let (n_, k) = (y_true.nrows(), y_true.ncols());
//     let (p_, k_) = (b_hat.nrows(), b_hat.ncols());
//     if p != p_ {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The X matrix is incompatible with b_hat.",
//         ));
//     }
//     if n != n_ {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The X matrix is incompatible with y.",
//         ));
//     }
//     if k != k_ {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The y matrix/vector is incompatible with b_hat.",
//         ));
//     }
//     let mut error_index: Vec<f64> = Vec::with_capacity(k);
//     for j in 0..k {
//         let y_pred: Array2<f64> = x * b_hat;
//         // Assumes y is a column-vector //TODO: Make it matrix-compatible
//         let n = y_true.len();
//         let min = y_true.min();
//         let max = y_true.max();
//         let (cor, _pval) = pearsons_correlation(
//             &y_true.column(j),
//             &y_pred.column(j),
//         )
//         .unwrap();
//         // let mbe = (y_true - &y_pred).mean() / (max - min);vec![0.0]
//         let mae = (y_true - &y_pred).iter().fold(0.0, |norm, &x| norm + x.abs()) / (max - min);
//         let mse = (y_true - &y_pred).iter().fold(0.0, |norm, &x| norm + x.powf(2.0)) / (max - min).powf(2.0);
//         let rmse = mse.sqrt() / (max - min);
//         error_index[j] = ((1.0 - cor.abs()) + mae + mse + rmse) / 4.0;
//     }
//     Ok(error_index)
// }

// fn k_split(x: &Array2<f64>, mut k: usize) -> io::Result<(Vec<usize>, usize, usize)> {
//     let n = x.nrows();
//     if (k >= n) | (n <= 2) {
//         return Err(Error::new(ErrorKind::Other, "The number of splits, i.e. k, needs to be less than the number of pools, n, and n > 2. We are aiming for fold sizes of 10 or greater."));
//     }
//     let mut s = (n as f64 / k as f64).floor() as usize;
//     while s < 10 {
//         if n < 20 {
//             println!("Warning: number of pools is less than 20, so we're using k=2.");
//             k = 2;
//             s = (n as f64 / k as f64).floor() as usize;
//             break;
//         }
//         k -= 1;
//         s = (n as f64 / k as f64).floor() as usize;
//     }
//     let mut g = (0..k)
//         .flat_map(|x| std::iter::repeat(x).take(s))
//         .collect::<Vec<usize>>();
//     if n - s > 0 {
//         for _i in 0..(n - s) {
//             g.push(k);
//         }
//     }
//     let mut rng = rand::thread_rng();
//     let shuffle = rand::seq::index::sample(&mut rng, n, n)
//         .into_iter()
//         .map(|x| x as usize)
//         .collect::<Vec<usize>>();
//     let mut out: Vec<usize> = Vec::new();
//     for i in 0..n {
//         out.push(g[shuffle[i]]);
//     }
//     Ok((out, k, s))
// }

// fn penalised_lambda_path_with_k_fold_cross_validation(
//     x: &Array2<f64>,
//     y: &Array2<f64>,
//     row_idx: &Vec<usize>,
//     iterative: bool,
//     norm: String,
//     lambda_step_size: f64,
//     r: usize,
// ) -> io::Result<(Array2<f64>, f64)> {
//     let (n, p) = (x.nrows(), x.ncols());
//     if n != y.nrows() {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "The X matrix is incompatible with y.",
//         ));
//     }

//     let max_usize: usize = (1.0 / lambda_step_size).round() as usize;
//     let lambda_path: Vec<f64> = (0..max_usize)
//         .into_iter()
//         .map(|x| (x as f64) / (max_usize as f64))
//         .collect::<Vec<f64>>();
//     let lambda_path = lambda_path.clone();
//     let mut performances: Vec<f64> = vec![];
//     let (_groupings, k, s) = k_split(x, 10).unwrap();
//     for _rep in 0..r {
//         let (groupings, k, s) = k_split(x, 10).unwrap();
//         for fold in 0..k {
//             let idx_validation: Vec<usize> = groupings
//                 .iter()
//                 .enumerate()
//                 .filter(|(_, x)| *x == &fold)
//                 .map(|(i, _)| i)
//                 .collect();
//             let idx_training: Vec<usize> = groupings
//                 .iter()
//                 .enumerate()
//                 .filter(|(_, x)| *x != &fold)
//                 .map(|(i, _)| i)
//                 .collect();
//             let x_training = x.select_rows(&idx_training);
//             let y_training = y.select_rows(&idx_training);
//             let x_validation = x.select_rows(&idx_validation);
//             let y_validation = y.select_rows(&idx_validation);

//             let (b_hat, _model_name) = if iterative {
//                 ols_iterative_for_penalisation(&x_training, &y_training).unwrap()
//             } else {
//                 ols(&x_training, &y_training).unwrap()
//             };

//             // Instantiate thread object for parallel execution
//             let mut thread_objects = Vec::new();
//             // Vector holding all returns from pileup2sync_chunk()
//             let thread_ouputs: Arc<Mutex<Vec<LambdaError>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker

//             for lambda in lambda_path.clone().into_iter() {
//                 // let b_hat_new: Array2<f64> =
//                 //     expand_and_contract(&b_hat, norm.clone(), lambda).unwrap();
//                 // performances
//                 //     .push(error_index(&b_hat_new, &x_validation, &y_validation).unwrap());
//                 let b_hat_clone: Array2<f64> = b_hat.clone();
//                 let norm_clone: String = norm.clone();
//                 let x_matrix_validation_clone: Array2<f64> = x_validation.clone();
//                 let y_matrix_validation_clone: Array2<f64> = y_validation.clone();
//                 let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
//                 let thread = std::thread::spawn(move || {
//                     let b_hat_new: Array2<f64> =
//                         expand_and_contract(&b_hat_clone, norm_clone, lambda).unwrap();
//                     let error = error_index(
//                         &b_hat_new,
//                         &x_matrix_validation_clone,
//                         &y_matrix_validation_clone,
//                     )
//                     .unwrap();
//                     let lambda_error = LambdaError {
//                         lambda: lambda,
//                         error: error,
//                     };
//                     thread_ouputs_clone.lock().unwrap().push(lambda_error);
//                 });
//                 thread_objects.push(thread);
//             }
//             // Waiting for all threads to finish
//             for thread in thread_objects {
//                 let _ = thread.join().expect("Unknown thread error occured.");
//             }
//             // Extract Vec<LambdaError>
//             let mut lambda_errors = vec![];
//             for l_e in thread_ouputs.lock().unwrap().iter() {
//                 lambda_errors.push(l_e.clone());
//             }
//             // Sort by lambda prior to pushing into the performance vector
//             // println!("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
//             // println!("BEFORE: lambda_errors={:?}", lambda_errors);
//             lambda_errors.sort_by(|a, b| a.lambda.partial_cmp(&b.lambda).unwrap());
//             // println!("AFTER: lambda_errors={:?}", lambda_errors);
//             for l_e in lambda_errors.into_iter() {
//                 performances.push(l_e.error);
//             }
//         }
//     }

//     let error_indices_across_folds_x_lambdas: Array2<f64> =
//         Array2::from_row_slice(r * k, lambda_path.len(), &performances);
//     // Array2::from_row_slice(k, lambda_path.len(), &(0..(k*lambda_path.len())).into_iter().map(|x| x as f64).collect::<Vec<f64>>());
//     // Find best lambda and estimate effects on the full dataset
//     let mean_error = error_indices_across_folds_x_lambdas.row_mean();
//     let idx = mean_error
//         .iter()
//         .position(|x| *x == mean_error.min())
//         .unwrap();
//     let lambda = lambda_path[idx];
//     // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
//     // println!("model={:?}; mean_error={:?}; lambda_path={:?}; lambda={:?}", norm, mean_error, lambda_path, lambda);
//     let (b_hat, _model_name) = ols(x, y).unwrap();
//     let b_hat_penalised = expand_and_contract(&b_hat, norm, lambda).unwrap();
//     Ok((b_hat_penalised, lambda))
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::gp::*;
//     use rand::prelude::*;
//     use statrs;
//     #[test]
//     fn test_penalised() {
//         let b: Array2<f64> =
//             Array2::from_column_slice(7, 1, &[5.0, -0.4, 0.0, 1.0, -0.1, 1.0, 0.0]);
//         let new_b: Array2<f64> = expand_and_contract(&b, "Lasso-like".to_owned(), 0.5).unwrap();
//         // let new_b = expand_and_contract(&b, "Ridge-like".to_owned(), 0.5);
//         println!("new_b={:?}", new_b);
//         let expected_output1: Array2<f64> =
//             Array2::from_column_slice(7, 1, &[4.5, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0]);
//         assert_eq!(expected_output1, new_b);
//     }
// }
