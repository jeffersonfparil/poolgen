use crate::gp::*;
use nalgebra::{self, DMatrix};
use std::io::{self, Error, ErrorKind};

#[function_name::named]
pub fn penalise(
    x: &DMatrix<f64>,
    y: &DMatrix<f64>,
    params: &Vec<f64>,
) -> io::Result<(DMatrix<f64>, String)> {
    let norm = match params[0] {
        1.0 => "Lasso-like",
        2.0 => "Ridge-like",
        _ => return Err(Error::new(ErrorKind::Other, "Invalid penalise params, i.e. for params[0] please use 1.0 or 2.0 for Lasso-like and Ridge-like penalisation, respectively."))
    };
    let lambda = if (params[1] >= 0.0) & (params[1] <= 1.0) {
        params[1]
    } else {
        return Err(Error::new(ErrorKind::Other, "Invalid penalise params, i.e. for params[1] please use a lambda parameter between the inclusive range of 0.0 to 1.0."));
    };
    let (mut b_hat, _model_name) = ols(x, y, params).unwrap();
    let (p, k) = b_hat.shape();
    // Norm 1 or norm 2
    let normed: DMatrix<f64> = if norm == "Lasso-like" {
        DMatrix::from_iterator(p, k, b_hat.iter().map(|x| x.abs()).collect::<Vec<f64>>())
    } else {
        DMatrix::from_iterator(
            p,
            k,
            b_hat.iter().map(|x| x.powf(2.0)).collect::<Vec<f64>>(),
        )
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
        if b_hat[i] >= 0.0 {
            b_hat[i] -= normed[(i, 0)];
            subtracted_penalised += normed[(i, 0)];
        } else {
            b_hat[i] += normed[(i, 0)];
            added_penalised += normed[(i, 0)];
        }
    }
    // Find total depenalised values
    let mut subtracted_depenalised = 0.0;
    let mut added_depenalised = 0.0;
    for i in idx_depenalised.clone().into_iter() {
        if b_hat[i] >= 0.0 {
            subtracted_depenalised += normed[(i, 0)];
        } else {
            added_depenalised += normed[(i, 0)];

        }
    }
    // Depenalise: expand
    for i in idx_depenalised.into_iter() {
        if b_hat[i] >= 0.0 {
            b_hat[i] += subtracted_penalised * (normed[(i, 0)]/subtracted_depenalised);
        } else {
            b_hat[i] -= added_penalised * (normed[(i, 0)]/added_depenalised);
        }
    }
    Ok((b_hat, function_name!().to_owned()))
}
