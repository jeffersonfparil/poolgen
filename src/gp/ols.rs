use nalgebra::{self, DMatrix};
use std::io::{self, Error, ErrorKind};

#[function_name::named]
pub fn ols(
    x: &DMatrix<f64>,
    y: &DMatrix<f64>,
    _other_params: &Vec<f64>,
) -> io::Result<(DMatrix<f64>, String)> {
    let (n, p) = x.shape();
    let (n_, _) = y.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let b_hat: DMatrix<f64> = if n < p {
        x.transpose() * (x * x.transpose()).try_inverse().unwrap() * y
    } else {
        (x.transpose() * x).try_inverse().unwrap() * x.transpose() * y
    };
    Ok((b_hat, function_name!().to_owned()))
}
