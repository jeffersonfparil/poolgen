use nalgebra::{self, DMatrix};
use std::io::{self, Error, ErrorKind};

#[function_name::named]
pub fn ols(x: &DMatrix<f64>, y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)> {
    let (n, p) = x.shape();
    let (n_, _) = y.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    if x.column(0).sum() < n as f64 {
        return Err(Error::new(
            ErrorKind::Other,
            "Please add the intercept in the X matrix.",
        ));
    }
    let b_hat: DMatrix<f64> = if n < p {
        x.transpose() * (x * x.transpose()).try_inverse().unwrap() * y
    } else {
        (x.transpose() * x).try_inverse().unwrap() * x.transpose() * y
    };
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((b_hat, function_name!().to_owned()))
}
