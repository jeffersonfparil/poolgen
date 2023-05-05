use crate::base::*;
use nalgebra::{self, DMatrix};
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, RwLock};

#[function_name::named]
pub fn ols(genotypes_and_phenotypes: Arc<RwLock<GenotypesAndPhenotypes>>) -> io::Result<(DMatrix<f64>, String)> {
    let data = genotypes_and_phenotypes.read().unwrap();
    let x = &data.intercept_and_allele_frequencies;
    let y = &data.phenotypes;
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
    let lambda = 1.0;
    let b_hat: DMatrix<f64> = if n < p {
        x.transpose()
            * ((x * x.transpose()).add_scalar(lambda))
                .pseudo_inverse(f64::EPSILON)
                .unwrap()
            * y
    } else {
        ((x.transpose() * x).add_scalar(lambda))
            .pseudo_inverse(f64::EPSILON)
            .unwrap()
            * x.transpose()
            * y
    };
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((b_hat, function_name!().to_owned()))
}
