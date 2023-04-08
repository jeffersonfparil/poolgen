use crate::base::*;
use crate::gwas::*;
use argmin::core::{self, CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use nalgebra::{self, DMatrix, DVector};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::f64::consts::PI;
use std::io::{self, Error, ErrorKind};

fn negative_likelihood_normal_distribution_sigma_and_beta(
    params: &Vec<f64>,
    x: &DMatrix<f64>,
    y: &DVector<f64>,
) -> f64 {
    let (n, p) = x.shape();
    assert_eq!(n, y.len());
    assert_eq!(p + 1, params.len()); // including sigma or error variance in the list of parameters
                                     // bound sigma with logit
    let sigma = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
    let betas = DVector::from_vec((&params[1..(p + 1)]).to_owned());
    (n as f64 / 2.00) * f64::ln(2.00 * PI * sigma) + (1.00 / sigma) * (y - x * betas).norm_squared()
}

impl CostFunction for UnivariateMaximumLikelihoodEstimation {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(negative_likelihood_normal_distribution_sigma_and_beta(
            &p, &self.x, &self.y,
        ))
    }
}

fn prepare_solver(p: f64, h: f64) -> NelderMead<Vec<f64>, f64> {
    let mut init_param: Vec<Vec<f64>> = Vec::new();
    for i in 0..(p as usize + 1) {
        init_param.push(vec![]);
        for j in 0..p as usize {
            if i == j {
                init_param[i].push(1.5 * h)
            } else {
                init_param[i].push(1.0 * h)
            }
        }
    }
    NelderMead::new(init_param)
}

impl Regression for UnivariateMaximumLikelihoodEstimation {
    fn new() -> Self {
        UnivariateMaximumLikelihoodEstimation {
            x: DMatrix::from_element(1, 1, f64::NAN),
            y: DVector::from_element(1, f64::NAN),
            b: DVector::from_element(1, f64::NAN),
            se: f64::NAN,
            v_b: DVector::from_element(1, f64::NAN),
            t: DVector::from_element(1, f64::NAN),
            pval: DVector::from_element(1, f64::NAN),
        }
    }

    fn remove_collinearities_in_x(&mut self) -> &mut Self {
        if self.x.ncols() == 2 {
            return self;
        }
        let mut i: usize = 1; // exclude the intercept
        let mut j: usize;
        let mut cor: f64;
        while i < self.x.ncols() {
            j = i + 1;
            while j < self.x.ncols() {
                (cor, _) = match pearsons_correlation(
                    &self.x.column(i).into_owned(),
                    &self.x.column(j).into_owned(),
                ) {
                    Ok(x) => x,
                    Err(_) => (0.0, f64::NAN),
                };
                if cor.abs() >= 0.99 {
                    self.x = self.x.clone().remove_column(j);
                    i -= 1;
                    j -= 1;
                }
                j += 1;
            }
            i += 1;
        }
        self
    }

    fn estimate_effects(&mut self) -> io::Result<&mut Self> {
        let (n, p) = self.x.shape();
        let (n_, _) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        self.remove_collinearities_in_x();
        let mut cost = UnivariateMaximumLikelihoodEstimation::new();
        cost.x = self.x.clone();
        cost.y = self.y.clone();
        let solver = prepare_solver(p as f64 + 1.0, 1.0);
        let res = match Executor::new(cost, solver)
            .configure(|state| state.max_iters(1_000))
            // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
            .run()
        {
            Ok(x) => x,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::Other,
                    "T_T Did not converge or something went terribly wrong!",
                ))
            } // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
        };
        // println!("CONVERGENCE: {:?}", res.state());
        let params = res.state().param.clone().unwrap();
        self.se = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
        self.b = DVector::from_vec((&params[1..(p + 1)]).to_owned());
        println!("MLE: {:?}", self);
        Ok(self)
    }

    fn estimate_variances(&mut self) -> io::Result<&mut Self> {
        let (n, p) = self.x.shape();
        let (n_, _) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let xt = self.x.transpose();
        let vcv: DMatrix<f64>;
        if n < p {
            let inv_xxt = match (&self.x * &xt).try_inverse() {
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xxt.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"));
            }
            vcv = self.se * &xt * &inv_xxt * &inv_xxt * &self.x;
        } else {
            let inv_xtx = match (&xt * &self.x).try_inverse() {
                Some(x) => x,
                None => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xtx.determinant() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"));
            }
            vcv = self.se * &inv_xtx;
        }
        self.v_b = DVector::from_element(p, f64::NAN);
        for i in 0..p {
            self.v_b[i] = vcv[(i, i)];
        }
        Ok(self)
    }

    fn estimate_significance(&mut self) -> io::Result<&mut Self> {
        if self.b[0].is_nan() {
            match self.estimate_effects() {
                Ok(x) => x,
                Err(y) => return Err(y),
            };
        }
        if self.v_b[0].is_nan() {
            match self.estimate_variances() {
                Ok(x) => x,
                Err(y) => return Err(y),
            };
        }
        let (n, p) = self.x.shape();
        let (n_, _k) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let d = StudentsT::new(0.0, 1.0, p as f64 - 1.0).unwrap();
        self.t = DVector::from_element(p, f64::NAN);
        self.pval = DVector::from_element(p, f64::NAN);
        for i in 0..p {
            self.t[i] = self.b[i] / self.v_b[i];
            if self.t[i].is_infinite() {
                self.pval[i] = 0.0
            } else if self.t[i].is_nan() {
                self.pval[i] = 1.0
            } else {
                self.pval[i] = 2.00 * (1.00 - d.cdf(self.t[i].abs()));
            }
        }
        Ok(self)
    }
}

// fn mle(x_matrix: &DMatrix<f64>, y_matrix: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, DMatrix<f64>, DMatrix<f64>)> {
//     let (n, mut p) = x_matrix.shape();
//     let (n_, k) = y_matrix.shape();
//     if n != n_ {
//         return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
//     }
//     let mut beta = DMatrix::from_element(p, k, 0.0);
//     let mut var_beta = DMatrix::from_element(p, k, 0.0);
//     let mut pval = DMatrix::from_element(p, k, 0.0);
//     for j in 0..k {
//         let mut ols_regression = UnivariateOrdinaryLeastSquares::new();
//         ols_regression.x = x_matrix.clone();
//         ols_regression.y = DVector::from_columns(&[y_matrix.column(j)]);
//         match ols_regression.estimate_significance() {
//             Ok(x) => x,
//             Err(_) => return Err(Error::new(ErrorKind::Other, "Regression failed."))
//         };
//         p = ols_regression.x.ncols();
//         for i in 0..p {
//             beta[(i,j)]     = ols_regression.b[i];
//             var_beta[(i,j)] = ols_regression.v_b[i];
//             pval[(i,j)]     = ols_regression.pval[i];
//         }
//     }
//     Ok((beta, var_beta, pval))
// }

// pub fn mle_iterate(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
//     // Filter and extract the allele frequencies
//     let locus_counts = match locus_counts_and_phenotypes
//                                                             .locus_counts
//                                                             .filter(filter_stats) {
//         Ok(x) => x,
//         Err(_) => return None
//     };
//     let mut locus_frequencies = match locus_counts.to_frequencies() {
//         Ok(x) => x,
//         Err(_) => return None
//     };
//     // Sort before we remove the major allele
//     let mut locus_frequencies = match locus_frequencies.sort_by_allele_freq(true) {
//         Ok(x) => x,
//         Err(_) => return None
//     };
//     // Keep p-1 alleles if p >= 2 so we have enough degrees of freedom to fit the intercept
//     if locus_frequencies.matrix.ncols() >= 2 {
//         locus_frequencies.matrix = locus_frequencies.matrix.clone().remove_columns(0, 1);
//         locus_frequencies.alleles_vector.remove(0);
//     }
//     // Extract the genotype and phenotypes
//     let mut x_matrix = locus_frequencies.matrix.clone();
//     let y_matrix = locus_counts_and_phenotypes.phenotypes.clone();
//     // Check if we have a compatible allele frequency and phenotype matrix or vector
//     let (n, mut p) =  x_matrix.shape();
//     let (m, k) = y_matrix.shape();
//     if n != m {
//         return None
//     }
//     // Insert intercept
//     x_matrix = x_matrix.clone().insert_column(0, 1.0);
//     p += 1;
//     // OLS and compute the p-values associated with each estimate
//     let (beta, _var_beta, pval) = match ridge(&x_matrix, &y_matrix) {
//         Ok(x) => x,
//         Err(_) => return None,
//     };
//     // Iterate across alleles
//     let first_2_col = vec![locus_frequencies.chromosome.clone(), locus_frequencies.position.to_string()];
//     let mut line: Vec<String> = vec![];
//     for i in 1..p {
//         // excluding the intercept
//         for j in 0..k {
//             line.append(&mut first_2_col.clone());
//             line.push(locus_frequencies.alleles_vector[i-1].clone());
//             line.push(parse_f64_roundup_and_own(x_matrix.column(i).mean(), 8));
//             line.push("Pheno_".to_string() + &(j.to_string())[..]);
//             line.push(parse_f64_roundup_and_own(beta[(i,j)], 8));
//             line.push(pval[(i,j)].to_string() + "\n");
//         }
//     }
//     let out = line.join(",").replace("\n,", "\n");
//     Some(out)
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_ols() {
        // Expected
        let expected_output1 = DVector::from_column_slice(&[
            -0.7317073161994628,
            5.528455290972765,
            6.422764230005505,
        ]);
        // Inputs
        let x: DMatrix<f64> = DMatrix::from_row_slice(
            5,
            3,
            &[
                1.0, 0.4, 0.1, 1.0, 0.2, 0.1, 1.0, 0.3, 0.2, 1.0, 0.4, 0.3, 1.0, 0.5, 0.5,
            ],
        );
        let y: DVector<f64> = DVector::from_column_slice(&[2.0, 1.0, 2.0, 4.0, 5.0]);
        let mut test = UnivariateMaximumLikelihoodEstimation::new();
        test.x = x;
        test.y = y;
        // Outputs
        test.estimate_effects().unwrap();
        // Assertions
        assert_eq!(test.b, expected_output1);
    }
}
