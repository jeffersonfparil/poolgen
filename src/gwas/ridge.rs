use crate::base::*;
use crate::gwas::*;
use argmin::core::observers;
use argmin::core::{self, CostFunction, Executor};
use nalgebra::{self, DMatrix, DVector};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::io::{self, Error, ErrorKind};

fn ridge_objective_function_lambda_and_beta(
    params: &Vec<f64>,
    x: &DMatrix<f64>,
    y: &DVector<f64>,
    xt: &DMatrix<f64>,
    xxt_or_xtx: &DMatrix<f64>,
) -> f64 {
    let (n, p) = x.shape();
    assert_eq!(n, y.len());
    let sigma2 = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
    let tau2 = bound_parameters_with_logit(&vec![params[1]], f64::EPSILON, 1e9)[0];
    let lambda = sigma2 / tau2;
    // let b = DVector::from_vec((&params[2..(p + 2)]).to_owned());
    // let xt = x.transpose();
    let b = if n < p {
        // &xt * ((x * &xt).add_scalar(lambda)).try_inverse().unwrap() * y
        xt * (xxt_or_xtx.add_scalar(lambda)).try_inverse().unwrap() * y
    } else {
        // ((&xt * x).add_scalar(lambda)).try_inverse().unwrap() * &xt * y
        (xxt_or_xtx.add_scalar(lambda)).try_inverse().unwrap() * xt * y
    };
    (y - (x * &b)).norm_squared() + (lambda * &b).norm_squared()
}

impl CostFunction for UnivariateRidgeRegression {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(ridge_objective_function_lambda_and_beta(
            &p,
            &self.x,
            &self.y,
            &self.xt,
            &self.xxt_or_xtx,
        ))
    }
}

impl Regression for UnivariateRidgeRegression {
    fn new() -> Self {
        UnivariateRidgeRegression {
            x: DMatrix::from_element(1, 1, f64::NAN),
            y: DVector::from_element(1, f64::NAN),
            xt: DMatrix::from_element(1, 1, f64::NAN),
            xxt_or_xtx: DMatrix::from_element(1, 1, f64::NAN),
            b: DVector::from_element(1, f64::NAN),
            sigma2: f64::NAN,
            tau2: f64::NAN,
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
        let cost = self.clone();
        // let mut cost = UnivariateRidgeRegression::new();
        // cost.x = self.x.clone();
        // cost.y = self.y.clone();
        let solver = prepare_solver_neldermead(2.0, 1.0); // estimte simga2 and tau2
        // let res = match Executor::new(cost, solver)
        let res = match Executor::new(cost, solver)
            .configure(|state| state.max_iters(1_000))
            .add_observer(
                observers::SlogLogger::term(),
                observers::ObserverMode::Always,
            )
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
        let params = res.state().param.clone().unwrap();
        self.sigma2 = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
        self.tau2 = bound_parameters_with_logit(&vec![params[1]], f64::EPSILON, 1e9)[0];
        let lambda = self.sigma2 / self.tau2;
        self.xt = self.x.transpose();
        self.b = if n < p {
            &self.xt
                * ((&self.xxt_or_xtx).add_scalar(lambda))
                    .try_inverse()
                    .unwrap()
                * &self.y
        } else {
            ((&self.xxt_or_xtx).add_scalar(lambda))
                .try_inverse()
                .unwrap()
                * &self.xt
                * &self.y
        };
        Ok(self)
    }

    fn estimate_variances(&mut self) -> io::Result<&mut Self> {
        if self.b[0].is_nan() {
            match self.estimate_effects() {
                Ok(x) => x,
                Err(y) => return Err(y),
            };
        }
        let (n, p) = self.x.shape();
        let (n_, _) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let inv_xxt_or_xtx_lambda = self
            .xxt_or_xtx
            .clone()
            .add_scalar(self.sigma2 / self.tau2)
            .try_inverse()
            .unwrap();
        println!("inv_xxt_or_xtx_lambda={:?}", inv_xxt_or_xtx_lambda);
        let vcv = if n < p {
            self.sigma2 * &self.xt * &inv_xxt_or_xtx_lambda * &inv_xxt_or_xtx_lambda * &self.x
        } else {
            self.sigma2 * &inv_xxt_or_xtx_lambda
        };
        self.v_b = DVector::from_element(p, f64::NAN);
        for i in 0..p {
            self.v_b[i] = vcv[(i, i)];
        }
        Ok(self)
    }

    fn estimate_significance(&mut self) -> io::Result<&mut Self> {
        if self.t[0].is_nan() {
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

fn ridge(
    x_matrix: &DMatrix<f64>,
    y_matrix: &DMatrix<f64>,
) -> io::Result<(DMatrix<f64>, DMatrix<f64>, DMatrix<f64>)> {
    let (n, mut p) = x_matrix.shape();
    let (n_, k) = y_matrix.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let mut beta = DMatrix::from_element(p, k, 0.0);
    let mut var_beta = DMatrix::from_element(p, k, 0.0);
    let mut pval = DMatrix::from_element(p, k, 0.0);
    for j in 0..k {
        let mut ridge_regression = UnivariateRidgeRegression::new();
        ridge_regression.x = x_matrix.clone();
        ridge_regression.y = DVector::from_columns(&[y_matrix.column(j)]);
        if p <= 6 {
            // Remove collinearities if we're performing iterative regression
            ridge_regression.remove_collinearities_in_x();
        }
        ridge_regression.xt = ridge_regression.x.transpose();
        if n < p {
            ridge_regression.xxt_or_xtx = &ridge_regression.x * &ridge_regression.xt;
        } else {
            ridge_regression.xxt_or_xtx = &ridge_regression.xt * &ridge_regression.x;
        }
        match ridge_regression.estimate_significance() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Regression failed.")),
        };
        p = ridge_regression.x.ncols();
        for i in 0..p {
            beta[(i, j)] = ridge_regression.b[i];
            var_beta[(i, j)] = ridge_regression.tau2;
            pval[(i, j)] = ridge_regression.pval[i];
        }
    }
    Ok((beta, var_beta, pval))
}

pub fn ridge_iterate(
    locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes,
    filter_stats: &FilterStats,
) -> Option<String> {
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
        .locus_counts
        .filter(filter_stats)
    {
        Ok(x) => x,
        Err(_) => return None,
    };
    let mut locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Sort before we remove the major allele
    let mut locus_frequencies = match locus_frequencies.sort_by_allele_freq(true) {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Keep p-1 alleles if p >= 2 so we have enough degrees of freedom to fit the intercept
    if locus_frequencies.matrix.ncols() >= 2 {
        locus_frequencies.matrix = locus_frequencies.matrix.clone().remove_columns(0, 1);
        locus_frequencies.alleles_vector.remove(0);
    }
    // Extract the genotype and phenotypes
    let mut x_matrix = locus_frequencies.matrix.clone();
    let y_matrix = locus_counts_and_phenotypes.phenotypes.clone();
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, mut p) = x_matrix.shape();
    let (m, k) = y_matrix.shape();
    if n != m {
        return None;
    }
    // Insert intercept
    x_matrix = x_matrix.clone().insert_column(0, 1.0);
    p += 1;
    // OLS and compute the p-values associated with each estimate
    let (beta, _var_beta, pval) = match ridge(&x_matrix, &y_matrix) {
        Ok(x) => x,
        Err(_) => return None,
    };
    // Iterate across alleles
    let first_2_col = vec![
        locus_frequencies.chromosome.clone(),
        locus_frequencies.position.to_string(),
    ];
    let mut line: Vec<String> = vec![];
    for i in 1..p {
        // excluding the intercept
        for j in 0..k {
            line.append(&mut first_2_col.clone());
            line.push(locus_frequencies.alleles_vector[i - 1].clone());
            line.push(parse_f64_roundup_and_own(x_matrix.column(i).mean(), 6));
            line.push("Pheno_".to_string() + &(j.to_string())[..]);
            line.push(parse_f64_roundup_and_own(beta[(i, j)], 6));
            line.push(pval[(i, j)].to_string() + "\n");
        }
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use statrs;
    #[test]
    fn test_ridge() {
        // Expected
        let expected_output1: String = "Chromosome1,12345,A,0.36,Pheno_0,5.528455,0.0000047446528681494016\nChromosome1,12345,A,0.36,Pheno_1,0.99187,0.0000002715235205563715\nChromosome1,12345,T,0.24,Pheno_0,6.422764,0.0000007577768363908888\nChromosome1,12345,T,0.24,Pheno_1,-0.406504,0.0000003484638257944539\n".to_owned();
        // Inputs
        let y: DMatrix<f64> =
            DMatrix::from_row_slice(5, 2, &[2.0, 0.5, 1.0, 0.2, 2.0, 0.5, 4.0, 0.0, 5.0, 0.5]);
        let counts: DMatrix<u64> =
            DMatrix::from_row_slice(5, 3, &[4, 1, 5, 2, 1, 7, 3, 2, 5, 4, 3, 3, 5, 5, 0]);
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![0.2, 0.2, 0.2, 0.2, 0.2],
        };
        let locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned(), "T".to_owned(), "D".to_owned()],
            matrix: counts,
        };
        let phenotypes: DMatrix<f64> = y.clone();
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: phenotypes,
            pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        };
        let mut rng = rand::thread_rng();
        let n = 100;
        let p = 1_000;
        let k = 10;
        let big_x_matrix = DMatrix::from_fn(n, p, |_, _| rng.gen_range(0.0..1.0));
        let mut big_b_true = DMatrix::from_element(p, 1, 0.0);
        let idx: DMatrix<usize> = DMatrix::from_fn(k, 1, |_, _| rng.gen_range(0..p));
        for i in idx.iter() {
            big_b_true[(*i, 0)] = 1.00;
        }
        let gauss_dist = statrs::distribution::Normal::new(0.0, 0.01).unwrap();
        let error = DMatrix::from_vec(
            n,
            1,
            gauss_dist
                .sample_iter(&mut rng)
                .take(n)
                .collect::<Vec<f64>>(),
        );
        let big_y_vector = (&big_x_matrix * &big_b_true) + &error;
        // Outputs
        let ridge_line = ridge_iterate(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        let (beta, _var_beta, _pval) = ridge(&big_x_matrix, &big_y_vector).unwrap();
        let y_hat = &big_x_matrix * &beta;
        let e = (&big_y_vector - &y_hat).sum();
        println!("corr={:?}", e);
        // Assertions
        assert_eq!(expected_output1, ridge_line);
        assert_eq!((e * 1e9).round(), 0.0);
        // assert_eq!(expectd_output1, "".to_owned());
    }
}
