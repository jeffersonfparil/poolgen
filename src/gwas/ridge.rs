use crate::base::*;
use crate::gwas::*;
use argmin::core::observers;
use argmin::core::{self, CostFunction, Executor};
use nalgebra::{self, DMatrix, DVector};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, Mutex};

fn multithread_matrix_multiplication_per_row(
    a: &DMatrix<f64>,
    b: &DMatrix<f64>,
    i: usize,
    j: usize,
) -> f64 {
    (a.row(i) * b.column(j))[(0, 0)]
}

fn multithread_matrix_multiplication(a: &DMatrix<f64>, b: &DMatrix<f64>) -> DMatrix<f64> {
    let (n, p) = a.shape();
    let (_p, m) = b.shape();
    // Vector holding all returns from pileup2sync_chunk()
    // let thread_ouputs: Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
    let thread_ouputs: Arc<Mutex<Vec<f64>>> = Arc::new(Mutex::new(Vec::with_capacity(p))); // Mutated within each thread worker
                                                                                // Making four separate threads calling the `search_for_word` function
std::thread::scope(|scope| {
    let mut thread_objects = Vec::new();
    for i in 0..n {
        for j in 0..m {
            let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
            let thread = scope.spawn(move || {
                let z = multithread_matrix_multiplication_per_row(a, b, i, j);
                thread_ouputs_clone.lock().unwrap().push(z);
            });
            thread_objects.push(thread);
        }
    }
    // Waiting for all threads to finish
    for thread in thread_objects {
        let _ = thread.join().expect("Unknown thread error occured.");
    }
});
    
    let mut out = DMatrix::from_element(n, m, f64::NAN);
    let mut i = 0;
    let mut j = 0;
    for z in thread_ouputs.lock().unwrap().iter() {
        println!("z={:?}", z);
        out[(i, j)] = *z;
        if (i+1) < n {
            if (j+1) < m {
                j += 1;
            } else {
                j = 0;
                i += 1;
            }
        } else {
            i = 0;
        }
    }
    out
}

fn ridge_objective_function_lambda_and_beta(
    params: &Vec<f64>,
    x: &DMatrix<f64>,
    y: &DVector<f64>,
) -> f64 {
    let (n, p) = x.shape();
    assert_eq!(n, y.len());
    let sigma2 = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
    let tau2 = bound_parameters_with_logit(&vec![params[1]], f64::EPSILON, 1e9)[0];
    let lambda = sigma2 / tau2;
    // let b = DVector::from_vec((&params[2..(p + 2)]).to_owned());
    let xt = x.transpose();
    let b = if n < p {
        &xt * ((x * &xt).add_scalar(lambda)).try_inverse().unwrap() * y
    } else {
        ((&xt * x).add_scalar(lambda)).try_inverse().unwrap() * &xt * y
    };
    (y - (x * &b)).norm_squared() + (lambda * &b).norm_squared()
}

impl CostFunction for UnivariateRidgeRegression {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(ridge_objective_function_lambda_and_beta(
            &p, &self.x, &self.y,
        ))
    }
}

impl Regression for UnivariateRidgeRegression {
    fn new() -> Self {
        UnivariateRidgeRegression {
            x: DMatrix::from_element(1, 1, f64::NAN),
            y: DVector::from_element(1, f64::NAN),
            b: DVector::from_element(1, f64::NAN),
            sigma2: f64::NAN,
            tau2: f64::NAN,
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
        let (n, _) = self.x.shape();
        let (n_, _) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        self.remove_collinearities_in_x();
        let (_, p) = self.x.shape();
        let mut cost = UnivariateRidgeRegression::new();
        cost.x = self.x.clone();
        cost.y = self.y.clone();
        let solver = prepare_solver_neldermead(2.0, 1.0);
        let res = match Executor::new(cost, solver)
            .configure(|state| state.max_iters(1_000))
            .add_observer(
                observers::SlogLogger::term(),
                observers::ObserverMode::NewBest,
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
        let xt = self.x.transpose();
        self.b = if n < p {
            &xt * ((&self.x * &xt).add_scalar(lambda)).try_inverse().unwrap() * &self.y
        } else {
            ((&xt * &self.x).add_scalar(lambda)).try_inverse().unwrap() * &xt * &self.y
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
        let (n, _) = self.x.shape();
        let (n_, _) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
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
        let (n, p) = self.x.shape();
        let (n_, _k) = self.y.shape();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let d = StudentsT::new(0.0, 1.0, p as f64 - 1.0).unwrap();
        self.t = DVector::from_element(p, f64::NAN);
        self.pval = DVector::from_element(p, f64::NAN);
        for i in 0..p {
            self.t[i] = self.b[i] / self.tau2;
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
    #[test]
    fn test_ridge() {
        // Expected
        let expected_output1: String = "Chromosome1,12345,A,0.36,Pheno_0,5.528455,1\nChromosome1,12345,A,0.36,Pheno_1,0.99187,1\nChromosome1,12345,T,0.24,Pheno_0,6.422764,1\nChromosome1,12345,T,0.24,Pheno_1,-0.406504,1\n".to_owned();
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
        let big_x_matrix = DMatrix::from_fn(100, 32_000, |_, _| rng.gen::<f64>());
        let big_y_vector = DMatrix::from_fn(100, 1, |_, _| rng.gen::<f64>());
        // let big_x_matrix: Array2<f64> = Array::from_shape_fn((100, 1_000), |_| rng.gen_range(0.0..1.0));
        // let big_y_vector: Array2<f64> = Array::from_shape_fn((100, 1), |_| rng.gen_range(-1.0..1.0));
        // println!("test={:?}", big_x_matrix.reversed_axes().shape());
        // println!("test={:?}", big_y_vector.shape());
        // let b = big_x_matrix.transpose() * big_y_vector;
        let b = multithread_matrix_multiplication(&big_x_matrix.transpose(), &big_y_vector);
        // println!("test={:?}", b);
        assert_eq!(b[(0,0)], b[(0,0)]);
        // Outputs
        let ols_line = ridge_iterate(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        // let (beta, _var_beta, _pval) = ridge(&big_x_matrix, &big_y_vector).unwrap();
        // println!("big_beta={:?}", beta);
        // Assertions
        assert_eq!(expected_output1, ols_line);
    }
}
