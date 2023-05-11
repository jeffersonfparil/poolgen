use crate::base::*;
use crate::gwas::*;
use argmin::core::{self, CostFunction, Executor};
use ndarray::prelude::*;
use ndarray_linalg::solve::Determinant;
use ndarray_linalg::Inverse;
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::f64::consts::PI;
use std::io::{self, Error, ErrorKind};

fn negative_likelihood_normal_distribution_sigma_and_beta(
    params: &Vec<f64>,
    x: &Array2<f64>,
    y: &Array1<f64>,
) -> f64 {
    let n = x.nrows();
    let p = x.ncols();
    assert_eq!(n, y.len());
    assert_eq!(p + 1, params.len()); // including sigma or error variance in the list of parameters
                                     // bound sigma with logit
    let sigma = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
    let betas = Array1::from_vec((&params[1..(p + 1)]).to_owned());
    (n as f64 / 2.00) * (2.00 * PI * sigma).ln()
        + (1.00 / sigma)
            * (y - x.dot(&betas))
                .iter()
                .fold(0.0, |sum_squared, &x| sum_squared + x.powf(2.0))
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

impl Regression for UnivariateMaximumLikelihoodEstimation {
    fn new() -> Self {
        UnivariateMaximumLikelihoodEstimation {
            x: Array2::from_elem((1, 1), f64::NAN),
            y: Array1::from_elem(1, f64::NAN),
            b: Array1::from_elem(1, f64::NAN),
            se: f64::NAN,
            v_b: Array1::from_elem(1, f64::NAN),
            t: Array1::from_elem(1, f64::NAN),
            pval: Array1::from_elem(1, f64::NAN),
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
                (cor, _) = match pearsons_correlation(&self.x.column(i), &self.x.column(j)) {
                    Ok(x) => x,
                    Err(_) => (0.0, f64::NAN),
                };
                if cor.abs() >= 0.99 {
                    self.x.remove_index(Axis(1), j);
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
        let n = self.x.nrows();
        let n_ = self.y.len();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        // self.remove_collinearities_in_x();
        let p = self.x.ncols();
        let mut cost = UnivariateMaximumLikelihoodEstimation::new();
        cost.x = self.x.clone();
        cost.y = self.y.clone();
        let solver = prepare_solver_neldermead(p as f64 + 1.0, 1.0);
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
        let params = res.state().param.clone().unwrap();
        self.se = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
        self.b = Array1::from_vec((&params[1..(p + 1)]).to_owned());
        Ok(self)
    }

    fn estimate_variances(&mut self) -> io::Result<&mut Self> {
        if self.b[0].is_nan() {
            match self.estimate_effects() {
                Ok(x) => x,
                Err(y) => return Err(y),
            };
        }
        let n = self.x.nrows();
        let p = self.x.ncols();
        let n_ = self.y.len();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let xt = self.x.t();
        let vcv: Array2<f64>;
        if n < p {
            let inv_xxt = match (self.x.dot(&xt)).inv() {
                Ok(x) => x,
                Err(_) => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xxt.det().unwrap() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"));
            }
            vcv = self.se * (&xt.dot(&inv_xxt).dot(&inv_xxt).dot(&self.x));
        } else {
            let inv_xtx = match (xt.dot(&self.x)).inv() {
                Ok(x) => x,
                Err(_) => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xtx.det().unwrap() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"));
            }
            vcv = self.se * &inv_xtx;
        }
        self.v_b = Array1::from_elem(p, f64::NAN);
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
        let n = self.x.nrows();
        let p = self.x.ncols();
        let n_ = self.y.len();
        if n != n_ {
            return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
        }
        let d = StudentsT::new(0.0, 1.0, p as f64 - 1.0).unwrap();
        self.t = Array1::from_elem(p, f64::NAN);
        self.pval = Array1::from_elem(p, f64::NAN);
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
        // println!("MLE: {:?}", self);
        Ok(self)
    }
}

fn mle(
    x_matrix: &Array2<f64>,
    y_matrix: &Array2<f64>,
) -> io::Result<(Array2<f64>, Array2<f64>, Array2<f64>)> {
    let n = x_matrix.nrows();
    let mut p = x_matrix.ncols();
    let n_ = y_matrix.nrows();
    let k = y_matrix.ncols();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    let mut beta = Array2::zeros((p, k));
    let mut var_beta = Array2::zeros((p, k));
    let mut pval = Array2::zeros((p, k));
    for j in 0..k {
        let mut mle_regression = UnivariateMaximumLikelihoodEstimation::new();
        mle_regression.x = x_matrix.clone();
        mle_regression.y = y_matrix.column(j).to_owned();
        if p <= 6 {
            // Remove collinearities if we're performing iterative regression
            mle_regression.remove_collinearities_in_x();
        }
        match mle_regression.estimate_significance() {
            Ok(x) => x,
            Err(_) => return Err(Error::new(ErrorKind::Other, "Regression failed.")),
        };
        p = mle_regression.x.ncols();
        for i in 0..p {
            beta[(i, j)] = mle_regression.b[i];
            var_beta[(i, j)] = mle_regression.v_b[i];
            pval[(i, j)] = mle_regression.pval[i];
        }
    }
    Ok((beta, var_beta, pval))
}

pub fn mle_iterate(
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
    match locus_frequencies.sort_by_allele_freq(true) {
        Ok(x) => x,
        Err(_) => return None,
    };
    // // Remove minimum allele, i.e. keep p-1 alleles if p >= 2 so we have enough degrees of freedom to fit the intercept
    if locus_frequencies.matrix.ncols() >= 2 {
        locus_frequencies.matrix.remove_index(Axis(1), 0);
        locus_frequencies.alleles_vector.remove(0);
    }
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let n = locus_frequencies.matrix.nrows();
    let mut p = locus_frequencies.matrix.ncols();
    let m = locus_counts_and_phenotypes.phenotypes.nrows();
    let k = locus_counts_and_phenotypes.phenotypes.ncols();
    if n != m {
        return None;
    }
    // Extract the genotype and phenotypes
    p += 1; // Include the intercept
    let mut x_matrix: Array2<f64> = Array2::ones((n, p));
    for j in 1..p {
        x_matrix
            .column_mut(j)
            .assign(&locus_frequencies.matrix.column(j - 1));
    }
    let y_matrix = locus_counts_and_phenotypes.phenotypes.clone();
    // OLS and compute the p-values associated with each estimate
    let (beta, _var_beta, pval) = match mle(&x_matrix, &y_matrix) {
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
            line.push(parse_f64_roundup_and_own(
                x_matrix.column(i).mean().unwrap(),
                8,
            ));
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
    #[test]
    fn test_ols() {
        // Expected
        let expected_output1: Array1<f64> = Array1::from_vec(vec![-0.73, 5.53, 6.42]);
        let expected_output2: Array1<f64> = Array1::from_vec(vec![0.61, 8.59, 3.99]);
        let expected_output3: Array1<f64> = Array1::from_vec(vec![0.35, 0.59, 0.25]);
        let expected_output4: Array1<f64> = Array1::from_vec(vec![0.08, 0.99, -0.41]);
        let expected_output5: Array1<f64> = Array1::from_vec(vec![0.25, 3.44, 1.60]);
        let expected_output6: Array1<f64> = Array1::from_vec(vec![0.77, 0.80, 0.82]);
        let expected_output7: String = "Chromosome1,12345,A,0.36,Pheno_0,5.528455,0.5856869833119951\nChromosome1,12345,A,0.36,Pheno_1,0.99187,0.8004426037481633\nChromosome1,12345,T,0.24,Pheno_0,6.422764,0.2485036431073504\nChromosome1,12345,T,0.24,Pheno_1,-0.406504,0.8230663350210885\n".to_owned();
        // Inputs
        let x: Array2<f64> = Array2::from_shape_vec(
            (5, 3),
            vec![
                1.0, 0.4, 0.1, 1.0, 0.2, 0.1, 1.0, 0.3, 0.2, 1.0, 0.4, 0.3, 1.0, 0.5, 0.5,
            ],
        )
        .unwrap();
        let y: Array2<f64> = Array2::from_shape_vec(
            (5, 2),
            vec![2.0, 0.5, 1.0, 0.2, 2.0, 0.5, 4.0, 0.0, 5.0, 0.5],
        )
        .unwrap();
        let counts: Array2<u64> =
            Array2::from_shape_vec((5, 3), vec![4, 1, 5, 2, 1, 7, 3, 2, 5, 4, 3, 3, 5, 5, 0])
                .unwrap();
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
        let phenotypes: Array2<f64> = y.clone();
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: phenotypes,
            pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        };
        // Outputs
        println!("x={:?}", x);
        println!("y={:?}", y);
        let (beta, var_beta, pval) = mle(&x, &y).unwrap();
        println!("beta={:?}", beta);
        let p1 = beta.nrows();
        let k1 = beta.ncols();
        let p2 = var_beta.nrows();
        let k2 = var_beta.ncols();
        let ols_line = mle_iterate(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        // Assertions
        assert_eq!(p1, 3); // Output dimensions
        assert_eq!(p1, p2);
        assert_eq!(k1, 2);
        assert_eq!(k1, k2);
        // estimated betas for the first phenotype
        assert_eq!(
            (expected_output1[0] * 100.0).round(),
            (beta[(0, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output1[1] * 100.0).round(),
            (beta[(1, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output1[2] * 100.0).round(),
            (beta[(2, 0)] * 100.0).round()
        );
        // estimated beta variances for the first phenotype
        assert_eq!(
            (expected_output2[0] * 100.0).round(),
            (var_beta[(0, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output2[1] * 100.0).round(),
            (var_beta[(1, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output2[2] * 100.0).round(),
            (var_beta[(2, 0)] * 100.0).round()
        );
        // estimated pvalues for the first phenotype
        assert_eq!(
            (expected_output3[0] * 100.0).round(),
            (pval[(0, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output3[1] * 100.0).round(),
            (pval[(1, 0)] * 100.0).round()
        );
        assert_eq!(
            (expected_output3[2] * 100.0).round(),
            (pval[(2, 0)] * 100.0).round()
        );
        // estimated betas for the second phenotype
        assert_eq!(
            (expected_output4[0] * 100.0).round(),
            (beta[(0, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output4[1] * 100.0).round(),
            (beta[(1, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output4[2] * 100.0).round(),
            (beta[(2, 1)] * 100.0).round()
        );
        // estimated beta variances for the second phenotype
        assert_eq!(
            (expected_output5[0] * 100.0).round(),
            (var_beta[(0, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output5[1] * 100.0).round(),
            (var_beta[(1, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output5[2] * 100.0).round(),
            (var_beta[(2, 1)] * 100.0).round()
        );
        // estimated pvalues for the second phenotype
        assert_eq!(
            (expected_output6[0] * 100.0).round(),
            (pval[(0, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output6[1] * 100.0).round(),
            (pval[(1, 1)] * 100.0).round()
        );
        assert_eq!(
            (expected_output6[2] * 100.0).round(),
            (pval[(2, 1)] * 100.0).round()
        );
        assert_eq!(expected_output7, ols_line);
    }
}
