use crate::base::*;
use crate::gwas::*;
use argmin::core::{self, CostFunction, Executor};
use ndarray::{prelude::*, Zip};
use ndarray_linalg::*;
use ndarray_linalg::{solve::Determinant, Inverse};
use statrs::distribution::{ContinuousCDF, StudentsT};
use std::f64::consts::PI;
use std::fs::OpenOptions;
use std::io::{self, prelude::*, Error, ErrorKind};
use std::time::{SystemTime, UNIX_EPOCH};

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
            ve: f64::NAN,
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
        self.ve = bound_parameters_with_logit(&vec![params[0]], f64::EPSILON, 1e9)[0];
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
            vcv = self.ve * (&xt.dot(&inv_xxt).dot(&inv_xxt).dot(&self.x));
        } else {
            let inv_xtx = match (xt.dot(&self.x)).inv() {
                Ok(x) => x,
                Err(_) => return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix")),
            };
            if inv_xtx.det().unwrap() == 0.0 {
                return Err(Error::new(ErrorKind::Other, "Non-invertible x_matrix"));
            }
            vcv = self.ve * &inv_xtx;
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
    remove_collinearities: bool,
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
        if remove_collinearities {
            // Remove collinearities if we're performing iterative regression
            // Note: does not account for the identities of the removed columns - alleles or covariates
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
    // Check struct
    locus_counts_and_phenotypes.check().unwrap();
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
    let (beta, _var_beta, pval): (Array2<f64>, Array2<f64>, Array2<f64>) =
        match mle(&x_matrix, &y_matrix, true) {
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

pub fn mle_with_covariate(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    xxt_eigen_variance_explained: f64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    // Check struct
    genotypes_and_phenotypes.check().unwrap();
    // Generate the covariate
    let g = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .slice(s![0.., 1..]);
    let (n, p) = (g.nrows(), g.ncols());
    let kinship: Array2<f64> = g.dot(&g.t()) / (p as f64);
    let (eigvals, eigvecs) = kinship.eig().unwrap(); // eigenvalues are sorted from high to low
    let sum_eigs = eigvals.iter().map(|&x| x.re).fold(0.0, |sum, x| sum + x);
    let mut n_eigenvecs: usize = n;
    let mut cummulative_variance_explained = eigvals
        .iter()
        .map(|x| x.re / sum_eigs)
        .collect::<Vec<f64>>();
    for i in 1..cummulative_variance_explained.len() {
        cummulative_variance_explained[i] =
            cummulative_variance_explained[i - 1] + cummulative_variance_explained[i];
        if (cummulative_variance_explained[i - 1] >= xxt_eigen_variance_explained)
            & (i - 1 < n_eigenvecs)
        {
            n_eigenvecs = i - 1;
        }
    }
    let covariate: Array2<f64> = eigvecs
        .map(|x| x.re)
        .slice(s![.., 0..n_eigenvecs])
        .to_owned();
    // Phenotype
    let y: Array2<f64> = genotypes_and_phenotypes.phenotypes.to_owned();
    let k = y.ncols();
    // Prepare arrays for parallel computations
    let mut beta: Array2<f64> = Array2::from_elem((p, k), f64::NAN);
    let mut varb: Array2<f64> = Array2::from_elem((p, k), f64::NAN);
    let mut pval: Array2<f64> = Array2::from_elem((p, k), f64::NAN);
    let allele_idx: Array2<usize> = Array2::from_shape_vec(
        (p, k),
        (0..p)
            .flat_map(|x| std::iter::repeat(x).take(k))
            .collect::<Vec<usize>>(),
    )
    .unwrap();
    let phenotype_idx: Array2<usize> = Array2::from_shape_vec(
        (p, k),
        std::iter::repeat((0..k).collect::<Vec<usize>>())
            .take(p)
            .flat_map(|x| x)
            .collect(),
    )
    .unwrap();

    println!("allele_idx={:?}", allele_idx);
    println!("phenotype_idx={:?}", phenotype_idx);
    println!("covariate={:?}", covariate);
    // Parallel OLS
    Zip::from(&mut beta)
        .and(&mut varb)
        .and(&mut pval)
        .and(&allele_idx)
        .and(&phenotype_idx)
        .par_for_each(|effect, variance, significance, &i, &j| {
            let mut x_matrix: Array2<f64> = Array2::ones((n, 2 + n_eigenvecs));
            for i_ in 0..n {
                for j_ in 1..(n_eigenvecs + 1) {
                    x_matrix[(i_, j_)] = covariate[(i_, j_ - 1)];
                }
                x_matrix[(i_, n_eigenvecs + 1)] = g[(i_, i)];
            }
            let y_matrix: Array2<f64> = Array2::from_shape_vec(
                (n, 1),
                y.column(j).iter().map(|x| x.clone()).collect::<Vec<f64>>(),
            )
            .unwrap();
            let (beta_, var_beta_, pval_): (Array2<f64>, Array2<f64>, Array2<f64>) =
                match mle(&x_matrix, &y_matrix, false) {
                    Ok(x) => x,
                    Err(_) => (
                        Array2::from_elem((p, 1), f64::NAN),
                        Array2::from_elem((p, 1), f64::NAN),
                        Array2::from_elem((p, 1), f64::NAN),
                    ),
                };
            *effect = beta_[(n_eigenvecs + 1, 0)];
            *variance = var_beta_[(n_eigenvecs + 1, 0)];
            *significance = pval_[(n_eigenvecs + 1, 0)];
        });

    println!("y={:?}", y);
    println!("g={:?}", g);
    println!("beta={:?}", beta);
    println!("pval={:?}", pval);

    // Write output
    let mut fname_output = fname_output.to_owned();
    if fname_output == "".to_owned() {
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        let bname = fname_input
            .split(".")
            .collect::<Vec<&str>>()
            .into_iter()
            .map(|a| a.to_owned())
            .collect::<Vec<String>>()
            .into_iter()
            .rev()
            .collect::<Vec<String>>()[1..]
            .to_owned()
            .into_iter()
            .rev()
            .collect::<Vec<String>>()
            .join(".");
        fname_output = bname.to_owned()
            + "-mle_iterative_xxt_"
            + &(n_eigenvecs + 1).to_string()
            + "_eigens-"
            + &time.to_string()
            + ".csv";
    }
    // Instatiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
    let mut file_out = OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output)
        .expect(&error_writing_file);
    file_out
        .write_all(("#chr,pos,alleles,phenotype,statistic,pvalue\n").as_bytes())
        .unwrap();
    for j in 0..k {
        for i in 0..p {
            let line = vec![
                genotypes_and_phenotypes.chromosome[i].to_string(),
                genotypes_and_phenotypes.position[i].to_string(),
                genotypes_and_phenotypes.allele[i].clone(),
                j.to_string(),
                beta[(i, j)].to_string(),
                pval[(i, j)].to_string(),
            ]
            .join(",")
                + "\n";
            file_out.write_all(line.as_bytes()).unwrap();
        }
    }

    Ok(fname_output)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_mle() {
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
        let mut genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: vec!["intercept".to_owned(), "X".to_owned(), "Y".to_owned()],
            position: vec![0, 123, 987],
            allele: vec!["intercept".to_owned(), "a".to_owned(), "g".to_owned()],
            intercept_and_allele_frequencies: x.clone(),
            phenotypes: y.clone(),
            pool_names: (0..5)
                .map(|x| "pool-".to_owned() + &x.to_string()[..])
                .collect(),
            coverages: Array2::from_elem((5, 2), 100.0),
        };
        genotypes_and_phenotypes.intercept_and_allele_frequencies[(0, 2)] = 10.0;
        genotypes_and_phenotypes.intercept_and_allele_frequencies[(1, 2)] = 8.0;
        genotypes_and_phenotypes.intercept_and_allele_frequencies[(2, 2)] = 5.0;
        genotypes_and_phenotypes.intercept_and_allele_frequencies[(3, 2)] = 2.0;
        genotypes_and_phenotypes.intercept_and_allele_frequencies[(4, 2)] = 1.0;
        let mle_iterate_with_covariate = mle_with_covariate(
            &genotypes_and_phenotypes,
            0.5,
            &"test-iterative_mle_with_xxt_eigens.sync".to_owned(),
            &"test-iterative_mle_with_xxt_eigens.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(
            mle_iterate_with_covariate,
            "test-iterative_mle_with_xxt_eigens.csv".to_owned()
        );
        // Outputs
        println!("x={:?}", x);
        println!("y={:?}", y);
        let (beta, var_beta, pval) = mle(&x, &y, true).unwrap();
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
        // assert_eq!(expected_output7, ols_line);
    }
}
