use nalgebra::{self, DMatrix};
use argmin::core::{self, CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use crate::base::*;

use statrs::distribution::{Beta, ContinuousCDF};

const PARAMETER_UPPER_LIMIT: f64 = 10.00;
const PARAMETER_LOWER_LIMIT: f64 = f64::EPSILON;

fn bound_parameters_with_logit(params: &Vec<f64>) -> Vec<f64> {
    // Map parameters with a logistic regression to bound them between 0 and PARAMETER_UPPER_LIMIT
    params.into_iter()
        .map(|x| PARAMETER_LOWER_LIMIT + 
                          ( (PARAMETER_UPPER_LIMIT-PARAMETER_LOWER_LIMIT) / 
                            (1.00 + (-x).exp()) 
                          ) 
            )
        .collect::<Vec<f64>>()
}

fn least_squares_beta(params: &Vec<f64>, percs_a: &DMatrix<f64>, percs_b: &DMatrix<f64>, q_prime: &DMatrix<f64>) -> f64 {
    let shapes = bound_parameters_with_logit(params);
    // println!("shapes={:?}", shapes);
    let a_dist = Beta::new(shapes[0], shapes[1]).expect(&shapes.clone().into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join("-")[..]);
    let b_dist = Beta::new(shapes[2], shapes[3]).unwrap();
    let n = percs_a.len();
    // if (n != percs_b.len()) | (n != q_prime.len()) {
    //     return Err(Error::new(ErrorKind::Other, "percs_a, percs_b, and q_prime are not the same sizes"))
    // }
    let mut a_sum_of_squares: f64 = 0.0;
    let mut b_sum_of_squares: f64 = 0.0;
    for i in 0..n {
        a_sum_of_squares += f64::powf(percs_a[i] - a_dist.cdf(q_prime[i]), 2.0);
        b_sum_of_squares += f64::powf(percs_b[i] - b_dist.cdf(q_prime[i]), 2.0);
    }
    let out = a_sum_of_squares + b_sum_of_squares;
    // println!("a1={}; a2={}; b1={}; b2={}; out={}", a_shape1, a_shape2, b_shape1, b_shape2, out);
    out
}

fn maximum_likelihood_beta(params: &Vec<f64>, percs_a: &DMatrix<f64>, percs_b: &DMatrix<f64>, percs_a0: &DMatrix<f64>, percs_b0: &DMatrix<f64>) -> f64 {
    let shapes = bound_parameters_with_logit(params);
    let a_dist = Beta::new(shapes[0], shapes[1]).unwrap();
    let b_dist = Beta::new(shapes[2], shapes[3]).unwrap();
    let n = percs_a.len();
    let mut a_log_likelihood: f64 = 0.0;
    let mut b_log_likelihood: f64 = 0.0;
    let mut diff_a: f64;
    let mut diff_b: f64;
    for i in 0..n {
        diff_a = a_dist.cdf(percs_a[i]) - a_dist.cdf(percs_a0[i]);
        diff_b = b_dist.cdf(percs_b[i]) - b_dist.cdf(percs_b0[i]);
        diff_a = match diff_a < f64::EPSILON {
            true => f64::EPSILON,
            false => diff_a
        };
        diff_b = match diff_b < f64::EPSILON {
            true => f64::EPSILON,
            false => diff_b
        };
        a_log_likelihood += f64::log10(diff_a);
        b_log_likelihood += f64::log10(diff_b);
        // println!("percs_a[i]={:?}; a_dist.ln_pdf(percs_a[i])={:?}", a_dist.ln_pdf(percs_a[i]), percs_a[i]);
        // println!("percs_b[i]={:?}; b_dist.ln_pdf(percs_b[i])={:?}", b_dist.ln_pdf(percs_b[i]), percs_b[i]);
        // a_log_likelihood += a_dist.ln_pdf(percs_a[i]);
        // b_log_likelihood += b_dist.ln_pdf(percs_b[i]);
    }
    let out = -a_log_likelihood - b_log_likelihood;
    // println!("a1={}; a2={}; b1={}; b2={}; a_log_likelihood={}; b_log_likelihood={}; out={}", shapes[0], shapes[1], shapes[2], shapes[3], a_log_likelihood, b_log_likelihood, out);
    out
}

impl CostFunction for LeastSquaresBeta {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(least_squares_beta(&p, &self.percs_a, &self.percs_b, &self.q_prime))
    }
}

impl CostFunction for MaximumLikelihoodBeta {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(maximum_likelihood_beta(&p, &self.percs_a, &self.percs_b, &self.percs_a0, &self.percs_b0))
    }
}

fn prepare_solver(h: f64) -> NelderMead<Vec<f64>, f64> {
    // let h = 1.0;
    let init_param: Vec<Vec<f64>> = vec![vec![1.00, 1.00, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                         vec![1.05, 1.00, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                         vec![1.00, 1.05, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                         vec![1.00, 1.00, 1.05, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                         vec![1.00, 1.00, 1.00, 1.05].into_iter().map(|x| x * h).collect::<Vec<f64>>()];
    NelderMead::new(init_param)
}

fn gwalpha_minimise_ls(solver: NelderMead<Vec<f64>, f64>, q_prime: DMatrix<f64>, percs_a: DMatrix<f64>, percs_b: DMatrix<f64>) -> Option<Vec<f64>> {
    let cost = LeastSquaresBeta { q_prime: q_prime.clone(),
                                                    percs_a: percs_a.clone(),
                                                    percs_b: percs_b.clone() };
    let res = match Executor::new(cost, solver)
            .configure(|state| {
                state
                    .max_iters(1_000)
            })
            // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
            .run() {
                Ok(x) => x,
                Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
            };
    // println!("CONVERGENCE: {:?}", res.state());
    let params = res.state().param.clone().unwrap();
    let solution = bound_parameters_with_logit(&params);
    Some(solution)
}

fn gwalpha_minimise_ml(solver: NelderMead<Vec<f64>, f64>, percs_a: DMatrix<f64>, percs_a0: DMatrix<f64>, percs_b: DMatrix<f64>, percs_b0: DMatrix<f64>) -> Option<Vec<f64>> {
    let cost = MaximumLikelihoodBeta { percs_a: percs_a,
                                                              percs_b: percs_b,
                                                              percs_a0: percs_a0,
                                                              percs_b0: percs_b0 };
    let res = match Executor::new(cost, solver)
            .configure(|state| {
                state
                    .max_iters(1_000)
            })
            // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
            .run() {
                Ok(x) => x,
                Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
            };
    // println!("CONVERGENCE: {:?}", res.state());
    let params = res.state().param.clone().unwrap();
    let solution = bound_parameters_with_logit(&params);
    Some(solution)
}

fn prepare_geno_and_pheno_stats(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<(LocusFrequencies, DMatrix<f64>, DMatrix<f64>, f64, f64, f64, usize, usize)> {
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
                                                            .locus_counts
                                                            .filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    };
    let mut locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None
    };
    // Sort before we remove the major allele
    let mut locus_frequencies = match locus_frequencies.sort_by_allele_freq(true) {
        Ok(x) => x,
        Err(_) => return None
    };
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if locus_frequencies.matrix.ncols() >= 2 {
        locus_frequencies.matrix = locus_frequencies.matrix.clone().remove_columns(0, 1);
        locus_frequencies.alleles_vector.remove(0);
    }
    // Extract phenotype information (Note: removing NEG_INFINITY from the bins and q columns if we have less than 3 pools, each corresponds to sig, MIN, and MAX rows with 3 as the minimum number of rows)
    let bins_tmp = locus_counts_and_phenotypes.phenotypes.column(0)
                                                                .into_iter()
                                                                .filter(|x| **x != f64::NEG_INFINITY)
                                                                .map(|x| x.to_owned())
                                                                .collect::<Vec<f64>>();
    let q_tmp = locus_counts_and_phenotypes.phenotypes.column(1)
                                                             .into_iter()
                                                             .filter(|x| **x != f64::NEG_INFINITY)
                                                             .map(|x| x.to_owned())
                                                             .collect::<Vec<f64>>();
    let bins = DMatrix::from_vec(bins_tmp.len(), 1, bins_tmp);
    let q = DMatrix::from_vec(q_tmp.len(), 1, q_tmp);
    let sig = locus_counts_and_phenotypes.phenotypes[(0,2)];
    let min = locus_counts_and_phenotypes.phenotypes[(1,2)];
    let max = locus_counts_and_phenotypes.phenotypes[(2,2)];
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let (n, p) =  locus_frequencies.matrix.shape();
    let m = bins.len();
    if n != m {
        return None
    }
    Some((locus_frequencies.clone(),
          bins,
          q,
          sig,
          min,
          max,
          n,
          p
        ))
}

fn prepare_freqs_and_qprime(locus_frequencies: &LocusFrequencies, bins: &DMatrix<f64>, q: &DMatrix<f64>, min: f64, max: f64, n: usize, j: usize) -> (f64, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    // let freqs_a: DMatrix<f64> = DMatrix::from_columns(&[locus_frequencies.matrix.column(j)]).add_scalar(1e-6); // Add a small value so we don't get negative log-likelihoods --> NO NEED HERE BECAUSE WE ARE BOUNDING OUR PARAMETER VALUES WITH LOGISTIC REGRESSION AND BOUNDING THE LOG-DIFFERENCE SO WE DON'T GET INFINITIES
    let freqs_a: DMatrix<f64> = DMatrix::from_columns(&[locus_frequencies.matrix.column(j)]);
    // let freqs_a_sum = freqs_a.sum();
    // for i in 0..freqs_a.len() {
    //     freqs_a[i] = freqs_a[i]/freqs_a_sum;
    // }
    let p_a = (freqs_a.clone().transpose() * bins.clone())[(0,0)]; // mean allele frequency across pools
    // println!("p_a={:?}", p_a);
    // Quantiles per pool (for least squares estimation)
    let mut q_prime: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
    for i in 1..n {
        q_prime[i] = (q[i] - min) / (max-min);
    }
    // println!("q_prime={:?}", q_prime);
    // Bins (sums up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
    let mut bins_a = DMatrix::from_element(n, 1, 0.0);
    let mut bins_b = DMatrix::from_element(n, 1, 0.0);
    for i in 0..n {
        bins_a[i] = (freqs_a[i]) * bins[i] / (p_a);
        bins_b[i] = (1.0 - freqs_a[i]) * bins[i] / (1.0 - p_a);
    }
    // println!("bins_a={:?}", bins_a);
    // println!("bins_b={:?}", bins_b);
    // Percentiles (cummulative bins summing up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
    let mut percs_a = bins_a.clone();
    let mut percs_b = bins_b.clone();
    for i in 1..bins_a.nrows() {
        percs_a[i] = bins_a.view((0,0), (i+1,1)).sum();
        percs_b[i] = bins_b.view((0,0), (i+1,1)).sum();
    }
    // println!("percs_a={:?}", percs_a);
    // println!("percs_b={:?}", percs_b);
    // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
    let mut percs_a0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
    let mut percs_b0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
    for i in 0..n-1 {
        percs_a0[i+1] = percs_a[i];
        percs_b0[i+1] = percs_b[i];
    }
    (p_a, q_prime, percs_a, percs_a0, percs_b, percs_b0)
}

pub fn gwalpha_ls(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // Prepare phenotype and genotype statistics
    let (locus_frequencies, bins, q, sig, min, max, n, p) = 
        match prepare_geno_and_pheno_stats(locus_counts_and_phenotypes, filter_stats) {
            Some(x) => x,
            None => return None
    };
    // Prepare output line
    let first_2_col = vec![locus_frequencies.chromosome.to_owned(), locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    // Instiate input optimisation data (NOTE: exluding percs_a0 and percs_b0)
    let (mut p_a, mut q_prime, mut percs_a, mut percs_b);
    // Instiate solver
    let mut solver: NelderMead<Vec<f64>, f64>;
    // Iterate across alleles
    for j in 0..p {
        // Prepare allele frequecies, quantiles and percentiles
        (p_a, q_prime, percs_a, _, percs_b, _) = prepare_freqs_and_qprime(&locus_frequencies, &bins, &q, min, max, n, j);
        // Optimise
        solver = prepare_solver(1.0);
        let solution = match gwalpha_minimise_ls(solver, q_prime, percs_a, percs_b) {
            Some(x) => x,
            None => return None
        };
        let a_mu_hat = min + (max-min) * (solution[0]/(solution[0]+solution[1]));
        let b_mu_hat = min + (max-min) * (solution[2]/(solution[2]+solution[3]));
        let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(locus_frequencies.matrix.column(j).mean().to_string());
        line.push("Pheno_0".to_string());
        line.push(alpha.to_string() + &",Unknown\n");
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

pub fn gwalpha_ml(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // Prepare phenotype and genotype statistics
    let (locus_frequencies, bins, q, sig, min, max, n, p) = 
        match prepare_geno_and_pheno_stats(locus_counts_and_phenotypes, filter_stats) {
        Some(x) => x,
        None => return None
    };
    // Prepare output line
    let first_2_col = vec![locus_frequencies.chromosome.to_owned(), locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    // Instiate input optimisation data (NOTE: exluding q_prime)
    let (mut p_a, mut percs_a, mut percs_a0, mut percs_b, mut percs_b0);
    // Instiate solver
    let mut solver: NelderMead<Vec<f64>, f64>;
    // Iterate across alleles
    for j in 0..p {
        // Prepare allele frequecies, quantiles and percentiles
        (p_a, _, percs_a, percs_a0, percs_b, percs_b0) = prepare_freqs_and_qprime(&locus_frequencies, &bins, &q, min, max, n, j);
        // Optimise
        solver = prepare_solver(1.0);
        let solution = match gwalpha_minimise_ml(solver, percs_a, percs_a0, percs_b, percs_b0) {
            Some(x) => x,
            None => return None
        };
        let a_mu_hat = min + (max-min) * (solution[0]/(solution[0]+solution[1]));
        let b_mu_hat = min + (max-min) * (solution[2]/(solution[2]+solution[3]));
        let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // let alpha =  (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(locus_frequencies.matrix.column(j).mean().to_string());
        line.push("Pheno_0".to_string());
        line.push(alpha.to_string() + &",Unknown\n");

    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}
