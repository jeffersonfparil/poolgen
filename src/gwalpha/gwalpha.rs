use std::f32::EPSILON;

// use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix};

use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{self, CostFunction, Gradient, Executor, Solver, State};
use finitediff::FiniteDiff;
use argmin::solver::neldermead::NelderMead;
use argmin::solver::linesearch::{MoreThuenteLineSearch, BacktrackingLineSearch, HagerZhangLineSearch};
use argmin::solver::quasinewton::BFGS;
use argmin::solver::quasinewton::LBFGS;
use argmin::solver::particleswarm::ParticleSwarm;
use crate::base::*;

use statrs::distribution::{Beta, ContinuousCDF};

struct LeastSquaresBeta {
    percs_a: DMatrix<f64>,
    percs_b: DMatrix<f64>,
    q_prime: DMatrix<f64>,
}

struct MaximumLikelihoodBeta {
    percs_a: DMatrix<f64>,
    percs_b: DMatrix<f64>,
    percs_a0: DMatrix<f64>,
    percs_b0: DMatrix<f64>,
}

fn search_params_to_shapes(params: &Vec<f64>) -> Vec<f64> {
    const UPPER_LIMIT: f64 = 10.00;
    let shapes = params.into_iter()
                                 .map(|x| 1e-8 + (UPPER_LIMIT / (1.00 + (-x).exp())).abs() )
                                 .collect::<Vec<f64>>();
    shapes
}

fn least_squares_beta(params: &Vec<f64>, percs_a: &DMatrix<f64>, percs_b: &DMatrix<f64>, q_prime: &DMatrix<f64>) -> f64 {
    let shapes = search_params_to_shapes(params);
    // println!("shapes={:?}", shapes);
    let a_dist = Beta::new(shapes[0], shapes[1]).expect(&params.clone().into_iter().map(|x| x.to_string()).collect::<Vec<String>>().join("-")[..]);
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
    let shapes = search_params_to_shapes(params);
    let a_dist = Beta::new(shapes[0], shapes[1]).unwrap();
    let b_dist = Beta::new(shapes[2], shapes[3]).unwrap();
    let n = percs_a.len();
    let mut a_log_likelihood: f64 = 0.0;
    let mut b_log_likelihood: f64 = 0.0;
    for i in 0..n {
        a_log_likelihood += f64::log10(a_dist.cdf(percs_a[i]) - a_dist.cdf(percs_a0[i]));
        b_log_likelihood += f64::log10(b_dist.cdf(percs_b[i]) - b_dist.cdf(percs_b0[i]));
    }
    let out = -a_log_likelihood - b_log_likelihood;
    // println!("a1={}; a2={}; b1={}; b2={}; a_log_likelihood={}; b_log_likelihood={}; out={}", a_shape1, a_shape2, b_shape1, b_shape2, a_log_likelihood, b_log_likelihood, out);
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

impl Gradient for LeastSquaresBeta {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, core::Error> {
        Ok((*p).forward_diff(&|x| least_squares_beta(&x, &self.percs_a, &self.percs_b, &self.q_prime)))
    }
}

impl Gradient for MaximumLikelihoodBeta {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Self::Param) -> Result<Self::Gradient, core::Error> {
        Ok((*p).forward_diff(&|x| maximum_likelihood_beta(&x, &self.percs_a, &self.percs_b, &self.percs_a0, &self.percs_b0)))
    }
}

pub fn gwalpha_ls(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
                                                            .locus_counts
                                                            .filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None
    };
    // Extract the genotype matrix
    let mut x_matrix = (*locus_frequencies).matrix.clone();
    let (_n, p) = x_matrix.shape();
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if p >= 2 {
        x_matrix = x_matrix.clone().remove_columns(p-1, 1);
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
    let (n, p) =  x_matrix.shape();
    let m = bins.len();
    if n != m {
        return None
    }
    // Prepare output line
    let first_2_col = vec![locus_frequencies.chromosome, locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    // Iterate across alleles
    for j in 0..p {
        let mut freqs_a: DMatrix<f64> = DMatrix::from_columns(&[x_matrix.column(j)]).add_scalar(1e-4); // Add a small value so we don't get negative log-likelihoods
        let freqs_a_sum = freqs_a.sum();
        for i in 0..freqs_a.len() {
            freqs_a[i] = freqs_a[i]/freqs_a_sum;
        }
        let p_a = (freqs_a.clone().transpose() * bins.clone())[(0,0)]; // mean allele frequency across pools
        // Quantiles per pool (for least squares estimation)
        let mut q_prime: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 1..n {
            q_prime[i] = (q[i] - min) / (max-min);
        }
        // Bins (sums up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut bins_a = DMatrix::from_element(n, 1, 0.0);
        let mut bins_b = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n {
            bins_a[i] = (freqs_a[i]) * bins[i] / (p_a);
            bins_b[i] = (1.0 - freqs_a[i]) * bins[i] / (1.0 - p_a);
        }
        // Percentiles (cummulative bins summing up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut percs_a = bins_a.clone();
        let mut percs_b = bins_b.clone();
        for i in 1..bins_a.nrows() {
            percs_a[i] = bins_a.view((0,0), (i+1,1)).sum();
            percs_b[i] = bins_b.view((0,0), (i+1,1)).sum();
        }
        // Remap 
        // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
        let mut percs_a0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        let mut percs_b0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n-1 {
            percs_a0[i+1] = percs_a[i];
            percs_b0[i+1] = percs_b[i];
        }
        // Define cost function
        let cost = LeastSquaresBeta { percs_a: percs_a.clone(),
                                                        percs_b: percs_b.clone(),
                                                        q_prime: q_prime.clone() };
        // Define initial parameter vector
        let init_param: Vec<f64> = vec![1.0, 1.0, 1.0, 1.0];
        // set up a line search
        let linesearch = MoreThuenteLineSearch::new();
        // Set up solver
        let solver = LBFGS::new(linesearch, 5);
        // Run solver
        let res = match Executor::new(cost, solver)
                .configure(|state| {
                    state
                        .param(init_param)
                        .max_iters(100)
                })
                // .add_observer(SlogLogger::term(), ObserverMode::Never)
                .run() {
                    Ok(x) => x,
                    Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
                };
        // Calculate the means of the best fit beta distributions describing the distribution of the focal and alternative allele frequences across pools
        let solution = res.state().param.clone().unwrap();
        let a_mu_hat = min + (max-min) * (solution[0]/(solution[0]+solution[1]));
        let b_mu_hat = min + (max-min) * (solution[2]/(solution[2]+solution[3]));
        let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(x_matrix.column(j).mean().to_string());
        line.push("Pheno_0".to_string());
        line.push(alpha.to_string() + &",Unknown\n");

    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}


pub fn gwalpha_ml(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // Filter and extract the allele frequencies
    let locus_counts = match locus_counts_and_phenotypes
                                                            .locus_counts
                                                            .filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    }; // filter by coverage in addition to MAF
    let mut locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None
    }; // Convert counts to frequencies
    match locus_frequencies.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None
    }; // Refilter by MAF with recomputed frequencies
    // Extract the genotype matrix
    locus_frequencies.sort_by_allele_freq(false).unwrap();
    let mut x_matrix = (*locus_frequencies).matrix.clone();
    let (_n, p) = x_matrix.shape();
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if p >= 2 {
        x_matrix = x_matrix.clone().remove_columns(p-1, 1);
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
    let (n, p) =  x_matrix.shape();
    let m = bins.len();
    if n != m {
        return None
    }
    // Prepare output line
    let first_2_col = vec![locus_frequencies.chromosome, locus_frequencies.position.to_string()];
    let mut line: Vec<String> = vec![];
    // Iterate across alleles
    for j in 0..p {
        let mut freqs_a: DMatrix<f64> = DMatrix::from_columns(&[x_matrix.column(j)]).add_scalar(1e-4); // Add a small value so we don't get negative log-likelihoods
        let freqs_a_sum = freqs_a.sum();
        for i in 0..freqs_a.len() {
            freqs_a[i] = freqs_a[i]/freqs_a_sum;
        }
        let p_a = (freqs_a.clone().transpose() * bins.clone())[(0,0)]; // mean allele frequency across pools
        // Quantiles per pool (for least squares estimation)
        let mut q_prime: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 1..n {
            q_prime[i] = (q[i] - min) / (max-min);
        }
        // Bins (sums up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut bins_a = DMatrix::from_element(n, 1, 0.0);
        let mut bins_b = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n {
            bins_a[i] = (freqs_a[i]) * bins[i] / (p_a);
            bins_b[i] = (1.0 - freqs_a[i]) * bins[i] / (1.0 - p_a);
        }
        // Percentiles (cummulative bins summing up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
        let mut percs_a = bins_a.clone();
        let mut percs_b = bins_b.clone();
        for i in 1..bins_a.nrows() {
            percs_a[i] = bins_a.view((0,0), (i+1,1)).sum();
            percs_b[i] = bins_b.view((0,0), (i+1,1)).sum();
        }

        // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
        let mut percs_a0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        let mut percs_b0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n-1 {
            percs_a0[i+1] = percs_a[i];
            percs_b0[i+1] = percs_b[i];
        }
        // // Define cost function
        // let cost = LeastSquaresBeta { percs_a: percs_a.clone(),
        //                                                 percs_b: percs_b.clone(),
        //                                                 q_prime: q_prime.clone() };
        let cost = MaximumLikelihoodBeta { percs_a: percs_a,
                                                                  percs_b: percs_b,
                                                                  percs_a0: percs_a0,
                                                                  percs_b0: percs_b0 };
        // Define initial parameter vector
        // let init_param: Vec<f64> = vec![1.0, 1.0, 1.0, 1.0];
        let h = 1e-9;
        let init_param: Vec<Vec<f64>> = vec![vec![1.00, 1.00, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                             vec![1.05, 1.00, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                             vec![1.00, 1.05, 1.00, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                             vec![1.00, 1.00, 1.05, 1.00].into_iter().map(|x| x * h).collect::<Vec<f64>>(),
                                             vec![1.00, 1.00, 1.00, 1.05].into_iter().map(|x| x * h).collect::<Vec<f64>>()];
        // // let init_param: Vec<Vec<f64>> = vec![vec![1.00, 1.00, 1.00, 1.00],
        // //                                      vec![-1.00, 1.00, 1.00, 1.00],
        // //                                      vec![-1.00, 10.00, 1.00, 1.00],
        // //                                      vec![-1.00, 10.00, 20.00, 1.00],
        // //                                      vec![-1.00, 10.00, 20.00, 30.00]];
        
        // let init_hessian: Vec<Vec<f64>> = vec![vec![1.0, 0.0, 0.0, 0.0],
        // vec![0.0, 1.0, 0.0, 0.0],
        // vec![0.0, 0.0, 1.0, 0.0],
        // vec![0.0, 0.0, 0.0, 1.0]];
        // set up a line search
        // let linesearch = MoreThuenteLineSearch::new();
        // .with_width_tolerance(1e-4)
        // .unwrap();
        // let linesearch = HagerZhangLineSearch::new();
        // Set up solver
        let solver = NelderMead::new(init_param);
        // .with_alpha(1.0).unwrap()
        // .with_gamma(2.0).unwrap()
        // .with_rho(0.5).unwrap()
        // .with_sigma(0.5).unwrap()
        // .with_sd_tolerance(1e-4)
        // .unwrap();
        // // let solver = BFGS::new(linesearch);
        // let solver = LBFGS::new(linesearch, 4);

        // Run solver
        let res = match Executor::new(cost, solver)
            .configure(|state| {
                state
                    // .param(init_param)
                    // .inv_hessian(init_hessian)
                    .max_iters(100)
            })
            // .add_observer(SlogLogger::term(), ObserverMode::Never)
            .run() {
                Ok(x) => x,
                Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
            };
        // let cost_function = cost;

        // let solver = ParticleSwarm::new((vec![0.00, 1e-8, 1e-8, 1e-8],
        //                                                                       vec![100.0, 100.0, 100.0, 100.0]), 
        //                                                               100);

        // let res = match Executor::new(cost_function, solver)
        //     .configure(|state| state.max_iters(1_000))
        //     .add_observer(SlogLogger::term(), ObserverMode::Never)
        //     .run() {
        //         Ok(x) => x,
        //         Err(_) => return None,
        //     };
        // Print result
        // let mut solution = res.state().get_best_param().unwrap().clone().position;
        let params = res.state().param.clone().unwrap();
        let solution = search_params_to_shapes(&params);
        // for i in 0..solution.len() {
        //     if solution[i] <= 0.0 {
        //         solution[i] = 0.00001
        //     }
        // }

        // let (solution, fx) = minimize_unbounded(
        //     |args| maximum_likelihood_beta(args, &percs_a, &percs_b, &percs_a0, &percs_b0),
        //     vec![10.0,10.0,10.0,10.0],
        //     1.0,
        //     Params::default(),
        //     1000);
         

        println!("LOCUS: {:?}", first_2_col);
        println!("SOL={:?}", solution);
        let a_mu_hat = min + (max-min) * (solution[0]/(solution[0]+solution[1]));
        let b_mu_hat = min + (max-min) * (solution[2]/(solution[2]+solution[3]));
        // let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(x_matrix.column(j).mean().to_string());
        line.push("Pheno_0".to_string());
        line.push(alpha.to_string() + &",Unknown\n");

    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}
