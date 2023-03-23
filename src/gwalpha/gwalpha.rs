// use std::io::{self, Error, ErrorKind};
use nalgebra::{self, DMatrix};

use argmin::core::observers::{ObserverMode, SlogLogger};
use argmin::core::{self, CostFunction, Gradient, Executor};
use finitediff::FiniteDiff;
// use argmin::solver::neldermead::NelderMead;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::BFGS;
use argmin::solver::quasinewton::LBFGS;
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

fn least_squares_beta(params: &Vec<f64>, percs_a: &DMatrix<f64>, percs_b: &DMatrix<f64>, q_prime: &DMatrix<f64>) -> f64 {
    let (mut a_shape1, mut a_shape2) = (params[0], params[1]);
    let (mut b_shape1, mut b_shape2) = (params[2], params[3]);
    // println!("a1={}; a2={}; b1={}; b2={}", a_shape1, a_shape2, b_shape1, b_shape2);
    if a_shape1 <= 0.0 {
        a_shape1 = a_shape1.abs() + 1e-9;
    }
    if a_shape2 <= 0.0 {
        a_shape2 = a_shape2.abs() + 1e-9;
    }
    if b_shape1 <= 0.0 {
        b_shape1 = b_shape1.abs() + 1e-9;
    }
    if b_shape2 <= 0.0 {
        b_shape2 = b_shape2.abs() + 1e-9;
    }
    let a_dist = match Beta::new(a_shape1, a_shape2) {
        Ok(x) => x,
        Err(_) => return f64::INFINITY
    };
    let b_dist = match Beta::new(b_shape1, b_shape2) {
        Ok(x) => x,
        Err(_) => return f64::INFINITY
    };
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
    let (mut a_shape1, mut a_shape2) = (params[0], params[1]);
    let (mut b_shape1, mut b_shape2) = (params[2], params[3]);
    if a_shape1 <= 0.0 {
        a_shape1 = a_shape1.abs() + 1e-9;
    }
    if a_shape2 <= 0.0 {
        a_shape2 = a_shape2.abs() + 1e-9;
    }
    if b_shape1 <= 0.0 {
        b_shape1 = b_shape1.abs() + 1e-9;
    }
    if b_shape2 <= 0.0 {
        b_shape2 = b_shape2.abs() + 1e-9;
    }
    let a_dist = match Beta::new(a_shape1, a_shape2) {
        Ok(x) => x,
        Err(_) => return f64::INFINITY
    };
    let b_dist = match Beta::new(b_shape1, b_shape2) {
        Ok(x) => x,
        Err(_) => return f64::INFINITY
    };
    // let a_dist = Beta::new(a_shape1, a_shape2).unwrap();
    // let b_dist = Beta::new(b_shape1, b_shape2).unwrap();
    // println!("percs_a={:?}", a_dist.cdf(percs_a[2]));
    // println!("percs_a0={:?}", percs_a0);
    // println!("percs_b={:?}", percs_b);
    // println!("percs_b0={:?}", percs_b0);
    let n = percs_a.len();
    let mut a_log_likelihood: f64 = 0.0;
    let mut b_log_likelihood: f64 = 0.0;
    for i in 0..n {
        a_log_likelihood += f64::log10(a_dist.cdf(percs_a[i]) - a_dist.cdf(percs_a0[i]));
        b_log_likelihood += f64::log10(b_dist.cdf(percs_b[i]) - b_dist.cdf(percs_b0[i]));
    }
    let out = -a_log_likelihood - b_log_likelihood;
    println!("a1={}; a2={}; b1={}; b2={}; a_log_likelihood={}; b_log_likelihood={}; out={}", a_shape1, a_shape2, b_shape1, b_shape2, a_log_likelihood, b_log_likelihood, out);
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

pub fn gwalpha(locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes, filter_stats: &FilterStats) -> Option<String> {
    // println!("locus_counts_and_phenotypes={:?}", locus_counts_and_phenotypes);
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
    let x_matrix = (*locus_frequencies).matrix.clone();
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
    // println!("x_matrix={:?}", x_matrix);
    // println!("bins={:?}", bins);
    // println!("q={:?}", q);
    // println!("sig={:?}", sig);
    // println!("min={:?}", min);
    // println!("max={:?}", max);
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
        let freqs_a: DMatrix<f64> = DMatrix::from_columns(&[x_matrix.column(j)]);
        // println!("bins={:?}", bins);
        // println!("freqs_a={:?}", freqs_a);
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
        let mut percs_a = DMatrix::from_element(bins_a.nrows(), 1, 0.0);
        let mut percs_b = DMatrix::from_element(bins_b.nrows(), 1, 0.0);
        for i in 1..bins_a.nrows() {
            for k in 0..i {
                percs_a[i] += bins_a[k];
                percs_b[i] += bins_b[k];
            }
        }
        // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
        let mut percs_a0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        let mut percs_b0: DMatrix<f64> = DMatrix::from_element(n, 1, 0.0);
        for i in 0..n-1 {
            percs_a0[i+1] = percs_a[i];
            percs_b0[i+1] = percs_b[i];
        }
        // Test
        println!("min={:?}", min);
        println!("max={:?}", max);
        println!("sig={:?}", sig);
        println!("p_a={:?}", p_a);
        println!("bins={:?}", bins);
        println!("q={:?}", q);
        println!("bins_a={:?}", bins_a);
        println!("bins_b={:?}", bins_b);
        println!("freqs_a={:?}", freqs_a);
        println!("q_prime={:?}", q_prime);
        println!("percs_a={:?}", percs_a);
        println!("percs_b={:?}", percs_b);
        println!("percs_a0={:?}", percs_a0);
        println!("percs_b0={:?}", percs_b0);
        // Define cost function
        // let cost = LeastSquaresBeta { percs_a: percs_a.clone(),
        //                                                 percs_b: percs_b.clone(),
        //                                                 q_prime: q_prime.clone() };
        let cost = MaximumLikelihoodBeta { percs_a: percs_a,
                                                                  percs_b: percs_b,
                                                                  percs_a0: percs_a0,
                                                                  percs_b0: percs_b0 };
        // Define initial parameter vector
        let init_param: Vec<f64> = vec![1.0, 1.0, 1.0, 1.0];
        let init_hessian: Vec<Vec<f64>> = vec![vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0],
        vec![0.0, 0.0, 1.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0]];
        // set up a line search
        let linesearch = MoreThuenteLineSearch::new()
        .with_bounds(1e-5, 100.0).unwrap();

        // Set up solver
        // let solver = BFGS::new(linesearch);
        let solver = LBFGS::new(linesearch, 7);

        // Run solver
        let res = match Executor::new(cost, solver)
            .configure(|state| {
                state
                    .param(init_param)
                    // .inv_hessian(init_hessian)
                    .max_iters(50)
            })
            // .add_observer(SlogLogger::term(), ObserverMode::Always)
            .add_observer(SlogLogger::term(), ObserverMode::Never)
            .run() {
                Ok(x) => x,
                Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
            };
        // Print result
        // println!("{res}");
        let mut solution = res.state().param.clone().unwrap();
        for i in 0..solution.len() {
            if solution[i] <= 0.0 {
                solution[i] = 0.00001
            }
        }
        let a_mu_hat = min + (max-min) * (solution[0]/(solution[0]+solution[1]));
        let b_mu_hat = min + (max-min) * (solution[2]/(solution[2]+solution[3]));
        let alpha = (2.00*f64::sqrt(p_a*(1.0-p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push("Pheno_0".to_string());
        line.push(alpha.to_string() + &",Unknown\n");

    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}
