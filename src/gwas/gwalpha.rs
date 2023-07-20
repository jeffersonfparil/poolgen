use crate::base::*;
use argmin::core::{self, CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use ndarray::prelude::*;

use statrs::distribution::{Beta, ContinuousCDF};

const PARAMETER_LOWER_LIMIT: f64 = f64::EPSILON;
const PARAMETER_UPPER_LIMIT: f64 = 10.00;

fn least_squares_beta(
    params: &Vec<f64>,
    percs_a: &Array1<f64>,
    percs_b: &Array1<f64>,
    q_prime: &Array1<f64>,
) -> f64 {
    let shapes = bound_parameters_with_logit(params, PARAMETER_LOWER_LIMIT, PARAMETER_UPPER_LIMIT);
    // println!("shapes={:?}", shapes);
    let a_dist = Beta::new(shapes[0], shapes[1]).expect(
        &shapes
            .clone()
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("-")[..],
    );
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

fn maximum_likelihood_beta(
    params: &Vec<f64>,
    percs_a: &Array1<f64>,
    percs_b: &Array1<f64>,
    percs_a0: &Array1<f64>,
    percs_b0: &Array1<f64>,
) -> f64 {
    let shapes = bound_parameters_with_logit(params, PARAMETER_LOWER_LIMIT, PARAMETER_UPPER_LIMIT);
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
            false => diff_a,
        };
        diff_b = match diff_b < f64::EPSILON {
            true => f64::EPSILON,
            false => diff_b,
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
        Ok(least_squares_beta(
            &p,
            &self.percs_a,
            &self.percs_b,
            &self.q_prime,
        ))
    }
}

impl CostFunction for MaximumLikelihoodBeta {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(maximum_likelihood_beta(
            &p,
            &self.percs_a,
            &self.percs_b,
            &self.percs_a0,
            &self.percs_b0,
        ))
    }
}

fn gwalpha_minimise_ls(
    solver: NelderMead<Vec<f64>, f64>,
    q_prime: Array1<f64>,
    percs_a: Array1<f64>,
    percs_b: Array1<f64>,
) -> Option<Vec<f64>> {
    let cost = LeastSquaresBeta {
        q_prime: q_prime,
        percs_a: percs_a,
        percs_b: percs_b,
    };
    let res = match Executor::new(cost, solver)
        .configure(|state| state.max_iters(1_000))
        // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
        .run()
    {
        Ok(x) => x,
        Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
    };
    // println!("CONVERGENCE: {:?}", res.state());
    let params = res.state().param.clone().unwrap();
    let solution =
        bound_parameters_with_logit(&params, PARAMETER_LOWER_LIMIT, PARAMETER_UPPER_LIMIT);
    Some(solution)
}

fn gwalpha_minimise_ml(
    solver: NelderMead<Vec<f64>, f64>,
    percs_a: Array1<f64>,
    percs_a0: Array1<f64>,
    percs_b: Array1<f64>,
    percs_b0: Array1<f64>,
) -> Option<Vec<f64>> {
    let cost = MaximumLikelihoodBeta {
        percs_a: percs_a,
        percs_b: percs_b,
        percs_a0: percs_a0,
        percs_b0: percs_b0,
    };
    let res = match Executor::new(cost, solver)
        .configure(|state| state.max_iters(1_000))
        // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
        .run()
    {
        Ok(x) => x,
        Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
    };
    // println!("CONVERGENCE: {:?}", res.state());
    let params = res.state().param.clone().unwrap();
    let solution =
        bound_parameters_with_logit(&params, PARAMETER_LOWER_LIMIT, PARAMETER_UPPER_LIMIT);
    Some(solution)
}

fn prepare_geno_and_pheno_stats(
    locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes,
    filter_stats: &FilterStats,
) -> Option<(
    LocusFrequencies,
    Array1<f64>,
    Array1<f64>,
    f64,
    f64,
    f64,
    usize,
    usize,
)> {
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
    // Keep p-1 alleles if p >= 2 so we have degrees of freedom to fit the intercept
    if locus_frequencies.matrix.ncols() >= 2 {
        locus_frequencies.matrix.remove_index(Axis(1), 0);
        locus_frequencies.alleles_vector.remove(0);
    }
    // Extract phenotype information (Note: removing NEG_INFINITY from the bins and q columns if we have less than 3 pools, each corresponds to sig, MIN, and MAX rows with 3 as the minimum number of rows)
    let bins_tmp = locus_counts_and_phenotypes
        .phenotypes
        .column(0)
        .into_iter()
        .filter(|x| **x != f64::NEG_INFINITY)
        .map(|x| x.to_owned())
        .collect::<Vec<f64>>();
    let q_tmp = locus_counts_and_phenotypes
        .phenotypes
        .column(1)
        .into_iter()
        .filter(|x| **x != f64::NEG_INFINITY)
        .map(|x| x.to_owned())
        .collect::<Vec<f64>>();
    let bins = Array1::from_vec(bins_tmp);
    let q = Array1::from_vec(q_tmp);
    let sig = locus_counts_and_phenotypes.phenotypes[(0, 2)];
    let min = locus_counts_and_phenotypes.phenotypes[(1, 2)];
    let max = locus_counts_and_phenotypes.phenotypes[(2, 2)];
    // Check if we have a compatible allele frequency and phenotype matrix or vector
    let n = locus_frequencies.matrix.nrows();
    let p = locus_frequencies.matrix.ncols();
    let m = bins.len();
    if n != m {
        return None;
    }
    Some((*locus_frequencies.clone(), bins, q, sig, min, max, n, p))
}

fn prepare_freqs_and_qprime(
    locus_frequencies: &LocusFrequencies,
    bins: &Array1<f64>,
    q: &Array1<f64>,
    min: f64,
    max: f64,
    n: usize,
    j: usize,
) -> (
    f64,
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
    Array1<f64>,
) {
    let freqs_a: ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>> =
        locus_frequencies.matrix.column(j);
    let p_a = freqs_a.t().dot(bins); // mean allele frequency across pools
                                     // println!("p_a={:?}", p_a);
                                     // Quantiles per pool (for least squares estimation)
    let mut q_prime: Array1<f64> = Array1::zeros(n);
    for i in 1..n {
        q_prime[i] = (q[i] - min) / (max - min);
    }
    // println!("q_prime={:?}", q_prime);
    // Bins (sums up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
    let mut bins_a = Array1::zeros(n);
    let mut bins_b = Array1::zeros(n);
    for i in 0..n {
        bins_a[i] = (freqs_a[i]) * bins[i] / (p_a);
        bins_b[i] = (1.0 - freqs_a[i]) * bins[i] / (1.0 - p_a);
    }
    // println!("bins_a={:?}", bins_a);
    // println!("bins_b={:?}", bins_b);
    // Percentiles (cummulative bins summing up to 1.0) of the current allele and its additive inverse representing the rest of the alleles
    let mut percs_a: Array1<f64> = bins_a.clone();
    let mut percs_b: Array1<f64> = bins_b.clone();
    for i in 1..bins_a.len() {
        percs_a[i] = bins_a.slice(s![0..(i + 1)]).sum();
        percs_b[i] = bins_b.slice(s![0..(i + 1)]).sum();
    }
    // println!("percs_a={:?}", percs_a);
    // println!("percs_b={:?}", percs_b);
    // Percentiles of the current allele and its additive inverse for modelling their distrbutions across pools
    let mut percs_a0: Array1<f64> = Array1::zeros(n);
    let mut percs_b0: Array1<f64> = Array1::zeros(n);
    for i in 0..n - 1 {
        percs_a0[i + 1] = percs_a[i];
        percs_b0[i + 1] = percs_b[i];
    }
    (p_a, q_prime, percs_a, percs_a0, percs_b, percs_b0)
}

pub fn gwalpha_ls(
    locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes,
    filter_stats: &FilterStats,
) -> Option<String> {
    // Check struct
    locus_counts_and_phenotypes.check().unwrap();
    // Prepare phenotype and genotype statistics
    let (locus_frequencies, bins, q, sig, min, max, n, p) =
        match prepare_geno_and_pheno_stats(locus_counts_and_phenotypes, filter_stats) {
            Some(x) => x,
            None => return None,
        };
    // Prepare output line
    let first_2_col = vec![
        locus_frequencies.chromosome.to_owned(),
        locus_frequencies.position.to_string(),
    ];
    let mut line: Vec<String> = vec![];
    // Instiate input optimisation data (NOTE: exluding percs_a0 and percs_b0)
    let (mut p_a, mut q_prime, mut percs_a, mut percs_b);
    // Instiate solver
    let mut solver: NelderMead<Vec<f64>, f64>;
    // Iterate across alleles
    for j in 0..p {
        // Prepare allele frequecies, quantiles and percentiles
        (p_a, q_prime, percs_a, _, percs_b, _) =
            prepare_freqs_and_qprime(&locus_frequencies, &bins, &q, min, max, n, j);
        // Optimise
        solver = prepare_solver_neldermead(4.0, 1.0);
        let solution = match gwalpha_minimise_ls(solver, q_prime, percs_a, percs_b) {
            Some(x) => x,
            None => return None,
        };
        let a_mu_hat = min + (max - min) * (solution[0] / (solution[0] + solution[1]));
        let b_mu_hat = min + (max - min) * (solution[2] / (solution[2] + solution[3]));
        let alpha = (2.00 * f64::sqrt(p_a * (1.0 - p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(parse_f64_roundup_and_own(
            locus_frequencies.matrix.column(j).mean().unwrap(),
            6,
        ));
        line.push("Pheno_0".to_string());
        line.push(parse_f64_roundup_and_own(alpha, 6) + &",Unknown\n");
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

pub fn gwalpha_ml(
    locus_counts_and_phenotypes: &mut LocusCountsAndPhenotypes,
    filter_stats: &FilterStats,
) -> Option<String> {
    // Check struct
    locus_counts_and_phenotypes.check().unwrap();
    // Prepare phenotype and genotype statistics
    let (locus_frequencies, bins, q, sig, min, max, n, p) =
        match prepare_geno_and_pheno_stats(locus_counts_and_phenotypes, filter_stats) {
            Some(x) => x,
            None => return None,
        };
    // Prepare output line
    let first_2_col = vec![
        locus_frequencies.chromosome.to_owned(),
        locus_frequencies.position.to_string(),
    ];
    let mut line: Vec<String> = vec![];
    // Instiate input optimisation data (NOTE: exluding q_prime)
    let (mut p_a, mut percs_a, mut percs_a0, mut percs_b, mut percs_b0);
    // Instiate solver
    let mut solver: NelderMead<Vec<f64>, f64>;
    // Iterate across alleles
    for j in 0..p {
        // Prepare allele frequecies, quantiles and percentiles
        (p_a, _, percs_a, percs_a0, percs_b, percs_b0) =
            prepare_freqs_and_qprime(&locus_frequencies, &bins, &q, min, max, n, j);
        // Optimise
        solver = prepare_solver_neldermead(4.0, 1.0);
        let solution = match gwalpha_minimise_ml(solver, percs_a, percs_a0, percs_b, percs_b0) {
            Some(x) => x,
            None => return None,
        };
        let a_mu_hat = min + (max - min) * (solution[0] / (solution[0] + solution[1]));
        let b_mu_hat = min + (max - min) * (solution[2] / (solution[2] + solution[3]));
        let alpha = (2.00 * f64::sqrt(p_a * (1.0 - p_a))) * (a_mu_hat - b_mu_hat) / sig;
        // let alpha =  (a_mu_hat - b_mu_hat) / sig;
        // Fill in output line
        line.append(&mut first_2_col.clone());
        line.push(locus_frequencies.alleles_vector[j].clone());
        line.push(parse_f64_roundup_and_own(
            locus_frequencies.matrix.column(j).mean().unwrap(),
            6,
        ));
        line.push("Pheno_0".to_string());
        line.push(parse_f64_roundup_and_own(alpha, 6) + &",Unknown\n");
    }
    let out = line.join(",").replace("\n,", "\n");
    Some(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_gwalpha() {
        // Expected
        let expected_output1: String = "Chromosome1,12345,A,0.353287,Pheno_0,5.816067,Unknown\nChromosome1,12345,T,0.267133,Pheno_0,9.176892,Unknown\n".to_owned();
        let expected_output2: String = "Chromosome1,12345,A,0.353287,Pheno_0,-3.293261,Unknown\nChromosome1,12345,T,0.267133,Pheno_0,-7.098985,Unknown\n".to_owned();
        // Inputs
        let counts: Array2<u64> =
            Array2::from_shape_vec((5, 3), vec![5, 2, 6, 2, 2, 7, 3, 2, 5, 4, 3, 3, 5, 5, 0])
                .unwrap();
        let gwalpha_fmt: Array2<f64> = Array2::from_shape_vec(
            (3, 5),
            vec![
                0.2,
                0.2,
                0.2,
                0.2,
                0.2,
                0.0,
                0.1,
                0.4,
                0.7,
                0.9,
                0.02,
                0.0,
                0.9,
                f64::NEG_INFINITY,
                f64::NEG_INFINITY,
            ],
        )
        .unwrap()
        .reversed_axes();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            pool_sizes: vec![20.0, 20.0, 20.0, 20.0, 20.0],
        };
        let locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            alleles_vector: vec!["A".to_owned(), "T".to_owned(), "D".to_owned()],
            matrix: counts,
        };
        let mut locus_counts_and_phenotypes = LocusCountsAndPhenotypes {
            locus_counts: locus_counts,
            phenotypes: gwalpha_fmt.clone(),
            pool_names: vec!["pool1", "pool2", "pool3", "pool4", "pool5"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
        };
        // Outputs
        let ls_line = gwalpha_ls(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        let ml_line = gwalpha_ml(&mut locus_counts_and_phenotypes, &filter_stats).unwrap();
        // Assertions
        assert_eq!(expected_output1, ls_line);
        assert_eq!(expected_output2, ml_line);
    }
}
