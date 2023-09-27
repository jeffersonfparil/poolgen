use crate::base::*;
use crate::popgen::*;
use ndarray::{prelude::*, Zip};
use statrs::distribution::Continuous;
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};
use argmin::core::{self, CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use statrs::distribution::Normal;

const PARAMETER_LOWER_LIMIT: f64 = f64::EPSILON;
const PARAMETER_UPPER_LIMIT: f64 = 1e24;

fn maximum_likelihood_normal(
    params: &Vec<f64>,
    q: &Array1<f64>,
) -> f64 {
    let sigma = bound_parameters_with_logit(&vec![params[1]], PARAMETER_LOWER_LIMIT, PARAMETER_UPPER_LIMIT)[0];
    // println!("shapes={:?}", shapes);
    let distribution = Normal::new(params[0], sigma).expect(
        ("mu=".to_owned() + &params[0].to_string()[..] + "; sigma=" + &sigma.to_string()[..]).as_str(),
    );
    q.iter()
        .map(|&x| -1.00 * distribution.ln_pdf(x))
        .fold(0.0, |sum, x| sum + x)
}

impl CostFunction for MaximumLikelihoodNormal {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(maximum_likelihood_normal(
            &p,
            &self.q,
        ))
    }
}

fn ml_normal_1d(
    solver: NelderMead<Vec<f64>, f64>,
    q: Array1<f64>,
) -> Option<Vec<f64>> {
    let cost = MaximumLikelihoodNormal {
        q: q,
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

/// gudmc: **g**enomewide **u**nbiased **d**etermination of the **m**odes of **c**onvergent evolution
/// using Tajima's D to detect peaks (balancing selection) and troughs (directional selection) as locations
pub fn gudmc(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    recombination_rate: &f64,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<i32> {
    // Calculate Tajima's D
    let fname_tajima = tajima_d(
        genotypes_and_phenotypes,
        pool_sizes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
        fname_input,
        &"gudmc_intermediate_file_tajimasD.tmp".to_owned(),
    )
    .unwrap();
    let (tajima_row_labels, tajima_col_labels, tajima) = load_sliding_window_tables(&fname_tajima, 
        &",".to_owned(), // rows are the populations, and columns are the windows
        &vec![0], // population ID
        &2, // skip the population ID and mean Tajima's D across the whole genome
        &10_000_000_000_000_000_000).unwrap();
        println!("tajima_row_labels={:?}", tajima_row_labels);
        println!("tajima_col_labels={:?}", tajima_col_labels);
        println!("tajima={:?}", tajima);
    // Calculate pairwise Fst (all pairwise combinations)
    let (_, fname_fst) = fst(
        genotypes_and_phenotypes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
        fname_input,
        &"gudmc_intermediate_file_Fst.tmp".to_owned(),
    )
    .unwrap();
    let (fst_row_labels, fst_col_labels, fst) = load_sliding_window_tables(&fname_fst, 
        &",".to_owned(), // rows are the windows, and columns are population pairs
        &vec![0,1,2], // chr, pos_ini, and pos_fin
        &3, // start of pairwise Fst per window
        &10_000_000_000_000_000_000).unwrap(); // set the end column to a very large number so we default to the last column
        println!("fst_row_labels={:?}", fst_row_labels);
        println!("fst_col_labels={:?}", fst_col_labels);
        println!("fst={:?}", fst);
    // Sanity checks
    let n = tajima.len();       // number of populations
    let w = tajima[0].len();    // number of windows
    let nxn = fst[0].len();     // number of population pairs
    let w_ = fst.len();         // number of windows according to Fst matrix - compare with that of Tajima's D matrix
    assert!(n*n == nxn, "Tajima's D and Fst calculations are not matching.");
    assert!(w == w_, "Tajima's D and Fst calculations are not matching.");
    // PER POPULATION: find significant troughs (selective sweeps) and peaks (balancing selection) and meausre their widths
    let mut tajima_pop: Vec<String> = vec![];
    let mut tajima_chr: Vec<String> = vec![];
    let mut tajima_pos_ini: Vec<u64> = vec![];
    let mut tajima_pos_fin: Vec<u64> = vec![];
    let mut tajima_d: Vec<f64> = vec![];
    let mut tajima_width: Vec<f64> = vec![];
    for i in 0..n {
        let d = Array1::from_vec(tajima[i].clone());
        println!("d={:?}", d);
        // Fit a normal distribution to d and find the location of d_i at p-value <= 5% (two-tailed)
        let solver = prepare_solver_neldermead(2.0, 1.0);
        let solution = match ml_normal_1d(solver, d) {
            Some(x) => x,
            None => vec![f64::NAN, f64::NAN],
        };
        println!("solution={:?}", solution);
    }

    // PER PAIR OF POPULATIONS: find the significant deviations in Fst within the above-identified Tajima's D peaks and troughs

    Ok(0)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_gudmc() {
        let x: Array2<f64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.5, 0.5, 0.0,
                0.5, 0.5, 1.0, 0.7, 0.2, 0.1, 0.7, 0.3, 1.0, 0.7, 0.2, 0.1, 0.7, 0.3,
            ],
        )
        .unwrap();
        let y: Array2<f64> = Array2::from_shape_vec(
            (5, 2),
            vec![2.0, 0.5, 1.0, 0.2, 2.0, 0.5, 4.0, 0.0, 5.0, 0.5],
        )
        .unwrap();
        let genotypes_and_phenotypes = GenotypesAndPhenotypes {
            chromosome: vec![
                "Intercept".to_owned(),
                "X".to_owned(),
                "X".to_owned(),
                "X".to_owned(),
                "Y".to_owned(),
                "Y".to_owned(),
            ],
            position: vec![0, 123, 123, 123, 456, 456],
            allele: vec![
                "Intercept".to_owned(),
                "a".to_string(),
                "g".to_string(),
                "d".to_string(),
                "c".to_string(),
                "t".to_string(),
            ],
            intercept_and_allele_frequencies: x.clone(),
            phenotypes: y.clone(),
            pool_names: vec![
                "Pop1".to_owned(),
                "Pop2".to_owned(),
                "Pop3".to_owned(),
                "Pop4".to_owned(),
                "Pop5".to_owned(),
            ],
            coverages: Array2::from_shape_vec(
                (5, 2),
                vec![
                    10.0, 10.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0,
                ],
            )
            .unwrap(),
        };
        // Outputs
        let out = gudmc(
            &genotypes_and_phenotypes,
            &vec![42.0, 42.0, 42.0, 42.0, 42.0],
            &3.14e-6,
            &100,
            &50,
            &1,
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();

        // Assertions
        assert_eq!(out, 0);
        // assert_eq!(0, 1);
    }
}
