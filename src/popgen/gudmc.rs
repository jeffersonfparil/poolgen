use crate::base::*;
use crate::popgen::*;
use argmin::core::{self, CostFunction, Executor};
use argmin::solver::neldermead::NelderMead;
use ndarray::prelude::*;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

use std::fs;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

const PARAMETER_LOWER_LIMIT: f64 = f64::EPSILON;
const PARAMETER_UPPER_LIMIT: f64 = 1e24;

fn maximum_likelihood_normal(params: &Vec<f64>, q: &Array1<f64>) -> f64 {
    let sigma = bound_parameters_with_logit(
        &vec![params[1]],
        PARAMETER_LOWER_LIMIT,
        PARAMETER_UPPER_LIMIT,
    )[0];
    // println!("shapes={:?}", shapes);
    let distribution = Normal::new(params[0], sigma).expect(
        ("mu=".to_owned() + &params[0].to_string()[..] + "; sigma=" + &sigma.to_string()[..])
            .as_str(),
    );
    q.iter()
        .map(|&x| -1.00 * distribution.ln_pdf(x))
        .fold(0.0, |sum, x| sum + x)
}

impl CostFunction for MaximumLikelihoodNormal {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Self::Param) -> Result<Self::Output, core::Error> {
        Ok(maximum_likelihood_normal(&p, &self.q))
    }
}

fn ml_normal_1d(solver: NelderMead<Vec<f64>, f64>, q: &Array1<f64>) -> Option<Vec<f64>> {
    let cost = MaximumLikelihoodNormal { q: q.clone() };
    let res = match Executor::new(cost, solver)
        .configure(|state| state.max_iters(10_000))
        // .add_observer(SlogLogger::term(), ObserverMode::NewBest)
        .run()
    {
        Ok(x) => x,
        Err(_) => return None, // Error occurs when the optimiser MoreThuenteLineSearch moves in the wrong direction
    };
    // println!("CONVERGENCE: {:?}", res.state());
    let params = res.state().param.clone().unwrap();
    let solution = vec![
        params[0],
        bound_parameters_with_logit(
            &vec![params[1]],
            PARAMETER_LOWER_LIMIT,
            PARAMETER_UPPER_LIMIT,
        )[0],
    ];
    Some(solution)
}

/// gudmc: **g**enomewide **u**nbiased **d**etermination of the **m**odes of **c**onvergent evolution
/// using Tajima's D to detect peaks (balancing selection) and troughs (directional selection) as locations
pub fn gudmc(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    sigma_threshold: &f64,
    recombination_rate_cM_per_Mb: &f64,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    /////////////////////////////////////////////////////////
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
    let (tajima_row_labels, tajima_col_labels, tajima) = load_table(
        &fname_tajima,
        &",".to_owned(), // rows are the populations, and columns are the windows
        &vec![0],        // population ID
        &2,              // skip the population ID and mean Tajima's D across the whole genome
        &(PARAMETER_UPPER_LIMIT as usize),
    )
    .unwrap();
    // println!("tajima_col_labels={:?}", tajima_col_labels);
    // println!("tajima_col_labels.len()={:?}", tajima_col_labels.len());
    // println!("tajima_row_labels={:?}", tajima_row_labels);
    // println!("tajima={:?}", tajima);
    /////////////////////////////////////////////////////////
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
    let (fst_row_labels, fst_col_labels, fst) = load_table(
        &fname_fst,
        &",".to_owned(), // rows are the windows, and columns are population pairs
        &vec![0, 1, 2],  // chr, pos_ini, and pos_fin
        &3,              // start of pairwise Fst per window
        &(PARAMETER_UPPER_LIMIT as usize),
    )
    .unwrap();
    // set the end column to a very large number so we default to the last column
    // println!("fst_row_labels={:?}", fst_row_labels);
    // println!("fst_row_labels.len()={:?}", fst_row_labels.len());
    // println!("fst_col_labels={:?}", fst_col_labels);
    // println!("fst={:?}", fst);
    // Sanity checks
    let n = tajima.len(); // number of populations
    let w = tajima[0].len(); // number of windows
    let nxn = fst[0].len(); // number of population pairs
    let w_ = fst.len(); // number of windows according to Fst matrix - compare with that of Tajima's D matrix
                        // println!("n={}; w={}; nxn={}; w_={}", n, w, nxn, w_);
    assert!(
        n * n == nxn,
        "Tajima's D and Fst calculations are not matching."
    );
    assert!(w == w_, "Tajima's D and Fst calculations are not matching.");
    /////////////////////////////////////////////////////////
    // PER POPULATION: find significant troughs (selective sweeps) and peaks (balancing selection) and measure their widths
    let mut tajima_pop: Vec<String> = vec![];
    let mut tajima_chr: Vec<Vec<String>> = vec![]; // each sub-vector represents a population
    let mut tajima_pos_ini: Vec<Vec<u64>> = vec![]; // each sub-vector represents a population
    let mut tajima_pos_fin: Vec<Vec<u64>> = vec![]; // each sub-vector represents a population
    let mut tajima_d: Vec<Vec<f64>> = vec![]; // each sub-vector represents a population
    let mut tajima_d_mean: Vec<Vec<f64>> = vec![]; // each sub-vector represents a population
    let mut tajima_d_sd: Vec<Vec<f64>> = vec![]; // each sub-vector represents a population
    let mut tajima_width: Vec<Vec<u64>> = vec![]; // each sub-vector represents a population
    for i in 0..n {
        // Instatntiate the information for the current population
        tajima_pop.push(tajima_row_labels[i].clone());
        tajima_chr.push(vec![]);
        tajima_pos_ini.push(vec![]);
        tajima_pos_fin.push(vec![]);
        tajima_d.push(vec![]);
        tajima_d_mean.push(vec![]);
        tajima_d_sd.push(vec![]);
        tajima_width.push(vec![]);
        let d = Array1::from_vec(
            tajima[i]
                .iter()
                .filter(|&&x| x.is_nan() == false)
                .map(|&x| x.clone())
                .collect::<Vec<f64>>(),
        );
        // Fit a normal distribution to d
        let solver = prepare_solver_neldermead(2.0, 1.0);
        let solution = match ml_normal_1d(solver, &d) {
            Some(x) => x,
            None => vec![f64::NAN, f64::NAN],
        };
        // Find troughs and peaks above the `sigma_threshold`
        for j in 0..d.len() {
            let window_id = tajima_col_labels[j].split("-").collect::<Vec<&str>>()[1]
                .to_owned()
                .split("_")
                .map(|x| x.to_owned())
                .collect::<Vec<String>>();
            tajima_chr[i].push(window_id[0..(window_id.len() - 2)].join("_").to_owned());
            tajima_pos_ini[i].push(window_id[window_id.len() - 2].parse::<u64>().unwrap());
            tajima_pos_fin[i].push(window_id[window_id.len() - 1].parse::<u64>().unwrap());
            tajima_d[i].push(d[j]);
            tajima_d_mean[i].push(solution[0]);
            tajima_d_sd[i].push(solution[1]);
            if (d[j] - solution[0]).abs() >= *sigma_threshold {
                // println!("tajima_col_labels[j]={:?}", tajima_col_labels[j]);
                tajima_width[i]
                    .push(tajima_pos_fin[i].last().unwrap() - tajima_pos_ini[i].last().unwrap());
                // Estimate trough and peak widths
                if tajima_chr[i].len() > 1 {
                    let idx_current = tajima_chr[i].len() - 1;
                    let idx_previous = idx_current - 1;
                    if (tajima_chr[i][idx_current] == tajima_chr[i][idx_previous])
                        & (tajima_pos_ini[i][idx_current] <= tajima_pos_fin[i][idx_previous])
                    {
                        tajima_width[i][idx_current] += tajima_width[i][idx_previous];
                    }
                }
                // // Skip if if have smaller than expected window sizes (Note: we expect at least half the window size in size, we can get smaller window sizes as function of the minimum coverage per window and we are noting window sizes based on the coordinates of the loci covered)
                // if tajima_width[i].last().unwrap() < &((*window_size_bp as f64 / 2.0).ceil() as u64)
                // {
                //     tajima_chr[i].pop();
                //     tajima_pos_ini[i].pop();
                //     tajima_pos_fin[i].pop();
                //     tajima_d[i].pop();
                //     tajima_d_mean[i].pop();
                //     tajima_d_sd[i].pop();
                //     tajima_width[i].pop();
                // }
            } else {
                // Non-significant troughs and peaks
                tajima_width[i].push(0);
            }
        }
    }
    // println!("tajima_row_labels={:?}", tajima_row_labels);
    // println!("tajima_row_labels.len()={:?}", tajima_row_labels.len());
    // println!("tajima_pop={:?}", tajima_pop);
    // println!("tajima_chr={:?}", tajima_chr);
    // println!("tajima_pos_ini={:?}", tajima_pos_ini);
    // println!("tajima_pos_fin={:?}", tajima_pos_fin);
    // println!("tajima_d={:?}", tajima_d);
    // println!("tajima_width={:?}", tajima_width);
    /////////////////////////////////////////////////////////
    // Extract pairwise Fst per window per population pair
    let mut fst_pop_a: Vec<String> = vec![];
    let mut fst_pop_b: Vec<String> = vec![];
    let mut fst_chr: Vec<Vec<String>> = vec![]; // each sub-vector represents a population
    let mut fst_pos_ini: Vec<Vec<u64>> = vec![]; // each sub-vector represents a population
    let mut fst_pos_fin: Vec<Vec<u64>> = vec![]; // each sub-vector represents a population
    let mut fst_f: Vec<Vec<f64>> = vec![]; // each sub-vector represents a population
    let mut fst_f_mean: Vec<f64> = vec![];
    let mut fst_f_sd: Vec<f64> = vec![];
    // let mut test: Vec<f64> = vec![];
    for j in 0..fst_col_labels.len() {
        let pops = fst_col_labels[j].split("_vs_").collect::<Vec<&str>>();
        fst_pop_a.push(pops[0].to_owned());
        fst_pop_b.push(pops[1].to_owned());
        fst_chr.push(vec![]);
        fst_pos_ini.push(vec![]);
        fst_pos_fin.push(vec![]);
        fst_f.push(vec![]);
        for i in 0..fst_row_labels.len() {
            let window = fst_row_labels[i].split("__-__").collect::<Vec<&str>>();
            fst_chr[j].push(window[0].to_owned());
            fst_pos_ini[j].push(window[1].parse::<u64>().unwrap());
            fst_pos_fin[j].push(window[2].parse::<u64>().unwrap());
            fst_f[j].push(fst[i][j]);
        }
        let f: Array1<f64> = Array1::from_vec(
            fst_f[j]
                .clone()
                .into_iter()
                .filter(|&x| x.is_nan() == false)
                .collect::<Vec<f64>>(),
        );
        // Fit a normal distribution to fst
        let solver = prepare_solver_neldermead(2.0, 1.0);
        let solution = match ml_normal_1d(solver, &f) {
            Some(x) => x,
            None => vec![f64::NAN, f64::NAN],
        };
        fst_f_mean.push(solution[0]);
        fst_f_sd.push(solution[1]);
    }
    // println!("fst_pop_a={:?}", fst_pop_a);
    // println!("fst_pop_b={:?}", fst_pop_b);
    // println!("fst_chr={:?}", fst_chr);
    // println!("fst_pos_ini={:?}", fst_pos_ini);
    // println!("fst_pos_fin={:?}", fst_pos_fin);
    // // println!("fst_f={:?}", fst_f);
    // // println!("fst_f_mean={:?}", fst_f_mean);
    // // println!("fst_f_sd={:?}", fst_f_sd);
    // println!("tajima_chr.len()={:?}", tajima_chr.len());
    // println!("fst_chr.len()={:?}", fst_chr.len());
    /////////////////////////////////////////////////////////
    // PER PAIR OF POPULATIONS: find the significant deviations in Fst
    // within the above-identified Tajima's D peaks and troughs
    let mut pop_a: Vec<String> = vec![];
    let mut pop_b: Vec<String> = vec![];
    let mut chr: Vec<Vec<String>> = vec![];
    let mut pos_ini: Vec<Vec<u64>> = vec![];
    let mut pos_fin: Vec<Vec<u64>> = vec![];
    let mut mean_tajima_d_pop_b: Vec<Vec<f64>> = vec![];
    let mut mean_fst: Vec<Vec<f64>> = vec![];
    let mut sd_tajima_d_pop_b: Vec<Vec<f64>> = vec![];
    let mut sd_fst: Vec<Vec<f64>> = vec![];
    let mut tajima_d_pop_b: Vec<Vec<f64>> = vec![];
    let mut tajima_width_pop_b: Vec<Vec<f64>> = vec![];
    let mut tajima_width_deviation_from_r_pop_b: Vec<Vec<f64>> = vec![];
    let mut tajima_width_one_tail_pval_pop_b: Vec<Vec<f64>> = vec![];
    let mut fst_delta: Vec<Vec<f64>> = vec![];
    let mut fst_delta_one_tail_pval: Vec<Vec<f64>> = vec![];
    // The recombination width in bp below is the estimated minimum width of Tajima's D troughs and peaks, however
    //      this is likely an overestimation as a function of the sparsity of the genotyping coverage.
    //      Hence, we are modeling the widths of Tajima's D as a normal distribution ,
    //      and attaching their respective one-tailed p-values.
    let recombination_width_bp = (recombination_rate_cM_per_Mb / 100.0) * 1.0e6;
    'outer: for i in 0..fst_pop_a.len() {
        let a = fst_pop_a[i].clone();
        let b = fst_pop_b[i].clone();
        // println!("a={:?}; b={:?}", a, b);
        let idx_tajima = match tajima_pop.iter().position(|x| x == &b) {
            Some(x) => x,
            None => continue 'outer,
        };
        // println!("idx_tajima={}", idx_tajima);
        // println!("tajima_d[idx_tajima]={:?}", tajima_d[idx_tajima]);
        // println!("tajima_d[idx_tajima].len()={}", tajima_d[idx_tajima].len());
        pop_a.push(a);
        pop_b.push(b);
        chr.push(vec![]);
        pos_ini.push(vec![]);
        pos_fin.push(vec![]);
        mean_tajima_d_pop_b.push(vec![]);
        mean_fst.push(vec![]);
        sd_tajima_d_pop_b.push(vec![]);
        sd_fst.push(vec![]);
        tajima_d_pop_b.push(vec![]);
        tajima_width_pop_b.push(vec![]);
        tajima_width_deviation_from_r_pop_b.push(vec![]);
        tajima_width_one_tail_pval_pop_b.push(vec![]);
        fst_delta.push(vec![]);
        fst_delta_one_tail_pval.push(vec![]);
        'inner: for j in 0..tajima_d[idx_tajima].len() {
            let tajima_window_id = tajima_chr[idx_tajima][j].to_owned()
                + ":"
                + &tajima_pos_ini[idx_tajima][j].to_string()
                + "-"
                + &tajima_pos_fin[idx_tajima][j].to_string();
            // println!("tajima_window_id={:?}", tajima_window_id);
            let idx_fst = match (0..fst_chr[i].len()).position(|idx| {
                tajima_window_id
                    == fst_chr[i][idx].to_owned()
                        + ":"
                        + &fst_pos_ini[i][idx].to_string()
                        + "-"
                        + &fst_pos_fin[i][idx].to_string()
            }) {
                Some(x) => x,
                None => continue 'inner,
            };
            // println!("idx_fst={:?}", idx_fst);
            chr[i].push(tajima_chr[idx_tajima][j].to_owned());
            pos_ini[i].push(tajima_pos_ini[idx_tajima][j]);
            pos_fin[i].push(tajima_pos_fin[idx_tajima][j]);
            mean_tajima_d_pop_b[i].push(tajima_d_mean[idx_tajima][j]);
            sd_tajima_d_pop_b[i].push(tajima_d_sd[idx_tajima][j]);
            tajima_d_pop_b[i].push(tajima_d[idx_tajima][j]);
            let width = tajima_width[idx_tajima][j] as f64;
            // Calculate Fst comparisons may be irrelevant if we have a significant Tajima's D trough/peak, i.e. Tajima's D trough/peak with is greater than zero
            tajima_width_pop_b[i].push(width);
            tajima_width_deviation_from_r_pop_b[i].push(width - recombination_width_bp);
            fst_delta[i].push(fst_f[i][idx_fst] - fst_f_mean[i]);
            mean_fst[i].push(fst_f_mean[i]); // redundant copying but I just want the all output vectors to be generated here - will need to refactor!
            sd_fst[i].push(fst_f_sd[i]); // redundant copying but I just want the all output vectors to be generated here - will need to refactor!
            let dist = Normal::new(fst_f_mean[i], fst_f_sd[i]).unwrap();
            let pval = if fst_f[i][idx_fst] < fst_f_mean[i] {
                dist.cdf(fst_f[i][idx_fst])
            } else {
                1.0 - dist.cdf(fst_f[i][idx_fst])
            };
            fst_delta_one_tail_pval[i].push(pval);
        }
        // Model the widths of the significantly deviating troughs and peaks of Tajima's D as a normal distribution
        //  and attach their respective one-tailed p-values
        let solver = prepare_solver_neldermead(2.0, 1.0);
        let solution = match ml_normal_1d(solver, &Array1::from_vec(tajima_width_pop_b[i].clone()))
        {
            Some(x) => x,
            None => vec![f64::NAN, f64::NAN],
        };
        let dist = Normal::new(solution[0], solution[1]).unwrap();
        for j in 0..tajima_width_pop_b[i].len() {
            let pval = if tajima_width_pop_b[i][j] < solution[0] {
                dist.cdf(tajima_width_pop_b[i][j])
            } else {
                1.0 - dist.cdf(tajima_width_pop_b[i][j])
            };
            tajima_width_one_tail_pval_pop_b[i].push(pval);
            // println!("pval={:?}", pval);
        }
    }
    // println!("tajima_pop={:?}", tajima_pop);

    // Write output
    let mut fname_output = fname_output.to_owned();
    if fname_output == "".to_owned() {
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        let bname = fname_input
            .split(".")
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
        fname_output = bname.to_owned() + "-gudmc-" + &time.to_string() + ".csv";
    }
    // Instantiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
    let mut file_out = fs::OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output)
        .expect(&error_writing_file);
    // Header
    let line: Vec<String> = vec![
        "pop_a",
        "pop_b",
        "chr",
        "pos_ini",
        "pos_fin",
        "mean_tajima_d_pop_b",
        "mean_fst",
        "sd_tajima_d_pop_b",
        "sd_fst",
        "tajima_d_pop_b",
        "tajima_width_pop_b",
        "tajima_width_deviation_from_r_pop_b",
        "tajima_width_one_tail_pval_pop_b",
        "fst_delta",
        "fst_delta_one_tail_pval",
    ]
    .iter()
    .map(|&x| x.to_owned())
    .collect();
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    for i in 0..pop_a.len() {
        for j in 0..chr[i].len() {
            let line = vec![
                pop_a[i].clone(),
                pop_b[i].clone(),
                chr[i][j].clone(),
                pos_ini[i][j].to_string(),
                pos_fin[i][j].to_string(),
                parse_f64_roundup_and_own(mean_tajima_d_pop_b[i][j], 7),
                parse_f64_roundup_and_own(mean_fst[i][j], 7),
                parse_f64_roundup_and_own(sd_tajima_d_pop_b[i][j], 7),
                parse_f64_roundup_and_own(sd_fst[i][j], 7),
                tajima_d_pop_b[i][j].to_string(),
                tajima_width_pop_b[i][j].to_string(),
                tajima_width_deviation_from_r_pop_b[i][j].to_string(),
                parse_f64_roundup_and_own(tajima_width_one_tail_pval_pop_b[i][j], 7),
                parse_f64_roundup_and_own(fst_delta[i][j], 7),
                parse_f64_roundup_and_own(fst_delta_one_tail_pval[i][j], 7),
            ]
            .join(",")
                + "\n";
            file_out.write_all(line.as_bytes()).unwrap();
        }
    }
    // Cleanup
    let _ = fs::remove_file("gudmc_intermediate_file_tajimasD.tmp");
    let _ = fs::remove_file("gudmc_intermediate_file_Fst.tmp");
    let _ = fs::remove_file(fname_fst);
    Ok(fname_output)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_gudmc() {
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "gudmc".to_owned(),
        };

        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        let filter_stats = FilterStats {
            remove_ns: false,
            keep_lowercase_reference: false,
            max_base_error_rate: 0.01,
            min_coverage_depth: 1,
            min_coverage_breadth: 1.0,
            min_allele_frequency: 0.001,
            max_missingness_rate: 0.0,
            pool_sizes: vec![42.0, 42.0, 42.0, 42.0, 42.0],
        };
        let genotypes_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &2)
            .unwrap(); // we need all alleles in each locus

        // Outputs
        let out = gudmc(
            &genotypes_and_phenotypes,
            &vec![42.0, 42.0, 42.0, 42.0, 42.0],
            &2.0,
            &0.73, // cM/Mb estimate in maize from [here](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-9-r103#Sec7)
            &100,
            &50,
            &20,
            &"test.something".to_owned(),
            &"test-gudmc.csv".to_owned(),
        )
        .unwrap();

        // Assertions
        assert_eq!(out, "test-gudmc.csv".to_owned());
        // assert_eq!(0, 1);
    }
}
