use crate::base::*;
use crate::gwas::pearsons_correlation;
use crate::tables::theta_pi;
use ndarray::{prelude::*, Zip};
use statrs::distribution::{Binomial, Discrete};
use std::f64::consts::PI;
use std::fs::OpenOptions;
use std::io::{self, prelude::*, Error, ErrorKind};
use std::time::{SystemTime, UNIX_EPOCH};

/// XP-CLR from https://www.ncbi.nlm.nih.gov/pubmed/20086244
/// Implementation and modifications (improvements, I hope) of https://github.com/hardingnj/xpclr/blob/master/xpclr/methods.py

// Probability density functions which will be integrated between 0 and 1 as allele frequencies are restricted between 0 and 1
fn pdf_base(x: f64, other_params: &Vec<f64>) -> io::Result<f64> {
    let q0 = other_params[0]; // allele frequency in the base population
    let s2 = other_params[1]; // estimated variance of the allele frequency across the base and focal populations
    let c = other_params[2]; // selection-recombination factor which is negatively related to selection intensity
    Ok(if (x < c) & (x < (1.0 - c)) {
        (1.0 / (2.0 * PI * s2).sqrt())
            * ((c - x) / c.powf(2.0))
            * f64::exp(-(((x - q0).powf(2.0)) / (2.0 * c.powf(2.0) * s2)))
    } else if x > (1.0 - c) {
        (1.0 / (2.0 * PI * s2).sqrt())
            * ((x + c - 1.0) / c.powf(2.0))
            * f64::exp(-(((x + c - 1.0 - (c * q0)).powf(2.0)) / (2.0 * c.powf(2.0) * s2)))
    } else {
        0.0
    })
}

fn pdf_focal(x: f64, other_params: &Vec<f64>) -> io::Result<f64> {
    let q0 = other_params[0]; // allele frequency in the base population
    let s2 = other_params[1]; // estimated variance of the allele frequency across the base and focal populations
    let c = other_params[2]; // selection-recombination factor which is negatively related to selection intensity
    let k = other_params[3] as u64; // allele count in the focal population
    let m = other_params[4] as u64; // total coverage (total number of all the alleles in the current locus) in the focal population
    let mut dist_binomial = Binomial::new(x, m).unwrap();
    Ok(pdf_base(x, &vec![q0, s2, c]).unwrap() * dist_binomial.pmf(k))
}

// Modification of https://codereview.stackexchange.com/a/161366
fn simple_reimann_sum_integration<F>(
    a: f64,
    b: f64,
    func: F,
    precision: u64,
    other_params: &Vec<f64>,
) -> io::Result<f64>
where
    F: Fn(f64, &Vec<f64>) -> io::Result<f64>,
{
    let delta: f64 = (b - a) / precision as f64;
    let sum: f64 = (0..precision)
        .map(|trapezoid| {
            let left_side = a + (delta * trapezoid as f64);
            let right_size = left_side + delta;
            0.5 * (func(left_side, other_params).unwrap() + func(right_size, other_params).unwrap())
                * delta
        })
        .sum();
    Ok(if a > b { -sum } else { sum })
}

// Compute XP-CLR across all pairs of pools per window
pub fn xpclr_per_window(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &usize,
    min_loci_per_window: &usize,
    selection_coefficient: &f64,
    recombination_rate: &f64,
    min_recombination_rate: &f64,
    effective_population_sizes: &Array1<f64>,
    integration_precision: &u64,
    correlation_threshold_between_loci: &f64,
) -> io::Result<Array2<f64>> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let (loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    let l = loci_idx.len() - 1; // the number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus

    // Prepare 3-dimensional (l x n x n) input arrays for parallel computations
    let loci: Array3<usize> = Array3::from_shape_vec(
        (l, n, n),
        (0..l)
            .flat_map(|x| std::iter::repeat(x).take(n * n))
            .collect(),
    )
    .unwrap();
    let pop1: Array3<usize> = Array3::from_shape_vec(
        (l, n, n),
        std::iter::repeat(
            (0..n)
                .flat_map(|x| std::iter::repeat(x).take(n))
                .collect::<Vec<usize>>(),
        )
        .take(l)
        .flat_map(|x| x)
        .collect::<Vec<usize>>(),
    )
    .unwrap();
    let pop2: Array3<usize> = Array3::from_shape_vec(
        (l, n, n),
        std::iter::repeat(
            std::iter::repeat((0..n).collect::<Vec<usize>>())
                .take(n)
                .flat_map(|x| x)
                .collect::<Vec<usize>>(),
        )
        .take(l)
        .flat_map(|x| x)
        .collect::<Vec<usize>>(),
    )
    .unwrap();
    // Estimate the effect of population histories between pairs of populations, i.e. omega (representing demography, splits, drift, bottleneck...)
    // where if the allele frequency at the base population, q0 is zero then we skip as we assume that the base population is the ancestral population
    // Generate an l x n x n array which is non-symmetric at nxn because of the reason above, i.e. allele frequencies cannot be fixed in the base (ancestral) population
    let mut squared_differences: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN);
    let mut denom_binomial_var: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN);
    Zip::from(&mut squared_differences)
        .and(&mut denom_binomial_var)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|d, v, &i, &j, &k| {
            if j==k {
                // If the base and focal populations are the same
                (*d, *v) = (f64::NAN, f64::NAN);
            } else {
                let idx_start = loci_idx[i];
                let idx_end = loci_idx[i + 1];
                let g = genotypes_and_phenotypes
                    .intercept_and_allele_frequencies
                    .slice(s![.., idx_start..idx_end]);
                // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
                let (idx_allele, _) = g.sum_axis(Axis(0)).iter().enumerate().fold(
                    (0, 0.0),
                    |(idx_max, max), (idx, &x)| if x > max { (idx, x) } else { (idx_max, max) },
                );
                let q0 = g[(j, idx_allele)];
                let q1 = g[(k, idx_allele)];
                (*d, *v) = if (q0 < f64::EPSILON) | (q0 == 1.00) {
                    // Catch where the allele frequency in the base population is zero and skip
                    // where if the denominator is zero then the locus is fixed
                    (f64::NAN, f64::NAN)
                } else {
                    ((q1 - q0).powf(2.0), (q0 * (1.00 - q0)))
                };
            }
        });
    println!("squared_differences={:?}", squared_differences);
    println!("denom_binomial_var={:?}", denom_binomial_var);
    println!(
        "&squared_differences / &denom_binomial_var = {:?}",
        &squared_differences / &denom_binomial_var
    );
    let omega: Array2<f64> =
        mean_axis_ignore_nan(&squared_differences / &denom_binomial_var, 0).unwrap();
    // Neutral allele frequency variance across time (sigma2) per locus (l x n x n)
    let mut sigma2: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN);
    for i in 0..l {
        for j in 0..n {
            for k in 0..n {
                sigma2[(i, j, k)] = omega[(j, k)] * denom_binomial_var[(i, j, k)];
            }
        }
    }
    println!("omega={:?}", omega);
    println!("sigma2={:?}", sigma2);
    // println!("loci_idx={:?}", loci_idx);
    // Below is the parallel computations that takes the most time! 20230906
    // Compute likelihood ratios per SNP
    // Instantiate an l x n x n array
    let mut xpclr: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN); // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    Zip::from(&mut xpclr)
        .and(&sigma2)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|ratio, &s2, &i, &j, &k| {
            if j==k {
                // If the base and focal populations are the same
                *ratio = f64::NAN;
            } else  {
                let idx_start = loci_idx[i];
                let idx_end = loci_idx[i + 1];
                let g = genotypes_and_phenotypes
                    .intercept_and_allele_frequencies
                    .slice(s![.., idx_start..idx_end]);
                // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
                let (idx_allele, _) = g.sum_axis(Axis(0)).iter().enumerate().fold(
                    (0, 0.0),
                    |(idx_max, max), (idx, &x)| if x > max { (idx, x) } else { (idx_max, max) },
                );
                let q0 = g[(j, idx_allele)];
                let q1 = g[(k, idx_allele)];

                let ne_mean = (effective_population_sizes[j] + effective_population_sizes[k]) / 2.0;
                let c = if *selection_coefficient <= f64::EPSILON {
                    0.0
                } else {
                    1.00 - f64::exp(
                        -(2.00 * ne_mean).ln()
                            * vec![recombination_rate, min_recombination_rate]
                                .iter()
                                .fold(recombination_rate, |max, &x| if x > max { x } else { max })
                            / selection_coefficient,
                    )
                };

                let m = genotypes_and_phenotypes.coverages[(k, i)];
                let k = (q1 * m).round();

                let params_of_pdf_p1 = vec![q0, s2, c];
                let params_of_cdf_p1 = vec![q0, s2, c, k, m];
                let likelihood_base = simple_reimann_sum_integration(
                    0.001,
                    0.999,
                    pdf_base,
                    *integration_precision,
                    &params_of_pdf_p1,
                )
                .unwrap();
                let likelihood_focal = simple_reimann_sum_integration(
                    0.001,
                    0.999,
                    pdf_focal,
                    *integration_precision,
                    &params_of_cdf_p1,
                )
                .unwrap();
                *ratio = if (likelihood_focal < f64::EPSILON) | (likelihood_base < f64::EPSILON) {
                    f64::EPSILON.ln()
                } else {
                    likelihood_focal.ln() - likelihood_base.ln()
                };
            }
        });
    println!("xpclr={:?}", xpclr);
    ///////////////////////////////////////////////////////////////////////
    // Output: XP-CLR per window per population pair
    // Summarize per non-overlapping window
    // Find window indices making sure we respect chromosomal boundaries
    // while filtering out windows with less than min_loci_per_window SNPs
    let m = loci_idx.len() - 1; // total number of loci, we subtract 1 as the last index refer to the last allele of the last locus and serves as an end marker
    let mut windows_idx: Vec<usize> = vec![0]; // indices in terms of the number of loci not in terms of genome coordinates - just to make it simpler
    let mut windows_chr: Vec<String> = vec![loci_chr[0].to_owned()];
    let mut windows_pos: Vec<u64> = vec![loci_pos[0] as u64];
    let mut windows_n_sites: Vec<usize> = vec![0];
    let mut j = windows_n_sites.len() - 1; // number of sites per window whose length is used to count the current number of windows
    for i in 0..m {
        let chr = loci_chr[i].to_owned(); // skipping the intercept at position 0
        let pos = loci_pos[i]; // skipping the intercept at position 0
        if (chr != windows_chr.last().unwrap().to_owned())
            | ((chr == windows_chr.last().unwrap().to_owned())
                & (pos > windows_pos.last().unwrap() + &(*window_size_bp as u64)))
        {
            if windows_n_sites[j] < *min_loci_per_window {
                windows_idx[j] = i;
                windows_chr[j] = chr.to_owned();
                windows_pos[j] = pos;
                windows_n_sites[j] = 1;
            } else {
                windows_idx.push(i);
                windows_chr.push(chr.to_owned());
                windows_pos.push(pos);
                windows_n_sites.push(1);
            }
        } else {
            windows_n_sites[j] += 1;
        }
        j = windows_n_sites.len() - 1;
    }
    // Add the last index of the final position
    windows_idx.push(m);
    windows_chr.push(windows_chr.last().unwrap().to_owned());
    windows_pos.push(*loci_pos.last().unwrap());
    if windows_n_sites.last().unwrap() < min_loci_per_window {
        windows_idx.pop();
        windows_chr.pop();
        windows_pos.pop();
        windows_n_sites.pop();
    }
    if windows_n_sites.len() < 1 {
        let error_message =
            "No window with at least ".to_owned() + &min_loci_per_window.to_string() + " SNPs.";
        return Err(Error::new(ErrorKind::Other, error_message));
    }
    // println!("loci_chr={:?}", loci_chr);
    // println!("loci_pos={:?}", loci_pos);
    // println!("m={:?}", m);
    // println!("windows_idx={:?}", windows_idx);
    // println!("windows_chr={:?}", windows_chr);
    // println!("windows_pos={:?}", windows_pos);
    println!("windows_n_sites={:?}", windows_n_sites);
    // Estimate weights per locus
    // Computing the correlation between SNPs across pools instead of between 2 populations and using individual genotypes because we have Pool-seq data, i.e. individual genotypes are not available
    let mut weights: Array1<f64> = Array1::from_elem(l, 0.0);
    for i in 0..l {
        let idx_start = loci_idx[i];
        let idx_end = loci_idx[i + 1];
        let g = genotypes_and_phenotypes
            .intercept_and_allele_frequencies
            .slice(s![.., idx_start..idx_end]);
        // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
        let (idx_allele, _) =
            g.sum_axis(Axis(0))
                .iter()
                .enumerate()
                .fold(
                    (0, 0.0),
                    |(idx_max, max), (idx, &x)| if x > max { (idx, x) } else { (idx_max, max) },
                );
        let qi = g.column(idx_allele);
        for j in 0..l {
            let idx_start = loci_idx[j];
            let idx_end = loci_idx[j + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
            let (idx_allele, _) = g.sum_axis(Axis(0)).iter().enumerate().fold(
                (0, 0.0),
                |(idx_max, max), (idx, &x)| if x > max { (idx, x) } else { (idx_max, max) },
            );
            let qj = g.column(idx_allele);
            let (corr, _) = pearsons_correlation(&qi, &qj).unwrap();
            // println!("qi={:?}", qi);
            // println!("qj={:?}", qj);
            // println!("corr={:?}", corr);
            weights[i] += if corr < *correlation_threshold_between_loci {
                0.0
            } else {
                1.0
            };
        }
        // println!("weights={:?}", weights);
        weights[i] = 1.0 / weights[i];
    }
    // Compute the cross-population composite likelihood ratio weighted by the inverse of the level of correlation between loci per window
    let n_windows = windows_idx.len() - 1;
    let mut xpclr_per_window_per_pop_pair: Array2<f64> = Array2::from_elem((n_windows, n * n), 0.0);
    for i in 0..n_windows {
        let idx_start = windows_idx[i];
        let idx_end = windows_idx[i + 1];
        for j in 0..n {
            for k in 0..n {
                let idx = (j * n) + k;
                let mut n_non_nan = 0;
                for locus_idx in idx_start..idx_end {
                    xpclr_per_window_per_pop_pair[(i, idx)] +=
                        if xpclr[(locus_idx, j, k)].is_nan() == false {
                            n_non_nan += 1;
                            weights[locus_idx] * xpclr[(locus_idx, j, k)]
                        } else {
                            0.0
                        };
                }
                if n_non_nan == 0 {
                    xpclr_per_window_per_pop_pair[(i, idx)] = f64::NAN;
                }
            }
        }
    }
    println!("weights={:?}", weights);
    println!(
        "xpclr_per_window_per_pop_pair={:?}",
        xpclr_per_window_per_pop_pair
    );
    Ok(xpclr_per_window_per_pop_pair)
}

// Compute and print out in a single csv file the XP-CLR across all pairs of pools per window
// across a range of selection coefficients and recombination rates
pub fn xpclr(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &usize,
    min_loci_per_window: &usize,
    integration_precision: &u64,
    correlation_threshold_between_loci: &f64,
    selection_coefficient_min: &f64,
    selection_coefficient_max: &f64,
    selection_coefficient_n_steps: &u64,
    recombination_rate_min: &f64,
    recombination_rate_max: &f64,
    recombination_rate_n_steps: &u64,
    mutation_rate: &f64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    let s_step_size = (selection_coefficient_max - selection_coefficient_min)
        / (*selection_coefficient_n_steps as f64);
    let r_step_size =
        (recombination_rate_max - recombination_rate_min) / (*recombination_rate_n_steps as f64);
    let range_of_selection_coefficients: Vec<f64> = (0..*selection_coefficient_n_steps)
        .map(|x| ((x as f64) * s_step_size) + selection_coefficient_min)
        .collect();
    let range_of_recombination_rates: Vec<f64> = (0..*recombination_rate_n_steps)
        .map(|x| ((x as f64) * r_step_size) + recombination_rate_min)
        .collect();
    let min_recombination_rate = recombination_rate_min;
    // println!("range_of_selection_coefficients={:?}", range_of_selection_coefficients);
    // println!("range_of_recombination_rates={:?}", range_of_recombination_rates);
    // Calculate heterozygosities to estimate the effective population sizes (theta = 4*Ne*mu)
    let (pi_per_pool_per_window, windows_chr, windows_pos) = theta_pi(
        genotypes_and_phenotypes,
        window_size_bp,
        min_loci_per_window,
    )
    .unwrap();
    // Estimate effective population sizes per pool
    let n = pi_per_pool_per_window.ncols();
    let n_windows = pi_per_pool_per_window.nrows();
    let pi_per_pool: Array1<f64> = pi_per_pool_per_window.mean_axis(Axis(0)).unwrap();
    let effective_population_sizes: Array1<f64> = pi_per_pool / (4.00 * mutation_rate);
    // Define the output filename
    let fname_output = if fname_output.to_owned() == "".to_owned() {
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
        let time = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64();
        bname + "-xpclr-" + &time.to_string() + ".csv"
    } else {
        fname_output.to_owned()
    };
    // Estimate XP-CLR using the mean effective population sizes between pairs of pools
    let mut vec_s: Vec<f64> = vec![];
    let mut vec_r: Vec<f64> = vec![];
    let mut vec_output: Vec<Array2<f64>> = vec![];
    for selection_coefficient in &range_of_selection_coefficients {
        for recombination_rate in &range_of_recombination_rates {
            vec_s.push(*selection_coefficient);
            vec_r.push(*recombination_rate);
            vec_output.push(
                xpclr_per_window(
                    genotypes_and_phenotypes,
                    window_size_bp,
                    min_loci_per_window,
                    &selection_coefficient,
                    &recombination_rate,
                    min_recombination_rate,
                    &effective_population_sizes,
                    integration_precision,
                    correlation_threshold_between_loci,
                )
                .unwrap(),
            );
        }
    }
    // println!("vec_s={:?}", vec_s);
    // println!("vec_r={:?}", vec_r);
    // println!("vec_output={:?}", vec_output);
    // Instantiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
    let mut file_out = OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output)
        .expect(&error_writing_file);
    // Header
    let mut line: Vec<String> = vec![
        "selection_coefficient".to_owned(),
        "recombination_rate".to_owned(),
        "chr".to_owned(),
        "pos_ini".to_owned(),
        "pos_fin".to_owned(),
    ];
    for j in 0..n {
        for k in 0..n {
            line.push(
                genotypes_and_phenotypes.pool_names[j].to_owned()
                    + "_into_"
                    + &genotypes_and_phenotypes.pool_names[k],
            );
        }
    }
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    // Write per line
    for j in 0..vec_output.len() {
        let xpclr_per_window_per_pop_pair = &vec_output[j];
        // let selection_recombination = sensible_round(vec_s[j], 4).to_string() + "," + &sensible_round(vec_r[j], 4).to_string() + ",";
        let selection_recombination = vec_s[j].to_string() + "," + &vec_r[j].to_string() + ",";
        for i in 0..n_windows {
            let window_pos_ini = windows_pos[i];
            let window_pos_fin = window_pos_ini + (*window_size_bp as u64);
            let coordinate = windows_chr[i].clone()
                + ","
                + &window_pos_ini.to_string()
                + ","
                + &window_pos_fin.to_string()
                + ",";
            let xpclr_string = xpclr_per_window_per_pop_pair
                .slice(s![i, ..])
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
                .join(",");
            let line = selection_recombination.to_owned() + &coordinate + &xpclr_string[..] + "\n";
            file_out.write_all(line.as_bytes()).unwrap();
        }
    }
    Ok(fname_output)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_xpclr() {
        let x: Array2<f64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.6, 0.4, 0.0,
                0.9, 0.1, 1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0,
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
        let window_size_bp = &100;
        let min_loci_per_window = &1;
        let selection_coefficient = &0.05;
        let recombination_rate = &1.0e-5;
        let min_recombination_rate = &1.0e-9;
        let effective_population_size =
            &Array1::from_vec(vec![9_000., 9_500., 10_000., 11_000., 12_000.]);
        let integration_precision = &100_000;
        let correlation_threshold_between_loci = &0.5;
        // Outputs
        let out = xpclr_per_window(
            &genotypes_and_phenotypes,
            window_size_bp,
            min_loci_per_window,
            selection_coefficient,
            recombination_rate,
            min_recombination_rate,
            effective_population_size,
            integration_precision,
            correlation_threshold_between_loci,
        )
        .unwrap();

        let selection_coefficient_min = 0.0;
        let selection_coefficient_max = 1.0;
        let selection_coefficient_n_steps = 10;
        let recombination_rate_min = 1.0e-8;
        let recombination_rate_max = 1.0e-3;
        let recombination_rate_n_steps = 10;
        // let selection_coefficient_min = 0.5;
        // let selection_coefficient_max = 0.5;
        // let selection_coefficient_n_steps = 1;
        // let recombination_rate_min = 1.0e-3;
        // let recombination_rate_max = 1.0e-3;
        // let recombination_rate_n_steps = 1;
        let mutation_rate = 9.01e-9; // Mutation rate: 9.01 x 10^(-9) - estimates from Oryza sativa (not using Arabipsis thaliana estimates of 4.35 x 10(-9)) from Wang et a, 2019
        let out = xpclr(
            &genotypes_and_phenotypes,
            window_size_bp,
            min_loci_per_window,
            integration_precision,
            correlation_threshold_between_loci,
            &selection_coefficient_min,
            &selection_coefficient_max,
            &selection_coefficient_n_steps,
            &recombination_rate_min,
            &recombination_rate_max,
            &recombination_rate_n_steps,
            &mutation_rate,
            &"test.sync".to_owned(),
            &"".to_owned(),
        )
        .unwrap();
        // assert_eq!(0, 1);
    }
}

// // Using SLiM-simulated data to test/validate XP-CLR
// cd poolgen/tests
// echo '
// initialize() {
//     // Define simulation variables
//     defineGlobal("N", 1e4);
//     defineGlobal("CHRPOS_INI",     0);
//     defineGlobal("CHRPOS_END", 9999);
//     defineGlobal("CHRPOS_SEL", 5000);
//     defineGlobal("SEL_COEF", 1.0);
//     // defineGlobal("MU", 1e-7);
//     defineGlobal("MU", 1e-5);
//     defineGlobal("R", 1e-8);
//     // Define the genomes and the initial population
//     initializeMutationRate(MU);
//     initializeMutationType("m1", 0.5, "f", 0.0);
//     initializeGenomicElementType("g1", m1, 1.0);
//     initializeGenomicElement(g1, CHRPOS_INI, CHRPOS_END);
//     initializeRecombinationRate(R);
//     // Define the mutation with selective advantage which will be introduced later into population 2
//     initializeMutationType("m2", 1.0, "f", SEL_COEF); // introduced mutation (initially neutral)
// }
// 1 early() {
//     sim.addSubpop("p1", N);
// }
// 1000 late() {
//     // sim.outputFull();
//     // p1.outputSample(500);
//     sim.addSubpopSplit("p2", N, p1);
//     target = sample(p2.genomes, 100);
//     target.addNewDrawnMutation(m2, CHRPOS_SEL);
// }
// 1100 late() {
//     p1.outputSample(1000, filePath="test-p1.txt"); // for 500 individuals or 2 x 500 haploid genomes
//     p2.outputSample(1000, filePath="test-p2.txt"); // for 500 individuals or 2 x 500 haploid genomes
// }
// ' > test.slim
// time slim test.slim > /dev/null
// //

// // Parsing in julia because python is too slow
// using ProgressMeter
// using UnicodePlots

// function parse_slim_output(;fname::String="test-p1.txt", n::Int64=500, l::Int64=10000)::Vector{Float64}
//     # fname = "test-p1.txt"
//     # n = 500
//     # l = 10000

//     i = 1
//     mut_id = []
//     mut_type = []
//     mut_pos = []
//     mutations_data_bool = false
//     genotype_data_bool = false

//     X = zeros(Float64, n, l)
//     f = open(fname)

//     ini = 0
//     while !eof(f)
//         x = readline(f)
//         if x == "Genomes:"
//             ini = ini + position(f)
//             break
//         end
//     end
//     seekend(f)
//     fin = position(f) - ini
//     seekstart(f)
//     p = Progress(fin, dt=1.0)

//     while !eof(f)
//         x = readline(f)
//         if x == "Mutations:"
//             mutations_data_bool = true
//             continue
//         end
//         if x == "Genomes:"
//             mutations_data_bool = false
//             genotype_data_bool = true
//             continue
//         end
//         if mutations_data_bool
//             xsplit = split(x, " ")
//             push!(mut_id, parse(Int64, xsplit[1]))
//             push!(mut_type, xsplit[3])
//             push!(mut_pos, parse(Int64, xsplit[4]))
//         end
//         if genotype_data_bool
//             xsplit = split(x, " ")
//             m = length(xsplit)
//             for _j = parse.(Int64, xsplit[3:m])
//                 # _j = parse(Int64, xsplit[3])
//                 j = mut_pos[mut_id .== _j][1] + 1 ### correct for julia starting with 1 instead of zero index
//                 X[i, j] = 1.0
//             end
//             i += 1
//             update!(p, position(f) - ini)
//         end
//     end
//     close(f)
//     q = sum(X, dims=1)[1,:] ./ n
//     return(q)
// end

// n = 500
// l = 10000
// q1 = parse_slim_output(fname="test-p1.txt", n=n, l=l)
// q2 = parse_slim_output(fname="test-p2.txt", n=n, l=l)

// x = collect(1:l)
// lineplot(x, q1)
// lineplot(x, q2)
