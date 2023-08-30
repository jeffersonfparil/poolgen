use crate::base::*;
use ndarray::{prelude::*, Zip};
use std::f64::consts::PI;
use statrs::distribution::{Binomial, Discrete};
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};
use crate::gwas::*;

/// XP-CLR from https://www.ncbi.nlm.nih.gov/pubmed/20086244
/// Implementation and modifications (improvements, I hope) of https://github.com/hardingnj/xpclr/blob/master/xpclr/methods.py

fn pdf_p1(p1: f64, other_params: &Vec<f64>) -> io::Result<f64> {
    let q0 = other_params[0];
    let s2 = other_params[1];
    let c = other_params[2];
    Ok(
        if (p1 < c) & (p1 < (1.0-c)) {
            (1.0 / (2.0*PI*s2).sqrt()) 
            * ((c-p1)/c.powf(2.0)) 
            * f64::exp(-( ((p1-q0).powf(2.0)) / (2.0*c.powf(2.0)*s2) ))
        } else if p1 > (1.0-c) {
            (1.0 / (2.0*PI*s2).sqrt()) 
            * ((p1+c-1.0)/c.powf(2.0)) 
            * f64::exp(-( ((p1+c-1.0-(c*q0)).powf(2.0)) / (2.0*c.powf(2.0)*s2) ))
        } else {
            0.0
        }
    )
}

fn cdf_p1(p1: f64, other_params: &Vec<f64>) -> io::Result<f64> {
    let q0 = other_params[0];
    let s2 = other_params[1];
    let c = other_params[2];
    let k = other_params[3] as u64;
    let m = other_params[4] as u64;
    let mut dist_binomial = Binomial::new(p1, m).unwrap();
    Ok(pdf_p1(p1, &vec![q0, s2, c]).unwrap() * dist_binomial.pmf(k))
}

// Modification of https://codereview.stackexchange.com/a/161366
fn simple_reimann_sum_integration<F>(a: f64, b: f64, func: F, precision: u64, other_params: &Vec<f64>) -> io::Result<f64>
    where F: Fn(f64, &Vec<f64>) -> io::Result<f64>
{
    let delta: f64 = (b - a) / precision as f64;
    let sum: f64 = (0..precision).map(|trapezoid| {
        let left_side = a + (delta * trapezoid as f64);
        let right_size = left_side + delta;
        0.5 * (func(left_side, other_params).unwrap() + func(right_size, other_params).unwrap()) * delta
        }).sum();
    Ok (
        if a > b {
            -sum
        } else {
            sum
        }
    )
}

pub fn xpclr_per_window(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &usize,
    min_snps_per_window: &usize,
    fname_output: &String,
) -> io::Result<String> {
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
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
            let (idx_allele, _) = g.sum_axis(Axis(0))
                .iter()
                .enumerate()
                .fold((0, 0.0), |(idx_max, max), (idx, &x)| if x>max{(idx, x)}else{(idx_max, max)});
            let q0 = g[(j, idx_allele)];
            let q1 = g[(k, idx_allele)];
            (*d, *v) = if (q0<f64::EPSILON) | (q0==1.00) {
                // Catch where the allele frequency in the base population is zero and skip
                // where if the denominator is zero then the locus is fixed
                (f64::NAN, f64::NAN)
            } else {
                ((q1-q0).powf(2.0), (q0*(1.00-q0)))
            };
        }
    );
    println!("squared_differences={:?}", squared_differences);
    println!("denom_binomial_var={:?}", denom_binomial_var);
    // Population history factor (omega) across all loci (n x n)
    let omega: Array2<f64> = (&squared_differences / &denom_binomial_var).mean_axis(Axis(0)).unwrap();
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
    // Joint effects of selection and recombination (c) across all loci and populations (scalar)
    let s = 0.05;
    let r = 1.0e-5;
    // let r = 1.0e-7;
    let r_min = 1.0e-9;
    let ne = 10_000;
    let precision = 100_000;
    let corr_cutoff = 0.5;
    let c = if s <= f64::EPSILON {
        0.0
    } else {
        1.00 - f64::exp(-(2.00 * ne as f64).ln() * vec![r, r_min]
            .iter()
            .fold(r, |max, &x| if x>max{x}else{max}) / s
        )
    };
    // let c = 0.50;
    println!("c={:?}", c);
    let params_of_pdf_p1 = vec![0.5, 0.1, c];
    let q = simple_reimann_sum_integration(0.001, 0.999, pdf_p1, precision, &params_of_pdf_p1).unwrap();
    println!("q={:?}", q);
    let q_001 = pdf_p1(0.001, &params_of_pdf_p1).unwrap();
    let q_500 = pdf_p1(0.500, &params_of_pdf_p1).unwrap();
    let q_999 = pdf_p1(0.999, &params_of_pdf_p1).unwrap();
    println!("q_001={:?}", q_001);
    println!("q_500={:?}", q_500);
    println!("q_999={:?}", q_999);
    println!("f64::EPSILON={:?}", f64::EPSILON);
    // Compute likelihood ratios per SNP
    // Instantiate an l x n x n array
    let mut xpclr: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN); // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    Zip::from(&mut xpclr)
        .and(&sigma2)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|ratio, &s2, &i, &j, &k| {
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
            let (idx_allele, _) = g.sum_axis(Axis(0))
                .iter()
                .enumerate()
                .fold((0, 0.0), |(idx_max, max), (idx, &x)| if x>max{(idx, x)}else{(idx_max, max)});
            let q0 = g[(j, idx_allele)];
            let q1 = g[(k, idx_allele)];

            let m = genotypes_and_phenotypes.coverages[(k, i)];
            let k = (q1 * m).round();

            let params_of_pdf_p1 = vec![q0, s2, c];
            let params_of_cdf_p1 = vec![q0, s2, c, k, m];
            let likelihood_b = simple_reimann_sum_integration(0.001, 0.999, pdf_p1, precision, &params_of_pdf_p1).unwrap();
            let likelihood_a = simple_reimann_sum_integration(0.001, 0.999, cdf_p1, precision, &params_of_cdf_p1).unwrap();
            *ratio = if (likelihood_a < f64::EPSILON) | (likelihood_b < f64::EPSILON) {
                f64::EPSILON.ln()
            } else {
                likelihood_a.ln() - likelihood_b.ln()
            };
        }
    );
    println!("xpclr={:?}", xpclr);
    ///////////////////////////////////////////////////////////////////////
    // Output: XP-CLR per window per population pair
    // Summarize per non-overlapping window
    // Find window indices making sure we respect chromosomal boundaries
    // while filtering out windows with less than min_snps_per_window SNPs
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
            if windows_n_sites[j] < *min_snps_per_window {
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
    if windows_n_sites.last().unwrap() < min_snps_per_window {
        windows_idx.pop();
        windows_chr.pop();
        windows_pos.pop();
        windows_n_sites.pop();
    }
    if windows_n_sites.len() < 1 {
        let error_message =
            "No window with at least ".to_owned() + &min_snps_per_window.to_string() + " SNPs.";
        return Ok(error_message);
    }
    println!("loci_chr={:?}", loci_chr);
    println!("loci_pos={:?}", loci_pos);
    println!("m={:?}", m);
    println!("windows_idx={:?}", windows_idx);
    println!("windows_chr={:?}", windows_chr);
    println!("windows_pos={:?}", windows_pos);
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
        let (idx_allele, _) = g.sum_axis(Axis(0))
            .iter()
            .enumerate()
            .fold((0, 0.0), |(idx_max, max), (idx, &x)| if x>max{(idx, x)}else{(idx_max, max)});
        let qi = g.column(idx_allele);
        for j in 0..l {
            let idx_start = loci_idx[j];
            let idx_end = loci_idx[j + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            // Using the major allele as the focal allele (NOTE: will need to account for multiallelic loci in a more satisfying way!)
            let (idx_allele, _) = g.sum_axis(Axis(0))
                .iter()
                .enumerate()
                .fold((0, 0.0), |(idx_max, max), (idx, &x)| if x>max{(idx, x)}else{(idx_max, max)});
            let qj = g.column(idx_allele);    
            let (corr, _) = pearsons_correlation(&qi, &qj).unwrap();
            // println!("qi={:?}", qi);
            // println!("qj={:?}", qj);
            // println!("corr={:?}", corr);
            weights[i] += if corr < corr_cutoff {
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
    let mut xpclr_per_window_per_pop_pair: Array2<f64> =
        Array2::from_elem((n_windows, n * n), 0.0);
    for i in 0..n_windows {
        let idx_start = windows_idx[i];
        let idx_end = windows_idx[i + 1];
        for j in 0..n {
            for k in 0..n {
                let idx = (j * n) + k;
                let mut n_non_nan = 0;
                for locus_idx in idx_start..idx_end {
                    xpclr_per_window_per_pop_pair[(i, idx)] += if xpclr[(locus_idx, j, k)].is_nan() == false {
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
    println!("xpclr_per_window_per_pop_pair={:?}", xpclr_per_window_per_pop_pair);
    // Instantiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
    let mut file_out = OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output)
        .expect(&error_writing_file);
    // Header
    let mut line: Vec<String> = vec!["chr".to_owned(), "pos_ini".to_owned(), "pos_fin".to_owned()];
    for j in 0..n {
        for k in 0..n {
            line.push(
                genotypes_and_phenotypes.pool_names[j].to_owned()
                    + "_vs_"
                    + &genotypes_and_phenotypes.pool_names[k],
            );
        }
    }
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    // Per line
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
        let line = coordinate + &xpclr_string[..] + "\n";
        file_out.write_all(line.as_bytes()).unwrap();
    }
    ///////////////////////////////////////////////////////////////////////

    Ok(fname_output.to_owned())
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
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4,
                1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 
                1.0, 0.6, 0.4, 0.0, 0.9, 0.1, 
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 
                1.0, 1.0, 0.0, 0.0, 0.0, 1.0,
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
        let out = xpclr_per_window(
            &genotypes_and_phenotypes,
            &100,
            &1,
            &"test-something.csv".to_owned(),
        )
        .unwrap();
    assert_eq!(0, 1);
        // let file = std::fs::File::open(&out_genomewide).unwrap();
        // let reader = std::io::BufReader::new(file);
        // let mut header: Vec<String> = vec![];
        // let mut pop: Vec<String> = vec![];
        // let mut xpclr: Vec<f64> = vec![];
        // for line in reader.lines() {
        //     let split = line
        //         .unwrap()
        //         .split(",")
        //         .map(|x| x.to_owned())
        //         .collect::<Vec<String>>();
        //     if header.len() == 0 {
        //         header = split;
        //     } else {
        //         pop.push(split[0].clone());
        //         for f in split[1..]
        //             .iter()
        //             .map(|x| x.parse::<f64>().unwrap())
        //             .collect::<Vec<f64>>()
        //         {
        //             xpclr.push(f);
        //         }
        //     }
        // }
        // let xpclr: Array2<f64> = Array2::from_shape_vec((pop.len(), pop.len()), xpclr).unwrap();
        // let diag: Array1<f64> = xpclr.diag().to_owned(); // itself, i.e. xpclr=0.0
        // let pop1_4 = xpclr[(0, 3)]; // the same populations, i.e. xpclr=0.0
        // let pop4_1 = xpclr[(3, 0)]; // the same populations, i.e. xpclr=0.0
        // let pop2_5 = xpclr[(1, 4)]; // totally different populations, i.e. xpclr=0.5, the same locus 1 and different locus 2
        // let pop5_2 = xpclr[(4, 1)]; // totally different populations, i.e. xpclr=0.5, the same locus 1 and different locus 2
        // println!("pop={:?}", pop);
        // println!("xpclr={:?}", xpclr);
        // println!("out_per_window={:?}", out_per_window);
        // // Assertions
        // assert_eq!(diag, Array1::from_elem(pop.len(), 0.0));
        // assert_eq!(pop1_4, 0.0);
        // assert_eq!(pop4_1, 0.0);
        // assert_eq!(pop2_5, 0.5);
        // assert_eq!(pop5_2, 0.5);
        // assert_eq!(0, 1);
    }
}
