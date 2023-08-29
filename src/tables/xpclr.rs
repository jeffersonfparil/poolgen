use crate::base::*;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::Scalar;
use std::f32::consts::LOG10_2;
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

/// XP-CLR from https://www.ncbi.nlm.nih.gov/pubmed/20086244
/// Implementation and modifications (improvements, I hope) of https://github.com/hardingnj/xpclr/blob/master/xpclr/methods.py
pub fn xpclr_per_window(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &usize,
    min_snps_per_window: &usize,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<(String, String)> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let (loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    let l = loci_idx.len() - 1; // the number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus

    // Prepare 3-dimensional arrays for parallel computations
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
    let mut squared_differences: Array3<f64> = Array3::from_elem((n, n, l), f64::NAN);
    let mut sum_half_heterozygotes: Array3<f64> = Array3::from_elem((n, n, l), f64::NAN);
    Zip::from(&mut squared_differences)
        .and(&mut sum_half_heterozygotes)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|d, h, &i, &j, &k| {
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            let (_, a) = g.dim();
            let mut denom = 0.0;
            for allele_a in 0..(a-1){
                for allele_b in 1..a {
                    denom += g[(i, allele_a)] * g[(i, allele_b)];
                }
            }
            (*d, *h) = if denom < f64::EPSILON {
                // Catch where the allele frequency in the base population is zero and skip
                // where if the denominator is zero then the locus is fixed
                (f64::NAN, f64::NAN)
            } else {
                // Compute CV-ish thingy: squared diff / allele freq var in base pop
                let mut squared_diff: f64 = 0.0;
                for allele in 0..a {
                    let q0 = g[(j, allele)];
                    let q1 = g[(k, allele)];
                    squared_diff += (q1-q0).powf(2.0);
                }
                (squared_diff, denom) // just use these as sums across alleles or average so we have means per locus not sums across alleles per locus????
            };
        });

    // Selection-recombination
    let s = 0.5;
    let r = 1.0e-7;
    let r_min = 1.0e-9;
    let ne = 100_000;
    let c = if s <= f64::EPSILON {
        0.0
    } else {
        1.00 - f64::exp(-(2.00 * ne as f64).ln() * vec![r, r_min]
            .iter()
            .fold(r, |max, &x| if x>max{x}else{max})
        )
    };

    // Compute likelihood ratios per SNP
    // Instantiate an l x n x n array
    let mut xpclr: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN); // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    Zip::from(&mut xpclr)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|r, &i, &j, &k| {
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            let omega = 
            *r = 0.0;
        });
    // Write output (Fst averaged across all loci)
    let mut fname_output = fname_output.to_owned();
    let mut fname_output_per_window = fname_output
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
    fname_output_per_window =
        fname_output_per_window + "-xpclr-" + &window_size_bp.to_string() + "_bp_windows.csv";
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
        fname_output =
            bname.to_owned() + "-xpclr-averaged_across_genome-" + &time.to_string() + ".csv";
        fname_output_per_window = bname.to_owned()
            + "-xpclr-"
            + &window_size_bp.to_string()
            + "_bp_windows-"
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
    // Header
    let mut line: Vec<String> = vec!["".to_owned()];
    for pool in &genotypes_and_phenotypes.pool_names {
        line.push(pool.to_owned());
    }
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    // Write the mean Fst across loci
    let fst_means: Array2<f64> = xpclr.mean_axis(Axis(0)).unwrap();
    for i in 0..n {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &fst_means
                .row(i)
                .iter()
                .map(|x| parse_f64_roundup_and_own(*x, 8))
                .collect::<Vec<String>>()
                .join(",")
            + "\n";
        file_out.write_all(line.as_bytes()).unwrap();
    }
    ///////////////////////////////////////////////////////////////////////
    // Additional output: Fst per window per population
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
        return Ok((fname_output, error_message));
    }
    // println!("loci_chr={:?}", loci_chr);
    // println!("loci_pos={:?}", loci_pos);
    // println!("m={:?}", m);
    // println!("windows_idx={:?}", windows_idx);
    // println!("windows_chr={:?}", windows_chr);
    // println!("windows_pos={:?}", windows_pos);
    // println!("windows_n_sites={:?}", windows_n_sites);
    // Take the means per window
    let n_windows = windows_idx.len() - 1;
    let mut fst_per_pool_x_pool_per_window: Array2<f64> =
        Array2::from_elem((n_windows, n * n), f64::NAN);
    for i in 0..n_windows {
        let idx_start = windows_idx[i];
        let idx_end = windows_idx[i + 1];
        for j in 0..n {
            for k in 0..n {
                let idx = (j * n) + k;
                fst_per_pool_x_pool_per_window[(i, idx)] = xpclr
                    .slice(s![idx_start..idx_end, j, k])
                    .mean_axis(Axis(0))
                    .unwrap()
                    .fold(0.0, |_, &x| x);
            }
        }
    }
    // Write output (Fst averaged per window)
    // Instatiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output_per_window;
    let mut file_out = OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output_per_window)
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
        let fst_string = fst_per_pool_x_pool_per_window
            .slice(s![i, ..])
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");
        let line = coordinate + &fst_string[..] + "\n";
        file_out.write_all(line.as_bytes()).unwrap();
    }
    ///////////////////////////////////////////////////////////////////////

    Ok((fname_output, fname_output_per_window))
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
        // // Outputs
        // let (out_genomewide, out_per_window) = xpclr_per_window(
        //     &genotypes_and_phenotypes,
        //     &100,
        //     &1,
        //     &"test.something".to_owned(),
        //     &"".to_owned(),
        // )
        // .unwrap();
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
