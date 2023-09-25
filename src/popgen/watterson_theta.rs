use crate::base::*;
use ndarray::prelude::*;
use std::fs::OpenOptions;
use std::io::{self, prelude::*, Error, ErrorKind};
use std::time::{SystemTime, UNIX_EPOCH};

/// Count the number of polymorphic sites per pool
fn polymorphic_loci_per_pool(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    loci_idx: &Vec<usize>,
    idx: usize,
) -> io::Result<Vec<u64>> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let mut out: Vec<u64> = std::iter::repeat(0).take(n).collect();
    for k in 0..n {
        let idx_allele_ini = loci_idx[idx + 0];
        let idx_allele_fin = loci_idx[idx + 1];
        let freq_max = genotypes_and_phenotypes
            .intercept_and_allele_frequencies
            .slice(s![k, idx_allele_ini..idx_allele_fin])
            .fold(0.0, |max, &x| if x > max { x } else { max });
        if freq_max < 1.0 {
            // assumes we are keeping all the alleles per locus
            out[k] += 1;
        }
    }
    Ok(out)
}

/// Watterson's estimator of theta
/// Simply (naively?) defined as $theta_w = number of segregating sites / \Sigma^{n-1}_{i=1}(1/i)$
/// For details see [Feretti et al, 2013](https://doi.org/10.1111/mec.12522)
pub fn theta_watterson(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
) -> io::Result<(Array2<f64>, Vec<usize>, Vec<usize>)> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let (loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    // Remove redundant trailing loci from `genotypes_and_phenotypes.count_loci()`
    // let mut loci_idx = loci_idx.to_owned(); loci_idx.pop(); // Do not remove the trailing index for loci_idx as it is needed in `polymorphic_loci_per_pool` to define the slice in `genotypes_and_phenotypes.intercept_and_allele_frequencies`
    let mut loci_chr = loci_chr.to_owned();
    loci_chr.pop();
    let mut loci_pos = loci_pos.to_owned();
    loci_pos.pop();
    //////////////////////////////////////////////////////////////////////////////////////////////
    // The following code is similar to `src/base/helpers.rs define_sliding_windows()` function //
    //////////////////////////////////////////////////////////////////////////////////////////////
    assert_eq!(loci_chr.len(), loci_pos.len());
    let l = loci_chr.len();
    // Indices, chromosome names, and positions of the start and end of the window, respectively (will be filtered to remove redundant tails to remove windows which are complete subsets of a bigger window near the end of chromosomes or scaffolds)
    let mut idx_head: Vec<usize> = vec![0];
    let mut idx_tail = idx_head.clone();
    let mut chr_head: Vec<String> = vec![loci_chr[0].to_owned()];
    let mut chr_tail = chr_head.clone();
    let mut pos_head: Vec<u64> = vec![loci_pos[0]];
    let mut pos_tail = pos_head.clone();
    // Number of genotyped loci included per window
    let mut cov: Vec<u64> = vec![1];
    // A boolean to mark whether the start of the next window has been found before the end of the current window has been reached
    let mut marker_next_window_head: bool = false;
    // The index of the sart of the next window
    let mut idx_next_head: usize = 0;
    // Index of the current locus
    let mut i: usize = 1;
    // Counter for the number of polymorphic loci per window per pool (where we increment polymorphic site count per pool if the maximum allele frequency at a locus is less than 1 and assumes we are keeping all the alleles per locus)
    let mut polymorphic: Vec<Vec<u64>> =
        vec![polymorphic_loci_per_pool(genotypes_and_phenotypes, &loci_idx, 0).unwrap()];
    // We iterate across the genome until we reach the last locus (may not be consecutive as the start of the next window may be within the body of the current window)
    while i < l {
        let chr = loci_chr[i].to_owned();
        let pos = loci_pos[i];
        // println!("i={:?}", i);
        // println!("idx_head={:?}", idx_head);
        // println!("idx_tail={:?}", idx_tail);
        // Did we reach the end of the chromosome or the end of the window according to window size?
        if (&chr != chr_head.last().unwrap()) | (pos > (pos_head.last().unwrap() + window_size_bp))
        {
            // If we found the start of the next window in body of the current (ending) window then,
            //  we use the next window head as the start of the next slide not the end of the window,
            //  otherwise we use the end of the window.
            i = if marker_next_window_head {
                idx_next_head
            } else {
                i
            };
            let chr = loci_chr[i].to_owned();
            let pos = loci_pos[i];
            // Do we have the minimum number of required loci in the current window?
            if cov.last().unwrap() >= min_loci_per_window {
                // If we have enough loci covered in the current (ending) window:
                // We also add the details of the start of the next window
                idx_head.push(i);
                idx_tail.push(i);
                chr_head.push(chr.to_owned());
                chr_tail.push(chr.to_owned());
                pos_head.push(pos);
                pos_tail.push(pos);
                cov.push(1);
                polymorphic.push(
                    polymorphic_loci_per_pool(genotypes_and_phenotypes, &loci_idx, i).unwrap(),
                );
            } else {
                // If we did no have enough loci covered in the current (ending) window:
                // We ditch the current (ending) window and replace it with the start of the next window
                let i_ = idx_head.len() - 1;
                idx_head[i_] = i;
                chr_head[i_] = chr;
                pos_head[i_] = pos;
                cov[i_] = 1;
                polymorphic[i_] =
                    polymorphic_loci_per_pool(genotypes_and_phenotypes, &loci_idx, i_).unwrap();
            }
            // Reset the marker for the start of the next window
            marker_next_window_head = false;
        } else {
            // If we have yet to reach the end of the current window or the end of the chromosome,
            // then we just replace the tail of the current window with the current locus
            // and add the another locus to the coverage counter.
            let i_ = idx_tail.len() - 1;
            idx_tail[i_] = i;
            chr_tail[i_] = chr;
            pos_tail[i_] = pos;
            cov[i_] += 1;
            // Increment polymorphic site count per pool if the maximum allele frequency at a locus is less than 1 and assumes we are keeping all the alleles per locus
            let polymorphic_loci_per_pool =
                polymorphic_loci_per_pool(genotypes_and_phenotypes, &loci_idx, i_).unwrap();
            for k in 0..n {
                polymorphic[i_][k] += polymorphic_loci_per_pool[k];
            }
            // We also check if we have reached the start of the next window and note the index if we have
            if (marker_next_window_head == false)
                & (pos >= (pos_head.last().unwrap() + window_slide_size_bp))
            {
                marker_next_window_head = true;
                idx_next_head = i;
            }
        }
        // Move to the next locus
        i += 1;
    }
    // Remove redundant tails
    assert_eq!(idx_head.len(), idx_tail.len());
    let n_windows = idx_head.len();
    let mut out_idx_head: Vec<usize> = vec![idx_head[0]];
    let mut out_idx_tail: Vec<usize> = vec![idx_tail[0]];
    let mut out_cov: Vec<u64> = vec![cov[0]];
    let mut out_polymorphic: Vec<Vec<u64>> = vec![polymorphic[0].clone()];
    for i in 1..n_windows {
        // println!("out_idx_tail={:?}", out_idx_tail);
        if &idx_tail[i] != out_idx_tail.last().unwrap() {
            out_idx_head.push(idx_head[i]);
            out_idx_tail.push(idx_tail[i]);
            out_cov.push(cov[i]);
            out_polymorphic.push(polymorphic[i].clone());
        }
    }
    let n_windows = out_idx_head.len();
    // Check if we have no windows containing at least the minimum number of loci
    if out_cov.len() < 1 {
        let error_message =
            "No window with at least ".to_owned() + &min_loci_per_window.to_string() + " SNPs.";
        return Err(Error::new(ErrorKind::Other, error_message));
    }
    // Calculate Watterson's estimator per window
    let mut watterson_theta_per_pool_per_window: Array2<f64> =
        Array2::from_elem((n_windows, n), f64::NAN);
    for i in 0..n_windows {
        for j in 0..n {
            let n_segregating_sites = (out_polymorphic[i][j] as f64) / (out_cov[i] as f64);
            let correction_factor =
                (1..(pool_sizes[j] as usize)).fold(0.0, |sum, x| sum + 1.0 / (x as f64)); // Should we account for coverage instead of pool sizes?
            watterson_theta_per_pool_per_window[(i, j)] = n_segregating_sites / correction_factor;
        }
    }
    Ok((
        watterson_theta_per_pool_per_window,
        out_idx_head,
        out_idx_tail,
    ))
}

/// Estimate and save into a csv file
pub fn watterson_estimator(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    // Calculate Watterson's estimator
    let (watterson_theta_per_pool_per_window, windows_idx_head, windows_idx_tail) =
        theta_watterson(
            genotypes_and_phenotypes,
            pool_sizes,
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
        )
        .unwrap();
    // println!("watterson_theta_per_pool_per_window={:?}", watterson_theta_per_pool_per_window);
    let (n_windows, n) = watterson_theta_per_pool_per_window.dim();
    // println!("n={}; n_windows={}", n, n_windows);
    let vec_watterson_theta_across_windows = watterson_theta_per_pool_per_window
        .mean_axis(Axis(0))
        .unwrap();
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
            + "-watterson-"
            + &window_size_bp.to_string()
            + "_bp_windows-"
            + &time.to_string()
            + ".csv";
    }
    // Define the loci for writing the output
    let (_loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    // Instantiate output file
    let error_writing_file = "Unable to create file: ".to_owned() + &fname_output;
    let mut file_out = OpenOptions::new()
        .create_new(true)
        .write(true)
        .append(false)
        .open(&fname_output)
        .expect(&error_writing_file);
    // Header
    let mut line: Vec<String> = vec!["Pool".to_owned()];
    for i in 0..n_windows {
        let idx_ini = windows_idx_head[i];
        let idx_fin = windows_idx_tail[i];
        let window_chr = loci_chr[idx_ini].clone();
        let window_pos_ini = loci_pos[idx_ini];
        let window_pos_fin = loci_pos[idx_fin];
        line.push(
            "Window-".to_owned()
                + &window_chr
                + "_"
                + &window_pos_ini.to_string()
                + "_"
                + &window_pos_fin.to_string(),
        );
    }
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    // Write the nucleotide diversity per window
    for i in 0..n {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &vec_watterson_theta_across_windows[i].to_string()
            + ","
            + &watterson_theta_per_pool_per_window
                .column(i)
                .iter()
                .map(|x| parse_f64_roundup_and_own(*x, 8))
                .collect::<Vec<String>>()
                .join(",")
            + "\n";
        file_out.write_all(line.as_bytes()).unwrap();
    }
    Ok(fname_output)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_theta_w() {
        let x: Array2<f64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.6, 0.4, 0.0,
                0.9, 0.1, 1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5,
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
        println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
        // Outputs
        let out = watterson_estimator(
            &genotypes_and_phenotypes,
            &vec![42.0, 42.0, 42.0, 42.0, 42.0],
            &100, // 100-bp windows
            &50,  // 50-bp window slide
            &1,
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();
        let file = std::fs::File::open(&out).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut header: Vec<String> = vec![];
        let mut pop: Vec<String> = vec![];
        let mut watterson: Vec<f64> = vec![];
        for line in reader.lines() {
            let split = line
                .unwrap()
                .split(",")
                .map(|x| x.to_owned())
                .collect::<Vec<String>>();
            if header.len() == 0 {
                header = split;
            } else {
                pop.push(split[0].clone());
                for f in split[1..]
                    .iter()
                    .map(|x| x.parse::<f64>().unwrap())
                    .collect::<Vec<f64>>()
                {
                    watterson.push(f);
                }
            }
        }
        println!("watterson={:?}", watterson);
        let watterson: Array2<f64> = Array2::from_shape_vec((5, 3), watterson).unwrap();
        println!("watterson={:?}", watterson);
        let pop2_locus1 = watterson[(1, 1)]; // locus fixed, i.e. watterson=0.0
        let pop2_locus2 = watterson[(1, 2)]; // locus fixed, i.e. watterson=0.0
        let pop3_locus1 = watterson[(2, 1)]; // locus fixed, i.e. watterson=0.0
        let pop3_locus2 = watterson[(2, 2)]; // locus polymorphic with pool size of 42, i.e. watterson = 1 polymorphic locus / (42 pools - 1) = 0.2324
        assert_eq!(parse_f64_roundup_and_own(pop2_locus1, 4), "0".to_owned());
        assert_eq!(parse_f64_roundup_and_own(pop2_locus2, 4), "0".to_owned());
        assert_eq!(
            parse_f64_roundup_and_own(pop3_locus1, 4),
            "0.2324".to_owned()
        );
        assert_eq!(
            parse_f64_roundup_and_own(pop3_locus2, 4),
            "0.2324".to_owned()
        );
        // assert_eq!(0, 1);
    }
}
