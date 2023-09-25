use crate::base::*;
use crate::popgen::*;
use ndarray::prelude::*;
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

/// Tajima's D
pub fn tajima_d(
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
    let n_pools = watterson_theta_per_pool_per_window.ncols();
    let n_windows = watterson_theta_per_pool_per_window.nrows();
    // Calculate heterozygosities
    let (pi_per_pool_per_window, windows_idx_head_pi, windows_idx_tail_pi) = theta_pi(
        genotypes_and_phenotypes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
    )
    .unwrap();
    // println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
    // println!("watterson_theta_per_pool_per_window={:?}", watterson_theta_per_pool_per_window);
    // println!("pi_per_pool_per_window={:?}", pi_per_pool_per_window);
    // Sanity checks
    assert_eq!(n_pools, pi_per_pool_per_window.ncols(), "The number of pools extracted from estimating the heterozygosities and Watterson's estimators are incompatible. Please report a bug.");
    assert_eq!(n_windows, pi_per_pool_per_window.nrows(), "The number of windows extracted from estimating the heterozygosities and Watterson's estimators are incompatible. Please report a bug.");
    assert_eq!(windows_idx_head, windows_idx_head_pi, "The chromosome names per window extracted from estimating the heterozygosities and Watterson's estimators are incompatible. Please report a bug.");
    assert_eq!(windows_idx_tail, windows_idx_tail_pi, "The SNP positions per window extracted from estimating the heterozygosities and Watterson's estimators are incompatible. Please report a bug.");
    // Calculate Tajima's D
    let mut tajimas_d_per_pool_per_window: Array2<f64> =
        Array2::from_elem((n_windows, n_pools), f64::NAN);
    for j in 0..n_pools {
        let n = pool_sizes[j] as usize;
        let a1 = (1..n).fold(0.0, |sum, x| sum + (1.00 / (x as f64)));
        let a2 = (1..n).fold(0.0, |sum, x| sum + (1.00 / (x as f64).powf(2.00)));
        let n = n as f64;
        let b1 = (n + 1.0) / (3.0 * (n - 1.0));
        let b2 = (2.0 * (n.powf(2.0) + n + 3.0)) / (9.0 * n * (n - 1.0));
        let c1 = b1 - (1.0 / a1);
        let c2 = b2 - ((n + 2.0) / (a1 * n)) + (a2 / a1.powf(2.0));
        let e1 = c1 / a1;
        let e2 = c2 / (a1.powf(2.0) + a2);
        for i in 0..n_windows {
            let s = if watterson_theta_per_pool_per_window[(i, j)] <= f64::EPSILON {
                0.0
            } else {
                watterson_theta_per_pool_per_window[(i, j)] / a1
            };
            let vd = e1 * s + e2 * s * (s - 1.0);
            tajimas_d_per_pool_per_window[(i, j)] = if (pi_per_pool_per_window[(i, j)]
                - watterson_theta_per_pool_per_window[(i, j)])
                .abs()
                <= f64::EPSILON
            {
                0.0
            } else if vd <= f64::EPSILON {
                // If heterozygosity - expected and observed are very low in the window set Tajima's D to zero to avoid infinities
                0.0
            } else {
                (pi_per_pool_per_window[(i, j)] - watterson_theta_per_pool_per_window[(i, j)])
                    / vd.sqrt()
            };
            // if tajimas_d_per_pool_per_window[(i, j)].is_infinite() {
            //     println!("pi_per_pool_per_window[(i, j)]={:?}", pi_per_pool_per_window[(i, j)]);
            //     println!("watterson_theta_per_pool_per_window[(i, j)]={:?}", watterson_theta_per_pool_per_window[(i, j)]);
            //     println!("vd={:?}", vd);
            //     println!("a1={:?}", a1);
            //     println!("a2={:?}", a2);
            //     println!("b1={:?}", b1);
            //     println!("b2={:?}", b2);
            //     println!("c1={:?}", c1);
            //     println!("c2={:?}", c2);
            //     println!("e1={:?}", e1);
            //     println!("e2={:?}", e2);
            // }
        }
    }
    // Mean Tajima's D across all the windows
    let vec_tajimas_d_across_windows = tajimas_d_per_pool_per_window.mean_axis(Axis(0)).unwrap();
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
            + "-Tajimas_D-"
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
    let mut line: Vec<String> = vec!["Pool".to_owned(), "Mean_across_windows".to_owned()];
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
    for i in 0..n_pools {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &vec_tajimas_d_across_windows[i].to_string()
            + ","
            + &tajimas_d_per_pool_per_window
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
    fn test_tajima_d() {
        let x: Array2<f64> = Array2::from_shape_vec(
            (5, 6),
            vec![
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.6, 0.4, 0.0,
                0.9, 0.1, 1.0, 0.01, 0.01, 0.98, 0.6, 0.4, 1.0, 1.0, 0.0, 0.0, 0.5, 0.5,
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
        let out = tajima_d(
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
        let mut d: Vec<f64> = vec![];
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
                    d.push(f);
                }
            }
        }
        let d: Array2<f64> = Array2::from_shape_vec((5, 3), d).unwrap();
        println!("d={:?}", d);
        let pop2_locus1 = d[(1, 1)]; // locus fixed - neutral, i.e. d=0.0
        let pop2_locus2 = d[(1, 2)]; // locus fixed - neutral, i.e. d=0.0
        let pop4_locus1 = d[(3, 1)]; // excess rare alleles - selective sweep, i.e. d<0.0
        let pop4_locus2 = d[(3, 2)]; // scarce rare alleles - balancing selection, i.e. d>0.0
        assert_eq!(parse_f64_roundup_and_own(pop2_locus1, 4), "0".to_owned());
        assert_eq!(parse_f64_roundup_and_own(pop2_locus2, 4), "0".to_owned());
        assert_eq!(
            parse_f64_roundup_and_own(pop4_locus1, 4),
            "-5.3954".to_owned()
        );
        assert_eq!(
            parse_f64_roundup_and_own(pop4_locus2, 4),
            "7.072".to_owned()
        );
        // assert_eq!(0, 1);
    }
}
