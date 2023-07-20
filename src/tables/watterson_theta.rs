use crate::base::*;
use ndarray::prelude::*;
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

/// Watterson's estimator of theta
/// Simply (naively?) defined as $theta_w = number of segregating sites / \Sigma^{n-1}_{i=1}(1/i)$
/// For details see [Feretti et al, 2013](https://doi.org/10.1111/mec.12522)
pub fn theta_watterson(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    window_size_bp: &usize,
) -> io::Result<(Array2<f64>, Vec<String>, Vec<u64>)> {
    let (n, p) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    assert_eq!(p, genotypes_and_phenotypes.chromosome.len(), "The number of entries in the 'chromosome' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
    assert_eq!(p, genotypes_and_phenotypes.position.len(), "The number of entries in the 'position' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
    // Count the number of loci (Note: assumes the loci are sorted) and extract the loci coordinates
    let mut windows_chr: Vec<String> = vec![genotypes_and_phenotypes.chromosome[1].to_owned()];
    let mut windows_pos: Vec<u64> = vec![genotypes_and_phenotypes.position[1] as u64];
    let mut window_n_sites: Vec<u64> = vec![1 as u64];
    let mut window_n_polymorphic_sites: Vec<Vec<u64>> = vec![vec![]];
    let idx = windows_chr.len() - 1;
    for j in 0..n {
        window_n_polymorphic_sites[idx].push(0);
        let freq = genotypes_and_phenotypes.intercept_and_allele_frequencies[(j, 0 + 1)];
        if (freq > 0.0) & (freq < 1.0) {
            window_n_polymorphic_sites[idx][j] += 1;
        }
    }
    for i in 2..p {
        if (genotypes_and_phenotypes.chromosome[i - 1] != genotypes_and_phenotypes.chromosome[i])
            | (genotypes_and_phenotypes.position[i - 1] != genotypes_and_phenotypes.position[i])
        {
            let chr = &genotypes_and_phenotypes.chromosome[i];
            let pos = &genotypes_and_phenotypes.position[i];
            // println!("windows_chr={:?}", windows_chr);
            // println!("windows_pos={:?}", windows_pos);
            // println!("chr={:?}", chr);
            // println!("pos={:?}", pos);
            if (chr != windows_chr.last().unwrap())
                | ((chr == windows_chr.last().unwrap())
                    & (pos > &(windows_pos.last().unwrap() + &(*window_size_bp as u64))))
            {
                windows_chr.push(chr.to_owned());
                windows_pos.push(*pos);
                window_n_sites.push(0);
                window_n_polymorphic_sites.push(vec![]);
                let idx = windows_chr.len() - 1;
                window_n_sites[idx] += 1;
                for j in 0..n {
                    window_n_polymorphic_sites[idx].push(0);
                    let freq = genotypes_and_phenotypes.intercept_and_allele_frequencies[(j, i)];
                    if (freq > 0.0) & (freq < 1.0) {
                        window_n_polymorphic_sites[idx][j] += 1;
                    }
                }
            }
        }
        let idx = windows_chr.len() - 1;
        window_n_sites[idx] += 1;
        for j in 0..n {
            let freq = genotypes_and_phenotypes.intercept_and_allele_frequencies[(j, i)];
            if (freq > 0.0) & (freq < 1.0) {
                window_n_polymorphic_sites[idx][j] += 1;
            }
        }
    }
    // Add the last index of the final position
    windows_chr.push(
        genotypes_and_phenotypes
            .chromosome
            .last()
            .unwrap()
            .to_owned(),
    );
    windows_pos.push(*genotypes_and_phenotypes.position.last().unwrap() - 1);
    // Calculate Watterson's estimator per window
    // println!("window_n_sites[0]={:?}", window_n_sites[0]);
    // println!("window_n_polymorphic_sites[0]={:?}", window_n_polymorphic_sites[0]);
    // println!("pool_sizes={:?}", pool_sizes);
    let n_windows = windows_chr.len() - 1;
    let mut watterson_theta_per_pool_per_window: Array2<f64> =
        Array2::from_elem((n_windows, n), f64::NAN);
    for i in 0..n_windows {
        for j in 0..n {
            let n_segregating_sites =
                (window_n_polymorphic_sites[i][j] as f64) / (window_n_sites[i] as f64);

            // Note: Do we need to account for coverage instead of pool sizes?
            let correction_factor =
                (1..(pool_sizes[j] as usize)).fold(0.0, |sum, x| sum + 1.0 / (x as f64));

            watterson_theta_per_pool_per_window[(i, j)] = n_segregating_sites / correction_factor;
        }
    }
    Ok((
        watterson_theta_per_pool_per_window,
        windows_chr,
        windows_pos,
    ))
}

/// Estimate and save into a csv file
pub fn watterson_estimator(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    window_size_bp: &usize,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    // Calculate Watterson's estimator
    let (watterson_theta_per_pool_per_window, windows_chr, windows_pos) =
        theta_watterson(genotypes_and_phenotypes, pool_sizes, window_size_bp).unwrap();
    // println!("watterson_theta_per_pool_per_window={:?}", watterson_theta_per_pool_per_window);
    let n = watterson_theta_per_pool_per_window.ncols();
    let n_windows = watterson_theta_per_pool_per_window.nrows();
    let vec_watterson_theta_across_windows = watterson_theta_per_pool_per_window.mean_axis(Axis(0)).unwrap();
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
        fname_output = bname.to_owned() + "-watterson-" + &time.to_string() + ".csv";
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
    let mut line: Vec<String> = vec!["Pool".to_owned()];
    for i in 0..n_windows {
        let window_chr = windows_chr[i].clone();
        let window_pos_ini = windows_pos[i];
        let window_pos_fin = window_pos_ini + (*window_size_bp as u64);
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
                .map(|x| parse_f64_roundup_and_own(*x, 4))
                .collect::<Vec<String>>()
                .join(",")
            + "\n";
        println!("line={:?}", line);
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
        // Outputs
        let out = watterson_estimator(
            &genotypes_and_phenotypes,
            &vec![42.0, 42.0, 42.0, 42.0, 42.0],
            &100, // 100-kb windows
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
        let watterson: Array2<f64> = Array2::from_shape_vec((5, 2), watterson).unwrap();
        println!("watterson={:?}", watterson);
        let pop2_locus1 = watterson[(1, 0)]; // locus fixed, i.e. watterson=0.0
        let pop2_locus2 = watterson[(1, 1)]; // locus fixed, i.e. watterson=0.0
        let pop3_locus1 = watterson[(2, 0)]; // locus fixed, i.e. watterson=0.0
        let pop3_locus2 = watterson[(2, 1)]; // locus at 0.5, i.e. watterson = 50 / (100-1) = 0.5051
        assert_eq!(pop2_locus1, 0.0);
        assert_eq!(pop2_locus2, 0.0);
        assert_eq!(pop3_locus1, 0.1549);
        assert_eq!(pop3_locus2, 0.2324);
        // assert_eq!(0, 2);
    }
}
