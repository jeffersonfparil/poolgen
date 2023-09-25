use crate::base::*;
use ndarray::{prelude::*, Zip};
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::io::{Error, ErrorKind};
use std::time::{SystemTime, UNIX_EPOCH};

/// Unbiased multi-allelic nucleotide diversity per population ($\pi$ or $\theta_{\pi}=4N_{e}\mu$), which is similar to [Korunes & Samuk 2019](https://doi.org/10.1111/1755-0998.13326) which assumes biallelic loci
pub fn theta_pi(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
) -> io::Result<(Array2<f64>, Vec<usize>, Vec<usize>)> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let (loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    let l = loci_idx.len() - 1; // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
                                // Each pi across loci is oriented row-wise, i.e. each row is a single value across columns for each locus
    let mut pi: Array2<f64> = Array2::from_elem((l, n), f64::NAN);
    let loci: Array2<usize> = Array2::from_shape_vec(
        (l, n),
        (0..(l))
            .flat_map(|x| std::iter::repeat(x).take(n))
            .collect(),
    )
    .unwrap();
    let pop: Array2<usize> = Array2::from_shape_vec(
        (l, n),
        std::iter::repeat((0..n).collect::<Vec<usize>>())
            .take(l)
            .flat_map(|x| x)
            .collect(),
    )
    .unwrap();
    // Parallel computations
    Zip::from(&mut pi)
        .and(&loci)
        .and(&pop)
        .par_for_each(|pi_, &i, &j| {
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            let nj = genotypes_and_phenotypes.coverages[(j, i)];
            // Nucleotide diversity (~ heterozygosity), where population across rows which means each column is the same value
            *pi_ = ((g.slice(s![j, ..]).fold(0.0, |sum, &x| sum + x.powf(2.0))
                * (nj / (nj - 1.00 + f64::EPSILON)))
                - (nj / (nj - 1.00 + f64::EPSILON)))
                .abs(); // equivalent to (n/(n-1))*(1-sum(p^2)) with a n/(n-1) factor on the heteroygosity to make it unbiased
        });
    // Summarize per non-overlapping window
    // Remove redundant trailing loci from `genotypes_and_phenotypes.count_loci()`
    let mut loci_chr_no_redundant_tail = loci_chr.to_owned();
    loci_chr_no_redundant_tail.pop();
    let mut loci_pos_no_redundant_tail = loci_pos.to_owned();
    loci_pos_no_redundant_tail.pop();
    // println!("loci_idx={:?}", loci_idx);
    // println!("loci_chr_no_redundant_tail={:?}", loci_chr_no_redundant_tail);
    // println!("loci_pos_no_redundant_tail={:?}", loci_pos_no_redundant_tail);
    // Define sliding windows
    let (windows_idx_head, windows_idx_tail) = define_sliding_windows(
        &loci_chr_no_redundant_tail,
        &loci_pos_no_redundant_tail,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
    )
    .unwrap();
    // println!("fst={:?}", fst);
    // println!("l={:?}", l);
    // println!("windows_idx_head={:?}", windows_idx_head);
    // println!("windows_idx_tail={:?}", windows_idx_tail);
    // Take the means per window
    let n_windows = windows_idx_head.len();
    assert!(n_windows > 0, "There were no windows defined. Please check the sync file, the window size, slide size, and the minimum number of loci per window.");
    let mut pi_per_pool_per_window: Array2<f64> = Array2::from_elem((n_windows, n), f64::NAN);
    for i in 0..n_windows {
        let idx_start = windows_idx_head[i];
        let idx_end = windows_idx_tail[i] + 1; // add one so that we include the tail index as part of the window
        for j in 0..n {
            let idx = j;
            pi_per_pool_per_window[(i, idx)] =
                match pi.slice(s![idx_start..idx_end, j]).mean_axis(Axis(0)) {
                    Some(x) => x.fold(0.0, |_, &x| x),
                    None => f64::NAN,
                };
        }
    }
    // println!("windows_idx_head={:?}", windows_idx_head);
    // println!("windows_idx_tail={:?}", windows_idx_tail);
    // println!("pi_per_pool_per_window={:?}", pi_per_pool_per_window);
    Ok((pi_per_pool_per_window, windows_idx_head, windows_idx_tail))
}

pub fn pi(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    // Calculate heterozygosities
    let (pi_per_pool_per_window, windows_idx_head, windows_idx_tail) = theta_pi(
        genotypes_and_phenotypes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
    )
    .unwrap();
    let n = pi_per_pool_per_window.ncols();
    let n_windows = pi_per_pool_per_window.nrows();
    assert!(n_windows==windows_idx_head.len(), "Please check the number of windows in the pi estimates and the starting indices of each window.");
    assert!(n_windows==windows_idx_tail.len(), "Please check the number of windows in the pi estimates and the ending indices of each window.");
    let vec_pi_across_windows = pi_per_pool_per_window.mean_axis(Axis(0)).unwrap();
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
            + "-pi-"
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
    for i in 0..n {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &vec_pi_across_windows[i].to_string()
            + ","
            + &pi_per_pool_per_window
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
    fn test_pi() {
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
        let out = pi(
            &genotypes_and_phenotypes,
            &100, // 100-bp windows
            &50,  // 50-bp window slide
            &1,   // minimum of 1 SNP per window
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();
        let file = std::fs::File::open(&out).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut header: Vec<String> = vec![];
        let mut pop: Vec<String> = vec![];
        let mut pi: Vec<f64> = vec![];
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
                    pi.push(f);
                }
            }
        }
        let pi: Array2<f64> = Array2::from_shape_vec((5, 3), pi).unwrap();
        let pop2_locus1 = pi[(1, 1)]; // locus fixed, i.e. pi=0.0
        let pop2_locus2 = pi[(1, 2)]; // locus fixed, i.e. pi=0.0
        let pop5_locus1 = pi[(4, 1)]; // locus fixed, i.e. pi=0.0
        let pop5_locus2 = pi[(4, 2)]; // locus at 0.5, i.e. pi = 50 / (100-1) = 0.5051
        assert_eq!(parse_f64_roundup_and_own(pop2_locus1, 4), "0".to_owned());
        assert_eq!(parse_f64_roundup_and_own(pop2_locus2, 4), "0".to_owned());
        assert_eq!(parse_f64_roundup_and_own(pop5_locus1, 4), "0".to_owned());
        assert_eq!(
            parse_f64_roundup_and_own(pop5_locus2, 4),
            "0.5051".to_owned()
        );
        // assert_eq!(0, 1);
    }
}
