use crate::base::*;
use ndarray::{prelude::*, Zip};
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

/// Unbiased multi-allelic version of Fst similar to [Gautier et al, 2019](https://doi.org/10.1111/1755-0998.13557) which assumes biallelic loci
/// Note: Window sizes smaller than the requested window size can be expected because the actual window sizes noted are dependent on the minimum coverage per window, such that window sizes are based on the coordinates of the loci covered
pub fn fst(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<(String, String)> {
    let (n, _) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    let (loci_idx, loci_chr, loci_pos) = genotypes_and_phenotypes.count_loci().unwrap();
    let l = loci_idx.len() - 1; // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    let mut fst: Array3<f64> = Array3::from_elem((l, n, n), f64::NAN);
    let loci: Array3<usize> = Array3::from_shape_vec(
        (l, n, n),
        (0..(l))
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
    // Parallel computations (NOTE: Not efficient yet. Compute only the upper or lower triangular in the future.)
    Zip::from(&mut fst)
        .and(&loci)
        .and(&pop1)
        .and(&pop2)
        .par_for_each(|f, &i, &j, &k| {
            let idx_start = loci_idx[i];
            let idx_end = loci_idx[i + 1];
            let g = genotypes_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![.., idx_start..idx_end]);
            // Make sure that allele frequencies per locus sums up to one
            assert!(g.sum_axis(Axis(1)).sum() as usize == n);
            let nj = genotypes_and_phenotypes.coverages[(j, i)];
            let nk = genotypes_and_phenotypes.coverages[(k, i)];
            let q1_j = (g.slice(s![j, ..]).fold(0.0, |sum, &x| sum + x.powf(2.0))
                * (nj / (nj - 1.00 + f64::EPSILON)))
                + (1.00 - (nj / (nj - 1.00 + f64::EPSILON))); // with a n/(n-1) factor on the heteroygosity to make it unbiased
            let q1_k = (g.slice(s![k, ..]).fold(0.0, |sum, &x| sum + x.powf(2.0))
                * (nk / (nk - 1.00 + f64::EPSILON)))
                + (1.00 - (nk / (nk - 1.00 + f64::EPSILON))); // with a n/(n-1) factor on the heteroygosity to make it unbiased
            let q2_jk = g
                .slice(s![j, ..])
                .iter()
                .zip(&g.slice(s![k, ..]))
                .fold(0.0, |sum, (&x, &y)| sum + (x * y));
            let f_unbiased = (0.5 * (q1_j + q1_k) - q2_jk) / (1.00 - q2_jk + f64::EPSILON); // The reason why we're getting NaN is that q2_jk==1.0 because the 2 populations are both fixed at the locus, i.e. the same allele is at 1.00.
            *f = if f_unbiased < 0.0 {
                // fst[(i, j, k)] = if f_unbiased < 0.0 {
                0.0
            } else if f_unbiased > 1.0 {
                1.0
            } else {
                f_unbiased
            };
        });
    // Write output (Fst averaged across all loci)
    let mut fname_output = fname_output.to_owned();
    let mut fname_output_per_window = fname_output
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
    fname_output_per_window =
        fname_output_per_window + "-fst-" + &window_size_bp.to_string() + "_bp_windows.csv";
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
        fname_output =
            bname.to_owned() + "-fst-averaged_across_genome-" + &time.to_string() + ".csv";
        fname_output_per_window = bname.to_owned()
            + "-fst-"
            + &window_size_bp.to_string()
            + "_bp_windows-"
            + &time.to_string()
            + ".csv";
    }
    // Instantiate output file
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
    let fst_means: Array2<f64> = fst.mean_axis(Axis(0)).unwrap();
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
    // Define sliding windows
    let mut loci_chr_no_redundant_tail = loci_chr.to_owned();
    loci_chr_no_redundant_tail.pop();
    let mut loci_pos_no_redundant_tail = loci_pos.to_owned();
    loci_pos_no_redundant_tail.pop();
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
    let mut fst_per_pool_x_pool_per_window: Array2<f64> =
        Array2::from_elem((n_windows, n * n), f64::NAN);
    for i in 0..n_windows {
        let idx_start = windows_idx_head[i];
        let idx_end = windows_idx_tail[i] + 1; // add one so that we include the tail index as part of the window
        for j in 0..n {
            for k in 0..n {
                let idx = (j * n) + k;
                // println!("#############################");
                // println!("start={}; end={}; j={}; k={}", idx_start, idx_end, j, j);
                // println!("fst.slice(s![idx_start..idx_end, j, k])={:?}", fst.slice(s![idx_start..idx_end, j, k]));
                fst_per_pool_x_pool_per_window[(i, idx)] =
                    match fst.slice(s![idx_start..idx_end, j, k]).mean_axis(Axis(0)) {
                        Some(x) => x.fold(0.0, |_, &x| x),
                        None => f64::NAN,
                    };
            }
        }
    }
    // println!(
    //     "fst_per_pool_x_pool_per_window={:?}",
    //     fst_per_pool_x_pool_per_window
    // );
    // Write output (Fst averaged per window)
    // Instantiate output file
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
        let idx_head = windows_idx_head[i];
        let idx_tail = windows_idx_tail[i];
        let chr = loci_chr[idx_head].to_owned();
        let pos_head = loci_pos[idx_head];
        let pos_tail = loci_pos[idx_tail];
        let coordinate = chr + "," + &pos_head.to_string() + "," + &pos_tail.to_string() + ",";
        let fst_string = fst_per_pool_x_pool_per_window
            .slice(s![i, ..])
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join(",");
        let line = coordinate + &fst_string[..] + "\n";
        file_out.write_all(line.as_bytes()).unwrap();
    }
    Ok((fname_output, fname_output_per_window))
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_fst() {
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
            start_index_of_each_locus: vec![0, 1, 4],
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
        let (out_genomewide, out_per_window) = fst(
            &genotypes_and_phenotypes,
            &100,
            &50,
            &1,
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();
        let file = std::fs::File::open(&out_genomewide).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut header: Vec<String> = vec![];
        let mut pop: Vec<String> = vec![];
        let mut fst: Vec<f64> = vec![];
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
                    fst.push(f);
                }
            }
        }
        let fst: Array2<f64> = Array2::from_shape_vec((pop.len(), pop.len()), fst).unwrap();
        let diag: Array1<f64> = fst.diag().to_owned(); // itself, i.e. fst=0.0
        let pop1_2 = fst[(0, 1)]; // the same populations, i.e. fst=0.0
        let pop2_1 = fst[(1, 0)]; // the same populations, i.e. fst=0.0
        let pop4_5 = fst[(3, 4)]; // totally different populations, i.e. fst=1.0
        let pop5_4 = fst[(4, 3)]; // totally different populations, i.e. fst=1.0
        let pop1_3 = fst[(0, 2)]; // fst ~ 0.5
        let pop3_1 = fst[(2, 1)]; // fst ~ 0.5
        println!("genotypes_and_phenotypes={:?}", genotypes_and_phenotypes);
        println!("pop={:?}", pop);
        println!("fst={:?}", fst);
        println!("out_per_window={:?}", out_per_window);
        // Assertions
        assert_eq!(diag, Array1::from_elem(pop.len(), 0.0));
        assert_eq!(pop1_2, 1.0);
        assert_eq!(pop2_1, 1.0);
        assert_eq!(pop4_5, 0.0);
        assert_eq!(pop5_4, 0.0);
        assert!(f64::abs(pop1_3 - 0.5) < 0.1);
        assert!(f64::abs(pop3_1 - 0.5) < 0.1);
        // assert_eq!(0, 1);
    }
}
