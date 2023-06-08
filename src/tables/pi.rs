use crate::base::*;
use ndarray::{prelude::*, Zip};
use std::fs::OpenOptions;
use std::io::{self, prelude::*, Error, ErrorKind};
use std::time::{SystemTime, UNIX_EPOCH};

pub fn pi(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    window_size_kb: &usize,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    let (n, p) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    assert_eq!(p, genotypes_and_phenotypes.chromosome.len(), "The number of entries in the 'chromosome' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
    assert_eq!(p, genotypes_and_phenotypes.position.len(), "The number of entries in the 'position' field and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'chromosome' fields of 'GenotypesAndPhenotypes' struct.");
    // Count the number of loci (Note: assumes the loci are sorted) and extract the loci coordinates
    let mut loci_idx: Vec<usize> = vec![];
    let mut loci_chr: Vec<&String> = vec![];
    let mut loci_pos: Vec<&u64> = vec![];
    for i in 1..p {
        // includes the intercept
        if (genotypes_and_phenotypes.chromosome[i - 1] != genotypes_and_phenotypes.chromosome[i])
            | (genotypes_and_phenotypes.position[i - 1] != genotypes_and_phenotypes.position[i])
        {
            loci_idx.push(i);
            loci_chr.push(&genotypes_and_phenotypes.chromosome[i]);
            loci_pos.push(&genotypes_and_phenotypes.position[i]);
        }
    }
    loci_idx.push(p); // last allele of the last locus
    let l = loci_idx.len();
    assert_eq!(l-1, genotypes_and_phenotypes.coverages.ncols(), "The number of loci with coverage information and the total number of loci are incompatible. Please check the 'intercept_and_allele_frequencies' and 'coverages' fields of 'GenotypesAndPhenotypes' struct.");
    // pi is the nucleotide diversity per population, and each pi across loci is oriented row-wise, i.e. each row is a single value across columns for each locus
    let mut pi: Array2<f64> = Array2::from_elem((l - 1, n), f64::NAN); // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    let loci: Array2<usize> = Array2::from_shape_vec(
        (l - 1, n),
        (0..(l - 1))
            .flat_map(|x| std::iter::repeat(x).take(n))
            .collect(),
    )
    .unwrap();
    let pop: Array2<usize> = Array2::from_shape_vec(
        (l - 1, n),
        std::iter::repeat((0..n)
        .collect::<Vec<usize>>())
        .take(l-1)
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
            *pi_ = ( (g.slice(s![j, ..]).fold(0.0, |sum, &x| sum + x.powf(2.0))
                * (nj / (nj - 1.00)))
                - (nj / (nj - 1.00)) ).abs(); // with a n/(n-1) factor on the heteroygosity to make it unbiased
        });
    // Summarize per non-overlapping window
    // Find window indices making sure we respect chromosomal boundaries
    let m = loci_idx.len() - 1; // total number of loci, we subtract 1 as the last index refer to the last allele of the last locus and serves as an end marker
    let mut windows_idx: Vec<usize> = vec![0]; // indices in terms of the number of loci not in terms of genome coordinates - just to make it simpler
    let mut windows_chr: Vec<String> = vec![loci_chr[0].to_owned()];
    let mut windows_pos: Vec<u64> = vec![loci_pos[0] + (*window_size_kb as u64)];
    for i in 1..m {
        let chr = loci_chr[i];
        let pos = loci_pos[i];
        if (chr != windows_chr.last().unwrap()) | ((chr == windows_chr.last().unwrap()) & (pos > windows_pos.last().unwrap())) {
            windows_idx.push(i);
            windows_chr.push(chr.to_owned());
            windows_pos.push(*pos);
        }
    }
    // Add the last index of the final position
    windows_idx.push(m);
    windows_chr.push(windows_chr.last().unwrap().to_owned());
    windows_pos.push(*loci_pos.last().unwrap() - 1);
    // Take the means per window
    let n_windows = windows_idx.len() - 1;
    let mut pi_per_pool_across_windows: Array2<f64> = Array2::from_elem((n_windows, n), f64::NAN);
    // let mut d_per_pool_across_windows: Array3<f64> = Array3::from_elem((n_windows, n, n), f64::NAN);
    for i in 0..n_windows {
        let idx_start = windows_idx[i];
        let idx_end = windows_idx[i + 1];
        for j in 0..n {
            pi_per_pool_across_windows[(i, j)] = pi
                .slice(s![idx_start..idx_end, j])
                .mean_axis(Axis(0))
                .unwrap()
                .fold(0.0, |_, &x| x);
        }
    }
    // println!(
    //     "pi_per_pool_across_windows={:?}",
    //     pi_per_pool_across_windows
    // );
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
        fname_output = bname.to_owned() + "-pi-" + &time.to_string() + ".csv";
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
        let window_pos_fin = windows_pos[i+1] + 1;
        line.push("Window-".to_owned() + &window_chr + "_" + &window_pos_ini.to_string() + "_" + &window_pos_fin.to_string());
    }
    let line = line.join(",") + "\n";
    file_out.write_all(line.as_bytes()).unwrap();
    // Write the nucleotide diversity per window
    for i in 0..n {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &pi
                .column(i)
                .iter()
                .map(|x| parse_f64_roundup_and_own(*x, 4))
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
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 
                1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 
                1.0, 0.6, 0.4, 0.0, 0.9, 0.1, 
                1.0, 0.4, 0.5, 0.1, 0.6, 0.4, 
                1.0, 1.0, 0.0, 0.0, 0.5, 0.5,
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
        let out = pi(
            &genotypes_and_phenotypes,
            &100, // 100-kb windows
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
        let pi: Array2<f64> = Array2::from_shape_vec((5, 2), pi).unwrap();
        let pop2_locus1 = pi[(1, 0)]; // locus fixed, i.e. pi=0.0
        let pop2_locus2 = pi[(1, 1)]; // locus fixed, i.e. pi=0.0
        let pop5_locus1 = pi[(4, 0)]; // locus fixed, i.e. pi=0.0
        let pop5_locus2 = pi[(4, 1)]; // locus at 0.5, i.e. pi = 50 / (100-1) = 0.5051
        assert_eq!(pop2_locus1, 0.0);
        assert_eq!(pop2_locus2, 0.0);
        assert_eq!(pop5_locus1, 0.0);
        assert_eq!(pop5_locus2, 0.5051);
        // assert_eq!(0, 1);
    }
}
