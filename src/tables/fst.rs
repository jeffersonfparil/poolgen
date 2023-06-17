use crate::base::*;
use ndarray::prelude::*;
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

pub fn fst(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<String> {
    let (n, p) = genotypes_and_phenotypes
        .intercept_and_allele_frequencies
        .dim();
    // Count the number of loci (Note: assumes the loci are sorted)
    let mut loci_idx: Vec<usize> = vec![];
    for i in 1..p {
        // includes the intercept
        if (genotypes_and_phenotypes.chromosome[i - 1] != genotypes_and_phenotypes.chromosome[i])
            | (genotypes_and_phenotypes.position[i - 1] != genotypes_and_phenotypes.position[i])
        {
            loci_idx.push(i);
        }
    }
    loci_idx.push(p); // last allele of the last locus
    let l = loci_idx.len();
    assert_eq!(l - 1, genotypes_and_phenotypes.coverages.ncols());
    // A probably not so unbiased and multi-allelic version of Gautier et al, 2019 (assumes biallelic loci)
    let mut fst: Array3<f64> = Array3::from_elem((l - 1, n, n), f64::NAN); // number of loci is loci_idx.len() - 1, i.e. less the last index - index of the last allele of the last locus
    let loci: Array3<usize> = Array3::from_shape_vec(
        (l - 1, n, n),
        (0..(l - 1))
            .flat_map(|x| std::iter::repeat(x).take(n * n))
            .collect(),
    )
    .unwrap();
    let pop1: Array3<usize> = Array3::from_shape_vec(
        (l - 1, n, n),
        std::iter::repeat(
            (0..n)
                .flat_map(|x| std::iter::repeat(x).take(n))
                .collect::<Vec<usize>>(),
        )
        .take(l - 1)
        .flat_map(|x| x)
        .collect::<Vec<usize>>(),
    )
    .unwrap();
    let pop2: Array3<usize> = Array3::from_shape_vec(
        (l - 1, n, n),
        std::iter::repeat(
            std::iter::repeat((0..n).collect::<Vec<usize>>())
                .take(n)
                .flat_map(|x| x)
                .collect::<Vec<usize>>(),
        )
        .take(l - 1)
        .flat_map(|x| x)
        .collect::<Vec<usize>>(),
    )
    .unwrap();
    // Parallel computations
    // Zip::from(&mut fst)
    //     .and(&loci)
    //     .and(&pop1)
    //     .and(&pop2)
    //     .par_for_each(|f, &i, &j, &k| {
    for i in 0..(l - 1) {
        for j in 0..n {
            for k in 0..n {
                if j != k {
                    let idx_start = loci_idx[i];
                    let idx_end = loci_idx[i + 1];
                    let g = genotypes_and_phenotypes
                        .intercept_and_allele_frequencies
                        .slice(s![.., idx_start..idx_end]);
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
                                                                                                    // *f = if f_unbiased < 0.0 {
                    fst[(i, j, k)] = if f_unbiased < 0.0 {
                        0.0
                    } else if f_unbiased > 1.0 {
                        1.0
                    } else {
                        f_unbiased
                    };
                } else {
                    // *f = 0.0; // Fst within itself is zero
                    fst[(i, j, k)] = 0.0; // Fst within itself is zero
                }
            }
        }
    }
    // });
    println!("loci={:?}", loci);
    println!("pop1={:?}", pop1);
    println!("pop2={:?}", pop2);
    println!("fst={:?}", fst);
    // println!("fst.mean_axis(Axis(0))={:?}", fst.mean_axis(Axis(0)));
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
        fname_output = bname.to_owned() + "-fst-" + &time.to_string() + ".csv";
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
    let fst_means: Array2<f64> = fst.mean_axis(Axis(0)).unwrap();
    for i in 0..n {
        let line = genotypes_and_phenotypes.pool_names[i].to_owned()
            + ","
            + &fst_means
                .row(i)
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
    fn test_fst() {
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
        // Outputs
        let out = fst(
            &genotypes_and_phenotypes,
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();
        let file = std::fs::File::open(&out).unwrap();
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
        let pop1_4 = fst[(0, 3)]; // the same populations, i.e. fst=0.0
        let pop4_1 = fst[(3, 0)]; // the same populations, i.e. fst=0.0
        let pop2_5 = fst[(1, 4)]; // totally different populations, i.e. fst=0.5, the same locus 1 and different locus 2
        let pop5_2 = fst[(4, 1)]; // totally different populations, i.e. fst=0.5, the same locus 1 and different locus 2
        println!("pop={:?}", pop);
        println!("fst={:?}", fst);
        // Assertions
        assert_eq!(diag, Array1::from_elem(pop.len(), 0.0));
        assert_eq!(pop1_4, 0.0);
        assert_eq!(pop4_1, 0.0);
        assert_eq!(pop2_5, 0.5);
        assert_eq!(pop5_2, 0.5);
        // assert_eq!(0, 1);
    }
}
