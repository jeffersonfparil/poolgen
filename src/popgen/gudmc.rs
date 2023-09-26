use crate::base::*;
use crate::popgen::*;
use ndarray::{prelude::*, Zip};
use std::fs::OpenOptions;
use std::io::{self, prelude::*};
use std::time::{SystemTime, UNIX_EPOCH};

/// gudmc: **g**enomewide **u**nbiased **d**etermination of the **m**odes of **c**onvergent evolution
/// using Tajima's D to detect peaks (balancing selection) and troughs (directional selection) as locations
pub fn gudmc(
    genotypes_and_phenotypes: &GenotypesAndPhenotypes,
    pool_sizes: &Vec<f64>,
    recombination_rate: &f64,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    fname_input: &String,
    fname_output: &String,
) -> io::Result<i32> {
    // Calculate Tajima's D
    let fname_tajima = tajima_d(
        genotypes_and_phenotypes,
        pool_sizes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
        fname_input,
        &"gudmc_intermediate_file_tajimasD.tmp".to_owned(),
    )
    .unwrap();

    // Calculate pairwise Fst (all pairwise combinations)
    let (_, fname_fst) = fst(
        genotypes_and_phenotypes,
        window_size_bp,
        window_slide_size_bp,
        min_loci_per_window,
        fname_input,
        &"gudmc_intermediate_file_Fst.tmp".to_owned(),
    )
    .unwrap();

    // Find significant troughs (selective sweeps) and peaks (balancing selection)
    //

    Ok(0)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_gudmc() {
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
        let out = gudmc(
            &genotypes_and_phenotypes,
            &vec![42.0, 42.0, 42.0, 42.0, 42.0],
            &3.14e-6,
            &100,
            &50,
            &1,
            &"test.something".to_owned(),
            &"".to_owned(),
        )
        .unwrap();

        // Assertions
        assert_eq!(out, 0);
        // assert_eq!(0, 1);
    }
}
