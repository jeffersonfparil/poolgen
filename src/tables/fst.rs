use crate::base::*;
use ndarray::prelude::*;
use statrs::distribution::{ChiSquared, ContinuousCDF};

pub fn fst(locus_counts: &mut LocusCounts, filter_stats: &FilterStats) -> Option<String> {
    println!("locus_counts={:?}", locus_counts);
    let locus_counts = match locus_counts.filter(filter_stats) {
        Ok(x) => x,
        Err(_) => return None,
    };
    let locus_frequencies = match locus_counts.to_frequencies() {
        Ok(x) => x,
        Err(_) => return None,
    };

    // Using Gautier et al, 2021 estimates (https://doi.org/10.1111/1755-0998.13557)
    println!("locus_frequencies={:?}", locus_frequencies);

    // Let's start with pairwise Fst
    let i_a: usize = 0; // pool a
    let i_b: usize = 1; // pool b
    let j: usize = 0; // allele 0
    let q1_a = locus_frequencies.matrix[(i_a, j)]; // probability of sampling two genes (or alleles) identical in state (IIS) within population a
    let q1_b = locus_frequencies.matrix[(i_b, j)]; // probability of sampling two genes (or alleles) IIS within population b
    let q1 = 0.5 * (q1_a + q1_b);
    let q2_ab = q1_a * q1_b; // probability of sampling two IIS genes from a and b
    let fst_ab = (q1 - q2_ab) / (1.00 - q2_ab);

    /// Apparently unbiasedestimates of q1 and q2
    let r: Array2<f64> = locus_counts.matrix.map(|&x| x as f64);
    let c: Array1<f64> = r.sum_axis(Axis(0));

    let p1 = 0; // population 1 index
    let p2 = 1; // population 2 index
    let a = 0; // allele index
    let q1 = 1.00
        - 2.00
            * ((c[a] / 2.00) / (c[a] / 2.00 - 1.00))
            * (r[(p1, a)] / c[a])
            * ((c[a] - r[(p1, a)]) / (c[a] - 1.00));
    
    let q2 = (r[(p1,a)]*r[(p2,a)] + (c[a]-r[(p1,a)]) * (c[a]-r[(p2,a)])) / (c[a]*c[a]);

    let r_ij: Array2<f64> = locus_counts.matrix.map(|&x| x as f64);
    let c_ij: Array1<f64> = r_ij.sum_axis(Axis(1));
    let n_j = r_ij.sum();

    println!("r_ij={:?}", r_ij);
    println!("c_ij={:?}", c_ij);
    println!("fst_ab={:?}", fst_ab);

    // Output line
    // let out = vec![
    //     locus_frequencies.chromosome.clone(),
    //     locus_frequencies.position.to_string(),
    //     locus_frequencies.alleles_vector.join(""),
    //     parse_f64_roundup_and_own(chi2, 6),
    //     pval.to_string(),
    // ]
    // .join(",")
    //     + "\n";
    // Some(out)
    Some("".to_owned())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_fst() {
        // Expected
        let expected_line = "Chromosome1,12345,AT,4,0.7797774084757156\n".to_owned(); // where df=7, then the pvalue is calculated as the lower tail because if chi2 < df
                                                                                      // Inputs
        let mut locus_counts = LocusCounts {
            chromosome: "Chromosome1".to_owned(),
            position: 12345,
            // alleles_vector: vec!["A", "T"]
            alleles_vector: vec!["A", "T", "C", "D"]
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>(),
            // matrix: Array2::from_shape_vec((4, 2), vec![0, 20, 20, 0, 0, 20, 20, 0]).unwrap(),
            matrix: Array2::from_shape_vec((2, 4), vec![0, 1, 10, 0, 0, 10, 1, 0]).unwrap(),
        };
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.01,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            // pool_sizes: vec![0.2, 0.2, 0.2, 0.2],
            pool_sizes: vec![0.5, 0.5],
        };
        // Outputs
        let out = fst(&mut locus_counts, &filter_stats).unwrap();
        // Assertions
        assert_eq!("0".to_owned(), out);
    }
}
