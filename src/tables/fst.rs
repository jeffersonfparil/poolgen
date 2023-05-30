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

    let counts_per_pool_per_allele: Array2<f64> = locus_counts.matrix.map(|&y| y as f64);
    let total_counts_per_pool: Array1<f64> = counts_per_pool_per_allele.sum_axis(Axis(1));
    let (n, p) = counts_per_pool_per_allele.dim();

    let i = 213; // locus index
    let a = 0; // arbitrary reference allele index
    let j = 0; // pool i index
    let k = 1; // pool k index

    let q_1i_j = 1.00 - (2.00 * 
        (counts_per_pool_per_allele[(j,a)]/total_counts_per_pool[j]) * 
        ((total_counts_per_pool[j]-counts_per_pool_per_allele[(j,a)]) / (total_counts_per_pool[j] - 1.00))
    );
    
    let q_1i_k = 1.00 - (2.00 * 
        (counts_per_pool_per_allele[(k,a)]/total_counts_per_pool[k]) * 
        ((total_counts_per_pool[k]-counts_per_pool_per_allele[(k,a)]) / (total_counts_per_pool[k] - 1.00))
    );

    let q_2i_j = (1.00 / (total_counts_per_pool[j]*total_counts_per_pool[k])) * 
    ( counts_per_pool_per_allele[(j,a)]*counts_per_pool_per_allele[(k,a)] + 
        (total_counts_per_pool[j]-counts_per_pool_per_allele[(j,a)]) 
            * (total_counts_per_pool[k]-counts_per_pool_per_allele[(k,a)])
    );

    let fst_i_jk = ((0.5*(q_1i_j+q_1i_k)) - q_2i_j) / (1.00 + q_2i_j);
    println!("fst_i_jk = {:?}", fst_i_jk);

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
    Some(fst_i_jk.to_string())
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
            matrix: Array2::from_shape_vec((2, 4), vec![0, 0, 10, 0, 0, 10, 0, 0]).unwrap(),
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
        assert_eq!("1".to_owned(), out);
    }
}
