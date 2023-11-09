use crate::base::*;
use ndarray::prelude::*;
use std::io;
use statrs::distribution::{Uniform, MultivariateNormal, Continuous};
use rand::distributions::Distribution;
use rand::{seq::IteratorRandom, thread_rng};

// Simulate genotypes with some linkage disequilibrium

pub fn simulate_genotypes(n: usize, p: usize, n_chr: usize, max_bp: usize, r2_50_perc_bp: usize) -> io::Result<Array2<f64>> {
    // Define approximately equally sized chromosomes and coordinates of the characterised loci or markers or SNPs
    let mut chromosome_sizes = vec![p/n_chr; n_chr];
    let p_ = chromosome_sizes.iter().fold(0, |sum, &x| sum + x);
    // If the we have less or more than the required number of loci we add or subtract loci from the last chromosome
    if p_ < p {
        chromosome_sizes[n_chr-1] += p - p_;
    } else if p_ > p {
        chromosome_sizes[n_chr-1] -= p_ - p;
    }
    let max_bp_per_chr = max_bp/n_chr;
    println!("Chromosome sizes: {:?}", chromosome_sizes);
    // Sample uniformly distributed loci per chromosome
    let dist_unif = Uniform::new(0.0, max_bp_per_chr as f64).unwrap();
    println!("Initialised the uniform distribution");
    let mut rng = rand::thread_rng();
    let mut idx: usize = 0;
    let mut chromosomes: Array1<String> = Array1::from_elem(p, "".to_owned());
    let mut positions: Array1<usize> = Array1::from_elem(p, 0);
    let mut alleles: Array1<String> = Array1::from_elem(p, "".to_owned());
    let atcgd = vec!["A", "T", "C", "G", "D"].iter().map(|&x| x.to_owned()).collect::<Vec<String>>();
    for i in 0..n_chr {
        // Total number of sites we want to sample in the chromosome
        let s = chromosome_sizes[i];
        let mut tmp_positions: Vec<usize> = vec![];
        for _ in 0..s {
            tmp_positions.push(dist_unif.sample(&mut rng).floor() as usize);
        }
        tmp_positions.sort();
        // println!("tmp_positions={:?}", tmp_positions);
        for j in 0..s {
            chromosomes[idx] = "chr".to_owned() + &i.to_string();
            positions[idx] = tmp_positions[j];
            alleles[idx] = atcgd.iter().choose_multiple(&mut rng, 1)[0].to_owned();
            idx += 1;
        }
    }
    // println!("chromosomes:\n{:?}", chromosomes);
    // println!("positions:\n{:?}", positions);
    // println!("alleles:\n{:?}", alleles);

    let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();


    return Ok(Array2::from_elem((1, 1), f64::NAN))
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_simulate_genotypes() {
        let n = 100;
        let p = 10_000;
        let n_chr = 7;
        let max_bp = 2.2e9 as usize;
        let r2_50_perc_bp = 10e6 as usize;
        let q = simulate_genotypes(n, p, n_chr, max_bp, r2_50_perc_bp).unwrap();
        assert_eq!(0, 1);
    }
}
