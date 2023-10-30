use crate::base::*;
use crate::gwas::pearsons_correlation;
use ndarray::{prelude::*, Zip};
use std::io::{self, Error, ErrorKind};

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation(
        &mut self,
        window_size_bp: &u64,
        window_slide_size_bp: &u64,
        min_loci_per_window: &u64,
        min_correlation: &f64,
        k_neighbours: &u64,
    ) -> io::Result<&mut Self> {
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        let (n_, l) = self.coverages.dim();
        let l_ = self.start_index_of_each_locus.len();
        let (idx_window_head, idx_window_tail) = define_sliding_windows(
            &self.chromosome,
            &self.position,
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
        )
        .unwrap();
        // println!("idx_window_head={:?}; idx_window_tail={:?}", idx_window_head, idx_window_tail);
        let w = idx_window_head.len();
        let w_ = idx_window_tail.len();
        assert_eq!(n, n_);
        assert_eq!(l, l_);
        assert_eq!(w, w_);

        // Parallel processing per window
        let mut vec_windows_freqs: Vec<Array2<f64>> = vec![Array2::from_elem((1, 1), f64::NAN); w];
        let idx_windows: Vec<usize> = (0..w).collect();
        // Zip::from(&mut vec_windows_freqs)
        //     .and(&idx_windows)
        //     .par_for_each(|window_freqs, &a| {
        for a in 0..idx_windows.len() {
                let idx_ini = idx_window_head[a];
                let idx_fin = idx_window_tail[a];
                let p = idx_fin - idx_ini;
                // *window_freqs = self
                let mut window_freqs = self
                    .intercept_and_allele_frequencies
                    .slice(s![.., idx_ini..idx_fin])
                    .to_owned();
                // Calculate correlations between alleles across all loci within the window (calculate only for the upper tringular and just mirror the results into the lower tringular for efficiency)
                let mut corr: Array2<f64> = Array2::from_elem((p, p), f64::NAN);
                for j0 in 0..p {
                    for j1 in j0..p {
                        match pearsons_correlation(
                            &window_freqs.column(j0),
                            &window_freqs.column(j1),
                            &"sensible_corr".to_owned(),
                        ) {
                            Ok(x) => {
                                corr[(j0, j1)] = x.0;
                                corr[(j1, j0)] = x.0
                            }
                            Err(_) => continue,
                        };
                    }
                }
                println!("corr:\n{:?}", corr);
                for j in 0..p {
                    if window_freqs
                        .column(j)
                        .fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum })
                        > 0
                    {
                        // Find the indexes of the linked (positively correlated) alleles to be used in calculating the distances between pools to find the nearest neighbours
                        let idx_linked_alleles = corr
                            .column(j)
                            .iter()
                            .enumerate()
                            .filter(|&(_i, x)| !x.is_nan())
                            .filter(|&(_i, x)| x >= min_correlation)
                            .map(|(i, _x)| i)
                            .collect::<Vec<usize>>();
                        // Use all the alleles to estimate distances between pools to find the nearest neighbours
                        let idx_linked_alleles = if idx_linked_alleles.len() < 2 {
                            (0..p).collect()
                        } else {
                            idx_linked_alleles
                        };
                        // Calculate the Euclidean distances between pools using the most positively correlated alleles (or all the alleles if none passed the minimum correlation threshold or if there is less that 2 loci that are most correlated - these minimum 2 is the allele itself and another allele)
                        let mut dist: Array2<f64> = Array2::from_elem((n, n), f64::NAN);
                        let window_freqs_linked_alleles = window_freqs.select(Axis(1), &idx_linked_alleles);
                        for i0 in 0..n {
                            let pool0 = window_freqs_linked_alleles
                                .row(i0);
                            // println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
                            // println!("idx_linked_alleles:\n{:?}", &idx_linked_alleles);
                            // println!("pool0:\n{:?}", &pool0);
                            for i1 in i0..n {
                                let pool1 = window_freqs_linked_alleles
                                    .row(i1)
                                    .to_owned();
                                // Keep only the loci present in both pools
                                let (pool0, pool1): (Vec<f64>, Vec<f64>) = pool0
                                    .iter()
                                    .zip(pool1.iter())
                                    .filter(|&(&x, &y)| (!x.is_nan()) & (!y.is_nan()))
                                    .unzip();
                                if pool0.len() == 0 {
                                    continue;
                                } else {
                                    let d = (&Array1::from_vec(pool0) - &Array1::from_vec(pool1))
                                        .map(|&x| x.powf(2.0))
                                        .sum()
                                        .sqrt();
                                    match d.is_nan() {
                                        true => continue,
                                        false => {
                                            dist[(i0, i1)] = d;
                                            dist[(i1, i0)] = d
                                        }
                                    };
                                }
                            }
                        }
                        println!("dist={:?}", dist);
                        // Now, let's find the pools needing imputation and impute them using the k-nearest neighbours
                        for i in 0..n {
                            let mut k = *k_neighbours as usize; // define the adaptive k for use later
                            if window_freqs[(i, j)].is_nan() {
                                // Find the k-nearest neighbours
                                let mut idx_pools: Vec<usize> = (0..n).collect();
                                idx_pools.sort_by(|&j0, &j1| {
                                    let a = match dist[(i, j0)].is_nan() {
                                        true => f64::INFINITY,
                                        false => dist[(i, j0)],
                                    };
                                    let b = match dist[(i, j1)].is_nan() {
                                        true => f64::INFINITY,
                                        false => dist[(i, j1)],
                                    };
                                    a.partial_cmp(&b).unwrap()
                                });
                                // println!("dist.column(i):\n{:?}", dist.column(i));
                                // println!("idx_pools:\n{:?}", idx_pools);
                                let freqs_sorted_neighbours = window_freqs
                                    .select(Axis(0), &idx_pools)
                                    .column(j)
                                    .to_owned();
                                // println!("freqs_sorted_neighbours:\n{:?}", freqs_sorted_neighbours);
                                // println!("freqs_unsorted_neighbours:\n{:?}", self.intercept_and_allele_frequencies.column(idx_ini + j));
                                // Extract the allele frequencies and distances of the k-neighbours from sorted lists using idx_pools whose indexes are sorted
                                let mut freqs_k_neighbours = freqs_sorted_neighbours
                                    .select(Axis(0), &(0..k).collect::<Vec<usize>>());
                                let mut dist_k_neighbours = dist
                                    .column(i)
                                    .select(Axis(0), &idx_pools)
                                    .select(Axis(0), &(0..k).collect::<Vec<usize>>());
                                // println!("freqs_k_neighbours:\n{:?}", freqs_k_neighbours);
                                // println!("dist_k_neighbours:\n{:?}", dist_k_neighbours);
                                while freqs_k_neighbours.fold(0, |sum, &x| {
                                    if x.is_nan() {
                                        sum + 1
                                    } else {
                                        sum
                                    }
                                }) == k
                                {
                                    k += 1;
                                    if k > n {
                                        break;
                                    }
                                    freqs_k_neighbours = freqs_sorted_neighbours
                                        .select(Axis(0), &(0..k).collect::<Vec<usize>>());
                                    dist_k_neighbours = dist
                                        .column(i)
                                        .select(Axis(0), &idx_pools)
                                        .select(Axis(0), &(0..k).collect::<Vec<usize>>());
                                }
                                // Keep only non-missing frequencies and distances among the k-neighbours
                                let (freqs_k_neighbours, dist_k_neighbours): (Vec<f64>, Vec<f64>) =
                                    freqs_k_neighbours
                                        .iter()
                                        .zip(dist_k_neighbours.iter())
                                        .filter(|&(&x, &y)| (!x.is_nan()) & (!y.is_nan()))
                                        .unzip();
                                if freqs_k_neighbours.len() == 0 {
                                    continue;
                                }
                                println!("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
                                println!("freqs_k_neighbours:\n{:?}", freqs_k_neighbours);
                                println!("dist_k_neighbours:\n{:?}", dist_k_neighbours);
                                // Impute the missing data using the weighted allele frequencies from the k-nearest neighbours
                                let dist_sum =
                                    dist_k_neighbours.iter().fold(0.0, |sum, &x| sum + x);
                                // println!("dist_sum:\n{:?}", dist_sum);
                                let weights = dist_k_neighbours
                                    .iter()
                                    .map(|&x| 1.0 - (x / dist_sum) + f64::EPSILON)
                                    .collect::<Vec<f64>>(); // add epsilon that we keep the weight of non-zero distance neighbours when zero distance neighbours exist
                                                            // println!("weights:\n{:?}", weights);
                                let weights_sum = weights.iter().fold(0.0, |sum, &x| sum + x);
                                let weights = Array1::from_vec(
                                    weights
                                        .iter()
                                        .map(|&x| x / weights_sum)
                                        .collect::<Vec<f64>>(),
                                ); // re-weigh after adding epsilon so that the weights still sum up to one
                                let values = Array1::from_vec(freqs_k_neighbours);
                                println!("weights:\n{:?}", weights);
                                println!("(values * weights).sum():\n{:?}", (&values * &weights).sum());
                                window_freqs[(i, j)] = (values * weights).sum();
                                // Need to correct for when the imputed allele frequencies do not add up to one!
                                if j > 0 {
                                    // Find the index of the index of each locus' end position
                                    let idx_idx_locus_fin = self
                                        .start_index_of_each_locus
                                        .iter()
                                        .enumerate()
                                        .filter(|&(_i, &x)| {
                                            (x > 0) & ((x as usize - 1) == j + idx_ini)
                                        })
                                        .map(|(i, _x)| i)
                                        .collect::<Vec<usize>>();
                                    let idx_idx_locus_fin = if idx_idx_locus_fin.len() > 0 {
                                        idx_idx_locus_fin[0]
                                    } else {
                                        continue;
                                    };
                                    let idx_idx_locus_ini = idx_idx_locus_fin - 1;
                                    let idx_locus_ini =
                                        self.start_index_of_each_locus[idx_idx_locus_ini] as usize
                                            - idx_ini;
                                    let idx_locus_fin =
                                        self.start_index_of_each_locus[idx_idx_locus_fin] as usize
                                            - idx_ini;
                                    let locus_freqs_imputed = window_freqs
                                        .slice(s![i, idx_locus_ini..idx_locus_fin])
                                        .to_owned();
                                    let locus_sum = locus_freqs_imputed.sum();
                                    let locus_freqs_imputed =
                                        locus_freqs_imputed.map(|&x| x / locus_sum);
                                    for j_ in idx_locus_ini..idx_locus_fin {
                                        window_freqs[(i, j_)] =
                                            locus_freqs_imputed[j_ - idx_locus_ini];
                                    }
                                }
                            }
                        }
                    }
                }
                // println!("*window_freqs:\n{:?}", *window_freqs);
                println!("*window_freqs:\n{:?}", window_freqs);
            // });
            }
        // println!("vec_windows_freqs[0]:\n{:?}", vec_windows_freqs[0]);
        // Write-out the imputed data
        for a in 0..w {
            // Use the indexes of each locus
            let idx_ini = idx_window_head[a];
            let idx_fin = idx_window_tail[a];
            let p = idx_fin - idx_ini;
            for i in 0..n {
                for j in 0..p {
                    self.intercept_and_allele_frequencies[(i, idx_ini + j)] =
                        vec_windows_freqs[a][(i, j)];
                }
            }
        }
        Ok(self)
    }
}

// Impute using adaptive linkage disequillibrium (estimated using correlations within a window) k-nearest neighbour weighted allele frequencies imputation
pub fn impute_aLDkNN(
    file_sync_phen: &FileSyncPhen,
    filter_stats: &FilterStats,
    min_depth_set_to_missing: &f64,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    min_correlation: &f64,
    k_neighbours: &u64,
    n_threads: &usize,
    out: &String,
) -> io::Result<String> {
    // All non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
    let keep_p_minus_1 = false;
    let start = std::time::SystemTime::now();
    let mut genotypes_and_phenotypes = file_sync_phen
        .into_genotypes_and_phenotypes(filter_stats, keep_p_minus_1, n_threads)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!("Parsing the sync into allele frequncies: {} seconds", duration.as_secs());
    println!("genotypes_and_phenotypes.intercept_and_allele_frequencies:\n{:?}", genotypes_and_phenotypes.intercept_and_allele_frequencies);
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(min_depth_set_to_missing)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!("Set missing loci below the minimum depth: {} seconds", duration.as_secs());
    println!("genotypes_and_phenotypes.intercept_and_allele_frequencies:\n{:?}", genotypes_and_phenotypes.intercept_and_allele_frequencies);
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .adaptive_ld_knn_imputation(
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
            min_correlation,
            k_neighbours,
        )
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!("Imputation: {} seconds", duration.as_secs());
    let out = genotypes_and_phenotypes
        .write_csv(filter_stats, keep_p_minus_1, out, n_threads)
        .unwrap();
    Ok(out)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_impute_aLDkNN() {
        let file_sync = FileSync {
            filename: "./tests/test.sync".to_owned(),
            test: "load".to_owned(),
        };
        let file_phen = FilePhen {
            filename: "./tests/test.csv".to_owned(),
            delim: ",".to_owned(),
            names_column_id: 0,
            sizes_column_id: 1,
            trait_values_column_ids: vec![2, 3],
            format: "default".to_owned(),
        };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        let filter_stats = FilterStats {
            remove_ns: true,
            min_quality: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            min_missingness_rate: 0.0,
            pool_sizes: vec![20., 20., 20., 20., 20.],
        };
        let n_threads = 2;
        let mut frequencies_and_phenotypes = file_sync_phen
            .into_genotypes_and_phenotypes(&filter_stats, false, &n_threads)
            .unwrap();
        println!(
            "Before simulating missing data:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let min_depth_set_to_missing = 5.0;
        let _ = frequencies_and_phenotypes
            .set_missing_by_depth(&min_depth_set_to_missing)
            .unwrap();
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 2)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 5)] = f64::NAN;
        // frequencies_and_phenotypes.intercept_and_allele_frequencies[(3, 6)] = f64::NAN;
        println!(
            "Before imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );
        let window_size_bp = 1e6 as u64;
        let window_slide_size_bp = window_size_bp;
        let min_loci_per_window = 2;
        let min_correlation = 0.8;
        let k_neighbours = 3;
        let _ = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &window_size_bp,
                &window_slide_size_bp,
                &min_loci_per_window,
                &min_correlation,
                &k_neighbours,
            )
            .unwrap();
        println!(
            "After imputation:\n{:?}",
            frequencies_and_phenotypes.intercept_and_allele_frequencies
        );

        let outname = impute_aLDkNN(
            &file_sync_phen,
            &filter_stats,
            &min_depth_set_to_missing,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
            &min_correlation,
            &k_neighbours,
            &n_threads,
            &"test-impute_aLDkNN.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_aLDkNN.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 1)],
            0.19663299601190515
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 2)],
            0.8033670039880948
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 13923)],
            0.02733470369319045
        );
        assert_eq!(
            frequencies_and_phenotypes.intercept_and_allele_frequencies[(0, 13924)],
            0.9726652963068096
        );
        // assert_eq!(0, 1);
    }
}


// SYNC=/home/jeff/poolgen/tests/test-10k_loci.sync
// PHEN=/home/jeff/poolgen/tests/test-10k_loci.csv
// NCORES=7
// OUT=/home/jeff/poolgen/tests/test-impute.10k_loci.csv
// time cargo run -- impute \
//     -f ${SYNC} \
//     -p ${PHEN} \
//     --phen-delim , \
//     --phen-name-col 0 \
//     --phen-pool-size-col 1 \
//     --phen-value-col 2 \
//     --min-allele-frequency 0.0001 \
//     --min-coverage 0 \
//     --min-quality 0.01 \
//     --min-missingness-rate 0.2 \
//     --min-depth-set-to-missing 5 \
//     --window-size-bp 1000000 \
//     --window-slide-size-bp 1000000 \
//     --min-loci-per-window 1 \
//     --min-correlation 0.5 \
//     --k-neighbours 5 \
//     --n-threads ${NCORES} \
//     -o ${OUT}