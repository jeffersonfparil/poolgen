use crate::base::*;
use crate::gwas::pearsons_correlation;
use ndarray::{prelude::*, Zip};
use std::io;

// // For benchmarking
// use std::time::Instant;

fn calculate_euclidean_distances(
    window_freqs: ArrayView2<f64>,
    corr: ArrayView1<f64>,
    n_loci_to_estimate_distance: &u64,
) -> io::Result<(Array2<f64>, bool)> {
    // println!("@@@@@@@@@@@@@@@@@@@calculate_euclidean_distances@@@@@@@@@@@@@@@@@@@");
    // let start = Instant::now();
    let (n, p) = window_freqs.dim();
    let n_loci_to_estimate_distance = if p < *n_loci_to_estimate_distance as usize {
        p
    } else {
        *n_loci_to_estimate_distance as usize
    };
    // Sort the correlations and select n_loci_to_estimate_distance loci to compute pairwise Euclidean distances between pools from
    let mut idx_linked_alleles = (0..p).collect::<Vec<usize>>();
    idx_linked_alleles.sort_by(|&j0, &j1| {
        let a = match corr[j0].is_nan() {
            true => f64::INFINITY,
            false => corr[j0],
        };
        let b = match corr[j1].is_nan() {
            true => f64::INFINITY,
            false => corr[j1],
        };
        // Sort from large to small values
        b.partial_cmp(&a).unwrap()
    });
    let idx_linked_alleles = idx_linked_alleles[0..n_loci_to_estimate_distance].to_vec();
    // Calculate the Euclidean distances between pools using the most positively correlated alleles (or all the alleles if none passed the minimum correlation threshold or if there is less that 2 loci that are most correlated - these minimum 2 is the allele itself and another allele)
    let mut dist: Array2<f64> = Array2::from_elem((n, n), f64::NAN);
    let window_freqs_linked_alleles = window_freqs.select(Axis(1), &idx_linked_alleles);
    let mut all_missing = true;
    for i0 in 0..n {
        let pool0_tmp = window_freqs_linked_alleles.row(i0);
        for i1 in i0..n {
            // Keep only the loci present in both pools
            let (pool0, pool1): (Vec<f64>, Vec<f64>) = pool0_tmp
                .iter()
                .zip(window_freqs_linked_alleles.row(i1))
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
                        dist[(i1, i0)] = d;
                        all_missing = match all_missing {
                            true => false,
                            false => false,
                        };
                    }
                };
            }
        }
    }
    // let duration = start.elapsed();
    // println!("duration = {:?}", duration);
    Ok((dist, all_missing))
}

fn mean_value_imputation(freqs: ArrayView1<f64>) -> io::Result<f64> {
    Ok(freqs
        .iter()
        .filter(|&x| !x.is_nan())
        .fold(0.0, |sum, &x| sum + 1.0)
        / (freqs.len() as f64))
}

fn find_k_nearest_neighbours(
    k: &mut usize,
    window_freqs: ArrayView1<f64>,
    dist: ArrayView1<f64>,
) -> io::Result<(Vec<f64>, Vec<f64>, Array1<f64>)> {
    // println!("@@@@@@@@@@@@@@@@@@@find_k_nearest_neighbours@@@@@@@@@@@@@@@@@@@");
    // let start = Instant::now();
    let n = window_freqs.len();
    // Find the k-nearest neighbours
    let mut idx_pools: Vec<usize> = (0..n).collect();
    idx_pools.sort_by(|&j0, &j1| {
        let a = match dist[j0].is_nan() {
            true => f64::INFINITY,
            false => dist[j0],
        };
        let b = match dist[j1].is_nan() {
            true => f64::INFINITY,
            false => dist[j1],
        };
        // Sort from small to large values
        a.partial_cmp(&b).unwrap()
    });
    let freqs_sorted_neighbours = window_freqs.select(Axis(0), &idx_pools).to_owned();
    // Extract the allele frequencies and distances of the k-neighbours from sorted lists using idx_pools whose indexes are sorted
    let mut freqs_k_neighbours =
        freqs_sorted_neighbours.select(Axis(0), &(0..*k).collect::<Vec<usize>>());
    while *k < n {
        if freqs_k_neighbours.fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum }) > 0 {
            break;
        }
        freqs_k_neighbours =
            freqs_sorted_neighbours.select(Axis(0), &(0..*k).collect::<Vec<usize>>());
        *k += 1;
    }
    let dist_k_neighbours = dist
        .select(Axis(0), &idx_pools)
        .select(Axis(0), &(0..*k).collect::<Vec<usize>>());
    // Keep only non-missing frequencies and distances among the k-neighbours
    let (freqs_k_neighbours, dist_k_neighbours): (Vec<f64>, Vec<f64>) = freqs_k_neighbours
        .iter()
        .zip(dist_k_neighbours.iter())
        .filter(|&(&x, &y)| (!x.is_nan()) & (!y.is_nan()))
        .unzip();
    // let duration = start.elapsed();
    // println!("duration = {:?}", duration);
    Ok((
        freqs_k_neighbours,
        dist_k_neighbours,
        freqs_sorted_neighbours,
    ))
}

fn correct_allele_frequencies_per_locus<'w>(
    a: usize,
    i: usize,
    j: usize,
    idx_ini: usize,
    window_freqs: &'w mut Array2<f64>,
    idx_window_head: &Vec<usize>,
    idx_window_tail: &Vec<usize>,
    loci_idx: &Vec<usize>,
) -> io::Result<&'w mut Array2<f64>> {
    // println!("@@@@@@@@@@@@@@@@@@@correct_allele_frequencies_per_locus@@@@@@@@@@@@@@@@@@@");
    // let start = Instant::now();
    // Include the start of the next window, i.e. the marker for the end of the last locus in the current window
    let loci_start_indexes_within_the_current_window =
        loci_idx[idx_window_head[a]..(idx_window_tail[a] + 2)].to_vec();
    for j_ in 1..loci_start_indexes_within_the_current_window.len() {
        // Are we at the last allele of the locus?
        if (loci_start_indexes_within_the_current_window[j_] - 1) == (idx_ini + j) {
            // If we are then we find the start of this locus, i.e. its local index
            let j_ini = loci_start_indexes_within_the_current_window[j_ - 1] - idx_ini;
            let freqs_sum = window_freqs
                .slice(s![i, j_ini..(j + 1)])
                .fold(0.0, |sum, &x| if !x.is_nan() { sum + x } else { sum })
                + f64::EPSILON;
            if freqs_sum != 1.0 {
                for j_ in j_ini..(j + 1) {
                    window_freqs[(i, j_)] = window_freqs[(i, j_)] / freqs_sum;
                }
            }
            break;
        }
    }
    // let duration = start.elapsed();
    // println!("duration = {:?}", duration);
    Ok(window_freqs)
}

impl GenotypesAndPhenotypes {
    pub fn adaptive_ld_knn_imputation(
        &mut self,
        window_size_bp: &u64,
        window_slide_size_bp: &u64,
        min_loci_per_window: &u64,
        // min_correlation: &f64,
        n_loci_to_estimate_distance: &u64,
        k_neighbours: &u64,
    ) -> io::Result<&mut Self> {
        self.check().unwrap();
        // We are assuming that all non-zero alleles across pools are kept, i.e. biallelic loci have 2 columns, triallelic have 3, and so on.
        let (n, _p) = self.intercept_and_allele_frequencies.dim();
        // Define sliding windows
        let (loci_idx, loci_chr, loci_pos) = self.count_loci().unwrap();
        let mut loci_chr_no_redundant_tail = loci_chr.to_owned();
        loci_chr_no_redundant_tail.pop();
        let mut loci_pos_no_redundant_tail = loci_pos.to_owned();
        loci_pos_no_redundant_tail.pop();
        let (idx_window_head, idx_window_tail) = define_sliding_windows(
            &loci_chr_no_redundant_tail,
            &loci_pos_no_redundant_tail,
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
        )
        .unwrap();
        let w = idx_window_head.len();
        // Parallel processing per window
        let mut vec_windows_freqs: Vec<Array2<f64>> = vec![Array2::from_elem((1, 1), f64::NAN); w];
        let idx_windows: Vec<usize> = (0..w).collect();
        Zip::from(&mut vec_windows_freqs)
            .and(&idx_windows)
            .par_for_each(
                |window_freqs, &a| {
                    // for a in 0..idx_windows.len() {
                    let idx_ini = loci_idx[idx_window_head[a]];
                    let idx_fin = loci_idx[idx_window_tail[a] + 1]; // add one so that we include the last part of the window!
                    let p = idx_fin - idx_ini;
                    *window_freqs = self
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
                    for j in 0..p {
                        if window_freqs
                            .column(j)
                            .fold(0, |sum, &x| if x.is_nan() { sum + 1 } else { sum })
                            == 0
                        {
                            // No missing data on the locus
                            continue;
                        } else {
                            let (dist, all_missing) = calculate_euclidean_distances(
                                window_freqs.view(),
                                corr.column(j),
                                n_loci_to_estimate_distance,
                            )
                            .unwrap();
                            // Now, let's find the pools needing imputation and impute them using the k-nearest neighbours
                            for i in 0..n {
                                let mut k = *k_neighbours as usize; // define the adaptive k for use later
                                if window_freqs[(i, j)].is_nan() == false {
                                    continue;
                                } else {
                                    // If the data is so sparse such that all pairwise distances are impossible to compute across n_loci_to_estimate_distance loci, then we simply impute using the mean allele frequency at the locus across pools
                                    if all_missing == true {
                                        window_freqs[(i, j)] =
                                            mean_value_imputation(window_freqs.column(j)).unwrap();
                                    } else {
                                        let (
                                            freqs_k_neighbours,
                                            dist_k_neighbours,
                                            freqs_sorted_neighbours,
                                        ) = find_k_nearest_neighbours(
                                            &mut k,
                                            window_freqs.column(j),
                                            dist.column(i),
                                        )
                                        .unwrap();
                                        if (freqs_k_neighbours.len() == 0) | (all_missing == true) {
                                            // If the pool freqs and distances do not correspond to non-missing info then we simply use the mean across all non-missing pools
                                            window_freqs[(i, j)] = mean_value_imputation(
                                                freqs_sorted_neighbours.view(),
                                            )
                                            .unwrap();
                                        } else {
                                            // Impute the missing data using the weighted allele frequencies from the k-nearest neighbours
                                            let dist_sum = dist_k_neighbours
                                                .iter()
                                                .fold(0.0, |sum, &x| sum + x)
                                                + f64::EPSILON; // Add a very small value so that we don;t get undefined values when all distances are zero
                                                                // println!("dist_sum:\n{:?}", dist_sum);
                                            let weights = dist_k_neighbours
                                                .iter()
                                                .map(|&x| 1.0 - (x / dist_sum) + f64::EPSILON)
                                                .collect::<Vec<f64>>(); // add epsilon that we keep the weight of non-zero distance neighbours when zero distance neighbours exist
                                                                        // println!("weights:\n{:?}", weights);
                                            let weights_sum =
                                                weights.iter().fold(0.0, |sum, &x| sum + x);
                                            let weights = Array1::from_vec(
                                                weights
                                                    .iter()
                                                    .map(|&x| x / weights_sum)
                                                    .collect::<Vec<f64>>(),
                                            ); // re-weigh after adding epsilon so that the weights still sum up to one
                                            let values = Array1::from_vec(freqs_k_neighbours);
                                            // if (&values * &weights).sum().is_nan() {
                                            //     println!("dist_sum:\n{:?}", dist_sum);
                                            //     println!("weights:\n{:?}", weights);
                                            //     println!("values:\n{:?}", values);
                                            // }
                                            window_freqs[(i, j)] = (&values * &weights).sum();
                                        }
                                    } // If distance matrix is missing due to extreme sparsity
                                }
                                // Need to correct for when the imputed allele frequencies do not add up to one!
                                if j > 0 {
                                    correct_allele_frequencies_per_locus(
                                        a,
                                        i,
                                        j,
                                        idx_ini,
                                        window_freqs,
                                        &idx_window_head,
                                        &idx_window_tail,
                                        &loci_idx,
                                    )
                                    .unwrap();
                                }
                            } // Impute across pools with missing data
                        } // Impute if we have missing data
                    } // Iterate across alleles across loci within the window
                      // println!("window_freqs:\n{:?}", &window_freqs);
                }, // Parallel processing across windows
            );
        // }
        // println!("vec_windows_freqs[0]:\n{:?}", vec_windows_freqs[0]);
        // Write-out the imputed data
        // println!("@@@@@@@@@@@@@@@@@@@Writing imputed data@@@@@@@@@@@@@@@@@@@");
        // let start = Instant::now();
        for a in 0..w {
            // Use the indexes of each locus
            let idx_ini = loci_idx[idx_window_head[a]];
            let idx_fin = loci_idx[idx_window_tail[a] + 1]; // add one so that we include the last part of the window!
            let p = idx_fin - idx_ini;
            for i in 0..n {
                for j in 0..p {
                    self.intercept_and_allele_frequencies[(i, idx_ini + j)] =
                        vec_windows_freqs[a][(i, j)];
                }
            }
        }
        // let duration = start.elapsed();
        // println!("duration = {:?}", duration);
        // Set missing coverages to infinity to mark imputed data
        // println!("@@@@@@@@@@@@@@@@@@@Set missing coverages to infinity to mark imputed data@@@@@@@@@@@@@@@@@@@");
        // let start = Instant::now();
        for j in 0..self.coverages.ncols() {
            // Mark only the imputed loci, i.e. loci which were not completely missing across all pools
            let n_non_missing = self.coverages.select(Axis(1), &[j]).fold(0, |sum, &x| {
                if x.is_nan() == false {
                    sum + 1
                } else {
                    sum
                }
            });
            if n_non_missing > 0 {
                for i in 0..self.coverages.nrows() {
                    if self.coverages[(i, j)].is_nan() {
                        self.coverages[(i, j)] = f64::INFINITY
                    };
                }
            }
        }
        // let duration = start.elapsed();
        // println!("duration = {:?}", duration);
        Ok(self)
    }
}

// Impute using adaptive linkage disequillibrium (estimated using correlations within a window) k-nearest neighbour weighted allele frequencies imputation
pub fn impute_aLDkNN(
    file_sync_phen: &FileSyncPhen,
    filter_stats: &FilterStats,
    min_depth_set_to_missing: &f64,
    frac_top_missing_pools: &f64,
    frac_top_missing_loci: &f64,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
    n_loci_to_estimate_distance: &u64,
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
    println!(
        "Parsed the sync file into allele frequncies: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .set_missing_by_depth(min_depth_set_to_missing)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Set missing loci below the minimum depth: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_pools(frac_top_missing_pools)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest pools: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(frac_top_missing_loci)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Filtered out sparsest loci: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .adaptive_ld_knn_imputation(
            window_size_bp,
            window_slide_size_bp,
            min_loci_per_window,
            n_loci_to_estimate_distance,
            k_neighbours,
        )
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Adaptive LD-kNN imputation: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Remove 100% of the loci with missing data
    let start = std::time::SystemTime::now();
    genotypes_and_phenotypes
        .filter_out_top_missing_loci(&1.00)
        .unwrap();
    let end = std::time::SystemTime::now();
    let duration = end.duration_since(start).unwrap();
    println!(
        "Missing data removed, i.e. loci which cannot be imputed because of extreme sparsity: {} pools x {} loci | Missingness: {}% | Duration: {} seconds",
        genotypes_and_phenotypes.coverages.nrows(),
        genotypes_and_phenotypes.coverages.ncols(),
        genotypes_and_phenotypes.missing_rate().unwrap(),
        duration.as_secs()
    );
    // Output
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
            max_bases_error_rate: 0.005,
            min_coverage: 1,
            min_allele_frequency: 0.005,
            max_missingness_rate: 0.0,
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
        let frac_top_missing_pools = 0.0;
        let frac_top_missing_loci = 0.2;
        let window_size_bp = 1e6 as u64;
        let window_slide_size_bp = window_size_bp;
        let min_loci_per_window = 1;
        let n_loci_to_estimate_distance = 10;
        let k_neighbours = 3;
        let _ = frequencies_and_phenotypes
            .adaptive_ld_knn_imputation(
                &window_size_bp,
                &window_slide_size_bp,
                &min_loci_per_window,
                &n_loci_to_estimate_distance,
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
            &frac_top_missing_pools,
            &frac_top_missing_loci,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
            &n_loci_to_estimate_distance,
            &k_neighbours,
            &n_threads,
            &"test-impute_aLDkNN.csv".to_owned(),
        )
        .unwrap();
        assert_eq!(outname, "test-impute_aLDkNN.csv".to_owned()); // Do better!!! Load data - thus working on improving load_table()

        println!("frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42])={:?}", frequencies_and_phenotypes.intercept_and_allele_frequencies.slice(s![0..5, 39..42]));
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 1..3])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 39..42])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 119..121])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
        );
        assert_eq!(
            frequencies_and_phenotypes
                .intercept_and_allele_frequencies
                .slice(s![0..5, 400..402])
                .sum_axis(Axis(1))
                .map(|x| sensible_round(*x, 2)),
            Array1::from_elem(5, 1.0)
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
//     --max-missingness-rate 0.2 \
//     --min-depth-set-to-missing 5 \
//     --window-size-bp 1000000 \
//     --window-slide-size-bp 1000000 \
//     --min-loci-per-window 1 \
//     --n-loci-to-estimate-distance 10 \
//     --k-neighbours 5 \
//     --n-threads ${NCORES} \
//     -o ${OUT}
