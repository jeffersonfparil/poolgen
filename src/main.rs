//! **poolgen**: quantitative and population genetics on pool sequencing (Pool-seq) data

use clap::Parser;
use ndarray::prelude::*;
#[allow(warnings)]
use std::io;
mod base;
mod gp;
mod gwas;
mod tables;
use base::{ChunkyReadAnalyseWrite, CrossValidation, LoadAll, Parse};
use gp::{
    ols, penalise_glmnet, penalise_lasso_like, penalise_lasso_like_with_iterative_proxy_norms,
    penalise_ridge_like, penalise_ridge_like_with_iterative_proxy_norms,
};
use gwas::{mle_with_covariate, ols_with_covariate};
use tables::{fst, pi, tajima_d, watterson_estimator};

#[derive(Parser, Debug)]
#[clap(
    author = "Jeff Paril",
    version = "0.1.0",
    about = "Quantitative and population genetics analyses using pool sequencing data.",
    long_about = "Quantitative and population genetics analyses using pool sequencing data: trying to continue the legacy of the now unmaintained popoolation2 package with the memory safety of Rust."
)]
struct Args {
    /// Analysis to perform (i.e. "pileup2sync", "sync2csv", "fisher_exact_test", "chisq_test", "pearson_corr", "ols_iter", "ols_iter_with_kinship", "mle_iter", "mle_iter_with_kinship", "gwalpha", "ridge_iter", "genomic_prediction_cross_validation", "fst", "heterozygosity", "watterson_estimator", "tajima_d")
    analysis: String,
    /// Filename of the input pileup or synchronised pileup file (i.e. *.pileup, *.sync, *.syncf, or *.syncx)
    #[clap(short, long)]
    fname: String,
    /// Output filename
    #[clap(short, long, default_value = "")]
    output: String,
    /// Minimum base quality
    #[clap(long, default_value_t = 0.01)]
    min_quality: f64,
    /// Minimum depth of coverage
    #[clap(long, default_value_t = 1)]
    min_coverage: u64,
    /// Minimum allele frequency for associating the genotypes with the phenotype/s
    #[clap(long, default_value_t = 0.001)]
    min_allele_frequency: f64,
    /// Keep ambiguous reads during SNP filtering, i.e. keep them coded as Ns
    #[clap(long, action)]
    keep_ns: bool,
    /// Input phenotype file: csv or tsv or any delimited file
    #[clap(short, long)]
    phen_fname: String,
    /// Delimiter of the input phenotype file: comma, tab, etc...
    #[clap(long, default_value = ",")]
    phen_delim: String,
    /// Column index containing the names or IDs of the indivudals in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value_t = 0)]
    phen_name_col: usize,
    /// Column index containing the the sizes of each pool or population: 0, 1, 2, ...
    #[clap(long, default_value_t = 1)]
    phen_pool_size_col: usize,
    /// Column indexes containing the phenotype values in the input phenotype file, e.g. 1 or 1,2,3 or 1,2,3,4 etc ...
    #[clap(
        long,
        use_value_delimiter = true,
        value_delimiter = ',',
        default_value = "2"
    )]
    phen_value_col: Vec<String>,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t = 1)]
    n_threads: usize,
    ////////////////////////////////////////////////////
    ////// Additional parameters
    //////////////////////////////////////////////////////
    /// GWAS iterative OLS using some of the kinship matrix's PCs as covariate
    #[clap(short, long, default_value_t = 0.75)]
    xxt_eigen_variance_explained: f64,
    /// GWAlpha inference method to use: "LS" for least squares or "ML" for maximum likelihood estimation
    #[clap(long, default_value = "ML")]
    gwalpha_method: String,
    /// Sync to csv file conversion to include all alleles or just p-1 excluding the minimum allele
    #[clap(long, action)]
    keep_p_minus_1: bool,
    /// Genomic prediction cross-validation: number of k-fold validations, i.e. number of time the data will be partitioned for training and testing each model
    #[clap(long, default_value_t = 10)]
    k_folds: usize,
    /// Genomic prediction cross-validation: number of replicates of k-fold cross-validation
    #[clap(long, default_value_t = 3)]
    n_reps: usize,
    /// Estimation of population genetics parameters per window, i.e. fst, pi, Watterson's theta, and Tajima's D per population per window: window size in terms of number of bases
    #[clap(long, default_value_t = 100)]
    window_size_bp: u64,
    /// Number of bases to slide the window (a good start will be half the window size)
    #[clap(long, default_value_t = 50)]
    window_slide_size_bp: u64,
    /// Estimation of population genetics parameters per window, i.e. fst, pi, Watterson's theta, and Tajima's D per population per window: minimum number of loci per window
    #[clap(long, default_value_t = 10)]
    min_loci_per_window: u64,
}

/// # poolgen: quantitative and population genetics on pool sequencing (Pool-seq) data
/// - *pileup2sync* - convert a pileup file into a synchronised pileup file with a header
/// - *sync2csv* - convert a synchronised pileup file (sync format) into a smple comma-separated file (csv format)
/// - *fisher_exact_test* - perform Fisher's exact test per locus (sync format), i.e. using allele counts matrix of n-pools x p-alleles (Note: phenotype data is only required for the pool sizes and/or pool names)
/// - *chisq_test* - perform $\Chi^2$ test per locus (sync format), i.e. using allele counts matrix of n-pools x p-alleles (Note: likewise phenotype data is only required for the pool sizes and/or pool names)
/// - *fst* - find the pairwise differentiation between populations using multiple/genome-wide allele frequencies (sync format) with unbiased Fst (fixation index) estimates from multiallelic loci (similar to [Gautier et al, 2019](https://doi.org/10.1111/1755-0998.13557))
/// - *heterozygosity* - find the heterozygosity or nucleotide diversity ($\pi$) of each population using multiple/genome-wide allele frequencies (sync format) with unbiased Fst (fixation index) estimates from multiallelic loci (similar to [Korunes & Samuk, 2021](https://doi.org/10.1111/1755-0998.13326))
/// - *pearson_corr* - find the correlation between phenotypes (csv format) and allele frequencies (sync format) per locus
/// - *ols_iter* - find the strength of association between phenotypes (csv format) and allele frequencies (sync format) per locus using ordinary least squares (OLS) regression
/// - *ols_iter_with_kinship* - similar to *ols_iter* but controls for the effects of kinship
/// - *mle_iter* - similar to *ols_iter* but uses maximum likelihood estimation
/// - *mle_iter_with_kinship* - similar to *ols_iter_with_kinship* but uses maximum likelihood estimation
/// - *gwalpha* - parametric allele effect estimation using Pool-seq data in cases where the number of pools is small, e.g. 5 pools genotyped across thousands to hundred of thousands of loci. See [Fournier-Level et al, 2017](https://doi.org/10.1093/bioinformatics/btw805) for details.
/// - *genomic_prediction_cross_validation* - perform *k*-fold cross-validation with *r* replicates using genomic prediction models (i.e. OLS and various penalised regression models)
///  
/// Please refer to the documentation of each module for more details.
///
/// ## Examples
/// ```shell
/// cargo run -- pileup2sync -f ./tests/test.pileup -p ./tests/test.csv
/// cargo run -- fisher_exact_test -f ./tests/test.sync -p ./tests/test.csv --n-threads 32 --min-coverage 10 --min-allele-frequency 0.01
/// cargo run -- chisq_test -f ./tests/test.sync -p ./tests/test.csv --n-threads 32 --min-coverage 10 --min-allele-frequency 0.01
/// cargo run -- pearson_corr -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage 10 --min-allele-frequency 0.01
/// cargo run -- fst -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// cargo run -- heterozygosity -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// cargo run -- ols_iter -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage 10 --min-allele-frequency 0.01
/// cargo run -- mle_iter -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --min-coverage 10 --min-allele-frequency 0.01
/// cargo run -- gwalpha  -f ./tests/test.sync -p ./tests/test.py --n-threads 32 --gwalpha-method ML
/// cargo run -- sync2csv -f ./tests/test.sync -p ./tests/test.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32 --keep-p-minus-1
/// # cargo run -- genomic_prediction_cross_validation -f ./tests/test_MORE_POOLS.sync -p ./tests/test_MORE_POOLS.csv --phen-delim , --phen-name-col 0 --phen-value-col 2,3  --n-threads 32
/// ```
fn main() {
    let args = Args::parse();
    let mut output: String = String::from("");
    // Prepare the mandatory inputs
    let mut phen_format = "default".to_string();
    if args.analysis == String::from("gwalpha") {
        phen_format = "gwalpha_fmt".to_string()
    }
    let phen_col = args
        .phen_value_col
        .into_iter()
        .map(|x| {
            x.parse::<usize>()
                .expect("Invalid integer input for the phenotype column/s (--phen-value-col).")
        })
        .collect::<Vec<usize>>();
    let file_phen = base::FilePhen {
        filename: args.phen_fname.clone(),
        delim: args.phen_delim.clone(),
        names_column_id: args.phen_name_col,
        sizes_column_id: args.phen_pool_size_col,
        trait_values_column_ids: phen_col.clone(),
        format: phen_format,
    };
    let phen = file_phen.lparse().unwrap();
    let filter_stats = base::FilterStats {
        remove_ns: !args.keep_ns,
        min_quality: args.min_quality,
        min_coverage: args.min_coverage,
        min_allele_frequency: args.min_allele_frequency,
        pool_sizes: phen.pool_sizes.clone(),
    };
    if args.analysis == String::from("pileup2sync") {
        // PILEUP INPUT
        let file_pileup = base::FilePileup {
            filename: args.fname,
            pool_names: phen.pool_names,
        };
        output = file_pileup
            .read_analyse_write(
                &filter_stats,
                &args.output,
                &args.n_threads,
                base::pileup_to_sync,
            )
            .unwrap();
    } else {
        // SYNC INPUT
        let file_sync = base::FileSync {
            filename: args.fname.clone(),
            test: args.analysis.clone(),
        };
        if args.analysis == String::from("fisher_exact_test") {
            output = file_sync
                .read_analyse_write(&filter_stats, &args.output, &args.n_threads, tables::fisher)
                .unwrap();
        } else if args.analysis == String::from("chisq_test") {
            output = file_sync
                .read_analyse_write(&filter_stats, &args.output, &args.n_threads, tables::chisq)
                .unwrap();
        } else if args.analysis == String::from("pearson_corr") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &args.output,
                    &args.n_threads,
                    gwas::correlation,
                )
                .unwrap();
        } else if args.analysis == String::from("ols_iter") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &args.output,
                    &args.n_threads,
                    gwas::ols_iterate,
                )
                .unwrap();
        } else if args.analysis == String::from("ols_iter_with_kinship") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let mut genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, args.keep_p_minus_1, &args.n_threads)
                .unwrap();
            output = ols_with_covariate(
                &mut genotypes_and_phenotypes,
                args.xxt_eigen_variance_explained,
                &args.fname.clone(),
                &args.output,
            )
            .unwrap();
        } else if args.analysis == String::from("mle_iter") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = file_sync_phen
                .read_analyse_write(
                    &filter_stats,
                    &args.output,
                    &args.n_threads,
                    gwas::mle_iterate,
                )
                .unwrap();
        } else if args.analysis == String::from("mle_iter_with_kinship") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, args.keep_p_minus_1, &args.n_threads)
                .unwrap();
            output = mle_with_covariate(
                &genotypes_and_phenotypes,
                args.xxt_eigen_variance_explained,
                &args.fname.clone(),
                &args.output,
            )
            .unwrap();
        } else if args.analysis == String::from("gwalpha") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            if args.gwalpha_method == "LS".to_owned() {
                output = file_sync_phen
                    .read_analyse_write(
                        &filter_stats,
                        &args.output,
                        &args.n_threads,
                        gwas::gwalpha_ls,
                    )
                    .unwrap()
            } else {
                // Defaut is ML, i.e. maximum likelihood estimation
                output = file_sync_phen
                    .read_analyse_write(
                        &filter_stats,
                        &args.output,
                        &args.n_threads,
                        gwas::gwalpha_ml,
                    )
                    .unwrap()
            }
        } else if args.analysis == String::from("sync2csv") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            output = file_sync_phen
                .write_csv(
                    &filter_stats,
                    args.keep_p_minus_1,
                    &args.output,
                    &args.n_threads,
                )
                .unwrap();
        } else if args.analysis == String::from("genomic_prediction_cross_validation") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, args.keep_p_minus_1, &args.n_threads)
                .unwrap();
            let functions: Vec<
                fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>,
            > = vec![
                ols,
                penalise_glmnet,
                penalise_lasso_like,
                penalise_ridge_like,
                penalise_lasso_like_with_iterative_proxy_norms,
                penalise_ridge_like_with_iterative_proxy_norms,
            ];
            let prediction_performances = genotypes_and_phenotypes
                .cross_validate(args.k_folds, args.n_reps, functions.clone())
                .unwrap();
            let (tabulated, _pred_v_expe, predictor_files) = genotypes_and_phenotypes
                .tabulate_predict_and_output(
                    &prediction_performances,
                    functions,
                    &args.fname,
                    &args.output,
                )
                .unwrap();
            output = tabulated;
            let message = "Predictors for each model are here:\n-".to_owned()
                + &predictor_files.join("\n-")[..];
            println!("{:?}", message);
        } else if args.analysis == String::from("fst") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, args.keep_p_minus_1, &args.n_threads)
                .unwrap();
            let (genome_wide, per_window) = fst(
                &genotypes_and_phenotypes,
                &args.window_size_bp,
                &args.window_slide_size_bp,
                &args.min_loci_per_window,
                &args.fname,
                &args.output,
            )
            .unwrap();
            output = genome_wide + " and " + &per_window[..];
        } else if args.analysis == String::from("heterozygosity") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &args.n_threads)
                .unwrap(); // we need all alleles in each locus
            output = pi(
                &genotypes_and_phenotypes,
                &args.window_size_bp,
                &args.window_slide_size_bp,
                &args.min_loci_per_window,
                &args.fname,
                &args.output,
            )
            .unwrap();
        } else if args.analysis == String::from("watterson_estimator") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &args.n_threads)
                .unwrap(); // we need all alleles in each locus
            output = watterson_estimator(
                &genotypes_and_phenotypes,
                &file_sync_phen.pool_sizes,
                &args.window_size_bp,
                &args.window_slide_size_bp,
                &args.min_loci_per_window,
                &args.fname,
                &args.output,
            )
            .unwrap();
        } else if args.analysis == String::from("tajima_d") {
            let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
            let genotypes_and_phenotypes = file_sync_phen
                .into_genotypes_and_phenotypes(&filter_stats, false, &args.n_threads)
                .unwrap(); // we need all alleles in each locus
            output = tajima_d(
                &genotypes_and_phenotypes,
                &file_sync_phen.pool_sizes,
                &args.window_size_bp,
                &args.window_slide_size_bp,
                &args.min_loci_per_window,
                &args.fname,
                &args.output,
            )
            .unwrap();
        } else if args.analysis == String::from("test") {
            let output = 0;
            println!("TEST={:?}", output);
        }
    }
    println!("{}", output);
}
