//! Structs and Traits

use ndarray::prelude::*;
use std::io;

/// The entry point genotype file struct, i.e. describing the genotype data in pileup format
/// - `filename` - filename of the pileup file (`*.pileup`)
/// - `pool_names` - names of the pools from the phenotype file
#[derive(Debug, Clone)]
pub struct FilePileup {
    pub filename: String,
    pub pool_names: Vec<String>,
}

/// The alternative entry point genotype file struct, i.e. describing the genotype data in vcf format
/// Note: Make sure that the annotate the vcf file with allele depths, e.g. `bcftools mpileup -a AD ...`
/// - `filename` - filename of the vcf file (`*.vcf` or `*.vcf.gz`)
#[derive(Debug, Clone)]
pub struct FileVcf {
    pub filename: String,
}

/// The main genotype file struct used for most of the analyses
/// - `filename` - filename of the synchronised pileup file (`*.sync`)
/// - `test` - name of statistical test, i.e. "sync2csv", "fisher_exact_test", "chisq_test", "fst", "heterozygosity, "pearson_corr", "ols_iter", "ols_iter_with_kinship", "mle_iter", "mle_iter_with_kinship", "gwalpha", "genomic_prediction_cross_validation""
#[derive(Debug, Clone)]
pub struct FileSync {
    pub filename: String,
    pub test: String,
}

/// Filename of the phenotype file which can be a simple delimited file (e.g. csv and tsv) or a specialised GWAlpha phenotype infomation file in a python file.
/// - `filename` - filename of the phenotype file (e.g. `*.csv`, `*.txt`, `*.tsv`, or `*.py`)
/// - `delim` - string delimiter of the phenotype file (e.g. `","` or `"\t"`)
/// - `names_column_id` - index of the column containing the names of the pools or populations
/// - `sizes_column_id` - index of the column contating the sizes of each pool or population
/// - `trait_values_column_ids` - vector of indexes corresponding to the column containing the trait values to be included in the analyses (Note that multi-trait analyese may not be available to all analyses types)
/// - `format` - string defining the format of the phenotype file as `default` for simple delimited file or `gwalpha_fmt` for GWAlpha-required format (for back-compatibility with github.com/aflevel/GWAlpha/GWAlpha.py)
#[derive(Debug, Clone)]
pub struct FilePhen {
    pub filename: String,
    pub delim: String,
    pub names_column_id: usize,
    pub sizes_column_id: usize,
    pub trait_values_column_ids: Vec<usize>,
    pub format: String,
}

/// Phenotype data including the names of the pools, the size of each pool, and the trait values
#[derive(Debug, Clone, PartialEq)]
pub struct Phen {
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: Array2<f64>,
}

/// Filename of the synchronised pileup file and its corresponding phenotype data
#[derive(Debug, Clone, PartialEq)]
pub struct FileSyncPhen {
    pub filename_sync: String,
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: Array2<f64>,
    pub test: String,
}

/// Locus and allele filtering settings
#[derive(Debug, Clone)]
pub struct FilterStats {
    pub remove_ns: bool,
    pub remove_monoallelic: bool,
    pub keep_lowercase_reference: bool,
    pub max_base_error_rate: f64,
    pub min_coverage_depth: u64,
    pub min_coverage_breadth: f64,
    pub min_allele_frequency: f64,
    pub max_missingness_rate: f64,
    pub pool_sizes: Vec<f64>,
}

/// A line of a pileup file corresponding to a single locus across all the pools
#[derive(Debug, Clone, PartialEq)]
pub struct PileupLine {
    pub chromosome: String,           // chromosome or scaffold name
    pub position: u64,                // position in number of bases
    pub reference_allele: char,       // reference allele
    pub coverages: Vec<u64>,          // number of times the locus was covered
    pub read_codes: Vec<Vec<u8>>, // utf8 read codes corresponding to 'A', 'T', 'C', or 'G' (252 other alleles can be accommodated)
    pub read_qualities: Vec<Vec<u8>>, // utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

/// A line of a vcf file corresponding to a single locus across all the pools
/// We are interested in extracting allele counts.
/// We are not interested in the genotype calls and their corresponding likelihoods.
#[derive(Debug, Clone, PartialEq)]
pub struct VcfLine {
    pub chromosome: String,             // chromosome or scaffold name
    pub position: u64,                  // position in number of bases
    pub reference_allele: char,         // reference allele
    pub alternative_alleles: Vec<char>, // vector of alternative alleles
    pub allele_depths: Vec<Vec<u64>>, // across samples average utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

/// Allele counts at a locus across pools
#[derive(Debug, Clone, PartialEq)]
pub struct LocusCounts {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: Array2<u64>, // n pools x p alleles
}

// Struct of allele frequencies to convert reads into syncx
// #[derive(Debug, Clone, PartialEq, PartialOrd)]
#[derive(Debug, Clone, PartialEq)]
pub struct LocusFrequencies {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: Array2<f64>, // n pools x p alleles
}

// Struct for tracking the insertions and deletions prepended with +/i which we want to skip as they refer to the next base position and not the current locus
#[derive(Debug, Clone)]
pub struct IndelMarker {
    pub indel: bool,  // insertion or deletion, i.e. +/- codes
    pub count: usize, // size of the indel, i.e. how many bases
    pub left: usize, // how many indel bases are left to be removed (initialised with the maximum possible value of 4294967295)
}

// Struct of allele counts and phenotypes per pool
#[derive(Debug, Clone)]
pub struct LocusCountsAndPhenotypes {
    pub locus_counts: LocusCounts,
    pub phenotypes: Array2<f64>, // n pools x k traits
    pub pool_names: Vec<String>,
}

// Struct of allele frequencies and phenotypes for genomic prediction
#[derive(Debug, Clone)]
pub struct GenotypesAndPhenotypes {
    pub chromosome: Vec<String>,                       // 1 + p
    pub position: Vec<u64>,                            // 1 + p
    pub allele: Vec<String>,                           // 1 + p
    pub intercept_and_allele_frequencies: Array2<f64>, // n pools x 1 + p alleles across loci
    pub phenotypes: Array2<f64>,                       // n pools x k traits
    pub pool_names: Vec<String>,                       // n
    pub coverages: Array2<f64>,                        // n pools x m loci
}

// Struct for GWAlpha's least squares cost function minimisation
#[derive(Debug, Clone)]
pub struct LeastSquaresBeta {
    pub q_prime: Array1<f64>,
    pub percs_a: Array1<f64>,
    pub percs_b: Array1<f64>,
}

// Struct for GWAlpha's maximum likelihood estimation
#[derive(Debug, Clone)]
pub struct MaximumLikelihoodBeta {
    pub percs_a: Array1<f64>,
    pub percs_a0: Array1<f64>,
    pub percs_b: Array1<f64>,
    pub percs_b0: Array1<f64>,
}

// Struct for gudmc's maximum likelihood estimation for the distribution of Tajima's D which is assumed to be normally distributed centred on zero
#[derive(Debug, Clone)]
pub struct MaximumLikelihoodNormal {
    pub q: Array1<f64>,
}

// Struct for regression objects
#[derive(Debug, Clone)]
pub struct UnivariateOrdinaryLeastSquares {
    pub x: Array2<f64>,
    pub y: Array1<f64>,
    pub b: Array1<f64>,
    pub xt: Array2<f64>,
    pub inv_xtx: Array2<f64>,
    pub inv_xxt: Array2<f64>,
    pub e: Array1<f64>,
    pub ve: f64,
    pub v_b: Array1<f64>,
    pub t: Array1<f64>,
    pub pval: Array1<f64>,
}

#[derive(Debug, Clone)]
pub struct UnivariateMaximumLikelihoodEstimation {
    pub x: Array2<f64>,
    pub y: Array1<f64>,
    pub b: Array1<f64>,
    pub ve: f64,
    pub v_b: Array1<f64>,
    pub t: Array1<f64>,
    pub pval: Array1<f64>,
}

#[derive(Debug, Clone)]
pub struct PredictionPerformance {
    pub n: usize,                                // number of observations
    pub p: usize,                                // number of predictors
    pub k: usize,                                // number cross-validation folds
    pub r: usize, // number replications where each cross-validation results in random groupings
    pub models: Vec<String>, // genomic prediction model used
    pub y_validation_and_predicted: Array4<f64>, // reps x n x models x traits+traits
    pub cor: Array4<f64>, // reps x folds x models x traits
    pub mbe: Array4<f64>, // reps x folds x models x traits
    pub mae: Array4<f64>, // reps x folds x models x traits
    pub mse: Array4<f64>, // reps x folds x models x traits
    pub rmse: Array4<f64>, // reps x folds x models x traits
}

////////////////////////////////////////////////////////////////////////////////
/// # TRAITS
////////////////////////////////////////////////////////////////////////////////

pub trait CheckStruct {
    fn check(&self) -> io::Result<()>;
}

pub trait Count {
    fn count_loci(&self) -> io::Result<(Vec<usize>, Vec<String>, Vec<u64>)>;
}

pub trait Parse<T> {
    fn lparse(&self) -> io::Result<Box<T>>;
}

pub trait Filter {
    fn to_counts(&self) -> io::Result<Box<LocusCounts>>;
    fn to_frequencies(&self) -> io::Result<Box<LocusFrequencies>>;
    fn filter(&mut self, filter_stats: &FilterStats) -> io::Result<&mut Self>;
}

pub trait Sort {
    fn sort_by_allele_freq(&mut self, decreasing: bool) -> io::Result<&mut Self>;
}

pub trait RemoveMissing {
    fn remove_missing(&mut self) -> io::Result<&mut Self>;
}

pub trait ChunkyReadAnalyseWrite<T, F> {
    fn per_chunk(
        &self,
        start: &u64,
        end: &u64,
        outname_ndigits: &usize,
        filter_stats: &FilterStats,
        function: F,
    ) -> io::Result<String>
    where
        F: Fn(&mut T, &FilterStats) -> Option<String>;
    fn read_analyse_write(
        &self,
        filter_stats: &FilterStats,
        out: &String,
        n_threads: &usize,
        function: F,
    ) -> io::Result<String>
    where
        F: Fn(&mut T, &FilterStats) -> Option<String>;
}

pub trait LoadAll {
    fn per_chunk_load(
        &self,
        start: &u64,
        end: &u64,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)>; // Allele frequencies and counts across pools and alleles per locus
    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<(Vec<LocusFrequencies>, Vec<LocusCounts>)>; // Allele frequencies and counts across pools and alleles per locus
    fn into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes>;
}

pub trait SaveCsv {
    fn write_csv(
        &self,
        filter_stats: &FilterStats,
        keep_p_minus_1: bool,
        out: &String,
        n_threads: &usize,
    ) -> io::Result<String>;
}

pub trait Regression {
    fn new() -> Self;
    fn remove_collinearities_in_x(&mut self) -> &mut Self;
    fn estimate_effects(&mut self) -> io::Result<&mut Self>;
    fn estimate_variances(&mut self) -> io::Result<&mut Self>;
    fn estimate_significance(&mut self) -> io::Result<&mut Self>;
}

pub trait MoorePenrosePseudoInverse {
    fn pinv(&self) -> io::Result<Array2<f64>>;
}

pub trait CrossValidation<F> {
    fn k_split(&self, k: usize) -> io::Result<(Vec<usize>, usize, usize)>;
    fn performance(
        &self,
        y_true: &Array2<f64>,
        y_hat: &Array2<f64>,
    ) -> io::Result<Vec<Array1<f64>>>;
    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        functions: Vec<F>,
    ) -> io::Result<PredictionPerformance>
    where
        F: Fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>;
    fn tabulate_predict_and_output(
        &self,
        prediction_performance: &PredictionPerformance,
        functions: Vec<F>,
        fname_input: &String,
        fname_output: &String,
    ) -> io::Result<(String, String, Vec<String>)>
    where
        F: Fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>;
}
