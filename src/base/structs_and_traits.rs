use nalgebra::{DMatrix, DVector};
use std::io;

// STRUCTS
#[derive(Debug, Clone)]
pub struct FilePileup {
    pub filename: String,
    pub pool_names: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct FileSync {
    pub filename: String,
    pub test: String,
}

#[derive(Debug, Clone)]
pub struct FilePhen {
    pub filename: String,
    pub delim: String,
    pub names_column_id: usize,
    pub sizes_column_id: usize,
    pub trait_values_column_ids: Vec<usize>,
    pub format: String,
}

#[derive(Debug, Clone)]
pub struct Phen {
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: DMatrix<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FileSyncPhen {
    pub filename_sync: String,
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: DMatrix<f64>,
    pub test: String,
}

#[derive(Debug, Clone)]
pub struct FilterStats {
    pub remove_ns: bool,
    pub min_quality: f64,
    pub min_coverage: u64,
    pub min_allele_frequency: f64,
    pub pool_sizes: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PileupLine {
    pub chromosome: String,           // chromosome or scaffold name
    pub position: u64,                // position in number of bases
    pub reference_allele: char,       // reference allele
    pub coverages: Vec<u64>,          // number of times the locus was covered
    pub read_codes: Vec<Vec<u8>>, // utf8 read codes corresponding to 'A', 'T', 'C', or 'G' (252 other alleles can be accommodated)
    pub read_qualities: Vec<Vec<u8>>, // utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

// Struct of allele counts to convert reads into sync
#[derive(Debug, Clone, PartialEq)]
pub struct LocusCounts {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: DMatrix<u64>, // n pools x p alleles
}

// Struct of allele frequencies to convert reads into syncx
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct LocusFrequencies {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: DMatrix<f64>, // n pools x p alleles
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
    pub phenotypes: DMatrix<f64>, // n pools x k traits
    pub pool_names: Vec<String>,
}

// Struct of allele frequencies and phenotypes for genomic prediction
#[derive(Debug, Clone)]
pub struct FrequenciesAndPhenotypes {
    pub chromosome: Vec<String>,
    pub position: Vec<u64>,
    pub intercept_and_frequencies: DMatrix<f64>, // n pools x 1 + p alleles across loci
    pub phenotypes: DMatrix<f64>,                // n pools x k traits
    pub pool_names: Vec<String>,
}

// Struct for GWAlpha's least squares cost function minimisation
#[derive(Debug, Clone)]
pub struct LeastSquaresBeta {
    pub q_prime: DMatrix<f64>,
    pub percs_a: DMatrix<f64>,
    pub percs_b: DMatrix<f64>,
}

// Struct for GWAlpha's maximum likelihood estimation
#[derive(Debug, Clone)]
pub struct MaximumLikelihoodBeta {
    pub percs_a: DMatrix<f64>,
    pub percs_a0: DMatrix<f64>,
    pub percs_b: DMatrix<f64>,
    pub percs_b0: DMatrix<f64>,
}

// Struct for regression objects
#[derive(Debug, Clone)]
pub struct UnivariateOrdinaryLeastSquares {
    pub x: DMatrix<f64>,
    pub y: DVector<f64>,
    pub b: DVector<f64>,
    pub xt: DMatrix<f64>,
    pub inv_xtx: DMatrix<f64>,
    pub inv_xxt: DMatrix<f64>,
    pub e: DVector<f64>,
    pub se: f64,
    pub v_b: DVector<f64>,
    pub t: DVector<f64>,
    pub pval: DVector<f64>,
}

#[derive(Debug, Clone)]
pub struct UnivariateMaximumLikelihoodEstimation {
    pub x: DMatrix<f64>,
    pub y: DVector<f64>,
    pub b: DVector<f64>,
    pub se: f64,
    pub v_b: DVector<f64>,
    pub t: DVector<f64>,
    pub pval: DVector<f64>,
}

#[derive(Debug, Clone)]
pub struct PredictionPerformance {
    pub n: usize,            // number of observations
    pub p: usize,            // number of predictors
    pub k: usize,            // number cross-validation folds
    pub r: usize, // number replications where each cross-validation results in random groupings
    pub models: Vec<String>, // genomic prediction model used
    pub cor: Vec<f64>,
    pub mbe: Vec<f64>,
    pub mae: Vec<f64>,
    pub mse: Vec<f64>,
    pub rmse: Vec<f64>,
}

// Struct for penalised-like regression lambda path search with cross-validation
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct LambdaError {
    pub lambda: f64,
    pub error: f64,
}

// TRAITS
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
        n_threads: &u64,
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
    ) -> io::Result<Vec<LocusFrequencies>>;
    fn load(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &u64,
    ) -> io::Result<Vec<LocusFrequencies>>;
    fn write_csv(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &u64,
    ) -> io::Result<String>;
    fn into_frequencies_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &u64,
    ) -> io::Result<FrequenciesAndPhenotypes>;
}

pub trait Regression {
    fn new() -> Self;
    fn remove_collinearities_in_x(&mut self) -> &mut Self;
    fn estimate_effects(&mut self) -> io::Result<&mut Self>;
    fn estimate_variances(&mut self) -> io::Result<&mut Self>;
    fn estimate_significance(&mut self) -> io::Result<&mut Self>;
}

pub trait EstmateAndPredict<F> {
    fn estimate_effects(&self, function: F) -> io::Result<(DMatrix<f64>, String)>
    where
        F: Fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>;
    fn predict_phenotypes(&self) -> io::Result<DMatrix<f64>>;
}

pub trait CrossValidate<F> {
    fn k_split(&self, k: usize) -> io::Result<(Vec<usize>, usize, usize)>;
    fn performance(&self, y_true: &DMatrix<f64>, y_hat: &DMatrix<f64>) -> io::Result<Vec<f64>>;
    fn cross_validate(
        &self,
        k: usize,
        r: usize,
        functions: Vec<F>,
    ) -> io::Result<PredictionPerformance>
    where
        F: Fn(&DMatrix<f64>, &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)>;
}
