use ndarray::prelude::*;
use std::io;

///////////////////////////////////////////////////////////////////////////////
// STRUCTS
///////////////////////////////////////////////////////////////////////////////
///
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
    pub phen_matrix: Array2<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FileSyncPhen {
    pub filename_sync: String,
    pub pool_names: Vec<String>,
    pub pool_sizes: Vec<f64>,
    pub phen_matrix: Array2<f64>,
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

// Struct of allele frequencies and phenotypes for genomic prediction
#[derive(Debug, Clone)]
pub struct GenotypesAndPhenotypes {
    pub chromosome: Vec<String>,                       // 1 + p
    pub position: Vec<u64>,                            // 1 + p
    pub allele: Vec<String>,                           // 1 + p
    pub intercept_and_allele_frequencies: Array2<f64>, // n pools x 1 + p alleles across loci
    pub phenotypes: Array2<f64>,                       // n pools x k traits
    pub pool_names: Vec<String>,                       // n
    pub coverages: Array2<f64>,                        // n pools x f(p) loci
}

#[derive(Debug, Clone)]
pub struct PredictionPerformance {
    pub n: usize,                                           // number of observations
    pub p: usize,                                           // number of predictors
    pub k: usize,                                           // number cross-validation folds
    pub r: usize, // number replications where each cross-validation results in random groupings
    pub models: Vec<String>, // genomic prediction model used
    pub y_validation_and_predicted: Array4<f64>, // reps x n x models x traits+traits
    pub b_histogram: Vec<(Vec<f64>, Vec<f64>, Vec<usize>)>, // bin_start (inclusive), bin_end (exclusive), counts or frequency
    pub cor: Array4<f64>,                                   // reps x folds x models x traits
    pub mbe: Array4<f64>,
    pub mae: Array4<f64>,
    pub mse: Array4<f64>,
    pub rmse: Array4<f64>,
}

///////////////////////////////////////////////////////////////////////////////
// TRAITS
///////////////////////////////////////////////////////////////////////////////

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
    fn write_csv(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        out: &String,
        n_threads: &usize,
    ) -> io::Result<String>;
    fn into_genotypes_and_phenotypes(
        &self,
        filter_stats: &FilterStats,
        keep_n_minus_1: bool,
        n_threads: &usize,
    ) -> io::Result<GenotypesAndPhenotypes>;
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
    ) -> io::Result<(String, Vec<String>)>
    where
        F: Fn(&Array2<f64>, &Array2<f64>, &Vec<usize>) -> io::Result<(Array2<f64>, String)>;
}
