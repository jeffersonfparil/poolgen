use std::io;
use nalgebra::DMatrix;

// STRUCTS
#[derive(Debug, Clone)]
pub struct FilePileup {
    pub filename: String,
    pub pool_names: String,
}

#[derive(Debug, Clone)]
pub struct FileSync {
    pub filename: String,
    pub test: String,
}

#[derive(Debug, Clone)]
pub struct FilePhen {
    pub filename: String,
    pub phen_delim: String,
    pub phen_name_col: usize,
    pub phen_value_col: Vec<usize>,
    pub format: String,
}

#[derive(Debug, Clone)]
pub struct FileSyncPhen {
    pub filename_sync: String,
    pub pool_names: Vec<String>,
    pub phen_matrix: DMatrix<f64>,
    pub test: String,
}

#[derive(Debug, Clone)]
pub struct FilterStats {
    pub remove_ns: bool,
    pub min_quality: f64,
    pub min_coverage: u64,
    pub min_allele_frequency: f64,
}

#[derive(Debug, Clone)]
pub struct PileupLine {
    pub chromosome: String,             // chromosome or scaffold name
    pub position: u64,                  // position in number of bases
    pub reference_allele: char,         // reference allele
    pub coverages: Vec<u64>,            // number of times the locus was covered
    pub read_codes: Vec<Vec<u8>>,       // utf8 read codes corresponding to 'A', 'T', 'C', or 'G' (252 other alleles can be accommodated)
    pub read_qualities: Vec<Vec<u8>>,   // utf8 base quality codes which can be transformed into bases error rate as 10^(-(u8 - 33)/10)
}

// Struct of allele counts to convert reads into sync
#[derive(Debug, Clone)]
pub struct LocusCounts {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: DMatrix<u64>, // n pools x p alleles
}

// Struct of allele frequencies to convert reads into syncx
#[derive(Debug, Clone)]
pub struct LocusFrequencies {
    pub chromosome: String,
    pub position: u64,
    pub alleles_vector: Vec<String>,
    pub matrix: DMatrix<f64>, // n pools x p alleles
}

// Struct for tracking the insertions and deletions prepended with +/i which we want to skip as they refer to the next base position and not the current locus
#[derive(Debug, Clone)]
pub struct IndelMarker {
    pub indel: bool,    // insertion or deletion, i.e. +/- codes
    pub count: usize,   // size of the indel, i.e. how many bases
    pub left: usize,    // how many indel bases are left to be removed (initialised with the maximum possible value of 4294967295)
}

// Struct of allele counts and phenotypes per pool
#[derive(Debug, Clone)]
pub struct LocusCountsAndPhenotypes {
    pub locus_counts: LocusCounts,
    pub phenotypes: DMatrix<f64>, // n pools x k traits
    pub pool_names: Vec<String>,
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

pub trait ChunkyReadAnalyseWrite<T, F> {
    fn per_chunk(&self, start: &u64, end: &u64, outname_ndigits: &usize, filter_stats: &FilterStats, function: F) -> io::Result<String>
        where F: Fn(&mut T, &FilterStats) -> Option<String>;
    fn read_analyse_write(&self, filter_stats: &FilterStats, out: &String, n_threads: &u64, function: F) -> io::Result<String>
        where F: Fn(&mut T, &FilterStats) -> Option<String>;
}