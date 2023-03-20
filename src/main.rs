// use std::env;
use clap::Parser;
mod base;
use base::{Parse, ChunkyReadAnalyseWrite};
mod io;
mod tables;
mod regression;
mod gwalpha;

// Instatiate arguments struct
#[derive(Parser, Debug)]
#[clap(author="Jeff Paril",
       version="0.1.0",
       about="Quantitative and population genetics analyses using pool sequencing data.",
       long_about="Quantitative and population genetics analyses using pool sequencing data: trying to continue the legacy of the now unmaintained popoolation2 package with the memory safety of Rust.")]
struct Args {
    /// Analysis to perform (i.e. "pileup2sync", "sync2syncx", "fisher_exact_test", "chisq_test", "pearson_corr", "ols_iter")
    analysis: String,
    /// Filename of the input pileup or synchronised pileup file (i.e. *.pileup, *.sync, *.syncf, or *.syncx)
    #[clap(short, long)]
    fname: String,
    /// Output filename
    #[clap(short, long, default_value="")]
    output: String,
    /// Text file containing the names of the pools, one name per line
    #[clap(long, default_value="")]
    pool_names: String,
    /// Minimum base quality
    #[clap(long, default_value_t=0.01)]
    min_quality: f64,
    /// Minimum depth of coverage
    #[clap(long, default_value_t=1)]
    min_coverage: u64,
    /// Minimum allele frequency for associating the genotypes with the phenotype/s
    #[clap(long, default_value_t=0.001)]
    min_allele_frequency: f64,
    /// Remove ambiguous reads during SNP filetering or keep them coded as Ns
    #[clap(long, default_value_t=true)]
    remove_ns: bool,
    /// Input phenotype file: csv or tsv or any delimited file
    #[clap(long, default_value="")]
    phen_fname: String,
    /// Delimiter of the input phenotype file: comma, tab, etc...
    #[clap(long, default_value=",")]
    phen_delim: String,
    /// Column index containing the names or IDs of the indivudals in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value_t=0)]
    phen_name_col: usize,
    /// Column indexes containing the phenotype values in the input phenotype file, e.g. 1 or 1,2,3 or 1,2,3,4 etc ...
    #[clap(long, use_value_delimiter=true, value_delimiter=',', default_value="1")]
    phen_value_col: Vec<String>,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t=1)]
    n_threads: u64,
}

fn main() {
    let args = Args::parse();
    let mut output: String = String::from("");
    let filter_stats = base::FilterStats{ remove_ns: args.remove_ns, min_quality: args.min_quality, min_coverage: args.min_coverage, min_allele_frequency: args.min_allele_frequency };
    if args.analysis == String::from("pileup2sync") {
        let file_pileup = base::FilePileup{ filename: args.fname, pool_names: args.pool_names };
        output = file_pileup.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        base::pileup_to_sync).unwrap();
    } else if args.analysis == String::from("fisher_exact_test") {
        let file_sync = base::FileSync{ filename: args.fname, test: String::from("fisher_exact_test") };
        output = file_sync.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        tables::fisher).unwrap();
    } else if args.analysis == String::from("chisq_test") {
        let file_sync = base::FileSync{ filename: args.fname, test: String::from("chisq_test") };
        output = file_sync.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        tables::chisq).unwrap();
    } else if args.analysis == String::from("pearson_corr") {
        let file_sync = base::FileSync{ filename: args.fname, test: String::from("pearson_corr") };
        let phen_col = args.phen_value_col.into_iter().map(|x| x.parse::<usize>().expect("Invalid integer input for the phenotype column/s (--phen-value-col).")).collect::<Vec<usize>>();
        let file_phen = base::FilePhen{ filename: args.phen_fname,
                                                  phen_delim: args.phen_delim,
                                                  phen_name_col: args.phen_name_col,
                                                  phen_value_col: phen_col };
        let file_sync_phen = (file_sync, file_phen).lparse().unwrap();
        output = file_sync_phen.read_analyse_write(&filter_stats,
                                                            &args.output,
                                                            &args.n_threads,
                                                            regression::correlation).unwrap();
    } else if args.analysis == String::from("ols_iter") {
        let file_sync = base::FileSync{ filename: args.fname, test: String::from("ols_iter") };
        let phen_col = args.phen_value_col.into_iter().map(|x| x.parse::<usize>().expect("Invalid integer input for the phenotype column/s (--phen-value-col).")).collect::<Vec<usize>>();
        let file_phen = base::FilePhen{ filename: args.phen_fname,
                                                  phen_delim: args.phen_delim,
                                                  phen_name_col: args.phen_name_col,
                                                  phen_value_col: phen_col };
        let file_sync_phen = (file_sync, file_phen).lparse().unwrap();
        output = file_sync_phen.read_analyse_write(&filter_stats,
                                                            &args.output,
                                                            &args.n_threads,
                                                            regression::ols_iterate).unwrap();
    } else if args.analysis == String::from("test") {
        let output = 0;
        println!("TEST={:?}", output);
        
    }
    println!("{}", output);
}
