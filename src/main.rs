// use std::env;
use clap::Parser;
mod io;
mod tables;
mod regression;

// Instatiate arguments struct
#[derive(Parser, Debug)]
#[clap(author="Jeff Paril",
       version="0.1.0",
       about="Quantitative and population genetics analyses using pool sequencing data.",
       long_about="Quantitative and population genetics analyses using pool sequencing data: trying to continue the legacy of the now unmaintained popoolation2 package with the memory safety of Rust.")]
struct Args {
    /// Analysis to perform (i.e. "pileup2sync", "sync2syncx", "load", "fisher_exact_test")
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
    min_qual: f64,
    /// Minimum depth of coverage
    #[clap(long, default_value_t=1)]
    min_cov: u64,
    /// Remove ambiguous reads during SNP filetering or keep them coded as Ns
    #[clap(long, default_value_t=true)]
    remove_ns: bool,
    /// Format of the output file: sync or syncx
    #[clap(long, default_value="sync")]
    file_format: String,
    /// Input phenotype file: csv or tsv or any delimited file
    #[clap(long, default_value="")]
    phen_fname: String,
    /// Delimiter of the input phenotype file: comma, tab, etc...
    #[clap(long, default_value=",")]
    phen_delim: String,
    /// Does the input phenotype file have a header?
    #[clap(long, default_value_t=true)]
    phen_header: bool,
    /// Column index containing the names or IDs of the indivudals in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value_t=0)]
    phen_name_col: usize,
    /// Column indexes containing the phenotype values in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value="[1]")]
    phen_phen_col: String,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t=1)]
    n_threads: u64,
}

fn main() {
    let args = Args::parse();
    if args.analysis == String::from("pileup2sync") {
        let out: String = io::pileup2sync(&args.fname,
                                          &args.output,
                                          &args.pool_names,
                                          &args.min_qual,
                                          &args.min_cov,
                                          &args.remove_ns,
                                          &args.file_format,
                                          &args.n_threads).unwrap();
        println!("{:?}", out);

    } else if args.analysis == String::from("sync2syncx") {
        let out: String = io::sync2syncx(&args.fname,
                                        &args.output,
                                         &args.min_cov,
                                         &args.n_threads).unwrap();
        println!("{:?}", out);
    } else if args.analysis == String::from("fisher_exact_test") {
        let out = tables::fisher(&args.fname, &args.output, &args.n_threads).unwrap();
    } else if args.analysis == String::from("chisq_test") {
        let out = tables::chisq(&args.fname, &args.output, &args.n_threads).unwrap();
    } else if args.analysis == String::from("test") {
        // let test = io::load_phen(&args.fname, &delim, &);
        // let test = loader_with_fun(&args.fname, &String::From(".sync"), );

        let phen_col: Vec<usize> = vec![1];
        let test = regression::correlation(&args.fname, &args.phen_fname, &args.phen_delim, &args.phen_header, &args.phen_name_col,
                &phen_col, &args.output, &args.n_threads);
        println!("TEST: {:?}", test);
    }
}
