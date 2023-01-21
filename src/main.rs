// use std::env;
use clap::Parser;
mod io;

// Instatiate arguments struct
#[derive(Parser, Debug)]
#[clap(author="Jeff Paril",
       version="0.1.0",
       about="Quantitative and population genetics analyses using pool sequencing data.",
       long_about="Quantitative and population genetics analyses using pool sequencing data: trying to continue the legacy of the now unmaintained popoolation2 package with the memory safety of Rust.")]
struct Args {
    /// Analysis to perform (i.e. "pileup2sync", "load")
    analysis: String,
    /// Filename of the input pileup or synchronised pileup file (i.e. *.pileup, *.sync, *.syncf, or *.syncx)
    #[clap(short, long)]
    fname: String,
    /// Text file containing the names of the pools, one name per line
    #[clap(long, default_value="")]
    pool_names: String,
    /// Minimum base quality
    #[clap(long, default_value_t=0.0001)]
    min_qual: f64,
    /// Minimum depth of coverage
    #[clap(long, default_value_t=1)]
    min_cov: u64,
    /// Format of the output file: sync or syncf
    #[clap(long, default_value="sync")]
    file_format: String,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t=1)]
    n_threads: u64,
}

fn main() {
    let args = Args::parse();
    if args.analysis == String::from("pileup2sync") {
        let p = io::read(&args.fname,
                                 &args.pool_names,
                                 &args.min_qual,
                                 &args.min_cov,
                                 &args.file_format,
                                 &args.n_threads).unwrap();
        println!("{:?}", p);

    } else if args.analysis == String::from("load") {
        let x = io::load(&args.fname).unwrap();
        println!("{:?}", x);
    }
}
