// use std::env;
use clap::Parser;
mod io;

// Instatiate arguments struct
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
   /// Input name of the pileup file
   #[clap(short, long)]
   fname: String,
   /// Text file containing the names of the pools, one name per line
   #[clap(short, long)]
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
    let p = io::read(&args.fname,
                     &args.pool_names,
                     &args.min_qual,
                     &args.min_cov,
                     &args.file_format,
                     &args.n_threads).unwrap();
    println!("{:?}", p);
}
