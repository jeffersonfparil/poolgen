// use std::env;
use clap::Parser;
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
    min_qual: f64,
    /// Minimum depth of coverage
    #[clap(long, default_value_t=1)]
    min_cov: u64,
    /// Minimum allele frequency for associating the genotypes with the phenotype/s
    #[clap(long, default_value_t=0.001)]
    maf: f64,
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
    } else if args.analysis == String::from("pearson_corr") {
        let phen_col = args.phen_value_col.into_iter().map(|x| x.parse::<usize>().expect("Invalid integer input for the phenotype column/s (--phen-value-col).")).collect::<Vec<usize>>();
        let out = regression::correlation(&args.fname,
                                                  &args.maf,
                                                  &args.phen_fname,
                                                  &args.phen_delim,
                                                  &args.phen_name_col,
                                                  &phen_col,
                                                  &args.output,
                                                  &args.n_threads).unwrap();
    } else if args.analysis == String::from("ols_iter") {
        let phen_col = args.phen_value_col.into_iter().map(|x| x.parse::<usize>().expect("Invalid integer input for the phenotype column/s (--phen-value-col).")).collect::<Vec<usize>>();
        let out = regression::ols_iterate(&args.fname,
                                          &args.maf,
                                          &args.phen_fname,
                                          &args.phen_delim,
                                          &args.phen_name_col,
                                          &phen_col,
                                          &args.output,
                                          &args.n_threads).unwrap();
    }else if args.analysis == String::from("test") {
        // let test = io::load_phen(&args.fname, &delim, &);
        // let test = loader_with_fun(&args.fname, &String::From(".sync"), );
        let phen_col: Vec<usize> = vec![1];
        
    }
}
