// For quick tests
#[allow(warnings)]

use clap::Parser;
mod base;
use base::{Parse, ChunkyReadAnalyseWrite};
mod tables;
mod gwas;

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
    #[clap(short, long)]
    phen_fname: String,
    /// Delimiter of the input phenotype file: comma, tab, etc...
    #[clap(long, default_value=",")]
    phen_delim: String,
    /// Column index containing the names or IDs of the indivudals in the input phenotype file: 0, 1, 2, ...
    #[clap(long, default_value_t=0)]
    phen_name_col: usize,
    /// Column index containing the the sizes of each pool or population: 0, 1, 2, ...
    #[clap(long, default_value_t=1)]
    phen_pool_size_col: usize,
    /// Column indexes containing the phenotype values in the input phenotype file, e.g. 1 or 1,2,3 or 1,2,3,4 etc ...
    #[clap(long, use_value_delimiter=true, value_delimiter=',', default_value="2")]
    phen_value_col: Vec<String>,
    /// Number of threads to use for parallel processing
    #[clap(long, default_value_t=1)]
    n_threads: u64,
    ////////////////////////////////////////////////////
    ////// Additional parameters
    ////////////////////////////////////////////////////
    /// GWAlpha inference method to use: "LS" for least squares or "ML" for maximum likelihood estimation
    #[clap(long, default_value="ML")]
    gwalpha_method: String,
}

fn main() {
    let args = Args::parse();
    let mut output: String = String::from("");
    // Prepare the mandatory inputs
    let mut phen_format = "default".to_string();
    if  args.analysis == String::from("gwalpha") {
        phen_format = "gwalpha_fmt".to_string()
    }
    let phen_col = args.phen_value_col.into_iter().map(|x| x.parse::<usize>().expect("Invalid integer input for the phenotype column/s (--phen-value-col).")).collect::<Vec<usize>>();
    let file_phen = base::FilePhen{ filename: args.phen_fname.clone(),
                                                delim: args.phen_delim.clone(),
                                                names_column_id: args.phen_name_col,
                                                sizes_column_id: args.phen_pool_size_col,
                                                trait_values_column_ids: phen_col.clone(),
                                                format: phen_format };
    
    // let file_sync = base::FileSync{ filename: args.fname.clone(), test: args.analysis.clone() }; // Note: file may be pileup instead of sync, but this does not matter, becuase if the input file is pileup then we only need to convert the pileup into sync and we only need the pool sizes from that
    // let file_sync_phen = *(file_sync.clone(), file_phen).lparse().unwrap(); // Note: file may be pileup instead of sync, but this does not matter, becuase if the input file is pileup then we only need to convert the pileup into sync and we only need the pool sizes from that

    let phen = file_phen.lparse().unwrap();

    let filter_stats = base::FilterStats{ remove_ns: args.remove_ns,
                                                       min_quality: args.min_quality,
                                                       min_coverage: args.min_coverage,
                                                       min_allele_frequency: args.min_allele_frequency,
                                                       pool_sizes: phen.pool_sizes.clone() };
    if args.analysis == String::from("pileup2sync") {
        let file_pileup = base::FilePileup{ filename: args.fname, pool_names: phen.pool_names };
        output = file_pileup.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        base::pileup_to_sync).unwrap();
    } else if args.analysis == String::from("fisher_exact_test") {
        let file_sync = base::FileSync{ filename: args.fname, test: args.analysis };
        output = file_sync.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        tables::fisher).unwrap();
    } else if args.analysis == String::from("chisq_test") {
        let file_sync = base::FileSync{ filename: args.fname, test: args.analysis };
        output = file_sync.read_analyse_write(&filter_stats,
                                                         &args.output,
                                                         &args.n_threads,
                                                        tables::chisq).unwrap();
    } else if args.analysis == String::from("pearson_corr") {
        let file_sync = base::FileSync{ filename: args.fname, test: args.analysis };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        output = file_sync_phen.read_analyse_write(&filter_stats,
                                                            &args.output,
                                                            &args.n_threads,
                                                            gwas::correlation).unwrap();
    } else if args.analysis == String::from("ols_iter") {
        let file_sync = base::FileSync{ filename: args.fname, test: args.analysis };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        output = file_sync_phen.read_analyse_write(&filter_stats,
                                                            &args.output,
                                                            &args.n_threads,
                                                            gwas::ols_iterate).unwrap();
    } else if args.analysis == String::from("gwalpha") {
        // Redefine combined sync and phenotype struct under GWAlpha analysis
        let file_sync = base::FileSync{ filename: args.fname, test: args.analysis };
        let file_sync_phen = *(file_sync, file_phen).lparse().unwrap();
        if args.gwalpha_method == "LS".to_owned() {
            output = file_sync_phen.read_analyse_write(&filter_stats,
                                                        &args.output,
                                                        &args.n_threads,
                                                        gwas::gwalpha_ls).unwrap()
        } else {
            // Defaut is ML, i.e. maximum likelihood estimation
            output = file_sync_phen.read_analyse_write(&filter_stats,
                                                        &args.output,
                                                        &args.n_threads,
                                                        gwas::gwalpha_ml).unwrap()
        }
    } else if args.analysis == String::from("test") {
        let output = 0;
        println!("TEST={:?}", output);
        
    }
    println!("{}", output);
}
