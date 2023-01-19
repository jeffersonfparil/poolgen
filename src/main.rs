use std::env;
mod io;

fn main() {
    let args: Vec<String> = env::args().collect();
    let fname: String = args[1].to_owned();
    let min_qual: f64 = 0.001;
    let min_cov: u64 = 1;
    // let file_format: String = "sync".to_owned();
    let file_format: String = "syncf".to_owned();
    let n_threads: u64 = 3;
    let p = io::read(&fname, &min_qual, &min_cov, &file_format, &n_threads).unwrap();
    // for i in p.into_iter() {
    //     println!("{:?}", i);
    // }
}
