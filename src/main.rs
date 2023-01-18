use std::env;
mod io;

fn main() {
    let args: Vec<String> = env::args().collect();
    let fname: String = args[1].to_owned();
    let min_qual: f64 = 0.001;
    let min_cov: u64 = 1;
    let p = io::read(&fname, &min_qual, &min_cov).unwrap();
    // for i in p.into_iter() {
    //     println!("{:?}", i);
    // }
}
