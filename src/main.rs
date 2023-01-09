mod io;

fn main() {
    let fname: &str = "./tests/test.pileup";
    let min_qual: f64 = 0.001;
    let p = io::read(fname, &min_qual).unwrap();
    for i in p.into_iter() {
        println!("{:?}", i);
    }
}
