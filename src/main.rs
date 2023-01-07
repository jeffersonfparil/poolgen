mod io;

fn main() {
    let fname: &str = "./tests/test.pileup";
    let p = io::read(fname);
    println!("{:?}", p);
}
