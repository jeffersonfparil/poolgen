mod io;

fn main() {
    let fname: &str = "/home/jeffersonfparil/Documents/poolgen/tests/test.pileup";
    let out = io::read(fname);
    println!("{:?}", out);
}
