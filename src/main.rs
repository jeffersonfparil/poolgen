mod io;

fn main() {
    let fname: &str = "/home/jeffersonfparil/Documents/poolgen/tests/test.pileup";
    let p = io::read(fname);
    println!("{:?}", p);
}
