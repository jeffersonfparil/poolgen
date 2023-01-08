mod io;

fn main() {
    let fname: &str = "./tests/test.pileup";
    let p = io::read(fname).unwrap();
    for i in p.into_iter() {
        println!("{:?}", i);
    }
}
