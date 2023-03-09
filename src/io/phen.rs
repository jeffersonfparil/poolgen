use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use nalgebra::{self, DMatrix};

#[derive(Debug, Clone)]
pub struct Phenotypes <T, R: nalgebra::Dim, C: nalgebra::Dim> {
    pub name: Vec<String>,
    pub phen: nalgebra::Matrix<T, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<T, R, C>>,
}

pub fn load_phen(fname: &String, delim: &String, header: &bool, name_col: &usize, phen_col: &Vec<usize>) -> io::Result<Phenotypes <f64, nalgebra::Dyn, nalgebra::Dyn>> {
    // // Other parameters
    // let delim: &String = &String::from(",");
    // let header: &bool = &true;
    // let name_col: &usize = &0;
    // let phen_col: &Vec<usize> = &vec![1, 2];

    // Determine the format of the input file
    let file = File::open(fname).unwrap();
    let all_lines: std::iter::Skip<std::io::Lines<BufReader<File>>>;
    if *header {
        all_lines = BufReader::new(file).lines().skip(1);
    } else {
        all_lines = BufReader::new(file).lines().skip(0);
    }
    let mut line: Vec<String>;
    let mut name: Vec<String> = vec![];
    let mut phen: Vec<f64> = vec![];
    for l in all_lines {
        line = l.unwrap().split(delim).collect::<Vec<&str>>().into_iter().map(|x| x.to_owned()).collect::<Vec<String>>();
        // println!("PHENO LINE: {:?}", line);
        name.push(line[*name_col].clone());
        for c in phen_col.iter() {
            let x = line[*c].parse::<f64>().expect(&("Check if column ".to_owned() + &(c + 1).to_string() + " is a valid number or if there are leading or trailing whitespace in line: " + &line.join(",")));
            // println!("PHENO x={:?}", x);
            phen.push(x);
        }
    }
    let mut Y = DMatrix::from_row_slice(name.len(), phen_col.len(), &phen); // generates the matrix row-wise
    let out = Phenotypes{ name: name, phen:  Y};
    // println!("PHENO ROWSUMS = {:?}", out.phen.column_sum());
    // println!("PHENO COLSUMS = {:?}", out.phen.row_sum());
    Ok(out)
}
