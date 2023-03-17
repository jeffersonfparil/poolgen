use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader};
use nalgebra::{self, DMatrix};

#[derive(Debug, Clone)]
pub struct Phenotypes <T, R: nalgebra::Dim, C: nalgebra::Dim> {
    pub name: Vec<String>,
    pub min: f64,
    pub max: f64,
    pub sdev: f64,
    pub phen: nalgebra::Matrix<T, nalgebra::Dyn, nalgebra::Dyn, nalgebra::VecStorage<T, R, C>>,
}

pub fn load_phen(fname: &String, delim: &String, summstats: &bool, name_col: &usize, phen_col: &Vec<usize>) -> io::Result<Phenotypes <f64, nalgebra::Dyn, nalgebra::Dyn>> {
    // // Other parameters
    // let delim: &String = &String::from(",");
    // let header: &bool = &true;
    // let name_col: &usize = &0;
    // let phen_col: &Vec<usize> = &vec![1, 2];

    // Determine the format of the input file
    let file = File::open(fname).unwrap();
    let all_lines = BufReader::new(file).lines();
    let mut line: Vec<String>;
    let mut name: Vec<String> = vec![];
    let mut stat: Vec<f64> = vec![];
    let mut phen: Vec<f64> = vec![];
    'lines_loop: for l in all_lines {
        line = l.unwrap()
                .trim()
                .split(delim)
                .collect::<Vec<&str>>()
                .into_iter()
                .map(|x| x.to_owned())
                .collect::<Vec<String>>();
        if line[0].as_bytes()[0] == 35 as u8 {
            if *summstats {
                line[0].remove(0);
                let _ = match line[0].parse::<f64>() {
                    Ok(x) => x,
                    Err(_) => continue 'lines_loop,
                };
                stat.append(&mut line.into_iter()
                           .map(|a| a.parse::<f64>().expect("Error in the line for minimum, maximum , and standard deviaton values of the phenotype"))
                           .collect::<Vec<f64>>());
            }
        } else {
            // println!("PHENO LINE: {:?}", line);
            name.push(line[*name_col].clone());
            for c in phen_col.iter() {
                let x = line[*c].parse::<f64>().expect(&("Check if column ".to_owned() + &(c + 1).to_string() + " is a valid number or if there are leading or trailing whitespace in line: " + &line.join(",")));
                // println!("PHENO x={:?}", x);
                phen.push(x);
            }
        }
    }
    if stat.len() == 0 {
        stat.append(&mut vec![420.69, 420.69, 420.69]);
    }
    let out = Phenotypes{ name: name.to_owned(),
                                                     min: stat[0],
                                                     max: stat[1],
                                                     sdev: stat[2],
                                                     phen: DMatrix::from_row_slice(name.len(),
                                                                                   phen_col.len(), 
                                                                                   &phen) // generates the matrix row-wise
                                                    };
    // println!("PHENO ROWSUMS = {:?}", out.phen.column_sum());
    // println!("PHENO COLSUMS = {:?}", out.phen.row_sum());
    Ok(out)
}
