use crate::base::*;
use argmin::solver::neldermead::NelderMead;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::svd::*;

use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, SeekFrom};
use std::io::{Error, ErrorKind};

// File splitting for thread allocation for parallele computation
fn find_start_of_next_line(fname: &String, pos: u64) -> u64 {
    let mut out = pos.clone();
    if out > 0 {
        let mut file = File::open(fname).unwrap();
        let _ = file.seek(SeekFrom::Start(out));
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let _ = reader.read_line(&mut line).unwrap();
        out = reader.seek(SeekFrom::Current(0)).unwrap();
    }
    return out;
}

pub fn find_file_splits(fname: &String, n_threads: &usize) -> io::Result<Vec<u64>> {
    let mut file = match File::open(fname) {
        Ok(x) => x,
        Err(_) => return Err(Error::new(ErrorKind::Other, "The input file: ".to_owned() + fname + " does not exist. Please make sure you are entering the correct filename and/or the correct path.")),
    };
    let _ = file.seek(SeekFrom::End(0));
    let mut reader = BufReader::new(file);
    let end = reader.seek(SeekFrom::Current(0)).unwrap();
    let mut out = (0..end)
        .step_by((end as usize) / n_threads)
        .collect::<Vec<u64>>();
    out.push(end);
    for i in 0..out.len() {
        out[i] = find_start_of_next_line(fname, out[i]);
    }
    out.dedup();
    return Ok(out);
}

pub fn sensible_round(x: f64, n_digits: usize) -> f64 {
    let factor = ("1e".to_owned() + &n_digits.to_string())
        .parse::<f64>()
        .unwrap();
    (x * factor).round() / factor
}

pub fn parse_f64_roundup_and_own(x: f64, n_digits: usize) -> String {
    let s = x.to_string();
    if s.len() < n_digits {
        return s;
    }
    sensible_round(x, n_digits).to_string()
}

pub fn bound_parameters_with_logit(
    params: &Vec<f64>,
    lower_limit: f64,
    upper_limit: f64,
) -> Vec<f64> {
    // Map parameters with a logistic regression to bound them between 0 and PARAMETER_UPPER_LIMIT
    params
        .into_iter()
        .map(|x| lower_limit + ((upper_limit - lower_limit) / (1.00 + (-x).exp())))
        .collect::<Vec<f64>>()
}

pub fn prepare_solver_neldermead(p: f64, h: f64) -> NelderMead<Vec<f64>, f64> {
    let mut init_param: Vec<Vec<f64>> = Vec::new();
    for i in 0..(p as usize + 1) {
        init_param.push(vec![]);
        for j in 0..p as usize {
            if i == j {
                // init_param[i].push(1.5 * h)
                init_param[i].push(h + 0.5)
            } else {
                init_param[i].push(h)
            }
        }
    }
    NelderMead::new(init_param)
}

pub fn multiply_views_xx(
    a: &Array2<f64>,
    b: &Array2<f64>,
    a_rows: &Vec<usize>,
    a_cols: &Vec<usize>,
    b_rows: &Vec<usize>,
    b_cols: &Vec<usize>,
) -> io::Result<Array2<f64>> {
    let n = a_rows.len();
    let m = b_cols.len();
    if a_cols.len() != b_rows.len() {
        return Err(Error::new(
            ErrorKind::Other,
            "The two matrices are incompatible.",
        ));
    }
    let mut out: Array2<f64> = Array2::zeros((n, m));
    let a_rows_mat = Array2::from_shape_vec((m, n), a_rows.repeat(m))
        .unwrap()
        .reversed_axes();
    let b_cols_mat = Array2::from_shape_vec((n, m), b_cols.repeat(n)).unwrap();
    Zip::from(&mut out)
        .and(&a_rows_mat)
        .and(&b_cols_mat)
        .par_for_each(|x, &a_i, &b_j| {
            for k in 0..a_cols.len() {
                let a_j = a_cols[k];
                let b_i = b_rows[k];
                *x += a[(a_i, a_j)] * b[(b_i, b_j)];
            }
        });
    Ok(out)
}

pub fn multiply_views_xtx(
    a: &Array2<f64>,
    b: &Array2<f64>,
    a_rows: &Vec<usize>,
    a_cols: &Vec<usize>,
    b_rows: &Vec<usize>,
    b_cols: &Vec<usize>,
) -> io::Result<Array2<f64>> {
    let n = a_cols.len(); // reversed a
    let m = b_cols.len();
    if a_rows.len() != b_rows.len() {
        // reversed a
        return Err(Error::new(
            ErrorKind::Other,
            "The two matrices are incompatible.",
        ));
    }
    let mut out: Array2<f64> = Array2::zeros((n, m));
    let a_cols_mat = Array2::from_shape_vec((m, n), a_cols.repeat(m))
        .unwrap()
        .reversed_axes();
    let b_cols_mat = Array2::from_shape_vec((n, m), b_cols.repeat(n)).unwrap();
    Zip::from(&mut out)
        .and(&a_cols_mat)
        .and(&b_cols_mat)
        .par_for_each(|x, &a_j, &b_j| {
            for k in 0..a_rows.len() {
                let a_i = a_rows[k];
                let b_i = b_rows[k];
                *x += a[(a_i, a_j)] * b[(b_i, b_j)];
            }
        });
    Ok(out)
}

pub fn multiply_views_xxt(
    a: &Array2<f64>,
    b: &Array2<f64>,
    a_rows: &Vec<usize>,
    a_cols: &Vec<usize>,
    b_rows: &Vec<usize>,
    b_cols: &Vec<usize>,
) -> io::Result<Array2<f64>> {
    let n = a_rows.len();
    let m = b_rows.len(); // reversed b
    if a_cols.len() != b_cols.len() {
        // reversed b
        return Err(Error::new(
            ErrorKind::Other,
            "The two matrices are incompatible.",
        ));
    }
    let mut out: Array2<f64> = Array2::zeros((n, m));
    let a_rows_mat = Array2::from_shape_vec((m, n), a_rows.repeat(m))
        .unwrap()
        .reversed_axes();
    let b_rows_mat = Array2::from_shape_vec((n, m), b_rows.repeat(n)).unwrap();
    Zip::from(&mut out)
        .and(&a_rows_mat)
        .and(&b_rows_mat)
        .par_for_each(|x, &a_i, &b_i| {
            for k in 0..a_cols.len() {
                let a_j = a_cols[k];
                let b_j = b_cols[k];
                *x += a[(a_i, a_j)] * b[(b_i, b_j)];
            }
        });
    Ok(out)
}

impl MoorePenrosePseudoInverse for Array2<f64> {
    fn pinv(&self) -> io::Result<Array2<f64>> {
        let n: usize = self.nrows();
        let p: usize = self.ncols();
        let svd = self.svd(true, true).unwrap();
        let u: Array2<f64> = svd.0.unwrap();
        let s: Array1<f64> = svd.1;
        let vt: Array2<f64> = svd.2.unwrap();
        let mut s_inv: Array2<f64> = Array2::zeros((n, p));
        let tolerance: f64 =
            f64::EPSILON * (s.len() as f64) * s.fold(s[0], |max, &x| if x > max { x } else { max });
        for j in 0..p {
            if s[j] > tolerance {
                s_inv[(j, j)] = 1.0 / s[j];
            }
        }
        let self_inv: Array2<f64> = vt.t().dot(&s_inv.t()).dot(&u.t());
        Ok(self_inv)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;
    #[test]
    fn test_helpers() {
        // Inputs
        let fname: &String = &"./tests/test.pileup".to_owned();
        let n_threads = 2;
        let number: f64 = 0.420000012435;
        let _betas = vec![
            0.0, 0.2, 0.3, 0.4, 0.5, 0.55, 0.6, 0.7, 0.75, 0.77, 0.8, 0.85, 0.86, 0.89, 1.0,
        ];
        let a: Array2<f64> =
            Array2::from_shape_vec((5, 3), (0..15).map(|x| x as f64).collect::<Vec<f64>>())
                .unwrap();
        let b: Array2<f64> = Array2::from_shape_vec(
            (5, 3),
            (0..15).map(|x| x as f64 / 2.0).collect::<Vec<f64>>(),
        )
        .unwrap();
        let idx_w3: Vec<usize> = vec![1, 3, 4];
        let idx_x2: Vec<usize> = vec![0, 2];
        let idx_y2: Vec<usize> = vec![1, 3];
        let idx_z2: Vec<usize> = vec![0, 1];
        // Outputs
        let splits = find_file_splits(fname, &n_threads).unwrap();
        let string_f64 = parse_f64_roundup_and_own(number, 4);
        let a_x_b = multiply_views_xx(&a, &b, &idx_w3, &idx_x2, &idx_y2, &idx_z2).unwrap();
        let at_x_b = multiply_views_xtx(&a, &b, &idx_w3, &idx_x2, &idx_w3, &idx_z2).unwrap();
        let a_x_bt = multiply_views_xxt(&a, &b, &idx_w3, &idx_x2, &idx_w3, &idx_z2).unwrap();

        // let (n, p) = (100, 50_000);
        // let a: Array2<f64> =
        //     Array2::from_shape_vec((n, p), (0..(n * p)).map(|x| x as f64).collect::<Vec<f64>>())
        //         .unwrap();
        // let b: Array2<f64> = Array2::from_shape_vec(
        //     (p, n),
        //     (0..(p * n)).map(|x| x as f64 / 2.0).collect::<Vec<f64>>(),
        // )
        // .unwrap();
        // let idx_w3: Vec<usize> = (0..n).collect();
        // let idx_x2: Vec<usize> = (0..p).collect();
        // let idx_y2: Vec<usize> = (0..p).collect();
        // let idx_z2: Vec<usize> = (0..n).collect();

        // let a_x_b_multhread =
        //     MULTITHREAD_multiply_views_xx(&a, &b, &idx_w3, &idx_x2, &idx_y2, &idx_z2).unwrap();
        // // println!("a_x_b_multhread={:?}", a_x_b_multhread);
        // let a_x_b = multiply_views_xx(&a, &b, &idx_w3, &idx_x2, &idx_y2, &idx_z2).unwrap();
        // // println!("a_x_b={:?}", a_x_b);
        // assert_eq!(a_x_b_multhread, a_x_b);

        // Assertion
        assert_eq!(splits.len(), 3);
        assert_eq!(string_f64, "0.42".to_owned());
        assert_eq!(
            a_x_b,
            Array2::from_shape_vec((3, 2), vec![27.0, 31.0, 63.0, 73.0, 81.0, 94.0]).unwrap()
        );
        assert_eq!(
            at_x_b,
            Array2::from_shape_vec((2, 2), vec![117.0, 129.0, 141.0, 156.0]).unwrap()
        );
        assert_eq!(
            a_x_bt,
            Array2::from_shape_vec(
                (3, 3),
                vec![14.5, 38.5, 50.5, 35.5, 95.5, 125.5, 46.0, 124.0, 163.0]
            )
            .unwrap()
        );
    }
}
