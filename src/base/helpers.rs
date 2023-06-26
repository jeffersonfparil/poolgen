//! Helper functions including parallel selective matrix multiplication and pseudo-inverse via singular value decomposition

use crate::base::*;
use argmin::solver::neldermead::NelderMead;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::svd::*;

use std;
use std::fs::File;
use std::io::{self, prelude::*, BufReader, SeekFrom};
use std::io::{Error, ErrorKind};

/// Find the start position of the next line given the current position,`pos`.
/// Positions are coded as the nth UTF8 character count in the file counting the newline characters at the end of each line.
/// This is used in file splitting to allocate a chunk of the file to a single thread for parallel processing.
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

/// Detect the cursor positions across the input file corresponding to the splits for parallel computation
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

/// Round-up an `f64` to `n_digits` decimal points
pub fn sensible_round(x: f64, n_digits: usize) -> f64 {
    let factor = ("1e".to_owned() + &n_digits.to_string())
        .parse::<f64>()
        .unwrap();
    (x * factor).round() / factor
}

/// Round-up an `f64` to `n_digits` decimal points and cast into a `String`
pub fn parse_f64_roundup_and_own(x: f64, n_digits: usize) -> String {
    let s = x.to_string();
    if s.len() < n_digits {
        return s;
    }
    sensible_round(x, n_digits).to_string()
}

/// Map a vector of `f64` with a logistic regression, so that they are bound between `lower_limit` and `upper_limit`
pub fn bound_parameters_with_logit(
    params: &Vec<f64>,
    lower_limit: f64,
    upper_limit: f64,
) -> Vec<f64> {
    params
        .into_iter()
        .map(|x| lower_limit + ((upper_limit - lower_limit) / (1.00 + (-x).exp())))
        .collect::<Vec<f64>>()
}

/// Instantiate the Nelder-Mead optimiser with `p` parameters an initial values equal defined by a `p+1 x p+1` kernel with `h` in the off-diagonals and `h+0.5 in the diagonals
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

/// Multi-threaded selective matrix multiplication variant: $AB$
/// Matrix multiplication of two 2-dimensional arrays where only the specified rows and columns are used
/// This is an attempt at multi-threadded selective matrix multiplication while minimising memory allocation by using pointers to arrays and a subset of its data, instead of coping data.
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

/// Multi-threaded selective matrix multiplication variant: $A^{T}B$
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

/// Multi-threaded selective matrix multiplication variant: $AB^{T}$
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
    /// Implementation of pseudo-inverse using singlar-value decomposition where $\Sigma_{j,j} < \epsilon_{machine}$ are set to zero
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
        if cfg!(windows) {
            assert_eq!(
                find_start_of_next_line(&"./tests/test.pileup".to_owned(), 56),
                58
            ); // line 1 has a total of 57 characters in windows and only 56 in unix, so the next line should start at the 58th character position
        } else {
            assert_eq!(
                find_start_of_next_line(&"./tests/test.pileup".to_owned(), 56),
                57
            ); // 56th character is the newline character as line 1 has a total of 56 characters, so the next line should start at the 57th character position
        }
        assert_eq!(
            find_file_splits(&"./tests/test.pileup".to_owned(), &2)
                .unwrap()
                .len(),
            3
        );
        assert_eq!(sensible_round(0.420000012435, 4), 0.42);
        assert_eq!(
            parse_f64_roundup_and_own(0.690000012435, 4),
            "0.69".to_owned()
        );
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
        assert_eq!(
            multiply_views_xx(&a, &b, &idx_w3, &idx_x2, &idx_y2, &idx_z2).unwrap(),
            Array2::from_shape_vec((3, 2), vec![27.0, 31.0, 63.0, 73.0, 81.0, 94.0]).unwrap()
        );
        assert_eq!(
            multiply_views_xtx(&a, &b, &idx_w3, &idx_x2, &idx_w3, &idx_z2).unwrap(),
            Array2::from_shape_vec((2, 2), vec![117.0, 129.0, 141.0, 156.0]).unwrap()
        );
        assert_eq!(
            multiply_views_xxt(&a, &b, &idx_w3, &idx_x2, &idx_w3, &idx_z2).unwrap(),
            Array2::from_shape_vec(
                (3, 3),
                vec![14.5, 38.5, 50.5, 35.5, 95.5, 125.5, 46.0, 124.0, 163.0]
            )
            .unwrap()
        );
    }
}
