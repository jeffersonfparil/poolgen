//! Helper functions including parallel selective matrix multiplication and pseudo-inverse via singular value decomposition

use crate::base::*;
use argmin::solver::neldermead::NelderMead;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::svd::*;

use std::fs::{self, File};
use std::io::{self, prelude::*, BufReader, SeekFrom, Error, ErrorKind};
use std::path::PathBuf;
use std::process::Command;

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

/// Run python scripts and append output file names to output string
pub fn run_python_and_append(output: &str, script_names: &[&str]) -> String {
    let abs_output = fs::canonicalize(output).expect("Failed to canonicalize output path");
    let scripts_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/python");

    let outputs: Vec<String> = std::iter::once(output.to_string())
        .chain(script_names.iter().map(|script_name| {
            let abs_script = fs::canonicalize(scripts_dir.join(script_name))
                .unwrap_or_else(|_| panic!("Failed to find script {}", script_name));

            let output = Command::new("python3")
                .arg(&abs_script)
                .arg(&abs_output)
                .output()
                .unwrap_or_else(|_| panic!("Failed to run {}", script_name));

            if !output.status.success() { panic!( "{} failed:\n{}", script_name, String::from_utf8_lossy(&output.stderr)) }

            String::from_utf8_lossy(&output.stdout).to_string()
        }))
        .collect();

    outputs.join("\n")
}

/// Run general python script
pub fn run_python(output: &str, script_names: &[&str]) {
    let abs_output = fs::canonicalize(output).expect("Failed to canonicalize output path");
    let scripts_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/python");

    for script_name in script_names {
        let abs_script = fs::canonicalize(scripts_dir.join(script_name))
            .unwrap_or_else(|_| panic!("Failed to find script {}", script_name));

        let output = Command::new("python3")
            .arg(&abs_script)
            .arg(&abs_output)
            .output()
            .unwrap_or_else(|_| panic!("Failed to run {}", script_name));

        if !output.status.success() { panic!( "{} failed:\n{}", script_name, String::from_utf8_lossy(&output.stderr)) }
    }
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

pub fn parse_valid_freq(value: &str) -> Result<f64, String> {
    let parsed: f64 = value.parse().map_err(|_| format!("`{}` isn't a valid number", value))?;
    if parsed < 0.0 || parsed > 1.0 {
        Err(format!("Value must be between 0.0 and 1.0, got `{}`", parsed))
    } else {
        Ok(parsed)
    }
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

/// Calculate the mean of a 1D array ignoring NaN
pub fn mean_array1_ignore_nan(
    x: &ArrayBase<ndarray::ViewRepr<&f64>, Dim<[usize; 1]>>,
) -> io::Result<f64> {
    let sum = x.fold(0.0, |sum, &a| if a.is_nan() { sum } else { sum + a });
    let counts = x.iter().filter(|&a| !a.is_nan()).count() as f64;
    Ok(sum / counts)
}

/// Calculate the axis-wise means of an array while ignoring NaN
pub fn mean_axis_ignore_nan<D>(
    a: &Array<f64, D>,
    axis: usize,
) -> io::Result<Array<f64, <D>::Smaller>>
where
    D: ndarray::Dimension + ndarray::RemoveAxis,
{
    let sum: Array<f64, <D>::Smaller> =
        a.fold_axis(
            Axis(axis),
            0.0,
            |&sum, &x| {
                if x.is_nan() {
                    sum
                } else {
                    sum + x
                }
            },
        );
    let counts: Array<f64, <D>::Smaller> = a.map_axis(Axis(axis), |x| {
        x.iter().filter(|&&y| !y.is_nan()).count() as f64
    });
    let out: Array<f64, <D>::Smaller> = sum / counts;
    Ok(out)
}

/// Extract the coordinates of each sliding window (can accommodate redundant and non-redundant loci)
pub fn define_sliding_windows(
    loci_chr: &Vec<String>,
    loci_pos: &Vec<u64>,
    window_size_bp: &u64,
    window_slide_size_bp: &u64,
    min_loci_per_window: &u64,
) -> io::Result<(Vec<usize>, Vec<usize>)> {
    assert_eq!(loci_chr.len(), loci_pos.len());
    let l = loci_chr.len();
    // Indices, chromosome names, and positions of the start and end of the window, respectively (will be filtered to remove redundant tails to remove windows which are complete subsets of a bigger window near the end of chromosomes or scaffolds)
    let mut idx_head: Vec<usize> = vec![0];
    let mut idx_tail = idx_head.clone();
    let mut chr_head: Vec<String> = vec![loci_chr[0].to_owned()];
    let mut chr_tail = chr_head.clone();
    let mut pos_head: Vec<u64> = vec![loci_pos[0]];
    let mut pos_tail = pos_head.clone();
    // Number of genotyped loci included per window
    let mut cov: Vec<u64> = vec![1];
    // A boolean to mark whether the start of the next window has been found before the end of the current window has been reached
    let mut marker_next_window_head: bool = false;
    // The index of the sart of the next window
    let mut idx_next_head: usize = 0;
    // Index of the current locus
    let mut i: usize = 1;
    // We iterate across the genome until we reach the last locus (may not be consecutive as the start of the next window may be within the body of the current window)
    while i < l {
        let chr = loci_chr[i].to_owned();
        let pos = loci_pos[i];
        // println!("i={:?}", i);
        // println!("idx_head={:?}", idx_head);
        // println!("idx_tail={:?}", idx_tail);
        // Did we reach the end of the chromosome or the end of the window according to window size?
        if (&chr != chr_head.last().unwrap()) | (pos > (pos_head.last().unwrap() + window_size_bp))
        {
            // If we found the start of the next window in body of the current (ending) window then,
            //  we use the next window head as the start of the next slide not the end of the window,
            //  otherwise we use the end of the window.
            i = if marker_next_window_head {
                idx_next_head
            } else {
                i
            };
            let chr = loci_chr[i].to_owned();
            let pos = loci_pos[i];
            // Do we have the minimum number of required loci in the current window?
            if cov.last().unwrap() >= min_loci_per_window {
                // If we have enough loci covered in the current (ending) window:
                // We also add the details of the start of the next window
                idx_head.push(i);
                idx_tail.push(i);
                chr_head.push(chr.to_owned());
                chr_tail.push(chr.to_owned());
                pos_head.push(pos);
                pos_tail.push(pos);
                cov.push(1);
            } else {
                // If we did no have enough loci covered in the current (ending) window:
                // We ditch the current (ending) window and replace it with the start of the next window
                let i_ = idx_head.len() - 1;
                idx_head[i_] = i;
                chr_head[i_] = chr;
                pos_head[i_] = pos;
                cov[i_] = 1;
            }
            // Reset the marker for the start of the next window
            marker_next_window_head = false;
        } else {
            // If we have yet to reach the end of the current window or the end of the chromosome,
            // then we just replace the tail of the current window with the current locus
            // and add the another locus to the coverage counter.
            let i_ = idx_tail.len() - 1;
            idx_tail[i_] = i;
            chr_tail[i_] = chr;
            pos_tail[i_] = pos;
            cov[i_] += 1;
            // We also check if we have reached the start of the next window and note the index if we have
            if (marker_next_window_head == false)
                & (pos >= (pos_head.last().unwrap() + window_slide_size_bp))
            {
                marker_next_window_head = true;
                idx_next_head = i;
            }
        }
        // Move to the next locus
        i += 1;
    }
    // Remove redundant tails
    assert_eq!(idx_head.len(), idx_tail.len());
    let n = idx_head.len();
    let mut out_idx_head: Vec<usize> = vec![idx_head[0]];
    let mut out_idx_tail: Vec<usize> = vec![idx_tail[0]];
    for i in 1..n {
        // println!("out_idx_tail={:?}", out_idx_tail);
        if &idx_tail[i] != out_idx_tail.last().unwrap() {
            out_idx_head.push(idx_head[i]);
            out_idx_tail.push(idx_tail[i]);
        }
    }
    // println!("#################################");
    // println!("loci_chr={:?}", loci_chr);
    // println!("loci_pos={:?}", loci_pos);
    // println!("window_size_bp={:?}", window_size_bp);
    // println!("min_loci_per_window={:?}", min_loci_per_window);
    // println!("cov={:?}", cov);
    // println!("idx_head={:?}", idx_head);
    // println!("idx_tail={:?}", idx_tail);
    // println!("out_idx_head={:?}", out_idx_head);
    // println!("out_idx_tail={:?}", out_idx_tail);
    Ok((out_idx_head, out_idx_tail))
}

/// Load table from a delimited text file
pub fn load_table(
    fname: &String,
    delimiter: &String,
    idx_row_labels: &Vec<usize>,
    data_start_col: &usize,
    data_end_col: &usize,
) -> io::Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>)> {
    let file = File::open(fname).unwrap();
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let column_labels = match lines.next() {
        Some(x) => x
            .unwrap()
            .split(delimiter.as_str())
            .map(|x| x.to_owned())
            .collect::<Vec<String>>(),
        None => return Err(Error::new(ErrorKind::Other, "No lines found.")),
    };
    let data_end_col = if column_labels.len() < *data_end_col {
        column_labels.len()
    } else {
        *data_end_col
    };
    let column_labels = column_labels[*data_start_col..data_end_col]
        .into_iter()
        .map(|x| x.to_owned())
        .collect::<Vec<String>>();
    let mut row_labels: Vec<String> = vec![];
    let mut data: Vec<Vec<f64>> = vec![];
    while let Some(line) = lines.next() {
        let mut line = line.unwrap();
        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }
        let line = line.split(delimiter.as_str()).collect::<Vec<&str>>();
        let mut lab = vec![];
        for i in idx_row_labels {
            lab.push(line[*i].to_owned())
        }
        row_labels.push(lab.join("__-__"));
        data.push(
            line[*data_start_col..data_end_col]
                .into_iter()
                .map(|x| match x.parse::<f64>() {
                    Ok(x) => x,
                    Err(_) => f64::NAN,
                })
                .collect::<Vec<f64>>(),
        );
    }
    Ok((row_labels, column_labels, data))
}

/// Implementation of pseudo-inverse using singlar-value decomposition where $\Sigma_{j,j} < \epsilon_{machine}$ are set to zero
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
        // Define sliding windows (non-redundant loci, i.e. per locus list with alleles ID removed)
        let loci_chr: Vec<String> = vec![
            "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr2", "chr3",
        ]
        .iter()
        .map(|&x| x.to_owned())
        .collect();
        let loci_pos: Vec<u64> = vec![123, 174, 220, 254, 55, 56, 100, 500, 765];
        let window_size_bp: u64 = 100;
        let window_slide_size_bp: u64 = 50;
        let min_loci_per_window: u64 = 1;
        let (windows_idx_head, windows_idx_tail) = define_sliding_windows(
            &loci_chr,
            &loci_pos,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
        )
        .unwrap();
        println!("windows_idx_head={:?}", windows_idx_head);
        println!("windows_idx_tail={:?}", windows_idx_tail);
        assert_eq!(windows_idx_head, vec![0, 1, 4, 7, 8]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
        assert_eq!(windows_idx_tail, vec![2, 3, 6, 7, 8]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
                                                           // Define sliding windows (redundant loci, i.e. per allele per locus)
        let loci_chr: Vec<String> = vec!["X", "X", "X", "Y", "Y"]
            .iter()
            .map(|&x| x.to_owned())
            .collect();
        let loci_pos: Vec<u64> = vec![123, 123, 123, 456, 456];
        let window_size_bp: u64 = 100;
        let window_slide_size_bp: u64 = 50;
        let min_loci_per_window: u64 = 1;
        let (windows_idx_head, windows_idx_tail) = define_sliding_windows(
            &loci_chr,
            &loci_pos,
            &window_size_bp,
            &window_slide_size_bp,
            &min_loci_per_window,
        )
        .unwrap();
        println!("windows_idx_head={:?}", windows_idx_head);
        println!("windows_idx_tail={:?}", windows_idx_tail);
        assert_eq!(windows_idx_head, vec![0, 3]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
        assert_eq!(windows_idx_tail, vec![2, 4]); // filtered out window start:2-end:3 which is a complete subset of window start:1-end:3
        let array1d: Array1<f64> = Array1::from_vec(vec![0.1, 0.2, 0.3, f64::NAN, 0.5]);
        assert_eq!(mean_array1_ignore_nan(&array1d.view()).unwrap(), 0.275);
        let mut array2d: Array2<f64> =
            Array2::from_shape_vec((2, 5), (0..10).map(|x| x as f64).collect::<Vec<f64>>())
                .unwrap();
        array2d[(0, 0)] = f64::NAN;
        println!("array2d={:?}", array2d);
        assert_eq!(
            Array1::from_shape_vec(5, vec![5.0, 3.5, 4.5, 5.5, 6.5]).unwrap(),
            mean_axis_ignore_nan(&array2d, 0).unwrap()
        );
        assert_eq!(
            Array1::from_shape_vec(2, vec![2.5, 7.0]).unwrap(),
            mean_axis_ignore_nan(&array2d, 1).unwrap()
        );
        let _fname = "./tests/test.csv".to_owned();
        let _delimiter = ",".to_owned();
        let _chr_col = 0;
        let _pos_start_col = 1;
        let _pos_end_col = 1;
        let _data_start_col = 2;
        let _data_end_col = 6;
        let (row_labels, column_labels, data) = load_table(
            &"./tests/test.csv".to_owned(),
            &",".to_owned(),
            &vec![0, 1],
            &2,
            &6,
        )
        .unwrap();
        println!("row_labels={:?}", row_labels);
        println!("column_labels={:?}", column_labels);
        println!("data={:?}", data);
        assert_eq!(row_labels.len(), 5);
        assert_eq!(column_labels.len(), 4);
        assert_eq!(data[0], vec![0.1, 83.2, 0.2, 0.0]);
        assert_eq!(data[4], vec![0.9, 12.0, 1.0, 0.0]);
    }
}
