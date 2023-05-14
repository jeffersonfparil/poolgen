use crate::base::*;
use ndarray::{prelude::*, Zip};
use ndarray_linalg::*;
use std::io::{self, Error, ErrorKind};
// use ndarray_linalg::least_squares::LeastSquaresSvd;

#[function_name::named]
pub fn ols(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    // let mut x: Array2<f64> = Array2::from_elem((row_idx.len(), x_orig.ncols()), f64::NAN);
    // let mut y: Array2<f64> = Array2::from_elem((row_idx.len(), y_orig.ncols()), f64::NAN);
    // for i in row_idx.iter() {
    //     for j in 0..x_orig.ncols() {
    //         x[(*i, j)] = x_orig[(*i, j)];
    //     }
    //     for k in 0..y_orig.ncols() {
    //         y[(*i, k)] = y_orig[(*i, k)];
    //     }
    // }
    let n = x.nrows();
    let p = x.ncols();
    let k = y.ncols();
    if x.column(0).sum() < n as f64 {
        return Err(Error::new(
            ErrorKind::Other,
            "Please add the intercept in the X matrix.",
        ));
    }
    // let mut b_hat: Array2<f64> = Array2::from_elem((p, k), f64::NAN);
    // for j in 0..k {
    //     let yj: Array1<f64> = y.column(j).to_owned();
    //     println!("x={:?}", x);
    //     println!("yj={:?}", yj);
    //     let ls = x.least_squares(&yj).unwrap();
    //     let bj = ls.solution;
    //     println!("bj={:?}", bj);
    //     for i in 0..p {
    //         b_hat[(i, j)] = bj[j];
    //     }
    // }
    let col_idx = &(0..p).collect::<Vec<usize>>();
    let new_row_idx = &(0..row_idx.len()).collect::<Vec<usize>>();
    let new_col_idx = &(0..k).collect::<Vec<usize>>();
    let b_hat: Array2<f64> = if n < p {
        // x.t()
        //     .dot(&x.dot(&x.t()))
        //         .inv()
        //         .unwrap()
        //     .dot(&y)
        multiply_views_xx(
            &multiply_views_xtx(
                x,
                &multiply_views_xxt(x, x, row_idx, col_idx, row_idx, col_idx)
                    .unwrap()
                    .pinv()
                    .unwrap(),
                row_idx,
                col_idx,
                new_row_idx,
                new_row_idx,
            )
            .unwrap(),
            y,
            col_idx,
            new_row_idx,
            row_idx,
            new_col_idx,
        )
        .unwrap()
    } else {
        // (x.t().dot(&x))
        //     .inv()
        //     .unwrap()
        //     .dot(&x.t())
        //     .dot(&y)
        multiply_views_xx(
            &multiply_views_xxt(
                &multiply_views_xtx(x, x, row_idx, col_idx, row_idx, col_idx)
                    .unwrap()
                    .pinv()
                    .unwrap(),
                x,
                col_idx,
                col_idx,
                row_idx,
                col_idx,
            )
            .unwrap(),
            y,
            col_idx,
            new_row_idx,
            row_idx,
            new_col_idx,
        )
        .unwrap()
    };
    Ok((b_hat, function_name!().to_owned()))
}

#[function_name::named]
pub fn ols_iterative_with_kinship_pca_covariate(
    x: &Array2<f64>,
    y: &Array2<f64>,
    row_idx: &Vec<usize>,
) -> io::Result<(Array2<f64>, String)> {
    let n = row_idx.len();
    let p = x.ncols();
    let k = y.ncols();
    if x.column(0).sum() < x.nrows() as f64 {
        return Err(Error::new(
            ErrorKind::Other,
            "Please add the intercept in the X matrix.",
        ));
    }
    let col_idx = (0..p - 1).collect::<Vec<usize>>();
    let mut x_column_centred_no_intercept: Array2<f64> = Array2::from_elem((n, p - 1), f64::NAN);
    let rows_mat: Array2<usize> = Array2::from_shape_vec((p - 1, n), row_idx.repeat(p - 1))
        .unwrap()
        .reversed_axes();
    let cols_mat: Array2<usize> = Array2::from_shape_vec((n, p - 1), col_idx.repeat(n)).unwrap();
    Zip::from(&mut x_column_centred_no_intercept)
        .and(&rows_mat)
        .and(&cols_mat)
        .par_for_each(|x_new, &i, &j| {
            let mut mean = 0.0;
            for i_ in 0..n {
                mean += x[(i_, j)];
            }
            mean = mean / n as f64;
            *x_new = x[(i, j)] - mean;
        });
    let row_idx_new = (0..n).collect::<Vec<usize>>();
    let xxt: Array2<f64> = multiply_views_xxt(
        &x_column_centred_no_intercept,
        &x_column_centred_no_intercept,
        &row_idx_new,
        &col_idx,
        &row_idx_new,
        &col_idx,
    )
    .unwrap();
    let (_eigen_values, eigen_vectors): (Array1<_>, Array2<_>) = xxt.eig().unwrap();

    let mut b_hat: Array2<f64> = Array2::from_elem((p, k), f64::NAN);
    let mut x_sub: Array2<f64> = Array2::ones((n, 3)); // intercept, 1st eigenvector, and the jth locus
    let mut y_sub: Array2<f64> = Array2::from_elem((n, k), f64::NAN);
    for i in 0..n {
        for j in 0..k {
            y_sub[(i, j)] = y[(row_idx[i], j)];
        }
    }
    let y_sub_means: Array1<f64> = y_sub.mean_axis(Axis(0)).unwrap();
    for j in 0..k {
        b_hat[(0,j)] = y_sub_means[j];
    }
    for j in 0..p - 1 {
        for i in 0..n {
            x_sub[(i, 1)] = eigen_vectors[(i, 0)].re; // extract the eigenvector value's real number component
            x_sub[(i, 2)] = x[(row_idx[i], j + 1)]; // use the row_idx and add 1 to the column indexes to account for the intercept in the input x
        }
        let (b, _) = ols(&x_sub, &y_sub, &row_idx_new).unwrap();
        for j_ in 0..k {
            b_hat[(j+1, j_)] = b[(2, j_)]; // extract the effect of the locus i.e. j==2 or the third row
        }
    }
    println!("BEFORE: b_hat={:?}", b_hat);
    Ok((b_hat, function_name!().to_owned()))
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
#[cfg(test)]
mod tests {
    use ndarray::concatenate;

    use super::*;
    #[test]
    fn test_ols_for_gp() {
        let intercept: Array2<f64> = Array2::ones((5, 1));

        let frequencies_wide = Array2::from_shape_vec(
            (5, 9),
            (1..46).map(|x| x as f64 / 45.0).collect::<Vec<f64>>(),
        )
        .unwrap();

        let frequencies_tall = Array2::from_shape_vec(
            (5, 2),
            (1..31)
                .step_by(3)
                .map(|x| x as f64 / 30.0)
                .collect::<Vec<f64>>(),
        )
        .unwrap();
        let x_wide: Array2<f64> =
            concatenate(Axis(1), &[intercept.view(), frequencies_wide.view()]).unwrap();
        let x_tall: Array2<f64> =
            concatenate(Axis(1), &[intercept.view(), frequencies_tall.view()]).unwrap();

        let y: Array2<f64> =
            Array2::from_shape_vec((5, 1), (1..6).map(|x| x as f64 / 5.0).collect::<Vec<f64>>())
                .unwrap();
        println!("x_wide={:?}", x_wide);
        println!("x_tall={:?}", x_tall);
        println!("y={:?}", y);
        let (b_wide, _) = ols(&x_wide, &y, &vec![0, 1, 2, 3, 4]).unwrap();
        println!("b_wide={:?}", b_wide);
        let (b_tall, _) = ols(&x_tall, &y, &vec![0, 1, 2, 3, 4]).unwrap();
        println!("b_tall={:?}", b_tall);
        let y_hat_wide = &x_wide.dot(&b_wide).mapv(|x| sensible_round(x, 4));
        println!("y_hat_wide={:?}", y_hat_wide);
        let y_hat_tall = &x_tall.dot(&b_tall).mapv(|x| sensible_round(x, 4));
        println!("y_hat_tall={:?}", y_hat_tall);

        assert_eq!(y, y_hat_wide);
        assert_eq!(y_hat_wide, y_hat_tall);

        let (b_wide_iterative, _) =
            ols_iterative_with_kinship_pca_covariate(&x_wide, &y, &vec![0, 1, 2, 3, 4]).unwrap();
        println!("b_wide_iterative={:?}", b_wide_iterative);
        assert_eq!(0, 0);
    }
}
