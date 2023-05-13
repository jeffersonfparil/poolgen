use crate::base::*;
use ndarray::prelude::*;
use std::io::{self, Error, ErrorKind};

use ndarray_linalg::least_squares::LeastSquaresSvd;

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
    }
}
