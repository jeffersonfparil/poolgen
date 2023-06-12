use ndarray::prelude::*;
use plotlib::page::Page;
use plotlib::repr::{Histogram, HistogramBins, Plot};
use plotlib::style::PointStyle;
use plotlib::view::ContinuousView;
use std::io::{self, Error, ErrorKind};

// Functions for tests and debugging mostly. This module will probably not be expanded/developed further.

pub fn plot_scatter_2d(
    x: &Array1<f64>,
    y: &Array1<f64>,
    xlab: &str,
    ylab: &str,
    fname_svg: &str,
) -> io::Result<String> {
    let n = x.len();
    if n != y.len() {
        return Err(Error::new(
            ErrorKind::Other,
            "The two vectors are not the same lengths.",
        ));
    }
    let data: Vec<(f64, f64)> = x.iter().zip(y).map(|(&x, &y)| (x, y)).collect();
    let plot = Plot::new(data).point_style(PointStyle::new().colour("#35C788"));
    let view = ContinuousView::new().add(plot).x_label(xlab).y_label(ylab);
    Page::single(&view).save(fname_svg).unwrap();
    Ok(fname_svg.to_owned())
}

pub fn plot_histogram_1d(
    x: &Array1<f64>,
    nbins: usize,
    lab: &str,
    fname_svg: &str,
) -> io::Result<String> {
    let h = Histogram::from_slice(
        &x.iter().map(|x| *x).collect::<Vec<f64>>(),
        HistogramBins::Count(nbins),
    );
    let view = ContinuousView::new()
        .add(h)
        .x_label(lab)
        .y_label("Frequency");
    Page::single(&view).save(fname_svg).unwrap();
    Ok(fname_svg.to_owned())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use statrs;
    #[test]
    fn test_plot() {
        let n = 100;
        let mut rng = rand::thread_rng();
        let d = statrs::distribution::Normal::new(0.0, 1.0).unwrap();
        let x: Array1<f64> =
            Array1::from_shape_vec(n, d.sample_iter(&mut rng).take(n).collect::<Vec<f64>>())
                .unwrap();
        let y = 2.0 * &x;
        let fname_scatter_plot =
            plot_scatter_2d(&x, &y, "xlab", "ylab", "test-scatterplot.svg").unwrap();
        let fname_histogram =
            plot_histogram_1d(&x, 17, "Some normal variable", "test-histogram.svg").unwrap();
        assert_eq!(fname_scatter_plot, "test-scatterplot.svg".to_owned());
        assert_eq!(fname_histogram, "test-histogram.svg".to_owned());
    }
}
