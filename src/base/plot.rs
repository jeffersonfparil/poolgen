use ndarray::prelude::*;
use plotlib::page::Page;
use plotlib::repr::Plot;
use plotlib::style::{PointMarker, PointStyle};
use plotlib::view::ContinuousView;
use std::io::{self, Error, ErrorKind};

pub fn plot_scatter_2d(x: &Array1<f64>, y: &Array1<f64>, fname_svg: String) -> io::Result<String> {
    let n = x.len();
    if n != y.len() {
        return Err(Error::new(
            ErrorKind::Other,
            "The two vectors are not the same lengths.",
        ));
    }
    let data: Vec<(f64, f64)> = x.iter().zip(y).map(|(&x, &y)| (x, y)).collect();
    let plot = Plot::new(data).point_style(PointStyle::new().colour("#35C788"));
    let view = ContinuousView::new().add(plot).x_label("x").y_label("y");
    Page::single(&view).save(&fname_svg[..]).unwrap();
    Ok(fname_svg)
}
