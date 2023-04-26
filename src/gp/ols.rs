use nalgebra::{self, DMatrix, DVector};
use nalgebra_lapack::{self, SVD};
use std::io::{self, Error, ErrorKind};
use std::sync::{Arc, Mutex};

#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct OlsBeta {
    j: usize,
    b: f64
}

fn ols_per_chunk(xjk: &DMatrix<f64>, y: &DVector<f64>, j: usize) -> io::Result<Vec<OlsBeta>> {
    let p = xjk.ncols();
    let mut out = vec![];
    for j_ in 0..p {
        let x = xjk.column(j_);
        if x.column(0).variance() < f64::EPSILON {
            out.push(OlsBeta{j:j+j_, b:0.0});
        }
        let x_with_intercept = x.clone().insert_column(0, 1.0);
        let ols_beta: OlsBeta = match (&x_with_intercept.transpose() * &x_with_intercept).try_inverse() {
            Some(inv_xtx) => OlsBeta{j:j+j_, b:(inv_xtx * &x_with_intercept.transpose() * y)[1]},
            None => OlsBeta{j:j+j_, b:0.0}
        };
        out.push(ols_beta);
    }
    Ok(out)
}

#[function_name::named]
pub fn ols_new(x: &DMatrix<f64>, y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)> {
    let (n, p) = x.shape();
    let (n_, _) = y.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    if x.column(0).sum() < n as f64 {
        return Err(Error::new(
            ErrorKind::Other,
            "Please add the intercept in the X matrix.",
        ));
    }
    
        let mut thread_objects = Vec::new();
        // Vector holding all returns from pileup2sync_chunk()
        let thread_ouputs: Arc<Mutex<Vec<OlsBeta>>> = Arc::new(Mutex::new(Vec::new())); // Mutated within each thread worker
                                                                                       // Making four separate threads calling the `search_for_word` function
        let max_threads = match p > 30_000 {
            true => 30_000,
            false => p
        };
        let mut chunk_size = ((p - 1) as f64 / max_threads as f64).ceil() as usize;
        for j_ in 0..(max_threads-1) {
            let j = (j_ * chunk_size) + 1;
            if (j + chunk_size) > p {
                chunk_size = p - j;
            }
            // Clone pileup2sync_chunk parameters
            let thread_ouputs_clone = thread_ouputs.clone(); // Mutated within the current thread worker
            let xjk = x.columns(j, chunk_size).clone_owned();
            let y = y.column(0).clone_owned();
            let thread = std::thread::spawn(move || {
                let ols_beta = ols_per_chunk(&xjk, &y, j).unwrap();
                for b in ols_beta.into_iter() {
                    thread_ouputs_clone
                        .lock()
                        .unwrap()
                        .push(b);
                }
            });
            thread_objects.push(thread);
        }
        // Waiting for all threads to finish
        for thread in thread_objects {
            let _ = thread.join().expect("Unknown thread error occured.");
        }
        let mut vec_b: Vec<OlsBeta> = Vec::new();
        for b in thread_ouputs.lock().unwrap().iter() {
            vec_b.push(b.to_owned());
        }
        vec_b.sort_by(|a, b| a.j.partial_cmp(&b.j).unwrap());
        let mut b_hat = DMatrix::from_element(p, 1, f64::NAN);
        b_hat[(0,0)] = y.mean();
        for j in 0..(p-1) {
            b_hat[(j+1,0)] = vec_b[j].b / (p-1) as f64;
        }
        

    Ok((b_hat, function_name!().to_owned()))
}


#[function_name::named]
pub fn ols(x: &DMatrix<f64>, y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)> {
    let (n, p) = x.shape();
    let (n_, _) = y.shape();
    if n != n_ {
        return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
    }
    if x.column(0).sum() < n as f64 {
        return Err(Error::new(
            ErrorKind::Other,
            "Please add the intercept in the X matrix.",
        ));
    }
    let lambda = 1.0;
    let b_hat: DMatrix<f64> = if n < p {
        x.transpose() * ((x * x.transpose()).add_scalar(lambda)).try_inverse().unwrap() * y
    } else {
        ((x.transpose() * x).add_scalar(lambda)).try_inverse().unwrap() * x.transpose() * y
    };
    // println!("##############################");
    // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
    Ok((b_hat, function_name!().to_owned()))
}


// #[function_name::named]
// pub fn ols2(x: &DMatrix<f64>, y: &DMatrix<f64>) -> io::Result<(DMatrix<f64>, String)> {
//     let (n, p) = x.shape();
//     let (n_, _) = y.shape();
//     if n != n_ {
//         return Err(Error::new(ErrorKind::Other, "The number of samples in the dependent and independent variables are not the same size."));
//     }
//     if x.column(0).sum() < n as f64 {
//         return Err(Error::new(
//             ErrorKind::Other,
//             "Please add the intercept in the X matrix.",
//         ));
//     }
//     let lambda = 0.5;
//     let b_hat: DMatrix<f64> = if n < p {
//         x.transpose() * ((x * x.transpose()).add_scalar(lambda)).try_inverse().unwrap() * y
//     } else {
//         ((x.transpose() * x).add_scalar(lambda)).try_inverse().unwrap() * x.transpose() * y
//     };
//     // println!("##############################");
//     // println!("{:?}: {:?}", function_name!().to_owned(), b_hat);
//     Ok((b_hat, function_name!().to_owned()))
// }
