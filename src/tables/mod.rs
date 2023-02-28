pub use self::{
    fisher_exact_test::fisher,
    chisq_test::chisq,
    cmh_test::cmh,
};

mod fisher_exact_test;
mod chisq_test;
mod cmh_test;
