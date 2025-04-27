// y'' + q(x)y' + p(x)y = f(x)
pub struct LinearBoundaryProblem {
    pub q: fn(x: f64) -> f64,
    pub p: fn(x: f64) -> f64,
    pub f: fn(x: f64) -> f64,
    pub x_0: f64,
    pub x_n: f64,
    pub a_1: f64,
    pub b_1: f64,
    pub u_1: f64,
    pub a_2: f64,
    pub b_2: f64,
    pub u_2: f64,
}

pub struct LinearBoundarySolution {
    pub x: std::vec::Vec<f64>,
    pub y: std::vec::Vec<f64>,
    pub method_name: String,
}
pub struct NonlinearBoundaryProblem {
    pub f: fn(x: f64, y: &[f64; 2]) -> [f64; 2],
    pub x_0: f64,
    pub y_0: f64,
    pub x_n: f64,
    pub y_n: f64,
}

pub struct NonlinearBoundarySolution {
    pub x: std::vec::Vec<f64>,
    pub y: std::vec::Vec<f64>,
    pub method_name: String,
}
