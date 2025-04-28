// y'' + q(x)y' + p(x)y = f(x)
// a_1 y'(x_0) + b_1 y(x_0) = u_1
// a_2 y'(x_n) + b_2 y(x_n) = u_2
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

// y'' = f(x, [y, y'])
// a_1 y'(x_0) + b_1 y(x_0) = u_1
// a_2 y'(x_n) + b_2 y(x_n) = u_2
pub struct NonlinearBoundaryProblem<F>
where
    F: FnMut(f64, &[f64; 2]) -> [f64; 2],
{
    pub f: F,
    pub x_0: f64,
    pub x_n: f64,
    pub a_1: f64,
    pub b_1: f64,
    pub u_1: f64,
    pub a_2: f64,
    pub b_2: f64,
    pub u_2: f64,
}

pub struct BoundarySolution {
    pub x: std::vec::Vec<f64>,
    pub y: std::vec::Vec<f64>,
    pub method_name: String,
}
