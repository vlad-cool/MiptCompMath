pub struct BoundaryProblem {
    pub f: fn(x: f64, y: &[f64; 2]) -> [f64; 2],
    pub x_0: f64,
    pub y_0: f64,
    pub x_n: f64,
    pub y_n: f64,
}

pub struct BoundarySolution {
    pub x: std::vec::Vec<f64>,
    pub y: std::vec::Vec<f64>,
    pub method_name: String,
}
