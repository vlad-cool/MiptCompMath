pub struct CauchyProblem<const N: usize> {
    pub f: fn(t: f64, &[f64; N]) -> [f64; N],
    pub start: f64,
    pub stop: f64,
    pub x_0: [f64; N],
}

pub struct CauchySolution<const N: usize> {
    pub t: std::vec::Vec<f64>,
    pub x: std::vec::Vec<[f64; N]>,
    pub method_name: String,
}

pub trait CauchySolver<const N: usize> {
    fn solve(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        print_progress: bool,
        save_every: Option<u32>,
    ) -> (CauchySolution<N>, Result<(), &'static str>);
    fn get_name(&self) -> String;
}

pub enum SolverType {
    Explicit,
    Implicit,
}

impl std::fmt::Display for SolverType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SolverType::Explicit => write!(f, "Explicit"),
            SolverType::Implicit => write!(f, "Implicit"),
        }
    }
}
