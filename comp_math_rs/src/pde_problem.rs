use crate::parameterized_function::ParameterizedFunction;

// du / dt + a du / dx = f(t, x)
pub struct PartialDerivativeEquation {
    pub psi: Box<dyn ParameterizedFunction<1, 1>>, // u(t, 0) = psi(t)
    pub phi: Box<dyn ParameterizedFunction<1, 1>>, // u(0, x) = phi(x)
    pub f: Box<dyn ParameterizedFunction<2, 1>>,
    pub a: f64,
    pub x_f: f64,
    pub t_f: f64,
    pub x_l: f64,
    pub t_l: f64,
}

pub struct PartialDerivativeSolution {
    pub solution: std::vec::Vec<std::vec::Vec<f64>>,
}

pub trait PartialDerivativeEquationSolver {
    fn solve(
        &self,
        problem: &mut PartialDerivativeEquation,
        tau: f64,
        h: f64,
    ) -> PartialDerivativeSolution;
}
