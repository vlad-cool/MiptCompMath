use crate::parameterized_function::ParameterizedFunction;

pub enum SturmLiuvilleBoundary {
    FixedValDer{
        x_f: f64,
        dx_f: f64,
        x_l: f64,
        dx_l: f64,
    },
}

pub struct SturmLiuvilleProblem {
    pub f: Box<dyn ParameterizedFunction<1, 1>>,
    pub boundary: SturmLiuvilleBoundary,
}

pub struct SturmLiuvilleSolution {
    pub lambda: f64,
    pub solution: std::vec::Vec::<[f64; 2]>,
}

pub trait SturmLiuvilleSolver {
    fn solve(&self, problem: &mut SturmLiuvilleProblem) -> SturmLiuvilleSolution;
}
