#![allow(dead_code)]

use std::io::Write;

use boundary_problem::NonlinearBoundaryProblem;
use boundary_problem::NonlinearBoundarySolution;
use cauchy_problem::CauchySolver;

mod adams_method;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod algebraic_equation_solvers;
mod runge_kutta_method;
mod shooting_method;
mod tridiagonal_method;
mod utils;

fn write_csv(path: String, solution: NonlinearBoundarySolution, step: f64, time: std::time::Duration) {
    let mut file = std::fs::File::create(format!("../task7_1_data/{}.csv", path))
        .expect("Failed to open file");
    file.write(format!("{}\n", solution.method_name).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", step).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    for i in 0..solution.x.len() {
        file.write(format!("{}", solution.x[i]).as_bytes())
            .expect("failed to write to file");
        file.write(format!(", {}", solution.y[i]).as_bytes())
            .expect("failed to write to file");
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn main() {
    let problem: NonlinearBoundaryProblem = NonlinearBoundaryProblem {
        f: |_: f64, x: &[f64; 2]| [-x[0].powi(2) / (2.0 - x[1]), x[0]],
        x_0: 0.0,
        y_0: 1.9,
        x_n: 1.0,
        y_n: 0.0,
    };

    let solver: backward_differentiation_method::BackwardDifferentiationMethod<2> =
        backward_differentiation_method::BackwardDifferentiationMethod::new(
            1,
            cauchy_problem::SolverType::Explicit,
        );
    let mut solver: Box<dyn CauchySolver<2>> = Box::new(solver);
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, _res) = crate::shooting_method::shooting_method(&problem, &mut solver, 0.001);
    let duration: std::time::Duration = start_time.elapsed();
    write_csv(format!("shooting_method"), solution, 0.001, duration);
}
