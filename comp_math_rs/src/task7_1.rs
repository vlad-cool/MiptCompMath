#![allow(dead_code)]

use std::io::Write;

use boundary_problem::BoundarySolution;
use boundary_problem::NonlinearBoundaryProblem;
use cauchy_problem::CauchySolver;

mod adams_method;
mod algebraic_equation_solvers;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod newton_boundary_method;
mod runge_kutta_method;
mod shooting_method;
mod tridiagonal_method;
mod utils;

fn write_csv(path: String, solution: BoundarySolution, step: f64, time: std::time::Duration) {
    let mut file = std::fs::File::create(format!("../task7_1_data/{}.csv", path))
        .expect("Failed to open file");
    file.write(format!("{}\n", solution.method_name).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", step).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    for i in 0..solution.x.len() {
        file.write(format!("{:.8}", solution.x[i]).as_bytes())
            .expect("failed to write to file");
        file.write(format!(", {:.8}", solution.y[i]).as_bytes())
            .expect("failed to write to file");
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn main() {
    let mut problem: NonlinearBoundaryProblem<_> = NonlinearBoundaryProblem {
        // f(x, [y, y'])
        f: |_: f64, x: &[f64; 2]| [x[1], -x[1].powi(2) / (2.0 - x[0])],
        x_0: 0.0,
        x_n: 1.0,
        a_1: 0.0,
        b_1: 1.0,
        u_1: 1.9,
        a_2: 0.0,
        b_2: 1.0,
        u_2: 0.0,
    };

    let solver: backward_differentiation_method::BackwardDifferentiationMethod<2> =
        backward_differentiation_method::BackwardDifferentiationMethod::new(
            1,
            cauchy_problem::SolverType::Explicit,
        );
    let mut solver: Box<dyn CauchySolver<2, _>> = Box::new(solver);
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, _res) =
        crate::shooting_method::shooting_method(&mut problem, &mut solver, 0.001);
    let duration: std::time::Duration = start_time.elapsed();
    write_csv(format!("shooting_method"), solution, 0.001, duration);

    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, _res) =
        crate::newton_boundary_method::newton_boundary_method(&mut problem, 0.0001, 1e-5);
    // let (solution, _res) = crate::newton_boundary_method::newton_boundary_method(&mut problem, 0.001, 1e-2);
    let duration: std::time::Duration = start_time.elapsed();
    write_csv(format!("newton_method"), solution, 0.001, duration);
}
