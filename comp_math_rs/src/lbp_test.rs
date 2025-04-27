#![allow(dead_code)]

use std::io::Write;

use boundary_problem::LinearBoundaryProblem;
use boundary_problem::BoundarySolution;

mod adams_method;
mod algebraic_equation_solvers;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod runge_kutta_method;
mod shooting_method;
mod tridiagonal_method;
mod utils;

fn write_csv(solution: BoundarySolution, step: f64, time: std::time::Duration) {
    let mut file = std::fs::File::create(format!("test.csv"))
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
    let problem: LinearBoundaryProblem = LinearBoundaryProblem {
        q: |_: f64| -1.0,
        p: |_: f64| 0.0,
        f: |_: f64| 0.0,
        x_0: 0.0,
        x_n: 5.0,
        a_1: 0.0,
        b_1: 1.0,
        u_1: 1.0,
        a_2: 0.0,
        b_2: 1.0,
        u_2: 10.0,
    };

    let start_time: std::time::Instant = std::time::Instant::now();
    let solution: BoundarySolution = crate::tridiagonal_method::tridiagonal_method(&problem, 0.001);
    let duration: std::time::Duration = start_time.elapsed();
    write_csv(solution, 0.001, duration);
}
