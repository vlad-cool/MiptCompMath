#![allow(dead_code)]

use std::io::Write;
use crate::cauchy_problem::*;
use crate::nordsieck_method::*;
use crate::adams_method::*;

mod adams_method;
mod algebraic_equation_solvers;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod runge_kutta_method;
mod nordsieck_method;
mod tridiagonal_method;
mod matrix_operations;
mod utils;

fn f(_t: f64, x: &[f64; 3]) -> [f64; 3] {
    // [x[0].powi(1)]
    // [x[0].powi(2)]
    [x[1], -x[0], 0.0]
}

fn write_csv<const N: usize>(
    group: String,
    path: String,
    solution: CauchySolution<N>,
    tau: f64,
    time: std::time::Duration,
) {
    let mut file = std::fs::File::create(format!("../task6_2_data/{}.csv", path))
        .expect("Failed to open file");
    file.write(format!("{}\n", group).as_bytes())
        .expect("failed to wrtite to file");
    file.write(format!("{}\n", solution.method_name).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", tau).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    for i in 0..solution.t.len() {
        file.write(format!("{}", solution.t[i]).as_bytes())
            .expect("failed to write to file");
        for j in 0..N {
            file.write(format!(", {}", solution.x[i][j]).as_bytes())
                .expect("failed to write to file");
        }
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn main() {
    let mut problem: CauchyProblem<3, _> = CauchyProblem {
        f: &mut f,
        start: 0.0,
        stop: 10.0,
        x_0: [0.0, 1.0, 0.0],
    };

    // let print_progress = false;

    let mut solver: AdamsMethod<3> =
        AdamsMethod::new(3, SolverType::Explicit);
    let mut solver: NordsieckMethod<3> =
        NordsieckMethod::new(NordsieckMethodType::ImplicitBackwardDifferentiation(4));

    let tau: f64 = 0.01;
    // let tau: f64 = 1.0;
    let save_every: u32 = 1;
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, res) = solver.solve(&mut problem, tau, true, Some(save_every));
    println!("{:?}", res);
    write_csv(
        "Nordsieck Method".to_string(),
        format!("nordsieck_test"),
        solution,
        tau,
        start_time.elapsed(),
    );
}
