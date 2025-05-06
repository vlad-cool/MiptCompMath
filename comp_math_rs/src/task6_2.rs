#![allow(dead_code)]

use std::io::Write;

use cauchy_problem::*;

mod adams_method;
mod algebraic_equation_solvers;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod matrix_operations;
mod nordsieck_method;
mod runge_kutta_method;
mod utils;

use crate::adams_method::AdamsMethod;
use crate::backward_differentiation_method::BackwardDifferentiationMethod;
use crate::nordsieck_method::NordsieckMethod;
use crate::nordsieck_method::NordsieckMethodType;
use crate::runge_kutta_method::RungeKuttaMethod;

fn f(_t: f64, x: &[f64; 3]) -> [f64; 3] {
    [
        77.27 * (x[1] + x[0] * (1.0 - 8.375e-6 * x[0] - x[1])),
        1.0 / 77.27 * (x[2] - (1.0 + x[0]) * x[1]),
        0.161 * (x[0] - x[2]),
    ]
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
        file.write(format!("{:.8}", solution.t[i]).as_bytes())
            .expect("failed to write to file");
        for j in 0..N {
            file.write(format!(", {:.8}", solution.x[i][j]).as_bytes())
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
        stop: 800.0,
        x_0: [0.5, 0.5, 0.5],
    };

    let print_progress: bool = false;

    let mut index: u32 = 0;
    let count_all: bool = std::env::args().any(|arg| arg == "--compute-all");

    if count_all {
        index += 1;
        let mut solver: RungeKuttaMethod<3, 3, 9> = RungeKuttaMethod::new(
            4,
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-1.0, 2.0, 0.0]],
            [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
            [0f64, 0.5f64, 1f64],
            "Kutta's third-order method (Explicit)".to_string(),
        );
        let tau: f64 = 0.00001;
        let save_every: u32 = (0.01f64 / tau).round() as u32;
        let start_time: std::time::Instant = std::time::Instant::now();
        let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
        println!("{:?}", res);
        write_csv(
            "Runge-Kutta".to_string(),
            format!("{}", index),
            solution,
            tau,
            start_time.elapsed(),
        );
    }

    // if count_all {
    //     index += 1;
    //     let mut solver: RungeKuttaMethod<3, 2, 6> = RungeKuttaMethod::new(
    //         2,
    //         [[0.0, 0.0], [0.5, 0.5]],
    //         [0.0f64, 1.0f64],
    //         [0.5f64, 0.5f64],
    //         "Crank-Nikolson method".to_string(),
    //     );
    //     let tau: f64 = 0.01;
    //     let save_every: u32 = (0.01f64 / tau).round() as u32;
    //     let start_time: std::time::Instant = std::time::Instant::now();
    //     let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
    //     println!("{:?}", res);
    //     write_csv(
    //         "Runge-Kutta".to_string(),
    //         format!("{}", index),
    //         solution,
    //         tau,
    //         start_time.elapsed(),
    //     );
    // }

    if count_all {
        index += 1;
        let mut solver: runge_kutta_method::RungeKuttaMethod<3, 4, 12> =
            runge_kutta_method::RungeKuttaMethod::new(
                3,
                [
                    [0.5, 0.0, 0.0, 0.0],
                    [1.0 / 6.0, 0.5, 0.0, 0.0],
                    [-0.5, 0.5, 0.5, 0.0],
                    [1.5, -1.5, 0.5, 0.5],
                ],
                [1.5, -1.5, 0.5, 0.5],
                [0.5, 2.0 / 3.0, 0.5, 1.0],
                "Four-stage, 3rd order, L-stable, Diagonally Implicit Runge-Kutta method"
                    .to_string(),
            );
        let tau: f64 = 0.01;
        let save_every: u32 = (0.01f64 / tau).round() as u32;
        let start_time: std::time::Instant = std::time::Instant::now();
        let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
        println!("{:?}", res);
        write_csv(
            "Runge-Kutta".to_string(),
            format!("{}", index),
            solution,
            tau,
            start_time.elapsed(),
        );
    }

    if count_all {
        for order in 1..5 {
            index += 1;
            let mut solver: AdamsMethod<3> = AdamsMethod::new(order, SolverType::Explicit);
            let tau: f64 = 0.000001;
            let save_every: u32 = (0.01f64 / tau).round() as u32;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Explicit Adams".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }
    }

    if count_all {
        for order in 1..5 {
            index += 1;
            let mut solver: AdamsMethod<3> = AdamsMethod::new(order, SolverType::Implicit);
            let tau: f64 = 0.01;
            let save_every: u32 = (0.01f64 / tau).round() as u32;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Implicit Adams".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }
    }

    if count_all {
        for order in 1..3 {
            index += 1;
            let mut solver: BackwardDifferentiationMethod<3> =
                BackwardDifferentiationMethod::new(order, SolverType::Explicit);
            let tau: f64 = 0.000001;
            let save_every: u32 = (0.01f64 / tau).round() as u32;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Explicit Backward Differentiation Method".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }
    }

    if count_all {
        for order in 1..4 {
            index += 1;
            let mut solver: BackwardDifferentiationMethod<3> =
                BackwardDifferentiationMethod::new(order, SolverType::Implicit);
            let tau: f64 = 0.01;
            let save_every: u32 = (0.01f64 / tau).round() as u32;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Implicit Backward Differentiation Method".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }
    }

    if count_all {
        for order in 1..5 {
            index += 1;
            let mut solver: NordsieckMethod<3> =
                NordsieckMethod::new(NordsieckMethodType::ImplicitBackwardDifferentiation(order));
            let tau: f64 = 0.01;
            let save_every: u32 = 1;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Nordsieck Method".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }

        for order in 1..5 {
            index += 1;
            let mut solver: NordsieckMethod<3> =
                NordsieckMethod::new(NordsieckMethodType::ImplicitAdams(order));
            let tau: f64 = 0.001;
            let save_every: u32 = 1;
            let start_time: std::time::Instant = std::time::Instant::now();
            let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
            println!("{:?}", res);
            write_csv(
                "Nordsieck Method".to_string(),
                format!("{}", index),
                solution,
                tau,
                start_time.elapsed(),
            );
        }
    }

    if !count_all {
        index += 1;
        let mut solver: NordsieckMethod<3> =
            NordsieckMethod::new(NordsieckMethodType::ImplicitBackwardDifferentiation(5));
        let tau: f64 = 0.001;
        let save_every: u32 = 1;
        let start_time: std::time::Instant = std::time::Instant::now();
        let (solution, res) = solver.solve(&mut problem, tau, print_progress, Some(save_every));
        println!("{:?}", res);
        write_csv(
            "Nordsieck Method".to_string(),
            format!("{}", index),
            solution,
            tau,
            start_time.elapsed(),
        );
    }
}
