#![allow(dead_code)]

use std::io::Write;

use algebraic_equation_solvers::solve_newton;
use cauchy_problem::CauchyProblem;
use cauchy_problem::CauchySolution;
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

fn write_csv<const N: usize>(
    group: String,
    solution: CauchySolution<N>,
    tau: f64,
    time: std::time::Duration,
) {
    let mut file =
        std::fs::File::create(format!("../comp_math_rs/test.csv")).expect("Failed to open file");
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
        file.write(format!(", {:.8}", solution.x[i][0]).as_bytes())
            .expect("failed to write to file");
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn main() {
    let mut cauchy_solver: runge_kutta_method::RungeKuttaMethod<2, 3, 6> =
        runge_kutta_method::RungeKuttaMethod::new(
            4,
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-1.0, 2.0, 0.0]],
            [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
            [0f64, 0.5f64, 1f64],
            "Kutta's third-order method (Explicit)".to_string(),
        );

    let equation = |v: &[f64; 1]| {
        let v: f64 = v[0];
        let a: f64 = v;
        let a: f64 = v.exp() + 10.0;

        let mut cauchy_problem: CauchyProblem<2, _> = CauchyProblem {
            f: &mut |_: f64, x: &[f64; 2]| [a * x[1], (1.0 - x[0].exp()) / a],
            start: 0.0,
            stop: 120.0,
            x_0: [0.0, 1.0],
        };
        let (solution, res) = cauchy_solver.solve(&mut cauchy_problem, 0.001, false, None);
        res.expect("Faield to solve differential equation");
        [solution.x.last().unwrap()[0].powi(2) * 1.0
            + (solution.x.last().unwrap()[1] - 1.0).powi(2) * 1.0]
    };

    // let solver: backward_differentiation_method::BackwardDifferentiationMethod<2> =
    //     backward_differentiation_method::BackwardDifferentiationMethod::new(
    //         1,
    //         cauchy_problem::SolverType::Explicit,
    //     );
    // // let mut solver: Box<dyn CauchySolver<2, _>> = Box::new(solver);
    // let start_time: std::time::Instant = std::time::Instant::now();
    // let (solution, _res) =
    //     crate::shooting_method::shooting_method(&mut problem, &mut solver, 0.001);
    // let duration: std::time::Duration = start_time.elapsed();
    let res = solve_newton(equation, &[0.0], None);
    let a: f64 = res.unwrap()[0];
    let a: f64 = a.exp() + 10.0;

    println!("a: {a}");

    let mut cauchy_problem: CauchyProblem<2, _> = CauchyProblem {
        f: &mut |_: f64, x: &[f64; 2]| [x[1], 1.0 - x[0].exp()],
        start: 0.0,
        stop: 120.0,
        x_0: [0.0, a],
    };

    let (solution, _res) = cauchy_solver.solve(&mut cauchy_problem, 0.001, false, None);
    println!("{:?}", solution.x.last().unwrap());

    write_csv(
        format!("shooting_method"),
        solution,
        0.001,
        std::time::Duration::from_micros(0),
    );
}
