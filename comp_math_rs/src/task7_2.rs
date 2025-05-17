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
mod phase_sl_method;
mod sturm_liuville_problem;
mod parameterized_function;

fn write_csv<const N: usize>(
    group: String,
    path: String,
    solution: CauchySolution<N>,
    tau: f64,
    self_value: f64,
    time: std::time::Duration,
) {
    let mut file =
        std::fs::File::create(format!("../task7_2_data/{path}.csv")).expect("Failed to open file");
    file.write(format!("{}\n", group).as_bytes())
        .expect("failed to wrtite to file");
    file.write(format!("{}\n", tau).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", self_value).as_bytes())
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

    for i in 5..30 {
        let i = i * 2;
        let start_time: std::time::Instant = std::time::Instant::now();
        let equation = |v: &[f64; 1]| {
            let v: f64 = v[0];
            let a: f64 = v;
            // let a: f64 = a.exp() + 0.0;

            let mut cauchy_problem: CauchyProblem<2, _> = CauchyProblem {
                f: &mut |_: f64, x: &[f64; 2]| [a * x[1], (1.0 - x[0].exp()) / a],
                start: 0.0,
                stop: 120.0,
                x_0: [0.0, 1.0],
            };
            let (solution, res) = cauchy_solver.solve(&mut cauchy_problem, 0.01, false, None);
            res.expect("Faield to solve differential equation");
            [(solution.x.last().unwrap()[0].powi(2) * 10.0
                + (solution.x.last().unwrap()[1] - 1.0).powi(2) * 1.0)]
        };

        let res: Result<[f64; 1], &'static str> = solve_newton(equation, &[i as f64], None);
        
        let a = match res {
            Err(_) => continue,
            Ok(val) => val,
        };
        
        let a: f64 = a[0];
        // let a: f64 = a.exp() + 10.0;

        // println!("a: {a}");

        let mut cauchy_problem: CauchyProblem<2, _> = CauchyProblem {
            f: &mut |_: f64, x: &[f64; 2]| [a * x[1], (1.0 - x[0].exp()) / a],
            start: 0.0,
            stop: 120.0,
            x_0: [0.0, 1.0],
        };

        let (solution, _res) = cauchy_solver.solve(&mut cauchy_problem, 0.001, false, None);
        // println!("{:?}", solution.x.last().unwrap());

        write_csv(
            format!("shooting_method"),
            format!("{i}"),
            solution,
            0.001,
            a,
            start_time.elapsed(),
        );
    }
}
