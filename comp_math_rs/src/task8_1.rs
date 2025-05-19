#![allow(dead_code)]

use std::io::Write;

mod adams_method;
mod algebraic_equation_solvers;
mod backward_differentiation_method;
mod boundary_problem;
mod cauchy_problem;
mod newton_boundary_method;
mod parameterized_function;
mod pde_problem;
mod pde_solvers;
mod phase_sl_method;
mod runge_kutta_method;
mod shooting_method;
mod sturm_liuville_problem;
mod tridiagonal_method;
mod utils;

use parameterized_function::ParameterizedFunction;
use pde_problem::*;
use pde_solvers::*;

fn write_csv(
    path: String,
    solution: crate::pde_problem::PartialDerivativeSolution,
    tau: f64,
    h: f64,
    error: f64,
    time: std::time::Duration,
) {
    let x_save_step: usize = std::cmp::max(1, (0.25 * 0.01 / h) as usize);
    let t_save_step: usize = std::cmp::max(1, (0.25 * 0.01 / tau) as usize);
    let mut file =
        std::fs::File::create(format!("../task8_1_data/{path}.csv")).expect("Failed to open file");
    file.write(format!("{}\n", h).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", tau).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", error).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    for t_i in (0..solution.solution.len()).step_by(t_save_step) {
        for x_i in (0..solution.solution[t_i].len()).step_by(x_save_step) {
            let mut avg: f64 = 0.0;
            let mut n: usize = 0;
            for t_j in 0..t_save_step {
                for x_j in 0..x_save_step {
                    if t_i + t_j < solution.solution.len()
                        && x_i + x_j < solution.solution[t_i + t_j].len()
                    {
                        avg += solution.solution[t_i + t_j][x_i + x_j];
                        n += 1;
                    }
                }
            }
            avg /= n as f64;
            file.write(format!("{:.8}, ", avg).as_bytes())
                .expect("failed to write to file");
        }
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

struct F {}

impl ParameterizedFunction<2, 1> for F {
    fn calc(&mut self, x: &[f64; 2]) -> [f64; 1] {
        [x[1].powi(2) + 4.0 * x[0] * x[1]]
        // [0.0]
    }

    fn get_parameter(
        &self,
        _id: usize,
    ) -> Result<f64, parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }

    fn set_parameter(
        &mut self,
        _id: usize,
        _parameter: f64,
    ) -> Result<(), parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }
}

struct Phi {}

impl ParameterizedFunction<1, 1> for Phi {
    fn calc(&mut self, x: &[f64; 1]) -> [f64; 1] {
        [(x[0] * 10.0 * std::f64::consts::PI).sin()]
    }

    fn get_parameter(
        &self,
        _id: usize,
    ) -> Result<f64, parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }

    fn set_parameter(
        &mut self,
        _id: usize,
        _parameter: f64,
    ) -> Result<(), parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }
}

struct Psi {}

impl ParameterizedFunction<1, 1> for Psi {
    fn calc(&mut self, x: &[f64; 1]) -> [f64; 1] {
        [x[0].powi(2)]
        // [(1.0 - (x[0] * 6.0 *std::f64::consts::PI).cos()) / 2.0]
    }

    fn get_parameter(
        &self,
        _id: usize,
    ) -> Result<f64, parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }

    fn set_parameter(
        &mut self,
        _id: usize,
        _parameter: f64,
    ) -> Result<(), parameterized_function::ParameterizedFunctionError> {
        Err(parameterized_function::ParameterizedFunctionError::WrongParameterId)
    }
}

fn main() {
    let mut problem: PartialDerivativeEquation = PartialDerivativeEquation {
        psi: Box::new(Psi {}),
        phi: Box::new(Phi {}),
        f: Box::new(F {}),
        a: 2.0,
        x_f: 0.0,
        x_l: 2.0,
        t_f: 0.0,
        t_l: 2.0,
    };

    for h_pow in 2..7 {
        for tau_pow in 2..7 {
            let h: f64 = 0.25f64.powi(h_pow);
            let tau: f64 = 0.25f64.powi(tau_pow);

            let t_n: usize = ((problem.t_l - problem.t_f) / tau).ceil() as usize;
            let x_n: usize = ((problem.x_l - problem.x_f) / h).ceil() as usize;

            let mut u: Vec<Vec<f64>> = vec![vec![0.0; x_n]; t_n];

            let start_time: std::time::Instant = std::time::Instant::now();
            for t_i in 0..t_n {
                for x_i in 0..x_n {
                    let t: f64 = t_i as f64 * tau + problem.t_f;
                    let x: f64 = x_i as f64 * h + problem.x_f;
                    u[t_i][x_i] = 0.0;
                    if x - t * problem.a >= 0.0 {
                        u[t_i][x_i] += problem.phi.calc(&[x - t * problem.a])[0];
                    }
                    if t - x / problem.a >= 0.0 {
                        u[t_i][x_i] += problem.psi.calc(&[t - x / problem.a])[0];
                    }
                    u[t_i][x_i] += x.powi(2) * t;
                }
            }

            write_csv(
                format!("analytical_{h_pow}_{tau_pow}"),
                PartialDerivativeSolution { solution: u.clone() },
                tau,
                h,
                0.0,
                start_time.elapsed(),
            );

            let u_analytical: Vec<Vec<f64>> = u;

            let solver: PDESolverRectangle = PDESolverRectangle {};

            let start_time: std::time::Instant = std::time::Instant::now();
            let solution: PartialDerivativeSolution = solver.solve(&mut problem, tau, h);

            let mut error: f64 = 0.0;
            let u: Vec<Vec<f64>> = solution.solution.clone();

            for i in 0..u.len() {
                for j in 0..u[i].len() {
                    error = if error < (u[i][j] - u_analytical[i][j]).abs() {
                        (u[i][j] - u_analytical[i][j]).abs()
                    } else {
                        error
                    }
                }
            }

            write_csv(
                format!("rectangle_{h_pow}_{tau_pow}"),
                solution,
                tau,
                h,
                error,
                start_time.elapsed(),
            );
        }
    }
}
