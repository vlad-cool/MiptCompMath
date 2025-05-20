use crate::pde_problem::*;

pub struct PDESolverRectangle {}

impl PartialDerivativeEquationSolver for PDESolverRectangle {
    fn solve(
        &self,
        problem: &mut PartialDerivativeEquation,
        tau: f64,
        h: f64,
    ) -> PartialDerivativeSolution {
        let t_n: usize = ((problem.t_l - problem.t_f) / tau).ceil() as usize;
        let x_n: usize = ((problem.x_l - problem.x_f) / h).ceil() as usize;

        let mut u: Vec<Vec<f64>> = vec![vec![0.0; x_n]; t_n];

        for x_i in 0..x_n {
            u[0][x_i] = problem.phi.calc(&[problem.x_f + x_i as f64 * h])[0];
        }

        for t_i in 0..t_n {
            u[t_i][0] = problem.psi.calc(&[problem.t_f + t_i as f64 * tau])[0];
        }

        for t_i in 1..t_n {
            for x_i in 1..x_n {
                u[t_i][x_i] = 0.0;

                // u[t_i][x_i] += problem.f.calc(&[
                //     problem.t_f + (t_i as f64 - 0.0) * tau,
                //     problem.x_f + (x_i as f64 - 0.0) * h,
                // ])[0];
                // u[t_i][x_i] += problem.f.calc(&[
                //     problem.t_f + (t_i as f64 - 1.0) * tau,
                //     problem.x_f + (x_i as f64 - 0.0) * h,
                // ])[0];
                // u[t_i][x_i] += problem.f.calc(&[
                //     problem.t_f + (t_i as f64 - 0.0) * tau,
                //     problem.x_f + (x_i as f64 - 1.0) * h,
                // ])[0];
                // u[t_i][x_i] += problem.f.calc(&[
                //     problem.t_f + (t_i as f64 - 1.0) * tau,
                //     problem.x_f + (x_i as f64 - 1.0) * h,
                // ])[0];
                // u[t_i][x_i] /= 4.0;

                u[t_i][x_i] += problem.f.calc(&[
                    problem.t_f + (t_i as f64 - 0.5) * tau,
                    problem.x_f + (x_i as f64 - 0.5) * h,
                ])[0];

                u[t_i][x_i] *= 2.0 * tau * h;

                u[t_i][x_i] += (u[t_i - 1][x_i - 1] + u[t_i - 1][x_i] - u[t_i][x_i - 1]) * h;
                u[t_i][x_i] +=
                    (u[t_i - 1][x_i - 1] - u[t_i - 1][x_i] + u[t_i][x_i - 1]) * problem.a * tau;
                u[t_i][x_i] /= problem.a * tau + h;
            }
        }

        PartialDerivativeSolution { solution: u }
    }
}

pub struct PDESolverAngleExplicitLeft {}

impl PartialDerivativeEquationSolver for PDESolverAngleExplicitLeft {
    fn solve(
        &self,
        problem: &mut PartialDerivativeEquation,
        tau: f64,
        h: f64,
    ) -> PartialDerivativeSolution {
        let t_n: usize = ((problem.t_l - problem.t_f) / tau).ceil() as usize;
        let x_n: usize = ((problem.x_l - problem.x_f) / h).ceil() as usize;

        let mut u: Vec<Vec<f64>> = vec![vec![0.0; x_n]; t_n];

        for x_i in 0..x_n {
            u[0][x_i] = problem.phi.calc(&[problem.x_f + x_i as f64 * h])[0];
        }

        for t_i in 0..t_n {
            u[t_i][0] = problem.psi.calc(&[problem.t_f + t_i as f64 * tau])[0];
        }

        for t_i in 1..t_n {
            for x_i in 1..x_n {
                u[t_i][x_i] = 0.0;

                u[t_i][x_i] += problem.f.calc(&[
                    problem.t_f + (t_i as f64 - 1.0) * tau,
                    problem.x_f + (x_i as f64 - 0.0) * h,
                ])[0];

                u[t_i][x_i] -= problem.a * (u[t_i - 1][x_i] - u[t_i - 1][x_i - 1]) / h;
                u[t_i][x_i] *= tau;
                u[t_i][x_i] += u[t_i - 1][x_i];
            }
        }

        PartialDerivativeSolution { solution: u }
    }
}

pub struct PDESolverLaxWendroff {}

impl PartialDerivativeEquationSolver for PDESolverLaxWendroff {
    fn solve(
        &self,
        problem: &mut PartialDerivativeEquation,
        tau: f64,
        h: f64,
    ) -> PartialDerivativeSolution {
        let t_n: usize = ((problem.t_l - problem.t_f) / tau).ceil() as usize;
        let x_n: usize = ((problem.x_l - problem.x_f) / h).ceil() as usize;

        let mut u: Vec<Vec<f64>> = vec![vec![0.0; x_n]; t_n];

        for x_i in 0..x_n {
            u[0][x_i] = problem.phi.calc(&[problem.x_f + x_i as f64 * h])[0];
        }

        for t_i in 0..t_n {
            u[t_i][0] = problem.psi.calc(&[problem.t_f + t_i as f64 * tau])[0];
        }

        for t_i in 1..t_n {
            for x_i in 1..(x_n - 1) {
                u[t_i][x_i] = 0.0;

                u[t_i][x_i] += problem.f.calc(&[
                    problem.t_f + (t_i as f64 - 1.0) * tau,
                    problem.x_f + (x_i as f64 - 0.0) * h,
                ])[0] * tau;

                u[t_i][x_i] -= ((problem.a * tau) / (2.0 * h)) * (u[t_i - 1][x_i + 1] - u[t_i - 1][x_i - 1]);
                u[t_i][x_i] += 0.5 * (u[t_i - 1][x_i + 1] + u[t_i - 1][x_i - 1]);
            }

            // Rectangle scheme for most right point
            let x_i: usize = x_n - 1;
            u[t_i][x_i] = 0.0;

            u[t_i][x_i] += problem.f.calc(&[
                problem.t_f + (t_i as f64 - 1.0) * tau,
                problem.x_f + (x_i as f64 - 0.0) * h,
            ])[0];

            u[t_i][x_i] -= problem.a * (u[t_i - 1][x_i] - u[t_i - 1][x_i - 1]) / h;
            u[t_i][x_i] *= tau;
            u[t_i][x_i] += u[t_i - 1][x_i];
        }

        PartialDerivativeSolution { solution: u }
    }
}
