use crate::cauchy_problem::*;
use crate::algebraic_equation_solvers;
use crate::utils;

pub struct RungeKuttaMethod<const N: usize, const M: usize, const MN: usize> {
    order: u32,
    solver_type: SolverType,
    a: [[f64; M]; M],
    b: [f64; M],
    c: [f64; M],
    name: String,
}

impl<const N: usize, const M: usize, const NM: usize> RungeKuttaMethod<N, M, NM> {
    pub fn new(order: u32, a: [[f64; M]; M], b: [f64; M], c: [f64; M], name: String) -> Self {
        let mut solver_type = SolverType::Explicit;
        for i in 0..M {
            for j in i..M {
                if a[i][j] != 0f64 {
                    solver_type = SolverType::Implicit;
                }
            }
        }

        Self {
            order,
            solver_type,
            a,
            b,
            c,
            name,
        }
    }

    fn step_explicit(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: f64,
        x: [f64; N],
    ) -> Result<(f64, [f64; N]), &'static str> {
        let mut k: [[f64; N]; M] = [[0f64; N]; M];
        for i in 0..M {
            let arg_1: f64 = t + tau * self.c[i];
            let mut arg_2: [f64; N] = [0f64; N];

            for a in 0..N {
                for j in 0..i {
                    arg_2[a] += self.a[i][j] * k[j][a];
                }
                arg_2[a] *= tau;
                arg_2[a] += x[a];
            }

            k[i] = (problem.f)(arg_1, &arg_2);
        }

        let mut res: [f64; N] = [0f64; N];
        for j in 0..N {
            for i in 0..M {
                res[j] += self.b[i] * k[i][j];
            }
            res[j] *= tau;
            res[j] += x[j];
        }

        Ok((t + tau, res))
    }

    fn step_implicit(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: f64,
        x: [f64; N],
    ) -> Result<(f64, [f64; N]), &'static str> {
        let equation = |k_0: &[f64; NM]| {
            let k_0: [[f64; N]; M] = utils::unflatten::<N, M>(k_0);
            let mut k_i: [[f64; N]; M] = [[0.0; N]; M];

            for i in 0..M {
                let arg_1: f64 = t + tau * self.c[i];
                let mut arg_2: [f64; N] = [0f64; N];

                for a in 0..N {
                    for j in 0..M {
                        arg_2[a] += self.a[i][j] * k_0[j][a];
                    }
                    arg_2[a] *= tau;
                    arg_2[a] += x[a];
                }
                for a in 0..N {
                    k_i[i][a] = (problem.f)(arg_1, &arg_2)[a] - k_0[i][a];
                }
            }

            let k_i: [f64; NM] = (*k_i.as_flattened()).try_into().unwrap();

            k_i
        };

        let k = match algebraic_equation_solvers::solve_newton(equation, &[0f64; NM], None) {
            Ok(x) => utils::unflatten::<N, M>(&x),
            Err(err) => {
                println!("Failed to solve, {}", err);
                return Err("Failed to solve");
            }
        };

        let mut res: [f64; N] = [0f64; N];
        for j in 0..N {
            for i in 0..M {
                res[j] += self.b[i] * k[i][j];
            }
            res[j] *= tau;
            res[j] += x[j];
        }

        Ok((t + tau, res))
    }

    fn step(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: f64,
        x: [f64; N],
    ) -> Result<(f64, [f64; N]), &'static str> {
        match self.solver_type {
            SolverType::Explicit => self.step_explicit(problem, tau, t, x),
            SolverType::Implicit => self.step_implicit(problem, tau, t, x),
        }
    }
}

impl<const N: usize, const M: usize, const NM: usize> CauchySolver<N>
    for RungeKuttaMethod<N, M, NM>
{
    fn solve(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        print_progress: bool,
        save_every: Option<u32>,
    ) -> (CauchySolution<N>, Result<(), &'static str>) {
        let mut t_i: f64 = problem.start;
        let mut x_i: [f64; N] = problem.x_0.clone();
        let mut iterations: u32 = 0u32;
        let save_every = save_every.unwrap_or(1);

        let mut solution: CauchySolution<N> = CauchySolution {
            t: vec![],
            x: vec![],
            method_name: self.get_name(),
        };

        solution.t.push(problem.start);
        solution.x.push(problem.x_0.clone());
        if print_progress {
            print!("t: {:>20.10}, iterations: {}", problem.start, iterations)
        }
        iterations += 1;

        while t_i < problem.stop {
            match self.step(&problem, tau, t_i, x_i) {
                Ok((t, x)) => {
                    t_i = t;
                    x_i = x.clone();
                    if iterations % save_every == 0 {
                        solution.t.push(t);
                        solution.x.push(x);
                        if print_progress {
                            print!("\rt: {:>10.6}, iterations: {:>16}", t, iterations)
                        }
                    }
                }
                Err(a) => {
                    if print_progress {
                        println!("\n\n");
                    }
                    println!("Failed to solve, reason: {}", a);
                    return (solution, Err("Failed to solve"));
                }
            };
            iterations += 1;
        }

        if print_progress {
            println!("");
        }

        (solution, Ok(()))
    }

    fn get_name(&self) -> String {
        return format!(
            "{}, order: {}, type: {}",
            self.name.clone(),
            self.order,
            self.solver_type
        );
    }
}
