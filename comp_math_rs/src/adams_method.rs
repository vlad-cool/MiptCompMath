use crate::cauchy_problem::*;
use crate::equation_solvers;

pub struct AdamsMethod<const N: usize> {
    order: usize,
    solver_type: SolverType,
}

impl<const N: usize> AdamsMethod<N> {
    pub fn new(order: usize, solver_type: SolverType) -> Self {
        Self { order, solver_type }
    }

    fn step_explicit(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: &std::vec::Vec<f64>,
        x: &std::vec::Vec<[f64; N]>,
        order: usize,
    ) -> Result<(f64, [f64; N]), &'static str> {
        if t.len() < order {
            self.step_explicit(problem, tau, t, x, t.len())
        } else {
            let n = t.len() - 1;
            match order {
                1 => {
                    let mut x_i: [f64; N] = [0.0; N];
                    for i in 0..N {
                        x_i[i] += 1.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                        x_i[i] /= 1.0;
                        x_i[i] *= tau;
                        x_i[i] += x[n][i];
                    }
                    Ok((t[n] + tau, x_i))
                }
                2 => {
                    let mut x_i: [f64; N] = [0.0; N];
                    for i in 0..N {
                        x_i[i] += 3.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                        x_i[i] -= 1.0 * (problem.f)(t[n - 1], &x[n - 1])[i];
                        x_i[i] /= 2.0;
                        x_i[i] *= tau;
                        x_i[i] += x[n][i];
                    }
                    Ok((t[n] + tau, x_i))
                }
                3 => {
                    let mut x_i: [f64; N] = [0.0; N];
                    for i in 0..N {
                        x_i[i] += 23.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                        x_i[i] -= 16.0 * (problem.f)(t[n - 1], &x[n - 1])[i];
                        x_i[i] += 05.0 * (problem.f)(t[n - 2], &x[n - 2])[i];
                        x_i[i] /= 12.0;
                        x_i[i] *= tau;
                        x_i[i] += x[n][i];
                    }
                    Ok((t[n] + tau, x_i))
                }
                4 => {
                    let mut x_i: [f64; N] = [0.0; N];
                    for i in 0..N {
                        x_i[i] += 55.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                        x_i[i] -= 59.0 * (problem.f)(t[n - 1], &x[n - 1])[i];
                        x_i[i] += 37.0 * (problem.f)(t[n - 2], &x[n - 2])[i];
                        x_i[i] -= 09.0 * (problem.f)(t[n - 3], &x[n - 3])[i];
                        x_i[i] /= 24.0;
                        x_i[i] *= tau;
                        x_i[i] += x[n][i];
                    }
                    Ok((t[n] + tau, x_i))
                }
                _ => Err("No method with such order"),
            }
        }
    }

    fn step_implicit(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: &std::vec::Vec<f64>,
        x: &std::vec::Vec<[f64; N]>,
        order: usize,
    ) -> Result<(f64, [f64; N]), &'static str> {
        if t.len() < order {
            self.step_implicit(problem, tau, t, x, t.len())
        } else {
            let n: usize = t.len() - 1;
            match order {
                1 => {
                    let equation = |x_next: &[f64; N]| {
                        let mut x_i: [f64; N] = [0.0; N];

                        for i in 0..N {
                            x_i[i] += 1.0 * (problem.f)(t[n - 0] + tau, x_next)[i];
                            x_i[i] /= 1.0;
                            x_i[i] *= tau;
                            x_i[i] += x[n][i];
                            x_i[i] -= x_next[i];
                        }

                        x_i
                    };

                    match equation_solvers::solve_newton(equation, &[0.0; N], None) {
                        Ok(x) => Ok((t[n] + tau, x)),
                        Err(err) => {
                            println!("Failed to solve, {}", err);
                            Err("Failed to solve")
                        }
                    }
                }
                2 => {
                    let equation = |x_next: &[f64; N]| {
                        let mut x_i: [f64; N] = [0.0; N];

                        for i in 0..N {
                            x_i[i] += 1.0 * (problem.f)(t[n - 0] + tau, x_next)[i];
                            x_i[i] += 1.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                            x_i[i] /= 2.0;
                            x_i[i] *= tau;
                            x_i[i] += x[n][i];
                            x_i[i] -= x_next[i];
                        }

                        x_i
                    };

                    match equation_solvers::solve_newton(equation, &x[n - 1], None) {
                        Ok(x) => Ok((t[n] + tau, x)),
                        Err(err) => {
                            println!("Failed to solve, {}", err);
                            Err("Failed to solve")
                        }
                    }
                }
                3 => {
                    let equation = |x_next: &[f64; N]| {
                        let mut x_i: [f64; N] = [0.0; N];

                        for i in 0..N {
                            x_i[i] += 5.0 * (problem.f)(t[n - 0] + tau, x_next)[i];
                            x_i[i] += 8.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                            x_i[i] -= 1.0 * (problem.f)(t[n - 1], &x[n - 1])[i];
                            x_i[i] /= 12.0;
                            x_i[i] *= tau;
                            x_i[i] += x[n][i];
                            x_i[i] -= x_next[i];
                        }

                        x_i
                    };

                    match equation_solvers::solve_newton(equation, &x[n - 1], None) {
                        Ok(x) => Ok((t[n] + tau, x)),
                        Err(err) => {
                            println!("Failed to solve, {}", err);
                            Err("Failed to solve")
                        }
                    }
                }
                4 => {
                    let equation = |x_next: &[f64; N]| {
                        let mut x_i: [f64; N] = [0.0; N];

                        for i in 0..N {
                            x_i[i] += 09.0 * (problem.f)(t[n - 0] + tau, x_next)[i];
                            x_i[i] += 19.0 * (problem.f)(t[n - 0], &x[n - 0])[i];
                            x_i[i] -= 05.0 * (problem.f)(t[n - 1], &x[n - 1])[i];
                            x_i[i] += 01.0 * (problem.f)(t[n - 2], &x[n - 2])[i];
                            x_i[i] /= 24.0;
                            x_i[i] *= tau;
                            x_i[i] += x[n][i];
                            x_i[i] -= x_next[i];
                        }

                        x_i
                    };

                    match equation_solvers::solve_newton(equation, &x[n - 1], None) {
                        Ok(x) => Ok((t[n] + tau, x)),
                        Err(err) => {
                            println!("Failed to solve, {}", err);
                            Err("Failed to solve")
                        }
                    }
                }
                _ => Err("No method with such order"),
            }
        }
    }
}

impl<const N: usize> CauchySolver<N> for AdamsMethod<N> {
    fn solve(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        print_progress: bool,
        save_every: Option<u32>,
    ) -> (CauchySolution<N>, Result<(), &'static str>) {
        let mut t_i: Vec<f64> = vec![problem.start];
        let mut x_i: Vec<[f64; N]> = vec![problem.x_0.clone()];
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

        while *t_i.last().unwrap() < problem.stop {
            let res: Result<(f64, [f64; N]), &str> = match self.solver_type {
                SolverType::Explicit => self.step_explicit(&problem, tau, &t_i, &x_i, self.order),
                SolverType::Implicit => self.step_implicit(&problem, tau, &t_i, &x_i, self.order),
            };

            match res {
                Ok((t, x)) => {
                    t_i.push(t);
                    x_i.push(x);
                    if t_i.len() > self.order {
                        t_i.remove(0);
                        x_i.remove(0);
                    }
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
        return format!("{} Adams method of order: {}", self.solver_type, self.order,);
    }
}
