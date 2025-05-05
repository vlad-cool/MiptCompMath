use crate::algebraic_equation_solvers::*;
use crate::cauchy_problem::*;
use crate::matrix_operations::*;
use crate::utils;

pub enum NordsieckMethodType {
    ImplicitAdams(u32),
    ImplicitBackwardDifferentiation(u32),
}

impl NordsieckMethodType {
    pub fn get_l(&self) -> Option<std::vec::Vec<f64>> {
        match self {
            Self::ImplicitAdams(n) => match n {
                1 => Some(vec![1.0 / 2.0, 1.0]),
                2 => Some(vec![5.0 / 12.0, 1.0, 1.0 / 2.0]),
                3 => Some(vec![3.0 / 8.0, 1.0, 3.0 / 4.0, 1.0 / 6.0]),
                4 => Some(vec![251.0 / 720.0, 1.0, 11.0 / 12.0, 1.0 / 3.0, 1.0 / 24.0]),
                5 => Some(vec![
                    95.0 / 288.0,
                    1.0,
                    25.0 / 24.0,
                    35.0 / 72.0,
                    5.0 / 48.0,
                    1.0 / 120.0,
                ]),
                6 => Some(vec![
                    19087.0 / 60480.0,
                    1.0,
                    137.0 / 120.0,
                    5.0 / 8.0,
                    17.0 / 96.0,
                    1.0 / 40.0,
                    1.0 / 720.0,
                ]),
                _ => None,
            },
            Self::ImplicitBackwardDifferentiation(n) => match n {
                1 => Some(vec![1.0, 1.0]),
                2 => Some(vec![2.0 / 3.0, 1.0, 1.0 / 3.0]),
                3 => Some(vec![6.0 / 11.0, 1.0, 6.0 / 11.0, 1.0 / 11.0]),
                4 => Some(vec![12.0 / 15.0, 1.0, 7.0 / 10.0, 1.0 / 5.0, 1.0 / 50.0]),
                5 => Some(vec![
                    60.0 / 137.0,
                    1.0,
                    225.0 / 274.0,
                    85.0 / 274.0,
                    15.0 / 274.0,
                    1.0 / 274.0,
                ]),
                6 => Some(vec![
                    20.0 / 49.0,
                    1.0,
                    58.0 / 63.0,
                    5.0 / 12.0,
                    25.0 / 252.0,
                    1.0 / 84.0,
                    1.0 / 1764.0,
                ]),
                _ => None,
            },
        }
    }

    pub fn get_name(&self) -> String {
        match self {
            Self::ImplicitAdams(n) => format!(
                "Nordsieck representation of implicit Adams method of order {}",
                n
            ),
            Self::ImplicitBackwardDifferentiation(n) => format!(
                "Nordsieck representation of implicit backward differentiation method of order {}",
                n
            ),
        }
    }
}

pub struct NordsieckMethod<const N: usize> {
    l: std::vec::Vec<f64>,
    name: String,
}

impl<const N: usize> NordsieckMethod<N> {
    pub fn new(method: NordsieckMethodType) -> Self {
        Self {
            l: method.get_l().expect("No such method is known"),
            name: method.get_name(),
        }
    }

    fn step<F>(
        &mut self,
        problem: &mut CauchyProblem<N, F>,
        tau: f64,
        t: f64,
        z: &std::vec::Vec<f64>,
    ) -> Result<(f64, std::vec::Vec<f64>), &'static str>
    where
        F: FnMut(f64, &[f64; N]) -> [f64; N],
    {
        let mut equation = |z_next: &std::vec::Vec<f64>| {
            let k: usize = self.l.len();
            let mut x_next: [f64; N] = [0.0; N];

            for i in 0..N {
                x_next[i] = z_next[i];
            }

            let mut p: Vec<Vec<f64>> = vec![vec![0.0; k]; k];

            for j in 0..k {
                for i in 0..(j + 1) {
                    p[i][j] = utils::c_n_k(j, i) as f64;
                }
            }

            // println!("p: {:?}", p);
            // println!("p: {p:?}");
            // let p: Vec<Vec<f64>> = matrix_transpose(&p);

            let mut identity: Vec<Vec<f64>> = vec![vec![0.0; N]; N];

            for i in 0..N {
                identity[i][i] = 1.0;
            }

            let mut e: Vec<Vec<f64>> = vec![vec![0.0; k]];

            // e[0][1] = 1.0;
            e[0][1] = 1.0;

            // let e: Vec<Vec<f64>> = matrix_transpose(&e);

            let z_col: Vec<Vec<f64>> = matrix_transpose(&vec![z.clone()]);

            // println!("P: {:?}", p);
            // println!("R: {:?}", identity);
            // println!("z_col: {:?}", z_col);
            // println!("E tim_cirx: {:?}", tensor_product(&matrix_product(&e, &p), &identity));
            // println!("e_1: {:?}", e);

            // println!("P: {}x{}", p.len(), p[0].len());
            // println!("E: {}x{}", identity.len(), identity[0].len());
            // println!("z_col: {}x{}", z_col.len(), z_col[0].len());
            // println!("E tim_cirx: {}x{}", tensor_product(&matrix_product(&e, &p), &identity).len(), tensor_product(&matrix_product(&e, &p), &identity)[0].len());
            // println!("e_1: {}x{}", e.len(), e[0].len());

            // std::io::stdout().flush().unwrap();

            // let a: Vec<Vec<f64>> = matrix_product(&tensor_product(&p, &identity), &z_col);
            // let b: Vec<Vec<f64>> =
            //     tensor_product(&matrix_transpose(&vec![self.l.clone()]), &identity);
            // let c: Vec<Vec<f64>> = matrix_transpose(&matrix_scale(
            //     &vec![(problem.f)(t + tau, &x_next).to_vec()],
            //     tau,
            // ));
            // let d: Vec<Vec<f64>> =
            //     matrix_product(&tensor_product(&matrix_product(&e, &p), &identity), &z_col);

            // // let b: Vec<Vec<f64>> = matrix_transpose(&b);
            // let e: Vec<Vec<f64>> = matrix_product(&b, &matrix_sub(&c, &d));
            // // let f: Vec<Vec<f64>> = matrix_sum(&a, &e);
            // // let g: Vec<Vec<f64>> = matrix_sub(&f, &z_col);

            // let q: Vec<Vec<f64>> = matrix_product(&p, &z_col);

            let l: Vec<f64> = self.l.clone();

            // println!("p: {:?}", p);
            // println!("z_col: {:?}", z_col);

            matrix_sub(
                &matrix_transpose(&matrix_sum(
                    &matrix_product(&tensor_product(&p, &identity), &z_col), // (P @ E) z_n
                    &matrix_product(
                        &tensor_product(&matrix_transpose(&vec![l]), &identity), // l @ E
                        &matrix_sub(
                            &matrix_scale(
                                &matrix_transpose(&vec![(problem.f)(t + tau, &x_next).to_vec()]),
                                tau,
                            ), // h f(t_(n+1), x_(n+1))
                            &matrix_product(
                                &tensor_product(&matrix_product(&e, &p), &identity),
                                &z_col,
                            ), // (e_1 P @ E) x_n
                        ), // h f(t_(n+1), x_(n+1)) - (e_1 P @ E) x_n
                    ), // (l @ E) (h f(t_(n+1), x_(n+1)) - (e_1 P @ E) x_n)
                )),
                &vec![z_next.to_vec()],
            )[0]
            .clone()
        };

        // println!("1: {:?}", equation(&vec![1.0, 0.0]));
        // println!("2: {:?}", equation(&vec![0.0, 1.0]));

        // let mut z_flat: Vec<f64> = vec![0.0; N * z.len()];

        // for i in 0..z.len() {
        //     for j in 0..N {
        //         z_flat[i * N + j] = z[i][j];
        //     }
        // }

        // let z_new_flat: Vec<f64> = match solve_newton_vec(&equation, &z_flat, None) {
        // println!("{:?}", z);
        // println!("{:?}", equation(&z));
        // println!("AAA BBB CCC {:?}", equation(&vec![0.0, 0.0]));
        // println!("AAA BBB CCC {:?}", equation(&vec![0.0, 1.0]));
        // println!("AAA BBB CCC {:?}", equation(&vec![0.0, -1.0]));

        // let mut z_start = z.to_vec();
        // for i in 0..z_start.len() {
        //     z_start[i] *= 0.0;
        // }

        let z_new: Vec<f64> = match solve_newton_vec(&mut equation, z, None) {
            Ok(x) => x,
            Err(err) => {
                println!("Failed to solve, {}", err);
                return Err("Failed to solve");
            }
        };

        Ok((t + tau, z_new))
    }
}

impl<const N: usize, F> CauchySolver<N, F> for NordsieckMethod<N>
where
    F: FnMut(f64, &[f64; N]) -> [f64; N],
{
    fn solve(
        &mut self,
        problem: &mut CauchyProblem<N, F>,
        tau: f64,
        print_progress: bool,
        save_every: Option<u32>,
    ) -> (CauchySolution<N>, Result<(), &'static str>) {
        let mut t_i: Vec<f64> = vec![problem.start];
        let mut iterations: u32 = 0u32;
        let save_every = save_every.unwrap_or(1);

        let mut z: Vec<f64> = vec![0.0; self.l.len() * N];

        for i in 0..N {
            z[N * 0 + i] = problem.x_0[i];
            z[N * 1 + i] = tau * (problem.f)(problem.start, &problem.x_0)[i];
        }

        let mut solution: CauchySolution<N> = CauchySolution {
            t: vec![],
            x: vec![],
            method_name: <NordsieckMethod<N> as CauchySolver<N, F>>::get_name(self),
        };

        solution.t.push(problem.start);
        solution.x.push(problem.x_0.clone());
        if print_progress {
            print!("t: {:>20.10}, iterations: {}", problem.start, iterations)
        }
        iterations += 1;

        while *t_i.last().unwrap() < problem.stop {
            let res: Result<(f64, Vec<f64>), &str> =
                self.step(problem, tau, *t_i.last().unwrap(), &z);

            match res {
                Ok((t, z_new)) => {
                    t_i.push(t);
                    z = z_new;
                    if iterations % save_every == 0 {
                        solution.t.push(t);
                        let mut x: [f64; N] = [0.0; N];
                        for i in 0..N {
                            x[i] = z[i]
                        }
                        solution.x.push(x);
                        if print_progress {
                            print!("\rt: {:>10.6}, iterations: {:>16}", t, iterations)
                        }
                    }
                }
                Err(a) => {
                    println!("Failed to solve, reason: {}", a);
                    if print_progress {
                        println!("\n\n");
                    }
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
        self.name.clone()
    }
}
