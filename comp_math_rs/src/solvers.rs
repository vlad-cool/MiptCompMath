// #![allow(dead_code)]
// use crate::cauchy_problem::*;
// use std::{f64, io::Write, vec};

// pub fn zero_enough_vec(a: &std::vec::Vec<f64>, epsilon: f64) -> bool {
//     let mut sum: f64 = 0f64;
//     for i in 0..a.len() {
//         sum += a[i].abs();
//     }

//     close_enough(sum, 0f64, epsilon)
// }

// pub fn partial_derivative_vec<F>(mut f: F, x: &vec::Vec<f64>, i: usize, j: usize, h: f64) -> f64
// where
//     F: FnMut(&vec::Vec<f64>) -> vec::Vec<f64>,
// {
//     let mut x_1: Vec<f64> = x.clone();
//     let mut x_2: Vec<f64> = x.clone();

//     x_1[j] += h;
//     x_2[j] -= h;

//     (f(&x_1)[i] - f(&x_2)[i]) / (2f64 * h)
// }

// pub fn solve_linear_system_vec(a: &vec::Vec<vec::Vec<f64>>, b: &vec::Vec<f64>) -> vec::Vec<f64> {
//     let n = b.len();

//     let mut l: Vec<Vec<f64>> = vec![vec![0f64; n]; n];
//     let mut u: Vec<Vec<f64>> = vec![vec![0f64; n]; n];

//     for i in 0..n {
//         l[i][i] = 1f64;
//         for j in i..n {
//             let mut sum: f64 = 0f64;
//             for k in 0..i {
//                 sum += l[i][k] * u[k][j];
//             }
//             u[i][j] = a[i][j] - sum;
//         }

//         for j in i..n {
//             let mut sum: f64 = 0f64;
//             for k in 0..i {
//                 sum += l[j][k] * u[k][i];
//             }
//             l[j][i] = (a[j][i] - sum) / u[i][i];
//         }
//     }

//     let mut v: Vec<f64> = vec![0f64; n];
//     for i in 0..n {
//         v[i] = b[i];
//         for j in 0..i {
//             v[i] -= l[i][j] * v[j];
//         }
//     }

//     let mut x: Vec<f64> = vec![0f64; n];

//     for i in (0..n).rev() {
//         x[i] = v[i] / u[i][i];
//         for j in i + 1..n {
//             x[i] -= u[i][j] * x[j] / u[i][i];
//         }
//     }

//     x
// }

// fn solve_newton_vec<F>(
//     f: &F,
//     x: &vec::Vec<f64>,
//     max_iterations: Option<u32>,
// ) -> Result<vec::Vec<f64>, &'static str>
// where
//     F: Fn(&vec::Vec<f64>) -> vec::Vec<f64>,
// {
//     let max_iterations: u32 = max_iterations.unwrap_or(100u32);
//     let n = x.len();

//     let mut x: Vec<f64> = x.clone();
//     let mut iterations: u32 = 0u32;

//     while !zero_enough_vec(&f(&x), 1e-6) {
//         if iterations == max_iterations {
//             return Err("Cannot solve system, too many iterations");
//         }

//         let mut jacobian: Vec<Vec<f64>> = vec![vec![0.0; n]; n];

//         for i in 0..n {
//             for j in 0..n {
//                 // jacobian[j][i] = -partial_derivative_vec(f, &x, i, j, 1e-5);
//                 jacobian[i][j] = -partial_derivative_vec(f, &x, i, j, 1e-5);
//             }
//         }

//         let delta_x: Vec<f64> = solve_linear_system_vec(&jacobian, &f(&x));
//         println!("jac: {:?}", jacobian);
//         println!("f: {:?}", f(&x));
//         println!("{:?}", delta_x);

//         println!("{:?}", jacobian);
//         println!("{:?}", f(&x));
//         println!("{:?}", delta_x);

//         for i in 0..n {
//             x[i] += delta_x[i];
//         }

//         iterations += 1;
//     }

//     Ok(x)
// }

// fn tensor_product(
//     l: &vec::Vec<vec::Vec<f64>>,
//     r: &vec::Vec<vec::Vec<f64>>,
// ) -> vec::Vec<vec::Vec<f64>> {
//     assert!(l.len() > 0);
//     assert!(r.len() > 0);
//     for i in 0..(l.len() - 1) {
//         assert!(l[i].len() == l[i + 1].len());
//     }
//     for i in 0..(r.len() - 1) {
//         assert!(r[i].len() == r[i + 1].len());
//     }

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len() * r[0].len()]; l.len() * r.len()];

//     for a in 0..r.len() {
//         for b in 0..r[0].len() {
//             for c in 0..l.len() {
//                 for d in 0..l[0].len() {
//                     res[a * l.len() + c][b * l[0].len() + d] = l[c][d] * r[a][b];
//                 }
//             }
//         }
//     }

//     res
// }

// fn matrix_product(
//     l: &vec::Vec<vec::Vec<f64>>,
//     r: &vec::Vec<vec::Vec<f64>>,
// ) -> vec::Vec<vec::Vec<f64>> {
//     assert!(l.len() > 0);
//     assert!(r.len() > 0);
//     for i in 0..(l.len() - 1) {
//         assert!(l[i].len() == l[i + 1].len());
//     }
//     for i in 0..(r.len() - 1) {
//         assert!(r[i].len() == r[i + 1].len());
//     }

//     assert!(l[0].len() == r.len());

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; r[0].len()]; l.len()];

//     for a in 0..l.len() {
//         for b in 0..r[0].len() {
//             for c in 0..r.len() {
//                 res[a][b] += r[c][b] * l[a][c];
//             }
//         }
//     }

//     res
// }

// fn matrix_sum(l: &vec::Vec<vec::Vec<f64>>, r: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
//     assert!(l.len() > 0);
//     assert!(r.len() > 0);
//     assert!(l.len() == r.len());
//     assert!(l[0].len() == r[0].len());
//     for i in 0..(l.len() - 1) {
//         assert!(l[i].len() == l[i + 1].len());
//     }
//     for i in 0..(r.len() - 1) {
//         assert!(r[i].len() == r[i + 1].len());
//     }

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

//     for a in 0..l.len() {
//         for b in 0..l[0].len() {
//             res[a][b] = l[a][b] + r[a][b];
//         }
//     }

//     res
// }

// fn matrix_sub(l: &vec::Vec<vec::Vec<f64>>, r: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
//     assert!(l.len() > 0);
//     assert!(r.len() > 0);
//     assert!(l.len() == r.len());
//     assert!(l[0].len() == r[0].len());
//     for i in 0..(l.len() - 1) {
//         assert!(l[i].len() == l[i + 1].len());
//     }
//     for i in 0..(r.len() - 1) {
//         assert!(r[i].len() == r[i + 1].len());
//     }

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

//     for a in 0..l.len() {
//         for b in 0..l[0].len() {
//             res[a][b] = l[a][b] - r[a][b];
//         }
//     }

//     res
// }

// fn matrix_transpose(a: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
//     for i in 0..(a.len() - 1) {
//         assert!(a[i].len() == a[i + 1].len());
//     }

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; a.len()]; a[0].len()];

//     for i in 0..a.len() {
//         for j in 0..a[0].len() {
//             res[j][i] = a[i][j];
//         }
//     }

//     res
// }

// fn matrix_scale(l: &vec::Vec<vec::Vec<f64>>, r: f64) -> vec::Vec<vec::Vec<f64>> {
//     for i in 0..(l.len() - 1) {
//         assert!(l[i].len() == l[i + 1].len());
//     }

//     let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

//     for i in 0..l.len() {
//         for j in 0..l[0].len() {
//             res[i][j] = l[i][j] * r;
//         }
//     }

//     res
// }

// fn c_n_k(n: usize, k: usize) -> usize {
//     if k > n {
//         return 0;
//     }

//     let k = std::cmp::min(k, n - k);
//     let mut res = 1;

//     for i in 1..=k {
//         res = res * (n - k + i) / i;
//     }

//     res
// }

// pub struct RosenbrockMethod<const N: usize, const M: usize, const MN: usize> {
//     order: u32,
//     solver_type: SolverType,
//     a: [[f64; M]; M],
//     b: [f64; M],
//     c: [f64; M],
//     name: String,
// }

// impl<const N: usize, const M: usize, const NM: usize> RosenbrockMethod<N, M, NM> {
//     pub fn new(order: u32, a: [[f64; M]; M], b: [f64; M], c: [f64; M], name: String) -> Self {
//         todo!()
//         // let mut solver_type = SolverType::Explicit;
//         // for i in 0..M {
//         //     for j in i..M {
//         //         if a[i][j] != 0f64 {
//         //             solver_type = SolverType::Implicit;
//         //         }
//         //     }
//         // }

//         // Self {
//         //     order,
//         //     solver_type,
//         //     a,
//         //     b,
//         //     c,
//         //     name,
//         // }
//     }

//     fn step(
//         &mut self,
//         problem: &CauchyProblem<N>,
//         tau: f64,
//         t: f64,
//         x: [f64; N],
//     ) -> Result<(f64, [f64; N]), &'static str> {
//         todo!()
//         // let equation = |k_0: &[f64; NM]| {
//         //     let k_0: [[f64; N]; M] = unflatten::<N, M>(k_0);
//         //     let mut k_i: [[f64; N]; M] = [[0.0; N]; M];

//         //     for i in 0..M {
//         //         let arg_1: f64 = t + tau * self.c[i];
//         //         let mut arg_2: [f64; N] = [0f64; N];

//         //         for a in 0..N {
//         //             for j in 0..M {
//         //                 arg_2[a] += self.a[i][j] * k_0[j][a];
//         //             }
//         //             arg_2[a] *= tau;
//         //             arg_2[a] += x[a];
//         //         }
//         //         for a in 0..N {
//         //             k_i[i][a] = (problem.f)(arg_1, &arg_2)[a] - k_0[i][a];
//         //         }
//         //     }

//         //     let k_i: [f64; NM] = (*k_i.as_flattened()).try_into().unwrap();

//         //     k_i
//         // };

//         // let k = match solve_newton(equation, &[0f64; NM], None) {
//         //     Ok(x) => unflatten::<N, M>(&x),
//         //     Err(err) => {
//         //         println!("Failed to solve, {}", err);
//         //         return Err("Failed to solve");
//         //     }
//         // };

//         // let mut res: [f64; N] = [0f64; N];
//         // for j in 0..N {
//         //     for i in 0..M {
//         //         res[j] += self.b[i] * k[i][j];
//         //     }
//         //     res[j] *= tau;
//         //     res[j] += x[j];
//         // }

//         // Ok((t + tau, res))
//     }
// }

// impl<const N: usize, const M: usize, const NM: usize> CauchySolver<N>
//     for RosenbrockMethod<N, M, NM>
// {
//     fn solve(
//         &mut self,
//         problem: &CauchyProblem<N>,
//         tau: f64,
//         print_progress: bool,
//         save_every: Option<u32>,
//     ) -> (CauchySolution<N>, Result<(), &'static str>) {
//         todo!()
//         // let mut t_i: f64 = problem.start;
//         // let mut x_i: [f64; N] = problem.x_0.clone();
//         // let mut iterations: u32 = 0u32;
//         // let save_every = save_every.unwrap_or(1);

//         // let mut solution: CauchySolution<N> = CauchySolution {
//         //     t: vec![],
//         //     x: vec![],
//         //     method_name: self.get_name(),
//         // };

//         // solution.t.push(problem.start);
//         // solution.x.push(problem.x_0.clone());
//         // if print_progress {
//         //     print!("t: {:>20.10}, iterations: {}", problem.start, iterations)
//         // }
//         // iterations += 1;

//         // while t_i < problem.stop {
//         //     match self.step(&problem, tau, t_i, x_i) {
//         //         Ok((t, x)) => {
//         //             t_i = t;
//         //             x_i = x.clone();
//         //             if iterations % save_every == 0 {
//         //                 solution.t.push(t);
//         //                 solution.x.push(x);
//         //                 if print_progress {
//         //                     print!("\rt: {:>10.6}, iterations: {:>16}", t, iterations)
//         //                 }
//         //             }
//         //         }
//         //         Err(a) => {
//         //             if print_progress {
//         //                 println!("\n\n");
//         //             }
//         //             println!("Failed to solve, reason: {}", a);
//         //             return (solution, Err("Failed to solve"));
//         //         }
//         //     };
//         //     iterations += 1;
//         // }

//         // if print_progress {
//         //     println!("");
//         // }

//         // (solution, Ok(()))
//     }

//     fn get_name(&self) -> String {
//         return format!(
//             "{}, order: {}, type: {}",
//             self.name.clone(),
//             self.order,
//             self.solver_type
//         );
//     }
// }

// pub enum NordsieckMethodType {
//     ImplicitAdams(u32),
//     ImplicitBackwardDifferentiation(u32),
// }

// impl NordsieckMethodType {
//     pub fn get_l(&self) -> Option<vec::Vec<f64>> {
//         match self {
//             Self::ImplicitAdams(n) => match n {
//                 1 => Some(vec![1.0 / 2.0, 1.0]),
//                 2 => Some(vec![5.0 / 12.0, 1.0, 1.0 / 2.0]),
//                 3 => Some(vec![3.0 / 8.0, 1.0, 3.0 / 4.0, 1.0 / 6.0]),
//                 4 => Some(vec![251.0 / 720.0, 1.0, 11.0 / 12.0, 1.0 / 3.0, 1.0 / 24.0]),
//                 5 => Some(vec![
//                     95.0 / 288.0,
//                     1.0,
//                     25.0 / 24.0,
//                     35.0 / 72.0,
//                     5.0 / 48.0,
//                     1.0 / 120.0,
//                 ]),
//                 6 => Some(vec![
//                     19087.0 / 60480.0,
//                     1.0,
//                     137.0 / 120.0,
//                     5.0 / 8.0,
//                     17.0 / 96.0,
//                     1.0 / 40.0,
//                     1.0 / 720.0,
//                 ]),
//                 _ => None,
//             },
//             Self::ImplicitBackwardDifferentiation(n) => match n {
//                 1 => Some(vec![1.0, 1.0]),
//                 2 => Some(vec![2.0 / 3.0, 1.0, 1.0 / 3.0]),
//                 3 => Some(vec![6.0 / 11.0, 1.0, 6.0 / 11.0, 1.0 / 11.0]),
//                 4 => Some(vec![12.0 / 15.0, 1.0, 7.0 / 10.0, 1.0 / 5.0, 1.0 / 50.0]),
//                 5 => Some(vec![
//                     60.0 / 137.0,
//                     1.0,
//                     225.0 / 274.0,
//                     85.0 / 274.0,
//                     15.0 / 274.0,
//                     1.0 / 274.0,
//                 ]),
//                 6 => Some(vec![
//                     20.0 / 49.0,
//                     1.0,
//                     58.0 / 63.0,
//                     5.0 / 12.0,
//                     25.0 / 252.0,
//                     1.0 / 84.0,
//                     1.0 / 1764.0,
//                 ]),
//                 _ => None,
//             },
//         }
//     }

//     pub fn get_name(&self) -> String {
//         match self {
//             Self::ImplicitAdams(n) => format!(
//                 "Nordsieck representation of implicit Adams method of order {}",
//                 n
//             ),
//             Self::ImplicitBackwardDifferentiation(n) => format!(
//                 "Nordsieck representation of implicit backward differentiation method of order {}",
//                 n
//             ),
//         }
//     }
// }

// pub struct NordsieckMethod<const N: usize> {
//     l: vec::Vec<f64>,
//     name: String,
// }

// impl<const N: usize> NordsieckMethod<N> {
//     pub fn new(method: NordsieckMethodType) -> Self {
//         Self {
//             l: method.get_l().expect("No such method is known"),
//             name: method.get_name(),
//         }
//     }

//     fn step(
//         &mut self,
//         problem: &CauchyProblem<N>,
//         tau: f64,
//         t: f64,
//         z: &vec::Vec<f64>,
//     ) -> Result<(f64, vec::Vec<f64>), &'static str> {
//         let equation = |z_next: &vec::Vec<f64>| {
//             let k: usize = self.l.len();
//             let mut x_next: [f64; N] = [0.0; N];

//             for i in 0..N {
//                 x_next[i] = z_next[i * N];
//             }

//             let mut p: Vec<Vec<f64>> = vec![vec![0.0; k]; k];

//             for j in 0..k {
//                 for i in 0..(j + 1) {
//                     p[i][j] = c_n_k(j, i) as f64;
//                 }
//             }

//             // println!("p: {:?}", p);
//             // println!("p: {p:?}");
//             // let p: Vec<Vec<f64>> = matrix_transpose(&p);

//             let mut identity: Vec<Vec<f64>> = vec![vec![0.0; N]; N];

//             for i in 0..N {
//                 identity[i][i] = 1.0;
//             }

//             let mut e: Vec<Vec<f64>> = vec![vec![0.0; k]];

//             // e[0][1] = 1.0;
//             e[0][1] = 1.0;

//             let z_col: Vec<Vec<f64>> = matrix_transpose(&vec![z.clone()]);

//             // println!("P: {:?}", p);
//             // println!("R: {:?}", identity);
//             // println!("z_col: {:?}", z_col);
//             // println!("E tim_cirx: {:?}", tensor_product(&matrix_product(&e, &p), &identity));
//             // println!("e_1: {:?}", e);

//             // println!("P: {}x{}", p.len(), p[0].len());
//             // println!("E: {}x{}", identity.len(), identity[0].len());
//             // println!("z_col: {}x{}", z_col.len(), z_col[0].len());
//             // println!("E tim_cirx: {}x{}", tensor_product(&matrix_product(&e, &p), &identity).len(), tensor_product(&matrix_product(&e, &p), &identity)[0].len());
//             // println!("e_1: {}x{}", e.len(), e[0].len());

//             // std::io::stdout().flush().unwrap();

//             // let a: Vec<Vec<f64>> = matrix_product(&tensor_product(&p, &identity), &z_col);
//             // let b: Vec<Vec<f64>> =
//             //     tensor_product(&matrix_transpose(&vec![self.l.clone()]), &identity);
//             // let c: Vec<Vec<f64>> = matrix_transpose(&matrix_scale(
//             //     &vec![(problem.f)(t + tau, &x_next).to_vec()],
//             //     tau,
//             // ));
//             // let d: Vec<Vec<f64>> =
//             //     matrix_product(&tensor_product(&matrix_product(&e, &p), &identity), &z_col);

//             // println!("###########################");
//             // println!("b: {:?}", b);
//             // println!("b: {}x{}", b.len(), b[0].len());
//             // println!("c: {:?}", c);
//             // println!("c: {}x{}", c.len(), c[0].len());
//             std::io::stdout().flush().unwrap();
//             let b: Vec<Vec<f64>> = matrix_transpose(&b);
//             // // println!("###########################");
//             // // println!("b: {:?}", b);
//             // // println!("b: {}x{}", b.len(), b[0].len());
//             // // println!("c: {:?}", c);
//             // // println!("c: {}x{}", c.len(), c[0].len());
//             // std::io::stdout().flush().unwrap();
//             // // let b: Vec<Vec<f64>> = matrix_transpose(&b);

//             // let e: Vec<Vec<f64>> = matrix_product(&b, &matrix_sub(&c, &d));
//             // let f: Vec<Vec<f64>> = matrix_sum(&a, &e);
//             // let g: Vec<Vec<f64>> = matrix_sub(&f, &z_col);
//             // matrix_transpose(&g)[0].clone()

//             let q: Vec<Vec<f64>> = matrix_product(&p, &z_col);

//             matrix_sub(
//                 &vec![z_next.clone()],
//                 &matrix_sum(&q, &matrix_scale(&vec![self.l.clone()], (tau * (problem.f)(t + tau, &[z[0]; N])[0] - matrix_product(&e, &q)[0][0]))),
//             )[0].clone()
//         };

//         let mut z_flat: Vec<f64> = vec![0.0; N * z.len()];

//         for i in 0..z.len() {
//             for j in 0..N {
//                 z_flat[i * N + j] = z[i][j]
//             }
//         }

//         let z_new_flat: Vec<f64> = match solve_newton_vec(&equation, &z_flat, None) {
//         println!("{:?}", z);
//         println!("{:?}", equation(&z));
//         println!("AAA BBB CCC {:?}", equation(&vec![0.0, 0.0]));
//         println!("AAA BBB CCC {:?}", equation(&vec![0.0, 1.0]));
//         println!("AAA BBB CCC {:?}", equation(&vec![0.0, -1.0]));

//         let z_new: Vec<f64> = match solve_newton_vec(&equation, &z, None) {
//             Ok(x) => x,
//             Err(err) => {
//                 println!("Failed to solve, {}", err);
//                 return Err("Failed to solve");
//             }
//         };

//         let mut z_new: Vec<[f64; N]> = vec![[0.0; N]; z.len()];

//         for i in 0..z.len() {
//             for j in 0..N {
//                 z_new[i][j] = z_new_flat[i * N + j]
//             }
//         }

//         Ok((t + tau, z_new))
//     }
// }

// impl<const N: usize> CauchySolver<N> for NordsieckMethod<N> {
//     fn solve(
//         &mut self,
//         problem: &CauchyProblem<N>,
//         tau: f64,
//         print_progress: bool,
//         save_every: Option<u32>,
//     ) -> (CauchySolution<N>, Result<(), &'static str>) {
//         let mut t_i: Vec<f64> = vec![problem.start];
//         let mut iterations: u32 = 0u32;
//         let save_every = save_every.unwrap_or(1);

//         let mut z: Vec<f64> = vec![0.0; self.l.len() * N];

//         for i in 0..N {
//             z[i * N + 0] = problem.x_0[i];
//             z[i * N + 1] = tau * (problem.f)(problem.start, &problem.x_0)[i];
//         }

//         let mut solution: CauchySolution<N> = CauchySolution {
//             t: vec![],
//             x: vec![],
//             method_name: self.get_name(),
//         };

//         solution.t.push(problem.start);
//         solution.x.push(problem.x_0.clone());
//         if print_progress {
//             print!("t: {:>20.10}, iterations: {}", problem.start, iterations)
//         }
//         iterations += 1;

//         while *t_i.last().unwrap() < problem.stop {
//             // let res: Result<(f64, [f64; N]), &str> = match self.solver_type {
//             //     SolverType::Explicit => self.step_explicit(&problem, tau, &t_i, &x_i, self.order),
//             //     SolverType::Implicit => self.step_implicit(&problem, tau, &t_i, &x_i, self.order),
//             // };

//             // let res = self.step(problem, tau, &t_i,&z);
//             let res: Result<(f64, Vec<f64>), &str> =
//                 self.step(problem, tau, *t_i.last().unwrap(), &z);

//             match res {
//                 Ok((t, z_new)) => {
//                     t_i.push(t);
//                     z = z_new;
//                     if iterations % save_every == 0 {
//                         solution.t.push(t);
//                         let mut x: [f64; N] = [0.0; N];
//                         for i in 0..N {
//                             x[i] = z[i * N]
//                         }
//                         solution.x.push(x);
//                         if print_progress {
//                             print!("\rt: {:>10.6}, iterations: {:>16}", t, iterations)
//                         }
//                     }
//                 }
//                 Err(a) => {
//                     println!("Failed to solve, reason: {}", a);
//                     if print_progress {
//                         println!("\n\n");
//                     }
//                     return (solution, Err("Failed to solve"));
//                 }
//             };
//             iterations += 1;
//         }

//         if print_progress {
//             println!("");
//         }

//         (solution, Ok(()))
//     }

//     fn get_name(&self) -> String {
//         self.name.clone()
//     }
// }

// /// Tests

// #[cfg(test)]
// mod tests {
//     use crate::solvers::*;

//     #[test]
//     fn test_solve_linear_system() {
//         const N: usize = 3;
//         let a: [[f64; 3]; 3] = [
//             [-4f64, 9f64, -4f64],
//             [-5f64, -5f64, 6f64],
//             [2f64, 5f64, -8f64],
//         ];

//         let b: [f64; 3] = [-64f64, 104.6f64, -85.2f64];

//         let x: [f64; 3] = crate::solvers::solve_linear_system(&a, &b);
//         for i in 0..N {
//             let mut sum: f64 = 0f64;
//             for j in 0..N {
//                 sum += a[i][j] * x[j];
//             }
//             assert!(close_enough(sum, b[i], 1e-6));
//         }
//     }

//     #[test]
//     fn test_derivative() {
//         fn f(x: f64) -> f64 {
//             x.sin()
//         }

//         assert!(close_enough(derivative(f, 0f64, 1e-6), 1f64, 1e-6));
//         assert!(close_enough(
//             derivative(f, std::f64::consts::FRAC_PI_2, 1e-6),
//             0f64,
//             1e-6
//         ));
//     }

//     #[test]
//     fn test_partial_derivative() {
//         fn f(x: &[f64; 2]) -> [f64; 2] {
//             [x[0] * x[0] + x[1] * x[1], 2f64 * x[0] * x[1]]
//         }

//         assert!(close_enough(
//             partial_derivative(f, &[0f64, 0f64], 0, 0, 1e-6),
//             0f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[0f64, 0f64], 0, 1, 1e-6),
//             0f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[0f64, 0f64], 1, 0, 1e-6),
//             0f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[0f64, 0f64], 1, 1, 1e-6),
//             0f64,
//             1e-6
//         ));

//         assert!(close_enough(
//             partial_derivative(f, &[1f64, 1f64], 0, 0, 1e-6),
//             2f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[1f64, 1f64], 0, 1, 1e-6),
//             2f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[1f64, 1f64], 1, 0, 1e-6),
//             2f64,
//             1e-6
//         ));
//         assert!(close_enough(
//             partial_derivative(f, &[1f64, 1f64], 1, 1, 1e-6),
//             2f64,
//             1e-6
//         ));
//     }

//     #[test]
//     fn test_solve_newton() {
//         fn f(x: &[f64; 2]) -> [f64; 2] {
//             [
//                 (x[0] + 1f64).sin() - x[1] - 1.2f64,
//                 2f64 * x[0] + x[1].cos() - 2f64,
//             ]
//         }

//         let solution: [f64; 2] = solve_newton(f, &[0f64, 0f64], None).unwrap();

//         assert!(close_enough_arr(&f(&solution), &[0f64, 0f64], 1e-6));
//     }

//     #[test]
//     fn test_solve_newton_vec() {
//         fn f(x: &vec::Vec<f64>) -> vec::Vec<f64> {
//             vec![
//                 (x[0] + 1f64).sin() - x[1] - 1.2f64,
//                 2f64 * x[0] + x[1].cos() - 2f64,
//             ]
//         }

//         let solution: Vec<f64> = solve_newton_vec(&f, &vec![0f64, 0f64], None).unwrap();

//         assert!(zero_enough_vec(&f(&solution), 1e-6));
//     }

//     #[test]
//     fn test_tensor_product() {
//         let a: Vec<Vec<f64>> = vec![
//             vec![0.0, 1.0, 0.0],
//             vec![2.0, 0.0, 4.0],
//             vec![3.0, 3.0, 0.0],
//         ];

//         let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

//         let result: Vec<Vec<f64>> = tensor_product(&a, &b);

//         let reference: Vec<Vec<f64>> = vec![
//             vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
//             vec![2.0, 0.0, 4.0, 4.0, 0.0, 8.0],
//             vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
//             vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
//             vec![3.0, 0.0, 6.0, 4.0, 0.0, 8.0],
//             vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
//         ];

//         println!("{:?}", result);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }

//         let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 3.0, 0.0]];

//         let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

//         let result: Vec<Vec<f64>> = tensor_product(&a, &b);

//         let reference: Vec<Vec<f64>> = vec![
//             vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
//             vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
//             vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
//             vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
//         ];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }
//     }

//     #[test]
//     fn test_matrix_product() {
//         let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 2.0, 0.0]];

//         let b: Vec<Vec<f64>> = vec![vec![1.0, 1.0], vec![1.5, 2.0], vec![2.0, 1.0]];

//         let result: Vec<Vec<f64>> = matrix_product(&a, &b);

//         let reference: Vec<Vec<f64>> = vec![vec![1.5, 2.0], vec![6.0, 7.0]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }

//         let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![1.0, 1.0, 1.0]];

//         let b: Vec<Vec<f64>> = vec![vec![1.0], vec![1.5], vec![2.0]];

//         let result: Vec<Vec<f64>> = matrix_product(&a, &b);

//         let reference: Vec<Vec<f64>> = vec![vec![1.5], vec![4.5]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }
//     }

//     #[test]
//     fn test_matrix_sum_sub() {
//         let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
//         let b: Vec<Vec<f64>> = vec![vec![0.5, 1.5, 2.5], vec![3.5, 4.5, 3.5]];

//         let result: Vec<Vec<f64>> = matrix_sum(&a, &b);
//         let reference: Vec<Vec<f64>> = vec![vec![1.5, 3.5, 5.5], vec![7.5, 9.5, 9.5]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }

//         let result: Vec<Vec<f64>> = matrix_sub(&a, &b);
//         let reference: Vec<Vec<f64>> = vec![vec![0.5, 0.5, 0.5], vec![0.5, 0.5, 2.5]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }
//     }

//     #[test]
//     fn test_matrix_transpose() {
//         let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

//         let result: Vec<Vec<f64>> = matrix_transpose(&a);
//         let reference: Vec<Vec<f64>> = vec![vec![1.0, 4.0], vec![2.0, 5.0], vec![3.0, 6.0]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }
//     }

//     #[test]
//     fn test_scale_transpose() {
//         let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

//         let result: Vec<Vec<f64>> = matrix_scale(&a, 0.5);
//         let reference: Vec<Vec<f64>> = vec![vec![0.5, 1.0, 1.5], vec![2.0, 2.5, 3.0]];

//         println!("{:?}", result);
//         println!("{:?}", reference);

//         assert!(result.len() == reference.len());
//         for i in 0..reference.len() {
//             assert!(result[i].len() == reference[i].len());
//             for j in 0..reference[i].len() {
//                 assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
//             }
//         }
//     }

//     #[test]
//     fn test_c_n_k() {
//         assert!(c_n_k(10, 0) == 1);
//         assert!(c_n_k(10, 1) == 10);
//         assert!(c_n_k(10, 2) == 45);
//         assert!(c_n_k(10, 9) == 10);
//         assert!(c_n_k(10, 10) == 1);
//     }
// }
