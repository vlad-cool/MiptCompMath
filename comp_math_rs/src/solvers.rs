#![allow(dead_code)]
use crate::cauchy_problem::*;
use std::{f64, io::Write, vec};


fn tensor_product(
    l: &vec::Vec<vec::Vec<f64>>,
    r: &vec::Vec<vec::Vec<f64>>,
) -> vec::Vec<vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len() * r[0].len()]; l.len() * r.len()];

    for a in 0..r.len() {
        for b in 0..r[0].len() {
            for c in 0..l.len() {
                for d in 0..l[0].len() {
                    res[a * l.len() + c][b * l[0].len() + d] = l[c][d] * r[a][b];
                }
            }
        }
    }

    res
}

fn matrix_product(
    l: &vec::Vec<vec::Vec<f64>>,
    r: &vec::Vec<vec::Vec<f64>>,
) -> vec::Vec<vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    assert!(l[0].len() == r.len());

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; r[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..r[0].len() {
            for c in 0..r.len() {
                res[a][b] += r[c][b] * l[a][c];
            }
        }
    }

    res
}

fn matrix_sum(l: &vec::Vec<vec::Vec<f64>>, r: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    assert!(l.len() == r.len());
    assert!(l[0].len() == r[0].len());
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..l[0].len() {
            res[a][b] = l[a][b] + r[a][b];
        }
    }

    res
}

fn matrix_sub(l: &vec::Vec<vec::Vec<f64>>, r: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
    assert!(l.len() > 0);
    assert!(r.len() > 0);
    assert!(l.len() == r.len());
    assert!(l[0].len() == r[0].len());
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }
    for i in 0..(r.len() - 1) {
        assert!(r[i].len() == r[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for a in 0..l.len() {
        for b in 0..l[0].len() {
            res[a][b] = l[a][b] - r[a][b];
        }
    }

    res
}

fn matrix_transpose(a: &vec::Vec<vec::Vec<f64>>) -> vec::Vec<vec::Vec<f64>> {
    for i in 0..(a.len() - 1) {
        assert!(a[i].len() == a[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; a.len()]; a[0].len()];

    for i in 0..a.len() {
        for j in 0..a[0].len() {
            res[j][i] = a[i][j];
        }
    }

    res
}

fn matrix_scale(l: &vec::Vec<vec::Vec<f64>>, r: f64) -> vec::Vec<vec::Vec<f64>> {
    for i in 0..(l.len() - 1) {
        assert!(l[i].len() == l[i + 1].len());
    }

    let mut res: Vec<Vec<f64>> = vec![vec![0.0; l[0].len()]; l.len()];

    for i in 0..l.len() {
        for j in 0..l[0].len() {
            res[i][j] = l[i][j] * r;
        }
    }

    res
}

pub struct RosenbrockMethod<const N: usize, const M: usize, const MN: usize> {
    order: u32,
    solver_type: SolverType,
    a: [[f64; M]; M],
    b: [f64; M],
    c: [f64; M],
    name: String,
}

impl<const N: usize, const M: usize, const NM: usize> RosenbrockMethod<N, M, NM> {
    pub fn new(order: u32, a: [[f64; M]; M], b: [f64; M], c: [f64; M], name: String) -> Self {
        todo!()
        // let mut solver_type = SolverType::Explicit;
        // for i in 0..M {
        //     for j in i..M {
        //         if a[i][j] != 0f64 {
        //             solver_type = SolverType::Implicit;
        //         }
        //     }
        // }

        // Self {
        //     order,
        //     solver_type,
        //     a,
        //     b,
        //     c,
        //     name,
        // }
    }

    fn step(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        t: f64,
        x: [f64; N],
    ) -> Result<(f64, [f64; N]), &'static str> {
        todo!()
        // let equation = |k_0: &[f64; NM]| {
        //     let k_0: [[f64; N]; M] = unflatten::<N, M>(k_0);
        //     let mut k_i: [[f64; N]; M] = [[0.0; N]; M];

        //     for i in 0..M {
        //         let arg_1: f64 = t + tau * self.c[i];
        //         let mut arg_2: [f64; N] = [0f64; N];

        //         for a in 0..N {
        //             for j in 0..M {
        //                 arg_2[a] += self.a[i][j] * k_0[j][a];
        //             }
        //             arg_2[a] *= tau;
        //             arg_2[a] += x[a];
        //         }
        //         for a in 0..N {
        //             k_i[i][a] = (problem.f)(arg_1, &arg_2)[a] - k_0[i][a];
        //         }
        //     }

        //     let k_i: [f64; NM] = (*k_i.as_flattened()).try_into().unwrap();

        //     k_i
        // };

        // let k = match solve_newton(equation, &[0f64; NM], None) {
        //     Ok(x) => unflatten::<N, M>(&x),
        //     Err(err) => {
        //         println!("Failed to solve, {}", err);
        //         return Err("Failed to solve");
        //     }
        // };

        // let mut res: [f64; N] = [0f64; N];
        // for j in 0..N {
        //     for i in 0..M {
        //         res[j] += self.b[i] * k[i][j];
        //     }
        //     res[j] *= tau;
        //     res[j] += x[j];
        // }

        // Ok((t + tau, res))
    }
}

impl<const N: usize, const M: usize, const NM: usize> CauchySolver<N>
    for RosenbrockMethod<N, M, NM>
{
    fn solve(
        &mut self,
        problem: &CauchyProblem<N>,
        tau: f64,
        print_progress: bool,
        save_every: Option<u32>,
    ) -> (CauchySolution<N>, Result<(), &'static str>) {
        todo!()
        // let mut t_i: f64 = problem.start;
        // let mut x_i: [f64; N] = problem.x_0.clone();
        // let mut iterations: u32 = 0u32;
        // let save_every = save_every.unwrap_or(1);

        // let mut solution: CauchySolution<N> = CauchySolution {
        //     t: vec![],
        //     x: vec![],
        //     method_name: self.get_name(),
        // };

        // solution.t.push(problem.start);
        // solution.x.push(problem.x_0.clone());
        // if print_progress {
        //     print!("t: {:>20.10}, iterations: {}", problem.start, iterations)
        // }
        // iterations += 1;

        // while t_i < problem.stop {
        //     match self.step(&problem, tau, t_i, x_i) {
        //         Ok((t, x)) => {
        //             t_i = t;
        //             x_i = x.clone();
        //             if iterations % save_every == 0 {
        //                 solution.t.push(t);
        //                 solution.x.push(x);
        //                 if print_progress {
        //                     print!("\rt: {:>10.6}, iterations: {:>16}", t, iterations)
        //                 }
        //             }
        //         }
        //         Err(a) => {
        //             if print_progress {
        //                 println!("\n\n");
        //             }
        //             println!("Failed to solve, reason: {}", a);
        //             return (solution, Err("Failed to solve"));
        //         }
        //     };
        //     iterations += 1;
        // }

        // if print_progress {
        //     println!("");
        // }

        // (solution, Ok(()))
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

/// Tests

#[cfg(test)]
mod tests {
    use crate::solvers::*;

    #[test]
    fn test_solve_linear_system() {
        const N: usize = 3;
        let a: [[f64; 3]; 3] = [
            [-4f64, 9f64, -4f64],
            [-5f64, -5f64, 6f64],
            [2f64, 5f64, -8f64],
        ];

        let b: [f64; 3] = [-64f64, 104.6f64, -85.2f64];

        let x: [f64; 3] = crate::solvers::solve_linear_system(&a, &b);
        for i in 0..N {
            let mut sum: f64 = 0f64;
            for j in 0..N {
                sum += a[i][j] * x[j];
            }
            assert!(close_enough(sum, b[i], 1e-6));
        }
    }

    #[test]
    fn test_derivative() {
        fn f(x: f64) -> f64 {
            x.sin()
        }

        assert!(close_enough(derivative(f, 0f64, 1e-6), 1f64, 1e-6));
        assert!(close_enough(
            derivative(f, std::f64::consts::FRAC_PI_2, 1e-6),
            0f64,
            1e-6
        ));
    }

    #[test]
    fn test_partial_derivative() {
        fn f(x: &[f64; 2]) -> [f64; 2] {
            [x[0] * x[0] + x[1] * x[1], 2f64 * x[0] * x[1]]
        }

        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 0, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 0, 1, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 1, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 1, 1, 1e-6),
            0f64,
            1e-6
        ));

        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 0, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 0, 1, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 1, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 1, 1, 1e-6),
            2f64,
            1e-6
        ));
    }

    #[test]
    fn test_solve_newton() {
        fn f(x: &[f64; 2]) -> [f64; 2] {
            [
                (x[0] + 1f64).sin() - x[1] - 1.2f64,
                2f64 * x[0] + x[1].cos() - 2f64,
            ]
        }

        let solution: [f64; 2] = solve_newton(f, &[0f64, 0f64], None).unwrap();

        assert!(close_enough_arr(&f(&solution), &[0f64, 0f64], 1e-6));
    }


    #[test]
    fn test_tensor_product() {
        let a: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0],
            vec![2.0, 0.0, 4.0],
            vec![3.0, 3.0, 0.0],
        ];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

        let result: Vec<Vec<f64>> = tensor_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
            vec![2.0, 0.0, 4.0, 4.0, 0.0, 8.0],
            vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
            vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
            vec![3.0, 0.0, 6.0, 4.0, 0.0, 8.0],
            vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
        ];

        println!("{:?}", result);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 3.0, 0.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 2.0], vec![1.5, 2.0]];

        let result: Vec<Vec<f64>> = tensor_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![
            vec![0.0, 1.0, 0.0, 0.0, 2.0, 0.0],
            vec![3.0, 3.0, 0.0, 6.0, 6.0, 0.0],
            vec![0.0, 1.5, 0.0, 0.0, 2.0, 0.0],
            vec![4.5, 4.5, 0.0, 6.0, 6.0, 0.0],
        ];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_product() {
        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![3.0, 2.0, 0.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0, 1.0], vec![1.5, 2.0], vec![2.0, 1.0]];

        let result: Vec<Vec<f64>> = matrix_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![vec![1.5, 2.0], vec![6.0, 7.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let a: Vec<Vec<f64>> = vec![vec![0.0, 1.0, 0.0], vec![1.0, 1.0, 1.0]];

        let b: Vec<Vec<f64>> = vec![vec![1.0], vec![1.5], vec![2.0]];

        let result: Vec<Vec<f64>> = matrix_product(&a, &b);

        let reference: Vec<Vec<f64>> = vec![vec![1.5], vec![4.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_sum_sub() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let b: Vec<Vec<f64>> = vec![vec![0.5, 1.5, 2.5], vec![3.5, 4.5, 3.5]];

        let result: Vec<Vec<f64>> = matrix_sum(&a, &b);
        let reference: Vec<Vec<f64>> = vec![vec![1.5, 3.5, 5.5], vec![7.5, 9.5, 9.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }

        let result: Vec<Vec<f64>> = matrix_sub(&a, &b);
        let reference: Vec<Vec<f64>> = vec![vec![0.5, 0.5, 0.5], vec![0.5, 0.5, 2.5]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_matrix_transpose() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let result: Vec<Vec<f64>> = matrix_transpose(&a);
        let reference: Vec<Vec<f64>> = vec![vec![1.0, 4.0], vec![2.0, 5.0], vec![3.0, 6.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

    #[test]
    fn test_scale_transpose() {
        let a: Vec<Vec<f64>> = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];

        let result: Vec<Vec<f64>> = matrix_scale(&a, 0.5);
        let reference: Vec<Vec<f64>> = vec![vec![0.5, 1.0, 1.5], vec![2.0, 2.5, 3.0]];

        println!("{:?}", result);
        println!("{:?}", reference);

        assert!(result.len() == reference.len());
        for i in 0..reference.len() {
            assert!(result[i].len() == reference[i].len());
            for j in 0..reference[i].len() {
                assert!((result[i][j] - reference[i][j]).abs() < 1e-5);
            }
        }
    }

}
