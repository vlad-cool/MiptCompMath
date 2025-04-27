use std::vec;

use crate::algebraic_equation_solvers::solve_newton;
use crate::boundary_problem::LinearBoundaryProblem;
use crate::boundary_problem::LinearBoundarySolution;
use crate::cauchy_problem::CauchyProblem;
use crate::cauchy_problem::CauchySolver;

fn solve_tridiagonal_linear_system(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n: usize = b.len();
    assert!(n > 0);
    assert!(a.len() == n);
    for i in 0..n {
        assert!(a[i].len() == n);
    }

    for i in 1..n {
        let c: f64 = -a[i][i - 1] / a[i - 1][i - 1];
        let end: usize = if i + 1 < n { i + 2 } else { n };
        for j in i - 1..end {
            a[i][j] += a[i - 1][j] * c;
        }
        b[i] += b[i - 1] * c;
    }

    println!("{:?}\n\n---\n", a);

    for i in (0..n - 1).rev() {
        let c: f64 = -a[i][i + 1] / a[i + 1][i + 1];
        let start: usize = if i == 0 { 0 } else { i - 1 };
        for j in start..i + 2 {
            a[i][j] += a[i + 1][j] * c;
        }
        b[i] += b[i + 1] * c;
    }

    println!("{:?}", a);

    for i in 0..n {
        b[i] /= a[i][i];
    }

    b
}

pub fn tridiagonal_method(problem: &LinearBoundaryProblem, step: f64) -> LinearBoundarySolution {
    let n: usize = ((problem.x_n - problem.x_0) / step).ceil() as usize + 1;
    let step: f64 = (problem.x_n - problem.x_0) / step as f64;
    let mut res_x = vec![0.0; n];

    let mut a: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut d: Vec<f64> = vec![0.0; n];
    for i in 1..n - 1 {
        let x: f64 = problem.x_0 + step * i as f64;
        res_x[i] = x;
        a[i][i - 1] = 1.0 - (problem.q)(x) * step / 2.0;
        a[i][i + 1] = 1.0 + (problem.q)(x) * step / 2.0;
        a[i][i] = step.powi(2) * (problem.p)(x) - a[i][i - 1] + a[i][i + 1];
        d[i] = (problem.f)(x) * step.powi(2);
    }
    a[0][0] = problem.b_1 - problem.a_1 / step;
    a[0][1] = problem.a_1 / step;
    a[n - 1][n - 2] = -problem.a_2 / step;
    a[n - 1][n - 1] = -problem.a_2 / step - problem.b_2;
    d[0] = problem.u_1;
    d[n - 1] = problem.u_2;
    let a: Vec<Vec<f64>> = a;
    let d: Vec<f64> = d;

    let res_y: Vec<f64> = solve_tridiagonal_linear_system(a, d);

    LinearBoundarySolution {
        x: res_x,
        y: res_y,
        method_name: "Tridiagonal method".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::vec;

    use crate::{tridiagonal_method::*, utils::close_enough};

    #[test]
    fn test_linear() {
        let a = vec![
            vec![1.0, 9.0, 0.0, 0.0],
            vec![0.0, 3.0, 3.0, 0.0],
            vec![0.0, 4.0, 5.0, 7.0],
            vec![0.0, 0.0, 2.0, 6.0],
        ];
        let b: Vec<f64> = vec![5.0, 3.0, 4.0, 1.0];

        let mut x: Vec<f64> = solve_tridiagonal_linear_system(a.clone(), b.clone());

        println!("{:?}", x);

        for i in 0..4 {
            let mut sum: f64 = 0.0;
            for j in 0..4 {
                sum += a[i][j] * x[j];
            }
            assert!(close_enough(sum - b[i], 0.0, 1e-5));
        }
    }
}
