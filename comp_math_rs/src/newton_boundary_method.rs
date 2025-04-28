use crate::algebraic_equation_solvers::solve_linear_system;
use crate::boundary_problem::BoundarySolution;
use crate::boundary_problem::NonlinearBoundaryProblem;
use crate::tridiagonal_method::solve_tridiagonal_linear_system;

fn get_first_derivative(y: &std::vec::Vec<[f64; 4]>, n: usize) -> f64 {
    let h: f64 = y[1][0] - y[0][0];
    if n == 0 {
        (y[n + 2][1] - 4.0 * y[n + 1][1] + 3.0 * y[n + 0][1]) / (-2.0 * h)
    } else if n + 1 == y.len() {
        (y[n - 2][1] - 4.0 * y[n - 1][1] + 3.0 * y[n - 0][1]) / (2.0 * h)
    } else {
        (y[n + 1][1] - y[n - 1][1]) / (2.0 * h)
    }
}

fn get_second_derivative(y: &std::vec::Vec<[f64; 4]>, n: usize) -> f64 {
    let h: f64 = y[1][0] - y[0][0];
    if n == 0 {
        (y[n + 2][1] - 2.0 * y[n + 1][1] + y[n + 0][1]) / h.powi(2)
    } else if n + 1 == y.len() {
        (y[n - 2][1] - 2.0 * y[n - 1][1] + y[n][1]) / h.powi(2)
    } else {
        (y[n - 1][1] + y[n + 1][1] - 2.0 * y[n][1]) / h.powi(2)
    }
}

fn calc_derivatives(res: &mut std::vec::Vec<[f64; 4]>) {
    let n = res.len();
    for i in 0..n {
        res[i][2] = get_first_derivative(&res, i);
        res[i][3] = get_second_derivative(&res, i);
    }
}

pub fn newton_boundary_method<F>(
    problem: &mut NonlinearBoundaryProblem<F>,
    step: f64,
    epsilon: f64,
) -> (BoundarySolution, Result<(), &'static str>)
where
    F: FnMut(f64, &[f64; 2]) -> [f64; 2],
{
    let n: usize = ((problem.x_n - problem.x_0) / step).ceil() as usize + 1;
    let step: f64 = (problem.x_n - problem.x_0) / n as f64;
    let mut res: Vec<[f64; 4]> = vec![[0.0; 4]; n]; // [[x, y, y', y'']]
    let h: f64 = step / 4.0;

    // System to find first approximation (p e^x + q)
    let a: [[f64; 2]; 2] = [
        [problem.x_0.exp() * (problem.a_1 + problem.b_1), problem.b_1],
        [problem.x_n.exp() * (problem.a_2 + problem.b_2), problem.b_2],
    ];
    let b: [f64; 2] = [problem.u_1, problem.u_2];
    let [p, q] = solve_linear_system(&a, &b);

    for i in 0..n {
        let x: f64 = problem.x_0 + step * i as f64;
        let y: f64 = p * x.exp() + q;
        res[i][0] = x;
        res[i][1] = y;
    }

    calc_derivatives(&mut res);

    loop {
        let mut a: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
        let mut d: Vec<f64> = vec![0.0; n];

        for i in 1..n - 1 {
            let q: f64 = -((problem.f)(res[i][0], &[res[i][1], res[i][2] + h])[1]
                - (problem.f)(res[i][0], &[res[i][1], res[i][2] - h])[1])
                / (2.0 * h);
            let p: f64 = -((problem.f)(res[i][0], &[res[i][1] + h, res[i][2]])[1]
                - (problem.f)(res[i][0], &[res[i][1] - h, res[i][2]])[1])
                / (2.0 * h);
            let r: f64 = (problem.f)(res[i][0], &[res[i][1], res[i][2]])[1] - res[i][3];

            a[i][i - 1] = 1.0 - q * step / 2.0;
            a[i][i + 1] = 1.0 + q * step / 2.0;
            a[i][i] = step.powi(2) * p - a[i][i - 1] - a[i][i + 1];
            d[i] = r * step.powi(2);
        }
        a[0][0] = problem.b_1 - problem.a_1 / step;
        a[0][1] = problem.a_1 / step;
        a[n - 1][n - 2] = -problem.a_2 / step;
        a[n - 1][n - 1] = -problem.a_2 / step - problem.b_2;
        d[0] = 0.0;
        d[n - 1] = 0.0;
        let a: Vec<Vec<f64>> = a;
        let d: Vec<f64> = d;

        let res_y: Vec<f64> = solve_tridiagonal_linear_system(a, d);

        for i in 0..n {
            res[i][1] += res_y[i];
        }
        calc_derivatives(&mut res);

        let mut discrepancy: f64 = 0.0;
        for i in 0..n {
            discrepancy = if res_y[i].abs() > discrepancy {
                res_y[i].abs()
            } else {
                discrepancy
            };
        }
        if discrepancy < epsilon {
            break;
        }
    }

    let mut res_x: Vec<f64> = vec![0.0; n];
    let mut res_y: Vec<f64> = vec![0.0; n];

    for i in 0..n {
        res_x[i] = res[i][0];
        res_y[i] = res[i][1];
    }

    (
        BoundarySolution {
            x: res_x,
            y: res_y,
            method_name: "Newton quazilinearization method".to_string(),
        },
        Ok(()),
    )
}
