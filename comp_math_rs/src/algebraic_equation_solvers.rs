use crate::utils;

pub fn solve_linear_system<const N: usize>(a: &[[f64; N]; N], b: &[f64; N]) -> [f64; N] {
    let mut l: [[f64; N]; N] = [[0f64; N]; N];
    let mut u: [[f64; N]; N] = [[0f64; N]; N];

    for i in 0..N {
        l[i][i] = 1f64;
        for j in i..N {
            let mut sum: f64 = 0f64;
            for k in 0..i {
                sum += l[i][k] * u[k][j];
            }
            u[i][j] = a[i][j] - sum;
        }

        for j in i..N {
            let mut sum: f64 = 0f64;
            for k in 0..i {
                sum += l[j][k] * u[k][i];
            }
            l[j][i] = (a[j][i] - sum) / u[i][i];
        }
    }

    let mut v: [f64; N] = [0f64; N];
    for i in 0..N {
        v[i] = b[i];
        for j in 0..i {
            v[i] -= l[i][j] * v[j];
        }
    }

    let mut x: [f64; N] = [0f64; N];

    for i in (0..N).rev() {
        x[i] = v[i] / u[i][i];
        for j in i + 1..N {
            x[i] -= u[i][j] * x[j] / u[i][i];
        }
    }

    x
}

fn get_discrepancy_vec(
    a: &std::vec::Vec<std::vec::Vec<f64>>,
    b: &std::vec::Vec<f64>,
    x: &std::vec::Vec<f64>,
) -> f64 {
    let n = a.len();
    let mut disc: f64 = 0.0;
    for i in 0..n {
        let mut sum: f64 = 0f64;
        for j in 0..n {
            sum += a[i][j] * x[j];
        }
        disc = if disc > (sum - b[i]).abs() {
            disc
        } else {
            (sum - b[i]).abs()
        };
    }
    disc
}

pub fn solve_linear_system_vec(
    a: &std::vec::Vec<std::vec::Vec<f64>>,
    b: &std::vec::Vec<f64>,
) -> std::vec::Vec<f64> {
    let n = b.len();

    let mut a_new: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    let mut b_new: Vec<f64> = vec![0.0; n];

    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                a_new[i][j] += a[i][k] * a[k][j];
            }
        }
    }

    for i in 0..n {
        for j in 0..n {
            b_new[i] += b[j] * a[i][j];
        }
    }

    let a: Vec<Vec<f64>> = a_new;
    let b: Vec<f64> = b_new;

    let mut x: Vec<f64> = vec![0f64; n];

    while get_discrepancy_vec(&a, &b, &x) > 1e-7 {
        for i in 0..n {
            x[i] = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                x[i] -= a[i][j] * x[j];
            }

            x[i] += b[i];
            x[i] /= a[i][i];
        }
    }
    x
}

pub fn solve_newton<F, const N: usize>(
    mut f: F,
    x: &[f64; N],
    max_iterations: Option<u32>,
) -> Result<[f64; N], &'static str>
where
    F: FnMut(&[f64; N]) -> [f64; N],
{
    let max_iterations: u32 = max_iterations.unwrap_or(100u32);

    let mut x: [f64; N] = x.clone();
    let mut iterations: u32 = 0u32;
    let mut discrepancy: [f64; N] = f(&x);

    while !utils::close_enough_arr(&discrepancy, &[0f64; N], 1e-6) {
        if iterations == max_iterations {
            return Err("Cannot solve system, too many iterations");
        }

        let mut jacobian: [[f64; N]; N] = [[0f64; N]; N];

        for i in 0..N {
            for j in 0..N {
                jacobian[i][j] = -utils::partial_derivative(&mut f, &x, i, j, 1e-5);
            }
        }

        let delta_x: [f64; N] = solve_linear_system(&jacobian, &discrepancy);

        for i in 0..N {
            x[i] += delta_x[i];
        }

        iterations += 1;
        discrepancy = f(&x);
    }

    Ok(x)
}

pub fn solve_newton_vec<F>(
    mut f: F,
    x: &std::vec::Vec<f64>,
    max_iterations: Option<u32>,
) -> Result<std::vec::Vec<f64>, &'static str>
where
    F: FnMut(&std::vec::Vec<f64>) -> std::vec::Vec<f64>,
{
    let max_iterations: u32 = max_iterations.unwrap_or(100u32);
    let n = x.len();

    let mut x: Vec<f64> = x.clone();
    let mut iterations: u32 = 0u32;

    while !utils::zero_enough_vec(&f(&x), 1e-6) {
        if iterations == max_iterations {
            return Err("Cannot solve system, too many iterations");
        }

        let mut jacobian: Vec<Vec<f64>> = vec![vec![0.0; n]; n];

        for i in 0..n {
            for j in 0..n {
                // jacobian[j][i] = -partial_derivative_vec(f, &x, i, j, 1e-5);
                jacobian[i][j] = -utils::partial_derivative_vec(&mut f, &x, i, j, 1e-5);
            }
        }

        let delta_x: Vec<f64> = solve_linear_system_vec(&jacobian, &f(&x));
        // println!("jac: {:?}", jacobian);
        // println!("f: {:?}", f(&x));
        // println!("{:?}", delta_x);

        // println!("{:?}", jacobian);
        // println!("{:?}", f(&x));
        // println!("{:?}", delta_x);

        for i in 0..n {
            x[i] += delta_x[i];
        }

        iterations += 1;
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use crate::algebraic_equation_solvers::*;
    use crate::utils::*;

    #[test]
    fn test_solve_linear_system() {
        const N: usize = 3;
        let a: [[f64; 3]; 3] = [
            [-4f64, 9f64, -4f64],
            [-5f64, -5f64, 6f64],
            [2f64, 5f64, -8f64],
        ];

        let b: [f64; 3] = [-64f64, 104.6f64, -85.2f64];

        let x: [f64; 3] = solve_linear_system(&a, &b);
        for i in 0..N {
            let mut sum: f64 = 0f64;
            for j in 0..N {
                sum += a[i][j] * x[j];
            }
            assert!(close_enough(sum, b[i], 1e-6));
        }
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
    fn test_solve_linear_system_vec() {
        const N: usize = 3;
        // let a: Vec<Vec<f64>> = vec![
        //     vec![-4f64, 9f64, -4f64],
        //     vec![-5f64, -5f64, 6f64],
        //     vec![2f64, 5f64, -8f64],
        // ];
        // let b: Vec<f64> = vec![1.0, 3.0, 2.0];
        let a: Vec<Vec<f64>> = vec![
            vec![1.0, 0.0, 3.0],
            vec![1.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.5],
        ];
        let b: Vec<f64> = vec![1.0, 1.0, 1.0];

        let x: Vec<f64> = solve_linear_system_vec(&a, &b);

        println!("{:?}", x);
        assert!(get_discrepancy_vec(&a, &b, &x) < 1e-6);
    }

    #[test]
    fn test_solve_newton_vec() {
        fn f(x: &std::vec::Vec<f64>) -> std::vec::Vec<f64> {
            vec![
                (x[0] + 1f64).sin() - x[1] - 1.2f64,
                2f64 * x[0] + x[1].cos() - 2f64,
            ]
        }

        let solution: Vec<f64> = solve_newton_vec(&f, &vec![0f64, 0f64], None).unwrap();

        assert!(zero_enough_vec(&f(&solution), 1e-6));
    }
}
