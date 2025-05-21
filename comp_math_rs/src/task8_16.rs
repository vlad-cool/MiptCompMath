#![allow(dead_code)]

use std::{io::Write, vec};

fn write_csv(
    path: String,
    u: Vec<Vec<f64>>,
    h_x: f64,
    h_y: f64,
    time: std::time::Duration,
    iterations: u32,
) {
    let mut file =
        std::fs::File::create(format!("../task8_16_data/{path}.csv")).expect("Failed to open file");
    file.write(format!("{}\n", h_x).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", h_y).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", iterations).as_bytes())
        .expect("failed to write to file");

    let x_step: usize = 1;
    let y_step: usize = 1;

    for y_i in (0..u.len()).step_by(y_step) {
        for x_i in (0..u[y_i].len()).step_by(x_step) {
            file.write(format!("{:.8}, ", u[y_i][x_i]).as_bytes())
                .expect("failed to write to file");
        }
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn solve_jacobi(mut u: Vec<Vec<f64>>, h_x: f64, h_y: f64, epsilon: f64) -> (u32, Vec<Vec<f64>>) {
    let mut u_new: &mut Vec<Vec<f64>> = &mut u.clone();
    let mut u: &mut Vec<Vec<f64>> = &mut u;
    let mut iterations: u32 = 0;

    let mut max_delta: f64 = epsilon;
    while max_delta >= epsilon {
        iterations += 1;
        max_delta = 0.0;
        for y_i in 1..(u.len() - 1) {
            for x_i in 1..(u[y_i].len() - 1) {
                u_new[y_i][x_i] = 0.0;
                u_new[y_i][x_i] += (u[y_i][x_i + 1] + u[y_i][x_i - 1]) * h_x.powi(2);
                u_new[y_i][x_i] += (u[y_i + 1][x_i] + u[y_i - 1][x_i]) * h_y.powi(2);
                u_new[y_i][x_i] /= (h_x.powi(2) + h_y.powi(2)) * 2.0;
                if (u[y_i][x_i] - u_new[y_i][x_i]).abs() > max_delta {
                    max_delta = (u[y_i][x_i] - u_new[y_i][x_i]).abs()
                }
            }
        }

        (u_new, u) = (u, u_new);
        // println!("{}", max_delta);
    }
    (iterations, (*u).clone())
}

fn solve_seidel(mut u: Vec<Vec<f64>>, h_x: f64, h_y: f64, epsilon: f64) -> (u32, Vec<Vec<f64>>) {
    let mut max_delta: f64 = epsilon;
    let mut iterations: u32 = 0;

    while max_delta >= epsilon {
        iterations += 1;
        max_delta = 0.0;
        for y_i in 1..(u.len() - 1) {
            for x_i in 1..(u[y_i].len() - 1) {
                let mut new_u: f64 = 0.0;
                new_u += (u[y_i][x_i + 1] + u[y_i][x_i - 1]) * h_x.powi(2);
                new_u += (u[y_i + 1][x_i] + u[y_i - 1][x_i]) * h_y.powi(2);
                new_u /= (h_x.powi(2) + h_y.powi(2)) * 2.0;
                if (u[y_i][x_i] - new_u).abs() > max_delta {
                    max_delta = (u[y_i][x_i] - new_u).abs()
                }
                u[y_i][x_i] = new_u;
            }
        }
        // for y_i in (1..(u.len() - 1)).rev() {
        //     for x_i in (1..(u[y_i].len() - 1)).rev() {
        //         let mut new_u: f64 = 0.0;
        //         new_u += (u[y_i][x_i + 1] + u[y_i][x_i - 1]) * h_x.powi(2);
        //         new_u += (u[y_i + 1][x_i] + u[y_i - 1][x_i]) * h_y.powi(2);
        //         new_u /= (h_x.powi(2) + h_y.powi(2)) * 2.0;
        //         if (u[y_i][x_i] - new_u).abs() > max_delta {
        //             max_delta = (u[y_i][x_i] - new_u).abs()
        //         }
        //         u[y_i][x_i] = new_u;
        //     }
        // }
        // println!("{}", max_delta);
    }
    (iterations, u)
}

fn solve_over_relaxation(
    mut u: Vec<Vec<f64>>,
    h_x: f64,
    h_y: f64,
    epsilon: f64,
    tau: f64,
) -> (u32, Vec<Vec<f64>>) {
    let mut max_delta: f64 = epsilon;
    let mut iterations: u32 = 0;

    while max_delta >= epsilon {
        iterations += 1;
        max_delta = 0.0;
        for y_i in 1..(u.len() - 1) {
            for x_i in 1..(u[y_i].len() - 1) {
                let mut new_u: f64 = 0.0;
                new_u += (u[y_i][x_i + 1] + u[y_i][x_i - 1]) * h_x.powi(2);
                new_u += (u[y_i + 1][x_i] + u[y_i - 1][x_i]) * h_y.powi(2);
                new_u /= (h_x.powi(2) + h_y.powi(2)) * 2.0;
                if (u[y_i][x_i] - new_u).abs() > max_delta {
                    max_delta = (u[y_i][x_i] - new_u).abs()
                }
                u[y_i][x_i] = new_u * tau + u[y_i][x_i] * (1.0 - tau);
            }
        }
        // println!("{}", max_delta);
    }
    (iterations, u)
}

fn main() {
    let l: f64 = 0.1;

    let t_ab: f64 = 1.0;
    let t_bc: f64 = 2.0;
    let t_cd: f64 = 3.0;
    let t_da: f64 = 4.0;

    for x_h_pow in 1..4 {
        for y_h_pow in 1..4 {
            let h_x: f64 = 0.05 * 0.1f64.powi(x_h_pow);
            let h_y: f64 = 0.05 * 0.1f64.powi(y_h_pow);

            let x_n: usize = (l / h_x) as usize;
            let y_n: usize = (l / h_y) as usize;

            let mut u_0: Vec<Vec<f64>> = vec![vec![2.5; x_n]; y_n];

            for y_i in 0..y_n {
                u_0[y_i][0] = t_ab;
                u_0[y_i][x_n - 1] = t_cd;
            }

            for x_i in 0..x_n {
                u_0[0][x_i] = t_bc;
                u_0[y_n - 1][x_i] = t_da;
            }
            
            //////////////////////////////////////////////

            let u: Vec<Vec<f64>> = u_0.clone();
            let start_time: std::time::Instant = std::time::Instant::now();
            let (iterations, u) = solve_jacobi(u, h_x, h_y, 1e-3);
            write_csv(
                format!("jacobi_{}_{}", x_h_pow, y_h_pow),
                u,
                h_x,
                h_y,
                start_time.elapsed(),
                iterations,
            );
            
            //////////////////////////////////////////////
            
            let u: Vec<Vec<f64>> = u_0.clone();
            let start_time: std::time::Instant = std::time::Instant::now();
            let (iterations, u) = solve_seidel(u, h_x, h_y, 1e-3);
            write_csv(
                format!("seidel_{}_{}", x_h_pow, y_h_pow),
                u,
                h_x,
                h_y,
                start_time.elapsed(),
                iterations,
            );
            
            //////////////////////////////////////////////
            
            for i in 1..5 {
                let tau: f64 = 1.0 + 0.2 * (i as f64);
                
                let u: Vec<Vec<f64>> = u_0.clone();
                let start_time: std::time::Instant = std::time::Instant::now();
                let (iterations, u) = solve_over_relaxation(u, h_x, h_y, 1e-3, tau);
                write_csv(
                    format!("or_1.{}_{}_{}", i * 2, x_h_pow, y_h_pow),
                    u,
                    h_x,
                    h_y,
                    start_time.elapsed(),
                    iterations,
                );
            }
        }
    }
}
