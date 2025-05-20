#![allow(dead_code)]

use std::{io::Write, vec};

fn write_csv(path: String, u: Vec<Vec<f64>>, h_x: f64, h_y: f64, time: std::time::Duration) {
    let mut file =
        std::fs::File::create(format!("../task8_16_data/{path}.csv")).expect("Failed to open file");
    file.write(format!("{}\n", h_x).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", h_y).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
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

fn solve_jacobi(mut u: Vec<Vec<f64>>, h_x: f64, h_y: f64, epsilon: f64) -> Vec<Vec<f64>> {
    let mut u_new: &mut Vec<Vec<f64>> = &mut u.clone();
    let mut u: &mut Vec<Vec<f64>> = &mut u;

    let mut max_delta: f64 = epsilon;
    while max_delta >= epsilon {
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
        println!("{}", max_delta);
    }
    (*u).clone()
}

fn solve_seidel(mut u: Vec<Vec<f64>>, h_x: f64, h_y: f64, epsilon: f64) -> Vec<Vec<f64>> {
    let mut max_delta: f64 = epsilon;
    while max_delta >= epsilon {
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
        println!("{}", max_delta);
    }
    u
}

fn main() {
    let l: f64 = 1.0;
    let h_x: f64 = 0.0005;
    let h_y: f64 = 0.0005;
    // let h_x: f64 = 1.0;
    // let h_y: f64 = 1.0;
    let x_n: usize = (l / h_x) as usize;
    let y_n: usize = (l / h_y) as usize;

    let t_ab: f64 = 1.0;
    let t_bc: f64 = 2.0;
    let t_cd: f64 = 3.0;
    let t_da: f64 = 4.0;

    let mut u: Vec<Vec<f64>> = vec![vec![2.5; x_n]; y_n];

    for y_i in 0..y_n {
        u[y_i][0] = t_ab;
        u[y_i][x_n - 1] = t_cd;
    }

    for x_i in 0..x_n {
        u[0][x_i] = t_bc;
        u[y_n - 1][x_i] = t_da;
    }

    // println!("{:?}", u);

    let start_time: std::time::Instant = std::time::Instant::now();
    // let u: Vec<Vec<f64>> = solve_seidel(u, h_x, h_y, 1e-3);
    let u: Vec<Vec<f64>> = solve_jacobi(u, h_x, h_y, 1e-3);
    write_csv(format!("test"), u, h_x, h_y, start_time.elapsed());
}
