use std::io::Write;

use solvers::*;

mod solvers;

fn f(_t: f64, x: &[f64; 3]) -> [f64; 3] {
    [
        77.27 * (x[1] + x[0] * (1.0 - 8.375e-6 * x[0] - x[1])),
        1.0 / 77.27 * (x[2] - (1.0 + x[0]) * x[1]),
        0.161 * (x[0] - x[2]),
    ]
}

fn write_csv<const N: usize>(
    group: String,
    solution: CauchySolution<N>,
    tau: f64,
    time: std::time::Duration,
) {
    let mut file = std::fs::File::create(format!("../task6_2_data/{}.csv", solution.method_name))
        .expect("Failed to open file");
    file.write(format!("{}\n", group).as_bytes())
        .expect("failed to wrtite to file");
    file.write(format!("{}\n", solution.method_name).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", tau).as_bytes())
        .expect("failed to write to file");
    file.write(format!("{}\n", time.as_secs_f64()).as_bytes())
        .expect("failed to write to file");
    for i in 0..solution.t.len() {
        file.write(format!("{}", solution.t[i]).as_bytes())
            .expect("failed to write to file");
        for j in 0..N {
            file.write(format!(", {}", solution.x[i][j]).as_bytes())
                .expect("failed to write to file");
        }
        file.write("\n".as_bytes())
            .expect("failed to write to file");
    }
}

fn main() {
    // let mut solver: solvers::RungeKuttaMethod<3, 3> = solvers::RungeKuttaMethod::new(
    //     f,
    //     0f64,
    //     800f64,
    //     0.000005f64,
    //     &[0.5f64, 0.5f64, 0.5f64],
    //     3,
    //     Some(10000u32),
    // );

    // solver.set_butcher_table(solvers::ButcherTable::<3>::new(
    //     [
    //         [0f64, 0f64, 0f64],
    //         [0.5f64, 0f64, 0f64],
    //         [-1f64, 2f64, 0f64],
    //     ],
    //     [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
    //     [0f64, 0.5f64, 1f64],
    //     "Kutta's third-order method",
    // ));

    // let mut solver: solvers::RungeKuttaMethod<3, 4> = solvers::RungeKuttaMethod::new(
    //     f,
    //     0f64,
    //     800f64,
    //     0.000005f64,
    //     &[0.5f64, 0.5f64, 0.5f64],
    //     4,
    //     Some(10000u32),
    // );

    // solver.set_butcher_table(solvers::ButcherTable::<4>::new(
    //     [
    //         [0f64, 0f64, 0f64, 0f64],
    //         [0.5f64, 0f64, 0f64, 0f64],
    //         [0f64, 0.5f64, 0f64, 0f64],
    //         [0f64, 0f64, 1f64, 0f64],
    //     ],
    //     [1f64 / 6f64, 1f64 / 3f64, 1f64 / 3f64, 1f64 / 6f64],
    //     [0f64, 1f64 / 2f64, 1f64 / 2f64, 1f64],
    //     "Classic fourth-order method",
    // ));

    let problem: CauchyProblem<3> = CauchyProblem {
        f,
        start: 0.0,
        stop: 800.0,
        x_0: [0.5, 0.5, 0.5],
    };

    let mut solver: solvers::RungeKuttaMethod<3, 3, 9> = solvers::RungeKuttaMethod::new(
        4,
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-1.0, 2.0, 0.0]],
        [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
        [0f64, 0.5f64, 1f64],
        "Kutta's third-order method".to_string(),
    );
    let tau: f64 = 0.00001;
    let save_every: u32 = (0.01f64 / tau).round() as u32;
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, res) = solver.solve(&problem, tau, false, Some(save_every));
    println!("{:?}", res);
    write_csv("Runge-Kutta".to_string(), solution, tau, start_time.elapsed());

    let mut solver: solvers::RungeKuttaMethod<3, 2, 6> = solvers::RungeKuttaMethod::new(
        2,
        [[0.0, 0.0], [0.5, 0.5]],
        [0.0f64, 1.0f64],
        [0.5f64, 0.5f64],
        "Crank-Nikolson method".to_string(),
    );
    let tau: f64 = 0.00001;
    let save_every: u32 = (0.01f64 / tau).round() as u32;
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, res) = solver.solve(&problem, tau, true, Some(save_every));
    println!("{:?}", res);
    write_csv("Runge-Kutta".to_string(), solution, tau, start_time.elapsed());

    // for order in 1..5 {
    //     let mut solver: solvers::AdamsMethod<3> =
    //         solvers::AdamsMethod::new(order, solvers::SolverType::Explicit);
    //     let tau: f64 = 0.000001;
    //     let save_every: u32 = (0.01f64 / tau).round() as u32;
    //     let start_time: std::time::Instant = std::time::Instant::now();
    //     let (solution, res) = solver.solve(&problem, tau, false, Some(save_every));
    //     println!("{:?}", res);
    //     write_csv("Explicit Adams".to_string(), solution, tau, start_time.elapsed());
    // }
    
    for order in 1..5 {
        let mut solver: solvers::AdamsMethod<3> =
            solvers::AdamsMethod::new(order, solvers::SolverType::Implicit);
        let tau: f64 = 0.01;
        let save_every: u32 = (0.01f64 / tau).round() as u32;
        let start_time: std::time::Instant = std::time::Instant::now();
        let (solution, res) = solver.solve(&problem, tau, false, Some(save_every));
        println!("{:?}", res);
        write_csv("Implicit Adams".to_string(), solution, tau, start_time.elapsed());
    }
}
