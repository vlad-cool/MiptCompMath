use std::io::Write;

use solvers::*;

mod solvers;

// fn f(t: f64, x: &[f64; 2]) -> [f64; 2] {
//     [x[1], std::f64::consts::E * (1f64 - x[0] * x[0]) * x[1] - x[0]]
// }

fn f(_t: f64, x: &[f64; 3]) -> [f64; 3] {
    [
        77.27 * (x[1] + x[0] * (1.0 - 8.375e-6 * x[0] - x[1])),
        1.0 / 77.27 * (x[2] - (1.0 + x[0]) * x[1]),
        0.161 * (x[0] - x[2]),
    ]
    // [x[1], std::f64::consts::E * (1f64 - x[0] * x[0]) * x[1] - x[0]]
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
        f: f,
        start: 0.0,
        stop: 800.0,
        x_0: [0.5, 0.5, 0.5],
    };

    let mut solver: solvers::RungeKuttaMethod<3, 3> = solvers::RungeKuttaMethod::new(
        3,
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-1.0, 2.0, 0.0]],
        [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
        [0f64, 0.5f64, 1f64],
        "Kutta's third-order method".to_string(),
    );

    let (solution, res) = solver.solve(problem, 0.00001, true, Some(10000));

    println!("{:?}", res);

    let mut file = std::fs::File::create("test.csv").unwrap();
    for i in 0..solution.t.len() {
        file.write(
            format!(
                "{}, {}, {}, {}\n",
                solution.t[i], solution.x[i][0], solution.x[i][1], solution.x[i][2]
            )
            .as_bytes(),
        )
        .expect("failed to write to file");
    }
}
