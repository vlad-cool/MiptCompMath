use std::io::Write;

use solvers::DifferentialEquationNumericMethod;

mod solvers;

// fn f(t: f64, x: &[f64; 2]) -> [f64; 2] {
//     [x[1], std::f64::consts::E * (1f64 - x[0] * x[0]) * x[1] - x[0]]
// }

fn f(t: f64, x: &[f64; 3]) -> [f64; 3] {
    [
        77.27 * (x[1] + x[0] * (1.0 - 8.375e-6 * x[0] - x[1])),
        1.0 / 77.27 * (x[2] - (1.0 + x[0]) * x[1]),
        0.161 * (x[0] - x[2]),
    ]
    // [x[1], std::f64::consts::E * (1f64 - x[0] * x[0]) * x[1] - x[0]]
}

fn main() {
    let mut solver: solvers::RungeKuttaMethod<3, 3> =
        solvers::RungeKuttaMethod::new(f, 0f64, 800f64, 0.000005f64, &[0.5f64, 0.5f64, 0.5f64], 3, Some(10000u32));

    solver.set_butcher_table(solvers::ButcherTable::<3>::new(
        [
            [0f64, 0f64, 0f64],
            [0.5f64, 0f64, 0f64],
            [-1f64, 2f64, 0f64],
        ],
        [0f64, 0.5f64, 1f64],
        [1f64 / 6f64, 2f64 / 3f64, 1f64 / 6f64],
    ));

    let mut file = std::fs::File::create("test.csv").unwrap();

    solver.solve();

    let (t, x) = solver.get_solution();

    for i in 0..t.len() {
        file.write(format!("{}, {}, {}, {}\n", t[i], x[i][0], x[i][1], x[i][2]).as_bytes());
    }
}
