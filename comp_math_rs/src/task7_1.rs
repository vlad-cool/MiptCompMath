use std::io::Write;

use solvers::*;

mod solvers;

fn write_csv<const N: usize>(
    path: String,
    solution: CauchySolution<N>,
    tau: f64,
    time: std::time::Duration,
) {
    let mut file = std::fs::File::create(format!("../task7_1_data/{}.csv", path))
        .expect("Failed to open file");
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
    let tau = 0.001;
    let equation = |v: &[f64; 1]| {
        let v: f64 = v[0];
        let problem: CauchyProblem<2> = CauchyProblem {
            f: |_: f64, x: &[f64; 2]| [-x[0].powi(2) / (2.0 - x[1]), x[0]],
            start: 0.0,
            stop: 1.0,
            x_0: [1.9, v],
        };
        let mut solver: solvers::BackwardDifferentiationMethod<2> =
            solvers::BackwardDifferentiationMethod::new(1, solvers::SolverType::Explicit);
        let (solution, res) = solver.solve(&problem, tau, false, None);
        res.expect("Faield to solve differential equation");
        println!("{}", solution.x.last().unwrap()[0]);
        [solution.x.last().unwrap()[0]]
    };

    let start_time: std::time::Instant = std::time::Instant::now();
    let res: Result<[f64; 1], &str> = solve_newton(equation, &[1.0], None);
    let duration: std::time::Duration = start_time.elapsed();
    let v: f64 = res.unwrap()[0];
    println!("v: {}", v);

    let problem: CauchyProblem<2> = CauchyProblem {
        f: |_: f64, x: &[f64; 2]| [-x[0].powi(2) / (2.0 - x[1]), x[0]],
        start: 0.0,
        stop: 1.0,
        x_0: [1.9, v],
    };
    let mut solver: solvers::BackwardDifferentiationMethod<2> =
        solvers::BackwardDifferentiationMethod::new(1, solvers::SolverType::Explicit);
    let (solution, res) = solver.solve(&problem, 0.001, false, None);
    res.expect("Faield to solve differential equation");

    write_csv(format!("shooting_method"), solution, tau, duration);
}
