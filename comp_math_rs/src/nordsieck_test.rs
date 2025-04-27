use std::io::Write;

use solvers::*;

mod solvers;

fn f(_t: f64, x: &[f64; 1]) -> [f64; 1] {
    [x[0].powi(1)]
}

fn write_csv<const N: usize>(
    group: String,
    path: String,
    solution: CauchySolution<N>,
    tau: f64,
    time: std::time::Duration,
) {
    let mut file = std::fs::File::create(format!("../task6_2_data/{}.csv", path))
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
    let problem: CauchyProblem<1> = CauchyProblem {
        f,
        start: 0.0,
        stop: 10.0,
        x_0: [1.0],
    };

    let print_progress = false;

    // let mut solver: solvers::AdamsMethod<1> =
    //     solvers::AdamsMethod::new(3, solvers::SolverType::Explicit);
    let mut solver: NordsieckMethod<1> =
        solvers::NordsieckMethod::new(solvers::NordsieckMethodType::ImplicitBackwardDifferentiation(5));

    // let tau: f64 = 0.01;
    let tau: f64 = 1.0;
    let save_every: u32 = 1;
    let start_time: std::time::Instant = std::time::Instant::now();
    let (solution, res) = solver.solve(&problem, tau, print_progress, Some(save_every));
    println!("{:?}", res);
    write_csv(
        "Nordsieck Method".to_string(),
        format!("nordsieck_test"),
        solution,
        tau,
        start_time.elapsed(),
    );
}
