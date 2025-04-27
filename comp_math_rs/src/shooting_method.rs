use crate::boundary_problem::BoundaryProblem;
use crate::boundary_problem::BoundarySolution;
use crate::cauchy_problem::CauchyProblem;
use crate::cauchy_problem::CauchySolver;
use crate::equation_solvers::solve_newton;

pub fn shooting_method(
    problem: &BoundaryProblem,
    cauchy_solver: &mut Box<dyn CauchySolver<2>>,
    step: f64,
) -> (BoundarySolution, Result<(), &'static str>) {
    let equation = |v: &[f64; 1]| {
        let v: f64 = v[0];
        let cauchy_problem: CauchyProblem<2> = CauchyProblem {
            f: problem.f,
            start: problem.x_0,
            stop: problem.x_n,
            x_0: [problem.y_0, v],
        };
        let (solution, res) = cauchy_solver.solve(&cauchy_problem, step, false, None);
        res.expect("Faield to solve differential equation");
        [solution.x.last().unwrap()[0] - problem.y_n]
    };

    let res: Result<[f64; 1], &str> = solve_newton(equation, &[1.0], None);
    let v: f64 = res.unwrap()[0];

    let problem: CauchyProblem<2> = CauchyProblem {
        f: problem.f,
        start: problem.x_0,
        stop: problem.x_n,
        x_0: [problem.y_0, v],
    };

    let (cauchy_solution, _res) = cauchy_solver.solve(&problem, step, false, None);

    let mut solution: BoundarySolution = BoundarySolution {
        x: std::vec::Vec::with_capacity(cauchy_solution.t.len()),
        y: std::vec::Vec::with_capacity(cauchy_solution.x.len()),
        method_name: "Shooting method".to_string(),
    };

    for i in 0..cauchy_solution.t.len() {
        solution.x.push(cauchy_solution.t[i]);
        solution.y.push(cauchy_solution.x[i][0]);
    }

    (solution, Ok(()))
}
