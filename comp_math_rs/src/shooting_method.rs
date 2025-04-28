use crate::algebraic_equation_solvers::solve_newton;
use crate::boundary_problem::BoundarySolution;
use crate::boundary_problem::NonlinearBoundaryProblem;
use crate::cauchy_problem::CauchyProblem;
use crate::cauchy_problem::CauchySolver;

pub fn shooting_method<F>(
    problem: &mut NonlinearBoundaryProblem<F>,
    cauchy_solver: &mut Box<dyn CauchySolver<2, F>>,
    step: f64,
) -> (BoundarySolution, Result<(), &'static str>)
where
    F: FnMut(f64, &[f64; 2]) -> [f64; 2],
{
    let equation =
        |v: &[f64; 1]| {
            let v: f64 = v[0];
            let (y, dy) = if problem.a_1 == 0.0 {
                (problem.u_1 / problem.b_1, v)
            } else if problem.b_1 == 0.0 {
                (v, problem.u_1 / problem.a_1)
            } else {
                ((problem.u_1 - problem.a_1 * v) / problem.b_1, v)
            };
            let mut cauchy_problem: CauchyProblem<2, F> = CauchyProblem {
                f: &mut problem.f,
                start: problem.x_0,
                stop: problem.x_n,
                x_0: [y, dy],
            };
            let (solution, res) = cauchy_solver.solve(&mut cauchy_problem, step, false, None);
            res.expect("Faield to solve differential equation");
            [solution.x.last().unwrap()[1] * problem.a_2
                + solution.x.last().unwrap()[0] * problem.b_2
                - problem.u_2]
        };

    let res: Result<[f64; 1], &str> = solve_newton(equation, &[0.0], None);
    let v: f64 = res.unwrap()[0];

    let (y, dy) = if problem.a_1 == 0.0 {
        (problem.u_1 / problem.b_1, v)
    } else if problem.b_1 == 0.0 {
        (v, problem.u_1 / problem.a_1)
    } else {
        ((problem.u_1 - problem.a_1 * v) / problem.b_1, v)
    };
    let mut problem: CauchyProblem<2, F> = CauchyProblem {
        f: &mut problem.f,
        start: problem.x_0,
        stop: problem.x_n,
        x_0: [y, dy],
    };

    let (cauchy_solution, _res) = cauchy_solver.solve(&mut problem, step, false, None);

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
