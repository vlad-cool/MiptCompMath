// use crate::sturm_liuville_problem::*;
// use crate::cauchy_problem::*;

// pub fn solve_phase_method(problem: &mut SturmLiuvilleProblem, n: usize) -> Result<SturmLiuvilleSolution, crate::utils::SolverError> {
//     fn get_u(rho: f64, phi: f64) -> f64 {
//         rho * phi.sin()
//     }
    
//     fn get_du(rho: f64, phi: f64) -> f64 {
//         rho * phi.cos()
//     }

//     let rho_0: f64 = (problem.x_0.powi(2) + problem.dx_0.powi(2)).sqrt();
//     let phi_0: f64 = (problem.x_0 / rho_0).asin();
    
// let mut cauchy_problem: CauchyProblem<2, _> = CauchyProblem {
//                 f: &mut |_: f64, x: &[f64; 2]| [problem.f.calc(x)],
//                 start: 0.0,
//                 stop: 120.0,
//                 x_0: [0.0, 1.0],
//             };
//             let (solution, res) = cauchy_solver.solve(&mut cauchy_problem, 0.01, false, None);
//             res.expect("Faield to solve differential equation");
//             [(solution.x.last().unwrap()[0].powi(2) * 10.0
//                 + (solution.x.last().unwrap()[1] - 1.0).powi(2) * 1.0)]

//     todo!()
// }