use crate::cauchy_problem::*;

struct ShootingSLMethod {
    cauchy_solver: Box<dyn CauchySolver>,
}

impl Sh
