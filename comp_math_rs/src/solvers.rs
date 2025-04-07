use core::arch;
use std::vec;

fn close_enough(a: f64, b: f64, epsilon: f64) -> bool {
    (a - b).abs() < epsilon
}

fn close_enough_arr<const N: usize>(a: &[f64; N], b: &[f64; N], epsilon: f64) -> bool {
    let mut sum: f64 = 0f64;
    for i in 0..N {
        sum += (a[i] - b[i]).abs();
    }

    close_enough(sum, 0f64, epsilon)
}

pub fn derivative(f: fn(f64) -> f64, x: f64, h: f64) -> f64 {
    (f(x + h) - f(x - h)) / (2f64 * h)
}

pub fn partial_derivative<const N: usize>(
    f: fn(&[f64; N]) -> [f64; N],
    x: &[f64; N],
    i: usize,
    j: usize,
    h: f64,
) -> f64 {
    let mut x_1: [f64; N] = x.clone();
    let mut x_2: [f64; N] = x.clone();

    x_1[j] += h;
    x_2[j] -= h;

    (f(&x_1)[i] - f(&x_2)[i]) / (2f64 * h)
}

pub fn solve_linear_system<const N: usize>(a: &[[f64; N]; N], b: &[f64; N]) -> [f64; N] {
    let mut l: [[f64; N]; N] = [[0f64; N]; N];
    let mut u: [[f64; N]; N] = [[0f64; N]; N];

    for i in 0..N {
        l[i][i] = 1f64;
        for j in i..N {
            let mut sum: f64 = 0f64;
            for k in 0..i {
                sum += l[i][k] * u[k][j];
            }
            u[i][j] = a[i][j] - sum;
        }

        for j in i..N {
            let mut sum: f64 = 0f64;
            for k in 0..i {
                sum += l[j][k] * u[k][i];
            }
            l[j][i] = (a[j][i] - sum) / u[i][i];
        }
    }

    let mut v: [f64; N] = [0f64; N];
    for i in 0..N {
        v[i] = b[i];
        for j in 0..i {
            v[i] -= l[i][j] * v[j];
        }
    }

    let mut x: [f64; N] = [0f64; N];

    for i in (0..N).rev() {
        x[i] = v[i] / u[i][i];
        for j in i + 1..N {
            x[i] -= u[i][j] * x[j] / u[i][i];
        }
    }

    x
}

fn solve_newton<const N: usize>(
    f: fn(&[f64; N]) -> [f64; N],
    x: &[f64; N],
    max_iterations: Option<u32>,
) -> Result<[f64; N], &'static str> {
    let max_iterations: u32 = max_iterations.unwrap_or(100u32);

    let mut x: [f64; N] = x.clone();
    let mut iterations: u32 = 0u32;

    while !close_enough_arr(&f(&x), &[0f64; N], 1e-6) {
        if iterations == max_iterations {
            return Err("Cannot solve system, too many iterations");
        }

        let mut jacobian: [[f64; N]; N] = [[0f64; N]; N];

        for i in 0..N {
            for j in 0..N {
                jacobian[i][j] = -partial_derivative(f, &x, i, j, 1e-5);
            }
        }

        let delta_x = solve_linear_system(&jacobian, &f(&x));

        for i in 0..N {
            x[i] += delta_x[i];
        }

        iterations += 1;
    }

    Ok(x)
}

pub trait DifferentialEquationNumericMethod<const N: usize> {
    fn new(
        f: fn(t: f64, &[f64; N]) -> [f64; N],
        start: f64,
        stop: f64,
        tau: f64,
        x_0: &[f64; N],
        order: u32,
        save_every: Option<u32>,
    ) -> Self;
    fn solve(&mut self) -> (f64, Result<(), &'static str>);
    fn get_solution(&self) -> (std::vec::Vec<f64>, std::vec::Vec<[f64; N]>);
}

enum SolverType {
    Explicit,
    Implicit,
}

pub struct ButcherTable<const N: usize> {
    solver_type: SolverType,
    a: [[f64; N]; N],
    b: [f64; N],
    c: [f64; N],
}

impl<const N: usize> ButcherTable<N> {
    pub fn new(a: [[f64; N]; N], b: [f64; N], c: [f64; N]) -> Self {
        let mut solver_type = SolverType::Explicit;
        for i in 0..N {
            for j in i..N {
                if a[i][j] != 0f64 {
                    solver_type = SolverType::Implicit;
                }
            }
        }

        Self {
            solver_type,
            a,
            b,
            c,
        }
    }
}

pub struct RungeKuttaMethod<const N: usize, const M: usize> {
    f: fn(t: f64, &[f64; N]) -> [f64; N],
    start: f64,
    stop: f64,
    tau: f64,
    x_0: [f64; N],
    save_every: u32,
    order: u32,

    butcher_table: Option<ButcherTable<M>>,
    t: std::vec::Vec<f64>,
    solution: std::vec::Vec<[f64; N]>,
}

impl<const N: usize, const M: usize> RungeKuttaMethod<N, M> {
    pub fn set_butcher_table(&mut self, butcher_table: ButcherTable<M>) {
        self.butcher_table = Some(butcher_table);
    }

    fn step_runge_kutta_explicit(&mut self, t: f64, x: [f64; N]) -> Result<(f64, [f64; N]), &'static str> {
        let tau: f64 = self.tau;

        let butcher_table: &ButcherTable<M> = self.butcher_table.as_ref().unwrap();
        let mut k: [[f64; N]; M] = [[0f64; N]; M];
        for i in 0..M {
            let arg_1: f64 = t + tau * butcher_table.c[i];
            let mut arg_2: [f64; N] = [0f64; N];

            for a in 0..N {
                for j in 0..i {
                    arg_2[a] += butcher_table.a[i][j] * k[j][a];
                }
                arg_2[a] *= tau;
                arg_2[a] += x[a];
            }

            k[i] = (self.f)(arg_1, &arg_2);
        }

        let mut res: [f64; N] = [0f64; N];
        for j in 0..N
        {
            for i in 0..M {
                res[j] += butcher_table.b[i] * k[i][j];
            }
            res[j] *= tau;
            res[j] += x[j];
        }

        // self.solution.push(res);
        // self.t.push(t + tau);

        Ok((t + tau, res))
    }

    //     ret_val = x - x
    //     for i in range(method.s):
    //         ret_val += method.b[i] * k[i]
    //     ret_val *= tau
    //     ret_val += x

    //     return ret_val, None

    // def step_runge_kutta_implicit(tau, t, x, f, method, epsilon=epsilon):
    //     def equation(k_0):
    //         k = []
    //         for i in range(method.s):
    //             arg_1 = t + method.s * tau
    //             arg_2 = 0
    //             for j in range(i):
    //                 arg_2 += method.a[i][j] * k_0[j]
    //             arg_2 *= tau
    //             arg_2 += x
    //             k.append(f(arg_1, arg_2))

    //         k = np.array(k)
    //         return k - k_0

    //     k = solve_newton(equation, np.array([f(t, x) for _ in range(method.s)]), epsilon=epsilon)

    //     ret_val = x - x
    //     for i in range(method.s):
    //         ret_val += method.b[i] * k[i]
    //     ret_val *= tau
    //     ret_val += x

    //     return ret_val, None
}

impl<const N: usize, const M: usize> DifferentialEquationNumericMethod<N>
    for RungeKuttaMethod<N, M>
{
    fn new(
        f: fn(t: f64, &[f64; N]) -> [f64; N],
        start: f64,
        stop: f64,
        tau: f64,
        x_0: &[f64; N],
        order: u32,
        save_every: Option<u32>,
    ) -> Self {
        Self {
            f,
            start,
            stop,
            tau,
            x_0: x_0.clone(),
            order,
            save_every: save_every.unwrap_or(1),
            butcher_table: None,
            t: vec![],
            solution: vec![],
        }
    }
    fn solve(&mut self) -> (f64, Result<(), &'static str>) {
        if self.butcher_table.is_none() {
            return (self.start, Err("No butcher table is set"));
        }

        let mut t_i = self.start;
        let mut x_i = self.x_0.clone();
        let mut iterations: u32 = 0u32;
        

        while t_i < self.stop {
            match self.step_runge_kutta_explicit(t_i, x_i) {
                Ok((t, x)) => {
                    t_i = t;
                    x_i = x.clone();
                    if iterations % self.save_every == 0
                    {
                        self.t.push(t_i);
                        self.solution.push(x_i);
                        println!("t: {}, iterations: {}", t_i, iterations)
                    }
                }
                Err(a) => {
                    return (t_i, Err("Failed to solve"));
                }
            };
            iterations += 1;
        }

        (t_i, Ok(()))
    }
    fn get_solution(&self) -> (std::vec::Vec<f64>, std::vec::Vec<[f64; N]>) {
        (self.t.clone(), self.solution.clone())
    }
}

// ################################################

// class ButcherTable:
//     def __init__(self, A):
//         if len(A) == len(A[0]):
//             self.a = []
//             self.c = []
//             self.b = []

//             self.type = "Explicit"
//             self.adaptive = False

//             self.s = len(A) - 1

//             for i in range(self.s):
//                 self.a.append([])
//                 for j in range(self.s):
//                     self.a[i].append(A[i][j + 1])

//                     if i > j and A[i][j + 1] != 0:
//                         self.type == "Implicit"

//                 self.b.append(A[self.s][i + 1])

//                 self.c.append(A[i][0])
//         elif len(A) == len(A[0]) + 1:
//             self.a = []
//             self.c = []
//             self.b_1 = []
//             self.b_2 = []
//             self.adaptive = True

//             self.type = "ExplicitAdaptive"

//             self.s = len(A) - 2

//             for i in range(self.s):
//                 self.a.append([])
//                 for j in range(self.s):
//                     self.a[i].append(A[i][j + 1])

//                     if i > j and A[i][j + 1] != 0:
//                         self.type == "ImplicitAdaptive"

//                 self.b_1.append(A[self.s + 0][i + 1])
//                 self.b_2.append(A[self.s + 1][i + 1])

//                 self.c.append(A[i][0])

// class Solution:
//     def __init__(self, t, x):
//         self.t = t
//         self.x = x

// class RungeKuttaMethod:
//     def __init__(self, name: str, order: int, table: ButcherTable, step=None):
//         self.name = name
//         self.order = order
//         self.table = table
//         self.solution = None
//         self.step = step

//     def set_solution(self, t, x):
//         self.solution = Solution(t, x)

// def step_runge_kutta_explicit(tau, t, x, f, method, epsilon=epsilon):
//     k = []
//     for i in range(method.s):
//         arg_1 = t + method.s * tau
//         arg_2 = 0
//         for j in range(i):
//             arg_2 += method.a[i][j] * k[j]
//         arg_2 *= tau
//         arg_2 += x
//         k.append(f(arg_1, arg_2))

//     k = np.array(k)

//     ret_val = x - x
//     for i in range(method.s):
//         ret_val += method.b[i] * k[i]
//     ret_val *= tau
//     ret_val += x

//     return ret_val, None

// def step_runge_kutta_implicit(tau, t, x, f, method, epsilon=epsilon):
//     def equation(k_0):
//         k = []
//         for i in range(method.s):
//             arg_1 = t + method.s * tau
//             arg_2 = 0
//             for j in range(i):
//                 arg_2 += method.a[i][j] * k_0[j]
//             arg_2 *= tau
//             arg_2 += x
//             k.append(f(arg_1, arg_2))

//         k = np.array(k)
//         return k - k_0

//     k = solve_newton(equation, np.array([f(t, x) for _ in range(method.s)]), epsilon=epsilon)

//     ret_val = x - x
//     for i in range(method.s):
//         ret_val += method.b[i] * k[i]
//     ret_val *= tau
//     ret_val += x

//     return ret_val, None

// def step_runge_kutta_explicit_adaptive(tau, t, x, f, method, epsilon=epsilon):
//     k = []
//     for i in range(method.s):
//         arg_1 = t + method.s * tau
//         arg_2 = 0
//         for j in range(i):
//             arg_2 += method.a[i][j] * k[j]
//         arg_2 *= tau
//         arg_2 += x
//         k.append(f(arg_1, arg_2))

//     k = np.array(k)

//     y_1 = x - x
//     y_2 = x - x
//     for i in range(method.s):
//         y_1 += method.b_1[i] * k[i]
//         y_2 += method.b_2[i] * k[i]
//     y_1 *= tau
//     y_2 *= tau
//     y_1 += x
//     y_2 += x

//     err = np.sum(np.abs(y_2 - y_1))

//     return y_2, err

// def step_runge_kutta_implicit_adaptive(tau, t, x, f, method, epsilon=epsilon):
//     def equation(k_0):
//         k = []
//         for i in range(method.s):
//             arg_1 = t + method.s * tau
//             arg_2 = 0
//             for j in range(i):
//                 arg_2 += method.a[i][j] * k_0[j]
//             arg_2 *= tau
//             arg_2 += x
//             k.append(f(arg_1, arg_2))

//         k = np.array(k)
//         return k - k_0

//     k = solve_newton(equation, np.array([f(t, x) for _ in range(method.s)]), epsilon=epsilon)

//     y_1 = x - x
//     y_2 = x - x
//     for i in range(method.s):
//         y_1 += method.b_1[i] * k[i]
//         y_2 += method.b_2[i] * k[i]
//     y_1 *= tau
//     y_2 *= tau
//     y_1 += x
//     y_2 += x

//     err = np.sum(np.abs(y_2 - y_1))

//     return y_2, err

// def step_runge_kutta(tau, t, x, f, method, epsilon=epsilon):
//     match method.type:
//         case "Explicit":
//             return step_runge_kutta_explicit(tau, t, x, f, method, epsilon=epsilon)
//         case "Implicit":
//             return step_runge_kutta_implicit(tau, t, x, f, method, epsilon=epsilon)
//         case "ExplicitAdaptive":
//             return step_runge_kutta_explicit_adaptive(tau, t, x, f, method, epsilon=epsilon)
//         case "ImplicitAdaptive":
//             return step_runge_kutta_implicit_adaptive(tau, t, x, f, method, epsilon=epsilon)

// def solve_runge_kutta(f, start, stop, tau, x_0, method, epsilon=epsilon, max_err=0.1, include_every_n=1, print_progress=False):
//     t = [start]

//     t_i = start
//     x_i = np.copy(x_0)

//     i = 0
//     max_tau = tau

//     x = [np.array(x_0)]

//     last_output = -1
//     start_time = time.time()

//     while t_i <= stop:

//         try:
//             x_i, err = step_runge_kutta(tau, t_i, x_i, f, method.table, epsilon=epsilon)

//             if err is not None and err != 0:
//                 tau *= (max_err / err) ** (1 / (method.order))
//                 if tau > max_tau:
//                     tau = max_tau

//             if err is None or err < max_err * 1.1:
//                 if print_progress and time.time() - last_output > 1:
//                     if method.table.adaptive:
//                         print(f"\rtau: {tau:.8f}, t: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
//                     else:
//                         print(f"\rt: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
//                     last_output = time.time()

//                 i += 1
//                 t_i += tau

//                 if i % include_every_n == 0:
//                     x.append(x_i)
//                     t.append(t_i)

//         except Exception as e:
//             print(e)
//             print(f"Failed to solve at t = {t_i}")
//             break

//     method.set_solution(
//         t,
//         np.array(x),
//     )

//     return method

// ################################################

// class AdamsMethod:
//     def __init__(self, order: int, type, step=None):
//         self.name = f"Adams method ({type})"
//         self.order = order
//         self.type = type
//         self.solution = None
//         self.step = step

//     def set_solution(self, t, x):
//         self.solution = Solution(t, x)

// def step_adams_explicit(tau, t, x, f, N):
//     if N > len(x):
//         return step_adams_explicit(tau, t, x, f, len(x))
//     else:
//         n = len(x) - 1
//         match N:
//             case 1:
//                 return x[n] + tau * f(t, x[n])
//             case 2:
//                 return x[n] + tau * (3 * f(t, x[n]) - 1 * f(t, x[n - 1])) / 2
//             case 3:
//                 return x[n] + tau * (23 * f(t, x[n]) - 16 * f(t, x[n - 1]) + 5 * f(t, x[n - 2])) / 12
//             case 4:
//                 return x[n] + tau * (55 * f(t, x[n]) - 59 * f(t, x[n - 1]) + 37 * f(t, x[n - 2]) - 9 * f(t, x[n - 3])) / 24

// def step_adams_implicit(tau, t, x, f, N, epsilon=epsilon):
//     if N > len(x):
//         return step_adams_implicit(tau, t, x, f, len(x))
//     else:
//         n = len(x) - 1
//         match N:
//             case 0:
//                 return solve_newton(lambda a: x[n] + tau * f(t, a) - a, x[n], epsilon=epsilon)
//             case 1:
//                 return solve_newton(lambda a: x[n] + tau * (f(t, a) + f(t, x[n])) / 2 - a, x[n], epsilon=epsilon)
//             case 2:
//                 return solve_newton(lambda a: x[n] + tau * (5 * f(t, a) + 8 * f(t, x[n]) - f(t, x[n - 1])) / 12 - a, x[n], epsilon=epsilon)
//             case 3:
//                 return solve_newton(lambda a: x[n] + tau * (9 * f(t, a) + 19 * f(t, x[n]) - 5 * f(t, x[n - 1]) + f(t, x[n - 2])) / 24 - a, x[n], epsilon=epsilon)

// def solve_adams(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
//     t = [start]
//     i = 1

//     x = [
//         np.copy(x_0)
//     ]

//     print(method.name, method.order)

//     while tau * i <= stop:
//         try:
//             x.append(step_adams_implicit(tau, t[i - 1], x, f, method.order))
//         except Exception as e:
//             print(e)
//             print(f"Failed to solve at t = {tau * i}")
//             break

//         t.append(tau * i)

//         i += 1

//         if verbose:
//             print(tau * i)

//     method.set_solution(
//         t,
//         np.array(x),
//     )

//     return method

// ################################################

// class BackwardDiffenetiationMethod:
//     def __init__(self, order: int, type, step=None):
//         self.name = f"Backward differentiation method method ({type})"
//         self.order = order
//         self.type = type
//         self.solution = None
//         self.step = step

//     def set_solution(self, t, x):
//         self.solution = Solution(t, x)

// def step_backward_differentiation_explicit(tau, t, x, f, N):
//     if N > len(x):
//         return step_backward_differentiation_explicit(tau, t, x, f, len(x))
//     else:
//         n = len(x) - 1
//         match N:
//             case 1:
//                 return x[n] + tau * f(t, x[n])
//             case 2:
//                 return 2 / 3 * (tau * f(t, x[n]) - 1 / 2 * x[n - 1] + 2 * x[n])
//             case 3:
//                 return 6 / 11 * (tau * f(t, x[n]) + 1 / 3 * x[n - 2] - 3 / 2 * x[n - 1] + 3 * x[n])
//             case 4:
//                 return 12 / 25 * (tau * f(t, x[n]) - 1 / 4 * x[n - 3] + 4 / 3 * x[n - 2] - 3 * x[n - 1] + 4 * x[n])

// def step_backward_differentiation_implicit(tau, t, x, f, N, epsilon=epsilon):
//     if N > len(x):
//         return step_backward_differentiation_implicit(tau, t, x, f, len(x), epsilon=epsilon)
//     else:
//         n = len(x) - 1
//         match N:
//             case 1:
//                 return solve_newton(lambda a: x[n] + tau * f(t, a) - a, x[n], epsilon=epsilon)
//             case 2:
//                 return solve_newton(lambda a: 2 / 3 * (tau * f(t, a) - 1 / 2 * x[n - 1] + 2 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)
//             case 3:
//                 return solve_newton(lambda a: 6 / 11 * (tau * f(t, a) + 1 / 3 * x[n - 2] - 3 / 2 * x[n - 1] + 3 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)
//             case 4:
//                 return solve_newton(lambda a: 12 / 25 * (tau * f(t, a) - 1 / 4 * x[n - 3] + 4 / 3 * x[n - 2] - 3 * x[n - 1] + 4 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)

// def solve_backward_differentiation(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
//     t = [start]
//     i = 0

//     x = [
//         np.array(x_0)
//     ]

//     t_last= [start]
//     x_last = [np.array(x)]

//     last_output = -1
//     start_time = time.time()

//     while t[-1] <= stop:
//         if print_progress and time.time() - last_output > 1:
//             print(f"\rt: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
//             last_output = time.time()

//         if method.type.lower() == "explicit":
//             x_last.append(step_backward_differentiation_explicit(tau, t[-1], x, f, method.order))
//         elif method.type.lower() == "implicit":
//             x_last.append(step_backward_differentiation_implicit(tau, t[-1], x, f, method.order))
//         t_last.append(tau * i)

//         i += 1

//         if len(t_last) > method.order:
//             t_last = t_last[len(t_last) - method.order:]
//             x_last = x_last[len(t_last) - method.order:]

//         if i % include_every_n == 0:
//             x.append(x_last[-1])
//             t.append(t_last[-1])

//     method.set_solution(
//         t,
//         np.array(x),
//     )

//     return method

// ################################################

// class NordsieckMethod:
//     def __init__(self, type, order: int, step=None):
//         self.name = f"Nordesieck representation of {type} method of order {order}"
//         self.order = order
//         self.type = type
//         self.solution = None
//         self.step = step
//         self.l = None

//         match type.lower():
//             case "implicit_adams":
//                 match order:
//                     case 1:
//                         self.l = [1/2, 1]
//                     case 2:
//                         self.l = [5/12, 1, 1/2]
//                     case 3:
//                         self.l = [3/8, 1, 3/4, 1/6]
//                     case 4:
//                         self.l = [251/720, 1, 11/12, 1/3, 1/24]
//                     case 5:
//                         self.l = [95/288, 1, 25/24, 35/72, 5/48, 1/120]
//                     case 6:
//                         self.l = [19087/60480, 1, 137/120, 5/8, 17/96, 1/40, 1/720]
//             case "implicit_backward_differentiation":
//                 match order:
//                     case 1:
//                         self.l = [1, 1]
//                     case 2:
//                         self.l = [2/3, 1, 1/3]
//                     case 3:
//                         self.l = [6/11, 1, 6/11, 1/11]
//                     case 4:
//                         self.l = [12/15, 1, 7/10, 1/5, 1/50]
//                     case 5:
//                         self.l = [60/137, 1, 225/274, 85/274, 15/274, 1/274]
//                     case 6:
//                         self.l = [20/49, 1, 58/63, 5/12, 25/252, 1/84, 1/1764]

//         if self.l is None:
//             raise ValueError("Unknown type or order")

//         self.l = np.array(self.l)

//     def set_solution(self, t, x):
//         self.solution = Solution(t, x)

// def step_nordsieck(tau, t, f, z, l):
//     # print(f"z: {z}", flush=True)
//     k = len(l)
//     m = len(z[0])
//     # E = np.eye(m)

//     P = np.zeros((k, k))
//     for i in range(k):
//         for j in range(k):
//             if i <= j:
//                 P[i][j] = math.factorial(j) // (math.factorial(i) * math.factorial(j - i))

//     e_1 = np.zeros(k)
//     e_1[1] = 1
//     # e_1[0] = 1

//     def equation(a):
//         res = []
//         for i in range(len(a[0])):
//             z_i = np.zeros(k)
//             a_i = np.zeros(k)
//             # print(z[i], flush=True)
//             for j in range(k):
//                 z_i[j] = z[j][i]
//                 a_i[j] = a[j][i]
//             # print(P)
//             # print(z_i)
//             # print("", flush=True)
//             # print(P @ z_i)
//             # print(f(t + tau, a[0])[i])
//             # print(e_1 @ P @ z_i)
//             # print("", flush=True)
//             res.append(P @ z_i + l * (tau * f(t + tau, a[0])[i] - e_1 @ P @ z_i) - a_i)

//         # print("res", np.array(res), flush=True)
//         return np.array(res)

//     # print("AAA", equation(z), flush=True)

//     # z_flat = z.flatten()
//     return solve_newton(equation, z)
//     # return solve_newton(lambda a: np.kron(P, E) * z_flat + np.kron(l, E) * (tau * f(t + tau, a[0]) - (e_1 @ np.kron(P, E)) @ z_flat) - a, z) # .reshape(z.shape)

// def solve_nordsieck(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
//     # t = [start]
//     # i = 0

//     # x = [
//     #     np.array(x_0)
//     # ]

//     # t_last= [start]
//     # x_last = [np.array(x)]

//     z_0 = np.zeros((len(method.l), len(x_0)))
//     z_0[0] = x_0
//     z_0[1] = f(start, x_0)

//     z = [z_0]

//     # z_0 =

//     # last_output = -1
//     # start_time = time.time()

//     last_output = -1
//     start_time = time.time()

//     t = [start]
//     # i = 1

//     print(method.name, method.order)
//     # def step_nordsieck(tau, t, f, z, l):
//     while t[-1] <= stop:
//         try:
//         # if True:
//             z.append(step_nordsieck(tau, t[-1], f, z[-1], method.l))
//         except Exception as e:
//             print(e)
//             print(f"Failed to solve at t = {t[-1] + tau}")
//             break

//         t.append(t[-1] + tau)

//     if print_progress and time.time() - last_output > 1:
//         if method.table.adaptive:
//             print(f"\rtau: {tau:.8f}, t: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
//         else:
//             print(f"\rt: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
//         last_output = time.time()

//     method.set_solution(
//         t,
//         np.array([i[0] for i in z]),
//     )

//     return method

#[cfg(test)]
mod tests {
    use crate::solvers::*;

    #[test]
    fn test_solve_linear_system() {
        const N: usize = 3;
        let a = [
            [-4f64, 9f64, -4f64],
            [-5f64, -5f64, 6f64],
            [2f64, 5f64, -8f64],
        ];

        let b = [-64f64, 104.6f64, -85.2f64];

        let x = crate::solvers::solve_linear_system(&a, &b);
        for i in 0..N {
            let mut sum: f64 = 0f64;
            for j in 0..N {
                sum += a[i][j] * x[j];
            }
            assert!(close_enough(sum, b[i], 1e-6));
        }
    }

    #[test]
    fn test_derivative() {
        fn f(x: f64) -> f64 {
            x.sin()
        }

        assert!(close_enough(derivative(f, 0f64, 1e-6), 1f64, 1e-6));
        assert!(close_enough(
            derivative(f, std::f64::consts::FRAC_PI_2, 1e-6),
            0f64,
            1e-6
        ));
    }

    #[test]
    fn test_partial_derivative() {
        fn f(x: &[f64; 2]) -> [f64; 2] {
            [x[0] * x[0] + x[1] * x[1], 2f64 * x[0] * x[1]]
        }

        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 0, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 0, 1, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 1, 0, 1e-6),
            0f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[0f64, 0f64], 1, 1, 1e-6),
            0f64,
            1e-6
        ));

        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 0, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 0, 1, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 1, 0, 1e-6),
            2f64,
            1e-6
        ));
        assert!(close_enough(
            partial_derivative(f, &[1f64, 1f64], 1, 1, 1e-6),
            2f64,
            1e-6
        ));
    }

    #[test]
    fn test_solve_newton() {
        fn f(x: &[f64; 2]) -> [f64; 2] {
            [
                (x[0] + 1f64).sin() - x[1] - 1.2f64,
                2f64 * x[0] + x[1].cos() - 2f64,
            ]
        }

        let solution = solve_newton(f, &[0f64, 0f64], None).unwrap();

        assert!(close_enough_arr(&f(&solution), &[0f64, 0f64], 1e-6));
    }
}
