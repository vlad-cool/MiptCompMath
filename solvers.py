import numpy as np

cheats = False
verbose = False

epsilon = 1e-5

def close_enough(A, B, epsilon=epsilon):
    return np.abs(A - B) < epsilon

def solve_linear_system(A, f, epsilon=epsilon): # Метод LU разложения
    N = len(f)
    
    L = np.zeros((N, N))
    U = np.zeros((N, N))

    for i in range(N):
        L[i][i] = 1

        for j in range(i, N):
            sum = 0
            for k in range(i):
                sum += L[i][k] * U[k][j]
            U[i][j] = A[i][j] - sum

        for j in range(i, N):
            sum = 0
            for k in range(i):
                sum += L[j][k] * U[k][i]
            L[j][i] = (A[j][i] - sum) / U[i][i]

    v = np.zeros(N)

    for i in range(N):
        v[i] = f[i]
        for j in range(i):
            v[i] -= L[i][j] * v[j]

    x = np.zeros(N)

    for i in range(N - 1, -1, -1):
        x[i] = v[i] / U[i][i]
        for j in range(i + 1, N):
            x[i] -= U[i][j] * x[j] / U[i][i]

    return x

def derivative(f, x, i, j, h=5e-5):
    x_1 = np.copy(x)
    x_1[j] += h
    x_2 = np.copy(x) 
    x_2[j] -= h
    return (f(x_1).ravel()[i] - f(x_2).ravel()[i]) / (2 * h)

## Cheats enabled
if cheats:
    from scipy.optimize import fsolve

def solve_newton(f, x, epsilon=epsilon, max_iterations=100):
    x = np.copy(x)
    original_shape = x.shape
    x = x.ravel()
    
        
    def f_flat(x):
        return f(x.ravel().reshape(original_shape)).ravel()
    
    if cheats:
        return fsolve(f_flat, x).reshape(original_shape)

    N = len(f_flat(x))
        
    iterations = 0
    while not close_enough(np.sum(np.abs(f_flat(x.reshape(original_shape)))), 0):
        if iterations == max_iterations:
            raise ArithmeticError("Cannot solve system, too many iterations")
        iterations += 1
        J = np.array([[0.0] * N] * N)

        for i in range(N):
            for j in range(N):
                J[i][j] = derivative(f_flat, x, i, j)

        if verbose:
            print(J)
        
        delta_x = solve_linear_system(J, -f_flat(x), epsilon=epsilon)
        
        if np.sum(delta_x) == 0:
            raise ArithmeticError("Cannot solve system, delta x is zero")
        if np.isnan(np.sum(delta_x)):
            raise ArithmeticError("Cannot solve system, got not a number value")
        
        x = x + delta_x
        
        if verbose:
            print(f"x = {x}\ndelta_x = {delta_x}\nnorm = {np.sum(np.abs(f_flat(x.reshape(original_shape))))}")
    
    return x.reshape(original_shape)

################################################

class ButcherTable:
    def __init__(self, A):
        if len(A) == len(A[0]):
            self.a = []
            self.c = []
            self.b = []
            
            self.type = "Explicit"
            
            self.s = len(A) - 1
            
            for i in range(self.s):
                self.a.append([])
                for j in range(self.s):
                    self.a[i].append(A[i][j + 1])

                    if i > j and A[i][j + 1] != 0:
                        self.type == "Implicit"

                self.b.append(A[self.s][i + 1])
                
                self.c.append(A[i][0])
        elif len(A) == len(A[0]) + 1:
            self.a = []
            self.c = []
            self.b = []
            
            self.type = "ExplicitAdaptive"
            
            self.s = len(A) - 1
            
            for i in range(self.s):
                self.a.append([])
                for j in range(self.s):
                    self.a[i].append(A[i][j + 1])

                    if i > j and A[i][j + 1] != 0:
                        self.type == "ImplicitAdaptive"

                self.b_1.append(A[self.s + 0][i + 1])
                self.b_2.append(A[self.s + 1][i + 1])

                self.c.append(A[i][0])

class Solution:
    def __init__(self, t, x):
        self.t = t
        self.x = x

class RungeKuttaMethod:
    def __init__(self, name: str, order: int, table: ButcherTable):
        self.name = name
        self.order = order
        self.table = table
        self.solution = None
    
    def set_solution(self, t, x):
        self.solution = Solution(t, x)

def step_runge_kutta_explicit(tau, t, x, f, method, epsilon=epsilon):
    k = []
    for i in range(method.s):
        arg_1 = t + method.s * tau
        arg_2 = 0
        for j in range(i):
            arg_2 += method.a[i][j] * k[j]
        arg_2 *= tau
        arg_2 += x
        k.append(f(arg_1, arg_2))
    
    k = np.array(k)
        
    ret_val = x - x
    for i in range(method.s):
        ret_val += method.b[i] * k[i]
    ret_val *= tau
    ret_val += x
        
    return ret_val, None

def step_runge_kutta_implicit(tau, t, x, f, method, epsilon=epsilon):
    def equation(k_0):
        k = []
        for i in range(method.s):
            arg_1 = t + method.s * tau
            arg_2 = 0
            for j in range(i):
                arg_2 += method.a[i][j] * k_0[j]
            arg_2 *= tau
            arg_2 += x
            k.append(f(arg_1, arg_2))
        
        k = np.array(k)
        return k - k_0
        
    k = solve_newton(equation, np.array([f(t, x) for _ in range(method.s)]), epsilon=epsilon)
        
    ret_val = x - x
    for i in range(method.s):
        ret_val += method.b[i] * k[i]
    ret_val *= tau
    ret_val += x
    
    return ret_val, None

def step_runge_kutta_explicit_adaptive(tau, t, x, f, method, epsilon=epsilon):
    k = []
    for i in range(method.s):
        arg_1 = t + method.s * tau
        arg_2 = 0
        for j in range(i):
            arg_2 += method.a[i][j] * k[j]
        arg_2 *= tau
        arg_2 += x
        k.append(f(arg_1, arg_2))
    
    k = np.array(k)
    
    y_1 = x - x
    y_2 = x - x
    for i in range(method.s):
        y_1 += method.b_1[i] * k[i]
        y_2 += method.b_2[i] * k[i]
    y_1 *= tau
    y_2 *= tau
    y_1 += x
    y_2 += x
    
    err = np.sum(np.abs(y_2 - y_1))
    
    return y_2, err

def step_runge_kutta_implicit_adaptive(tau, t, x, f, method, epsilon=epsilon):
    def equation(k_0):
        k = []
        for i in range(method.s):
            arg_1 = t + method.s * tau
            arg_2 = 0
            for j in range(i):
                arg_2 += method.a[i][j] * k_0[j]
            arg_2 *= tau
            arg_2 += x
            k.append(f(arg_1, arg_2))
        
        k = np.array(k)
        return k - k_0
        
    k = solve_newton(equation, np.array([f(t, x) for _ in range(method.s)]), epsilon=epsilon)
        
    y_1 = x - x
    y_2 = x - x
    for i in range(method.s):
        y_1 += method.b_1[i] * k[i]
        y_2 += method.b_2[i] * k[i]
    y_1 *= tau
    y_2 *= tau
    y_1 += x
    y_2 += x
    
    err = np.sum(np.abs(y_2 - y_1))
    
    return y_2, err

def step_runge_kutta(tau, t, x, f, method, epsilon=epsilon):
    match method.type:
        case "Explicit":
            return step_runge_kutta_explicit(tau, t, x, f, method, epsilon=epsilon)
        case "Implicit":
            return step_runge_kutta_implicit(tau, t, x, f, method, epsilon=epsilon)
        case "ExplicitAdaptive":
            return step_runge_kutta_explicit_adaptive(tau, t, x, f, method, epsilon=epsilon)
        case "ImplicitAdaptive":
            return step_runge_kutta_implicit_adaptive(tau, t, x, f, method, epsilon=epsilon)

def solve_runge_kutta(f, start, stop, tau, x_0, method, epsilon=epsilon, max_err=0.1, include_every_n=1):
    t = [start]
    
    t_i = start
    x_i = np.copy(x_0)
    
    i = 0
    max_tau = tau

    x = [np.array(x_0)]
        
    while t_i <= stop:
        try:
            x_i, err = step_runge_kutta(tau, t_i, x_i, f, method.table, epsilon=epsilon)
            
            if err is not None and err != 0:
                tau *= (max_err / err) ** (1 / (method.order))
                if tau > max_tau:
                    tau = max_tau
    
            if err is None or err < max_err * 1.1:
                i += 1
                t_i += tau
                
                if i % include_every_n == 0:
                    x.append(x_i)
                    t.append(t_i)
        
        except Exception as e:
            print(e)
            print(f"Failed to solve at t = {t_i}")
            break
        
    method.set_solution(
        t,
        np.array(x),
    )
    
    return method

################################################

class AdamsMethod:
    def __init__(self, order: int, type):
        self.name = f"Adams method ({type})"
        self.order = order
        self.type = type
        self.solution = None
    
    def set_solution(self, t, x):
        self.solution = Solution(t, x)

def step_adams_explicit(tau, t, x, f, N):
    if N > len(x):
        return step_adams_explicit(tau, t, x, f, len(x))
    else:
        n = len(x) - 1
        match N:
            case 1:
                return x[n] + tau * f(t, x[n])
            case 2:
                return x[n] + tau * (3 * f(t, x[n]) - 1 * f(t, x[n - 1])) / 2
            case 3:
                return x[n] + tau * (23 * f(t, x[n]) - 16 * f(t, x[n - 1]) + 5 * f(t, x[n - 2])) / 12
            case 4:
                return x[n] + tau * (55 * f(t, x[n]) - 59 * f(t, x[n - 1]) + 37 * f(t, x[n - 2]) - 9 * f(t, x[n - 3])) / 24

def step_adams_implicit(tau, t, x, f, N, epsilon=epsilon):
    if N > len(x):
        return step_adams_implicit(tau, t, x, f, len(x))
    else:
        n = len(x) - 1
        match N:
            case 0:
                return solve_newton(lambda a: x[n] + tau * f(t, a) - a, x[n], epsilon=epsilon)
            case 1:
                return solve_newton(lambda a: x[n] + tau * (f(t, a) + f(t, x[n])) / 2 - a, x[n], epsilon=epsilon)
            case 2:
                return solve_newton(lambda a: x[n] + tau * (5 * f(t, a) + 8 * f(t, x[n]) - f(t, x[n - 1])) / 12 - a, x[n], epsilon=epsilon)
            case 3:
                return solve_newton(lambda a: x[n] + tau * (9 * f(t, a) + 19 * f(t, x[n]) - 5 * f(t, x[n - 1]) + f(t, x[n - 2])) / 24 - a, x[n], epsilon=epsilon)

def solve_adams(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1):
    t = [start]
    i = 0
        
    x = [
        np.array(x_0)
    ]
    
    t_last= [start]
    x_last = [np.array(x)]

    while t[-1] <= stop:
        if method.type.lower() == "explicit":
            x_last.append(step_adams_explicit(tau, t[-1], x, f, method.order))
        elif method.type.lower() == "implicit":
            x_last.append(step_adams_implicit(tau, t[-1], x, f, method.order))
        t_last.append(tau * i)
        
        i += 1
        
        if len(t_last) > method.order:
            t_last = t_last[len(t_last) - method.order:]
            x_last = x_last[len(t_last) - method.order:]
        
        if i % include_every_n == 0:
            x.append(x_last[-1])
            t.append(t_last[-1])
        
    method.set_solution(
        t,
        np.array(x),
    )
    
    return method

################################################

class BackwardDiffenetiationMethod:
    def __init__(self, order: int, type):
        self.name = f"Backward differentiation method method ({type})"
        self.order = order
        self.type = type
        self.solution = None
    
    def set_solution(self, t, x):
        self.solution = Solution(t, x)

def step_backward_differentiation_explicit(tau, t, x, f, N):
    if N > len(x):
        return step_backward_differentiation_explicit(tau, t, x, f, len(x))
    else:
        n = len(x) - 1
        match N:
            case 1:
                return x[n] + tau * f(t, x[n])
            case 2:
                return 2 / 3 * (tau * f(t, x[n]) - 1 / 2 * x[n - 1] + 2 * x[n])
            case 3:
                return 6 / 11 * (tau * f(t, x[n]) + 1 / 3 * x[n - 2] - 3 / 2 * x[n - 1] + 3 * x[n])
            case 4:
                return 12 / 25 * (tau * f(t, x[n]) - 1 / 4 * x[n - 3] + 4 / 3 * x[n - 2] - 3 * x[n - 1] + 4 * x[n])

def step_backward_differentiation_implicit(tau, t, x, f, N, epsilon=epsilon):
    if N > len(x):
        return step_backward_differentiation_implicit(tau, t, x, f, len(x), epsilon=epsilon)
    else:
        n = len(x) - 1
        match N:
            case 1:
                return solve_newton(lambda a: x[n] + tau * f(t, a) - a, x[n], epsilon=epsilon)
            case 2:
                return solve_newton(lambda a: 2 / 3 * (tau * f(t, a) - 1 / 2 * x[n - 1] + 2 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)
            case 3:
                return solve_newton(lambda a: 6 / 11 * (tau * f(t, a) + 1 / 3 * x[n - 2] - 3 / 2 * x[n - 1] + 3 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)
            case 4:
                return solve_newton(lambda a: 12 / 25 * (tau * f(t, a) - 1 / 4 * x[n - 3] + 4 / 3 * x[n - 2] - 3 * x[n - 1] + 4 * x[n]) - a, x[n] + (x[n] - x[n - 1]) / 2, epsilon=epsilon)

def solve_backward_differentiation_implicit(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1):
    t = [start]
    i = 0
        
    x = [
        np.array(x_0)
    ]
    
    t_last= [start]
    x_last = [np.array(x)]

    while t[-1] <= stop:
        if method.type.lower() == "explicit":
            x_last.append(step_backward_differentiation_explicit(tau, t[-1], x, f, method.order))
        elif method.type.lower() == "implicit":
            x_last.append(step_backward_differentiation_implicit(tau, t[-1], x, f, method.order))
        t_last.append(tau * i)
        
        i += 1
        
        if len(t_last) > method.order:
            t_last = t_last[len(t_last) - method.order:]
            x_last = x_last[len(t_last) - method.order:]
        
        if i % include_every_n == 0:
            x.append(x_last[-1])
            t.append(t_last[-1])
        
    method.set_solution(
        t,
        np.array(x),
    )
    
    return method
