import numpy as np
import time
import math

verbose = False
cheats = False
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
            self.adaptive = False
            
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
            self.b_1 = []
            self.b_2 = []
            self.adaptive = True
            
            self.type = "ExplicitAdaptive"
            
            self.s = len(A) - 2
            
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
    def __init__(self, name: str, order: int, table: ButcherTable, step=None):
        self.name = name
        self.order = order
        self.table = table
        self.solution = None
        self.step = step
    
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

def solve_runge_kutta(f, start, stop, tau, x_0, method, epsilon=epsilon, max_err=0.1, include_every_n=1, print_progress=False):
    t = [start]
    
    t_i = start
    x_i = np.copy(x_0)
    
    i = 0
    max_tau = tau

    x = [np.array(x_0)]
        
    
    last_output = -1
    start_time = time.time()

    while t_i <= stop:
        
        try:
            x_i, err = step_runge_kutta(tau, t_i, x_i, f, method.table, epsilon=epsilon)
            
            if err is not None and err != 0:
                tau *= (max_err / err) ** (1 / (method.order))
                if tau > max_tau:
                    tau = max_tau
    
            if err is None or err < max_err * 1.1:
                if print_progress and time.time() - last_output > 1:
                    if method.table.adaptive:
                        print(f"\rtau: {tau:.8f}, t: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
                    else:
                        print(f"\rt: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
                    last_output = time.time()
                
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
    def __init__(self, order: int, type, step=None):
        self.name = f"Adams method ({type})"
        self.order = order
        self.type = type
        self.solution = None
        self.step = step
    
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


def solve_adams(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
    t = [start]
    i = 1

    x = [
        np.copy(x_0)
    ]

    print(method.name, method.order)
    
    while tau * i <= stop:
        try:
            x.append(step_adams_implicit(tau, t[i - 1], x, f, method.order))
        except Exception as e:
            print(e)
            print(f"Failed to solve at t = {tau * i}")
            break
        
        t.append(tau * i)
        
        i += 1
        
        if verbose:
            print(tau * i)
    
    method.set_solution(
        t, 
        np.array(x),
    )
    
    return method

################################################

class BackwardDiffenetiationMethod:
    def __init__(self, order: int, type, step=None):
        self.name = f"Backward differentiation method method ({type})"
        self.order = order
        self.type = type
        self.solution = None
        self.step = step
    
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

def solve_backward_differentiation(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
    t = [start]
    i = 0
        
    x = [
        np.array(x_0)
    ]
    
    t_last= [start]
    x_last = [np.array(x)]
    
    last_output = -1
    start_time = time.time()

    while t[-1] <= stop:
        if print_progress and time.time() - last_output > 1:
            print(f"\rt: {t[-1]:.8f}, len: {len(t):010}, {time.time() - start_time:010.8f}", end="")
            last_output = time.time()
        
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

################################################

class NordsieckMethod:
    def __init__(self, type, order: int, step=None):
        self.name = f"Nordesieck representation of {type} method of order {order}"
        self.order = order
        self.type = type
        self.solution = None
        self.step = step
        self.l = None
        
        match type.lower():
            case "implicit_adams":
                match order:
                    case 1:
                        self.l = [1/2, 1]
                    case 2:
                        self.l = [5/12, 1, 1/2]
                    case 3:
                        self.l = [3/8, 1, 3/4, 1/6]
                    case 4:
                        self.l = [251/720, 1, 11/12, 1/3, 1/24]
                    case 5:
                        self.l = [95/288, 1, 25/24, 35/72, 5/48, 1/120]
                    case 6:
                        self.l = [19087/60480, 1, 137/120, 5/8, 17/96, 1/40, 1/720]
            case "implicit_backward_differentiation":
                match order:
                    case 1:
                        self.l = [1, 1]
                    case 2:
                        self.l = [2/3, 1, 1/3]
                    case 3:
                        self.l = [6/11, 1, 6/11, 1/11]
                    case 4:
                        self.l = [12/15, 1, 7/10, 1/5, 1/50]
                    case 5:
                        self.l = [60/137, 1, 225/274, 85/274, 15/274, 1/274]
                    case 6:
                        self.l = [20/49, 1, 58/63, 5/12, 25/252, 1/84, 1/1764]
        
        if self.l is None:
            raise ValueError("Unknown type or order")

        self.l = np.array(self.l)
    
    def set_solution(self, t, x):
        self.solution = Solution(t, x)

def step_nordsieck(tau, t, f, z, l):
    # print(f"z: {z}", flush=True)
    k = len(l)
    m = len(z[0])
    # E = np.eye(m)
    
    P = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            if i <= j:
                P[i][j] = math.factorial(j) // (math.factorial(i) * math.factorial(j - i))
   
    e_1 = np.zeros(k)
    e_1[1] = 1
    # e_1[0] = 1
    
    def equation(a):
        res = []
        for i in range(len(a[0])):
            z_i = np.zeros(k)
            a_i = np.zeros(k)
            # print(z[i], flush=True)
            for j in range(k):
                z_i[j] = z[j][i]
                a_i[j] = a[j][i]
            # print(P)
            # print(z_i)
            # print("", flush=True)
            # print(P @ z_i)
            # print(f(t + tau, a[0])[i])
            # print(e_1 @ P @ z_i)
            # print("", flush=True)
            res.append(P @ z_i + l * (tau * f(t + tau, a[0])[i] - e_1 @ P @ z_i) - a_i)
        
        # print("res", np.array(res), flush=True)
        return np.array(res)
    
    # print("AAA", equation(z), flush=True)
    
    # z_flat = z.flatten()
    return solve_newton(equation, z)
    # return solve_newton(lambda a: np.kron(P, E) * z_flat + np.kron(l, E) * (tau * f(t + tau, a[0]) - (e_1 @ np.kron(P, E)) @ z_flat) - a, z) # .reshape(z.shape)
    
def solve_nordsieck(f, start, stop, tau, x_0, method, epsilon=epsilon, include_every_n=1, print_progress=False):
    # t = [start]
    # i = 0
        
    # x = [
    #     np.array(x_0)
    # ]
    
    # t_last= [start]
    # x_last = [np.array(x)]
    
    z_0 = np.zeros((len(method.l), len(x_0)))
    z_0[0] = x_0
    z_0[1] = f(start, x_0)
    
    z = [z_0]
    
    # z_0 = 
    
    # last_output = -1
    # start_time = time.time()

    t = [start]
    # i = 1

    print(method.name, method.order)
    # def step_nordsieck(tau, t, f, z, l):
    while t[-1] <= stop:
        try:
        # if True:
            z.append(step_nordsieck(tau, t[-1], f, z[-1], method.l))
        except Exception as e:
            print(e)
            print(f"Failed to solve at t = {t[-1] + tau}")
            break
        
        t.append(t[-1] + tau)
        
    method.set_solution(
        t, 
        np.array([i[0] for i in z]),
    )
    
    return method