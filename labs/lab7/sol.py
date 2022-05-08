import numpy as np
import matplotlib.pyplot as plt
import math as mt
from scipy.optimize import curve_fit
import random

# Просто линейная функция
def Linear(x, a, b):
        return a * x + b

# Вычисляет максимум отклонения сеточной функции x_arr -> y_arr от функции func на сетке x_arr
def CalculateError(x_arr, y_arr, func):
    err = []
    for i in range(0, len(x_arr)):
        err.append(np.abs(y_arr[i] - func(x_arr[i])))
    return np.max(err)


# Метод прогонки для упакованной трёхдиагональной матрицы вида
def ThomasCompact(a : np.array, f : np.array):
    n = f.shape[0]
    if a.shape[0] != n or a.shape[1] != 3:
        raise RuntimeError('Matrix and/or vector dimensions do not correspond')
    u = np.zeros(n)
    v = np.zeros(n)
    u[0] = - a[0][2] / a[0][1]
    v[0] = f[0] / a[0][1]
    for i in range(1, n - 1):
        denum = a[i][1] + a[i][0]*u[i-1]
        u[i] = -a[i][2]/denum
        v[i] = (f[i] - a[i][0]*v[i-1])/denum
    x = np.zeros(n)
    x[n-1] = (f[n-1] - a[n-1][0] * v[n-2]) / (a[n-1][1] + a[n-1][0] * u[n-2])
    for i in range(n-2, -1, -1):
        x[i] = u[i]*x[i+1] + v[i]
    return x

# Метод касательнызх
def Tangent(f, x0, epsilon):
    x = np.copy(x0)
    xp = x0 - 2*epsilon
    n = 0 # Счётчик итераций
    f_x = f(x)
    while(abs(f_x) > epsilon):
        df = (f_x - f(xp))/(x - xp)
        xn = x - f_x/df
        xp = np.copy(x)
        x = np.copy(xn)
        n += 1
        f_x = f(x)
    return x, n


# Явная схема Эйлера
def EulerExplicit(x_prev : np.array, t_prev, step, f):
    return x_prev + step * f(t_prev, x_prev)

# Р-К 4 порядка
def RungeKutta(x_prev, t_prev, step, f):
    k1 = f(t_prev, x_prev)
    k2 = f(t_prev + (step/2), x_prev + (k1 * (step/2)))
    k3 = f(t_prev + (step/2), x_prev + (k2 * (step/2)))
    k4 = f(t_prev + step, x_prev + (k3 * step))
    return x_prev + (step/6) * (k1 + 2 * k2 + 2 * k3 + k4)


def CauchyProblem(diff_scheme, f, x0 : np.array, t_0, t_max, tau):
    order = x0.shape[0] # Порядок системы
    n = mt.floor((t_max - t_0) / tau)
    last = t_max - t_0 - n * tau
    x = np.zeros([n + 2, order])
    t = np.full(n + 2, t_0)
    x[0] = np.copy(x0)
    for i in range(0, n + 1):
        t[i] += i * tau
    t[n + 1] = t_max
    for i in range(1, n + 1):
        x[i] = diff_scheme(x[i - 1], t[i - 1], tau, f)
    x[n + 1] = diff_scheme(x[n], t_max - last, last, f)
    return t, x

def CauchyTerminalPoint(diff_scheme, f, x0 : np.array, t_0, t_max, tau):
    order = x0.shape[0] # Порядок системы
    n = mt.floor((t_max - t_0) / tau)
    last = t_max - t_0 - n * tau
    x = np.zeros(order)
    x_prev = np.copy(x0)
    for i in range(1, n + 1):
        t = t_0 + tau * i
        x = diff_scheme(x_prev, t - tau, tau, f)
        x_prev = np.copy(x)
    x = diff_scheme(x_prev, t_max - last, last, f)
    return x

def SolveLinear(
    p, q, f, # Уравнение
    a, b, y_a, y_b, # Краевые условия
    int_step, # Шаг интегрирования решения
    method = EulerExplicit, # Разностная схема для интегрирования
   ):
    def f_vec(t, x : np.array):
        return np.array([x[1], -p(t) * x[1] - q(t) * x[0] + f(t)])
    def f_homo(t, x : np.array):
        return np.array([x[1], -p(t) * x[1] - q(t) * x[0]])
    
    x, no_homo_y = CauchyProblem(method, f_vec, np.array([y_a, 0.0]), a, b, int_step)
    x, homo_y_b = CauchyProblem(method, f_homo, np.array([0.0, random.uniform(-5, 5)]), a, b, int_step)
    while (abs(homo_y_b[-1, 0]) < int_step):
        x, homo_y_b = CauchyTerminalPoint(method, f_homo, np.array([0.0, random.uniform(-5, 5)]), a, b, int_step)
    c = (y_b - no_homo_y[-1, 0]) / homo_y_b[-1, 0]
    y = homo_y_b[:, 0] * c + no_homo_y[:, 0]
    return x, y

def Shooting(p, q, f, a, b, y_a, y_b, int_step, method = EulerExplicit, epsilon = 0.01):
    def f_vec(t, x : np.array):
        return np.array([x[1], -p(t) * x[1] - q(t) * x[0] + f(t)])
    
    def y_p(parameter):
        return CauchyTerminalPoint(method, f_vec, np.array([y_a, parameter]), a, b, epsilon)[0] - y_b
    actual_p = Tangent(y_p, 10.0, epsilon)[0]
    t_sol, x_sol = CauchyProblem(method, f_vec, np.array([y_a, actual_p]), a, b, int_step)
    return t_sol, x_sol[:,0]

def FinDiff(
    p, q, f, # Уравнение
    a, b, y_a, y_b, # Краевые условия
    int_step # Шаг интегрирования решения
    ):
    n = mt.ceil((b - a) / int_step)
    h = (b - a) / n #new_step
    x = np.zeros(n + 1)
    for i in range(0, n):
        x[i] = a + i * h
    x[n] = b
    h_sq = 1 / (h * h)
    net_matrix = np.zeros([n - 1, 3]) # [0], ..., [n - 2]
    net_column = np.zeros(n - 1)
    net_matrix[0] = np.array([
        0,
        q(x[1]) - 2 * h_sq,
        h_sq + p(x[1]) / (2 * h)
    ])
    net_column[0] = f(x[1]) - y_a * (h_sq - p(x[1]) / (2*h))

    net_matrix[n - 2] = np.array([
        h_sq - p(x[n - 1]) / (2 * h),
        q(x[n - 1]) - 2 * h_sq,
        0
    ])
    net_column[n - 2] = f(x[n - 1]) - y_b * (h_sq + p(x[n - 1]) / (2*h))
    
    for i in range(1, n - 2):
        net_matrix[i] = np.array([
            h_sq - p(x[i + 1]) / (2*h),
            q(x[i + 1]) - 2*h_sq,
            h_sq + p(x[i + 1]) / (2*h)
        ])
        net_column[i] = f(x[i + 1])
    y = np.array([y_a])
    y = np.append(y, ThomasCompact(net_matrix, net_column))
    y = np.append(y, y_b)
    return x, y

def FinDiffMixedCond(
    p, q, f, # Уравнение
    a, b, dy_a, y_b, # Краевые условия
    int_step # Шаг интегрирования решения
    ):
    n = mt.ceil((b - a) / int_step)
    h = (b - a) / n #new_step
    x = np.zeros(n + 1)
    for i in range(0, n):
        x[i] = a + i * h
    x[n] = b
    h_sq = 1 / (h * h)
    net_matrix = np.zeros([n, 3])
    net_column = np.zeros(n)
    net_matrix[0] = np.array([0, q(x[0]) - 2*h_sq, 2*h_sq])
    net_column[0] = f(x[0]) - p(x[0]) * dy_a + 2 * dy_a / h
    net_matrix[n - 1] = np.array([
        h_sq - p(x[n - 1]) / (2*h),
        q(x[n - 1]) - 2 * h_sq,
        0
    ])
    net_column[n - 1] = f(x[n - 1]) - y_b * (h_sq + p(x[n - 1]) / (2*h))
    for i in range(1, n - 1):
        net_matrix[i] = np.array([
            h_sq - p(x[i]) / (2*h),
            q(x[i]) - 2*h_sq,
            h_sq + p(x[i]) / (2*h)
        ])
        net_column[i] = f(x[i])
    y = ThomasCompact(net_matrix, net_column)
    y = np.append(y, y_b)
    return x, y

def paint_orig(x, y, name):
    y_orig = x * np.log(x) + x ** 3 / 4 + x
    plt.scatter(x, y, color = 'r', label = name + ' solution', linestyle = ':')
    plt.plot(x, y_orig, color = 'black', label = "original_solution", linestyle = ':')
    plt.legend()
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.grid()
    plt.show()

def PlotErrors(
    p, q, f,
    a, b, y_a, y_b, dy_a,
    solution
    ):
    err_lin = []
    err_lin_2_ord = []
    err_lin_rk = []
    err_sh = []
    err_sh_2ord = []
    err_fd = []
    err_fd_mx = []
    log_h = []
    for h in np.linspace(0.1, 0.001, 4):
        print('step = ', h)
        x1, y1 = Shooting(p, q, f, a, b, y_a, y_b, h)
        x2, y2 = FinDiff(p, q, f, a, b, y_a, y_b, h)
        x3, y3 = FinDiffMixedCond(p, q, f, a, b, dy_a, y_b, h)
        x7, y7 = SolveLinear(p, q, f, a, b, y_a, y_b, h, RungeKutta)
        err_lin_rk.append(np.log10(CalculateError(x7, y7, solution)))
        err_sh.append(np.log10(CalculateError(x1, y1, solution)))
        err_fd.append(np.log10(CalculateError(x2, y2, solution)))
        err_fd_mx.append(np.log10(CalculateError(x3, y3, solution)))
        log_h.append(np.log10(h))

    popt, pcov = curve_fit(Linear, log_h, err_lin_rk, p0 = [0, -1])
    lin_rk_ord = popt[0]
    popt, pcov = curve_fit(Linear, log_h, err_sh, p0 = [0, -1])
    sh_ord = popt[0]
    popt, pcov = curve_fit(Linear, log_h, err_fd, p0 = [0, -2])
    fd_ord = popt[0]
    popt, pcov = curve_fit(Linear, log_h, err_fd_mx, p0 = [0, -2])
    fdmx_ord = popt[0]

    paint_orig(x1, y1, "Meтод стрельбы")
    paint_orig(x2, y2, "Метод сеток")
    paint_orig(x3, y3, "Метод сеток с фиктивной точкой")
    paint_orig(x7, y7, "Линейность + схема Рунге-Кутта")
    plt.grid(True)
    plt.ylabel('lg(r)')
    plt.xlabel('lg(h)')
    fmt = 'Использование линейности + схема Рунге-Кутты, k = {order:.3f}'
    plt.plot(log_h, err_lin_rk, label = fmt.format(order = lin_rk_ord))
    fmt = 'Метод стрельбы + схема Эйлера 1 порядка, k = {order:.3f}'
    plt.plot(log_h, err_sh, '-.', label = fmt.format(order = 2 * sh_ord))
    fmt = 'Метод сеток, k = {order:.3f}'
    plt.plot(log_h, err_fd, '--', label = fmt.format(order = fd_ord))
    fmt = 'Метод сеток с фиктивной точкой, k = {order:.3f}'
    plt.plot(log_h, err_fd_mx, '--', label = fmt.format(order = fdmx_ord))
    plt.legend()
    plt.show()

def main():
    def p(x):
        return -1 / x
    def q(x):
        return 1 / ( x * x )
    def f(x):
        return x
    def Solution(x):
        return x * np.log(x) + x ** 3 / 4 + x
    def SolutionDer(x):
        return 2 + np.log(x) + 0.75 * x ** 2
    a = 1.0
    b = 2.0
    y_a = Solution(a)
    y_b = Solution(b)
    dy_a = SolutionDer(a)
    PlotErrors(p, q, f, a, b, y_a, y_b, dy_a, Solution)

if __name__ == '__main__':
    main() 
