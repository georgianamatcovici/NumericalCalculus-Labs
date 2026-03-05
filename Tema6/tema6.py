import numpy as np
import matplotlib.pyplot as plt

# functia
def f(x):
    return 1 / (1 + x**2)
# derivata
def df(x):
    return -2*x / (1 + x**2)**2

# Horner: evaluare polinom
# coef = [a0, a1, ..., am] (puteri crescatoare)

def horner_eval(coef: np.ndarray, x: float) -> float:
    # P(x) = a0 + a1 x + ... + am x^m
    # Horner: (((am)x + a_{m-1})x + ... + a0)
    p = 0.0
    for a in coef[::-1]:
        p = p * x + a
    return p

# Least Squares polinomial 
# B_ij = sum_k x_k^(i+j)
# f_i  = sum_k y_k x_k^i
def least_squares_poly_coeffs(x: np.ndarray, y: np.ndarray, m: int) -> np.ndarray:

    n = len(x) - 1

    S = np.array([np.sum(x ** p) for p in range(2 * m + 1)], dtype=float)

    B = np.zeros((m + 1, m + 1), dtype=float)
    for i in range(m + 1):
        for j in range(m + 1):
            B[i, j] = S[i + j] # construieste matricea B

    fvec = np.array([np.sum(y * (x ** i)) for i in range(m + 1)], dtype=float) # construieste vectorul f

    a = np.linalg.solve(B, fvec)  # rezolva sistemul
    return a

# Spline cubic clamped C^2 (derivate la capete date)
# - h_i = x_{i+1}-x_i
# - sistem HA = rhs (A0..An)
# - apoi b_i, c_i si evaluare pe interval
def clamped_cubic_spline_setup(x, y, da, db):
    n = len(x) - 1
    h = np.diff(x)

    H = np.zeros((n+1, n+1), float)
    rhs = np.zeros(n+1, float)

    # clamped left
    H[0,0] = 2*h[0]
    H[0,1] = h[0]
    rhs[0] = 6*((y[1]-y[0])/h[0] - da)

    # interior
    for i in range(1, n):
        H[i,i-1] = h[i-1]
        H[i,i]   = 2*(h[i-1] + h[i])
        H[i,i+1] = h[i]
        rhs[i]   = 6*((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])

    # clamped right
    H[n,n-1] = h[n-1]
    H[n,n]   = 2*h[n-1]
    rhs[n]   = 6*(db - (y[n]-y[n-1])/h[n-1])

    M = np.linalg.solve(H, rhs)   # M_i = S''(x_i)
    return h, M

def clamped_cubic_spline_eval(x, y, h, M, xbar):
    i = np.searchsorted(x, xbar) - 1
    i = max(0, min(i, len(x) - 2))

    xi, xip1 = x[i], x[i+1]
    hi = h[i]

    term1 = M[i]   * (xip1 - xbar)**3 / (6*hi)
    term2 = M[i+1] * (xbar - xi)**3   / (6*hi)
    term3 = (y[i]   - M[i]*hi*hi/6)   * (xip1 - xbar) / hi
    term4 = (y[i+1] - M[i+1]*hi*hi/6) * (xbar - xi)   / hi

    return term1 + term2 + term3 + term4

def main():
    np.random.seed(0)

    x0 = float(input("x0 = "))
    xn = float(input("xn = "))
    if x0 >= xn:
        raise ValueError("Trebuie x0 < xn!")

    n = int(input("n = "))
    if n < 2:
        raise ValueError("Ia n >= 2.")

    xbar = float(input("x_bar (sa NU fie egal cu vreun xi) = "))
    if not (x0 < xbar < xn):
        print("Atentie: ideal x_bar in (x0, xn)!")

    # m < 6
    m_list_str = input("Lista grade m (<6), ex: 1 2 3 4 5 = ").strip()
    m_list = [int(s) for s in m_list_str.split()]
    for m in m_list:
        if m >= 6 or m < 0:
            raise ValueError("Toate m trebuie sa fie in {0,1,2,3,4,5}.")

    # x1..x_{n-1} random in (x0,xn), sortate
    interior = np.random.uniform(x0, xn, size=n - 1)
    interior.sort()
    x = np.concatenate(([x0], interior, [xn]))
    y = f(x)

    # derivate la capete (caz clamped)
    da = df(x0)
    db = df(xn)

    true_f_xbar = float(f(np.array([xbar]))[0])

    print("\nNoduri (xi) si valori (yi):")
    for i in range(n + 1):
        print(f"i={i:2d}: x={x[i]: .6f}, y={y[i]: .6f}")

    # LS polinomial pentru fiecare m
    print("\nAproximare LS (cele mai mici patrate) + Horner")
    ls_results = {}
    for m in m_list:
        a = least_squares_poly_coeffs(x, y, m)  # a0..am
        P_xbar = horner_eval(a, xbar)
        err_xbar = abs(P_xbar - true_f_xbar)

        # suma |Pm(xi) - yi|
        P_at_nodes = np.array([horner_eval(a, xi) for xi in x])
        sum_abs = float(np.sum(np.abs(P_at_nodes - y)))

        ls_results[m] = (a, P_xbar)

        print(f"\nm = {m}")
        print("coef [a0..am] =", np.array2string(a, precision=6, suppress_small=True))
        print(f"Pm(x_bar) = {P_xbar:.10f}")
        print(f"|Pm(x_bar) - f(x_bar)| = {err_xbar:.10e}")
        print(f"Sum_i |Pm(xi) - yi| = {sum_abs:.10e}")

    # spline cubic C^2 clamped 
    print("\nSpline cubic C^2 (clamped)")
    h, M = clamped_cubic_spline_setup(x, y, da, db)
    S_xbar = clamped_cubic_spline_eval(x, y, h, M, xbar)
    err_s = abs(S_xbar - true_f_xbar)
    print(f"Sf(x_bar) = {S_xbar:.10f}")
    print(f"|Sf(x_bar) - f(x_bar)| = {err_s:.10e}")

    # grafice 
    grid = np.linspace(x0, xn, 600)
    f_grid = f(grid)

    plt.figure()
    plt.plot(grid, f_grid, label="f(x)")
    plt.scatter(x, y, marker="o", label="Noduri (xi, yi)")

    # plot pentru fiecare m
    for m in m_list:
        a, _ = ls_results[m]
        P_grid = np.array([horner_eval(a, t) for t in grid])
        plt.plot(grid, P_grid, label=f"P_m, m={m}")

    # spline
    S_grid = np.array([clamped_cubic_spline_eval(x, y, h, M, t) for t in grid])
    plt.plot(grid, S_grid, label="Spline C^2")

    plt.axvline(xbar, linestyle="--", label="x_bar")
    plt.title("f(x) vs aproximari (LS & Spline)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
