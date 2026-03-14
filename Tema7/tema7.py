import random
import math

def horner(coeff, x):
    b = coeff[0]
    for a in coeff[1:]:
        b = b * x + a
    return b

def derivative(coeff):
    n = len(coeff) - 1
    return [coeff[i] * (n - i) for i in range(n)]

def second_derivative(coeff):
    return derivative(derivative(coeff))

def compute_R(coeff):
    a0 = coeff[0]
    A = max(abs(a) for a in coeff[1:])
    R = (abs(a0) + A) / abs(a0)
    return R

def newton_method(coeff, x0, eps, kmax=1000):
    d1 = derivative(coeff)

    x = x0
    k = 0

    while k < kmax:
        Px = horner(coeff, x)
        Pdx = horner(d1, x)

        if abs(Pdx) < eps:
            return None, k

        dx = Px / Pdx
        x = x - dx

        if abs(dx) < eps:
            return x, k

        k += 1

    return None, k

def olver_method(coeff, x0, eps, kmax=1000):
    d1 = derivative(coeff)
    d2 = second_derivative(coeff)

    x = x0
    k = 0

    while k < kmax:

        Px = horner(coeff, x)
        Pdx = horner(d1, x)
        Pd2x = horner(d2, x)

        if abs(Pdx) < eps:
            return None, k

        c = (Px**2 * Pd2x) / (Pdx**3)

        dx = Px / Pdx + 0.5 * c

        x = x - dx

        if abs(dx) < eps:
            return x, k

        k += 1

    return None, k

def distinct_roots(roots, eps):
    unique = []

    for r in roots:
        ok = True
        for u in unique:
            if abs(r - u) <= eps:
                ok = False
                break
        if ok:
            unique.append(r)

    return unique

def main():

    coeff = [1, -6, 11, -6]

    eps = 1e-6

    R = compute_R(coeff)

    print("Interval radacini reale:", -R, R)

    roots_newton = []
    roots_olver = []

    steps_newton = []
    steps_olver = []

    for _ in range(50):

        x0 = random.uniform(-R, R)

        rN, kN = newton_method(coeff, x0, eps)
        rO, kO = olver_method(coeff, x0, eps)

        if rN is not None:
            roots_newton.append(rN)
            steps_newton.append(kN)

        if rO is not None:
            roots_olver.append(rO)
            steps_olver.append(kO)

    roots_newton = distinct_roots(roots_newton, eps)
    roots_olver = distinct_roots(roots_olver, eps)

    print("\nRadacini Newton:")
    for r in roots_newton:
        print(r)

    print("\nRadacini Olver:")
    for r in roots_olver:
        print(r)

    if steps_newton:
        print("\nPasi medii Newton:", sum(steps_newton)/len(steps_newton))

    if steps_olver:
        print("Pasi medii Olver:", sum(steps_olver)/len(steps_olver))

    with open("radacini.txt", "w") as f: #trebuie pusa toata calea

        f.write("Radacini distincte:\n")

        roots = distinct_roots(roots_newton + roots_olver, eps)

        for r in roots:
            f.write(str(r) + "\n")

main()
