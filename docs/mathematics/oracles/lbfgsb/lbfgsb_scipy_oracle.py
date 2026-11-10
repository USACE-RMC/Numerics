"""SciPy L-BFGS-B reference results for the Numerics LBFGSB tests.

Runs scipy.optimize.minimize(method="L-BFGS-B") with analytic gradients, maxcor=10, gtol=1e-5 and
ftol=2.220446049250313e-09 (the defaults, factr = 1e7) on the problems below, and writes the optimum,
objective, iteration and evaluation counts, active set and final projected-gradient norm to
lbfgsb-scipy.json next to this script.

Usage: python lbfgsb_scipy_oracle.py
"""
import json
import math
import os
import platform

import numpy as np
import scipy
from scipy.optimize import minimize

INF = math.inf


def rosenbrock(x):
    return sum(100.0 * (x[i + 1] - x[i] ** 2) ** 2 + (1.0 - x[i]) ** 2 for i in range(len(x) - 1))


def rosenbrock_grad(x):
    g = np.zeros(len(x))
    for i in range(len(x) - 1):
        a = x[i + 1] - x[i] ** 2
        g[i] += -400.0 * x[i] * a - 2.0 * (1.0 - x[i])
        g[i + 1] += 200.0 * a
    return g


def fxyz(p):
    x, y, z = p
    return (4 * x - 0.5) ** 2 + (3 * y - 0.6) ** 2 + (2 * z - 0.7) ** 2


def fxyz_grad(p):
    x, y, z = p
    return np.array([8 * (4 * x - 0.5), 6 * (3 * y - 0.6), 4 * (2 * z - 0.7)])


def dejong(x):
    return float(np.sum(np.asarray(x) ** 2))


def dejong_grad(x):
    return 2.0 * np.asarray(x)


def sum_of_powers(x):
    return sum(abs(x[i]) ** (i + 2) for i in range(len(x)))


def sum_of_powers_grad(x):
    return np.array([(i + 2) * abs(x[i]) ** (i + 1) * np.sign(x[i]) for i in range(len(x))])


def booth(p):
    x, y = p
    return (x + 2 * y - 7) ** 2 + (2 * x + y - 5) ** 2


def booth_grad(p):
    x, y = p
    a, b = x + 2 * y - 7, 2 * x + y - 5
    return np.array([2 * a + 4 * b, 4 * a + 2 * b])


def matyas(p):
    x, y = p
    return 0.26 * (x * x + y * y) - 0.48 * x * y


def matyas_grad(p):
    x, y = p
    return np.array([0.52 * x - 0.48 * y, 0.52 * y - 0.48 * x])


def mccormick(p):
    x, y = p
    return math.sin(x + y) + (x - y) ** 2 - 1.5 * x + 2.5 * y + 1.0


def mccormick_grad(p):
    x, y = p
    c = math.cos(x + y)
    return np.array([c + 2 * (x - y) - 1.5, c - 2 * (x - y) + 2.5])


def beale(p):
    x, y = p
    return (1.5 - x + x * y) ** 2 + (2.25 - x + x * y * y) ** 2 + (2.625 - x + x * y ** 3) ** 2


def beale_grad(p):
    x, y = p
    a, b, c = 1.5 - x + x * y, 2.25 - x + x * y * y, 2.625 - x + x * y ** 3
    return np.array([2 * a * (y - 1) + 2 * b * (y * y - 1) + 2 * c * (y ** 3 - 1),
                     2 * a * x + 4 * b * x * y + 6 * c * x * y * y])


def corner(p):
    return (p[0] - 2.0) ** 2 + (p[1] + 2.0) ** 2


def corner_grad(p):
    return np.array([2 * (p[0] - 2.0), 2 * (p[1] + 2.0)])


W = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
C = np.array([-1.0, 0.5, 2.0, 0.3, -0.2, 1.5])


def weighted_quadratic(x):
    return float(np.sum(W * (np.asarray(x) - C) ** 2))


def weighted_quadratic_grad(x):
    return 2.0 * W * (np.asarray(x) - C)


A = np.array([3.0, -2.0, 1.0, -4.0])


def tilted(x):
    x = np.asarray(x)
    return float(A @ x + 0.05 * x @ x)


def tilted_grad(x):
    return A + 0.1 * np.asarray(x)


def linear(x):
    return -x[0]


def linear_grad(x):
    return np.array([-1.0])


def rosen_start(n):
    return [-1.2 if i % 2 == 0 else 1.0 for i in range(n)]


PROBLEMS = [
    ("FXYZ", fxyz, fxyz_grad, [0.2, 0.5, 0.5], [0.0] * 3, [1.0] * 3),
    ("DeJong5", dejong, dejong_grad, [1.0, -1.0, 2.0, -2.0, 1.0], [-5.12] * 5, [5.12] * 5),
    ("SumOfPowers3", sum_of_powers, sum_of_powers_grad, [0.5, -0.5, 0.5], [-1.0] * 3, [1.0] * 3),
    ("Rosenbrock2Box", rosenbrock, rosenbrock_grad, [0.0, 0.0], [-2.048] * 2, [2.048] * 2),
    ("Booth", booth, booth_grad, [0.0, 0.0], [-10.0] * 2, [10.0] * 2),
    ("Matyas", matyas, matyas_grad, [1.0, -1.0], [-10.0] * 2, [10.0] * 2),
    ("McCormick", mccormick, mccormick_grad, [0.0, 0.0], [-1.5, -3.0], [4.0, 4.0]),
    ("Beale", beale, beale_grad, [0.0, 0.0], [-4.5] * 2, [4.5] * 2),
    ("Corner", corner, corner_grad, [0.5, 0.5], [0.0] * 2, [1.0] * 2),
    ("Linear", linear, linear_grad, [0.0], [0.0], [1.0]),
    ("Rosenbrock2", rosenbrock, rosenbrock_grad, rosen_start(2), [-INF] * 2, [INF] * 2),
    ("Rosenbrock5", rosenbrock, rosenbrock_grad, rosen_start(5), [-INF] * 5, [INF] * 5),
    ("Rosenbrock10", rosenbrock, rosenbrock_grad, rosen_start(10), [-INF] * 10, [INF] * 10),
    ("WeightedQuadraticActive", weighted_quadratic, weighted_quadratic_grad, [0.5] * 6, [0.0] * 6, [1.0] * 6),
    ("TiltedAllBounds", tilted, tilted_grad, [0.0] * 4, [-1.0] * 4, [1.0] * 4),
    ("Rosenbrock4LowerOnly", rosenbrock, rosenbrock_grad, rosen_start(4), [-INF, 1.2, -INF, -INF], [INF] * 4),
    ("Rosenbrock4UpperOnly", rosenbrock, rosenbrock_grad, rosen_start(4), [-INF] * 4, [INF, 0.8, INF, INF]),
    ("Rosenbrock4Fixed", rosenbrock, rosenbrock_grad, rosen_start(4), [-INF, 0.5, -INF, -INF], [INF, 0.5, INF, INF]),
]


def projected_gradient_norm(x, g, lower, upper):
    norm = 0.0
    for xi, gi, li, ui in zip(x, g, lower, upper):
        if gi < 0 and math.isfinite(ui):
            gi = max(xi - ui, gi)
        elif gi >= 0 and math.isfinite(li):
            gi = min(xi - li, gi)
        norm = max(norm, abs(gi))
    return norm


def main():
    results = []
    for name, f, g, x0, lower, upper in PROBLEMS:
        res = minimize(f, np.array(x0, dtype=float), jac=g, method="L-BFGS-B",
                       bounds=list(zip(lower, upper)),
                       options={"maxcor": 10, "gtol": 1e-5, "ftol": 2.220446049250313e-09,
                                "maxiter": 15000, "maxfun": 15000, "maxls": 20})
        x = [float(v) for v in res.x]
        active = [i for i in range(len(x))
                  if (math.isfinite(lower[i]) and x[i] == lower[i]) or (math.isfinite(upper[i]) and x[i] == upper[i])]
        results.append({
            "name": name,
            "x0": x0,
            "lower": [None if math.isinf(v) else v for v in lower],
            "upper": [None if math.isinf(v) else v for v in upper],
            "x": x,
            "fun": float(res.fun),
            "nit": int(res.nit),
            "nfev": int(res.nfev),
            "success": bool(res.success),
            "message": str(res.message),
            "active": active,
            "projectedGradientNorm": projected_gradient_norm(x, g(res.x), lower, upper),
        })
    out = {
        "reference": "scipy.optimize.minimize(method='L-BFGS-B', maxcor=10, gtol=1e-5, ftol=2.220446049250313e-09, maxls=20)",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "problems": results,
    }
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lbfgsb-scipy.json")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
