"""Wall-time reference for Benchmarks_Numerics/Optimization/LBFGSBBenchmarks.cs.

Times scipy.optimize.minimize(method="L-BFGS-B") with analytic NumPy gradients and the default options
(maxcor=10, gtol=1e-5, ftol=2.220446049250313e-09) on the same four problems as the C# benchmark, and prints
the median of the repetitions with the iteration and evaluation counts. The results are written to
lbfgsb-scipy-timings.json next to this script.

Usage: python lbfgsb_scipy_timing.py [repetitions]
"""
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
from scipy.optimize import minimize


def rosenbrock(x):
    a = x[1:] - x[:-1] ** 2
    return float(np.sum(100.0 * a * a + (1.0 - x[:-1]) ** 2))


def rosenbrock_grad(x):
    a = x[1:] - x[:-1] ** 2
    g = np.zeros_like(x)
    g[:-1] = -400.0 * x[:-1] * a - 2.0 * (1.0 - x[:-1])
    g[1:] += 200.0 * a
    return g


def problem(name):
    if name == "Quadratic1000":
        n = 1000
        i = np.arange(n)
        w = 1.0 + i % 10
        c = 1.5 * np.sin(i)
        f = lambda x: float(np.sum(w * (x - c) ** 2))
        g = lambda x: 2.0 * w * (x - c)
        return f, g, np.zeros(n), [(-1.0, 1.0)] * n
    n = 10 if name == "Rosenbrock10" else 100
    x0 = np.array([-1.2 if i % 2 == 0 else 1.0 for i in range(n)])
    bounds = [(-2.0, 0.9)] * n if name == "BoxRosenbrock100" else None
    return rosenbrock, rosenbrock_grad, x0, bounds


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    rows = []
    for name in ["Rosenbrock10", "Rosenbrock100", "BoxRosenbrock100", "Quadratic1000"]:
        f, g, x0, bounds = problem(name)
        times = []
        res = None
        for _ in range(reps + 1):
            start = time.perf_counter()
            res = minimize(f, x0, jac=g, method="L-BFGS-B", bounds=bounds)
            times.append(time.perf_counter() - start)
        median = float(np.median(times[1:]))
        rows.append({"problem": name, "seconds": median, "nit": int(res.nit), "nfev": int(res.nfev), "fun": float(res.fun)})
        print(f"{name:18s} {median * 1e3:9.3f} ms  nit {res.nit:4d}  nfev {res.nfev:4d}  f {res.fun:.6g}")
    out = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
           "repetitions": reps, "results": rows}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lbfgsb-scipy-timings.json")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
