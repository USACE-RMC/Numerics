"""Timing of scikit-learn GaussianProcessRegressor for the GaussianProcessRegression benchmark.

The data are deterministic and identical in the C#, Python and R benchmarks:
    x[i, j] = frac((i + 1) * sqrt(p_j)),  p_j the j-th prime (2, 3, 5, ...),
    y[i]    = sum_j sin(2 pi x[i, j]) / (j + 1) + 0.1 * (frac((i + 1) * sqrt(97)) - 0.5),
for training rows i = 0..n-1; the m = 1000 prediction points use rows n..n+m-1.

Fit: Matern 5/2 with one length scale per feature times a constant, plus a white-noise term,
one optimizer start (n_restarts_optimizer=0, L-BFGS-B), normalize_y=False.
Predict: mean and standard deviation at the m points.
The median wall time over the repetitions is recorded for each size.

Usage:
    python gpr_benchmark_sklearn.py <d> <out.json> n:reps [n:reps ...]
"""
import json
import os
import platform
import sys
import time

import numpy as np
import scipy
import sklearn
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]


def data(rows, d, offset=0):
    i = np.arange(offset + 1, offset + rows + 1, dtype=float)[:, None]
    v = i * np.sqrt(np.array(PRIMES[:d], dtype=float))[None, :]
    x = v - np.floor(v)
    e = i[:, 0] * np.sqrt(97.0)
    y = (np.sin(2 * np.pi * x) / np.arange(1, d + 1)).sum(axis=1) + 0.1 * (e - np.floor(e) - 0.5)
    return x, y


def main():
    d, out = int(sys.argv[1]), sys.argv[2]
    sizes = [tuple(int(v) for v in a.split(":")) for a in sys.argv[3:]]
    blas = np.show_config(mode="dicts")["Build Dependencies"]["blas"]
    result = {"library": "scikit-learn", "sklearn": sklearn.__version__, "numpy": np.__version__,
              "scipy": scipy.__version__, "python": platform.python_version(),
              "blas": blas.get("name", "") + " " + blas.get("version", ""),
              "processor": platform.processor(), "cpuCount": os.cpu_count(), "d": d, "m": 1000, "sizes": []}
    for n, reps in sizes:
        x, y = data(n, d)
        xt, _ = data(1000, d, n)
        fit, pred, lml = [], [], None
        for _ in range(reps):
            kernel = ConstantKernel(1.0, (1e-5, 1e5)) * Matern(np.ones(d), (1e-3, 1e3), nu=2.5) + WhiteKernel(0.01, (1e-8, 1e2))
            gp = GaussianProcessRegressor(kernel=kernel, alpha=0.0, n_restarts_optimizer=0)
            t0 = time.perf_counter()
            gp.fit(x, y)
            fit.append(time.perf_counter() - t0)
            t0 = time.perf_counter()
            gp.predict(xt, return_std=True)
            pred.append(time.perf_counter() - t0)
            lml = float(gp.log_marginal_likelihood_value_)
        result["sizes"].append({"n": n, "reps": reps, "fitSeconds": float(np.median(fit)),
                                "predictSeconds": float(np.median(pred)), "logMarginalLikelihood": lml})
        print(result["sizes"][-1], flush=True)
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=1)


if __name__ == "__main__":
    main()
