"""pycma side of the Numerics CMAES comparison.

Runs pycma (IPOP restarts, BoundTransform bounds) on the same problems as the
Numerics harness, in the same unit-cube scaling: the objective is evaluated at
lo + (hi - lo) * u with u in [0, 1]^n, x0 = box center for the first run and a
uniform random point for each restart, sigma0 = 0.3, popsize doubling.

For each seed it records the number of evaluations until the first objective
value below 1e-8 (or -1), the best value, and the evaluations used.
A timing mode measures the wall time per evaluation on a cheap objective with
restarts disabled.

Usage:
    python cmaes_pycma_comparison.py bench <problem> <n> <budget> <seeds> <out.csv>
    python cmaes_pycma_comparison.py time <n> <evals> <reps>
"""
import csv
import math
import sys
import time

import cma
import numpy as np


def rotation(n):
    """Deterministic orthogonal matrix: product of Givens rotations G(p, q, 0.5 + p*q), p < q."""
    a = np.eye(n)
    for p in range(n):
        for q in range(p + 1, n):
            t = 0.5 + p * q
            c, s = math.cos(t), math.sin(t)
            ap, aq = a[p].copy(), a[q].copy()
            a[p] = c * ap - s * aq
            a[q] = s * ap + c * aq
    return a


def problem(name, n):
    """Return (f, lo, hi) matching the Numerics harness."""
    if name == "sphere":
        return (lambda x: float(np.sum((x - 1) ** 2))), -5.0, 5.0
    if name in ("ellipsoid", "ellipsoid_axis"):
        r = rotation(n) if name == "ellipsoid" else np.eye(n)
        c = 10.0 ** (6.0 * np.arange(n) / (n - 1))
        return (lambda x: float(np.sum(c * (r @ (x - 1)) ** 2))), -5.0, 5.0
    if name == "rosenbrock":
        return (lambda x: float(np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2))), -5.0, 5.0
    if name == "rastrigin":
        return (lambda x: float(10 * n + np.sum((x - 1) ** 2 - 10 * np.cos(2 * np.pi * (x - 1))))), -5.12, 5.12
    raise ValueError(name)


def bench(name, n, budget, seeds, out):
    f, lo, hi = problem(name, n)
    rows = []
    for seed in range(1, seeds + 1):
        state = {"count": 0, "hit": -1, "best": math.inf}

        def g(u):
            state["count"] += 1
            v = f(lo + (hi - lo) * np.asarray(u))
            if state["hit"] < 0 and v < 1e-8:
                state["hit"] = state["count"]
            state["best"] = min(state["best"], v)
            return v

        rng = np.random.default_rng(seed)
        runs = {"k": 0}

        def x0():
            runs["k"] += 1
            return 0.5 * np.ones(n) if runs["k"] == 1 else rng.random(n)

        opts = {"bounds": [0, 1], "seed": seed, "verbose": -9, "verb_log": 0,
                "maxfevals": budget, "ftarget": 1e-8, "tolfun": 1e-12}
        cma.fmin(g, x0, 0.3, opts, restarts=9, incpopsize=2)
        rows.append([seed, state["hit"], repr(state["best"]), state["count"]])
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["seed", "hit", "fbest", "fe"])
        w.writerows(rows)


def timing(n, evals, reps):
    lo, hi = -1e3, 2e3
    weights = np.arange(1, n + 1)
    best = math.inf
    fe = 0
    for _ in range(reps):
        state = {"count": 0}

        def g(u):
            state["count"] += 1
            x = lo + (hi - lo) * np.asarray(u)
            return float(np.sum(weights * x * x))

        opts = {"bounds": [0, 1], "seed": 1, "verbose": -9, "verb_log": 0,
                "maxfevals": evals, "tolfun": 1e-15, "tolx": 1e-15}
        t = time.perf_counter()
        cma.fmin(g, 0.5 * np.ones(n), 0.3, opts)
        dt = time.perf_counter() - t
        if dt < best:
            best, fe = dt, state["count"]
    print(f"n={n},fe={fe},seconds={best:.4f},us_per_eval={1e6 * best / fe:.2f}")


if __name__ == "__main__":
    if sys.argv[1] == "bench":
        bench(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])
    else:
        timing(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
