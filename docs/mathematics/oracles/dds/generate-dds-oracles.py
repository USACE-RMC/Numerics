"""Dynamically dimensioned search (DDS) oracles for Test_Numerics/Mathematics/Optimization/Global.

Two independent checks of Numerics.Mathematics.Optimization.DynamicallyDimensionedSearch:

golden
    A line-by-line Python transcription of Figure 1 of Tolson and Shoemaker (2007),
    doi:10.1029/2005WR004723: initial solution or best of M = max(5, floor(0.005 m)) uniform
    random solutions (paragraph 18), P(i) = 1 - ln(i)/ln(m) with i the evaluations so far,
    N(0, 1) perturbations with sigma = r (x_max - x_min), fully reflective bounds (step 4),
    greedy acceptance F(x_new) <= F_best (step 5), stop at i = m (step 6).
    It draws from a pure-Python MT19937 that replicates Numerics.Sampling.MersenneTwister
    (init_genrand seeding, NextDouble = genrand_int32 / 2^32, Next(n) by rejection), and
    scipy.special.ndtri for N(0, 1) by inversion, consuming draws in the same order as the
    C# solver: D selection uniforms, one Next(D) when the neighborhood is empty, then one
    uniform per selected variable in index order. The objectives use products only, so the
    fixtures printed to stdout pin the C# trajectory to within inverse-normal rounding.

spotpy OUT.json
    Runs spotpy 1.6.7 DDS over 30 seeds on the benchmark problems and budgets and writes
    the best value of every run. Note that spotpy implements the 2008 revision of DDS
    (50% reflect / 50% absorb at a bound) and computes P over the post-initialization
    budget, so it is a statistical comparison only.

Run: python generate-dds-oracles.py golden > fixtures.cs.txt
     python generate-dds-oracles.py spotpy spotpy-results.json
"""
import io
import json
import math
import sys
import contextlib

import numpy as np
from scipy.special import ndtri


class MT19937:
    """Replica of Numerics.Sampling.MersenneTwister(int seed)."""

    def __init__(self, seed):
        self.mt = [0] * 624
        self.mt[0] = seed & 0xFFFFFFFF
        for i in range(1, 624):
            self.mt[i] = (1812433253 * (self.mt[i - 1] ^ (self.mt[i - 1] >> 30)) + i) & 0xFFFFFFFF
        self.mti = 624

    def int32(self):
        if self.mti >= 624:
            mt = self.mt
            for kk in range(624):
                y = (mt[kk] & 0x80000000) | (mt[(kk + 1) % 624] & 0x7FFFFFFF)
                mt[kk] = mt[(kk + 397) % 624] ^ (y >> 1) ^ (0x9908B0DF if y & 1 else 0)
            self.mti = 0
        y = self.mt[self.mti]
        self.mti += 1
        y ^= y >> 11
        y ^= (y << 7) & 0x9D2C5680
        y ^= (y << 15) & 0xEFC60000
        y ^= y >> 18
        return y & 0xFFFFFFFF

    def next_double(self):
        return self.int32() * (1.0 / 4294967296.0)

    def next_int(self, n):
        threshold = 0xFFFFFFFF - (0xFFFFFFFF % n)
        while True:
            r = self.int32()
            if r < threshold:
                return r % n


def standard_z(u):
    return float(ndtri(5e-324 if u == 0.0 else u))


def reflect(v, lo, hi):
    """Step 4 of Figure 1."""
    if v < lo:
        v = lo + (lo - v)
        if v > hi:
            v = lo
    elif v > hi:
        v = hi - (v - hi)
        if v < lo:
            v = hi
    return v


def dds(f, lo, hi, m, seed, r=0.2, x0=None):
    """Reference DDS. Returns (best x, best f, best-so-far f after every evaluation)."""
    rng = MT19937(seed)
    D = len(lo)
    trace = []
    if x0 is not None:
        xb = list(x0)
        fb = f(xb)
        trace.append(fb)
        i = 1
    else:
        M = max(5, int(math.floor(0.005 * m)))
        xb, fb = None, math.inf
        for k in range(M):
            x = [lo[j] + (hi[j] - lo[j]) * rng.next_double() for j in range(D)]
            fx = f(x)
            if k == 0 or fx <= fb:
                xb, fb = x, fx
            trace.append(min(fb, trace[-1]) if trace else fb)
        i = M
    sigma = [r * (hi[j] - lo[j]) for j in range(D)]
    while i < m:
        P = 1.0 - math.log(i) / math.log(m)
        sel = [rng.next_double() < P for _ in range(D)]
        if not any(sel):
            sel[rng.next_int(D)] = True
        x = list(xb)
        for j in range(D):
            if sel[j]:
                x[j] = reflect(xb[j] + sigma[j] * standard_z(rng.next_double()), lo[j], hi[j])
        fx = f(x)
        i += 1
        if fx <= fb:
            xb, fb = x, fx
        trace.append(fb)
    return xb, fb, trace


def rosenbrock(x):
    s = 0.0
    for i in range(len(x) - 1):
        a = x[i + 1] - x[i] * x[i]
        b = 1 - x[i]
        s += 100 * a * a + b * b
    return s


def six_hump_camel(v):
    x, y = v[0], v[1]
    x2, y2 = x * x, y * y
    return (4 - 2.1 * x2 + x2 * x2 / 3) * x2 + x * y + (-4 + 4 * y2) * y2


def corner(x):
    """Optimum at the upper corner of [0, 1]^D; exercises the reflection rule."""
    s = 0.0
    for v in x:
        d = v - 1.0
        s += d * d
    return s


def golden():
    cases = [
        ("RosenbrockRandomStart", rosenbrock, [-2.0] * 3, [2.0] * 3, 300, 12345, 0.2, None),
        ("CamelInitialPoint", six_hump_camel, [-3.0, -2.0], [3.0, 2.0], 200, 7, 0.2, [2.5, -1.5]),
        ("CornerWideSteps", corner, [0.0] * 4, [1.0] * 4, 100, 99, 0.5, None),
    ]
    for name, f, lo, hi, m, seed, r, x0 in cases:
        xb, fb, trace = dds(f, lo, hi, m, seed, r, x0)
        checkpoints = [trace[k - 1] for k in (10, m // 2, m)]
        print(f"        // {name}: m={m}, seed={seed}, r={r}, x0={x0}")
        print(f"        private static readonly double[] {name}X = {{ {', '.join(repr(v) for v in xb)} }};")
        print(f"        private const double {name}F = {fb!r};")
        print(f"        private static readonly double[] {name}Trace = {{ {', '.join(repr(v) for v in checkpoints)} }};")


# Benchmark problems shared with the C# harness: name -> (f, lower, upper, f*).
def rastrigin(x):
    return 10.0 * len(x) + sum(v * v - 10.0 * math.cos(2 * math.pi * v) for v in x)


def griewank(x):
    s = sum(v * v for v in x) / 4000.0
    p = 1.0
    for i, v in enumerate(x):
        p *= math.cos(v / math.sqrt(i + 1))
    return s - p + 1.0


PROBLEMS = {
    "Rosenbrock": (rosenbrock, -2.0, 2.0, 0.0),
    "Rastrigin": (rastrigin, -5.12, 5.12, 0.0),
    "Griewank": (griewank, -600.0, 600.0, 0.0),
    "SixHumpCamel": (six_hump_camel, None, None, -1.031628453489877),
}


def spotpy_runs(out_path):
    import spotpy
    from spotpy.parameter import Uniform

    class Setup:
        def __init__(self, f, lo, hi):
            self.f = f
            self.params = [Uniform(f"x{j}", lo[j], hi[j]) for j in range(len(lo))]

        def parameters(self):
            return spotpy.parameter.generate(self.params)

        def simulation(self, x):
            return [self.f(list(x))]

        def evaluation(self):
            return [0.0]

        def objectivefunction(self, simulation, evaluation):
            return -simulation[0]  # spotpy DDS maximizes

    results = []
    for name, (f, lo, hi, fstar) in PROBLEMS.items():
        dims = [2] if name == "SixHumpCamel" else [5, 10]
        for D in dims:
            if name == "SixHumpCamel":
                lob, hib = [-3.0, -2.0], [3.0, 2.0]
            else:
                lob, hib = [lo] * D, [hi] * D
            for m in (1000, 5000):
                for seed in range(1, 31):
                    sampler = spotpy.algorithms.dds(Setup(f, lob, hib), dbformat="ram", save_sim=False, random_state=seed)
                    with contextlib.redirect_stdout(io.StringIO()):
                        sampler.sample(m)
                    best = -float(np.max(sampler.getdata()["like1"]))
                    results.append({"problem": name, "D": D, "m": m, "seed": seed, "best": best, "fstar": fstar})
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump({"spotpy": spotpy.__version__, "numpy": np.__version__, "results": results}, fh, indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "golden":
        golden()
    elif sys.argv[1] == "spotpy":
        spotpy_runs(sys.argv[2])
    else:
        raise SystemExit(__doc__)
