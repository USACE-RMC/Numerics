"""Frozen generalized Pareto survival-quantile references for GeneralizedPareto.InverseCCDF.

Inputs: numpy.random.default_rng(141); per case, in order, location uniform(0, 5), scale
uniform(0.1, 3), hydrologic shape choice([0, uniform(-0.4, 0.6)]) and survival probability
10**uniform(-12, 0). The Hosking shape is kappa = -shape (zero stays +0).

Columns, all binary64 values as big-endian hexadecimal bit patterns:
- location, scale, kappa, q: the inputs;
- reference: the binary64 expression u + sigma*(expm1(xi*L)/xi) with L = -log(q) (u + sigma*L at
  xi = 0), using the algorithm of Numerics Tools.Expm1: the direct difference exp(x) - 1 for
  |x| > 1/2, and otherwise the compensated expm1 of Kahan, u = exp(x), (u - 1) x / log(u);
- exact: the 200-bit mpmath value of u + sigma*expm1(xi*L)/xi, rounded to binary64.
"""
import csv
import math
import struct
from pathlib import Path

import mpmath
import numpy as np

mpmath.mp.prec = 200
def expm1(x):
    """exp(x) - 1 as Tools.Expm1 evaluates it: the direct difference for |x| > 1/2, else Kahan's compensated form."""
    if abs(x) > 0.5:
        return (math.exp(x) if x < 709.78 else math.inf) - 1.0
    u = math.exp(x)
    if u == 1.0:
        return x
    return (u - 1.0) * x / math.log(u)


def bits(value):
    return struct.pack(">d", value).hex().upper()


rng = np.random.default_rng(141)
rows = []
for _ in range(300):
    u = float(rng.uniform(0, 5))
    sigma = float(rng.uniform(0.1, 3))
    xi = float(rng.choice([0, rng.uniform(-0.4, 0.6)]))
    q = float(10 ** rng.uniform(-12, 0))
    log = -math.log(q)
    reference = u + sigma * log if xi == 0 else u + sigma * (expm1(xi * log) / xi)
    big_log = -mpmath.log(mpmath.mpf(q))
    exact = mpmath.mpf(u) + mpmath.mpf(sigma) * (big_log if xi == 0 else mpmath.expm1(mpmath.mpf(xi) * big_log) / mpmath.mpf(xi))
    kappa = -xi if xi != 0 else 0.0
    rows.append([bits(u), bits(sigma), bits(kappa), bits(q), bits(reference), bits(float(exact))])

with open(Path(__file__).with_suffix(".csv"), "w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["location", "scale", "kappa", "q", "reference", "exact"])
    writer.writerows(rows)
