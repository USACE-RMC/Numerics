"""Frozen generalized extreme value survival-quantile references for GeneralizedExtremeValue.InverseCCDF.

Grid: location 10, scale 2, nine Hosking shapes kappa in {-0.5, -0.3, -0.1, -1e-6, 0, 1e-6, 0.1, 0.3,
0.5} and survival probabilities q = 10**(-e/20) for e = 1..240 plus 1e-20, 1e-50, 1e-100 and 1e-300.
Binary64 inputs are converted exactly. The exact quantile is the 400-bit mpmath value of
10 + 2 (1 - exp(kappa g)) / kappa (10 - 2 g at kappa = 0) with g = log(-log1p(-q)), rounded to binary64.

Columns, all binary64 values as big-endian hexadecimal bit patterns: kappa, q, exact.
"""
import csv
import struct
from pathlib import Path

import mpmath

mpmath.mp.prec = 400


def bits(value):
    return struct.pack(">d", value).hex().upper()


rows = []
for kappa in [-0.5, -0.3, -0.1, -1e-6, 0.0, 1e-6, 0.1, 0.3, 0.5]:
    probabilities = [10 ** (-e / 20) for e in range(1, 241)] + [1e-20, 1e-50, 1e-100, 1e-300]
    for q in probabilities:
        g = mpmath.log(-mpmath.log1p(-mpmath.mpf(q)))
        k = mpmath.mpf(kappa)
        exact = 10 - 2 * g if kappa == 0 else 10 + 2 * (1 - mpmath.exp(k * g)) / k
        rows.append([bits(kappa), bits(q), bits(float(exact))])

with open(Path(__file__).with_suffix(".csv"), "w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["kappa", "q", "exact"])
    writer.writerows(rows)
