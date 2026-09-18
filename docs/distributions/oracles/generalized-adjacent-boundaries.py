"""Frozen adjacent Kappa support references; independent 400-digit Decimal formulas.

This extends the existing kappa-four boundary fixture from the fourth interior
binary64 value to the first, second, third and fourth values. No Numerics helper
is loaded. Binary64 shapes and evaluation arguments are converted exactly.
"""
import csv
import math
from decimal import Decimal as D, localcontext
from pathlib import Path

rows = []
with localcontext() as context:
    context.prec = 400
    for kf in [-2., -.5, 0., .2, 2., 10.]:
        for hf in [.2, .5, 1., 2., 10.]:
            k, h = D.from_float(kf), D.from_float(hf)
            lower = h.ln() if not k else (1 - (-k * h.ln()).exp()) / k
            endpoints = [(lower, math.inf, "lower")]
            if k > 0:
                endpoints.append((1 / k, -math.inf, "upper"))
            for endpoint, direction, label in endpoints:
                x_float = float(endpoint)
                for step in range(1, 5):
                    x_float = math.nextafter(x_float, direction)
                    x = D.from_float(x_float)
                    logt = -x if not k else (1 - k * x).ln() / k
                    logcdf = (1 - h * logt.exp()).ln() / h
                    logpdf = (1 - k) * logt + (1 - h) * logcdf
                    rows.append([kf, hf, label, step, x_float, float(logcdf), float(logpdf)])

target = Path(__file__).with_suffix(".csv")
with target.open("w", newline="", encoding="ascii") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(["kappa", "hondo", "endpoint", "step", "x", "logcdf", "logpdf"])
    writer.writerows(rows)
print(f"Wrote {len(rows)} independent adjacent-boundary rows to {target.name}")
