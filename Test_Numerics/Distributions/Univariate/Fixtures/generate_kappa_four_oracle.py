"""Regenerate the frozen Kappa Four probability oracle using Python's 90-digit Decimal.
Formula provenance: SciPy scipy.stats.kappa4 and Hosking's Kappa quantile.
No Numerics code, NumPy, SciPy installation, or runtime network is used.
The represented binary64 shapes/probabilities/arguments are converted exactly to Decimal.
Rows whose rounded quantile is outside the mathematical open support are omitted.
"""
from decimal import Decimal as D, localcontext
from pathlib import Path
import csv
import math

target = Path(__file__).with_name("kappa-four-probabilities.csv")
shapes = [-2., -1., -.5, -.2, -1e-8, -1e-12, -1e-16, 0., 1e-16, 1e-12, 1e-8, .2, .5, 1., 2.]
probabilities = [1e-12, .01, .1, .5, .9, .99, .999999999999]
rows = []

def quantile(k, h, p):
    a = -p.ln() if not h else (1 - (h*p.ln()).exp()) / h
    return -a.ln() if not k else (1 - (k*a.ln()).exp()) / k

with localcontext() as context:
    context.prec = 90
    for kf in shapes:
        for hf in shapes:
            k, h = D.from_float(kf), D.from_float(hf)
            for pf in probabilities:
                p = D.from_float(pf)
                q = quantile(k, h, p)
                x = D.from_float(float(q))
                base = 1 - k*x
                if base <= 0:
                    continue
                logt = -x if not k else base.ln()/k
                t = logt.exp()
                b = 1 - h*t
                if b <= 0:
                    continue
                logf = -t if not h else b.ln()/h
                logpdf = (1-k)*logt + (1-h)*logf
                # High-precision symmetric differences independently check production derivatives.
                step = D('1e-25')
                dk = (quantile(k+step, h, p) - quantile(k-step, h, p))/(2*step)
                dh = (quantile(k, h+step, p) - quantile(k, h-step, p))/(2*step)
                rows.append([kf, hf, pf, float(x), float(q), float(logf.exp()), float(logpdf.exp()), float(dk), float(dh)])
with target.open("w", newline="", encoding="ascii") as output:
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["kappa", "hondo", "probability", "x", "quantile", "cdf", "pdf", "gradient_kappa", "gradient_hondo"])
    writer.writerows(rows)
print(f"Wrote {len(rows)} cases to {target.name}")

boundary_rows = []
with localcontext() as context:
    # Extra digits retain the subnormal distance above the exact zero boundary at h=1.
    context.prec = 400
    for kf in [-2., -.5, 0., .2, 2., 10.]:
        for hf in [.2, .5, 1., 2., 10.]:
            k, h = D.from_float(kf), D.from_float(hf)
            lower = h.ln() if not k else (1-(-k*h.ln()).exp())/k
            endpoints = [(lower, math.inf)]
            if k > 0:
                endpoints.append((1/k, -math.inf))
            for endpoint, direction in endpoints:
                xf = float(endpoint)
                for unused in range(4):
                    xf = math.nextafter(xf, direction)
                x = D.from_float(xf)
                logt = -x if not k else (1-k*x).ln()/k
                logcdf = (1-h*logt.exp()).ln()/h
                logpdf = (1-k)*logt+(1-h)*logcdf
                boundary_rows.append([kf, hf, xf, float(logcdf), float(logpdf)])
boundary_target = target.with_name("kappa-four-boundaries.csv")
with boundary_target.open("w", newline="", encoding="ascii") as output:
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["kappa", "hondo", "x", "logcdf", "logpdf"])
    writer.writerows(boundary_rows)
print(f"Wrote {len(boundary_rows)} boundary cases to {boundary_target.name}")
