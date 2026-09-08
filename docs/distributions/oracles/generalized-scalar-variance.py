"""Independent scalar references for generalized-family range regressions.

Run with Python 3.11+ from any directory. Uses only decimal and the previously
frozen R covariance fixture. It never imports Numerics or regenerates that
fixture. Binary64 inputs are converted exactly before 450-digit arithmetic.
"""

import csv
from decimal import Decimal as D, localcontext
from pathlib import Path


def covariance(rows, family, kappa, hondo):
    """Read the original R covariance without symmetrizing its last-bit output."""
    size = 4 if family == "K4" else 3
    result = [[D(0)] * size for _ in range(size)]
    for row in rows:
        if (row["family"] == family and row["quantity"] == "covariance"
                and float(row["kappa"]) == kappa
                and (row["hondo"] == "NA" or float(row["hondo"]) == hondo)):
            result[int(row["row"]) - 1][int(row["column"]) - 1] = D(row["value"])
    return result


def gno_covariance(k):
    """Exact closed-form GNO information inverse, evaluated directly."""
    v = k * k
    r = ((1 + v) * v.exp() - 1 - 2 * v) / v**2
    c = (1 - (-v / 2).exp()) / v
    return [[1 + c * c / r, -k, c / r],
            [-k, v + D(".5"), k / 2], [c / r, k / 2, v / 2 + 1 / r]]


def transformed_gradient(z, k):
    e = (-k * z).exp()
    return [D(1), (1 - e) / k, (e - 1 + k * z * e) / k**2]


def kappa_gradient(p, k, h):
    logp = p.ln()
    if h == 0:
        t, dt = -logp, -logp * logp / 2
    else:
        e = (h * logp).exp()
        t, dt = (1 - e) / h, (-h * e * logp - 1 + e) / h**2
    logt = t.ln()
    if k == 0:
        return [D(1), -logt, -logt * logt / 2, -dt / t]
    e = (k * logt).exp()
    return [D(1), (1 - e) / k, (e - 1 - k * logt * e) / k**2, -e * dt / t]


def quadratic(matrix, gradient, alpha=1.0, sample_size=100):
    return (sum(gradient[i] * matrix[i][j] * gradient[j]
                for i in range(len(gradient)) for j in range(len(gradient)))
            * D.from_float(alpha)**2 / sample_size)


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 450
        with Path(__file__).with_name("generalized-fisher.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        print("GLO_zero_median", float(covariance(rows, "GLO", 0, -1)[0][0] / 100))
        print("K4_zero_median", float(quadratic(covariance(rows, "K4", 0, 0),
                                               kappa_gradient(D(".5"), D(0), D(0)))))
        print("GNO_covariance_overflow", float(quadratic(gno_covariance(D(2)),
                                                        [D(1), D(0), D(0)], 1e155)))
        p, k = D.from_float(1e-300), D.from_float(.2)
        print("GLO_covariance_underflow", float(quadratic(covariance(rows, "GLO", .2, -1),
                              transformed_gradient(p.ln() - (1 - p).ln(), k), 1e-170)))
        print("K4_covariance_underflow", float(quadratic(covariance(rows, "K4", -1, -.2),
                kappa_gradient(D.from_float(1 - 1e-16), D(-1), D.from_float(-.2)), 1e-170)))
        # Independently frozen R qnorm values for the exact binary64 probabilities.
        for name, k, alpha, z in [
            ("GNO_unit_gradient_overflow", 20.0, 1e-200, -37.047096299361199),
            ("GNO_endpoint_cancellation", 30.0, 1e200, 8.209536151601387),
        ]:
            shape = D.from_float(k)
            print(name, float(quadratic(gno_covariance(shape),
                                  transformed_gradient(D.from_float(z), shape), alpha)))
