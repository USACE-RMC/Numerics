"""mpmath oracle for the special functions of Numerics.Mathematics.SpecialFunctions.

Every reference value is evaluated by mpmath at 50 significant digits and rounded once to the
nearest double. The arguments are doubles, written with 17 significant digits so that the C#
tests parse exactly the same inputs.

Columns: function, n (integer order, or 0), a (shape parameter, or 0), x (argument or
probability), value (the reference result).

Functions and their mpmath definitions:
    i0, i1, in        besseli(n, x)
    k0, k1, kn        besselk(n, x)
    j0, j1, jn        besselj(n, x)
    y0, y1, yn        bessely(n, x)
    debye3            3 / x^3 * integral_0^x t^3 / (e^t - 1) dt
    gammap            gammainc(a, 0, x, regularized=True)     (lower regularized P)
    gammaq            gammainc(a, x, inf, regularized=True)   (upper regularized Q)
    gammapinv         x with P(a, x) = p, p the double in column x
    gammaqinv         x with Q(a, x) = q, q the double in column x
    erfc              erfc(x)
    ndtr              ncdf(x), the standard normal CDF

Usage:
    python special_functions_oracle.py <out.csv>
"""
import sys

import mpmath as mp

mp.mp.dps = 50


def d(v):
    """Round an mpf to the nearest double."""
    return float(v)


def fmt(v):
    return repr(float(v))


def logspace(lo, hi, count):
    """Doubles spaced evenly in log10 from lo to hi (the arguments are rounded before evaluation)."""
    return [float(10.0 ** (float(lo) + (float(hi) - float(lo)) * i / (count - 1))) for i in range(count)]


def linspace(lo, hi, count):
    return [lo + (hi - lo) * i / (count - 1) for i in range(count)]


def debye3(x):
    x = mp.mpf(x)
    return 3 / x**3 * mp.quad(lambda t: t**3 / mp.expm1(t), [0, x])


def rows():
    out = []

    # Modified Bessel functions of the first kind.
    xs = logspace(-8, mp.log10(700), 120) + linspace(0.25, 30.0, 120)
    xs += [-x for x in xs[::7]]
    for x in xs:
        out.append(("i0", 0, 0.0, x, mp.besseli(0, x)))
        out.append(("i1", 0, 0.0, x, mp.besseli(1, x)))
    for n in (2, 3, 5, 10, 20, 50):
        for x in logspace(-3, mp.log10(600), 60) + [-0.5, -7.0, -40.0]:
            out.append(("in", n, 0.0, x, mp.besseli(n, x)))

    # Modified Bessel functions of the second kind.
    xs = logspace(-8, mp.log10(700), 120) + linspace(0.25, 30.0, 120)
    for x in xs:
        out.append(("k0", 0, 0.0, x, mp.besselk(0, x)))
        out.append(("k1", 0, 0.0, x, mp.besselk(1, x)))
    for n in (2, 3, 5, 10, 20, 50):
        for x in logspace(-2, mp.log10(600), 60):
            out.append(("kn", n, 0.0, x, mp.besselk(n, x)))

    # Bessel functions of the first and second kind: a uniform grid, a wide logarithmic grid and
    # points within 1E-3 and 1E-6 of the first twenty zeros.
    base = linspace(0.05, 60.0, 400) + logspace(-8, 4, 120)
    for name, order in (("j0", 0), ("j1", 1)):
        xs = list(base)
        for k in range(1, 21):
            z = d(mp.besseljzero(order, k))
            xs += [z, z - 1e-3, z + 1e-3, z - 1e-6, z + 1e-6]
        xs += [-x for x in base[::9]]
        for x in xs:
            out.append((name, 0, 0.0, x, mp.besselj(order, x)))
    for name, order in (("y0", 0), ("y1", 1)):
        xs = list(base)
        for k in range(1, 21):
            z = d(mp.besselyzero(order, k))
            xs += [z, z - 1e-3, z + 1e-3, z - 1e-6, z + 1e-6]
        for x in xs:
            out.append((name, 0, 0.0, x, mp.bessely(order, x)))
    for n in (2, 3, 5, 10, 20, 50):
        for x in linspace(0.1, 80.0, 120) + logspace(-2, 3, 30):
            out.append(("jn", n, 0.0, x, mp.besselj(n, x)))
            out.append(("yn", n, 0.0, x, mp.bessely(n, x)))
        out.append(("jn", n, 0.0, -3.5, mp.besselj(n, -3.5)))

    # Debye function of order 3.
    for x in logspace(-6, 3, 150) + linspace(0.05, 12.0, 120):
        out.append(("debye3", 0, 0.0, x, debye3(x)))

    # Regularized incomplete gamma functions over shapes from 1E-3 to 1E6 and arguments from
    # 1E-4 a to 30 a, plus the far tails.
    shapes = logspace(-3, 6, 28)
    ratios = [1e-4, 1e-2, 0.1, 0.3, 0.5, 0.8, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.25, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0]
    for a in shapes:
        sd = float(mp.sqrt(a))
        xs = [a * r for r in ratios] + [max(a - k * sd, a * 1e-3) for k in (1, 3, 6)] + [a + k * sd for k in (1, 3, 6)]
        for x in xs:
            if x <= 0:
                continue
            p, q = gamma_pq(a, x)
            out.append(("gammap", 0, a, x, p))
            out.append(("gammaq", 0, a, x, q))
    for a, x in ((0.5, 800.0), (10.0, 400.0), (2.0, 1e-160), (100.0, 0.5), (1000.0, 200.0), (1.5, 700.0)):
        p, q = gamma_pq(a, x)
        out.append(("gammap", 0, a, x, p))
        out.append(("gammaq", 0, a, x, q))

    # Complementary error function and standard normal CDF over the whole range where the results are
    # normal doubles, plus points in the subnormal tail.
    for x in linspace(-6.0, 26.5, 651) + logspace(-300, 0, 40) + [26.6, 27.0, 27.2]:
        out.append(("erfc", 0, 0.0, x, mp.erfc(x)))
    for x in linspace(-37.5, 8.5, 921) + [-38.0, -38.4, -1e-300]:
        out.append(("ndtr", 0, 0.0, x, mp.ncdf(x)))

    # Inverse incomplete gamma: the reference root is solved in mpmath for the double probability.
    probs = [1e-300, 1e-100, 1e-20, 1e-10, 1e-5, 0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999, 1 - 1e-10]
    for a in logspace(-2, 5, 15):
        for p in probs:
            for name, lower in (("gammapinv", True), ("gammaqinv", False)):
                root = invert(a, p, lower)
                if root is None:
                    continue
                out.append((name, 0, a, p, root))
    return out


def gamma_pq(a, x):
    """Return (P, Q); the smaller one is evaluated directly and the other as its complement."""
    a = mp.mpf(a)
    x = mp.mpf(x)
    if x < a:
        # P(a, x) = x^a e^-x / Gamma(a + 1) * 1F1(1; a + 1; x), summed with a raised term limit.
        p = mp.exp(a * mp.log(x) - x - mp.loggamma(a + 1)) * mp.hyp1f1(1, a + 1, x, maxterms=10**7)
        return p, 1 - p
    q = upper_cf(a, x)
    return 1 - q, q


def upper_cf(a, x):
    """Q(a, x) for x >= a by the Legendre continued fraction, summed by the modified Lentz method."""
    tiny = mp.mpf(10) ** -300
    b = x + 1 - a
    c = 1 / tiny
    dd = 1 / b
    h = dd
    i = 1
    while True:
        an = -i * (i - a)
        b += 2
        dd = an * dd + b
        if abs(dd) < tiny:
            dd = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        dd = 1 / dd
        delta = dd * c
        h *= delta
        if abs(delta - 1) < mp.mpf(10) ** -45:
            break
        i += 1
    return mp.exp(a * mp.log(x) - x - mp.loggamma(a)) * h


def invert(a, p, lower):
    """Solve P(a, x) = p (lower) or Q(a, x) = p in log space by bisection on log x."""
    a = mp.mpf(a)
    p = mp.mpf(p)

    def f(logx):
        x = mp.exp(logx)
        lower_p, upper_q = gamma_pq(a, x)
        if lower:
            return mp.log(lower_p) - mp.log(p)
        return mp.log(p) - mp.log(upper_q)

    lo, hi = mp.mpf(-800), mp.log(a) + 60 + mp.log(1 + a)
    flo, fhi = f(lo), f(hi)
    if flo > 0 or fhi < 0:
        return None
    for _ in range(400):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < mp.mpf(10) ** -35:
            break
    x = mp.exp((lo + hi) / 2)
    if x < mp.mpf("2.3e-308") or x > mp.mpf("1.7e308"):
        return None
    return x


def main(path):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write("function,n,a,x,value\n")
        for name, n, a, x, value in rows():
            fh.write(f"{name},{n},{fmt(a)},{fmt(x)},{fmt(value)}\n")


if __name__ == "__main__":
    main(sys.argv[1])
