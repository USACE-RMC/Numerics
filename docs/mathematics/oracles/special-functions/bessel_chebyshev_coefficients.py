"""Chebyshev coefficients of the Bessel function approximations in Numerics.Mathematics.SpecialFunctions.Bessel.

Each smooth auxiliary function is sampled by mpmath at 50 significant digits on the Chebyshev
points of its interval and expanded as f(t) = c0/2 + sum_k c_k T_k(t), t in [-1, 1] (Clenshaw, 1962).
The series is truncated where the remaining coefficients sum below 1E-18 of the largest value of f,
and the coefficients are printed order-zero term first, as Evaluate.Chebyshev expects (the Clenshaw
sum of the Cephes routine chbevl, Moshier, 1989). Running this script reproduces the arrays in Bessel.cs.

Usage:
    python bessel_chebyshev_coefficients.py
"""
import sys

import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")

mp.mp.dps = 50
NODES = 96
TOL = mp.mpf("1e-18")


def chebyshev(f, lo, hi):
    """Chebyshev coefficients of f on [lo, hi], truncated by TOL."""
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    theta = [mp.pi * (j + mp.mpf(1) / 2) / NODES for j in range(NODES)]
    xs = [(hi + lo) / 2 + (hi - lo) / 2 * mp.cos(t) for t in theta]
    fs = [f(x) for x in xs]
    scale = max(abs(v) for v in fs)
    c = [2 * mp.fsum(fs[j] * mp.cos(k * theta[j]) for j in range(NODES)) / NODES for k in range(NODES)]
    n = NODES
    while n > 1 and mp.fsum(abs(v) for v in c[n - 1:]) < TOL * scale:
        n -= 1
    if n > NODES - 8:
        raise ValueError("the series did not converge within the node count")
    return c[:n]


def emit(name, comment, coefficients):
    print(f"        /// <summary>{comment}</summary>")
    print(f"        private static readonly double[] {name} =")
    print("        {")
    values = [repr(float(v)) for v in coefficients]
    for i in range(0, len(values), 3):
        line = ", ".join(v.upper().replace("E", "E") for v in values[i:i + 3])
        tail = "," if i + 3 < len(values) else ""
        print(f"            {line}{tail}")
    print("        };")
    print()


def split(z):
    hi = float(z)
    lo = float(z - mp.mpf(hi))
    return hi, lo


def zeros():
    for name, fn in (("J0", lambda k: mp.besseljzero(0, k)), ("J1", lambda k: mp.besseljzero(1, k)),
                     ("Y0", lambda k: mp.besselyzero(0, k)), ("Y1", lambda k: mp.besselyzero(1, k))):
        for k in (1, 2, 3):
            hi, lo = split(fn(k))
            print(f"        // {name} zero {k}: {repr(hi)} + {repr(lo)}")
    print()


def pq(order, x):
    """Hankel's P and Q of the given order, recovered from J and Y."""
    chi = x - (2 * order + 1) * mp.pi / 4
    s = mp.sqrt(mp.pi * x / 2)
    j, y = mp.besselj(order, x), mp.bessely(order, x)
    return s * (j * mp.cos(chi) + y * mp.sin(chi)), s * (y * mp.cos(chi) - j * mp.sin(chi))


J01, J02 = mp.besseljzero(0, 1), mp.besseljzero(0, 2)
J11, J12 = mp.besseljzero(1, 1), mp.besseljzero(1, 2)
Y01, Y02, Y03 = mp.besselyzero(0, 1), mp.besselyzero(0, 2), mp.besselyzero(0, 3)
Y11, Y12 = mp.besselyzero(1, 1), mp.besselyzero(1, 2)
ASYMPTOTIC = 25


def main():
    zeros()
    sq = mp.sqrt
    emit("J0Small", "J₀(x) / ((x² - j₀,₁²)(x² - j₀,₂²)) in z = x² on [0, 64].",
         chebyshev(lambda z: mp.besselj(0, sq(z)) / ((z - J01**2) * (z - J02**2)), 0, 64))
    emit("J1Small", "J₁(x) / (x (x² - j₁,₁²)(x² - j₁,₂²)) in z = x² on [0, 64].",
         chebyshev(lambda z: mp.besselj(1, sq(z)) / (sq(z) * (z - J11**2) * (z - J12**2)), 0, 64))
    emit("Y0Log", "Y₀(x) - (2/π) ln(x) J₀(x) in z = x² on [0, 0.25].",
         chebyshev(lambda z: mp.bessely(0, sq(z)) - 2 / mp.pi * mp.log(sq(z)) * mp.besselj(0, sq(z)), 0, 0.25))
    emit("Y0Zero1", "Y₀(x) / (x - y₀,₁) on [0.5, 1.5].", chebyshev(lambda x: mp.bessely(0, x) / (x - Y01), 0.5, 1.5))
    emit("Y0Mid", "Y₀(x) on [1.5, 3].", chebyshev(lambda x: mp.bessely(0, x), 1.5, 3))
    emit("Y0Zero2", "Y₀(x) / (x - y₀,₂) on [3, 5.5].", chebyshev(lambda x: mp.bessely(0, x) / (x - Y02), 3, 5.5))
    emit("Y0Zero3", "Y₀(x) / (x - y₀,₃) on [5.5, 8].", chebyshev(lambda x: mp.bessely(0, x) / (x - Y03), 5.5, 8))
    emit("Y1Log", "(Y₁(x) - (2/π)(ln(x) J₁(x) - 1/x)) / x in z = x² on [0, 2.25].",
         chebyshev(lambda z: (mp.bessely(1, sq(z)) - 2 / mp.pi * (mp.log(sq(z)) * mp.besselj(1, sq(z)) - 1 / sq(z))) / sq(z), 0, 2.25))
    emit("Y1Zero1", "Y₁(x) / (x - y₁,₁) on [1.5, 3.5].", chebyshev(lambda x: mp.bessely(1, x) / (x - Y11), 1.5, 3.5))
    emit("Y1Zero2", "Y₁(x) / (x - y₁,₂) on [3.5, 8].", chebyshev(lambda x: mp.bessely(1, x) / (x - Y12), 3.5, 8))
    for order in (0, 1):
        emit(f"P{order}Mid", f"Hankel's P{order}(x) in u = 8/x on [8, {ASYMPTOTIC}].",
             chebyshev(lambda u: pq(order, 8 / u)[0], mp.mpf(8) / ASYMPTOTIC, 1))
        emit(f"Q{order}Mid", f"Hankel's Q{order}(x) / u in u = 8/x on [8, {ASYMPTOTIC}].",
             chebyshev(lambda u: pq(order, 8 / u)[1] / u, mp.mpf(8) / ASYMPTOTIC, 1))
    emit("K0Small", "K₀(x) + ln(x/2) I₀(x) in z = x² on [0, 4].",
         chebyshev(lambda z: mp.besselk(0, sq(z)) + mp.log(sq(z) / 2) * mp.besseli(0, sq(z)), 0, 4))
    emit("K0Large", "exp(x) sqrt(x) K₀(x) in t = 4/x - 1 on [2, ∞).",
         chebyshev(lambda t: (lambda x: mp.exp(x) * sq(x) * mp.besselk(0, x))(4 / (t + 1)), -1, 1))
    emit("K1Small", "x (K₁(x) - ln(x/2) I₁(x)) in z = x² on [0, 4].",
         chebyshev(lambda z: sq(z) * (mp.besselk(1, sq(z)) - mp.log(sq(z) / 2) * mp.besseli(1, sq(z))), 0, 4))
    emit("K1Large", "exp(x) sqrt(x) K₁(x) in t = 4/x - 1 on [2, ∞).",
         chebyshev(lambda t: (lambda x: mp.exp(x) * sq(x) * mp.besselk(1, x))(4 / (t + 1)), -1, 1))


if __name__ == "__main__":
    main()
