"""Coefficient tables of the scaled complementary error function in Erf.cs, recomputed with mpmath.

erfcx(x) = exp(x^2) erfc(x) is evaluated in two pieces:

* 0 <= x < 3: the Chebyshev series of Schonfelder (1978) in t = (8x - 30) / (4x + 15), summed by the
  Clenshaw recurrence b_k = t b_(k+1) - b_(k+2) + a_k with erfcx = (b_0 - b_2) / 2, that is
  erfcx = a_0 / 2 + sum_(k >= 1) a_k T_k(t / 2). The coefficients a_k are recomputed here by
  Chebyshev-Gauss quadrature and compared with the table (ScaledErfcCoefficients).
* x >= 3: g(s) = sqrt(pi) x erfcx(x) as a polynomial in s = 1 / x^2 on [0, 1/9], the near-minimax
  Chebyshev fit of mpmath.chebyfit (ScaledErfcAsymptoticCoefficients, highest degree first).

Both pieces are then evaluated in double precision, as Erf.cs does, and compared with mpmath.

Usage:
    python erfc_coefficients.py
"""
import mpmath as mp

mp.mp.dps = 50

# Schonfelder (1978) coefficients as tabulated in Erf.ScaledErfcCoefficients, kept as decimal strings.
TABLE = [
    "6.10143081923200417926465815756E-1", "-4.34841272712577471828182820888E-1",
    "1.76351193643605501125840298123E-1", "-6.0710795609249414860051215825E-2",
    "1.7712068995694114486147141191E-2", "-4.321119385567293818599864968E-3",
    "8.54216676887098678819832055E-4", "-1.27155090609162742628893940E-4",
    "1.1248167243671189468847072E-5", "3.13063885421820972630152E-7",
    "-2.70988068537762022009086E-7", "3.0737622701407688440959E-8",
    "2.515620384817622937314E-9", "-1.028929921320319127590E-9",
    "2.9944052119949939363E-11", "2.6051789687266936290E-11",
    "-2.634839924171969386E-12", "-6.43404509890636443E-13",
    "1.12457401801663447E-13", "1.7281533389986098E-14",
    "-4.264101694942375E-15", "-5.45371977880191E-16",
    "1.58697607761671E-16", "2.0899837844334E-17",
    "-5.900526869409E-18", "-9.41893387554E-19",
    "2.14977356470E-19", "4.6660985008E-20",
    "-7.243011862E-21", "-2.387966824E-21",
    "1.91177535E-22", "1.20482568E-22",
    "-6.72377E-25", "-5.747997E-24",
    "-4.28493E-25", "2.44856E-25",
    "4.3793E-26", "-8.151E-27",
    "-3.089E-27", "9.3E-29",
    "1.74E-28", "1.6E-29",
    "-8.0E-30", "-2.0E-30",
]
SCHONFELDER_TERMS = 27
SPLIT = 3.0
ASYMPTOTIC_TERMS = 22


def erfcx(x):
    x = mp.mpf(x)
    return mp.erfc(x) * mp.exp(x * x)


def schonfelder_coefficients(count, nodes=400):
    """a_k = (2 / M) sum_j erfcx(x(u_j)) T_k(u_j) at the M Chebyshev-Gauss nodes u_j, with t = 2u."""
    out = [mp.mpf(0)] * count
    for j in range(nodes):
        theta = mp.pi * (j + mp.mpf(1) / 2) / nodes
        t = 2 * mp.cos(theta)
        x = 15 * (t + 2) / (4 * (2 - t))
        f = erfcx(x)
        for k in range(count):
            out[k] += f * mp.cos(k * theta)
    return [2 * v / nodes for v in out]


def g(s):
    if s == 0:
        return mp.mpf(1)
    x = 1 / mp.sqrt(s)
    return mp.sqrt(mp.pi) * x * erfcx(x)


def clenshaw_double(x, coefficients):
    t = (8.0 * x - 30.0) / (4.0 * x + 15.0)
    bm = b = bp = 0.0
    for a in reversed(coefficients):
        bp, b = b, bm
        bm = t * b - bp + a
    return 0.5 * (bm - bp)


def horner_double(s, coefficients):
    v = 0.0
    for a in coefficients:
        v = v * s + a
    return v


def main():
    recomputed = schonfelder_coefficients(len(TABLE))
    worst = max(abs(r - mp.mpf(a)) for r, a in zip(recomputed, TABLE))
    print(f"Schonfelder table: max |recomputed - tabulated| = {mp.nstr(worst, 3)} over {len(TABLE)} coefficients")

    asymptotic, fit_error = mp.chebyfit(g, [0, 1 / mp.mpf(SPLIT) ** 2], ASYMPTOTIC_TERMS, error=True)
    print(f"Asymptotic fit on s in [0, 1/{SPLIT:g}^2], {ASYMPTOTIC_TERMS} terms: max error {mp.nstr(fit_error, 3)}")
    print("ScaledErfcAsymptoticCoefficients (highest degree first):")
    doubles = [float(c) for c in asymptotic]
    for c in doubles:
        print(f"    {c!r},")

    table = [float(a) for a in TABLE[:SCHONFELDER_TERMS]]
    w1 = max(abs(clenshaw_double(x, table) / erfcx(x) - 1) for x in (k / 400 for k in range(int(SPLIT * 400))))
    xs = [SPLIT + k * (40.0 - SPLIT) / 4000 for k in range(4001)]
    w2 = max(abs(horner_double(1.0 / (x * x), doubles) / g(1 / mp.mpf(x) ** 2) - 1) for x in xs)
    print(f"Double precision: Clenshaw (x < {SPLIT:g}) max rel {mp.nstr(w1, 3)}; "
          f"Horner g (x >= {SPLIT:g}) max rel {mp.nstr(w2, 3)}")


if __name__ == "__main__":
    main()
