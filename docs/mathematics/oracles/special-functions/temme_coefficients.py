"""Coefficients of Temme's uniform asymptotic expansion of the incomplete gamma ratios.

Temme (1979) writes, with lambda = x / a and eta = sign(lambda - 1) sqrt(2 (lambda - 1 - ln lambda)),

    Q(a, x) = erfc(eta sqrt(a/2)) / 2 + R,   P(a, x) = erfc(-eta sqrt(a/2)) / 2 - R,
    R ~ exp(-a eta^2 / 2) / sqrt(2 pi a) * sum_k c_k(eta) a^-k

(NIST DLMF 8.12.3, 8.12.4 and 8.12.8), where c_0(eta) = 1/(lambda - 1) - 1/eta and
c_k(eta) = (1/eta) d c_(k-1)/d eta + (-1)^k g_k / (lambda - 1) for k >= 1 (DLMF 8.12.9 to 8.12.11), with g_k the
coefficients of the Stirling series of Gamma*(a) (DLMF 5.11.3). Each c_k is regular at eta = 0, and its Maclaurin
coefficients d[k][n] (c_k(eta) = sum_n d[k][n] eta^n) form the table used by DiDonato and Morris (1986) and by
Cephes igam. This script computes them with mpmath power series at 80 digits, checks the expansion against
mpmath's incomplete gamma, and prints the table in C# form.

Usage:
    python temme_coefficients.py
"""
import sys

import mpmath as mp

mp.mp.dps = 80
K = 20          # number of c_k rows
N = 25          # Maclaurin terms kept per row
M = N + 2 * K + 4  # working order of the power series


def mul(a, b):
    out = [mp.mpf(0)] * M
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(M - i):
            out[i + j] += ai * b[j]
    return out


def inv(a):
    """Reciprocal of a power series with a[0] != 0."""
    out = [mp.mpf(0)] * M
    out[0] = 1 / a[0]
    for n in range(1, M):
        out[n] = -mp.fsum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0]
    return out


def sqrt1(a):
    """Square root of a power series with a[0] = 1."""
    out = [mp.mpf(0)] * M
    out[0] = mp.mpf(1)
    for n in range(1, M):
        out[n] = (a[n] - mp.fsum(out[k] * out[n - k] for k in range(1, n))) / 2
    return out


def compose(a, b):
    """a(b(t)) for b[0] = 0."""
    out = [mp.mpf(0)] * M
    power = [mp.mpf(0)] * M
    power[0] = mp.mpf(1)
    for coefficient in a:
        if coefficient != 0:
            out = [o + coefficient * p for o, p in zip(out, power)]
        power = mul(power, b)
    return out


def eta_over_mu():
    """eta / mu as a series in mu, from eta^2 / 2 = mu - ln(1 + mu)."""
    # 2 (mu - ln(1 + mu)) / mu^2 = sum_{k>=0} 2 (-1)^k mu^k / (k + 2)
    s = [mp.mpf(2) * (-1) ** k / (k + 2) for k in range(M)]
    return sqrt1(s)


def mu_of_eta():
    """mu as a series in eta, by Newton iteration on the composition eta(mu(t)) = t."""
    e = [mp.mpf(0)] + eta_over_mu()[:M - 1]
    de = [(i + 1) * e[i + 1] for i in range(M - 1)] + [mp.mpf(0)]
    mu = [mp.mpf(0)] * M
    mu[1] = mp.mpf(1)
    t = [mp.mpf(0)] * M
    t[1] = mp.mpf(1)
    for _ in range(8):
        residual = [a - b for a, b in zip(compose(e, mu), t)]
        mu = [a - b for a, b in zip(mu, mul(residual, inv(compose(de, mu))))]
    return mu


def gamma_star_coefficients(count):
    """g_k with Gamma*(a) ~ sum_k g_k a^-k, from ln Gamma*(a) = sum_m B_2m / (2m (2m - 1) a^(2m - 1))."""
    log = [mp.mpf(0)] * (count + 1)
    for m in range(1, count // 2 + 2):
        if 2 * m - 1 <= count:
            log[2 * m - 1] = mp.bernoulli(2 * m) / (2 * m * (2 * m - 1))
    out = [mp.mpf(0)] * (count + 1)
    out[0] = mp.mpf(1)
    for n in range(1, count + 1):
        out[n] = mp.fsum(k * log[k] * out[n - k] for k in range(1, n + 1)) / n
    return out


def table():
    mu = mu_of_eta()
    mu_over_eta = mu[1:] + [mp.mpf(0)]
    eta_over_mu_eta = inv(mu_over_eta)  # eta / mu as a series in eta
    g = gamma_star_coefficients(K + 1)
    # c_0 = 1/mu - 1/eta = (eta/mu - 1) / eta
    c = eta_over_mu_eta[1:] + [mp.mpf(0)]
    rows = [c]
    for k in range(1, K):
        derivative = [(n + 1) * c[n + 1] for n in range(M - 1)] + [mp.mpf(0)]
        numerator = [d + (-1) ** k * g[k] * e for d, e in zip(derivative, eta_over_mu_eta)]
        if abs(numerator[0]) > mp.mpf(10) ** -50:
            raise ValueError(f"the pole of c_{k} does not cancel: {numerator[0]}")
        c = numerator[1:] + [mp.mpf(0)]
        rows.append(c)
    return [row[:N] for row in rows]


def check(rows):
    """Worst relative error of the truncated expansion against mpmath for a few shapes and arguments."""
    worst = mp.mpf(0)
    for a in (25, 60, 250, 2000):
        for ratio in (0.75, 0.9, 0.99, 1.0, 1.01, 1.1, 1.25):
            a_ = mp.mpf(a)
            x = a_ * ratio
            lam = x / a_
            eta = mp.sign(lam - 1) * mp.sqrt(2 * (lam - 1 - mp.log(lam)))
            s = mp.fsum(mp.fsum(rows[k][n] * eta ** n for n in range(N)) * a_ ** -k for k in range(K))
            q = mp.erfc(eta * mp.sqrt(a_ / 2)) / 2 + mp.exp(-a_ * eta ** 2 / 2) / mp.sqrt(2 * mp.pi * a_) * s
            exact = mp.gammainc(a_, x, mp.inf, regularized=True)
            worst = max(worst, abs(q / exact - 1))
    return worst


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    rows = table()
    print(f"// worst relative error of the truncated expansion against mpmath: {mp.nstr(check(rows), 3)}")
    print("        private static readonly double[][] TemmeCoefficients =")
    print("        {")
    for k, row in enumerate(rows):
        values = [repr(float(v)).upper() for v in row]
        print("            new double[]")
        print("            {")
        for i in range(0, len(values), 4):
            tail = "," if i + 4 < len(values) else ""
            print("                " + ", ".join(values[i:i + 4]) + tail)
        print("            }" + ("," if k + 1 < len(rows) else ""))
    print("        };")


if __name__ == "__main__":
    main()
