"""mpmath 50-digit oracle for the AdaptiveGaussKronrod comparison cases (mpmath 1.3.0).

Usage: python agk_oracle.py compare.csv oracle.csv, where compare.csv comes from
dotnet run -c Release --project Benchmarks_Numerics -- agk-compare compare.csv.
"""
import sys
import mpmath as mp

mp.mp.dps = 50
F = mp.mpf
pi = mp.pi


def phi(z):
    return mp.npdf(z)


def Phi(z):
    return mp.ncdf(z)


def ndtri(p):
    return -mp.sqrt(2) * mp.erfinv(1 - 2 * F(p)) if F(p) != F(1) / 2 else F(0)


GXI, GAL, GKA = F(100), F(30), F(-0.1)


def gev_q(p):
    return GXI + GAL / GKA * (1 - (-mp.log(p)) ** GKA)


def gev_pdf(x):
    y = -mp.log(1 - GKA * (x - GXI) / GAL) / GKA
    return mp.exp(-(1 - GKA) * y - mp.exp(-y)) / GAL


def gamma_pdf(x):  # scale 10, shape 5
    return x ** 4 * mp.exp(-x / 10) / (mp.gamma(5) * F(10) ** 5)


LN_MU, LN_SIGMA = F("2.282974736417405"), F("0.198042200435365")


def ln_pdf(x):
    return mp.exp(-(mp.log(x) - LN_MU) ** 2 / (2 * LN_SIGMA ** 2)) / (x * LN_SIGMA * mp.sqrt(2 * pi))


def frag(x):
    return Phi(mp.log(x / 250) / F(0.25))


def damage(x):
    return F(10) ** 6 * max(F(0), x - 150) / (x + 50)


def peak(x):
    t = (x - F(0.3)) * 60
    return mp.exp(-t * t) + F(0.1) * x


B = mp.beta(F(2.2), F(3.8))
THR = [ndtri(F("0.01")), ndtri(F("0.02")), ndtri(F("0.05"))]


def union_z(z):
    s = F(0)
    for t in THR:
        s += mp.log1p(-Phi((t - mp.sqrt(F(0.5)) * z) / mp.sqrt(F(0.5))))
    return -mp.expm1(s) * phi(z)


P150 = mp.exp(-((1 - (F(150) - GXI) * GKA / GAL) ** (1 / GKA)))  # Q(P150) = 150

FUN = {
    "t_x3": lambda x: x ** 3,
    "t_cos": mp.cos,
    "t_sin": mp.sin,
    "t_fxx": lambda x: F(0.5) + 24 * x + 3 * x * x,
    "t_fxxx": lambda x: F(0.5) + 24 * x + 3 * x * x + 8 * x ** 3,
    "t_gamma_mean": lambda x: x * gamma_pdf(x),
    "t_cvar": lambda x: x * ln_pdf(x),
    "t_expcos": lambda x: mp.exp(-x) * mp.cos(3 * x),
    "t_sinc": lambda x: mp.sinc(x - F(0.5)),
    "t_peak_rec": peak, "t_peak_bins": peak, "t_peak": peak,
    "t_x2_d2": lambda x: x * x,
    "n_moment_gev_mean": lambda x: x * gev_pdf(x),
    "n_moment_gev_m2": lambda x: (x - 100) ** 2 * gev_pdf(x),
    "n_moment_gamma_m4": lambda x: ((x - 50) / mp.sqrt(500)) ** 4 * gamma_pdf(x),
    "n_nct": lambda x: x * x / 9 * mp.exp(-(x - 3) ** 2 / 2),
    "n_pert_mean": lambda x: x * x ** F(1.2) * (10 - x) ** F(2.8) / (B * F(10) ** 5),
    "n_pertz_mean": lambda z: Phi(z) * (z + 2) ** F(1.2) * (3 - z) ** F(2.8) / (B * F(5) ** 5),
    "n_truncated_m2": lambda t: t * t * phi(F(0.5) + F(1.5) * t) / (Phi(2) - Phi(-1)) * F(1.5),
    "n_vonmises_mid": lambda t: mp.exp(10 - 20 * mp.sin((t - F(0.5)) / 2) ** 2),
    "n_vonmises_tail": lambda t: mp.exp(10 - 20 * mp.sin((t - F(0.5)) / 2) ** 2),
    "tr_probe_pf": lambda p: frag(gev_q(p)),
    "tr_risk_mean": lambda p: frag(gev_q(p)) * damage(gev_q(p)),
    "tr_risk_ensemble": lambda p: frag(gev_q(p)) * damage(gev_q(p)),
    "tr_conditional": lambda z: phi(z) * Phi((z - F(1.5)) / F(0.3)),
    "h_lorentz_peak": lambda x: 1 / ((x - F(0.3)) ** 2 + F(10) ** -4),
    "h_cos100": lambda x: mp.cos(100 * x),
    "h_xsin50": lambda x: x * mp.sin(50 * x),
}

CLOSED = {
    "t_step_d3": F(7) / 10,
    "t_inv_sqrt_abs": F(4),
    "n_dmi_normal_k0": F(1) / 2, "n_dmi_normal_k2": F(1) / 2, "n_dmi_normal_k4": F(3) / 2,
    "n_dmi_t10_k4": F(25) / 8,
    "n_mvt_radial": 1 - mp.betainc(F(5) / 2, F(1) / 2, 0, 5 / (5 + F(1.3) ** 2), regularized=True) / 2,
    "n_kappa4_mean": mp.euler,
    "n_kappa4_m2": mp.euler ** 2 + pi ** 2 / 6,
    "h_inv_sqrt": F(2), "h_inv_sqrt_tight": F(2), "h_log": F(-1), "h_x_m09": F(10), "h_log_over_sqrt": F(-4),
    "h_abs_kink": F(5) / 18, "h_exp": mp.e - 1, "h_runge": pi / 4,
    "h_gauss_wide": mp.sqrt(pi) * mp.erf(10),
}

SPLITS = {"t_peak_rec": [F(0.3)], "t_peak_bins": [F(0.3)], "t_peak": [F(0.3)], "h_lorentz_peak": [F(0.3)],
          "tr_risk_mean": [P150], "tr_risk_ensemble": [P150]}


def integrate(name, lo, hi):
    f = FUN[name]
    pts = [lo] + [s for s in SPLITS.get(name, []) if lo < s < hi] + [hi]
    if name in ("h_cos100", "h_xsin50"):
        n = 400
        pts = [lo + (hi - lo) * k / n for k in range(n + 1)]
    return mp.quad(f, pts, maxdegree=10)


def main(path, out):
    rows = []
    with open(path) as fh:
        header = fh.readline().strip().split(";")
        for line in fh:
            rec = dict(zip(header, line.strip().split(";")))
            name = rec["name"]
            if name == "t_inv_x":
                rows.append((name, "divergent"))
                continue
            if name == "n_union_single_factor":
                lo, hi = (F(float(v)) for v in rec["ranges"].split(":"))
                val = mp.quad(union_z, [ndtri(lo), -3, 0, 3, ndtri(hi)], maxdegree=10)
            elif name in CLOSED:
                val = CLOSED[name]
            else:
                val = F(0)
                for rng in rec["ranges"].split(","):
                    lo, hi = (F(float(v)) for v in rng.split(":"))
                    if hi > lo:
                        val += integrate(name, lo, hi)
            rows.append((name, mp.nstr(val, 40)))
    with open(out, "w") as fh:
        for r in rows:
            fh.write(";".join(r) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
