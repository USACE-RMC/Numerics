"""Joins the agk-compare CSV with the agk_oracle.py output into the old-versus-new error table.

Usage: python agk_table.py <dir holding compare.csv and oracle.csv>
"""
import sys, mpmath as mp
mp.mp.dps = 50
s = sys.argv[1]
orc = dict(l.strip().split(";") for l in open(s + r"\oracle.csv"))
rows = [l.strip().split(";") for l in open(s + r"\compare.csv")]
h = rows[0]
print("name|tol|err_old|err_new|ev_old|ev_new|st_old|st_new|flag")
for r in rows[1:]:
    d = dict(zip(h, r)); n = d["name"]; o = orc[n]
    if o == "divergent":
        print(f"{n}|-|div|div|{d['old_evals']}|{d['new_evals']}|{d['old_status']}|{d['new_status']}|"); continue
    O = mp.mpf(o)
    def err(v):
        v = float(v)
        return float(abs(mp.mpf(v) - O)) if v == v else float("nan")
    eo, en = err(d["old_result"]), err(d["new_result"])
    tol = max(float(d["abs"]), float(d["rel"]) * abs(float(O)))
    flag = []
    if not (en <= eo) and en > 0: flag.append("LESS_ACC" + ("(within tol)" if en <= tol else "(OVER TOL)"))
    if int(d["new_evals"]) > int(d["old_evals"]): flag.append("MORE_EVALS")
    comp = ""
    if d.get("old_composite"):
        comp = f" composite_err old {err(d['old_composite']):.2e} new {err(d['new_composite']):.2e}"
    print(f"{n}|{tol:.1e}|{eo:.2e}|{en:.2e}|{d['old_evals']}|{d['new_evals']}|{d['old_status'][:12]}|{d['new_status'][:12]}|{' '.join(flag)}{comp}")
