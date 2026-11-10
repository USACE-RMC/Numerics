"""Independent numpy oracle for Numerics.Mathematics.Optimization.DynamicProgramming and Knapsack.

The fixtures are defined by integer formulas, so the C# tests rebuild exactly the same inputs.
Floating-point operations are written in the same order as the C# solver (sequential sums over
outcomes), so the optimal values agree to the last bit.

Usage:
    python dynamic_programming_oracle.py golden <out.csv>   # small golden fixtures for the tests
    python dynamic_programming_oracle.py bench              # median wall times of the benchmark sizes
"""
import sys
import time

import numpy as np


# ---------------------------------------------------------------------------
# Reservoir model (shared with the C# tests and benchmark)
#   storage s in 0..S-1, release r in 0..A-1
#   x = s + q - r; infeasible if x < 0; spill = max(x - (S-1), 0); s' = min(x, S-1)
#   cost = d*d + 0.5*spill + 0.1*(S-1-s'), d = r - demand(t)
#   terminal value 0.25*(S-1-s)
# ---------------------------------------------------------------------------

def inflow(t, A):
    return (3 + 5 * t) % A


def demand(t, A):
    return (1 + 3 * t) % A


def stage_cost(s, r, q, D, S):
    """Vectorized stage cost and next storage; returns (cost, next, feasible)."""
    x = s + q - r
    feasible = x >= 0
    spill = np.maximum(x - (S - 1), 0)
    nxt = np.minimum(x, S - 1)
    d = (r - D).astype(np.float64)
    cost = d * d + 0.5 * spill.astype(np.float64) + 0.1 * (S - 1 - nxt).astype(np.float64)
    return cost, np.where(feasible, nxt, 0), feasible


def solve_deterministic(T, S, A, gamma):
    s = np.arange(S)[:, None]
    r = np.arange(A)[None, :]
    V = 0.25 * (S - 1 - np.arange(S)).astype(np.float64)
    values = np.empty((T + 1, S))
    policy = np.empty((T, S), dtype=np.int64)
    values[T] = V
    for t in range(T - 1, -1, -1):
        cost, nxt, feasible = stage_cost(s, r, inflow(t, A), demand(t, A), S)
        cand = cost + gamma * V[nxt]
        cand = np.where(feasible & np.isfinite(cand), cand, np.inf)
        pol = np.argmin(cand, axis=1)
        V = cand[np.arange(S), pol]
        pol = np.where(np.isfinite(V), pol, -1)
        values[t] = V
        policy[t] = pol
    return values, policy


def markov_matrix(K):
    P = np.empty((K, K))
    for j in range(K):
        w = [1.0 / (1.0 + abs(j - jn)) for jn in range(K)]
        total = 0.0
        for x in w:
            total += x
        for jn in range(K):
            P[j, jn] = w[jn] / total
    return P


def solve_stochastic(T, S, A, K, q_step, gamma):
    """State index = storage * K + class; class j has inflow j * q_step; next class ~ P[j, :]."""
    P = markov_matrix(K)
    n = S * K
    storage = (np.arange(n) // K)[:, None]
    cls = (np.arange(n) % K)[:, None]
    r = np.arange(A)[None, :]
    V = np.empty(n)
    for idx in range(n):
        V[idx] = 0.25 * (S - 1 - idx // K)
    values = np.empty((T + 1, n))
    policy = np.empty((T, n), dtype=np.int64)
    values[T] = V
    Pj = P[cls[:, 0]]  # (n, K) transition row of each state's class
    for t in range(T - 1, -1, -1):
        cost, nxt, feasible = stage_cost(storage, r, cls * q_step, demand(t, A), S)
        expected = np.zeros((n, A))
        base = nxt * K
        for jn in range(K):
            p = Pj[:, jn][:, None]
            term = np.where(p > 0, p * V[base + jn], 0.0)
            expected = expected + term
        cand = cost + gamma * expected
        cand = np.where(feasible & np.isfinite(cand), cand, np.inf)
        pol = np.argmin(cand, axis=1)
        V = cand[np.arange(n), pol]
        pol = np.where(np.isfinite(V), pol, -1)
        values[t] = V
        policy[t] = pol
    return values, policy


# ---------------------------------------------------------------------------
# Knapsack (shared fixture): w_i = 1 + (53 i) mod 197, v_i = 1 + (37 i) mod 101 + 0.001 i
# ---------------------------------------------------------------------------

def knapsack_items(n):
    w = np.array([1 + (53 * i) % 197 for i in range(n)], dtype=np.int64)
    v = np.array([1.0 + float((37 * i) % 101) + 0.001 * i for i in range(n)])
    return v, w


def knapsack(v, w, W):
    n = len(v)
    best = np.zeros(W + 1)
    keep = np.zeros((n, W + 1), dtype=bool)
    for i in range(n):
        wi = int(w[i])
        if wi > W or not v[i] > 0:
            continue
        cand = best[: W + 1 - wi] + v[i]
        improve = cand > best[wi:]
        keep[i, wi:] = improve
        best[wi:] = np.where(improve, cand, best[wi:])
    sel = []
    c = W
    for i in range(n - 1, -1, -1):
        if keep[i, c]:
            sel.append(i)
            c -= int(w[i])
    return best[W], sorted(sel)


def golden(path):
    rows = ["case,index,value,policy"]
    vals, pol = solve_deterministic(12, 11, 6, 0.95)
    for s in range(11):
        rows.append(f"det,{s},{float(vals[0, s])!r},{pol[0, s]}")
    for t in range(12):
        rows.append(f"detpolicy,{t},0,{' '.join(map(str, pol[t]))}")
    vals, pol = solve_stochastic(12, 11, 6, 3, 2, 0.95)
    for s in range(33):
        rows.append(f"sto,{s},{float(vals[0, s])!r},{pol[0, s]}")
    for t in range(12):
        rows.append(f"stopolicy,{t},0,{' '.join(map(str, pol[t]))}")
    v, w = knapsack_items(30)
    total, sel = knapsack(v, w, 400)
    rows.append(f"knap,-1,{float(total)!r},{' '.join(map(str, sel))}")
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(rows) + "\n")


def median_time(fn, reps):
    times = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


def bench():
    print(f"numpy {np.__version__}, python {sys.version.split()[0]}")
    det = median_time(lambda: solve_deterministic(365, 200, 50, 0.99), 7)
    v, _ = solve_deterministic(365, 200, 50, 0.99)
    print(f"deterministic T=365 S=200 A=50: {det * 1e3:.2f} ms, V0[100]={float(v[0, 100])!r}")
    sto = median_time(lambda: solve_stochastic(365, 200, 50, 5, 10, 0.99), 5)
    v, _ = solve_stochastic(365, 200, 50, 5, 10, 0.99)
    print(f"stochastic T=365 S=200x5 A=50 K=5: {sto * 1e3:.2f} ms, V0[500]={float(v[0, 500])!r}")
    kv, kw = knapsack_items(200)
    kn = median_time(lambda: knapsack(kv, kw, 10000), 11)
    total, sel = knapsack(kv, kw, 10000)
    print(f"knapsack n=200 W=10000: {kn * 1e3:.2f} ms, value={float(total)!r}, items={len(sel)}")


if __name__ == "__main__":
    if sys.argv[1] == "golden":
        golden(sys.argv[2])
    else:
        bench()
