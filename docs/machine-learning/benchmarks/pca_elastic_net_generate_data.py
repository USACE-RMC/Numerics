"""Benchmark data for PrincipalComponentAnalysis and ElasticNetRegression.

Usage:  python pca_elastic_net_generate_data.py <output-folder>

Writes little-endian float64 binaries, row-major, read identically by the C#, Python and R timing
scripts. The C# benchmarks find the folder through the NUMERICS_ML_BENCHMARKS environment variable.

PCA:          pca_<n>x<p>.bin                  n x p correlated normal data
Elastic net:  enet_<n>x<p>_X.bin, _y.bin       AR(1)-correlated predictors (rho = 0.5), 20 nonzero
                                               coefficients, noise standard deviation 1
"""
import json
import pathlib
import sys

import numpy as np

PCA_SIZES = [(1_000, 20), (10_000, 50), (100_000, 100)]
ENET_SIZES = [(100, 200), (1_000, 50), (10_000, 100), (100_000, 500)]


def main(folder):
    out = pathlib.Path(folder)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261105)
    for n, p in PCA_SIZES:
        mix = rng.normal(size=(p, p)) / np.sqrt(p)
        X = rng.normal(size=(n, p)) @ mix + rng.normal(size=p) * 10
        X.astype("<f8").tofile(out / f"pca_{n}x{p}.bin")
    for n, p in ENET_SIZES:
        Z = rng.normal(size=(n, p))
        X = np.empty_like(Z)
        X[:, 0] = Z[:, 0]
        for j in range(1, p):
            X[:, j] = 0.5 * X[:, j - 1] + np.sqrt(0.75) * Z[:, j]
        beta = np.zeros(p)
        idx = rng.choice(p, size=20, replace=False)
        beta[idx] = rng.normal(scale=2.0, size=20)
        y = 3.0 + X @ beta + rng.normal(size=n)
        X.astype("<f8").tofile(out / f"enet_{n}x{p}_X.bin")
        y.astype("<f8").tofile(out / f"enet_{n}x{p}_y.bin")
    (out / "sizes.json").write_text(json.dumps({"pca": PCA_SIZES, "enet": ENET_SIZES}))


if __name__ == "__main__":
    main(sys.argv[1])
