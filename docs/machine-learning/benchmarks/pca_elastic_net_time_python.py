"""Times scikit-learn on the benchmark data written by pca_elastic_net_generate_data.py.

Usage:  python pca_elastic_net_time_python.py <data-folder> [repeats]

Writes pca_elastic_net_timings_python.json (median seconds per case) with library versions and the BLAS vendor.
Cases, matched to the C# benchmarks:
  pca_fit     PCA().fit(X)  (svd_solver="auto", which chooses covariance_eigh for n >= 10 p)
  enet_path   standardize X and y, then enet_path with l1_ratio = 0.5, 100 penalties,
              eps = 1e-4 (n > p) or 1e-2 (n <= p), tol = 1e-7, precompute = "auto"
  enet_cv     ElasticNetCV with the same path and 10-fold KFold; the faster of n_jobs = 1 and n_jobs = -1
"""
import json
import pathlib
import platform
import statistics
import sys
import time

import numpy as np
import scipy
import sklearn
from threadpoolctl import threadpool_info
from sklearn.decomposition import PCA
from sklearn.linear_model import ElasticNetCV, enet_path
from sklearn.model_selection import KFold


def median_time(fn, repeats):
    fn()  # warm-up
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return statistics.median(times)


def main(folder, repeats):
    root = pathlib.Path(folder)
    sizes = json.loads((root / "sizes.json").read_text())
    result = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "blas": [{k: d.get(k) for k in ("internal_api", "version", "num_threads", "architecture")} for d in threadpool_info()],
        "cpu": platform.processor(),
        "cases": {},
    }
    for n, p in sizes["pca"]:
        X = np.fromfile(root / f"pca_{n}x{p}.bin", dtype="<f8").reshape(n, p)
        result["cases"][f"pca_fit_{n}x{p}"] = median_time(lambda: PCA().fit(X), repeats)
        print("pca", n, p, result["cases"][f"pca_fit_{n}x{p}"], flush=True)
    for n, p in sizes["enet"]:
        X = np.fromfile(root / f"enet_{n}x{p}_X.bin", dtype="<f8").reshape(n, p)
        y = np.fromfile(root / f"enet_{n}x{p}_y.bin", dtype="<f8")
        eps = 1e-4 if n > p else 1e-2
        reps = repeats if n * p < 10_000_000 else max(3, repeats // 2)

        def path():
            Xs = (X - X.mean(axis=0)) / X.std(axis=0)
            ys = (y - y.mean()) / y.std()
            enet_path(Xs, ys, l1_ratio=0.5, eps=eps, n_alphas=100, tol=1e-7, max_iter=100_000)

        result["cases"][f"enet_path_{n}x{p}"] = median_time(path, reps)

        best = None
        for jobs in (None, -1):
            cv = ElasticNetCV(l1_ratio=0.5, eps=eps, alphas=100, tol=1e-7, max_iter=100_000,
                              cv=KFold(10), n_jobs=jobs)
            t = median_time(lambda: cv.fit(X, y), reps)
            best = t if best is None else min(best, t)
        result["cases"][f"enet_cv_{n}x{p}"] = best
        print("enet", n, p, result["cases"][f"enet_path_{n}x{p}"], best, flush=True)
    (root / "pca_elastic_net_timings_python.json").write_text(json.dumps(result, indent=1))
    print(json.dumps(result["cases"], indent=1))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 7)
