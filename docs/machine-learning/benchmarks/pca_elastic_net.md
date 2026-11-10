# Machine Learning Benchmarks

Timing comparison of `PrincipalComponentAnalysis` and `ElasticNetRegression` with scikit-learn and R on identical data.

## Running

```bash
python docs/machine-learning/benchmarks/pca_elastic_net_generate_data.py <folder>
python docs/machine-learning/benchmarks/pca_elastic_net_time_python.py <folder> 5
Rscript docs/machine-learning/benchmarks/pca_elastic_net_time_r.R <folder> 5
```

Then set `NUMERICS_ML_BENCHMARKS=<folder>` and run the comparison, which times the C# methods on the same files and fails unless each median is below both the Python and R medians. Run it on an otherwise idle machine:

```bash
dotnet run -c Release --project Benchmarks_Numerics -- compare pca-enet
```

`compare gpr` times `GaussianProcessRegression` fits and predictions against the scikit-learn and DiceKriging reference timings in `gpr-reference-timings.json`, recorded by `gpr_benchmark_sklearn.py` and `gpr_benchmark_dicekriging.R`; `compare` with no name runs both.

The `Benchmarks_Numerics` project has BenchmarkDotNet versions of the same cases (`PcaBenchmarks`, `ElasticNetBenchmarks`), which read the folder from the same environment variable:

```bash
dotnet run -c Release --project Benchmarks_Numerics -- --filter *MachineLearning*
```

## Cases

| Case | C# | scikit-learn | R |
|---|---|---|---|
| `pca_fit` | `new PrincipalComponentAnalysis(X).Train()` | `PCA().fit(X)` (auto solver; `covariance_eigh` for n ≥ 10p) | `prcomp(X, retx = FALSE)` |
| `enet_path` | `Fit()` with `Alpha = 0.5` | standardize, then `enet_path(l1_ratio=0.5, 100 penalties, tol=1e-7)` | `glmnet(X, y, alpha = 0.5)` |
| `enet_cv` | `CrossValidate(10, seed)` with `Alpha = 0.5` | `ElasticNetCV(l1_ratio=0.5, 100 penalties, KFold(10), tol=1e-7)`, faster of `n_jobs` = 1 and -1 | `cv.glmnet(X, y, alpha = 0.5, nfolds = 10)` |

All three use 100 penalties with the glmnet default ratio (10⁻⁴ for n > p, 10⁻² otherwise). C# and glmnet use the same convergence rule and threshold (10⁻⁷); scikit-learn's duality-gap tolerance is set to the same value. Every timing is the median of repeated runs after a warm-up.

## Results

Intel Core Ultra 9 185H (16 cores), Windows 11, .NET 10; Python 3.14.3, NumPy 2.4.2, SciPy 1.17.1, scikit-learn 1.8.0 with OpenBLAS 0.3.31 (scipy-openblas, multithreaded); R 4.4.3 with glmnet 4.1.10 and R's reference BLAS and LAPACK. The machine was shared with other work while these were recorded, so absolute times vary by up to a factor of two between runs; the ratios are stable.

| Case | n × p | C# | scikit-learn | R |
|---|---|---|---|---|
| PCA fit | 1,000 × 20 | 0.24 ms | 0.35 ms | 1.5 ms |
| PCA fit | 10,000 × 50 | 1.1 ms | 2.7 ms | 78 ms |
| PCA fit | 100,000 × 100 | 16 ms | 48 ms | 4.4 s |
| Elastic net path | 100 × 200 | 2.5 ms | 10.9 ms | 5.2 ms |
| Elastic net path | 1,000 × 50 | 1.3 ms | 10.5 ms | 4.1 ms |
| Elastic net path | 10,000 × 100 | 15 ms | 33 ms | 45 ms |
| Elastic net path | 100,000 × 500 | 0.82 s | 16.6 s | 6.7 s |
| Elastic net 10-fold CV | 100 × 200 | 27 ms | 125 ms | 110 ms |
| Elastic net 10-fold CV | 1,000 × 50 | 8.3 ms | 102 ms | 165 ms |
| Elastic net 10-fold CV | 10,000 × 100 | 46 ms | 271 ms | 890 ms |
| Elastic net 10-fold CV | 100,000 × 500 | 1.6 s | 203 s | 48 s |

glmnet uses its covariance updates for p < 500 and its naive updates for p = 500, as by default.

