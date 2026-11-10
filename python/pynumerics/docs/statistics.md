# Statistics, hypothesis tests, goodness of fit and probability

pynumerics wraps `Numerics.Data.Statistics` in four modules:

| Module | Numerics types |
| --- | --- |
| `pynumerics.statistics` | `Statistics`, `Correlation`, `Autocorrelation`, `BoxCox`, `YeoJohnson`, `PlottingPositions`, `RunningStatistics`, `RunningCovarianceMatrix`, `Histogram`, `Histogram.Bin`, `GlobalSensitivity` |
| `pynumerics.hypothesis_tests` | `HypothesisTests`, `MultipleGrubbsBeckTest` |
| `pynumerics.goodness_of_fit` | `GoodnessOfFit` |
| `pynumerics.probability` | `Probability`, `Probability.DependencyType`, `Probability.ExclusiveEnumerationStatus` |

Every function calls one Numerics member. Arrays cross into .NET with one bulk copy, so the cost
of a call does not grow with a per-element crossing, and results are the Numerics bits.

## Descriptive statistics

```python
import numpy as np
from pynumerics import statistics

x = np.random.default_rng(1).gamma(2.5, 4.0, 200)
print(statistics.mean(x), statistics.variance(x), statistics.standard_deviation(x))
print(statistics.skewness(x), statistics.kurtosis(x))          # bias corrected, as scipy bias=False
print(statistics.minimum(x), statistics.maximum(x), statistics.sum(x))
print(statistics.geometric_mean(x), statistics.harmonic_mean(x))
print(statistics.population_variance(x), statistics.coefficient_of_variation(x))
print(statistics.mean_variance(x))                              # (mean, sample variance) in one pass
print(statistics.product_moments(x))                            # mean, sd, skew, excess kurtosis
print(statistics.five_number_summary(x), statistics.seven_number_summary(x))
print(statistics.percentile(x, [0.05, 0.5, 0.95]))              # Type 7
print(statistics.linear_moments(x))                             # l1, l2, t3, t4
print(statistics.linear_moments(x, nmom=5))                     # l1, l2, t3, t4, t5
print(statistics.plotting_position_percentile(x, [0.9, 0.1], "hazen", weights=np.ones(200)))

w = np.random.default_rng(2).random(200)
print(statistics.mean(x, w), statistics.skewness(x, w, "reliability"))
print(statistics.effective_sample_size(w))                      # Kish (sum w)^2 / sum w^2
```

The jackknife takes the statistic as a Python callable, C# source or a .NET delegate. Numerics
evaluates the leave-one-out samples in parallel; C# source avoids a Python call per sample.

```python
import numpy as np
from pynumerics import statistics

x = np.random.default_rng(3).normal(size=50)
print(statistics.jackknife_standard_error(x, np.mean))
print(statistics.jackknife_sample(x, np.median)[:5])
```

## Correlation, autocorrelation and transformations

```python
import numpy as np
from pynumerics import statistics

rng = np.random.default_rng(4)
x = rng.normal(size=100)
y = 0.6 * x + rng.normal(size=100)
print(statistics.pearson(x, y), statistics.spearman(x, y), statistics.kendalls_tau(x, y))
print(statistics.covariance(x, y))

series = np.cumsum(rng.normal(size=120)) * 0.1 + rng.normal(size=120)
lags, acf = statistics.autocorrelation(series, 10)
lower, upper = statistics.autocorrelation_confidence_interval(series.size)
print(lags[acf > upper])
```

`autocorrelation(series, 10, centered=False)` takes the mean as zero, for a process whose mean is
known to be zero. `autocorrelation` refuses a series that is not one-dimensional, holds NaN or infinity, or has
fewer than two values, and a `max_lag` that is not an integer in [1, n - 1], each with a
`ValueError` naming the argument.

```python
import numpy as np
from pynumerics import statistics

flows = np.random.default_rng(5).lognormal(5.0, 0.6, 80)
ok, lam = statistics.box_cox_fit_lambda(flows)
z = statistics.box_cox_transform(flows, lam)
back = statistics.box_cox_inverse_transform(z, lam)
print(ok, lam, np.max(np.abs(back - flows)))

signed = np.random.default_rng(6).normal(1.0, 3.0, 80)
lam = statistics.yeo_johnson_fit_lambda(signed)
print(lam, statistics.yeo_johnson_transform(signed, lam)[:3])
```

## Streaming statistics and histograms

```python
import numpy as np
from pynumerics import statistics

rs = statistics.RunningStatistics()
rs.push(np.arange(10.0))
rs.push(42.0)
print(rs.count, rs.mean, rs.standard_deviation, rs.population_skewness)

cov = statistics.RunningCovarianceMatrix(2)
cov.push(np.random.default_rng(7).normal(size=(500, 2)))       # one Welford update per row
print(cov.n, cov.mean, cov.sample_covariance)                  # shrunk toward the identity

hist = statistics.Histogram(np.random.default_rng(8).normal(size=1000), 12)
lower, upper, counts = hist.bins()
print(counts, hist.mean, hist.mode, hist[0])
hist.add_data([5.0, -5.0])                                     # widens the first and last bins
```

## Given-data sensitivity

```python
import numpy as np
from pynumerics import statistics

rng = np.random.default_rng(9)
x1, x2 = rng.random(5000), rng.random(5000)
y = x1 + 2.0 * x2
print(statistics.first_order_sobol(x1, y), statistics.first_order_sobol(x2, y))
print(statistics.pawn_median(x1, y), statistics.borgonovo_delta(x2, y))
```

## Hypothesis tests

All functions return p-values except `hartigan_dip`, `two_sample_kolmogorov_smirnov` (the statistic
only), `f_test_models` (statistic and p-value) and the Grubbs-Beck tests.

```python
import numpy as np
from pynumerics import hypothesis_tests as ht

rng = np.random.default_rng(10)
a, b = rng.normal(10.0, 2.0, 30), rng.normal(11.0, 2.5, 40)
print(ht.one_sample_t_test(a, 10.0), ht.equal_variance_t_test(a, b), ht.unequal_variance_t_test(a, b))
print(ht.paired_t_test(a, a + rng.normal(0.3, 0.5, 30)), ht.f_test(a, b))
print(ht.jarque_bera_test(b), ht.mann_whitney_test(a, b))
print(ht.two_sample_kolmogorov_smirnov(a, b))                 # exact D = sup |F1 - F2|, no p-value

trend = np.arange(40.0) * 0.05 + rng.normal(size=40)
print(ht.mann_kendall_test(trend), ht.linear_trend_test(np.arange(40.0), trend))
print(ht.wald_wolfowitz_test(trend), ht.ljung_box_test(trend, 5))
print(ht.f_test_models(sse_restricted=120.0, sse_full=100.0, df_restricted=38, df_full=37))

mixed = np.concatenate([rng.normal(0.0, 1.0, 60), rng.normal(6.0, 1.0, 60)])
print(ht.hartigan_dip(mixed), ht.hartigan_dip_test(mixed), ht.unimodality_test(mixed))
```

Low-outlier screening for flood series (Bulletin 17C):

```python
import numpy as np
from pynumerics import hypothesis_tests as ht

peaks = np.concatenate([np.random.default_rng(11).lognormal(8.0, 0.4, 60), [150.0, 300.0]])
print(ht.multiple_grubbs_beck_test(peaks))       # number of potentially influential low floods
print(ht.grubbs_beck_test(peaks))                # (high, low) thresholds of the original test
```

## Goodness of fit

Paired metrics take the observed values first and the modeled values second. The
distribution statistics take any `pynumerics.distributions` object and sort the sample unless
`data_is_sorted=True`.

```python
import numpy as np
from pynumerics import distributions, statistics
from pynumerics import goodness_of_fit as gof

rng = np.random.default_rng(12)
observed = rng.gamma(3.0, 10.0, 60) + 5.0
modeled = observed * 0.9 + rng.normal(0.0, 4.0, 60) + 3.0
print(gof.nash_sutcliffe_efficiency(observed, modeled), gof.kling_gupta_efficiency(observed, modeled))
print(gof.rmse(observed, modeled), gof.mae(observed, modeled), gof.pbias(observed, modeled))  # pbias > 0: underestimates
print(gof.index_of_agreement(observed, modeled), gof.r_squared(observed, modeled))

model = distributions.Gumbel(100.0, 20.0)
sample = rng.gumbel(100.0, 20.0, 75)
print(gof.kolmogorov_smirnov(sample, model), gof.anderson_darling(sample, model),
      gof.chi_squared(sample, model))
print(gof.rmse(sample, model))                    # against quantiles at Weibull positions

ll = model.log_likelihood(sample)
print(gof.aic(2, ll), gof.aicc(sample.size, 2, ll), gof.bic(sample.size, 2, ll))
print(gof.aic_weights([250.0, 252.5, 260.0]))

u = np.sort(rng.random((40, 500)), axis=0)       # 500 samples of 40 probabilities, one per column
print(gof.column_kolmogorov_smirnov_uniform(u)[:3], gof.column_anderson_darling_uniform(u)[:3])

grid = rng.gamma(2.0, 30.0, (60, 2000))           # one row per year, one column per cell
variances, excluded = statistics.column_variances(grid)
moments, _ = statistics.column_linear_moments(grid)
print(variances[:3], moments[1, :3] / moments[0, :3])   # L-CV of the first three cells
```

The column members exclude NaN and infinite entries, return the number excluded from each
column, and give each column bitwise the one-sample result, evaluating all columns in one call.
When bitwise agreement with Numerics is not needed, NumPy can be faster on small matrices: up to
tens of thousands of entries for the variance and standard deviation, and up to a few hundred
thousand for the Kolmogorov-Smirnov statistic. The Numerics kernels themselves are at parity there:
on a 60 x 2000 matrix the .NET Kolmogorov-Smirnov computation takes about 0.15 ms, against 0.14 ms
for NumPy. Most of the remaining cost of the Python call is copying the array into a .NET matrix,
about 0.2 ms at that size. On larger matrices the parallel Numerics members are faster.

## Probability combination

```python
import numpy as np
from pynumerics import multivariate
from pynumerics import probability as pr

p = np.array([0.1, 0.25, 0.4])
corr = np.array([[1.0, 0.5, 0.3], [0.5, 1.0, 0.4], [0.3, 0.4, 1.0]])

print(pr.a_and_b(0.1, 0.25, 0.5), pr.a_or_b(0.1, 0.25, 0.5))
print(pr.independent_union(p), pr.union(p, "perfectly_positive"), pr.union_single_factor(p, 0.35))
print(pr.joint_probability(p), pr.joint_probability_hpcm(p, [1, 1, 1], corr))
print(pr.union_pcm(p, corr))
print(pr.union_mvn(p, multivariate.MultivariateNormal(np.zeros(3), corr)))

events, indicators, status = pr.independent_exclusive_lazy(p)
print(events.sum(), indicators.shape, status)
print(pr.common_cause_adjustment(p, corr))
```
