# Univariate distributions and uncertainty analysis

pynumerics wraps `Numerics.Distributions` univariate families in two modules:

| Module | Numerics types |
| --- | --- |
| `pynumerics.distributions` | `UnivariateDistributionBase` and every univariate family, `TruncatedDistribution`, `CompetingRisks`, `UnivariateDistributionFactory`, the estimation and standard-error interfaces |
| `pynumerics.bootstrap` | `BootstrapAnalysis`, `UncertaintyAnalysisResults` |

Every class wraps one .NET object (`.net`). Evaluation members take a scalar or any array: an
array crosses into .NET once as a `System.Double[]`, so no numerical work happens in Python and the
values are the Numerics bits. Parameters follow the Numerics conventions: the GEV, GPA, GLO, GNO and
kappa shapes use Hosking's sign (negative is heavy-tailed), `LogNormal` and `LogPearsonTypeIII` take
log-space parameters in a base (10 by default), and `LnNormal` takes the real-space mean and
standard deviation.

## Families

Each typed class takes its parameters by name, with the Numerics defaults, and exposes them as
properties under the same names.

| Class | Parameters |
| --- | --- |
| `Normal` | `mean`, `standard_deviation` |
| `LogNormal` | `mean_of_log`, `standard_deviation_of_log`, `base` |
| `LnNormal` | `mean`, `standard_deviation` (real space); `mean_of_log`, `standard_deviation_of_log` |
| `GeneralizedExtremeValue`, `GeneralizedPareto`, `GeneralizedLogistic`, `GeneralizedNormal` | `location`, `scale`, `shape` |
| `KappaFour` | `location`, `scale`, `shape`, `shape2` |
| `Gumbel`, `Exponential`, `Logistic`, `Cauchy` | `location`, `scale` |
| `PearsonTypeIII` | `mean`, `standard_deviation`, `skew`; gamma form `support_bound`, `scale`, `shape` |
| `LogPearsonTypeIII` | `mean_of_log`, `standard_deviation_of_log`, `skew_of_log`, `base`; `log_space()` |
| `GammaDistribution` | `scale`, `shape` (and `rate`) |
| `Weibull`, `Pareto`, `InverseGamma` | `scale`, `shape` |
| `Rayleigh` | `scale` |
| `StudentT` | `location`, `scale`, `degrees_of_freedom` |
| `NoncentralT` | `degrees_of_freedom`, `noncentrality` |
| `ChiSquared` | `degrees_of_freedom` |
| `InverseChiSquared` | `degrees_of_freedom`, `scale` (τ²) |
| `Beta` | `alpha`, `beta` |
| `GeneralizedBeta` | `alpha`, `beta`, `minimum`, `maximum`; `GeneralizedBeta.pert(...)` |
| `Uniform` | `minimum`, `maximum` |
| `Triangular`, `Pert` | `minimum`, `most_likely`, `maximum` |
| `PertPercentile`, `PertPercentileZ` | `percentile_5th`, `percentile_50th`, `percentile_95th` |
| `TruncatedNormal` | `mu`, `sigma`, `minimum`, `maximum` |
| `TruncatedDistribution` | `distribution`, `lower_bound`, `upper_bound` |
| `VonMises` | `mu`, `kappa` |
| `Deterministic` | `value` |
| `Poisson` | `rate` |
| `Bernoulli` | `probability` |
| `Binomial` | `probability`, `number_of_trials` |
| `Geometric` | `probability` |
| `UniformDiscrete` | `minimum`, `maximum` |
| `Mixture` | `weights`, `distributions` |
| `CompetingRisks` | `distributions` |
| `Empirical`, `KernelDensity` | value tables and samples |

```python
import numpy as np
from pynumerics import distributions

gev = distributions.GeneralizedExtremeValue(100.0, 15.0, -0.2)
p = np.array([0.5, 0.9, 0.99])
print(gev.inverse_cdf(p))                     # quantiles
print(gev.inverse_ccdf(1.0 - p))              # survival quantiles
print(gev.pdf(120.0), gev.cdf(120.0), gev.ccdf(120.0), gev.hf(120.0))
print(gev.log_pdf([110.0, 150.0]), gev.log_cdf(150.0), gev.log_ccdf(150.0))
print(gev.mean, gev.median, gev.mode, gev.variance, gev.skewness, gev.kurtosis)
print(gev.location, gev.scale, gev.shape, gev.parameters)

gamma = distributions.GammaDistribution(scale=2.0, shape=3.5)
binomial = distributions.Binomial(0.3, 20)
print(gamma.inverse_cdf(0.99), binomial.pdf(5.0), binomial.log_pdf(5.0))

natural = distributions.LogNormal(1.5, 0.4, base=np.e)
lp3 = distributions.LogPearsonTypeIII(2.0, 0.25, -0.3)
print(natural.base, lp3.support_bound, lp3.log_space())

print(distributions.standard_z(0.975), distributions.standard_cdf(1.96), distributions.standard_log_pdf(0.0))
```

`inverse_ccdf(q)` returns the x with `ccdf(x) = q`. The generalized Pareto distribution has a
closed form in ln q that never forms `1 - q`, which keeps the precision of small exceedance
probabilities; the other families return `inverse_cdf(1 - q)`.

## Fitting and estimation

`estimate` fits in place; `mle`, the moment methods and `bootstrap` return values and leave the
distribution unchanged. The moment methods also work on the class, without an instance.

```python
import numpy as np
from pynumerics import distributions, statistics

sample = distributions.Gumbel(50.0, 8.0).generate_random_values(200, seed=5)

fit = distributions.Gumbel()
fit.estimate(sample, "maximum_likelihood")
print(fit.parameters, fit.log_likelihood(sample))
print(distributions.Gumbel().mle(sample))
initial, lower, upper = distributions.Gumbel().parameter_constraints(sample)

moments = statistics.linear_moments(sample)
print(distributions.GeneralizedExtremeValue.parameters_from_linear_moments(moments))
print(distributions.GeneralizedPareto.parameters_from_linear_moments(moments))
print(distributions.KappaFour.linear_moments_from_parameters([10.0, 2.0, 0.1, 0.2]))
print(distributions.Gumbel.moments_from_parameters([50.0, 8.0]))

gpa = distributions.GeneralizedPareto.fit_known_location(sample[sample > 60.0], 60.0)
replicate = fit.bootstrap("maximum_likelihood", 200, seed=12)
print(gpa, replicate)
```

## Standard errors

Families that implement `IStandardError` give the asymptotic parameter covariance, the quantile
gradient and Jacobian, and the delta-method quantile variance.

```python
from pynumerics import distributions

gev = distributions.GeneralizedExtremeValue(100.0, 10.0, -0.1)
print(gev.parameter_covariance(80, "maximum_likelihood"))
print(gev.quantile_variance(0.99, 80, "maximum_likelihood"))
print(gev.quantile_gradient(0.99))
jacobian, determinant = gev.quantile_jacobian([0.1, 0.5, 0.9])
print(determinant, gev.log_abs_quantile_jacobian([0.1, 0.5, 0.9]))
print(gev.expected_information_matrix(80))
```

## Peaks over threshold

```python
from pynumerics import distributions

# Annual maximum of a Poisson(3.7) count of GPD(10, 2, -0.15) exceedances
gev = distributions.GeneralizedExtremeValue.from_poisson_generalized_pareto(3.7, 10.0, 2.0, -0.15)
gpd = distributions.GeneralizedPareto(10.0, 2.0, -0.15)
print(gev, gev.cdf(20.0), gpd.inverse_ccdf(1e-6))

# The inverse: the exceedance GPD of an annual-maximum GEV at 3.7 exceedances per year
back = distributions.GeneralizedPareto.from_poisson_generalized_extreme_value(3.7, gev.location, gev.scale, gev.shape)
print(back.location, back.scale, back.shape)   # 10.0, 2.0, -0.15 to rounding
```

## Many parameter sets

```python
import numpy as np
from pynumerics import distributions

gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
rows = np.array([[4.0, 1.0, -0.1], [5.0, 2.0, 0.1], [4.0, -1.0, 0.0]])
quantiles, valid = gev.parameter_set_quantile_table(rows, [0.5, 0.99])   # NaN rows where invalid
print(valid, gev.parameter_set_quantiles(rows[valid], 0.99))
try:
    gev.parameter_set_quantile_table(rows, [0.5, 0.99], validate=True)
except ValueError as error:
    print(error)                                # Parameter set 2: The scale parameter ... must be positive.
```

## Growth curves and per-element counts

`inverse_cdf` on an array is a growth curve: many probabilities under one parameter set, each bitwise
equal to the scalar call. `poisson_inverse_cdf`, `poisson_log_pmf` and `poisson_cdf` take one rate
per element, and `binomial_log_pmf` and `binomial_cdf` one success probability and number of trials
per element; the arguments broadcast, and each value is bitwise equal to the same call on a
`Poisson` or `Binomial` instance. Arrays of 4096 or more elements are evaluated in parallel in
.NET, and the Poisson quantile search starts from a Cornish-Fisher
estimate; neither changes a result.

```python
import numpy as np
from pynumerics import distributions

regional = distributions.GeneralizedExtremeValue(0.8376, 0.2967, 0.0312)
print(regional.inverse_cdf([0.5, 0.9, 0.99, 0.999]))      # the regional growth curve

rates = np.array([0.8, 2.5, 12.0])
print(distributions.poisson_inverse_cdf(rates, [0.5, 0.9, 0.99]))   # [ 1.  5. 21.]
print(distributions.poisson_log_pmf(rates, [0, 3, 12]), distributions.poisson_cdf(rates, 4))
print(distributions.binomial_log_pmf([0.3, 0.5], [10, 4], [3, 2]), distributions.binomial_cdf(0.3, 10, np.arange(11)))
```

## Truncation, competing risks, mixtures and tables

```python
import numpy as np
from pynumerics import distributions

truncated = distributions.TruncatedDistribution(distributions.Normal(1.0, 2.0), -1.0, 4.0)
risks = distributions.CompetingRisks([distributions.Weibull(10.0, 2.0), distributions.Exponential(0.0, 8.0)])
mixture = distributions.Mixture([0.3, 0.7], [distributions.Normal(), distributions.Normal(3.0, 1.0)])
print(truncated.cdf(2.0), risks.cdf(5.0), mixture.inverse_cdf(0.5))
print(len(risks.cumulative_incidence_functions()))

p = np.linspace(0.01, 0.99, 25)
table = distributions.Empirical(distributions.Normal(5.0, 1.0).inverse_cdf(p), p)
print(distributions.Empirical.convolve([table, table]).mean)
kde = distributions.KernelDensity(distributions.Normal().generate_random_values(100, seed=3))
print(kde.bandwidth, kde.pdf(0.0))
```

## Moments, conditional moments and plotting tables

```python
from pynumerics import distributions, sampling

gumbel = distributions.Gumbel(50.0, 8.0)
print(gumbel.central_moments(), gumbel.central_moments(steps=1000))
print(gumbel.conditional_mean(45.0, 70.0), gumbel.conditional_expected_value(0.9))
print(gumbel.log_likelihood_right_censored(70.0, 2))
print(gumbel.create_pdf_graph().shape, gumbel.create_cdf_graph(sampling.stratify_x_values(30.0, 90.0, 12)).shape)
```

## Values, families and serialization

Distributions compare by value; `frozen()` gives a read-only, hashable copy. `FAMILIES` maps each
Numerics type name to its class; `create_distribution`, `try_create_distribution` and `from_xml`
build distributions by name or from XML.

```python
from pynumerics import distributions

weibull = distributions.Weibull(10.0, 2.5)
print(weibull == distributions.Weibull(10.0, 2.5), hash(weibull.frozen()) is not None)
print(distributions.FAMILIES["GammaDistribution"].parameter_names)
print(distributions.create_distribution("Weibull", [10.0, 2.5]) == weibull)
print(distributions.try_create_distribution("Mixture"))     # None: needs components
print(distributions.from_xml(weibull.to_xml()) == weibull)
```

## Bootstrap uncertainty

```python
from pynumerics import bootstrap, distributions

gev = distributions.GeneralizedExtremeValue(100.0, 10.0, -0.1)
analysis = bootstrap.BootstrapAnalysis(gev, "maximum_likelihood", 50, replications=500, seed=7)
results = analysis.estimate([0.5, 0.9, 0.99])
print(results.mode_curve, results.confidence_intervals, results.mean_curve)
print(analysis.normal_quantile_ci([0.99]), analysis.bias_corrected_quantile_ci([0.99]))
values, fitness, weight = analysis.parameter_sets()
restored = bootstrap.UncertaintyAnalysisResults.from_bytes(results.to_bytes())
print(values.shape, restored.mean_curve)
```
