# Sampling, the bootstrap and MCMC

pynumerics wraps `Numerics.Sampling` in two modules:

| Module | Numerics types |
| --- | --- |
| `pynumerics.mcmc` | `RWMH`, `ARWMH`, `DEMCz`, `DEMCzs`, `HMC`, `NUTS`, `Gibbs`, `SNIS`, `MCMCSampler`, `MCMCResults`, `ParameterResults`, `ParameterStatistics`, `MCMCDiagnostics` |
| `pynumerics.sampling` | `Bootstrap<TData>` and its results, `SobolSequence`, `LatinHypercube`, `MersenneTwister`, `Stratify`, `StratificationBin`, `StratificationOptions` |

Every class wraps one .NET object (`.net`). Chains, draws and replicates come back as NumPy arrays.
Every run is reproducible: the samplers and the bootstrap draw from Mersenne Twister streams seeded
by `seed` (Numerics default 12345), one per chain or replicate, so a seeded run gives the same
numbers whether or not it runs on several threads.

## Functions you pass in

The samplers take a log-likelihood, and the bootstrap takes resample, fit and statistic functions.
Each can be:

- a Python callable. Arrays reach it as NumPy copies. Each call crosses from .NET into Python
  (about 2.5 microseconds plus the function itself) and holds the GIL, so the samplers set
  `ParallelizeChains = False` for you;
- C# source, compiled once with the bundled Roslyn compiler. It runs entirely in .NET and the
  chains stay parallel;
- a .NET delegate, for example from `pynumerics.mcmc.distribution_log_likelihood`.

## MCMC

The priors bound each parameter and seed the chains. The log-likelihood is the full log-posterior up
to a constant: Numerics adds no prior density, so with uniform priors the log-likelihood alone is
enough, and with other priors add their log-densities.

```python
import numpy as np
from scipy import stats
from pynumerics import distributions, mcmc

data = stats.norm(100.0, 15.0).rvs(size=50, random_state=1)

def log_likelihood(theta):
    mu, sigma = theta
    return float(stats.norm(mu, sigma).logpdf(data).sum())

priors = [distributions.Uniform(50.0, 150.0), distributions.Uniform(1.0, 50.0)]
sampler = mcmc.RWMH(priors, log_likelihood, initialize="MAP", iterations=2000, warmup_iterations=500,
                    thinning_interval=1, output_length=4000)
results = sampler.sample()
print(sampler.parallelize_chains)          # False: the callback is Python
for name, result in zip(["mu", "sigma"], results.parameter_results):
    s = result.summary_statistics
    print(name, round(s.mean, 2), round(s.lower_ci, 2), round(s.upper_ci, 2), round(s.rhat, 3), round(s.ess))
print(sampler.markov_chains.shape)         # (chains, iterations, parameters)
```

`markov_chains` (on a sampler or on `MCMCResults`) and `markov_chain_log_likelihoods` read the
chains from .NET once per run (about 3 microseconds per recorded state) and keep them; each later
access returns a fresh NumPy copy. Running the sampler again refreshes them.

`RWMH` takes a proposal covariance; the default zero matrix is replaced by the MAP covariance with
`initialize="MAP"`, as in the Numerics examples. The other samplers need no proposal:

```python
import numpy as np
from pynumerics import distributions, mcmc

def log_posterior(theta):
    # Beta(2, 3) prior and 7 successes in 20 trials: the posterior is Beta(9, 16)
    p = theta[0]
    return 8.0 * np.log(p) + 15.0 * np.log1p(-p) if 0.0 < p < 1.0 else -np.inf

settings = dict(iterations=1000, warmup_iterations=200, thinning_interval=1, output_length=2000)
prior = [distributions.Uniform(0.0, 1.0)]
for sampler in (mcmc.ARWMH(prior, log_posterior, **settings),
                mcmc.DEMCzs(prior, log_posterior, **settings),
                mcmc.HMC(prior, log_posterior, step_size=0.05, steps=5, **settings),
                mcmc.NUTS(prior, log_posterior, max_tree_depth=7, **settings)):
    summary = sampler.sample().parameter_results[0].summary_statistics
    print(type(sampler).__name__, round(summary.mean, 3))   # about 0.36 = 9 / 25
```

`HMC` and `NUTS` take an optional `gradient` (finite differences by default), and `NUTS` takes
`log_likelihood_and_gradient` returning `(value, gradient)` in one call. `Gibbs` takes a proposal
`proposal(values, prng)` that draws from the full conditionals with the chain's .NET generator, and
`SNIS` draws independent points from the priors (or an importance distribution) and resamples them;
like the other samplers it evaluates them in parallel only while `parallelize_chains` is True, so a
Python likelihood runs on the calling thread with the same bits.

### Likelihoods that run in .NET

`distribution_log_likelihood` builds the log-likelihood of data under any Numerics distribution, so
no Python runs per evaluation and the chains keep running in parallel. Parameters not sampled keep
the template's values, and `add_log_prior=True` adds the priors' log-densities:

```python
from scipy import stats
from pynumerics import distributions, mcmc

data = stats.norm(3.0, 2.0).rvs(size=30, random_state=2)
prior = distributions.Normal(0.0, 5.0)
model = mcmc.distribution_log_likelihood(distributions.Normal(0.0, 2.0), data, priors=[prior],
                                         parameter_indices=[0], add_log_prior=True)   # sigma = 2 is known
sampler = mcmc.ARWMH([prior], model, iterations=2000, warmup_iterations=500, thinning_interval=1,
                     output_length=4000)
print(sampler.parallelize_chains)          # True
print(round(sampler.sample().parameter_results[0].summary_statistics.mean, 2))
```

C# source works the same way; it can use any Numerics type:

```python
from pynumerics import distributions, mcmc

sampler = mcmc.DEMCz([distributions.Uniform(0.0, 1.0)],
                     "p => (p[0] <= 0 || p[0] >= 1) ? double.NegativeInfinity"
                     " : 8 * Math.Log(p[0]) + 15 * Math.Log(1 - p[0])",
                     iterations=1000, warmup_iterations=200, thinning_interval=1, output_length=2000)
print(round(sampler.sample().parameter_results[0].summary_statistics.mean, 3))
```

### Results and diagnostics

`MCMCResults` holds the chains, the output draws, the MAP, the posterior mean and one
`ParameterResults` per parameter (summary statistics, kernel density, histogram, autocorrelation).
R-hat and ESS are the rank-normalized forms of Vehtari et al. (2021), the values ArviZ reports:

```python
from pynumerics import distributions, mcmc

sampler = mcmc.DEMCz([distributions.Uniform(0.0, 1.0)], "p => 8 * Math.Log(p[0]) + 15 * Math.Log(1 - p[0])",
                     iterations=1000, warmup_iterations=200, thinning_interval=1, output_length=2000)
sampler.sample()
results = sampler.results(alpha=0.05)      # 95 percent intervals
print(results.map, results.posterior_mean)
print(mcmc.gelman_rubin(sampler, warmup_iterations=sampler.warmup_iterations))
ess, autocorrelation = mcmc.effective_sample_sizes(sampler)
print(ess, autocorrelation[0][:3])
print(mcmc.effective_sample_size(results.output[:, 0]))
print(mcmc.minimum_sample_size(0.975, 0.005, 0.95))   # Raftery-Lewis: 3700
saved = results.to_bytes()
print(mcmc.MCMCResults.from_bytes(saved).output.shape)
```

`gelman_rubin` and `effective_sample_sizes` also take any array of shape
`(chains, iterations, parameters)`, so chains from other tools can be checked the same way.

## The bootstrap

`Bootstrap` wraps `Bootstrap<TData>`. With NumPy data `TData` is `double[]`. `resample(data,
parameters, prng)` draws a data set with the replicate's .NET generator, `fit(data)` returns the
fitted parameters and `statistic(parameters)` the quantities of interest:

```python
import numpy as np
from scipy import stats
from pynumerics import sampling

x = stats.norm(10.0, 2.0).rvs(size=40, random_state=5)

def resample(data, parameters, prng):
    return data[[prng.Next(data.size) for _ in range(data.size)]]

boot = sampling.Bootstrap(x, [x.mean()], resample=resample, fit=lambda d: [d.mean()],
                          statistic=lambda p: p, jackknife=lambda d, i: np.delete(d, i),
                          sample_size=lambda d: d.size, replicates=2000)
boot.run()
print(boot.use_parallel)                    # False: the functions are Python
# statistic_results use the requested method; parameter_results are always percentile intervals
for method in ("percentile", "bias_corrected", "bca", "normal"):
    interval = boot.get_confidence_intervals(method, alpha=0.1).statistic_results[0]
    print(method, round(interval.lower_ci, 3), round(interval.upper_ci, 3), round(interval.standard_error, 3))
```

Python functions set `use_parallel = False`, so the replicates run in order on the calling thread
(each with its own seed, so the results are the bits of a parallel run). The leave-one-out jackknife
of the BCa intervals follows the same switch, with the same bits either way. C# source keeps the
replicates in parallel, much faster:

```python
from scipy import stats
from pynumerics import sampling

x = stats.norm(10.0, 2.0).rvs(size=40, random_state=5)
fast = sampling.Bootstrap(
    x, [x.mean()],
    resample="(d, p, r) => { var o = new double[d.Length]; for (int i = 0; i < d.Length; i++) "
             "o[i] = d[r.Next(d.Length)]; return o; }",
    fit="d => new Numerics.Mathematics.Optimization.ParameterSet(new[] { d.Average() }, double.NaN)",
    statistic="p => p.Values", replicates=2000)
print(fast.run().get_confidence_intervals("percentile").statistic_results[0])
```

`run_double_bootstrap`, `run_with_studentized_bootstrap` (for `"bootstrap_t"`) and the
covariance-aware `run_pivotal_bootstrap` (with `fit_with_covariance` returning
`(parameters, covariance)` and a parent `BootstrapFit`) follow the same pattern.

## Designs and generators

```python
from pynumerics import distributions, sampling

print(sampling.SobolSequence(2).next_doubles(4))              # unscrambled; scipy's points 1 to 4
print(sampling.SobolSequence(2, seed=7).next_doubles(2))      # scrambled and shifted, reproducible
print(sampling.latin_hypercube(5, 2, seed=12345))
twister = sampling.MersenneTwister(12345)
print(twister.next_doubles(3), twister.gen_rand_int32())

bins = sampling.stratify_x_values(0.0, 10.0, 5)
weights = sampling.x_to_probability(bins, distributions.Normal(5.0, 2.0).cdf).weight
print(weights, weights.sum())
probabilities = sampling.stratify_probabilities(0.001, 0.999, 4, importance="normal")
print(sampling.probability_to_x(probabilities, distributions.Normal(0.0, 1.0).inverse_cdf).upper)
print(sampling.StratificationBins.from_xml(probabilities.to_xml()).lower)
```
