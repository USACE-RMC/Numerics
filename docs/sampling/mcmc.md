# MCMC Sampling

[← Previous: Goodness-of-Fit](../statistics/goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](convergence-diagnostics.md)

Markov Chain Monte Carlo (MCMC) methods sample from complex posterior distributions that are difficult to sample directly. The ***Numerics*** library provides multiple MCMC samplers for Bayesian inference, uncertainty quantification, and parameter estimation with full posterior distributions [[1]](#1).

## Available MCMC Samplers

| Sampler | Full Name | Best For | Key Features |
|---------|-----------|----------|--------------|
| **RWMH** | Random Walk Metropolis-Hastings | General purpose, small dimensions | Simple, robust baseline |
| **ARWMH** | Adaptive Random Walk M-H | Medium dimensions (2-20) | Self-tuning proposal |
| **DEMCz** | Differential Evolution MCMC | High dimensions, multimodal | Population-based, efficient |
| **DEMCzs** | DE-MCMC with snooker update | Very high dimensions | Enhanced DE-MCMC |
| **HMC** | Hamiltonian Monte Carlo | Smooth posteriors | Uses gradient information |
| **Gibbs** | Gibbs Sampler | Conditional distributions available | No rejections |

## Common MCMC Interface

All samplers inherit from `MCMCSampler` base class with common properties:

```cs
// Configuration
int PRNGSeed               // Random seed (default: 12345)
int InitialIterations      // Initialization phase (default: 10)
int WarmupIterations       // Burn-in period (default: 1750)
int Iterations             // Main sampling (default: 3500)
int NumberOfChains         // Parallel chains (default: 4)
int ThinningInterval       // Keep every nth sample (default: 20)

// Inputs
List<IUnivariateDistribution> PriorDistributions
LogLikelihood LogLikelihoodFunction

// Outputs (after sampling)
ParameterSet[] ParameterSets      // All samples
double[] LogLikelihoods           // Log-likelihood values
double[] LogPosteriors            // Log-posterior values
int[] SampleCount                 // Samples per chain
```

## Defining the Model

### Step 1: Define Prior Distributions

```cs
using Numerics.Distributions;
using Numerics.Sampling.MCMC;

// Example: Linear regression y = a + b*x + ε
// Parameters: [a (intercept), b (slope), σ (noise)]

var priors = new List<IUnivariateDistribution>
{
    new Normal(0, 10),        // Intercept: N(0, 10)
    new Normal(0, 10),        // Slope: N(0, 10)
    new Uniform(0.1, 5.0)     // Noise std dev: Uniform(0.1, 5)
};

Console.WriteLine("Prior Distributions:");
for (int i = 0; i < priors.Count; i++)
{
    Console.WriteLine($"  θ{i}: {priors[i].DisplayName}");
}
```

### Step 2: Define Log-Likelihood Function

The log-likelihood function computes the log-probability of the data given parameters:

```cs
// Observed data
double[] xData = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
double[] yData = { 2.5, 4.8, 6.2, 8.5, 10.1, 12.8, 14.2, 16.5, 18.1, 20.3 };

// Define log-likelihood function
LogLikelihood logLikelihood = (parameters) =>
{
    double a = parameters[0];     // Intercept
    double b = parameters[1];     // Slope
    double sigma = parameters[2]; // Noise std dev
    
    // Prior log-likelihood
    double logPrior = priors[0].LogPDF(a) + 
                     priors[1].LogPDF(b) + 
                     priors[2].LogPDF(sigma);
    
    // Data log-likelihood: y ~ N(a + b*x, σ²)
    double logData = 0;
    for (int i = 0; i < xData.Length; i++)
    {
        double predicted = a + b * xData[i];
        double residual = yData[i] - predicted;
        logData += -0.5 * Math.Log(2 * Math.PI * sigma * sigma) - 
                   0.5 * residual * residual / (sigma * sigma);
    }
    
    return logPrior + logData;
};
```

**Important:** The log-likelihood function should return the sum of:
1. Log-prior probability: `Σ log(p(θᵢ))`
2. Log-data likelihood: `log(p(data|θ))`

This gives the log-posterior: `log(p(θ|data)) ∝ log(p(θ)) + log(p(data|θ))`

## Random Walk Metropolis-Hastings (RWMH)

The simplest and most robust MCMC algorithm [[2]](#2):

```cs
using Numerics.Sampling.MCMC;

// Create sampler
var rwmh = new RWMH(priors, logLikelihood);

// Configure sampling
rwmh.PRNGSeed = 12345;
rwmh.InitialIterations = 50;      // 5x parameters
rwmh.WarmupIterations = 2000;     // Burn-in
rwmh.Iterations = 5000;           // Main sampling
rwmh.NumberOfChains = 4;          // Parallel chains
rwmh.ThinningInterval = 10;       // Keep every 10th sample

// Run sampler
Console.WriteLine("Running RWMH sampler...");
rwmh.Sample();

// Access results
var samples = rwmh.ParameterSets;
Console.WriteLine($"Generated {samples.Length} samples");
Console.WriteLine($"Samples per chain: {string.Join(", ", rwmh.SampleCount)}");

// Posterior statistics
for (int i = 0; i < priors.Count; i++)
{
    var values = samples.Select(s => s.Values[i]).ToArray();
    double mean = values.Average();
    double std = Statistics.StandardDeviation(values);
    double q025 = Statistics.Percentile(values.OrderBy(x => x).ToArray(), 2.5);
    double q975 = Statistics.Percentile(values.OrderBy(x => x).ToArray(), 97.5);
    
    Console.WriteLine($"θ{i}: {mean:F3} ± {std:F3}, 95% CI: [{q025:F3}, {q975:F3}]");
}
```

**When to use RWMH:**
- General purpose baseline
- Low-dimensional problems (< 10 parameters)
- When simplicity and robustness are priorities
- As a reference for comparing other samplers

## Adaptive Random Walk M-H (ARWMH)

ARWMH automatically tunes the proposal distribution during warmup [[3]](#3):

```cs
var arwmh = new ARWMH(priors, logLikelihood);

// Configuration
arwmh.PRNGSeed = 12345;
arwmh.WarmupIterations = 2000;    // Adaptation happens here
arwmh.Iterations = 5000;
arwmh.NumberOfChains = 4;

Console.WriteLine("Running Adaptive RWMH sampler...");
arwmh.Sample();

var samples = arwmh.ParameterSets;
Console.WriteLine($"Generated {samples.Length} samples");

// ARWMH adapts proposal covariance to achieve ~23% acceptance rate
Console.WriteLine("ARWMH automatically tuned proposal during warmup");
```

**When to use ARWMH:**
- Medium-dimensional problems (2-20 parameters)
- When you don't want to manually tune proposals
- Correlated parameters
- Default choice for most applications

**Advantages:**
- No manual tuning required
- Adapts to parameter correlations
- Generally more efficient than fixed RWMH

## Differential Evolution MCMC (DEMCz)

Population-based sampler using differential evolution [[4]](#4):

```cs
var demcz = new DEMCz(priors, logLikelihood);

// Configuration
demcz.PRNGSeed = 12345;
demcz.NumberOfChains = 10;        // More chains for population diversity
demcz.WarmupIterations = 2000;
demcz.Iterations = 5000;

Console.WriteLine("Running DE-MCMC sampler...");
demcz.Sample();

var samples = demcz.ParameterSets;
Console.WriteLine($"Generated {samples.Length} samples from {demcz.NumberOfChains} chains");

// DEMCz is particularly effective for multimodal posteriors
```

**When to use DEMCz:**
- High-dimensional problems (20+ parameters)
- Multimodal posteriors
- Complex posterior geometry
- When ARWMH struggles with convergence

**Advantages:**
- Excellent for high dimensions
- Handles multimodal distributions
- Robust to initialization
- Self-tuning proposals from population

### DEMCz with Snooker Update (DEMCzs)

Enhanced version with improved mixing:

```cs
var demczs = new DEMCzs(priors, logLikelihood);

demczs.NumberOfChains = 12;      // Even more chains recommended
demczs.WarmupIterations = 2000;
demczs.Iterations = 5000;

Console.WriteLine("Running DE-MCMC with snooker update...");
demczs.Sample();

// Snooker update provides better exploration in very high dimensions
```

## Hamiltonian Monte Carlo (HMC)

Uses gradient information for efficient sampling [[5]](#5):

```cs
var hmc = new HMC(priors, logLikelihood);

// HMC-specific settings
hmc.NumberOfChains = 4;
hmc.WarmupIterations = 1000;      // HMC converges faster
hmc.Iterations = 2000;

Console.WriteLine("Running Hamiltonian Monte Carlo...");
hmc.Sample();

// HMC produces high-quality samples with lower autocorrelation
var samples = hmc.ParameterSets;
Console.WriteLine($"Generated {samples.Length} high-quality samples");
```

**When to use HMC:**
- Smooth, differentiable posteriors
- When gradient information is available
- Need for low autocorrelation
- Medium to high dimensions with smooth geometry

**Advantages:**
- Very efficient (low autocorrelation)
- Explores parameter space quickly
- Excellent for smooth posteriors

**Disadvantages:**
- Requires gradient computation
- Less robust to discontinuities
- More complex to tune

## Gibbs Sampler

Samples each parameter conditionally given others:

```cs
var gibbs = new Gibbs(priors, logLikelihood);

gibbs.NumberOfChains = 4;
gibbs.WarmupIterations = 1500;
gibbs.Iterations = 3500;

Console.WriteLine("Running Gibbs sampler...");
gibbs.Sample();

// Gibbs has no rejections - every proposal is accepted
Console.WriteLine("Gibbs sampler completed (no rejections)");
```

**When to use Gibbs:**
- Conditional distributions available in closed form
- Conjugate prior-likelihood pairs
- Hierarchical models

**Note:** Gibbs is most efficient when conditional distributions are easy to sample from. For general problems, ARWMH or DEMCz are often better choices.

## Complete Bayesian Inference Example

### Example 1: Linear Regression with Full Uncertainty

```cs
using Numerics.Distributions;
using Numerics.Sampling.MCMC;
using Numerics.Data.Statistics;

// Generate synthetic data
double trueIntercept = 2.0;
double trueSlope = 1.8;
double trueNoise = 0.5;

var random = new MersenneTwister(123);
int n = 20;
double[] x = Enumerable.Range(1, n).Select(i => (double)i).ToArray();
double[] yTrue = x.Select(xi => trueIntercept + trueSlope * xi).ToArray();
double[] y = yTrue.Select(yi => yi + new Normal(0, trueNoise).InverseCDF(random.NextDouble())).ToArray();

Console.WriteLine($"Generated {n} data points");
Console.WriteLine($"True parameters: a={trueIntercept}, b={trueSlope}, σ={trueNoise}");

// Define priors
var priors = new List<IUnivariateDistribution>
{
    new Normal(0, 10),        // Intercept
    new Normal(0, 10),        // Slope
    new Uniform(0.1, 5.0)     // Noise
};

// Define log-likelihood
LogLikelihood logLik = (theta) =>
{
    double a = theta[0];
    double b = theta[1];
    double sigma = theta[2];
    
    // Prior
    double logPrior = priors[0].LogPDF(a) + priors[1].LogPDF(b) + priors[2].LogPDF(sigma);
    
    // Likelihood
    double logData = 0;
    for (int i = 0; i < n; i++)
    {
        double mu = a + b * x[i];
        logData += new Normal(mu, sigma).LogPDF(y[i]);
    }
    
    return logPrior + logData;
};

// Run MCMC
var sampler = new ARWMH(priors, logLik);
sampler.WarmupIterations = 2000;
sampler.Iterations = 5000;
sampler.NumberOfChains = 4;
sampler.ThinningInterval = 5;

Console.WriteLine("\nRunning MCMC...");
sampler.Sample();

// Analyze results
var samples = sampler.ParameterSets;
Console.WriteLine($"\nPosterior Summary ({samples.Length} samples):");
Console.WriteLine("Parameter | True  | Post Mean | Post SD | 95% Credible Interval");
Console.WriteLine("---------------------------------------------------------------");

string[] names = { "Intercept", "Slope", "Noise SD" };
double[] trueVals = { trueIntercept, trueSlope, trueNoise };

for (int i = 0; i < 3; i++)
{
    var vals = samples.Select(s => s.Values[i]).OrderBy(v => v).ToArray();
    double mean = vals.Average();
    double std = Statistics.StandardDeviation(vals);
    double lower = Statistics.Percentile(vals, 2.5);
    double upper = Statistics.Percentile(vals, 97.5);
    
    Console.WriteLine($"{names[i],-9} | {trueVals[i],5:F2} | {mean,9:F3} | {std,7:F3} | [{lower:F3}, {upper:F3}]");
}

// Posterior predictive
Console.WriteLine("\nPosterior Predictive at x=15:");
double xNew = 15;
var predictions = samples.Select(s => s.Values[0] + s.Values[1] * xNew).ToArray();
double predMean = predictions.Average();
double predSD = Statistics.StandardDeviation(predictions);
double predLower = Statistics.Percentile(predictions.OrderBy(p => p).ToArray(), 2.5);
double predUpper = Statistics.Percentile(predictions.OrderBy(p => p).ToArray(), 97.5);

Console.WriteLine($"E[y|x={xNew}] = {predMean:F2} ± {predSD:F2}");
Console.WriteLine($"95% Credible Interval: [{predLower:F2}, {predUpper:F2}]");
```

### Example 2: Distribution Parameter Estimation

```cs
// Observed data from unknown GEV distribution
double[] annualMaxima = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300 };

Console.WriteLine("Bayesian Estimation of GEV Parameters");
Console.WriteLine("=" + new string('=', 50));

// Prior distributions for GEV parameters [ξ, α, κ]
var priors = new List<IUnivariateDistribution>
{
    new Normal(15000, 5000),      // Location: N(15000, 5000)
    new Uniform(100, 5000),       // Scale: U(100, 5000)
    new Uniform(-0.5, 0.5)        // Shape: U(-0.5, 0.5)
};

// Log-likelihood
LogLikelihood logLik = (theta) =>
{
    double xi = theta[0];
    double alpha = theta[1];
    double kappa = theta[2];
    
    // Check parameter validity
    if (alpha <= 0) return double.NegativeInfinity;
    
    // Prior
    double logPrior = priors[0].LogPDF(xi) + priors[1].LogPDF(alpha) + priors[2].LogPDF(kappa);
    
    // Likelihood
    var gev = new GeneralizedExtremeValue(xi, alpha, kappa);
    if (!gev.ParametersValid) return double.NegativeInfinity;
    
    double logData = annualMaxima.Sum(x => gev.LogPDF(x));
    
    return logPrior + logData;
};

// Sample posterior
var sampler = new ARWMH(priors, logLik);
sampler.WarmupIterations = 3000;
sampler.Iterations = 10000;
sampler.NumberOfChains = 4;
sampler.ThinningInterval = 10;

Console.WriteLine("Sampling posterior distribution...");
sampler.Sample();

var samples = sampler.ParameterSets;
Console.WriteLine($"Generated {samples.Length} posterior samples\n");

// Parameter estimates
string[] paramNames = { "Location (ξ)", "Scale (α)", "Shape (κ)" };
for (int i = 0; i < 3; i++)
{
    var vals = samples.Select(s => s.Values[i]).OrderBy(v => v).ToArray();
    Console.WriteLine($"{paramNames[i]}:");
    Console.WriteLine($"  Mean: {vals.Average():F2}");
    Console.WriteLine($"  Median: {Statistics.Percentile(vals, 50):F2}");
    Console.WriteLine($"  95% CI: [{Statistics.Percentile(vals, 2.5):F2}, {Statistics.Percentile(vals, 97.5):F2}]");
}

// Posterior predictive quantiles
Console.WriteLine("\nPosterior Predictive 100-year Flood:");
var q100 = samples.Select(s =>
{
    var dist = new GeneralizedExtremeValue(s.Values[0], s.Values[1], s.Values[2]);
    return dist.InverseCDF(0.99);
}).OrderBy(q => q).ToArray();

Console.WriteLine($"  Mean: {q100.Average():F0} cfs");
Console.WriteLine($"  Median: {Statistics.Percentile(q100, 50):F0} cfs");
Console.WriteLine($"  95% CI: [{Statistics.Percentile(q100, 2.5):F0}, {Statistics.Percentile(q100, 97.5):F0}] cfs");
```

## Thinning

Thinning reduces autocorrelation by keeping only every nth sample:

```cs
// Without thinning
sampler.ThinningInterval = 1;     // Keep all samples
sampler.Iterations = 10000;       // Need many iterations

// With thinning  
sampler.ThinningInterval = 20;    // Keep every 20th sample
sampler.Iterations = 10000;       // Total iterations
// Effective samples = 10000 / 20 = 500 per chain
```

**Thinning trade-offs:**
- **Reduces autocorrelation** in final samples
- **Saves memory** for long runs
- **Doesn't improve efficiency** (better to run longer without thinning)
- **Rule of thumb:** Keep thinning interval ≤ autocorrelation length

## Multiple Chains

Running multiple chains helps assess convergence:

```cs
sampler.NumberOfChains = 4;       // Standard choice

// Access chain-specific information
int[] samplesPerChain = sampler.SampleCount;

Console.WriteLine("Samples per chain:");
for (int i = 0; i < samplesPerChain.Length; i++)
{
    Console.WriteLine($"  Chain {i + 1}: {samplesPerChain[i]} samples");
}
```

**Benefits of multiple chains:**
1. Assess convergence via R-hat statistic
2. Detect multimodal posteriors
3. Parallelize computation
4. More robust inference

## Warmup (Burn-in)

Warmup iterations are discarded to allow chains to reach stationarity:

```cs
sampler.InitialIterations = 50;      // Quick initialization
sampler.WarmupIterations = 2000;     // Burn-in / adaptation
sampler.Iterations = 5000;           // Kept samples

// Total iterations = Initial + Warmup + Main
// Only main iterations are kept
```

**Warmup guidelines:**
- **RWMH**: 2000-5000 iterations
- **ARWMH**: 2000-3000 (adapts during warmup)
- **DEMCz**: 1500-3000 (converges faster)
- **HMC**: 1000-2000 (efficient exploration)
- **Rule of thumb**: Warmup ≥ 50% of main iterations

## Best Practices

### 1. Always Check Convergence

```cs
// Visual inspection of traces
// Check R-hat < 1.1 for all parameters
// Effective sample size > 100 per parameter
```

### 2. Run Multiple Chains

```cs
// Minimum 4 chains for convergence assessment
sampler.NumberOfChains = 4;
```

### 3. Start with Enough Iterations

```cs
// Conservative starting point
sampler.WarmupIterations = 2000;
sampler.Iterations = 5000;
// Can adjust based on convergence diagnostics
```

### 4. Choose Appropriate Sampler

```cs
// Low dimensions (< 5): RWMH or ARWMH
// Medium (5-20): ARWMH (default choice)
// High (20+): DEMCz or DEMCzs
// Smooth posteriors: HMC
// Conjugate models: Gibbs
```

### 5. Informative Priors

```cs
// Use reasonable prior distributions
// Too vague: var prior = new Uniform(-1e10, 1e10);  // Bad!
// Better: var prior = new Normal(expectedValue, reasonableSD);
```

### 6. Validate Log-Likelihood

```cs
// Test log-likelihood function before sampling
double[] testParams = { 2.0, 1.5, 0.5 };
double logLik = logLikelihoodFunction(testParams);

Console.WriteLine($"Test log-likelihood: {logLik}");
if (double.IsNaN(logLik) || double.IsInfinity(logLik))
{
    Console.WriteLine("Warning: Check log-likelihood function!");
}
```

## Comparing Samplers

```cs
var samplers = new (string Name, MCMCSampler Sampler)[]
{
    ("RWMH", new RWMH(priors, logLik)),
    ("ARWMH", new ARWMH(priors, logLik)),
    ("DEMCz", new DEMCz(priors, logLik))
};

foreach (var (name, sampler) in samplers)
{
    sampler.WarmupIterations = 2000;
    sampler.Iterations = 5000;
    
    var watch = System.Diagnostics.Stopwatch.StartNew();
    sampler.Sample();
    watch.Stop();
    
    Console.WriteLine($"{name}: {sampler.ParameterSets.Length} samples in {watch.ElapsedMilliseconds}ms");
}
```

---

## References

<a id="1">[1]</a> Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

<a id="2">[2]</a> Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of state calculations by fast computing machines. *The Journal of Chemical Physics*, 21(6), 1087-1092.

<a id="3">[3]</a> Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm. *Bernoulli*, 7(2), 223-242.

<a id="4">[4]</a> ter Braak, C. J., & Vrugt, J. A. (2008). Differential evolution Markov chain with snooker updater and fewer chains. *Statistics and Computing*, 18(4), 435-446.

<a id="5">[5]</a> Neal, R. M. (2011). MCMC using Hamiltonian dynamics. *Handbook of Markov Chain Monte Carlo*, 2(11), 2.

---

[← Previous: Goodness-of-Fit](../statistics/goodness-of-fit.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](convergence-diagnostics.md)
