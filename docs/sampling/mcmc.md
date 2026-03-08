# MCMC Sampling

[← Previous: Random Generation](random-generation.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](convergence-diagnostics.md)

Markov Chain Monte Carlo (MCMC) methods sample from complex posterior distributions that are difficult to sample directly. The ***Numerics*** library provides multiple MCMC samplers for Bayesian inference, uncertainty quantification, and parameter estimation with full posterior distributions [[1]](#1).

## Available MCMC Samplers

| Sampler | Full Name | Best For | Key Features |
|---------|-----------|----------|--------------|
| **RWMH** | Random Walk Metropolis-Hastings | General purpose, small dimensions | Simple, robust baseline |
| **ARWMH** | Adaptive Random Walk M-H | Medium dimensions (2-20) | Self-tuning proposal |
| **DEMCz** | Differential Evolution MCMC | High dimensions, multimodal | Population-based, efficient |
| **DEMCzs** | DE-MCMC with snooker update | Very high dimensions | Enhanced DE-MCMC |
| **HMC** | Hamiltonian Monte Carlo | Smooth posteriors | Uses gradient information |
| **NUTS** | No-U-Turn Sampler | General smooth posteriors | Auto-tuning HMC |
| **Gibbs** | Gibbs Sampler | Conditional distributions available | No rejections |

## MCMC Fundamentals

Before using any specific sampler, it helps to understand the core theory that underpins all MCMC methods.

### Target Distribution

The goal of MCMC is to draw samples from a posterior distribution that is known only up to a normalizing constant. By Bayes' theorem, the posterior is proportional to the product of the prior and the likelihood:

```math
\pi(\theta \mid y) \propto \pi(\theta) \cdot L(y \mid \theta)
```

Taking the logarithm, which is how the Numerics library works internally:

```math
\log \pi(\theta \mid y) = \log \pi(\theta) + \log L(y \mid \theta) + \text{const}
```

The `LogLikelihoodFunction` delegate in Numerics should return the sum $\log \pi(\theta) + \log L(y \mid \theta)$, i.e., the unnormalized log-posterior.

### Detailed Balance

A Markov chain with transition kernel $T(\theta \to \theta')$ satisfies **detailed balance** with respect to $\pi$ if:

```math
\pi(\theta) \, T(\theta \to \theta') = \pi(\theta') \, T(\theta' \to \theta)
```

This reversibility condition guarantees that $\pi$ is a stationary distribution of the chain. All Metropolis-Hastings-based samplers in Numerics (RWMH, ARWMH, DEMCz, HMC, NUTS) enforce detailed balance through the accept/reject step, while Gibbs sampling satisfies it by construction when sampling from exact conditional distributions.

### Ergodicity

For the chain to converge to the target distribution regardless of its starting point, it must be **ergodic** -- meaning it is irreducible (can reach any region of positive probability) and aperiodic (does not cycle deterministically). In practice, ergodicity requires that the proposal distribution has support that covers the full posterior. The ARWMH sampler explicitly ensures ergodicity by mixing in a small identity-covariance proposal with probability $\beta$.

### Mixing

Mixing describes how quickly the chain "forgets" its starting point and produces effectively independent samples. Good mixing means low autocorrelation between successive states. Gradient-based samplers (HMC, NUTS) typically mix much faster than random-walk samplers because they follow the geometry of the posterior rather than exploring by diffusion. The **effective sample size (ESS)** quantifies mixing: a well-mixing chain has ESS close to the actual number of post-warmup samples.

### Acceptance Rate Guidelines

The acceptance rate is the fraction of proposed moves that are accepted. Acceptance rates that are too high indicate overly timid proposals (small steps), while rates that are too low mean the proposals are too aggressive. The following table summarizes established optimal acceptance rates from the literature:

| Sampler | Optimal Acceptance Rate | Reference |
|---------|------------------------|-----------|
| RWMH ($d = 1$) | ~44% | Roberts et al. (1997) [[7]](#7) |
| RWMH ($d \to \infty$) | ~23.4% | Roberts et al. (1997) [[7]](#7) |
| HMC | ~65% | Neal (2011) [[5]](#5) |
| NUTS | ~80% (target $\delta$) | Hoffman & Gelman (2014) [[6]](#6) |
| DEMCz | Varies by dimension | ter Braak & Vrugt (2008) [[4]](#4) |

These are rules of thumb. In practice, monitor the acceptance rates reported by `MCMCResults.AcceptanceRates` and adjust sampler settings if rates fall outside the expected range.

## Common MCMC Interface

All samplers inherit from `MCMCSampler` base class with common properties:

```cs
// Configuration
int PRNGSeed               // Random seed (default: 12345)
int InitialIterations      // Initialization phase (base default: 10; most samplers override, e.g. RWMH/ARWMH/NUTS use 100 * NumberOfParameters)
int WarmupIterations       // Burn-in period (default: 1750)
int Iterations             // Main sampling (default: 3500)
int NumberOfChains         // Parallel chains (default: 4)
int ThinningInterval       // Keep every nth sample (default: 20)

// Inputs
List<IUnivariateDistribution> PriorDistributions
LogLikelihood LogLikelihoodFunction

// Outputs (after sampling)
List<ParameterSet>[] Output       // Posterior samples (one list per chain)
List<ParameterSet>[] MarkovChains // Raw MCMC chains
ParameterSet MAP                  // Maximum a posteriori estimate
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

The simplest and most robust MCMC algorithm [[2]](#2).

### Mathematical Foundation

At each iteration, the Metropolis-Hastings algorithm proposes a new state $\theta^*$ from a proposal distribution $q(\theta^* \mid \theta)$ and accepts it with probability:

```math
\alpha(\theta^* \mid \theta) = \min\left(1, \frac{\pi(\theta^*) \cdot q(\theta \mid \theta^*)}{\pi(\theta) \cdot q(\theta^* \mid \theta)}\right)
```

In the RWMH implementation (`RWMH.cs`), the proposal is a symmetric multivariate normal centered at the current state:

```math
\theta^* \sim \mathcal{N}(\theta, \Sigma)
```

Because the proposal is symmetric -- $q(\theta^* \mid \theta) = q(\theta \mid \theta^*)$ -- the proposal densities cancel in the acceptance ratio, simplifying to:

```math
\alpha = \min\left(1, \frac{\pi(\theta^*)}{\pi(\theta)}\right)
```

In log space, which is how the source code computes it:

```math
\log \alpha = \log \pi(\theta^*) - \log \pi(\theta)
```

The proposal is accepted if $\log U \leq \log \alpha$, where $U \sim \text{Uniform}(0,1)$. If any proposed parameter falls outside its prior bounds, the proposal is immediately rejected without evaluating the log-likelihood.

The proposal covariance matrix $\Sigma$ must be provided by the user. A common starting point is a scaled identity matrix, but better performance is achieved when $\Sigma$ approximates the shape of the posterior (e.g., from a preliminary optimization).

```cs
using Numerics.Sampling.MCMC;
using Numerics.Sampling;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Data.Statistics;
using System.Linq;

// Create proposal covariance matrix (identity scaled for initial exploration)
int nParams = priors.Count;
var proposalSigma = Matrix.Identity(nParams) * 0.1;

// Create sampler
var rwmh = new RWMH(priors, logLikelihood, proposalSigma);

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
var samples = rwmh.Output.SelectMany(chain => chain).ToList();
Console.WriteLine($"Generated {samples.Count} samples");
Console.WriteLine($"Samples per chain: {string.Join(", ", rwmh.SampleCount)}");

// Posterior statistics
for (int i = 0; i < priors.Count; i++)
{
    var values = samples.Select(s => s.Values[i]).ToArray();
    double mean = values.Average();
    double std = Statistics.StandardDeviation(values);
    double q025 = Statistics.Percentile(values, 0.025);
    double q975 = Statistics.Percentile(values, 0.975);
    
    Console.WriteLine($"θ{i}: {mean:F3} ± {std:F3}, 95% CI: [{q025:F3}, {q975:F3}]");
}
```

**When to use RWMH:**
- General purpose baseline
- Low-dimensional problems (< 10 parameters)
- When simplicity and robustness are priorities
- As a reference for comparing other samplers

## Adaptive Random Walk M-H (ARWMH)

ARWMH automatically tunes the proposal distribution during warmup [[3]](#3).

### Mathematical Foundation

The key idea behind adaptive MCMC is to learn the proposal covariance from the chain's own history, eliminating the need for manual tuning [[8]](#8). The optimal scaling theory of Roberts and Rosenthal (2001) shows that for Gaussian targets in $d$ dimensions, the optimal proposal covariance is:

```math
\Sigma_{\text{opt}} = \frac{2.38^2}{d} \, \Sigma_{\text{target}}
```

where $\Sigma_{\text{target}}$ is the covariance of the target distribution. Since the target covariance is unknown, ARWMH estimates it online from the chain history.

In the Numerics implementation (`ARWMH.cs`), the proposal at iteration $t$ uses a mixture:

```math
\theta^* \sim \begin{cases} \mathcal{N}\!\left(\theta, \, \frac{0.1^2}{d} \, I_d\right) & \text{with probability } \beta \\ \mathcal{N}\!\left(\theta, \, \frac{2.38^2}{d} \, \hat{\Sigma}_t\right) & \text{with probability } 1-\beta \end{cases}
```

where:

- $d$ is the number of parameters (`NumberOfParameters`)
- $\beta = 0.05$ by default (the `Beta` property)
- $\hat{\Sigma}_t$ is the empirical covariance matrix computed as a running covariance of accepted samples (and current states after warmup)
- $I_d$ is the $d$-dimensional identity matrix
- The scale factor $s = 2.38^2/d$ is the `Scale` property

The small identity component (used with probability $\beta$, and also for the first $100 \times d$ samples) ensures **ergodicity**: even if the adaptive covariance estimate is poor, the chain can still reach any region of the parameter space.

The acceptance criterion is identical to RWMH -- since both proposal components are symmetric multivariate normals centered at $\theta$, the Hastings ratio simplifies to the posterior ratio.

```cs
var arwmh = new ARWMH(priors, logLikelihood);

// Configuration
arwmh.PRNGSeed = 12345;
arwmh.WarmupIterations = 2000;    // Adaptation happens here
arwmh.Iterations = 5000;
arwmh.NumberOfChains = 4;

Console.WriteLine("Running Adaptive RWMH sampler...");
arwmh.Sample();

var samples = arwmh.Output.SelectMany(chain => chain).ToList();
Console.WriteLine($"Generated {samples.Count} samples");

// ARWMH adapts the proposal covariance matrix during sampling (it does not explicitly target a specific acceptance rate)
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

Population-based sampler using differential evolution [[4]](#4).

### Mathematical Foundation

DEMCz combines differential evolution (DE) with MCMC by using a population of chains and a history of past states to generate proposals. The mutation formula from the source code (`DEMCz.cs`) is:

```math
\theta^*_i = \theta_i + \gamma \, (z_{R_1} - z_{R_2}) + e
```

where:

- $\gamma = 2.38 / \sqrt{2d}$ is the default jump rate (`Jump` property), with $d$ the number of parameters
- $z_{R_1}$ and $z_{R_2}$ are two randomly selected states from the **population matrix** (a memory of past states from all chains)
- $e \sim \mathcal{N}(0, b^2)$ is a small noise perturbation with default $b = 10^{-3}$ (`Noise` property)
- $R_1$ and $R_2$ are drawn uniformly without replacement from $\{1, 2, \ldots, M\}$, where $M$ is the current size of the population matrix

The proposal is accepted using the standard Metropolis ratio in log space:

```math
\log \alpha = \log \pi(\theta^*) - \log \pi(\theta)
```

To enable **mode-jumping** in multimodal posteriors, the jump rate $\gamma$ is set to $1.0$ with probability equal to `JumpThreshold` (default 0.1). When $\gamma = 1$, the proposal jumps the full difference between two past states, which can bridge gaps between separated modes.

The key insight of DEMCz is that the population matrix $Z$ serves as a memory of past states from **all** chains, providing a rich set of difference vectors for generating proposals. This eliminates the need for a manually specified proposal covariance matrix -- the population automatically learns the scale and orientation of the posterior.

```cs
var demcz = new DEMCz(priors, logLikelihood);

// Configuration
demcz.PRNGSeed = 12345;
demcz.NumberOfChains = 10;        // More chains for population diversity
demcz.WarmupIterations = 2000;
demcz.Iterations = 5000;

Console.WriteLine("Running DE-MCMC sampler...");
demcz.Sample();

var samples = demcz.Output.SelectMany(chain => chain).ToList();
Console.WriteLine($"Generated {samples.Count} samples from {demcz.NumberOfChains} chains");

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

Uses gradient information for efficient sampling [[5]](#5) [[9]](#9).

### Mathematical Foundation

HMC augments the parameter space with auxiliary momentum variables $\phi$ and simulates Hamiltonian dynamics to generate distant, high-quality proposals. The Hamiltonian is defined as:

```math
H(\theta, \phi) = U(\theta) + K(\phi) = -\log \pi(\theta) + \frac{1}{2} \phi^T M^{-1} \phi
```

where $U(\theta) = -\log \pi(\theta)$ is the **potential energy** (negative log-posterior) and $K(\phi)$ is the **kinetic energy**. The mass matrix $M$ is diagonal in the Numerics implementation (the `Mass` vector property).

Hamilton's equations of motion govern the dynamics:

```math
\frac{d\theta}{dt} = \frac{\partial H}{\partial \phi} = M^{-1}\phi, \qquad \frac{d\phi}{dt} = -\frac{\partial H}{\partial \theta} = \nabla \log \pi(\theta)
```

These continuous dynamics are approximated using the **leapfrog integrator**, as implemented in `HMC.cs`:

1. Half-step momentum: $\phi \leftarrow \phi + \frac{\varepsilon}{2} \nabla \log \pi(\theta)$
2. Full-step position: $\theta \leftarrow \theta + \varepsilon \, M^{-1} \phi$
3. Half-step momentum: $\phi \leftarrow \phi + \frac{\varepsilon}{2} \nabla \log \pi(\theta)$

Steps 2-3 are repeated for $L$ leapfrog steps. After the trajectory, a Metropolis correction accounts for numerical integration error:

```math
\alpha = \min\left(1, \exp\!\left(-H(\theta^*, \phi^*) + H(\theta, \phi)\right)\right)
```

In log space, the source code computes this as:

```math
\log \alpha = \left[\log \pi(\theta^*) - \frac{1}{2}\phi^{*T} M^{-1} \phi^*\right] - \left[\log \pi(\theta) - \frac{1}{2}\phi^T M^{-1} \phi\right]
```

**Implementation details from `HMC.cs`:**

- The step size $\varepsilon$ is **jittered**: each iteration draws $\varepsilon \sim \text{Uniform}(0, \, 2\varepsilon_0)$ where $\varepsilon_0$ is the `StepSize` property. This avoids resonant trajectories.
- The number of leapfrog steps $L$ is **jittered**: each iteration draws $L \sim \text{UniformDiscrete}(1, \, 2L_0)$ where $L_0$ is the `Steps` property.
- The mass vector $M$ is diagonal (default: identity). Users can set it via the `mass` constructor parameter.
- If no gradient function is provided, numerical finite differences are used via `NumericalDerivative.Gradient`, with probes clamped to prior bounds.

```cs
var hmc = new HMC(priors, logLikelihood);

// HMC-specific settings
hmc.NumberOfChains = 4;
hmc.WarmupIterations = 1000;      // HMC converges faster
hmc.Iterations = 2000;

Console.WriteLine("Running Hamiltonian Monte Carlo...");
hmc.Sample();

// HMC produces high-quality samples with lower autocorrelation
var samples = hmc.Output.SelectMany(chain => chain).ToList();
Console.WriteLine($"Generated {samples.Count} high-quality samples");
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

## No-U-Turn Sampler (NUTS)

NUTS automatically tunes the trajectory length that HMC requires as a manual setting, making it the recommended gradient-based sampler for most problems [[6]](#6).

### Mathematical Foundation

NUTS eliminates HMC's most sensitive tuning parameter -- the number of leapfrog steps $L$ -- by automatically determining when to stop the trajectory. It builds a balanced binary tree of leapfrog states using the following procedure:

1. Sample momentum $\phi \sim \mathcal{N}(0, M)$
2. Recursively double the trajectory in a randomly chosen direction (forward or backward)
3. Stop when a **U-turn** is detected or the maximum tree depth is reached

The **U-turn criterion** (from `NUTS.cs`) checks whether the trajectory endpoints are moving apart or starting to return:

```math
(\theta^+ - \theta^-) \cdot \phi^- < 0 \quad \text{or} \quad (\theta^+ - \theta^-) \cdot \phi^+ < 0
```

where $\theta^+$ and $\theta^-$ are the forward and backward endpoints of the trajectory, and $\phi^+$ and $\phi^-$ are their corresponding momenta. When either dot product is negative, the trajectory has started to curve back on itself.

**Candidate selection** uses multinomial sampling weighted by $\exp(-H)$:

```math
P(\text{select } \theta_j) = \frac{\exp(-H(\theta_j, \phi_j))}{\sum_k \exp(-H(\theta_k, \phi_k))}
```

This is implemented via log-sum-exp arithmetic for numerical stability.

**Step size adaptation** uses the dual averaging scheme from Hoffman and Gelman (2014), Algorithm 5. The running statistic $\bar{H}_m$ is updated at each adaptation step $m$:

```math
\bar{H}_m = \left(1 - \frac{1}{m + t_0}\right)\bar{H}_{m-1} + \frac{1}{m + t_0}(\delta - \alpha_m)
```

```math
\log \varepsilon_m = \mu - \frac{\sqrt{m}}{\gamma}\bar{H}_m
```

where:

- $\delta = 0.80$ is the target acceptance rate (`DELTA_TARGET`)
- $\gamma = 0.05$ is the adaptation regularization (`GAMMA`)
- $t_0 = 10$ prevents early instability (`T0`)
- $\mu = \log(10 \cdot \varepsilon_0)$ is the bias point, with $\varepsilon_0$ the initial step size
- $\alpha_m$ is the average Metropolis acceptance probability from the current tree

The smoothed step size uses an exponential moving average with decay exponent $\kappa = 0.75$:

```math
\log \bar{\varepsilon}_m = m^{-\kappa} \log \varepsilon_m + (1 - m^{-\kappa}) \log \bar{\varepsilon}_{m-1}
```

After warmup, the step size is fixed to $\exp(\log \bar{\varepsilon})$.

**Implementation details from `NUTS.cs`:**

- `MaxTreeDepth` defaults to 10, capping trajectories at $2^{10} = 1024$ leapfrog steps
- Divergence threshold: if $H - H_0 > 1000$, the trajectory is considered divergent and tree-building stops
- NUTS always accepts a candidate from the tree (acceptance is built into the multinomial weighting), so `AcceptCount` increments every iteration
- Step size adaptation occurs only during the warmup phase, with step sizes clamped to $[10^{-10}, \, 10^{5}]$

```cs
var nuts = new NUTS(priors, logLikelihood);

// NUTS-specific settings
nuts.NumberOfChains = 4;
nuts.WarmupIterations = 1000;       // Step size adapts during warmup
nuts.Iterations = 2000;

// Optional: set step size and max tree depth
// nuts = new NUTS(priors, logLikelihood, stepSize: 0.5, maxTreeDepth: 10);

Console.WriteLine("Running No-U-Turn Sampler...");
nuts.Sample();

// Analyze results using MCMCResults
var results = new MCMCResults(nuts);

Console.WriteLine($"Parameter 0 - Mean: {results.ParameterResults[0].SummaryStatistics.Mean:F3}");
Console.WriteLine($"Parameter 0 - Median: {results.ParameterResults[0].SummaryStatistics.Median:F3}");
```

**When to use NUTS:**
- Smooth, differentiable posteriors (same as HMC)
- When you don't want to tune the number of leapfrog steps
- Default gradient-based sampler for most applications
- Medium to high dimensions

**Advantages over HMC:**
- No manual tuning of trajectory length
- Adapts step size automatically via dual averaging
- Eliminates wasteful U-turns in the trajectory
- Generally more efficient per computation

**Key settings:**
- `stepSize`: Initial leapfrog step size (adapted during warmup)
- `maxTreeDepth`: Maximum binary tree depth (default 10, caps trajectory at 2^10 = 1024 steps)
- `gradientFunction`: Optional analytical gradient, provided as a constructor parameter only (not a settable property). If not provided, numerical finite differences are used.

## Gibbs Sampler

### Mathematical Foundation

The Gibbs sampler updates each parameter in turn by sampling from its **full conditional distribution** -- the distribution of one parameter given all the others and the data:

```math
\theta_j^{(t+1)} \sim p\!\left(\theta_j \mid \theta_1^{(t+1)}, \ldots, \theta_{j-1}^{(t+1)}, \theta_{j+1}^{(t)}, \ldots, \theta_d^{(t)}, y\right)
```

The key property is that **there is no rejection step** -- every proposed sample is accepted. This makes Gibbs sampling highly efficient when the conditional distributions are available in closed form (e.g., conjugate prior-likelihood pairs).

In the Numerics implementation (`Gibbs.cs`), the user provides a `Proposal` delegate with the signature `double[] Proposal(double[] parameters, Random prng)`. This delegate takes the current parameter vector and a pseudo-random number generator, and returns a new parameter vector sampled from the conditional distributions. The sampler always accepts the returned values (no Metropolis step).

Default settings differ from other samplers: the Gibbs sampler uses 1 chain, no thinning (`ThinningInterval = 1`), minimal warmup (`WarmupIterations = 1`), and 100,000 iterations with an output buffer of 10,000 samples.

Samples each parameter conditionally given others:

```cs
// Define a proposal function that samples each parameter from its conditional distribution.
// The Gibbs.Proposal delegate signature is: double[] Proposal(double[] parameters, Random prng)
Gibbs.Proposal proposalFunction = (parameters, prng) =>
{
    var proposed = (double[])parameters.Clone();
    for (int i = 0; i < priors.Count; i++)
    {
        // Sample each parameter from a small perturbation around the current value
        proposed[i] = parameters[i] + new Normal(0, 0.1).InverseCDF(prng.NextDouble());
        // Clamp to prior bounds
        proposed[i] = Math.Max(priors[i].Minimum, Math.Min(priors[i].Maximum, proposed[i]));
    }
    return proposed;
};

var gibbs = new Gibbs(priors, logLikelihood, proposalFunction);

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
using Numerics.Sampling;
using Numerics.Data.Statistics;
using System.Linq;

// Generate synthetic data
double trueIntercept = 2.0;
double trueSlope = 1.8;
double trueNoise = 0.5;

var random = new MersenneTwister(123);
int n = 20;
double[] x = Enumerable.Range(1, n).Select(i => (double)i).ToArray();
double[] yTrue = x.Select(xi => trueIntercept + trueSlope * xi).ToArray();
var noiseDist = new Normal(0, trueNoise);
double[] y = yTrue.Select(yi => yi + noiseDist.InverseCDF(random.NextDouble())).ToArray();

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
var samples = sampler.Output.SelectMany(chain => chain).ToList();
Console.WriteLine($"\nPosterior Summary ({samples.Count} samples):");
Console.WriteLine("Parameter | True  | Post Mean | Post SD | 95% Credible Interval");
Console.WriteLine("---------------------------------------------------------------");

string[] names = { "Intercept", "Slope", "Noise SD" };
double[] trueVals = { trueIntercept, trueSlope, trueNoise };

for (int i = 0; i < 3; i++)
{
    var vals = samples.Select(s => s.Values[i]).ToArray();
    double mean = vals.Average();
    double std = Statistics.StandardDeviation(vals);
    double lower = Statistics.Percentile(vals, 0.025);
    double upper = Statistics.Percentile(vals, 0.975);
    
    Console.WriteLine($"{names[i],-9} | {trueVals[i],5:F2} | {mean,9:F3} | {std,7:F3} | [{lower:F3}, {upper:F3}]");
}

// Posterior predictive
Console.WriteLine("\nPosterior Predictive at x=15:");
double xNew = 15;
var predictions = samples.Select(s => s.Values[0] + s.Values[1] * xNew).ToArray();
double predMean = predictions.Average();
double predSD = Statistics.StandardDeviation(predictions);
double predLower = Statistics.Percentile(predictions, 0.025);
double predUpper = Statistics.Percentile(predictions, 0.975);

Console.WriteLine($"E[y|x={xNew}] = {predMean:F2} ± {predSD:F2}");
Console.WriteLine($"95% Credible Interval: [{predLower:F2}, {predUpper:F2}]");
```

### Example 2: Bayesian Distribution Fitting with Real Streamflow Data

This example fits a Normal distribution to Tippecanoe River streamflow data using ARWMH, following the same methodology used in the library's unit tests. Results are validated against rstan with comparable MCMC settings.

**Data source:** Rao, A. R. & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press, Table 5.1.1.
See also: [`example-data/tippecanoe-river-streamflow.csv`](../example-data/tippecanoe-river-streamflow.csv)

```cs
using Numerics.Distributions;
using Numerics.Sampling.MCMC;
using Numerics.Data.Statistics;
using System.Linq;

// Tippecanoe River near Delphi, IN — 48 years of annual peak streamflow (cfs)
// Source: Rao & Hamed (2000), Table 5.1.1
double[] annualPeaks = {
    6290, 2700, 13100, 16900, 14600, 9600, 7740, 8490, 8130, 12000,
    17200, 15000, 12400, 6960, 6500, 5840, 10400, 18800, 21400, 22600,
    14200, 11000, 12800, 15700, 4740, 6950, 11800, 12100, 20600, 14600,
    14600, 8900, 10600, 14200, 14100, 14100, 12500, 7530, 13400, 17600,
    13400, 19200, 16900, 15500, 14500, 21900, 10400, 7460
};

Console.WriteLine("Bayesian Estimation of Normal Distribution Parameters");
Console.WriteLine("Tippecanoe River near Delphi, IN (n=48)");

// Prior distributions for Normal(μ, σ)
// Use weakly informative priors centered on data scale
var priors = new List<IUnivariateDistribution>
{
    new Uniform(0, 50000),     // μ: Uniform over reasonable flow range
    new Uniform(0, 50000)      // σ: Uniform positive
};

// Log-likelihood function
LogLikelihood logLik = (theta) =>
{
    double mu = theta[0];
    double sigma = theta[1];
    if (sigma <= 0) return double.NegativeInfinity;

    var model = new Normal(mu, sigma);
    return model.LogLikelihood(annualPeaks);
};

// Run ARWMH sampler
var sampler = new ARWMH(priors, logLik);
sampler.PRNGSeed = 12345;
sampler.WarmupIterations = 2000;
sampler.Iterations = 5000;
sampler.NumberOfChains = 4;
sampler.ThinningInterval = 10;

sampler.Sample();

// Analyze posterior
var samples = sampler.Output.SelectMany(chain => chain).ToList();
string[] paramNames = { "μ (mean)", "σ (std dev)" };

Console.WriteLine($"\nPosterior Summary ({samples.Count} samples):");
for (int i = 0; i < 2; i++)
{
    var vals = samples.Select(s => s.Values[i]).ToArray();
    Console.WriteLine($"  {paramNames[i]}:");
    Console.WriteLine($"    Mean:   {vals.Average():F1}");
    Console.WriteLine($"    SD:     {Statistics.StandardDeviation(vals):F1}");
    Console.WriteLine($"    95% CI: [{Statistics.Percentile(vals, 0.025):F1}, " +
                     $"{Statistics.Percentile(vals, 0.975):F1}]");
}

// Expected results (validated against rstan):
//   μ:  mean ≈ 12664, sd ≈ 707
//   σ:  mean ≈ 4844,  sd ≈ 519

// MAP estimate
Console.WriteLine($"\nMAP estimate: μ={sampler.MAP.Values[0]:F1}, σ={sampler.MAP.Values[1]:F1}");
```

### Post-Processing with MCMCResults

The `MCMCResults` class provides structured post-processing of MCMC output, including summary statistics, convergence diagnostics, kernel density estimates, and histograms for each parameter:

```cs
using Numerics.Sampling.MCMC;

// After running the sampler, create results object
var results = new MCMCResults(sampler, alpha: 0.1);

// Access structured parameter summaries
Console.WriteLine("Parameter Summary:");
Console.WriteLine($"{"Param",-8} {"Mean",10} {"Median",10} {"SD",10} {"5%",10} {"95%",10} {"R-hat",8} {"ESS",8}");

for (int i = 0; i < results.ParameterResults.Length; i++)
{
    var stats = results.ParameterResults[i].SummaryStatistics;
    Console.WriteLine($"{paramNames[i],-8} {stats.Mean,10:F2} {stats.Median,10:F2} " +
                      $"{stats.StandardDeviation,10:F2} {stats.LowerCI,10:F2} " +
                      $"{stats.UpperCI,10:F2} {stats.Rhat,8:F4} {stats.ESS,8:F0}");
}

// MAP and posterior mean estimates
Console.WriteLine($"\nMAP: [{string.Join(", ", results.MAP.Values.Select(v => v.ToString("F2")))}]");
Console.WriteLine($"Posterior Mean: [{string.Join(", ", results.PosteriorMean.Values.Select(v => v.ToString("F2")))}]");

// Access all posterior samples
Console.WriteLine($"\nTotal posterior samples: {results.Output.Count}");

// Access individual chain results
if (results.MarkovChains != null)
{
    for (int c = 0; c < results.MarkovChains.Length; c++)
    {
        Console.WriteLine($"  Chain {c}: {results.MarkovChains[c].Count} samples");
    }
}
```

**`MCMCResults` properties:**

| Property | Type | Description |
|----------|------|-------------|
| `ParameterResults` | `ParameterResults[]` | Per-parameter statistics, KDE, and histograms |
| `MAP` | `ParameterSet` | Maximum a posteriori estimate |
| `PosteriorMean` | `ParameterSet` | Mean of posterior samples |
| `Output` | `List<ParameterSet>` | All posterior samples (aggregated across chains) |
| `MarkovChains` | `List<ParameterSet>[]` | Per-chain sample lists |
| `AcceptanceRates` | `double[]` | Acceptance rate for each chain |
| `MeanLogLikelihood` | `List<double>` | Average log-likelihood per iteration |

**`ParameterStatistics` properties** (via `ParameterResults[i].SummaryStatistics`):

| Property | Type | Description |
|----------|------|-------------|
| `Mean` | `double` | Posterior mean |
| `Median` | `double` | Posterior median |
| `StandardDeviation` | `double` | Posterior standard deviation |
| `LowerCI` | `double` | Lower confidence interval (default 5th percentile) |
| `UpperCI` | `double` | Upper confidence interval (default 95th percentile) |
| `Rhat` | `double` | Gelman-Rubin convergence diagnostic |
| `ESS` | `double` | Effective sample size |
| `N` | `int` | Total sample count |

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
- **NUTS**: 1000-2000 (step size adapts during warmup)
- **Rule of thumb**: Warmup ≥ 50% of main iterations

## Algorithm Selection Guide

Use the following decision tree to choose the most appropriate sampler for your problem:

```
Do you have gradient information (or a smooth, differentiable log-posterior)?
├── YES: Use NUTS (auto-tunes trajectory length)
│   └── If NUTS is too slow per iteration: Try HMC with manually tuned ε, L
├── NO: How many parameters?
│   ├── d ≤ 5: ARWMH (simple, robust, self-tuning)
│   ├── 5 < d ≤ 20: ARWMH or DEMCz
│   │   └── Multimodal posterior? → DEMCz (population-based exploration)
│   └── d > 20: DEMCz or DEMCzs
│       └── Very high dimensions → DEMCzs (snooker update for better mixing)
└── Conjugate model with known conditional distributions?
    └── YES: Gibbs (no rejections, maximum efficiency)
```

**Key trade-offs:**

- **RWMH** is the simplest algorithm and a useful baseline, but requires manual proposal tuning. Use it when you want full control or need a reference sampler.
- **ARWMH** is the recommended default for most problems without gradients. It eliminates manual tuning and adapts to posterior correlations automatically.
- **DEMCz / DEMCzs** excel at high-dimensional and multimodal problems. They require more chains (minimum 3, typically 10+) but need no proposal covariance specification.
- **HMC / NUTS** provide the best mixing per sample for smooth posteriors, but require gradient computation (analytical or numerical). NUTS is preferred over HMC because it eliminates manual trajectory-length tuning.
- **Gibbs** is maximally efficient when exact conditional distributions are available, but its applicability is limited to conjugate or conditionally tractable models.

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
// Smooth posteriors: NUTS (preferred) or HMC
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
using Numerics.Sampling.MCMC;
using Numerics.Mathematics.LinearAlgebra;
using System.Linq;

// RWMH requires a proposal covariance matrix
var rwmhProposal = Matrix.Identity(priors.Count) * 0.1;

var samplers = new (string Name, MCMCSampler Sampler)[]
{
    ("RWMH", new RWMH(priors, logLik, rwmhProposal)),
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

    int totalSamples = sampler.Output.Sum(chain => chain.Count);
    Console.WriteLine($"{name}: {totalSamples} samples in {watch.ElapsedMilliseconds}ms");
}
```

## Self-Normalizing Importance Sampling (SNIS)

The `SNIS` sampler performs self-normalizing importance sampling rather than iterative MCMC. It draws independent samples from a proposal distribution and re-weights them by the likelihood. This is useful when the posterior is close to the prior or when a good proposal distribution is available.

```cs
using Numerics.Sampling.MCMC;
using Numerics.Distributions;

// Observed data
double[] observations = { 3.2, 4.8, 5.1, 6.3, 4.9, 5.5, 6.1, 3.8, 5.0, 4.7 };

// Define priors
var priors = new List<IUnivariateDistribution>
{
    new Normal(0, 10),
    new Uniform(0, 50)
};

// Log-likelihood function
// Note: use a lowercase name to avoid shadowing the LogLikelihood delegate type.
LogLikelihood computeLogLikelihood = (double[] theta) =>
{
    double mu = theta[0], sigma = theta[1];
    if (sigma <= 0) return double.NegativeInfinity;
    var model = new Normal(mu, sigma);
    return model.LogLikelihood(observations);
};

// Naive Monte Carlo (sample from priors)
var snis = new SNIS(priors, computeLogLikelihood);
snis.Iterations = 100000;
snis.Sample();

// With a Multivariate Normal proposal distribution (importance sampling)
var proposal = new MultivariateNormal(
    new double[] { 5, 10 },   // mean vector
    new double[,] { { 4, 0 }, { 0, 25 } }  // covariance matrix
);
var snisIS = new SNIS(priors, computeLogLikelihood, proposal);
snisIS.Iterations = 100000;
snisIS.Sample();
```

**Key differences from MCMC samplers:**
- Single chain only (`NumberOfChains` = 1)
- No warmup or thinning — all samples contribute
- Produces weighted samples (weights stored in `ParameterSet.Weight`)
- No autocorrelation — samples are independent
- Best for low-dimensional problems with informative proposals

---

## References

<a id="1">[1]</a> Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

<a id="2">[2]</a> Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of state calculations by fast computing machines. *The Journal of Chemical Physics*, 21(6), 1087-1092.

<a id="3">[3]</a> Haario, H., Saksman, E., & Tamminen, J. (2001). An adaptive Metropolis algorithm. *Bernoulli*, 7(2), 223-242.

<a id="4">[4]</a> ter Braak, C. J., & Vrugt, J. A. (2008). Differential evolution Markov chain with snooker updater and fewer chains. *Statistics and Computing*, 18(4), 435-446.

<a id="5">[5]</a> Neal, R. M. (2011). MCMC using Hamiltonian dynamics. *Handbook of Markov Chain Monte Carlo*, 2(11), 2.

<a id="6">[6]</a> Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn Sampler: Adaptively setting path lengths in Hamiltonian Monte Carlo. *Journal of Machine Learning Research*, 15(47), 1593-1623.

<a id="7">[7]</a> Roberts, G. O., Gelman, A., & Gilks, W. R. (1997). Weak convergence and optimal scaling of random walk Metropolis algorithms. *Annals of Applied Probability*, 7(1), 110-120.

<a id="8">[8]</a> Roberts, G. O., & Rosenthal, J. S. (2001). Optimal scaling for various Metropolis-Hastings algorithms. *Statistical Science*, 16(4), 351-367.

<a id="9">[9]</a> Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint arXiv:1701.02434*.

---

[← Previous: Random Generation](random-generation.md) | [Back to Index](../index.md) | [Next: Convergence Diagnostics →](convergence-diagnostics.md)
