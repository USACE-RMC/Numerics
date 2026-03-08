# MCMC Convergence Diagnostics

[← Previous: MCMC Methods](mcmc.md) | [Back to Index](../index.md)

Convergence diagnostics assess whether MCMC samplers have reached their stationary distribution and provide reliable samples from the posterior. The ***Numerics*** library provides essential diagnostic tools including Gelman-Rubin statistic and Effective Sample Size.

## Why Convergence Matters

MCMC samplers:
1. Start from arbitrary initial values
2. Explore parameter space stochastically
3. Eventually converge to target distribution
4. Must discard "burn-in" samples

**Key questions:**
- Have chains converged to stationary distribution?
- How many independent samples do we have?
- Is the warmup period sufficient?

### Theoretical Foundation

The theoretical justification for MCMC rests on the **ergodic theorem**: under regularity conditions (irreducibility, aperiodicity, positive recurrence), the time average of a function along the Markov chain converges to its expectation under the stationary distribution:

```math
\frac{1}{N}\sum_{t=1}^{N}f(\theta_t) \;\xrightarrow{\;a.s.\;}\; \mathbb{E}_\pi[f(\theta)] \quad \text{as } N \to \infty
```

where $\pi$ is the target (posterior) distribution. This guarantee is asymptotic -- for any finite $N$, we need diagnostics to assess whether we are close enough to this limit.

**Burn-in (warmup)** refers to the initial transient phase where the chain has not yet reached the stationary distribution. Samples from this phase are drawn from a distribution that depends on the arbitrary starting point, not from $\pi$. Discarding these samples is essential for valid inference.

**Mixing** describes how quickly the chain "forgets" its current state and explores the full support of $\pi$. A well-mixing chain has rapidly decaying autocorrelation -- the correlation between $\theta_t$ and $\theta_{t+k}$ diminishes quickly with lag $k$. Poorly mixing chains remain in local regions for long periods, producing highly correlated samples and requiring far more iterations to achieve reliable estimates.

## Gelman-Rubin Statistic (R̂)

The Gelman-Rubin diagnostic compares within-chain and between-chain variance [[1]](#1). Values near 1.0 indicate convergence.

### Computing R̂

```cs
using Numerics.Sampling.MCMC;
using Numerics.Mathematics.Optimization;

// Run sampler with multiple chains
var sampler = new ARWMH(priors, logLikelihood);
sampler.NumberOfChains = 4;
sampler.WarmupIterations = 2000;
sampler.Iterations = 5000;
sampler.Sample();

// Get chains from sampler
var chains = sampler.MarkovChains.ToList();

// Compute Gelman-Rubin for each parameter
int warmup = sampler.WarmupIterations;
double[] rHat = MCMCDiagnostics.GelmanRubin(chains, warmup);

Console.WriteLine("Gelman-Rubin Diagnostics (R̂):");
for (int i = 0; i < rHat.Length; i++)
{
    Console.WriteLine($"  Parameter {i}: R̂ = {rHat[i]:F4}");
    
    if (rHat[i] < 1.1)
        Console.WriteLine("    ✓ Converged");
    else if (rHat[i] < 1.2)
        Console.WriteLine("    ⚠ Marginal - run longer");
    else
        Console.WriteLine("    ✗ Not converged - investigate");
}
```

### Interpretation

| R̂ Value | Interpretation | Action |
|---------|---------------|--------|
| R̂ < 1.01 | Excellent convergence | Proceed |
| R̂ < 1.1 | Good convergence | Safe to use |
| 1.1 ≤ R̂ < 1.2 | Marginal | Run longer |
| R̂ ≥ 1.2 | Poor convergence | Investigate |

**Formula:**

```math
\hat{R} = \sqrt{\frac{\hat{V}}{W}}
```
Where $W$ is the mean within-chain variance, $B$ is the between-chain variance, and:
```math
\hat{V} = \frac{n-1}{n}W + \frac{1}{n}B
```

### Mathematical Derivation

Consider $m$ chains, each of length $n$, sampling a scalar parameter $\theta$. Let $\theta_{jt}$ denote the $t$-th sample from chain $j$.

**Step 1: Chain means and overall mean.** Compute the mean of each chain and the grand mean across all chains:

```math
\bar{\theta}_j = \frac{1}{n}\sum_{t=1}^{n}\theta_{jt}, \qquad \bar{\theta}_{..} = \frac{1}{m}\sum_{j=1}^{m}\bar{\theta}_j
```

**Step 2: Between-chain variance $B$.** This measures how much the chain means differ from each other. Large $B$ relative to the within-chain variability indicates that chains have not converged to the same distribution:

```math
B = \frac{n}{m-1}\sum_{j=1}^{m}\left(\bar{\theta}_j - \bar{\theta}_{..}\right)^2
```

The factor $n/(m-1)$ scales $B$ so that it estimates the variance of $\theta$ under stationarity (the scaling by $n$ converts from variance-of-means to variance-of-individual-draws).

**Step 3: Within-chain variance $W$.** This is the average of the individual chain variances. Each chain's variance $s_j^2$ uses the Bessel-corrected estimator:

```math
s_j^2 = \frac{1}{n-1}\sum_{t=1}^{n}\left(\theta_{jt} - \bar{\theta}_j\right)^2
```

```math
W = \frac{1}{m}\sum_{j=1}^{m}s_j^2
```

**Step 4: Pooled variance estimate.** The pooled estimate combines within-chain and between-chain information:

```math
\hat{V} = \frac{n-1}{n}\,W + \frac{1}{n}\,B
```

This is a weighted average that has a key property: $\hat{V}$ **overestimates** the true target variance when chains have not converged, because $B$ captures the additional spread from chains being in different regions. Meanwhile, $W$ **underestimates** the true target variance because each finite chain has only explored a portion of the full parameter space. At convergence, the between-chain contribution vanishes ($B/n \to 0$ relative to $W$), and $\hat{V} \to W$.

**Step 5: The diagnostic ratio.** The potential scale reduction factor is:

```math
\hat{R} = \sqrt{\frac{\hat{V}}{W}}
```

Since $\hat{V} \geq W$ in general, we have $\hat{R} \geq 1$. At perfect convergence $\hat{R} = 1$; values substantially above 1 indicate that the chains have not mixed and further sampling is needed.

**Split-$\hat{R}$.** Modern practice [[4]](#4) recommends splitting each chain in half before computing $\hat{R}$, which doubles the number of chains from $m$ to $2m$. This helps detect non-stationarity *within* individual chains -- for example, a chain that drifted during the first half but settled during the second half. The ***Numerics*** implementation does not perform split-$\hat{R}$ automatically; to use this approach, split each chain manually before passing them to `GelmanRubin()`.

### Common Causes of High R̂

1. **Insufficient warmup** - Chains haven't reached stationarity
2. **Poor mixing** - Chains explore slowly
3. **Multimodal posterior** - Chains stuck in different modes
4. **Bad initialization** - Starting values too extreme
5. **Wrong sampler** - Algorithm not suited to problem

## Effective Sample Size (ESS)

ESS quantifies number of independent samples, accounting for autocorrelation [[2]](#2).

### Computing ESS

```cs
// Extract samples for parameter 0 from the first chain
double[] samples = sampler.MarkovChains[0]
    .Select(ps => ps.Values[0])
    .ToArray();

double ess = MCMCDiagnostics.EffectiveSampleSize(samples);

Console.WriteLine($"Effective Sample Size: {ess:F0}");
Console.WriteLine($"Actual samples: {samples.Length}");
Console.WriteLine($"Efficiency: {ess / samples.Length:P1}");

// Rule of thumb: ESS > 100 per parameter
if (ess < 100)
    Console.WriteLine("⚠ Warning: Low ESS - run longer or thin more");
```

### ESS Across All Parameters

```cs
// Compute ESS for all parameters across chains
double[] essValues = MCMCDiagnostics.EffectiveSampleSize(chains, out double[][,] avgACF);

Console.WriteLine("Effective Sample Size by Parameter:");
for (int i = 0; i < essValues.Length; i++)
{
    Console.WriteLine($"  θ{i}: ESS = {essValues[i]:F0}");
}

// Check minimum ESS
double minESS = essValues.Min();
Console.WriteLine($"\nMinimum ESS: {minESS:F0}");

if (minESS > 400)
    Console.WriteLine("✓ Excellent: ESS > 400");
else if (minESS > 100)
    Console.WriteLine("✓ Good: ESS > 100");
else
    Console.WriteLine("✗ Poor: ESS < 100 - need more samples");
```

### Interpretation

**Guidelines:**
- **ESS > 400**: Excellent - precise estimates
- **ESS > 100**: Good - adequate for most purposes
- **ESS < 100**: Poor - increase iterations or improve mixing

**Autocorrelation impact:**
- High autocorrelation → Low ESS → Need more iterations
- Good mixing → High ESS → Efficient sampling

### ESS Formula

```math
\text{ESS} = \frac{N}{1 + 2\sum_{k=1}^{K} \rho_k}
```
where $N$ is the number of samples, $\rho_k$ is the autocorrelation at lag $k$, and the sum is truncated when $\rho_k$ becomes negligible.

### Mathematical Derivation

The ESS formula arises from analyzing the variance of the sample mean of a correlated sequence. For a stationary process $\lbrace\theta_1, \theta_2, \ldots, \theta_N\rbrace$ with marginal variance $\sigma^2$ and autocorrelation function $\rho_k = \text{Corr}(\theta_t, \theta_{t+k})$, the variance of the sample mean $\bar{\theta} = \frac{1}{N}\sum_{t=1}^{N}\theta_t$ is:

```math
\text{Var}(\bar{\theta}) = \frac{\sigma^2}{N}\left(1 + 2\sum_{k=1}^{N-1}\left(1 - \frac{k}{N}\right)\rho_k\right)
```

For large $N$, the $(1 - k/N)$ correction becomes negligible for the lags that matter, and the expression simplifies to:

```math
\text{Var}(\bar{\theta}) \approx \frac{\sigma^2}{N}\left(1 + 2\sum_{k=1}^{\infty}\rho_k\right)
```

If the samples were independent, we would have $\rho_k = 0$ for all $k \geq 1$, giving $\text{Var}(\bar{\theta}) = \sigma^2/N$. The effective sample size is defined as the number of *independent* samples that would give the same variance for the sample mean:

```math
\frac{\sigma^2}{\text{ESS}} = \frac{\sigma^2}{N}\left(1 + 2\sum_{k=1}^{\infty}\rho_k\right)
```

Solving for ESS:

```math
\text{ESS} = \frac{N}{1 + 2\sum_{k=1}^{\infty}\rho_k}
```

The quantity $\tau = 1 + 2\sum_{k=1}^{\infty}\rho_k$ is called the **integrated autocorrelation time**. It represents how many MCMC iterations correspond to one independent draw: $\text{ESS} = N/\tau$.

**Truncation strategy.** In practice, the infinite sum must be truncated. The ***Numerics*** implementation uses a simple truncation rule: the sum is cut off at the first lag $k$ where $\rho_k < 0$. This works because for a well-behaved MCMC chain, the autocorrelation function decays monotonically toward zero and oscillations below zero represent noise rather than genuine correlation.

Geyer (1992) [[3]](#3) proposed a more robust alternative called the **initial positive sequence estimator**, which sums consecutive *pairs* of autocorrelations $(\rho_{2k} + \rho_{2k+1})$ and stops when a pair sum becomes negative. This approach is theoretically guaranteed to produce a non-negative variance estimate. The ***Numerics*** implementation uses the simpler first-negative truncation, which is adequate for chains with good mixing behavior.

**Multi-chain ESS.** When $M$ chains of length $N$ are available, the implementation computes the autocorrelation sum $\rho_m$ for each chain $m$ separately, then averages across chains:

```math
\bar{\rho} = \frac{1}{M}\sum_{m=1}^{M}\rho_m \qquad \text{where} \quad \rho_m = \sum_{k=1}^{K_m}\hat{\rho}_k^{(m)}
```

Here $K_m$ is the truncation point for chain $m$ (the first lag at which the autocorrelation is negative). The total effective sample size is then:

```math
\text{ESS} = \frac{N \cdot M}{1 + 2\bar{\rho}}
```

This is capped at $N \cdot M$ (the total number of samples) since the effective sample size cannot exceed the actual number of draws.

### ESS Requirements

The minimum ESS needed depends on what posterior summary you are estimating:

| Inference Goal | Minimum ESS | Rationale |
|---------------|------------|-----------|
| Posterior mean | 100 | MCSE is approximately 10% of posterior SD |
| Posterior standard deviation | 200 | Variance estimation requires more samples than mean estimation |
| 95% credible interval | 400 | Quantile estimation demands greater precision in the tails |
| Tail probabilities (e.g., $P(\theta > c)$) | 1000+ | Extreme quantiles are estimated from sparse tail samples |

These thresholds are guidelines, not strict rules. For life-safety applications, err on the side of larger ESS.

### ESS per Second

When comparing MCMC samplers, ESS alone is not sufficient -- the computational cost per iteration matters. The **ESS per second** metric accounts for this:

```math
\text{ESS/s} = \frac{\text{ESS}}{T_{\text{compute}}}
```

where $T_{\text{compute}}$ is the total wall-clock time for sampling. This metric is critical for sampler selection: Hamiltonian Monte Carlo (HMC) typically achieves higher ESS per iteration than Random Walk Metropolis-Hastings (RWMH) because HMC's proposals are guided by gradient information and produce less correlated samples. However, each HMC iteration requires evaluating the gradient of the log-posterior (and often multiple leapfrog steps), making it more expensive per iteration. The optimal sampler is the one that maximizes ESS/s for the problem at hand.

## Monte Carlo Standard Error (MCSE)

The Monte Carlo Standard Error quantifies the precision of posterior estimates due to finite sampling [[5]](#5). While the posterior standard deviation describes uncertainty about the parameter, the MCSE describes uncertainty about the *estimate itself*.

### Formula

For a posterior mean estimate $\bar{\theta}$ computed from MCMC output with posterior standard deviation $\text{SD}(\theta)$ and effective sample size ESS:

```math
\text{MCSE} = \frac{\text{SD}(\theta)}{\sqrt{\text{ESS}}}
```

This is the standard error of the Monte Carlo estimate of $\mathbb{E}_\pi[\theta]$. It tells you how much the posterior mean would vary if you repeated the entire MCMC run.

### Interpretation

A useful rule of thumb is that the MCSE should be less than 5% of the posterior standard deviation:

```math
\frac{\text{MCSE}}{\text{SD}(\theta)} = \frac{1}{\sqrt{\text{ESS}}} < 0.05 \quad \Longrightarrow \quad \text{ESS} > 400
```

This means your Monte Carlo error is small relative to the inherent posterior uncertainty. For life-safety applications, a more stringent threshold (e.g., MCSE/SD < 0.02, requiring ESS > 2500) may be appropriate.

### Computing MCSE

```cs
// Compute MCSE from ESS and posterior standard deviation
double[] values = samples.Select(s => s.Values[0]).ToArray();
double sd = Statistics.StandardDeviation(values);
double mcse = sd / Math.Sqrt(ess[0]);

Console.WriteLine($"Posterior SD:    {sd:F6}");
Console.WriteLine($"ESS:             {ess[0]:F0}");
Console.WriteLine($"MCSE:            {mcse:F6}");
Console.WriteLine($"MCSE/SD ratio:   {mcse / sd:P2}");

if (mcse / sd > 0.05)
    Console.WriteLine("⚠ MCSE is large relative to posterior SD - increase sampling");
else
    Console.WriteLine("✓ MCSE is acceptably small");
```

## Autocorrelation

### Visualizing Autocorrelation

```cs
// Compute autocorrelation function
// (avgACF from ESS calculation above)

Console.WriteLine("Autocorrelation at selected lags:");
Console.WriteLine("Lag  | ACF");
Console.WriteLine("-----|------");

for (int lag = 0; lag <= 20; lag += 5)
{
    // Access from avgACF[parameter][lag, column] where column 1 has ACF values
    double acf = avgACF[0][lag, 1];  // Parameter 0, average ACF at this lag
    Console.WriteLine($"{lag,4} | {acf,5:F3}");
}

// Ideal: ACF drops quickly to zero
// Problem: ACF remains high (slow decorrelation)
```

### Autocorrelation Guidelines

| Lag-1 ACF | Interpretation | ESS Impact |
|-----------|---------------|------------|
| < 0.1 | Excellent mixing | ESS ≈ N |
| 0.1-0.3 | Good mixing | ESS ≈ 0.5N |
| 0.3-0.6 | Moderate | ESS ≈ 0.2N |
| > 0.6 | Poor mixing | ESS << 0.1N |

## Minimum Sample Size

The Raftery-Lewis method provides a theoretical lower bound on the number of MCMC samples needed to estimate a particular posterior quantile with specified precision. Given a quantile of interest $q$, tolerance $r$, and desired probability $s$, the minimum sample size is:

```math
N_{\min} = \frac{q(1 - q)\left[z_{(1+s)/2}\right]^2}{r^2}
```

where $z_{(1+s)/2}$ is the standard normal quantile (inverse CDF) evaluated at $(1+s)/2$. This formula arises from the normal approximation to the binomial: the indicator $I(\theta \leq \theta_q)$ has variance $q(1-q)$ under the posterior, and $z_{(1+s)/2}$ ensures that the estimate falls within tolerance $r$ of the true quantile with probability $s$. The ***Numerics*** implementation rounds the result to the nearest hundred.

For example, to estimate the 99th percentile ($q = 0.99$) within $\pm 0.01$ ($r = 0.01$) with 95% confidence ($s = 0.95$):

```math
N_{\min} = \frac{0.99 \times 0.01 \times (1.96)^2}{0.01^2} = \frac{0.0099 \times 3.8416}{0.0001} \approx 380
```

Note that this is a lower bound assuming independent samples. With autocorrelated MCMC output, the actual number of iterations needed is $N_{\min} \times \tau$, where $\tau$ is the integrated autocorrelation time ($\tau = N/\text{ESS}$).

Determine required sample size for desired precision:

```cs
// For quantile estimation
double quantile = 0.99;       // 100-year event
double tolerance = 0.01;      // ±1% of true quantile
double probability = 0.95;    // 95% confidence

int minN = MCMCDiagnostics.MinimumSampleSize(quantile, tolerance, probability);

Console.WriteLine($"Minimum sample size needed:");
Console.WriteLine($"  For {quantile:P1} quantile");
Console.WriteLine($"  With ±{tolerance:P1} tolerance");
Console.WriteLine($"  At {probability:P0} confidence");
Console.WriteLine($"  Need N ≥ {minN}");

// Adjust MCMC iterations accordingly
```

## Practical Diagnostics Workflow

### Complete Convergence Check

```cs
using Numerics.Sampling.MCMC;
using Numerics.Data.Statistics;

// Step 1: Run sampler
var sampler = new ARWMH(priors, logLikelihood);
sampler.NumberOfChains = 4;
sampler.WarmupIterations = 2000;
sampler.Iterations = 5000;
sampler.ThinningInterval = 10;
sampler.Sample();

Console.WriteLine("MCMC Convergence Diagnostics");
Console.WriteLine("=" + new string('=', 60));

// Step 2: Extract samples (flatten all chains)
var samples = sampler.Output.SelectMany(chain => chain).ToList();
int nParams = samples[0].Values.Length;

Console.WriteLine($"\nSampling Summary:");
Console.WriteLine($"  Chains: {sampler.NumberOfChains}");
Console.WriteLine($"  Warmup: {sampler.WarmupIterations}");
Console.WriteLine($"  Iterations: {sampler.Iterations}");
Console.WriteLine($"  Thinning: {sampler.ThinningInterval}");
Console.WriteLine($"  Total samples: {samples.Count}");

// Step 3: Check Gelman-Rubin
Console.WriteLine($"\nGelman-Rubin Statistics:");
var chains = sampler.MarkovChains.ToList();
double[] rhat = MCMCDiagnostics.GelmanRubin(chains, sampler.WarmupIterations);

bool converged = true;
for (int i = 0; i < nParams; i++)
{
    string status = rhat[i] < 1.1 ? "✓" : "✗";
    Console.WriteLine($"  {status} θ{i}: R̂ = {rhat[i]:F4}");
    if (rhat[i] >= 1.1) converged = false;
}

// Step 4: Check ESS
Console.WriteLine($"\nEffective Sample Size:");
double[] ess = MCMCDiagnostics.EffectiveSampleSize(chains, out double[][,] acf);

int minESS = (int)ess.Min();
// Use MarkovChains total sample count for efficiency (ESS is computed from MarkovChains, not Output)
int chainSampleCount = chains.Sum(c => c.Count);
for (int i = 0; i < nParams; i++)
{
    double efficiency = ess[i] / chainSampleCount;
    Console.WriteLine($"  θ{i}: ESS = {ess[i]:F0} ({efficiency:P1} efficiency)");
}

// Step 5: Overall assessment
Console.WriteLine($"\nOverall Assessment:");
if (converged && minESS > 100)
    Console.WriteLine("  ✓ PASS: Chains converged, sufficient samples");
else if (!converged)
    Console.WriteLine("  ✗ FAIL: Chains not converged - run longer");
else
    Console.WriteLine("  ⚠ MARGINAL: Converged but low ESS - consider more iterations");

// Step 6: Recommendations
if (minESS < 100)
{
    int recommended = (int)(sampler.Iterations * 100.0 / minESS);
    Console.WriteLine($"\nRecommendation: Increase iterations to ~{recommended}");
}
```

## Visual Diagnostics

While the library doesn't provide plotting, these are essential checks:

### Trace Plots

Plot parameter values vs. iteration:

```cs
Console.WriteLine("Export trace data for plotting:");
Console.WriteLine("Iteration | Parameter Values");

for (int i = 0; i < Math.Min(samples.Count, 100); i++)
{
    Console.Write($"{i,9} | ");
    foreach (var val in samples[i].Values)
    {
        Console.Write($"{val,8:F3} ");
    }
    Console.WriteLine();
}

Console.WriteLine("\nGood traces: Stationary, well-mixed 'hairy caterpillar'");
Console.WriteLine("Bad traces: Trending, stuck, or oscillating patterns");
```

### Posterior Distributions

```cs
// Export posterior samples for histogram
for (int param = 0; param < nParams; param++)
{
    var values = samples.Select(s => s.Values[param]).ToArray();
    
    Console.WriteLine($"\nParameter {param} summary:");
    Console.WriteLine($"  Mean: {values.Average():F4}");
    Console.WriteLine($"  Median: {Statistics.Percentile(values, 0.50):F4}");
    Console.WriteLine($"  SD: {Statistics.StandardDeviation(values):F4}");
    Console.WriteLine($"  95% CI: [{Statistics.Percentile(values, 0.025):F4}, " +
                     $"{Statistics.Percentile(values, 0.975):F4}]");
}
```

## Troubleshooting Convergence Issues

### Problem: High R̂ (> 1.1)

**Diagnosis:**
```cs
if (rhat.Max() > 1.1)
{
    Console.WriteLine("Convergence issue detected");
    Console.WriteLine("Possible causes:");
    Console.WriteLine("  1. Insufficient warmup");
    Console.WriteLine("  2. Poor mixing");
    Console.WriteLine("  3. Multimodal posterior");
    Console.WriteLine("  4. Bad initialization");
}
```

**Solutions:**
1. **Increase warmup**
   ```cs
   sampler.WarmupIterations = 5000;  // Double it
   ```

2. **Try different sampler**
   ```cs
   // Switch from RWMH to DEMCz for better mixing
   var betterSampler = new DEMCz(priors, logLikelihood);
   ```

3. **Check initialization**
   ```cs
   // Ensure initial values are reasonable
   // Not too far from expected values
   ```

### Problem: Low ESS (< 100)

**Diagnosis:**
```cs
if (ess.Min() < 100)
{
    Console.WriteLine($"Low ESS detected: {ess.Min():F0}");
    Console.WriteLine("Chains are highly autocorrelated");
}
```

**Solutions:**
1. **Increase iterations**
   ```cs
   int factor = (int)Math.Ceiling(100.0 / ess.Min());
   sampler.Iterations *= factor;
   Console.WriteLine($"Suggested iterations: {sampler.Iterations}");
   ```

2. **Increase thinning**
   ```cs
   sampler.ThinningInterval = 20;  // Keep every 20th sample
   ```

3. **Improve sampler**
   ```cs
   // Use ARWMH instead of RWMH
   // Use DEMCz for high dimensions
   ```

### Problem: Multimodal Posterior

**Diagnosis:**
```cs
// Check if chains explore different modes
foreach (int param in new[] { 0, 1, 2 })
{
    var chainMeans = chains.Select(chain => 
        chain.Select(ps => ps.Values[param]).Average()).ToArray();
    
    double meanRange = chainMeans.Max() - chainMeans.Min();
    double overallSD = Statistics.StandardDeviation(
        samples.Select(s => s.Values[param]).ToArray());
    
    if (meanRange > 2 * overallSD)
    {
        Console.WriteLine($"Parameter {param} may have multiple modes");
        Console.WriteLine("  Chain means differ substantially");
    }
}
```

**Solutions:**
1. **Use population-based sampler**
   ```cs
   var demcz = new DEMCz(priors, logLikelihood);
   demcz.NumberOfChains = 10;  // More chains
   ```

2. **Check for model identification issues**

3. **Consider transforming parameters**

## Best Practices

### 1. Always Run Multiple Chains

```cs
// Minimum 4 chains
sampler.NumberOfChains = 4;

// Benefits:
// - Can compute R̂
// - Detect multimodality
// - More robust inference
```

### 2. Adequate Warmup

```cs
// Rule of thumb: Warmup ≥ 50% of iterations
sampler.WarmupIterations = Math.Max(2000, sampler.Iterations / 2);
```

### 3. Check Diagnostics BEFORE Inference

```cs
// Always check before using samples
bool ready = (rhat.Max() < 1.1) && (ess.Min() > 100);

if (!ready)
{
    Console.WriteLine("⚠ WARNING: Do not use these samples!");
    Console.WriteLine("Run diagnostics and extend sampling");
    return;
}

// Proceed with inference...
```

### 4. Document Settings

```cs
Console.WriteLine("MCMC Configuration:");
Console.WriteLine($"  Sampler: {sampler.GetType().Name}");
Console.WriteLine($"  Chains: {sampler.NumberOfChains}");
Console.WriteLine($"  Warmup: {sampler.WarmupIterations}");
Console.WriteLine($"  Iterations: {sampler.Iterations}");
Console.WriteLine($"  Thinning: {sampler.ThinningInterval}");
Console.WriteLine($"  Total samples: {samples.Count}");
Console.WriteLine($"  R̂ range: [{rhat.Min():F3}, {rhat.Max():F3}]");
Console.WriteLine($"  ESS range: [{ess.Min():F0}, {ess.Max():F0}]");
```

### 5. Iterative Improvement

**Note:** Setting any sampler property (e.g., `Iterations`) calls `Reset()` internally, which clears all chains and output. This means increasing `Iterations` and calling `Sample()` runs a completely new simulation with the updated iteration count -- it does not extend the previous run.

```cs
int iteration = 1;
int totalIterations = sampler.Iterations;
while (rhat.Max() > 1.05 || ess.Min() < 200)
{
    totalIterations += 5000;
    Console.WriteLine($"\nIteration {iteration}: Restarting with {totalIterations} iterations...");

    // Setting Iterations triggers Reset(), clearing all previous chains/output.
    // Sample() then runs a fresh simulation with the new iteration count.
    sampler.Iterations = totalIterations;
    sampler.Sample();

    // Recompute diagnostics
    chains = sampler.MarkovChains.ToList();
    rhat = MCMCDiagnostics.GelmanRubin(chains, sampler.WarmupIterations);
    ess = MCMCDiagnostics.EffectiveSampleSize(chains, out _);

    iteration++;

    if (iteration > 5)
    {
        Console.WriteLine("Convergence issues persist - check model/sampler");
        break;
    }
}
```

## Summary

| Diagnostic | Target | Action if Not Met |
|------------|--------|------------------|
| **R̂** | < 1.1 | Increase warmup, run longer |
| **ESS** | > 100 (per param) | Increase iterations, improve mixing |
| **Visual traces** | Stationary | Check initialization, try different sampler |
| **ACF** | Drops quickly | Increase thinning, better sampler |

---

## References

<a id="1">[1]</a> Gelman, A., & Rubin, D. B. (1992). Inference from iterative simulation using multiple sequences. *Statistical Science*, 7(4), 457-472.

<a id="2">[2]</a> Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B., Vehtari, A., & Rubin, D. B. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.

<a id="3">[3]</a> Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical Science*, 7(4), 473-483.

<a id="4">[4]</a> Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: An improved R̂ for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667-718.

<a id="5">[5]</a> Flegal, J. M., Haran, M., & Jones, G. L. (2008). Markov chain Monte Carlo: Can we trust the third significant figure? *Statistical Science*, 23(2), 250-260.

---

[← Previous: MCMC Methods](mcmc.md) | [Back to Index](../index.md)
