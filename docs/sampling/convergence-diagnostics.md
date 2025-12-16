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

// Get chains
var chains = new List<List<ParameterSet>>();
// Extract chains from sampler output
// (Implementation depends on sampler structure)

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
```
R̂ = √(Var_total / W)

Where:
- W = within-chain variance (average)
- B = between-chain variance
- Var_total = ((n-1)/n)W + (1/n)B
```

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
// For single parameter series
double[] samples = /* Extract parameter samples from chain */;

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

```
ESS = N / (1 + 2·Σ ρ_k)

Where:
- N = number of samples
- ρ_k = autocorrelation at lag k
- Sum until ρ_k becomes negligible
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
    // Access from avgACF[parameter][chain, lag]
    double acf = avgACF[0][0, lag];  // Parameter 0, Chain 0
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

// Step 2: Extract samples
var samples = sampler.ParameterSets;
int nParams = samples[0].Values.Length;

Console.WriteLine($"\nSampling Summary:");
Console.WriteLine($"  Chains: {sampler.NumberOfChains}");
Console.WriteLine($"  Warmup: {sampler.WarmupIterations}");
Console.WriteLine($"  Iterations: {sampler.Iterations}");
Console.WriteLine($"  Thinning: {sampler.ThinningInterval}");
Console.WriteLine($"  Total samples: {samples.Length}");

// Step 3: Check Gelman-Rubin
Console.WriteLine($"\nGelman-Rubin Statistics:");
var chains = ExtractChains(sampler);  // Helper function
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
for (int i = 0; i < nParams; i++)
{
    double efficiency = ess[i] / samples.Length;
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

for (int i = 0; i < Math.Min(samples.Length, 100); i++)
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
    Console.WriteLine($"  Median: {Statistics.Percentile(values.OrderBy(v => v).ToArray(), 50):F4}");
    Console.WriteLine($"  SD: {Statistics.StandardDeviation(values):F4}");
    Console.WriteLine($"  95% CI: [{Statistics.Percentile(values.OrderBy(v => v).ToArray(), 2.5):F4}, " +
                     $"{Statistics.Percentile(values.OrderBy(v => v).ToArray(), 97.5):F4}]");
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
Console.WriteLine($"  Total samples: {samples.Length}");
Console.WriteLine($"  R̂ range: [{rhat.Min():F3}, {rhat.Max():F3}]");
Console.WriteLine($"  ESS range: [{ess.Min():F0}, {ess.Max():F0}]");
```

### 5. Iterative Improvement

```cs
int iteration = 1;
while (rhat.Max() > 1.05 || ess.Min() < 200)
{
    Console.WriteLine($"\nIteration {iteration}: Extending sampling...");
    
    sampler.Iterations += 5000;
    sampler.Sample();
    
    // Recompute diagnostics
    chains = ExtractChains(sampler);
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

---

[← Previous: MCMC Methods](mcmc.md) | [Back to Index](../index.md)
