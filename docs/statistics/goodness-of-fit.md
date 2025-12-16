# Goodness-of-Fit Metrics

[← Previous: Descriptive Statistics](descriptive.md) | [Back to Index](../index.md) | [Next: MCMC Sampling →](../sampling/mcmc.md)

Goodness-of-fit (GOF) metrics evaluate how well a statistical model fits observed data. The ***Numerics*** library provides comprehensive metrics for model selection, distribution fitting validation, and hydrological model evaluation.

## Model Selection Criteria

### Information Criteria

Information criteria balance model fit with complexity, penalizing additional parameters [[1]](#1):

```cs
using Numerics.Data.Statistics;
using Numerics.Distributions;

double[] data = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Fit a distribution
var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

// Compute log-likelihood
double logLikelihood = 0;
foreach (var x in data)
{
    logLikelihood += gev.LogPDF(x);
}

int n = data.Length;
int k = gev.NumberOfParameters;

// Akaike Information Criterion
double aic = GoodnessOfFit.AIC(k, logLikelihood);

// Corrected AIC (for small samples)
double aicc = GoodnessOfFit.AICc(n, k, logLikelihood);

// Bayesian Information Criterion
double bic = GoodnessOfFit.BIC(n, k, logLikelihood);

Console.WriteLine($"Model Selection Criteria:");
Console.WriteLine($"  AIC:  {aic:F2} (lower is better)");
Console.WriteLine($"  AICc: {aicc:F2} (small sample correction)");
Console.WriteLine($"  BIC:  {bic:F2} (stronger penalty for parameters)");
```

**Formulas:**
```
AIC  = 2k - 2·ln(L)
AICc = AIC + 2k(k+1)/(n-k-1)
BIC  = k·ln(n) - 2·ln(L)
```

Where:
- `k` = number of parameters
- `n` = sample size
- `L` = likelihood

### Comparing Multiple Models

```cs
double[] data = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

// Fit multiple distributions
var models = new[]
{
    ("GEV", new GeneralizedExtremeValue()),
    ("Gumbel", new Gumbel()),
    ("LogNormal", new LogNormal()),
    ("Normal", new Normal())
};

var results = new List<(string Name, double AIC, double BIC, double LogLik)>();

foreach (var (name, dist) in models)
{
    dist.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
    
    double logLik = data.Sum(x => dist.LogPDF(x));
    double aic = GoodnessOfFit.AIC(dist.NumberOfParameters, logLik);
    double bic = GoodnessOfFit.BIC(data.Length, dist.NumberOfParameters, logLik);
    
    results.Add((name, aic, bic, logLik));
}

// Sort by AIC (lower is better)
results.Sort((a, b) => a.AIC.CompareTo(b.AIC));

Console.WriteLine("Model Comparison:");
Console.WriteLine("Model       | Params | AIC     | BIC     | Log-Lik");
Console.WriteLine("--------------------------------------------------------");
foreach (var (name, aic, bic, logLik) in results)
{
    var dist = models.First(m => m.Item1 == name).Item2;
    Console.WriteLine($"{name,-11} | {dist.NumberOfParameters,6} | {aic,7:F2} | {bic,7:F2} | {logLik,7:F2}");
}

Console.WriteLine($"\nBest model by AIC: {results[0].Name}");
```

### AIC Weights

Compute relative model probabilities:

```cs
var aicValues = results.Select(r => r.AIC).ToList();
double[] aicWeights = GoodnessOfFit.AICWeights(aicValues);

Console.WriteLine("\nModel Weights (relative probabilities):");
for (int i = 0; i < results.Count; i++)
{
    Console.WriteLine($"{results[i].Name,-11}: {aicWeights[i]:P1}");
}
```

## Distribution Fit Metrics

### Root Mean Square Error (RMSE)

```cs
double[] observed = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

var gev = new GeneralizedExtremeValue();
gev.Estimate(observed, ParameterEstimationMethod.MethodOfLinearMoments);

// RMSE using empirical plotting positions
double rmse = GoodnessOfFit.RMSE(observed, gev);

Console.WriteLine($"RMSE: {rmse:F2}");

// With custom plotting positions
var plottingPos = PlottingPositions.Weibull(observed.Length);
double rmse2 = GoodnessOfFit.RMSE(observed, plottingPos, gev);

// With parameter penalty
double rmse3 = GoodnessOfFit.RMSE(observed, gev.InverseCDF(plottingPos).ToArray(), k: gev.NumberOfParameters);

Console.WriteLine($"RMSE (Weibull plotting): {rmse2:F2}");
Console.WriteLine($"RMSE (with penalty): {rmse3:F2}");
```

### RMSE Weights

For model averaging:

```cs
var rmseValues = new List<double>();

foreach (var (name, dist) in models)
{
    dist.Estimate(observed, ParameterEstimationMethod.MethodOfLinearMoments);
    rmseValues.Add(GoodnessOfFit.RMSE(observed, dist));
}

double[] rmseWeights = GoodnessOfFit.RMSEWeights(rmseValues);

Console.WriteLine("RMSE-based Weights:");
for (int i = 0; i < models.Length; i++)
{
    Console.WriteLine($"{models[i].Item1,-11}: {rmseWeights[i]:P1}");
}
```

### Mean Square Error (MSE) and Mean Absolute Error (MAE)

```cs
double[] observed = { 100, 105, 98, 110, 95 };
double[] modeled = { 102, 104, 99, 108, 96 };

double mse = GoodnessOfFit.MSE(observed, modeled);
double mae = GoodnessOfFit.MAE(observed, modeled);
double rmse = Math.Sqrt(mse);

Console.WriteLine($"MSE:  {mse:F2}");
Console.WriteLine($"RMSE: {rmse:F2}");
Console.WriteLine($"MAE:  {mae:F2}");
```

## Hydrological Model Performance

### Nash-Sutcliffe Efficiency (NSE)

The NSE is widely used in hydrology [[2]](#2):

```cs
double[] observed = { 125, 135, 180, 220, 250, 280, 260, 230, 190, 150, 130, 120 };
double[] modeled = { 120, 140, 175, 225, 245, 275, 265, 225, 195, 145, 135, 115 };

double nse = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);

Console.WriteLine($"Nash-Sutcliffe Efficiency: {nse:F3}");

// Interpretation
if (nse >= 0.75)
    Console.WriteLine("Very good performance");
else if (nse >= 0.65)
    Console.WriteLine("Good performance");
else if (nse >= 0.50)
    Console.WriteLine("Satisfactory performance");
else if (nse >= 0.40)
    Console.WriteLine("Acceptable performance");
else
    Console.WriteLine("Unsatisfactory performance");
```

**Formula:**
```
NSE = 1 - Σ(O - M)² / Σ(O - Ō)²
```

Range: (-∞, 1], where 1 is perfect fit, 0 means model is as good as mean, <0 means worse than mean.

### Log Nash-Sutcliffe Efficiency

For better performance on low flows:

```cs
double logNSE = GoodnessOfFit.LogNashSutcliffeEfficiency(observed, modeled);

Console.WriteLine($"Log-NSE: {logNSE:F3}");
Console.WriteLine("Log-NSE emphasizes low flow performance");
```

### Kling-Gupta Efficiency (KGE)

Decomposes error into correlation, bias, and variability [[3]](#3):

```cs
double kge = GoodnessOfFit.KlingGuptaEfficiency(observed, modeled);

Console.WriteLine($"Kling-Gupta Efficiency: {kge:F3}");

// Modified KGE (variability ratio based on CV)
double kgeMod = GoodnessOfFit.KlingGuptaEfficiencyMod(observed, modeled);

Console.WriteLine($"Modified KGE: {kgeMod:F3}");

// KGE interpretation
if (kge >= 0.75)
    Console.WriteLine("Good model performance");
else if (kge >= 0.50)
    Console.WriteLine("Intermediate performance");
else
    Console.WriteLine("Poor performance");
```

**Formula:**
```
KGE = 1 - √[(r-1)² + (β-1)² + (γ-1)²]
```

Where:
- `r` = correlation coefficient
- `β` = bias ratio (μ_modeled / μ_observed)
- `γ` = variability ratio (CV_modeled / CV_observed)

### Percent Bias (PBIAS)

```cs
double pbias = GoodnessOfFit.PBIAS(observed, modeled);

Console.WriteLine($"Percent Bias: {pbias:F1}%");

// Interpretation
if (Math.Abs(pbias) < 10)
    Console.WriteLine("Very good (bias < ±10%)");
else if (Math.Abs(pbias) < 15)
    Console.WriteLine("Good (bias < ±15%)");
else if (Math.Abs(pbias) < 25)
    Console.WriteLine("Satisfactory (bias < ±25%)");
else
    Console.WriteLine("Unsatisfactory (bias ≥ ±25%)");

if (pbias > 0)
    Console.WriteLine("Model underestimates (positive bias)");
else if (pbias < 0)
    Console.WriteLine("Model overestimates (negative bias)");
```

### RMSE-Observations Standard Deviation Ratio (RSR)

```cs
double rsr = GoodnessOfFit.RSR(observed, modeled);

Console.WriteLine($"RSR: {rsr:F3}");

// Performance ratings
if (rsr <= 0.50)
    Console.WriteLine("Very good (RSR ≤ 0.50)");
else if (rsr <= 0.60)
    Console.WriteLine("Good (0.50 < RSR ≤ 0.60)");
else if (rsr <= 0.70)
    Console.WriteLine("Satisfactory (0.60 < RSR ≤ 0.70)");
else
    Console.WriteLine("Unsatisfactory (RSR > 0.70)");
```

### R-Squared (Coefficient of Determination)

```cs
double r2 = GoodnessOfFit.RSquared(observed, modeled);

Console.WriteLine($"R²: {r2:F3}");
Console.WriteLine($"Model explains {r2:P1} of variance");
```

### Index of Agreement

Willmott's index of agreement:

```cs
double ioa = GoodnessOfFit.IndexOfAgreement(observed, modeled);
double ioaMod = GoodnessOfFit.ModifiedIndexOfAgreement(observed, modeled);
double ioaRefined = GoodnessOfFit.RefinedIndexOfAgreement(observed, modeled);

Console.WriteLine($"Index of Agreement: {ioa:F3}");
Console.WriteLine($"Modified IoA: {ioaMod:F3}");
Console.WriteLine($"Refined IoA: {ioaRefined:F3}");

// Range: [0, 1] where 1 is perfect agreement
```

### Volumetric Efficiency

For water balance assessment:

```cs
double ve = GoodnessOfFit.VolumetricEfficiency(observed, modeled);

Console.WriteLine($"Volumetric Efficiency: {ve:F3}");

// VE = 1 - |Σ(O - M)| / Σ(O)
// Perfect when VE = 1 (volumes match exactly)
```

## Error Metrics

### Mean Absolute Percentage Error (MAPE)

```cs
double[] observed = { 100, 105, 98, 110, 95, 102 };
double[] modeled = { 102, 104, 99, 108, 96, 101 };

double mape = GoodnessOfFit.MAPE(observed, modeled);

Console.WriteLine($"MAPE: {mape:F2}%");

// Symmetric MAPE (better for asymmetric errors)
double smape = GoodnessOfFit.sMAPE(observed, modeled);

Console.WriteLine($"sMAPE: {smape:F2}%");
```

## Distribution-Specific Tests

### Kolmogorov-Smirnov Test

Tests if data comes from a specified distribution:

```cs
double[] data = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200 };

var gev = new GeneralizedExtremeValue();
gev.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

double ksStatistic = GoodnessOfFit.KolmogorovSmirnov(data, gev);

Console.WriteLine($"Kolmogorov-Smirnov D: {ksStatistic:F4}");
Console.WriteLine("Smaller D indicates better fit");

// Critical value at α=0.05 for n=8: approximately 0.457
// If D < critical value, fail to reject null hypothesis (good fit)
```

### Anderson-Darling Test

More sensitive to tail deviations:

```cs
double adStatistic = GoodnessOfFit.AndersonDarling(data, gev);

Console.WriteLine($"Anderson-Darling A²: {adStatistic:F4}");
Console.WriteLine("Smaller A² indicates better fit");
```

### Chi-Squared Test

For discrete or binned data:

```cs
double chiSq = GoodnessOfFit.ChiSquared(data, gev);

Console.WriteLine($"Chi-Squared χ²: {chiSq:F4}");
Console.WriteLine("Compare with χ² critical value from tables");
```

## Classification Metrics

For binary classification problems:

```cs
double[] observed = { 1, 0, 1, 1, 0, 1, 0, 0, 1, 1 };  // 1 = positive, 0 = negative
double[] predicted = { 1, 0, 1, 1, 0, 0, 1, 0, 1, 0 };

double accuracy = GoodnessOfFit.Accuracy(observed, predicted);
double precision = GoodnessOfFit.Precision(observed, predicted);
double recall = GoodnessOfFit.Recall(observed, predicted);
double f1 = GoodnessOfFit.F1Score(observed, predicted);
double specificity = GoodnessOfFit.Specificity(observed, predicted);
double balancedAcc = GoodnessOfFit.BalancedAccuracy(observed, predicted);

Console.WriteLine("Classification Metrics:");
Console.WriteLine($"  Accuracy: {accuracy:P1}");
Console.WriteLine($"  Precision: {precision:P1}");
Console.WriteLine($"  Recall (Sensitivity): {recall:P1}");
Console.WriteLine($"  F1 Score: {f1:F3}");
Console.WriteLine($"  Specificity: {specificity:P1}");
Console.WriteLine($"  Balanced Accuracy: {balancedAcc:P1}");
```

## Practical Examples

### Example 1: Complete Distribution Comparison

```cs
using Numerics.Data.Statistics;
using Numerics.Distributions;

double[] annualPeaks = { 12500, 15300, 11200, 18700, 14100, 16800, 13400, 17200, 10500, 19300 };

// Candidate distributions
var candidates = new (string Name, IUnivariateDistribution Dist)[]
{
    ("LP3", new LogPearsonTypeIII()),
    ("GEV", new GeneralizedExtremeValue()),
    ("Gumbel", new Gumbel()),
    ("LogNormal", new LogNormal())
};

Console.WriteLine("Distribution Comparison for Annual Peak Flows");
Console.WriteLine("=" + new string('=', 70));

var results = new List<(string Name, double AIC, double BIC, double RMSE, double KS, double AD)>();

foreach (var (name, dist) in candidates)
{
    // Fit distribution
    dist.Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);
    
    // Compute metrics
    double logLik = annualPeaks.Sum(x => dist.LogPDF(x));
    double aic = GoodnessOfFit.AIC(dist.NumberOfParameters, logLik);
    double bic = GoodnessOfFit.BIC(annualPeaks.Length, dist.NumberOfParameters, logLik);
    double rmse = GoodnessOfFit.RMSE(annualPeaks, dist);
    double ks = GoodnessOfFit.KolmogorovSmirnov(annualPeaks, dist);
    double ad = GoodnessOfFit.AndersonDarling(annualPeaks, dist);
    
    results.Add((name, aic, bic, rmse, ks, ad));
}

// Display results
Console.WriteLine("\nDistribution | Params | AIC     | BIC     | RMSE   | K-S    | A-D");
Console.WriteLine("-----------------------------------------------------------------------");

foreach (var r in results)
{
    var dist = candidates.First(c => c.Name == r.Name).Dist;
    Console.WriteLine($"{r.Name,-12} | {dist.NumberOfParameters,6} | {r.AIC,7:F1} | {r.BIC,7:F1} | " +
                     $"{r.RMSE,6:F0} | {r.KS,6:F4} | {r.AD,6:F4}");
}

// Recommendations
var bestAIC = results.OrderBy(r => r.AIC).First();
var bestBIC = results.OrderBy(r => r.BIC).First();
var bestRMSE = results.OrderBy(r => r.RMSE).First();

Console.WriteLine($"\nBest by AIC:  {bestAIC.Name}");
Console.WriteLine($"Best by BIC:  {bestBIC.Name}");
Console.WriteLine($"Best by RMSE: {bestRMSE.Name}");
```

### Example 2: Hydrological Model Evaluation

```cs
double[] observedFlow = { 125, 135, 180, 220, 250, 280, 260, 230, 190, 150, 130, 120 };
double[] modeledFlow = { 120, 140, 175, 225, 245, 275, 265, 225, 195, 145, 135, 115 };

Console.WriteLine("Hydrological Model Performance Evaluation");
Console.WriteLine("=" + new string('=', 50));

// Compute all metrics
double nse = GoodnessOfFit.NashSutcliffeEfficiency(observedFlow, modeledFlow);
double logNSE = GoodnessOfFit.LogNashSutcliffeEfficiency(observedFlow, modeledFlow);
double kge = GoodnessOfFit.KlingGuptaEfficiency(observedFlow, modeledFlow);
double pbias = GoodnessOfFit.PBIAS(observedFlow, modeledFlow);
double rsr = GoodnessOfFit.RSR(observedFlow, modeledFlow);
double r2 = GoodnessOfFit.RSquared(observedFlow, modeledFlow);
double rmse = GoodnessOfFit.RMSE(observedFlow, modeledFlow);
double mae = GoodnessOfFit.MAE(observedFlow, modeledFlow);

Console.WriteLine("\nPerformance Metrics:");
Console.WriteLine($"  NSE:     {nse,6:F3} {GetNSERating(nse)}");
Console.WriteLine($"  Log-NSE: {logNSE,6:F3}");
Console.WriteLine($"  KGE:     {kge,6:F3} {GetKGERating(kge)}");
Console.WriteLine($"  PBIAS:   {pbias,6:F1}% {GetPBIASRating(pbias)}");
Console.WriteLine($"  RSR:     {rsr,6:F3} {GetRSRRating(rsr)}");
Console.WriteLine($"  R²:      {r2,6:F3}");
Console.WriteLine($"  RMSE:    {rmse,6:F1} cfs");
Console.WriteLine($"  MAE:     {mae,6:F1} cfs");

// Helper functions for ratings
string GetNSERating(double nse)
{
    if (nse >= 0.75) return "(Very Good)";
    if (nse >= 0.65) return "(Good)";
    if (nse >= 0.50) return "(Satisfactory)";
    if (nse >= 0.40) return "(Acceptable)";
    return "(Unsatisfactory)";
}

string GetKGERating(double kge)
{
    if (kge >= 0.75) return "(Good)";
    if (kge >= 0.50) return "(Intermediate)";
    return "(Poor)";
}

string GetPBIASRating(double pbias)
{
    var abs = Math.Abs(pbias);
    if (abs < 10) return "(Very Good)";
    if (abs < 15) return "(Good)";
    if (abs < 25) return "(Satisfactory)";
    return "(Unsatisfactory)";
}

string GetRSRRating(double rsr)
{
    if (rsr <= 0.50) return "(Very Good)";
    if (rsr <= 0.60) return "(Good)";
    if (rsr <= 0.70) return "(Satisfactory)";
    return "(Unsatisfactory)";
}
```

### Example 3: Time Series Model Selection

```cs
double[] observed = { 100, 105, 98, 110, 95, 102, 108, 97, 103, 106 };
double[] model1 = { 102, 104, 99, 108, 96, 101, 107, 98, 102, 105 };
double[] model2 = { 101, 106, 97, 111, 94, 103, 109, 96, 104, 107 };
double[] model3 = { 100, 105, 98, 110, 95, 102, 108, 97, 103, 106 };

var models = new[] { ("Model 1", model1), ("Model 2", model2), ("Model 3", model3) };

Console.WriteLine("Time Series Model Comparison");
Console.WriteLine("=" + new string('=', 60));
Console.WriteLine("Model   | RMSE  | MAE   | MAPE  | NSE   | R²");
Console.WriteLine("-------------------------------------------------------");

foreach (var (name, modeled) in models)
{
    double rmse = Math.Sqrt(GoodnessOfFit.MSE(observed, modeled));
    double mae = GoodnessOfFit.MAE(observed, modeled);
    double mape = GoodnessOfFit.MAPE(observed, modeled);
    double nse = GoodnessOfFit.NashSutcliffeEfficiency(observed, modeled);
    double r2 = GoodnessOfFit.RSquared(observed, modeled);
    
    Console.WriteLine($"{name,-7} | {rmse,5:F2} | {mae,5:F2} | {mape,5:F2} | {nse,5:F3} | {r2,5:F3}");
}
```

## Best Practices

1. **Use multiple metrics** - No single metric captures all aspects of fit
2. **Information criteria** - Use AIC/BIC for model selection with different parameters
3. **Hydrological standards** - Follow Moriasi et al. [[2]](#2) criteria for hydrology
4. **Distribution tests** - Use K-S, A-D, or χ² for distribution validation
5. **Sample size** - Information criteria more reliable with n > 30
6. **Outliers** - Consider robust metrics like MAE instead of RMSE
7. **Context matters** - Different applications prioritize different metrics

---

## References

<a id="1">[1]</a> Akaike, H. (1974). A new look at the statistical model identification. *IEEE Transactions on Automatic Control*, 19(6), 716-723.

<a id="2">[2]</a> Moriasi, D. N., Arnold, J. G., Van Liew, M. W., Bingner, R. L., Harmel, R. D., & Veith, T. L. (2007). Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. *Transactions of the ASABE*, 50(3), 885-900.

<a id="3">[3]</a> Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling. *Journal of Hydrology*, 377(1-2), 80-91.

---

[← Previous: Descriptive Statistics](descriptive.md) | [Back to Index](../index.md) | [Next: MCMC Sampling →](../sampling/mcmc.md)
