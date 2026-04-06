# Goodness-of-Fit Metrics

[← Previous: Descriptive Statistics](descriptive.md) | [Back to Index](../index.md) | [Next: Hypothesis Tests →](hypothesis-tests.md)

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

```math
\text{AIC} = 2k - 2\ln(\hat{L})
```
```math
\text{AICc} = \text{AIC} + \frac{2k(k+1)}{n-k-1}
```
```math
\text{BIC} = k\ln(n) - 2\ln(\hat{L})
```

Where:
- `k` = number of parameters
- `n` = sample size
- `L` = likelihood

#### AIC and Kullback-Leibler Divergence

AIC is derived as an asymptotic approximation to the expected Kullback-Leibler (KL) divergence between the true data-generating model $f$ and the fitted candidate model $g$. The KL divergence measures the information lost when $g$ is used to approximate $f$:

```math
D_{KL}(f \| g) = \int f(x) \log \frac{f(x)}{g(x|\hat{\theta})} \, dx
```

Akaike (1974) showed that the expected KL divergence can be estimated, up to a constant that is the same for all candidate models, by:

```math
E[D_{KL}] \approx -2\log \hat{L} + 2k
```

Hence $\text{AIC} = 2k - 2\ln(\hat{L})$ estimates twice the expected information loss. The model with the lowest AIC is expected to be closest to the unknown true model in the KL sense.

#### AICc Small-Sample Correction

When the ratio $n/k < 40$, the standard AIC can exhibit substantial bias. The corrected AICc adds a finite-sample bias adjustment:

```math
\text{AICc} = \text{AIC} + \frac{2k(k+1)}{n-k-1}
```

This correction term becomes negligible as $n$ grows large, so AICc converges to AIC for large samples. It is generally recommended to use AICc by default, since the correction is harmless when the sample is large.

#### BIC and Bayesian Model Selection

BIC approximates the log marginal likelihood (or evidence) for a model under certain regularity conditions. It can be interpreted as an approximation to the Bayes factor between a candidate model $M_i$ and a saturated model:

```math
\text{BIC} \approx -2 \log P(\text{data} | M_i) + C
```

where $C$ is a constant common to all models. BIC penalizes model complexity more severely than AIC because the $\ln(n)$ term grows with sample size (compared to the constant penalty of 2 in AIC). As a result, BIC tends to prefer simpler models, especially for large $n$.

#### Interpreting ΔAIC

When comparing a set of candidate models, compute $\Delta_i = \text{AIC}_i - \text{AIC}_{\min}$ for each model. The following rules of thumb from Burnham and Anderson (2002) [[4]](#4) guide interpretation:

| ΔAIC | Evidence for model |
|------|-------------------|
| 0--2 | Substantial support; model is competitive with the best |
| 4--7 | Considerably less support |
| > 10 | Essentially no support; model is implausible |

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

// With custom plotting positions (observed data must be sorted in ascending order)
var sortedObserved = observed.OrderBy(x => x).ToArray();
var plottingPos = PlottingPositions.Weibull(sortedObserved.Length);
double rmse2 = GoodnessOfFit.RMSE(sortedObserved, plottingPos, gev);

// With parameter penalty (both arrays must be in the same order)
double rmse3 = GoodnessOfFit.RMSE(sortedObserved, gev.InverseCDF(plottingPos).ToArray(), k: gev.NumberOfParameters);

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

```math
\text{NSE} = 1 - \frac{\sum_{i=1}^{n}(O_i - M_i)^2}{\sum_{i=1}^{n}(O_i - \bar{O})^2}
```

Range: (-∞, 1], where 1 is perfect fit, 0 means model is as good as mean, <0 means worse than mean.

#### NSE Decomposition

NSE can be decomposed into three interpretable components (Murphy, 1988 [[7]](#7)):

```math
\text{NSE} = 2\alpha r - \alpha^2 - \beta^2
```

where:
- $r$ = Pearson correlation coefficient between observed and modeled values
- $\alpha = \sigma_M / \sigma_O$ (ratio of standard deviations, measuring variability bias)
- $\beta = (\mu_M - \mu_O) / \sigma_O$ (normalized bias)

This decomposition reveals whether poor NSE is caused by:
- **Poor correlation** ($r$ far from 1) -- indicates timing or phasing errors in the model
- **Wrong variability** ($\alpha$ far from 1) -- indicates the model over- or under-estimates the amplitude of fluctuations
- **Systematic bias** ($\beta$ far from 0) -- indicates consistent over- or under-prediction of the mean

### Log Nash-Sutcliffe Efficiency

The Log-NSE variant applies a logarithmic transformation to both observed and modeled values before computing NSE. This gives more weight to low-flow performance, since the logarithm compresses large values and expands small values. The transformation is:

```math
\text{Log-NSE} = 1 - \frac{\sum_{i=1}^{n}(\ln(O_i + \varepsilon) - \ln(M_i + \varepsilon))^2}{\sum_{i=1}^{n}(\ln(O_i + \varepsilon) - \overline{\ln(O + \varepsilon)})^2}
```

where $\varepsilon$ is a small constant (default: $\bar{O}/100$) added to prevent taking the logarithm of zero. Log-NSE is particularly useful for water quality assessment and ecological flow requirements, where accurate simulation of low-flow conditions is critical.

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

```math
\text{KGE} = 1 - \sqrt{(r-1)^2 + (\beta-1)^2 + (\gamma-1)^2}
```

Where:
- `r` = correlation coefficient
- `β` = bias ratio (μ_modeled / μ_observed)
- `γ` = variability ratio (CV_modeled / CV_observed)

#### KGE Component Interpretation

The KGE decomposes model error into three orthogonal components, each diagnosing a distinct type of model deficiency:

| Component | Meaning | Perfect Value | Diagnostic |
|-----------|---------|---------------|------------|
| $r$ | Pearson correlation | 1 | Timing and phasing of peaks and recessions |
| $\beta$ | Bias ratio ($\mu_M / \mu_O$) | 1 | Systematic over-estimation ($\beta > 1$) or under-estimation ($\beta < 1$) |
| Variability ratio | Spread of modeled vs. observed | 1 | Whether model reproduces the observed variability |

**Original KGE (Gupta et al., 2009)** [[3]](#3): The `KlingGuptaEfficiency` method uses the standard deviation ratio $\alpha = \sigma_M / \sigma_O$ as the variability component:

```math
\text{KGE} = 1 - \sqrt{(r-1)^2 + (\alpha - 1)^2 + (\beta - 1)^2}
```

**Modified KGE (Kling et al., 2012)** [[8]](#8): The `KlingGuptaEfficiencyMod` method replaces $\alpha$ with the coefficient of variation ratio $\gamma = \text{CV}_M / \text{CV}_O$, which avoids cross-correlation between the bias and variability components:

```math
\text{KGE'} = 1 - \sqrt{(r-1)^2 + (\gamma - 1)^2 + (\beta - 1)^2}
```

where $\gamma = (\sigma_M / \mu_M) / (\sigma_O / \mu_O)$. The modified version is generally preferred because the bias ratio $\beta$ and the variability ratio $\gamma$ are mathematically independent, making the decomposition cleaner.

**Diagnostic approach:** When KGE is low, examine which component dominates the Euclidean distance to guide model improvement:
- **Low $r$** -- improve model structure or parameterization to better capture temporal dynamics
- **$\beta \neq 1$** -- adjust bias correction or calibrate volume-related parameters
- **Variability ratio $\neq 1$** -- the model under- or over-estimates the spread; adjust parameters controlling flow variability

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
    Console.WriteLine("Model overestimates (positive bias)");
else if (pbias < 0)
    Console.WriteLine("Model underestimates (negative bias)");
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

The Kolmogorov-Smirnov (K-S) test statistic measures the maximum vertical distance between the empirical cumulative distribution function (ECDF) and the theoretical CDF:

```math
D_n = \sup_x |F_n(x) - F(x)|
```

where $F_n(x)$ is the empirical CDF defined as the proportion of observations less than or equal to $x$:

```math
F_n(x) = \frac{1}{n}\sum_{i=1}^{n} \mathbf{1}_{[x_i \leq x]}
```

and $F(x)$ is the CDF of the hypothesized distribution.

In practice, since the ECDF is a step function, the supremum reduces to checking each order statistic. The implementation computes:

```math
D_n = \max_{1 \leq i \leq n} \left\{ F(x_{i:n}) - \frac{i-1}{n}, \;\; \frac{i}{n} - F(x_{i:n}) \right\}
```

**Limitations:**
- Most sensitive near the center (median) of the distribution, less sensitive in the tails
- Critical values depend on whether parameters are estimated from data; when parameters are estimated, the Lilliefors correction is needed for valid p-values
- Conservative when used with estimated parameters (rejects too infrequently)
- Power decreases for heavy-tailed distributions
- The Anderson-Darling test is generally preferred when tail behavior is important

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

The Anderson-Darling (A-D) test statistic applies a weighting function that gives more emphasis to discrepancies in the tails of the distribution compared to the Kolmogorov-Smirnov test [[6]](#6). The test statistic is computed from the order statistics $x_{1:n} \leq x_{2:n} \leq \cdots \leq x_{n:n}$:

```math
A^2 = -n - \sum_{i=1}^{n} \frac{2i-1}{n}\left[\ln F(x_{i:n}) + \ln(1 - F(x_{n+1-i:n}))\right]
```

where $F$ is the CDF of the hypothesized distribution.

The implicit weighting function underlying the A-D statistic is $1/[F(x)(1-F(x))]$, which increases without bound as $F(x) \to 0$ or $F(x) \to 1$. This makes the test particularly sensitive to departures from the hypothesized distribution in the tails. For hydrological applications where tail behavior determines flood risk or drought severity, this tail sensitivity is a critical advantage over the K-S test.

**Key advantage over K-S:** The tail weighting makes the Anderson-Darling test more powerful for detecting misfit in the extreme quantiles that matter most for risk analysis.

More sensitive to tail deviations:

```cs
double adStatistic = GoodnessOfFit.AndersonDarling(data, gev);

Console.WriteLine($"Anderson-Darling A²: {adStatistic:F4}");
Console.WriteLine("Smaller A² indicates better fit");
```

### Chi-Squared Test

The Chi-Squared ($\chi^2$) goodness-of-fit test compares observed bin frequencies to the frequencies expected under the hypothesized distribution. The data are divided into $k$ bins, and the test statistic is:

```math
\chi^2 = \sum_{i=1}^{k} \frac{(O_i - E_i)^2}{E_i}
```

where $O_i$ is the observed frequency in bin $i$, and $E_i = n \cdot [F(b_i) - F(b_{i-1})]$ is the expected frequency computed from the hypothesized CDF $F$ evaluated at the bin boundaries $b_{i-1}$ and $b_i$.

Under the null hypothesis, $\chi^2$ follows approximately a chi-squared distribution with $\text{df} = k - 1 - p$ degrees of freedom, where $p$ is the number of parameters estimated from the data.

**Rules of thumb for bin selection:**
- Minimum expected frequency per bin should be at least 5 for the chi-squared approximation to be valid
- Too few bins reduces the power of the test; too many bins makes the test statistic unstable
- The implementation uses a default histogram binning; for formal hypothesis testing, consider the sensitivity of results to bin choice

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

This example compares candidate distributions for the White River near Nora, Indiana. Results can be cross-checked with the worked examples in Rao & Hamed (2000), Chapter 7.

**Data source:** Rao, A. R. & Hamed, K. H. (2000). *Flood Frequency Analysis*. CRC Press, Table 7.1.2.
See also: [`example-data/white-river-nora-floods.csv`](../example-data/white-river-nora-floods.csv)

```cs
using Numerics.Data.Statistics;
using Numerics.Distributions;

// White River near Nora, IN — 62 years of annual peak streamflow (cfs)
double[] annualPeaks = {
    23200, 2950, 10300, 23200, 4540, 9960, 10800, 26900, 23300, 20400,
    8480, 3150, 9380, 32400, 20800, 11100, 7270, 9600, 14600, 14300,
    22500, 14700, 12700, 9740, 3050, 8830, 12000, 30400, 27000, 15200,
    8040, 11700, 20300, 22700, 30400, 9180, 4870, 14700, 12800, 13700,
    7960, 9830, 12500, 10700, 13200, 14700, 14300, 4050, 14600, 14400,
    19200, 7160, 12100, 8650, 10600, 24500, 14400, 6300, 9560, 15800,
    14300, 28700
};

// Candidate distributions
var candidates = new (string Name, UnivariateDistributionBase Dist)[]
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
    ((IEstimation)dist).Estimate(annualPeaks, ParameterEstimationMethod.MethodOfLinearMoments);
    
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

## When to Use Which Metric

Selecting the right goodness-of-fit metric depends on the analysis goal. The following table summarizes recommended metrics for common tasks, all of which are available in the `GoodnessOfFit` class:

| Goal | Recommended Metric | Notes |
|------|-------------------|-------|
| Model selection (different number of parameters) | AIC / AICc / BIC | Information criteria; lower is better |
| Model selection (same number of parameters) | RMSE or NSE | Direct comparison of fit quality |
| Distribution fit validation | K-S, A-D, $\chi^2$ | A-D preferred for tail sensitivity |
| Hydrological model calibration | KGE | Decomposes error into correlation, bias, and variability |
| Hydrological model validation | NSE + PBIAS + RSR | Multiple complementary metrics recommended (Moriasi et al., 2007) [[2]](#2) |
| Forecast accuracy | MAPE / sMAPE | Scale-independent percentage errors |
| Water balance assessment | VE + PBIAS | Volume-focused metrics |
| Low-flow emphasis | Log-NSE | Log transform gives more weight to low values |
| Classification problems | F1 Score, Balanced Accuracy | When output is binary (0/1) |

### Moriasi et al. (2007) Performance Ratings

The following performance ratings from Moriasi et al. (2007) [[2]](#2) are widely used in watershed modeling to evaluate simulation results:

| Performance | NSE | RSR | PBIAS (streamflow) |
|------------|-----|-----|-------------------|
| Very Good | > 0.75 | 0.00--0.50 | < +/-10% |
| Good | 0.65--0.75 | 0.50--0.60 | +/-10--15% |
| Satisfactory | 0.50--0.65 | 0.60--0.70 | +/-15--25% |
| Unsatisfactory | < 0.50 | > 0.70 | >= +/-25% |

These thresholds apply to streamflow simulations at a monthly time step. For other constituents (e.g., sediment, nutrients) or sub-monthly time steps, the thresholds may be relaxed. No single metric fully characterizes model performance; using multiple complementary metrics provides a more complete assessment.

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

<a id="4">[4]</a> Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach* (2nd ed.). Springer.

<a id="5">[5]</a> Schwarz, G. (1978). Estimating the dimension of a model. *Annals of Statistics*, 6(2), 461-464.

<a id="6">[6]</a> Anderson, T. W., & Darling, D. A. (1954). A test of goodness of fit. *Journal of the American Statistical Association*, 49(268), 765-769.

<a id="7">[7]</a> Murphy, A. H. (1988). Skill scores based on the mean square error and their relationships to the correlation coefficient. *Monthly Weather Review*, 116(12), 2417-2424.

<a id="8">[8]</a> Kling, H., Fuchs, M., & Paulin, M. (2012). Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios. *Journal of Hydrology*, 424-425, 264-277.

---

[← Previous: Descriptive Statistics](descriptive.md) | [Back to Index](../index.md) | [Next: Hypothesis Tests →](hypothesis-tests.md)
