# Goodness-of-Fit Metrics

The `GoodnessOfFit` class provides comprehensive metrics for evaluating model performance, including information criteria, error metrics, efficiency coefficients, and statistical tests [[1]](#ref1) [[2]](#ref2).

## Overview

| Category | Metrics |
|----------|---------|
| Information Criteria | AIC, AICc, BIC |
| Error Metrics | RMSE, MAE, MSE, MBE |
| Efficiency Metrics | NSE, KGE, R², d |
| Bias Metrics | PBIAS, RSR |
| Statistical Tests | Kolmogorov-Smirnov, Anderson-Darling, Chi-Squared |

---

## Information Criteria

Information criteria balance model fit against complexity, penalizing additional parameters to prevent overfitting [[3]](#ref3).

### Akaike Information Criterion (AIC)

```math
\text{AIC} = -2\ln(L) + 2k
```

where $L$ is the maximum likelihood and $k$ is the number of parameters.

```cs
using Numerics.Data.Statistics;

int k = 3;  // Number of parameters
double logLikelihood = -125.5;

double aic = GoodnessOfFit.AIC(k, logLikelihood);
Console.WriteLine($"AIC: {aic:F2}");
```

### Corrected AIC (AICc)

For small sample sizes (when $n/k < 40$):

```math
\text{AICc} = \text{AIC} + \frac{2k^2 + 2k}{n - k - 1}
```

```cs
int n = 30;  // Sample size
double aicc = GoodnessOfFit.AICc(n, k, logLikelihood);
Console.WriteLine($"AICc: {aicc:F2}");
```

### Bayesian Information Criterion (BIC)

```math
\text{BIC} = -2\ln(L) + k\ln(n)
```

BIC penalizes complexity more heavily than AIC for large samples.

```cs
double bic = GoodnessOfFit.BIC(n, k, logLikelihood);
Console.WriteLine($"BIC: {bic:F2}");
```

### Model Comparison

**Lower values indicate better fit** for all information criteria.

```cs
// Compare multiple distributions
var distributions = new[]
{
    (Name: "Normal", Dist: new Normal()),
    (Name: "LogNormal", Dist: new LogNormal()),
    (Name: "GEV", Dist: new GeneralizedExtremeValue()),
    (Name: "Gumbel", Dist: new Gumbel())
};

Console.WriteLine("Distribution      AIC       BIC");
Console.WriteLine("------------      ---       ---");

foreach (var (name, dist) in distributions)
{
    dist.SetParameters(dist.ParametersFromLinearMoments(data));
    double ll = dist.LogLikelihood(data);
    int k = dist.NumberOfParameters;
    
    double aic = GoodnessOfFit.AIC(k, ll);
    double bic = GoodnessOfFit.BIC(data.Length, k, ll);
    
    Console.WriteLine($"{name,-15}   {aic,7:F1}   {bic,7:F1}");
}
```

### AIC Weights

Convert AIC values to probabilities:

```cs
double[] aicValues = { 125.3, 128.1, 130.5, 126.2 };
double minAIC = aicValues.Min();

// Delta AIC
double[] deltaAIC = aicValues.Select(a => a - minAIC).ToArray();

// Relative likelihoods
double[] relLikelihood = deltaAIC.Select(d => Math.Exp(-0.5 * d)).ToArray();

// AIC weights (sum to 1)
double sumRL = relLikelihood.Sum();
double[] aicWeights = relLikelihood.Select(r => r / sumRL).ToArray();

Console.WriteLine("Model   ΔAIC    Weight");
for (int i = 0; i < aicValues.Length; i++)
{
    Console.WriteLine($"  {i+1}    {deltaAIC[i],5:F1}    {aicWeights[i]:F3}");
}
```

---

## Error Metrics

### Root Mean Square Error (RMSE)

```math
\text{RMSE} = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(O_i - P_i)^2}
```

where $O_i$ are observed values and $P_i$ are predicted values.

```cs
double[] observed = { 100, 150, 120, 180, 140 };
double[] predicted = { 105, 145, 125, 175, 138 };

double rmse = GoodnessOfFit.RMSE(observed, predicted);
Console.WriteLine($"RMSE: {rmse:F2}");
```

**Interpretation**: RMSE is in the same units as the data. Lower is better, with 0 indicating perfect fit.

### Mean Absolute Error (MAE)

```math
\text{MAE} = \frac{1}{n}\sum_{i=1}^{n}|O_i - P_i|
```

```cs
double mae = GoodnessOfFit.MAE(observed, predicted);
Console.WriteLine($"MAE: {mae:F2}");
```

**Interpretation**: Less sensitive to outliers than RMSE. Lower is better.

### Mean Squared Error (MSE)

```math
\text{MSE} = \frac{1}{n}\sum_{i=1}^{n}(O_i - P_i)^2
```

```cs
double mse = GoodnessOfFit.MSE(observed, predicted);
Console.WriteLine($"MSE: {mse:F2}");
```

### Mean Bias Error (MBE)

```math
\text{MBE} = \frac{1}{n}\sum_{i=1}^{n}(P_i - O_i)
```

```cs
double mbe = GoodnessOfFit.MBE(observed, predicted);
Console.WriteLine($"MBE: {mbe:F2}");
```

**Interpretation**: Positive MBE indicates overprediction; negative indicates underprediction.

---

## Efficiency Metrics

### Nash-Sutcliffe Efficiency (NSE)

The most widely used efficiency metric in hydrology [[4]](#ref4):

```math
\text{NSE} = 1 - \frac{\sum_{i=1}^{n}(O_i - P_i)^2}{\sum_{i=1}^{n}(O_i - \bar{O})^2}
```

```cs
double nse = GoodnessOfFit.NSE(observed, predicted);
Console.WriteLine($"NSE: {nse:F3}");
```

**Interpretation**:

| NSE | Performance |
|-----|-------------|
| NSE = 1 | Perfect match |
| 0.75 < NSE ≤ 1.0 | Very good |
| 0.65 < NSE ≤ 0.75 | Good |
| 0.50 < NSE ≤ 0.65 | Satisfactory |
| NSE ≤ 0.50 | Unsatisfactory |
| NSE < 0 | Worse than mean |

### Kling-Gupta Efficiency (KGE)

Decomposes model performance into correlation, bias, and variability [[5]](#ref5):

```math
\text{KGE} = 1 - \sqrt{(r-1)^2 + (\beta-1)^2 + (\gamma-1)^2}
```

where:
- $r$ = Pearson correlation coefficient
- $\beta = \mu_P / \mu_O$ (bias ratio)
- $\gamma = \text{CV}_P / \text{CV}_O$ (variability ratio)

```cs
double kge = GoodnessOfFit.KGE(observed, predicted);
Console.WriteLine($"KGE: {kge:F3}");

// Get KGE components
var kgeResult = GoodnessOfFit.KGEComponents(observed, predicted);
Console.WriteLine($"  Correlation (r): {kgeResult.Correlation:F3}");
Console.WriteLine($"  Bias ratio (β):  {kgeResult.BiasRatio:F3}");
Console.WriteLine($"  Variability (γ): {kgeResult.VariabilityRatio:F3}");
```

**Interpretation**: Similar to NSE, KGE = 1 is perfect. KGE > 0.5 is generally satisfactory.

### Coefficient of Determination (R²)

```math
R^2 = \left(\frac{\sum_{i=1}^{n}(O_i - \bar{O})(P_i - \bar{P})}{\sqrt{\sum_{i=1}^{n}(O_i - \bar{O})^2}\sqrt{\sum_{i=1}^{n}(P_i - \bar{P})^2}}\right)^2
```

```cs
double r2 = GoodnessOfFit.RSquared(observed, predicted);
Console.WriteLine($"R²: {r2:F3}");
```

**Interpretation**: R² ranges from 0 to 1. Measures linear association only; insensitive to bias.

### Index of Agreement (d)

Willmott's index of agreement [[6]](#ref6):

```math
d = 1 - \frac{\sum_{i=1}^{n}(O_i - P_i)^2}{\sum_{i=1}^{n}(|P_i - \bar{O}| + |O_i - \bar{O}|)^2}
```

```cs
double d = GoodnessOfFit.IndexOfAgreement(observed, predicted);
Console.WriteLine($"Index of Agreement: {d:F3}");
```

**Interpretation**: d ranges from 0 to 1, with 1 indicating perfect agreement.

---

## Bias Metrics

### Percent Bias (PBIAS)

```math
\text{PBIAS} = \frac{\sum_{i=1}^{n}(O_i - P_i)}{\sum_{i=1}^{n}O_i} \times 100
```

```cs
double pbias = GoodnessOfFit.PBIAS(observed, predicted);
Console.WriteLine($"PBIAS: {pbias:F1}%");
```

**Interpretation**:

| |PBIAS| | Performance |
|---------|-------------|
| < 10% | Very good |
| 10-15% | Good |
| 15-25% | Satisfactory |
| > 25% | Unsatisfactory |

Positive PBIAS indicates model underestimation.

### RMSE-Observations Standard Deviation Ratio (RSR)

```math
\text{RSR} = \frac{\text{RMSE}}{\sigma_O}
```

```cs
double rsr = GoodnessOfFit.RSR(observed, predicted);
Console.WriteLine($"RSR: {rsr:F3}");
```

**Interpretation**:

| RSR | Performance |
|-----|-------------|
| 0.00 ≤ RSR ≤ 0.50 | Very good |
| 0.50 < RSR ≤ 0.60 | Good |
| 0.60 < RSR ≤ 0.70 | Satisfactory |
| RSR > 0.70 | Unsatisfactory |

---

## Statistical Tests

### Kolmogorov-Smirnov Test

Tests whether data follow a specified distribution [[7]](#ref7):

```math
D = \max_i |F_n(x_i) - F(x_i)|
```

```cs
double[] data = { 12.5, 15.2, 11.8, 18.9, 14.2, 16.5, 13.4, 17.8 };

var normal = new Normal();
normal.SetParameters(normal.ParametersFromMoments(data));

var ksResult = GoodnessOfFit.KolmogorovSmirnov(data, normal);

Console.WriteLine($"K-S Statistic: {ksResult.Statistic:F4}");
Console.WriteLine($"P-value: {ksResult.PValue:F4}");
Console.WriteLine($"Reject H0 at α=0.05? {ksResult.PValue < 0.05}");
```

### Anderson-Darling Test

More sensitive to tail deviations than K-S [[7]](#ref7):

```math
A^2 = -n - \frac{1}{n}\sum_{i=1}^{n}(2i-1)[\ln F(x_{(i)}) + \ln(1-F(x_{(n+1-i)}))]
```

```cs
var adResult = GoodnessOfFit.AndersonDarling(data, normal);

Console.WriteLine($"A-D Statistic: {adResult.Statistic:F4}");
Console.WriteLine($"P-value: {adResult.PValue:F4}");
```

### Chi-Squared Test

Tests goodness-of-fit using binned data:

```cs
var chiResult = GoodnessOfFit.ChiSquared(data, normal, numberOfBins: 5);

Console.WriteLine($"χ² Statistic: {chiResult.Statistic:F4}");
Console.WriteLine($"Degrees of freedom: {chiResult.DegreesOfFreedom}");
Console.WriteLine($"P-value: {chiResult.PValue:F4}");
```

### Probability Plot Correlation Coefficient (PPCC)

Measures linearity of probability plot:

```cs
double ppcc = GoodnessOfFit.PPCC(data, normal);
Console.WriteLine($"PPCC: {ppcc:F4}");
```

**Interpretation**: PPCC close to 1 indicates good fit. Critical values depend on sample size and distribution.

---

## Comprehensive Model Evaluation

### Multi-Criteria Assessment

```cs
public static void EvaluateModel(double[] observed, double[] predicted, string modelName)
{
    Console.WriteLine($"=== {modelName} Performance ===");
    Console.WriteLine();
    
    // Error metrics
    Console.WriteLine("Error Metrics:");
    Console.WriteLine($"  RMSE:  {GoodnessOfFit.RMSE(observed, predicted):F3}");
    Console.WriteLine($"  MAE:   {GoodnessOfFit.MAE(observed, predicted):F3}");
    Console.WriteLine($"  MBE:   {GoodnessOfFit.MBE(observed, predicted):F3}");
    Console.WriteLine();
    
    // Efficiency metrics
    Console.WriteLine("Efficiency Metrics:");
    Console.WriteLine($"  NSE:   {GoodnessOfFit.NSE(observed, predicted):F3}");
    Console.WriteLine($"  KGE:   {GoodnessOfFit.KGE(observed, predicted):F3}");
    Console.WriteLine($"  R²:    {GoodnessOfFit.RSquared(observed, predicted):F3}");
    Console.WriteLine($"  d:     {GoodnessOfFit.IndexOfAgreement(observed, predicted):F3}");
    Console.WriteLine();
    
    // Bias metrics
    Console.WriteLine("Bias Metrics:");
    Console.WriteLine($"  PBIAS: {GoodnessOfFit.PBIAS(observed, predicted):F1}%");
    Console.WriteLine($"  RSR:   {GoodnessOfFit.RSR(observed, predicted):F3}");
}
```

### Distribution Selection Example

```cs
double[] annualMax = { 12500, 15200, 11800, 18900, 14200, 16500, 
                       13400, 17800, 19200, 10500, 21000, 14800 };

var candidates = new (string Name, IUnivariateDistribution Dist)[]
{
    ("Normal", new Normal()),
    ("LogNormal", new LogNormal()),
    ("GEV", new GeneralizedExtremeValue()),
    ("Gumbel", new Gumbel()),
    ("Log-Pearson III", new LogPearsonTypeIII())
};

Console.WriteLine("Distribution       AIC      BIC    K-S p    A-D p");
Console.WriteLine("------------       ---      ---    -----    -----");

foreach (var (name, dist) in candidates)
{
    // Fit distribution
    if (dist is ILinearMomentEstimation lmom)
        dist.SetParameters(lmom.ParametersFromLinearMoments(annualMax));
    else if (dist is IMomentEstimation mom)
        dist.SetParameters(mom.ParametersFromMoments(annualMax));
    
    // Information criteria
    double ll = dist.LogLikelihood(annualMax);
    int k = dist.NumberOfParameters;
    int n = annualMax.Length;
    double aic = GoodnessOfFit.AIC(k, ll);
    double bic = GoodnessOfFit.BIC(n, k, ll);
    
    // Statistical tests
    var ks = GoodnessOfFit.KolmogorovSmirnov(annualMax, dist);
    var ad = GoodnessOfFit.AndersonDarling(annualMax, dist);
    
    Console.WriteLine($"{name,-17} {aic,7:F1} {bic,7:F1}  {ks.PValue,6:F3}  {ad.PValue,6:F3}");
}
```

---

## Performance Rating Guidelines

Summary of performance ratings from Moriasi et al. [[1]](#ref1) [[2]](#ref2):

| Performance | RSR | NSE | PBIAS (Streamflow) |
|-------------|-----|-----|-------------------|
| Very Good | 0.00 ≤ RSR ≤ 0.50 | 0.75 < NSE ≤ 1.00 | PBIAS < ±10% |
| Good | 0.50 < RSR ≤ 0.60 | 0.65 < NSE ≤ 0.75 | ±10% ≤ PBIAS < ±15% |
| Satisfactory | 0.60 < RSR ≤ 0.70 | 0.50 < NSE ≤ 0.65 | ±15% ≤ PBIAS < ±25% |
| Unsatisfactory | RSR > 0.70 | NSE ≤ 0.50 | PBIAS ≥ ±25% |

---

## Summary Table

| Metric | Optimal | Range | Interpretation |
|--------|---------|-------|----------------|
| AIC/BIC | Lower | $-\infty$ to $\infty$ | Lower is better |
| RMSE | 0 | 0 to $\infty$ | Same units as data |
| MAE | 0 | 0 to $\infty$ | Same units as data |
| NSE | 1 | $-\infty$ to 1 | >0.5 satisfactory |
| KGE | 1 | $-\infty$ to 1 | >0.5 satisfactory |
| R² | 1 | 0 to 1 | >0.5 moderate |
| d | 1 | 0 to 1 | >0.5 acceptable |
| PBIAS | 0 | $-\infty$ to $\infty$ | <±25% satisfactory |
| RSR | 0 | 0 to $\infty$ | <0.7 satisfactory |

---

## References

<a id="ref1">[1]</a> Moriasi, D. N., Arnold, J. G., Van Liew, M. W., Bingner, R. L., Harmel, R. D., & Veith, T. L. (2007). Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. *Transactions of the ASABE*, 50(3), 885-900.

<a id="ref2">[2]</a> Moriasi, D. N., Gitau, M. W., Pai, N., & Daggupati, P. (2015). Hydrologic and water quality models: Performance measures and evaluation criteria. *Transactions of the ASABE*, 58(6), 1763-1785.

<a id="ref3">[3]</a> Burnham, K. P., & Anderson, D. R. (2002). *Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach* (2nd ed.). Springer.

<a id="ref4">[4]</a> Nash, J. E., & Sutcliffe, J. V. (1970). River flow forecasting through conceptual models part I—A discussion of principles. *Journal of Hydrology*, 10(3), 282-290.

<a id="ref5">[5]</a> Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling. *Journal of Hydrology*, 377(1-2), 80-91.

<a id="ref6">[6]</a> Legates, D. R., & McCabe, G. J. (1999). Evaluating the use of "goodness-of-fit" measures in hydrologic and hydroclimatic model validation. *Water Resources Research*, 35(1), 233-241.

<a id="ref7">[7]</a> D'Agostino, R. B., & Stephens, M. A. (1986). *Goodness-of-Fit Techniques*. Marcel Dekker.
