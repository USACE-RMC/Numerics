# Regional Frequency Analysis

[← Previous: Multivariate Distributions](multivariate.md) | [Back to Index](../index.md) | [Next: Univariate Functions →](../functions/index.md)

The `Numerics.Distributions.Regional` namespace implements the regional frequency analysis of Hosking and Wallis [[1]](#1) [[2]](#2): the index-flood procedure with L-moments. Sites whose frequency distributions are the same apart from a site-specific scale factor (the index flood, here the site mean) pool their data, so a regional growth curve is estimated from many more years than any single site holds. The statistics are those of Hosking's `REGTST` Fortran routine and of R's `lmomRFA::regtst`.

| Class | Purpose |
|---|---|
| `RegionalData` | Site names, record lengths, means and L-moment ratios; regional averages and discordancy |
| `RegionalHomogeneityTest` | Heterogeneity measures H1, H2, H3 and goodness-of-fit measures Z by kappa simulation |
| `RegionalGrowthCurve` | Growth curve fitted to the regional L-moments, growth factors and site quantiles |

## Regional Data

A site carries its record length $n_i$, its sample mean $l_1^{(i)}$, its L-CV $t^{(i)} = l_2^{(i)} / l_1^{(i)}$ and its L-moment ratios $t_3^{(i)}$, $t_4^{(i)}$ and $t_5^{(i)}$. `RegionalData.FromSamples` computes them from each site's raw sample with `Statistics.LinearMoments(sample, 5)` (unbiased sample L-moments; NaN values are missing years and are left out), and the constructor takes them precomputed.

The regional average L-moment ratios weight each site by its record length [[2]](#2) (section 4.1, equation 4.3):

```math
t^R = \frac{\sum_{i=1}^N n_i t^{(i)}}{\sum_{i=1}^N n_i}
```

and likewise for $t_3^R$, $t_4^R$ and $t_5^R$. `RegionalAverage()` returns $\{1, t^R, t_3^R, t_4^R, t_5^R\}$, the regional L-moments scaled to $l_1 = 1$.

## Discordancy

With $\mathbf{u}_i = (t^{(i)}, t_3^{(i)}, t_4^{(i)})^T$, their unweighted mean $\bar{\mathbf{u}}$ and $\mathbf{A} = \sum_i (\mathbf{u}_i - \bar{\mathbf{u}})(\mathbf{u}_i - \bar{\mathbf{u}})^T$, the discordancy of site $i$ is [[2]](#2) (section 3.2.3):

```math
D_i = \frac{N}{3} (\mathbf{u}_i - \bar{\mathbf{u}})^T \mathbf{A}^{-1} (\mathbf{u}_i - \bar{\mathbf{u}})
```

The $D_i$ average 1 and cannot exceed $(N - 1)/3$. When $\mathbf{A}$ is singular (the ratios lie in a plane) every $D_i$ is NaN; this differs deliberately from `REGTST`, which leaves $D_i = 0$, and agrees with lmomRFA, which reports NA. `RegionalData.DiscordancyCriticalValue(N)` gives the 10% critical value of Table 3.1 (1.333 at 5 sites rising to 3 from 15 sites on); a 5% table is available as well. Regions of three or fewer sites give $D_i = 1$.

## Heterogeneity

A four-parameter kappa distribution [[3]](#3) is fitted to $\{1, t^R, t_3^R, t_4^R\}$; when the kappa fit fails, which happens when $(t_3^R, t_4^R)$ lies above the generalized logistic line, the generalized logistic (the kappa with $h = -1$) stands in, as in `REGTST`. Each simulated region draws, independently at every site, as many values from that parent as the site's record length. Three dispersion statistics are measured about the regional averages:

```math
V_1 = \left\{ \frac{\sum n_i (t^{(i)} - t^R)^2}{\sum n_i} \right\}^{1/2}, \quad
V_2 = \frac{\sum n_i \{(t^{(i)} - t^R)^2 + (t_3^{(i)} - t_3^R)^2\}^{1/2}}{\sum n_i}, \quad
V_3 = \frac{\sum n_i \{(t_3^{(i)} - t_3^R)^2 + (t_4^{(i)} - t_4^R)^2\}^{1/2}}{\sum n_i}
```

and $H_j = (V_{j,\text{obs}} - \mu_V)/\sigma_V$, with $\mu_V$ and $\sigma_V$ the mean and sample standard deviation of the simulated values [[1]](#1) [[2]](#2) (section 4.3.3). A region is acceptably homogeneous for $H < 1$, possibly heterogeneous for $1 \le H < 2$ and definitely heterogeneous for $H \ge 2$; $H_1$ is the recommended measure.

The simulations default to 500. Simulated region $m$ uses its own `MersenneTwister`, seeded by the $m$-th integer drawn from a generator seeded with the master seed, so a seed gives the same results whatever the number of threads. The results agree with Hosking's within Monte Carlo error, not exactly, because the random number generators differ. The sites are simulated as independent.

## Goodness of Fit

For each of the generalized logistic, generalized extreme value, generalized normal, Pearson type III and generalized Pareto distributions fitted to $\{1, t^R, t_3^R\}$ [[2]](#2) (section 5.2.3):

```math
Z^{\text{DIST}} = \frac{\tau_4^{\text{DIST}} - t_4^R + B_4}{\sigma_4}
```

where $\tau_4^{\text{DIST}}$ is the fitted distribution's L-kurtosis, $B_4$ is the mean simulated regional $t_4^R$ less the observed $t_4^R$, and $\sigma_4$ is the standard deviation of the simulated regional $t_4^R$, all from the same kappa simulations. A fit is acceptable at the 10% level when $|Z| \le 1.64$.

## Growth Curve and Site Quantiles

`RegionalGrowthCurve` fits the chosen distribution to the regional L-moments with $l_1 = 1$, so its quantile function is the growth curve $q(F)$ with mean 1. The quantile of site $i$ is $\hat{Q}_i(F) = l_1^{(i)} \hat{q}(F)$ [[2]](#2) (section 6.2). The supported families are GLO, GEV, GNO, PE3, GPA and the kappa (which also uses $t_4^R$).

## Example: North Cascades

The 19-site North Cascades annual precipitation region of Hosking's LMOMENTS package (CASCADES.DAT, also lmomRFA's `Cascades`) reproduces the published CASCADES.OUT:

```cs
using Numerics.Distributions;
using Numerics.Distributions.Regional;

string[] names = { "350304", "351433", "351862", "351897", "352997", "353445", "353770", "356907", "357169", "357331",
                   "357354", "358466", "450945", "451233", "453284", "454764", "454769", "457773", "458773" };
int[] n = { 98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82 };
double[] mean = { 19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932,
                  17.560, 69.518, 47.653, 102.501, 52.413, 79.696, 44.643, 58.655, 39.024 };
double[] lcv = { 0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179,
                 0.1308, 0.1119, 0.1018, 0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046 };
double[] t3 = { 0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492,
                0.0940, -0.0429, 0.0435, 0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128 };
double[] t4 = { 0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900,
                0.1273, 0.0927, 0.1446, 0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583 };

var data = new RegionalData(names, n, mean, lcv, t3, t4);
double[] d = data.Discordancy();                      // 0.60, 1.02, ..., 2.63 (site 6), ..., 0.38
double[] r = data.RegionalAverage();                  // {1, 0.1103, 0.0279, 0.1366, NaN}

var test = new RegionalHomogeneityTest(data);         // 500 simulations, seed 12345
Console.WriteLine(string.Join(", ", test.KappaParameters));   // 0.9542, 0.1533, 0.1236, -0.2955
Console.WriteLine(string.Join(", ", test.Heterogeneity));     // about 0.55, -1.45, -2.31
Console.WriteLine(test.GetZ(UnivariateDistributionType.GeneralizedNormal)); // about -1.50
foreach (var accepted in test.AcceptableDistributions())
    Console.WriteLine(accepted);                      // GeneralizedNormal, PearsonTypeIII

var curve = new RegionalGrowthCurve(UnivariateDistributionType.GeneralizedNormal, data);
double[] f = { 0.5, 0.9, 0.99 };
double[] q = curve.GrowthFactors(f);                  // 0.994, 1.254, 1.480
double[,] sites = curve.SiteQuantiles(data, f);       // site mean times growth factor
```

The fixed values (D, regional means, kappa parameters, observed V, growth factors) match CASCADES.OUT to its printed precision. CASCADES.OUT reports H = 0.62, -1.49, -2.37 and Z(GNO) = -1.51 from its own 500 simulations; other seeds move H1 by a few hundredths and Z by about 0.1 to 0.3.

## Limitations

- The simulated sites are independent. Intersite correlation narrows the spread of the observed L-CVs, so a correlated region reads more homogeneous than it is (H below its independent-site value) [[2]](#2) (section 4.3).
- A simulated site with a zero or negative sample mean (a parent with a large L-CV and short records) has an infinite or negative L-CV and corrupts H, as in `REGTST`; `NonpositiveMeanSimulations` counts the simulated regions affected. A zero spread of the simulated statistics ($\sigma_V$ or $\sigma_4$ equal to 0) gives an infinite or NaN H or Z without an error.
- The Wakeby distribution is not available, so the Wakeby growth curve of CASCADES.OUT is not reproduced.
- The regional quantile-accuracy simulation (section 6.4) is not included.

## References

<a id="1">[1]</a> Hosking, J. R. M., & Wallis, J. R. (1993). Some statistics useful in regional frequency analysis. *Water Resources Research*, 29(2), 271-281. https://doi.org/10.1029/92WR01980

<a id="2">[2]</a> Hosking, J. R. M., & Wallis, J. R. (1997). *Regional Frequency Analysis: An Approach Based on L-Moments*. Cambridge University Press. https://doi.org/10.1017/CBO9780511529443

<a id="3">[3]</a> Hosking, J. R. M. (1994). The four-parameter kappa distribution. *IBM Journal of Research and Development*, 38(3), 251-258. https://doi.org/10.1147/rd.383.0251

[← Previous: Multivariate Distributions](multivariate.md) | [Back to Index](../index.md) | [Next: Univariate Functions →](../functions/index.md)
