# Regional frequency analysis

`pynumerics.regional` wraps `Numerics.Distributions.Regional`, the regional frequency analysis of Hosking and Wallis
(1997): the index-flood procedure with L-moments, as in Hosking's `REGTST` Fortran routine and R's `lmomRFA::regtst`.
No Python package provides it.

- `RegionalData` holds each site's name, record length, mean (the index flood), L-CV and L-moment ratios t3, t4 and
  t5. Build it from raw samples with `RegionalData.from_samples` (NaN marks a missing year) or from precomputed
  L-moments. `regional_average()` is the record-length weighted regional average `(1, t^R, t3^R, t4^R, t5^R)`, and
  `discordancy()` the discordancy measure `D_i` of section 3.2.3, with `discordancy_critical_value` (Table 3.1: 3 from
  15 sites on).
- `RegionalHomogeneityTest` fits a four-parameter kappa to the regional L-moments (the generalized logistic stands in
  when the kappa fit fails), simulates 500 homogeneous regions with the sites' record lengths, and gives the
  heterogeneity measures `H1`, `H2`, `H3` (section 4.3.3) and the goodness-of-fit measures `Z` of GLO, GEV, GNO, PE3 and
  GPA with the bias correction `B4` and `sigma4` (section 5.2.3). Each simulated region takes its own seed drawn from the
  master seed, so a seed gives the same results whatever the number of threads. The sites are simulated as
  independent.
- `RegionalGrowthCurve` fits GLO, GEV, GNO, PE3, GPA or the kappa to the regional L-moments with mean 1, and gives growth
  factors at non-exceedance probabilities and site quantiles (index flood times growth factor).

Distributions are named by Hosking's abbreviation (`"GNO"`) or by the `UnivariateDistributionType` name
(`"GeneralizedNormal"`).

## The North Cascades example

The 19 sites of CASCADES.DAT in Hosking's LMOMENTS package reproduce CASCADES.OUT: the discordancy, regional means,
kappa parameters and growth curve to the printed precision, and H and Z within Monte Carlo error.

```python
import numpy as np

from pynumerics import regional

n = [98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82]
mean = [19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932, 17.560, 69.518, 47.653,
        102.501, 52.413, 79.696, 44.643, 58.655, 39.024]
l_cv = [0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179, 0.1308, 0.1119, 0.1018,
        0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046]
t3 = [0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492, 0.0940, -0.0429, 0.0435,
      0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128]
t4 = [0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900, 0.1273, 0.0927, 0.1446,
      0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583]
data = regional.RegionalData([f"site{i + 1}" for i in range(19)], n, mean, l_cv, t3, t4)

print(np.round(data.discordancy(), 2))            # 0.6 1.02 0.38 ... 2.63 (site 6) ... 0.38
print(np.round(data.regional_average()[:4], 4))   # [1. 0.1103 0.0279 0.1366]

test = regional.RegionalHomogeneityTest(data)     # 500 simulations, seed 12345
print(np.round(test.kappa_parameters, 4))         # [ 0.9542  0.1533  0.1236 -0.2955]
print(np.round(test.heterogeneity, 2))            # about [ 0.55 -1.45 -2.31]
print({k: round(v, 2) for k, v in test.z.items()})
print(test.acceptable_distributions())            # ['GeneralizedNormal', 'PearsonTypeIII']

curve = regional.RegionalGrowthCurve("GNO", data)
print(np.round(curve.growth_factors([0.5, 0.9, 0.99]), 3))   # [0.994 1.254 1.48 ]
print(curve.site_quantiles(data, [0.99]).shape)               # (19, 1)
```

## From raw samples

```python
import numpy as np

from pynumerics import regional

rng = np.random.default_rng(1)
samples = [rng.gumbel(100.0, 20.0, size=years) for years in (40, 55, 62, 48, 70, 35)]
samples[0][3] = np.nan                            # a missing year
data = regional.RegionalData.from_samples(samples)
print(data.record_lengths)                        # [39 55 62 48 70 35]
print(np.round(data.site_linear_moments(1), 4))   # (l1, l2, t3, t4, t5) of site 2

test = regional.RegionalHomogeneityTest(data, simulations=200, seed=7)
print(np.round(test.heterogeneity, 2))
gev = regional.RegionalGrowthCurve.from_ratios("GEV", *data.regional_average()[1:3])
print(np.round(gev.site_quantiles(data.mean[0], [0.99]), 1))
```
