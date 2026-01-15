# Numerics Python Demos - Project Scope

## Overview

This document outlines the scope for creating a PythonNet demo repository that showcases the Numerics library for Python users. The goal is to validate interest before investing in a full Python wrapper package.

---

## Repository

**Name:** `numerics-python-examples`

**Rationale:** Clear, searchable name that indicates these are examples/demos rather than a full wrapper package.

**Location:** Separate repository from the core Numerics library to:
- Keep demos independent from core library releases
- Make it easier for Python users to discover and clone
- Allow different contribution guidelines

---

## Repository Structure

```
numerics-python-examples/
├── README.md                      # Overview, installation, getting started
├── LICENSE                        # Same as Numerics (BSD-style USACE-RMC)
├── requirements.txt               # pythonnet, numpy, matplotlib, jupyter
├── environment.yml                # Conda environment (alternative)
├── SETUP.md                       # Detailed setup instructions
│
├── assets/
│   └── numerics.dll               # Pre-built Numerics DLL (or instructions)
│
├── src/
│   └── numerics_helper.py         # Shared helper functions for demos
│
├── notebooks/
│   ├── 00_getting_started.ipynb
│   ├── 01_distributions.ipynb
│   ├── 02_distribution_fitting.ipynb
│   ├── 03_mcmc_basics.ipynb
│   ├── 04_mcmc_bayesian_inference.ipynb
│   ├── 05_mcmc_advanced.ipynb
│   ├── 06_optimization.ipynb
│   ├── 07_statistics.ipynb
│   ├── 08_time_series.ipynb
│   └── 09_machine_learning.ipynb
│
├── examples/
│   ├── flood_frequency_analysis.py    # Real-world hydrology example
│   ├── reliability_analysis.py        # Engineering reliability
│   └── bayesian_regression.py         # Classic Bayesian example
│
└── docs/
    ├── api_cheatsheet.md          # Quick reference for Python users
    └── pythonnet_tips.md          # Common patterns and gotchas
```

---

## Jupyter Notebooks

### 00. Getting Started
**Purpose:** Environment setup and "Hello World"

**Content:**
- Installing .NET runtime (Windows/Linux/Mac)
- Installing PythonNet via pip
- Loading the Numerics DLL
- First example: create a Normal distribution, compute PDF/CDF
- Troubleshooting common issues

**Code Preview:**
```python
import clr
clr.AddReference("Numerics")
from Numerics.Distributions import Normal

dist = Normal(100, 15)
print(f"Mean: {dist.Mean}, Std: {dist.StandardDeviation}")
print(f"P(X < 120) = {dist.CDF(120):.4f}")
```

---

### 01. Distributions
**Purpose:** Tour of the 42+ univariate distributions

**Content:**
- Creating distributions (Normal, LogNormal, Gamma, Weibull, GEV, etc.)
- PDF, CDF, InverseCDF (quantile function)
- Generating random samples
- Statistical properties (mean, std, skew, kurtosis)
- Plotting distributions with matplotlib
- Discrete distributions (Poisson, Binomial)
- Special distributions (PERT, Triangular, Mixture)

**Real-world context:** "Which distribution for my data?"

---

### 02. Distribution Fitting
**Purpose:** Parameter estimation from data

**Content:**
- Method of Moments (MOM)
- Maximum Likelihood Estimation (MLE)
- L-Moments (Linear Moments)
- Comparing estimation methods
- Goodness-of-fit tests (KS, Anderson-Darling, Chi-squared)
- Plotting fitted vs empirical CDF
- Bootstrap confidence intervals on parameters
- Plotting positions (Weibull, Gringorten)

**Real-world context:** Flood frequency analysis - fitting annual peak flows

---

### 03. MCMC Basics
**Purpose:** Introduction to Bayesian inference with Numerics

**Content:**
- What is MCMC and why use it?
- Simple example: estimating mean of Normal distribution
- Setting up priors
- Defining log-likelihood function
- Running RWMH sampler
- Interpreting results (trace plots, posterior distributions)
- Burn-in and thinning

**Code Preview:**
```python
from Numerics.Sampling.MCMC import RWMH
from Numerics.Distributions import Normal, Uniform
from System import Func, Double
from System.Collections.Generic import List

# Priors
priors = List[UnivariateDistributionBase]()
priors.Add(Uniform(0, 200))   # Prior for mu
priors.Add(Uniform(1, 50))    # Prior for sigma

# Log-likelihood
def log_likelihood(params):
    mu, sigma = params[0], params[1]
    dist = Normal(mu, sigma)
    return sum(dist.LogPDF(x) for x in data)

# Create sampler
sampler = RWMH(priors, Func[List[Double], Double](log_likelihood))
sampler.Sample()
```

---

### 04. MCMC Bayesian Inference
**Purpose:** Practical Bayesian modeling workflows

**Content:**
- Linear regression with uncertainty
- Comparing different samplers (RWMH vs DEMCz vs HMC)
- **Performance benchmark: Numerics vs PyMC3** (key selling point)
- Posterior predictive checks
- Credible intervals
- Model comparison

**Real-world context:** Estimating a rating curve with uncertainty

---

### 05. MCMC Advanced
**Purpose:** Advanced MCMC techniques

**Content:**
- Multiple chains and convergence diagnostics (R-hat, ESS)
- Adaptive samplers (ARWMH)
- Hamiltonian Monte Carlo (HMC) with gradients
- Differential Evolution MCMC (DEMCz, DEMCzs)
- Hierarchical models
- Handling multimodal posteriors
- Tuning sampler parameters

---

### 06. Optimization
**Purpose:** Finding minima/maxima

**Content:**
- Local optimization
  - Nelder-Mead (gradient-free)
  - BFGS (quasi-Newton)
  - Powell's method
- Global optimization
  - Particle Swarm Optimization
  - Simulated Annealing
- Constrained optimization
- Practical tips: choosing initial points, handling bounds

**Real-world context:** Calibrating a model to observed data

---

### 07. Statistics
**Purpose:** Statistical analysis toolkit

**Content:**
- Descriptive statistics
- Correlation (Pearson, Spearman, Kendall)
- Hypothesis tests (t-test, Mann-Whitney)
- Bootstrap resampling
- Outlier detection (Multiple Grubbs-Beck)
- Box-Cox and Yeo-Johnson transformations
- Histogram and kernel density estimation

---

### 08. Time Series
**Purpose:** Working with temporal data

**Content:**
- Creating TimeSeries objects
- Smoothing (moving average, exponential, Gaussian)
- Block operations (annual max, monthly mean)
- Interpolation and resampling
- USGS data download integration
- Plotting time series

**Real-world context:** Analyzing streamflow records

---

### 09. Machine Learning
**Purpose:** ML capabilities in Numerics

**Content:**
- Random Forest (classification and regression)
- K-Means clustering
- Gaussian Mixture Models
- K-Nearest Neighbors
- Decision Trees
- When to use Numerics ML vs scikit-learn

---

## Standalone Python Examples

### flood_frequency_analysis.py
Complete workflow for hydrologic frequency analysis:
- Load annual peak flow data
- Fit multiple distributions (LogPearson III, GEV, LogNormal)
- Compare using goodness-of-fit
- Compute design flows with confidence intervals
- Generate flood frequency curve plot

### reliability_analysis.py
Engineering reliability example:
- Define limit state function
- Monte Carlo simulation
- MCMC for parameter uncertainty
- Compute probability of failure

### bayesian_regression.py
Classic statistics example:
- Linear regression with Bayesian inference
- Compare to frequentist results
- Posterior predictive distribution

---

## Helper Module

**File:** `src/numerics_helper.py`

Utility functions to reduce boilerplate in demos:

```python
"""Helper functions for using Numerics with Python."""

import clr
import numpy as np
from pathlib import Path

def load_numerics(dll_path=None):
    """Load the Numerics DLL and return the assembly."""
    if dll_path is None:
        dll_path = Path(__file__).parent.parent / "assets" / "Numerics.dll"
    clr.AddReference(str(dll_path))

def to_net_list(py_list, dtype=float):
    """Convert Python list to .NET List<T>."""
    from System.Collections.Generic import List
    from System import Double, Int32

    type_map = {float: Double, int: Int32}
    net_list = List[type_map[dtype]]()
    for item in py_list:
        net_list.Add(item)
    return net_list

def to_numpy(net_array):
    """Convert .NET array to numpy array."""
    return np.array(list(net_array))

def create_log_likelihood(data, dist_class):
    """Create a log-likelihood function for MCMC."""
    from System import Func, Double
    from System.Collections.Generic import List

    def log_lik(params):
        dist = dist_class(*[params[i] for i in range(params.Count)])
        return sum(dist.LogPDF(x) for x in data)

    return Func[List[Double], Double](log_lik)

def plot_distribution(dist, ax=None, **kwargs):
    """Plot PDF and CDF of a distribution."""
    import matplotlib.pyplot as plt
    # ... implementation
```

---

## Documentation Files

### README.md
- Project overview and purpose
- Quick start (5 minutes to first plot)
- Link to full Numerics documentation
- Performance comparison highlights
- Contributing guidelines

### SETUP.md
- Detailed installation for Windows/Linux/Mac
- .NET runtime installation
- Virtual environment setup
- Troubleshooting common PythonNet issues
- VS Code / PyCharm configuration

### docs/api_cheatsheet.md
Quick reference mapping common operations:

| Task | Numerics (C#) | Numerics via PythonNet |
|------|---------------|------------------------|
| Create Normal | `new Normal(0, 1)` | `Normal(0, 1)` |
| Evaluate PDF | `dist.PDF(x)` | `dist.PDF(x)` |
| Random samples | `dist.GenerateRandomValues(100)` | `list(dist.GenerateRandomValues(100))` |
| Fit to data | `Normal.MLE(data)` | `Normal.MLE(to_net_list(data))` |

### docs/pythonnet_tips.md
- Converting between Python and .NET types
- Working with delegates/Func<>
- Handling .NET exceptions in Python
- Performance tips
- Common pitfalls and solutions

---

## Deliverables Summary

| Item | Count | Description |
|------|-------|-------------|
| Jupyter notebooks | 10 | Core demo content |
| Python examples | 3 | Standalone real-world scripts |
| Helper module | 1 | Shared utilities |
| Documentation files | 4 | README, SETUP, cheatsheet, tips |
| Config files | 3 | requirements.txt, environment.yml, LICENSE |

**Total: ~21 files**

---

## Release Plan

### MVP Release (Phase 1)
Priority deliverables to validate interest:

1. `00_getting_started.ipynb`
2. `01_distributions.ipynb`
3. `03_mcmc_basics.ipynb`
4. `src/numerics_helper.py`
5. `README.md`
6. `SETUP.md`
7. `requirements.txt`

### Phase 2 Release
Fitting and advanced MCMC (key differentiators):

8. `02_distribution_fitting.ipynb`
9. `04_mcmc_bayesian_inference.ipynb` (includes performance benchmark)
10. `examples/flood_frequency_analysis.py`
11. `docs/api_cheatsheet.md`

### Phase 3 Release
Complete the suite:

12. `05_mcmc_advanced.ipynb`
13. `06_optimization.ipynb`
14. `07_statistics.ipynb`
15. `08_time_series.ipynb`
16. `09_machine_learning.ipynb`
17. `examples/reliability_analysis.py`
18. `examples/bayesian_regression.py`
19. `docs/pythonnet_tips.md`
20. `environment.yml`

---

## Success Metrics

To determine if a full Python wrapper is warranted, track:

1. **GitHub stars/forks** on the demo repository
2. **Issues/questions** from users (indicates engagement)
3. **Download/clone counts**
4. **Community feedback** on desired features
5. **Requests for pip-installable package**

If significant interest develops, proceed with official `rmc-numerics` Python package development.

---

## Dependencies

### Python Requirements
```
pythonnet>=3.0.0
numpy>=1.20.0
matplotlib>=3.5.0
jupyter>=1.0.0
pandas>=1.3.0  # Optional, for data handling examples
```

### .NET Requirements
- .NET 6.0+ Runtime (or .NET Framework 4.8.1 on Windows)
- Numerics.dll (from NuGet or built from source)

---

## Notes

- All demos should work on Windows, Linux, and macOS
- Test on both .NET 6/8 and .NET Framework where applicable
- Include sample datasets for reproducible examples
- Each notebook should be self-contained (can run independently)
- Code should be well-commented for educational purposes
