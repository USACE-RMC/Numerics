---
title: 'Numerics: A .NET Library for Numerical Computing, Statistical Analysis, and Risk Assessment'
tags:
  - C#
  - .NET
  - numerical methods
  - statistics
  - probability distributions
  - Bayesian inference
  - machine learning
  - optimization
authors:
  - name: C. Haden Smith
    orcid: 0000-0002-4651-9890
    affiliation: 1
    corresponding: true
  - name: Woodrow L. Fields
    affiliation: 1
  - name: Julian Gonzalez
    affiliation: 1
  - name: Sadie Niblett
    affiliation: 1
  - name: Brennan Beam
    affiliation: 2
  - name: Brian Skahill
    orcid: 0000-0002-2164-0301
    affiliation: 3
affiliations:
  - name: U.S. Army Corps of Engineers, Risk Management Center, Lakewood, Colorado, USA
    index: 1
  - name: U.S. Army Corps of Engineers, Hydrologic Engineering Center, Davis, California, USA
    index: 2
  - name: Fariborz Maseeh Department of Mathematics and Statistics, Portland State University, Portland, Oregon, USA
    index: 3
date: 8 March 2026
bibliography: paper.bib
---

# Summary

Numerics is a free and open-source numerical computing library for .NET, developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). The library provides over 40 univariate probability distributions with L-moment and maximum likelihood parameter estimation, eight Markov Chain Monte Carlo (MCMC) samplers for Bayesian inference, bootstrap uncertainty quantification, bivariate copulas, optimization algorithms, and numerical methods for integration, differentiation, root finding, and linear algebra. Numerics targets engineers, scientists, and analysts working in infrastructure risk assessment, flood frequency analysis, and Monte Carlo simulation within .NET-based enterprise systems. The library supports .NET 8.0 through 10.0 and .NET Framework 4.8.1 and is validated by over 1,000 unit tests against published reference values.

# Statement of Need

Infrastructure risk assessment for dams and levees requires the integration of extreme value statistics, uncertainty quantification, and Monte Carlo simulation. Practitioners face challenges including small sample sizes from historical records, regulatory requirements for specific statistical methods such as the Log-Pearson Type III distribution mandated by USGS Bulletin 17C [@england2018], robust outlier handling via the Multiple Grubbs-Beck Test [@cohn2013], and the need to incorporate expert judgment through Bayesian inference.

The .NET framework is a widely used platform for enterprise and government desktop applications in the United States, yet the ecosystem lacks comprehensive numerical libraries tailored for engineering risk analysis. Python offers mature tools such as NumPy [@harris2020] and SciPy [@virtanen2020], and R provides specialized hydrology packages such as `lmom` [@hosking2019lmom], `nsRFA` [@viglione2024nsrfa], and `SpatialExtremes` [@ribatet2022spatialextremes], but neither ecosystem integrates with .NET-based engineering workflows without introducing language interop complexity, performance overhead, and deployment challenges for regulated government systems. In addition, Numerics benefits from .NET's ahead-of-time compilation, absence of a global interpreter lock, and native `Parallel.For` support, delivering significantly higher throughput for computationally intensive methods such as bootstrap resampling and MCMC sampling compared to interpreted alternatives.

Numerics fills this gap by providing domain-specific capabilities within .NET:

- **Hydrology-specific distributions**: Log-Pearson Type III, Generalized Extreme Value, Pearson Type III, Generalized Pareto, and Kappa-Four, with L-moment parameter estimation [@hosking1990; @hosking1997] that outperforms conventional moments for the small samples typical of flood records.
- **Mixed-population and competing-risk models**: Mixture distributions [@waylen1982; @alila2002] for sites where floods arise from distinct causal mechanisms such as rainfall, snowmelt, or dam-regulated releases, and competing-risks models [@crowder2001; @bedford2001] for reliability analysis with multiple failure modes.
- **Bayesian inference**: Eight MCMC samplers, including Random Walk Metropolis-Hastings, Adaptive RWMH [@haario2001], Differential Evolution MCMC (DE-MCz and DE-MCzs) [@terbraak2008], Hamiltonian Monte Carlo [@neal2011], No-U-Turn Sampler (NUTS) [@hoffman2014], Gibbs, and Self-Normalizing Importance Sampling (SNIS), with improved Gelman-Rubin convergence diagnostics [@gelmanrubin1992; @vehtari2021].
- **Uncertainty quantification**: Bootstrap resampling methods [@efron1993] for confidence intervals on design estimates, a requirement in dam and levee safety risk assessments.
- **Global optimization**: Differential Evolution [@storn1997], Shuffled Complex Evolution [@duan1994], Particle Swarm Optimization [@kennedy1995], Multi-Level Single-Linkage (MLSL) [@rinnooy1987], and Nelder-Mead [@nelder1965] for calibrating complex, multi-modal objective functions.
- **Adaptive numerical integration**: Gauss-Kronrod quadrature [@piessens1983] and VEGAS adaptive Monte Carlo integration [@lepage1978] for efficiently evaluating complex risk integrals involving multiple system components and failure modes, a core computation in quantitative risk assessment.
- **Machine learning**: Supervised and unsupervised algorithms, including generalized linear models, decision trees, random forests, k-nearest neighbors, k-means clustering, and Gaussian mixture models for regression, classification, and clustering tasks in risk assessment workflows.

# State of the Field

General-purpose .NET numerical libraries exist, most notably Math.NET Numerics [@mathnetnumerics], which provides linear algebra, probability distributions, and basic statistics. However, Math.NET Numerics does not include L-moment estimation, hydrology-specific distributions (Log-Pearson Type III, Kappa-Four), mixture and competing-risk models, adaptive integration, bootstrap resampling, copulas, or global optimization algorithms. Math.NET Numerics includes basic MCMC samplers (Metropolis-Hastings, Hybrid Monte Carlo, Slice), but lacks the adaptive and ensemble methods (Adaptive RWMH, DE-MCzs, NUTS) commonly applied for Bayesian inference of hydrologic models. Contributing these features to Math.NET Numerics was not pursued because the scope of domain-specific functionality, the distinct API design requirements for risk assessment workflows, and the need for long-term maintenance by a domain-expert organization warranted an independent library.

In other language ecosystems, Python's SciPy [@virtanen2020] provides broad numerical capabilities but lacks specialized hydrological distributions and L-moment estimation. R packages such as `lmom` [@hosking2019lmom] and `evd` [@stephenson2002evd] offer these features individually, but integrating R into .NET production applications introduces runtime dependencies and complicates deployment in regulated government environments. No single package in any ecosystem consolidates mixture distributions, competing-risk models, adaptive integration, extreme-value statistics, global optimization, bootstrapping, machine learning, and MCMC into a unified library for engineering risk analysis.

Numerics consolidates these capabilities into a single, self-contained .NET library purpose-built for quantitative risk assessment in water resources engineering. Computationally intensive operations, including bootstrap resampling and MCMC chain evaluation, are parallelized using `Parallel.For` for high-throughput execution on modern multi-core hardware. Backed by long-term USACE-RMC maintenance, Numerics provides the ecosystem, performance, and institutional support needed for quantitative risk analysis in dam and levee safety.

# Research Impact

Numerics serves as the computational engine for six USACE-RMC production applications used for infrastructure safety decisions:

- **[RMC-BestFit](https://github.com/USACE-RMC/RMC-BestFit)**: Bayesian estimation and fitting software that supports flood frequency analysis for flood risk management and dam and levee safety studies.
- **[RMC-RFA](https://github.com/USACE-RMC/RMC-RFA)**: Reservoir frequency analysis software to facilitate flood hazard assessments for dam safety studies.
- **[RMC-TotalRisk](https://github.com/USACE-RMC/RMC-TotalRisk)**: Quantitative risk analysis software designed to support dam and levee safety investment decisions.
- **[LifeSim](https://github.com/USACE-RMC/LifeSim)**: A consequence estimation engine used to simulate population redistribution during an evacuation.
- **[Levee Screening Tool](https://lst2.sec.usace.army.mil/)**: A web-based portfolio-level levee risk screening tool that has been used to evaluate more than 7,500 levee segments across the national inventory.
- **[Dam Screening Tool](https://dst.sec.usace.army.mil/)**: A web-based portfolio-level dam risk screening tool that will be used to evaluate more than 90,000 dams across the national inventory.

These applications support risk-informed decisions for thousands of dams and levees managed by USACE, other federal agencies, and private dam and levee owners, directly affecting public safety for millions of Americans. The RMC software suite is also used internationally by organizations in Australia, the Netherlands, and Canada. By open-sourcing Numerics, USACE-RMC enables independent verification of the algorithms underpinning these critical assessments and invites contributions from the broader water resources engineering community.

# Software Design

The library is organized around abstract base classes that define consistent interfaces for each computational domain. `UnivariateDistributionBase` provides PDF, CDF, inverse CDF, log-likelihood evaluation, and random variate generation for all 43 univariate distributions. Parameter estimation is decoupled from distribution classes through segregated interfaces (`IMaximumLikelihoodEstimation`, `ILinearMomentEstimation`, `IMomentEstimation`), so each distribution declares which estimation methods it supports and analysts can swap candidate distributions without changing workflow code. `Optimizer` provides a common minimization interface shared by all local and global optimization algorithms, with built-in function evaluation tracking, convergence testing, and numerical Hessian computation. `MCMCSampler` manages chain initialization, parallel execution via `Parallel.For`, thinning, and posterior output collection for all eight samplers, while convergence diagnostics are handled by a separate `MCMCDiagnostics` utility class. Bivariate copulas follow the same pattern through `BivariateCopula` and its `ArchimedeanCopula` specialization.

The library has zero external runtime dependencies, relying only on .NET Base Class Libraries, which eliminates NuGet dependency conflicts in large government applications and simplifies deployment in regulated environments. Multi-targeting (.NET 8.0, 9.0, 10.0, and .NET Framework 4.8.1) ensures compatibility with both modern and legacy enterprise systems. Because .NET 8 and later are cross-platform, Numerics runs on Windows, Linux, and macOS, supporting both desktop engineering tools and cloud-hosted web applications. Python users can also access Numerics through PythonNet, making the library's performance and capabilities available to the broader scientific and academic community.

Continuous integration via GitHub Actions runs the test suite (over 1,000 tests validated against published references [@press2007; @rao2000; @england2018]) on every pull request, preventing regressions in numerical accuracy.

# Example

The following example illustrates a bootstrap uncertainty analysis [@efron1993] for flood frequency quantiles, a common workflow in dam and levee safety studies. A Log-Pearson Type III distribution is fit to a synthetic sample using L-moments, and bias-corrected and accelerated (BCa) confidence intervals are computed for a range of non-exceedance probabilities:

```csharp
using Numerics.Distributions;

// Generate a synthetic flood record from a Log-Pearson Type III distribution
int n = 80;
int B = 10000;
var lp3 = new LogPearsonTypeIII(meanOfLog: 3, standardDeviationOfLog: 0.5, skewOfLog: 0.25);
double[] sample = lp3.GenerateRandomValues(n, seed: 12345);

// Estimate distribution parameters from the sample using L-moments
lp3.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);

// Configure the bootstrap analysis
var boot = new BootstrapAnalysis(distribution: lp3,
                                 estimationMethod: ParameterEstimationMethod.MethodOfLinearMoments,
                                 sampleSize: n,
                                 replications: B,
                                 seed: 12345);

// Define non-exceedance probabilities for quantile estimation
var probabilities = new double[] { 0.999, 0.99, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };

// Compute 90% BCa confidence intervals for each quantile
var CIs = boot.BCaQuantileCI(sample, probabilities, alpha: 0.1);
```

# AI Usage Disclosure

Generative AI was used to assist with XML documentation comments, markdown documentation content, and code review during library development. All core library code, class architecture, numerical methods, and algorithms were designed and implemented by the authors. All AI-generated content was reviewed, edited, and validated by the authors, who made all design decisions and accept full responsibility for the work.

# Acknowledgements

The authors acknowledge the U.S. Army Corps of Engineers Risk Management Center for supporting the development and open-source release of this software. The Numerics software would not exist without the support of Risk Management Center leadership, in particular the RMC Director Bryant A. Robbins and RMC Lead Engineers David A. Margo (retired) and John F. England. The authors are grateful to all who have contributed to the development of Numerics and the content of this paper.

# References
