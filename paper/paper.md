---
title: 'Numerics: A .NET Library for Numerical Computing, Statistical Analysis, and Risk Assessment'
tags:
  - C#
  - .NET
  - numerical methods
  - statistics
  - probability distributions
  - Bayesian inference
  - hydrology
  - risk assessment
authors:
  - name: C. Haden Smith
    orcid: 0000-0001-7881-5814
    affiliation: 1
    corresponding: true
  - name: Woodrow L. Fields
    affiliation: 1
  - name: Julian Gonzalez
    affiliation: 1
  - name: Sadie Niblett
    affiliation: 1
  - name: Brennan Beam
    affiliation: 1
  - name: Brian Skahill
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: U.S. Army Corps of Engineers, Risk Management Center, Lakewood, Colorado, USA
    index: 1
date: 8 March 2026
bibliography: paper.bib
---

# Summary

Numerics is a free and open-source numerical computing library for .NET, developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). The library provides over 40 univariate probability distributions with L-moment and maximum likelihood parameter estimation, six Markov Chain Monte Carlo (MCMC) samplers for Bayesian inference, bootstrap uncertainty quantification, bivariate copulas, optimization algorithms, and numerical methods for integration, differentiation, root finding, and linear algebra. Numerics targets engineers, scientists, and analysts working in infrastructure risk assessment, flood frequency analysis, and Monte Carlo simulation within .NET-based enterprise systems. The library is self-contained with zero external runtime dependencies, supports .NET 8.0 through 10.0 and .NET Framework 4.8.1, and is validated by over 1,000 unit tests against published reference values.

# Statement of Need

Infrastructure risk assessment for dams and levees requires the integration of extreme value statistics, uncertainty quantification, and Monte Carlo simulation. Practitioners face challenges including small sample sizes from historical records, regulatory requirements for specific statistical methods such as the Log-Pearson Type III distribution mandated by USGS Bulletin 17C [@england2019], robust outlier handling via the Multiple Grubbs-Beck Test [@cohn2013], and the need to incorporate expert judgment through Bayesian inference.

The .NET framework is the dominant platform for enterprise and government desktop applications in the United States, yet the ecosystem lacks comprehensive numerical libraries tailored for engineering risk analysis. Python offers mature tools such as NumPy [@harris2020] and SciPy [@virtanen2020], and R provides specialized hydrology packages, but neither ecosystem integrates with .NET-based engineering workflows without introducing language interop complexity, performance overhead, and deployment challenges for regulated government systems.

Numerics fills this gap by providing domain-specific capabilities within .NET:

- **Hydrology-specific distributions**: Log-Pearson Type III, Generalized Extreme Value, Pearson Type III, Generalized Pareto, and Kappa-Four, with L-moment parameter estimation [@hosking1990; @hosking1997] that outperforms conventional moments for the small samples typical of flood records.
- **Bayesian inference**: Six MCMC samplers---Random Walk Metropolis-Hastings, Adaptive RWMH [@haario2001], Differential Evolution MCMC [@terbraak2008], Hamiltonian Monte Carlo [@neal2011], and Gibbs---with Gelman-Rubin convergence diagnostics [@vehtari2021].
- **Uncertainty quantification**: Bootstrap resampling methods [@efron1993] for confidence intervals on design estimates, a requirement in dam safety risk assessments.
- **Global optimization**: Differential Evolution [@storn1997], Shuffled Complex Evolution [@duan1994], Particle Swarm Optimization [@kennedy1995], and Nelder-Mead [@nelder1965] for calibrating complex, multi-modal objective functions.

# State of the Field

General-purpose .NET numerical libraries exist, most notably Math.NET Numerics [@mathnetnumerics], which provides linear algebra, probability distributions, and basic statistics. However, Math.NET Numerics does not include L-moment estimation, hydrology-specific distributions (Log-Pearson Type III, Kappa-Four), MCMC samplers, bootstrap resampling, copulas, or global optimization algorithms---capabilities that are essential for engineering risk analysis. Contributing these features to Math.NET Numerics was not pursued because the scope of domain-specific functionality, the distinct API design requirements for risk assessment workflows, and the need for long-term maintenance by a domain-expert organization warranted an independent library.

In other language ecosystems, Python's SciPy [@virtanen2020] provides broad numerical capabilities but lacks specialized hydrological distributions and L-moment estimation. R packages such as `lmom` [@hosking2019lmom] and `evd` [@stephenson2002evd] offer these features individually, but integrating R into .NET production applications introduces runtime dependencies and complicates deployment in regulated government environments.

Numerics consolidates these capabilities into a single, self-contained .NET library purpose-built for quantitative risk assessment in water resources engineering.

# Software Design

Numerics is organized into namespaces that mirror its functional domains: `Numerics.Distributions` for probability distributions and copulas, `Numerics.Data.Statistics` for statistical analysis, `Numerics.Mathematics` for numerical methods (integration, optimization, root finding, linear algebra, ODE solvers), `Numerics.Sampling.MCMC` for Bayesian samplers, and `Numerics.MachineLearning` for supervised and unsupervised algorithms.

All univariate distributions inherit from a common base class that enforces a consistent API: probability density, cumulative distribution, inverse CDF (quantile), random variate generation, and parameter estimation via pluggable methods (L-moments, maximum likelihood, method of moments). This design allows analysts to switch between distributions without changing workflow code---a critical feature when comparing candidate models during frequency analysis.

The library has zero external runtime dependencies, relying only on .NET Base Class Libraries. This eliminates NuGet dependency conflicts in large government applications and simplifies deployment. Multi-targeting (.NET 8.0, 9.0, 10.0, and .NET Framework 4.8.1) ensures compatibility with both modern and legacy enterprise systems. Continuous integration via GitHub Actions runs the full test suite---over 1,000 tests validated against published references [@press2007; @rao2000; @england2019]---on every pull request.

# Research Impact

Numerics serves as the computational engine for six USACE-RMC production applications used for infrastructure safety decisions:

- **RMC-BestFit**: Distribution fitting and flood frequency analysis with uncertainty quantification
- **RMC-RFA**: Regional frequency analysis for ungauged sites using L-moment-based index flood methods
- **RMC-TotalRisk**: Comprehensive dam and levee risk quantification combining hazard, fragility, and consequence analyses
- **LifeSim**: Life loss estimation for dam and levee failure scenarios
- **Levee Screening Tool**: Portfolio-level levee risk screening across the national inventory
- **Dam Screening Tool**: Portfolio-level dam risk screening across the national inventory

These applications support risk-informed decisions for thousands of dams and levees managed by USACE and other federal agencies, directly affecting public safety for millions of Americans. By open-sourcing Numerics, USACE-RMC enables independent verification of the algorithms underpinning these critical assessments and invites contributions from the broader water resources engineering community.

# AI Usage Disclosure

Generative AI (Anthropic Claude) was used to assist with XML documentation comments, markdown documentation content, and code review during library development. The JOSS paper was drafted with AI assistance. All core library code, class architecture, numerical methods, and algorithms were designed and implemented by the human authors. All AI-generated content was reviewed, edited, and validated by the authors, who made all design decisions and accept full responsibility for the work.

# Acknowledgements

The authors acknowledge the U.S. Army Corps of Engineers Risk Management Center for supporting the development and open-source release of this software.

# References
