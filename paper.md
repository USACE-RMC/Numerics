---
title: 'Numerics: A .NET Library for Numerical Computing and Statistical Analysis in Water Resources Engineering'
tags:
  - C#
  - .NET
  - statistics
  - probability distributions
  - MCMC
  - hydrology
  - flood frequency analysis
  - copulas
  - L-moments
  - optimization
authors:
  - name: Haden Smith
    orcid: 0000-0001-7881-5814
    affiliation: 1
affiliations:
  - name: U.S. Army Corps of Engineers, Risk Management Center, Lakewood, CO, USA
    index: 1
date: 3 March 2026
bibliography: paper.bib
---

# Summary

Numerics is a free and open-source .NET library providing numerical methods, probability distributions, statistical analysis, and Bayesian inference tools for quantitative risk assessment. Developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC), the library targets .NET 8, .NET 9, and .NET Framework 4.8.1. It is distributed via NuGet as `RMC.Numerics` and licensed under the BSD-3-Clause license.

The library implements 43 univariate probability distributions, 7 bivariate copulas, 8 Markov chain Monte Carlo (MCMC) samplers, 17 optimization algorithms, and supporting numerical infrastructure including integration, differentiation, root finding, linear algebra, and interpolation. Every univariate distribution supports probability density, cumulative distribution, quantile functions, and random variate generation. Parameter estimation is available via the method of moments, L-moments [@Hosking1990], and maximum likelihood estimation, with bootstrap confidence intervals for uncertainty quantification.

# Statement of Need

Flood frequency analysis, dam safety risk assessment, and water resources planning require statistical methods that are robust, well-tested, and accessible to practicing engineers. The USACE and other federal agencies rely on guidelines such as Bulletin 17C [@England2019] that prescribe specific distributional models (e.g., Log-Pearson Type III) and estimation procedures (e.g., Expected Moments Algorithm, multiple Grubbs-Beck test) not available in general-purpose numerical libraries.

Numerics was developed to fill this gap. It serves as the computational engine for several USACE-RMC software applications:

- **RMC-BestFit**: Bayesian mixture-model distribution fitting using MCMC with Dirichlet priors over mixture weights and copula-based bivariate models.
- **RMC-TotalRisk**: Dam and levee risk assessment integrating probabilistic hazard, fragility, and consequence models through Monte Carlo simulation.
- **RMC-RFA**: Regional flood frequency analysis using L-moments following @HoskingWallis1997.
- **LifeSim**: Life loss estimation for dam and levee failure scenarios using Monte Carlo simulation with competing risk models.

These applications require a unified library that combines classical frequentist estimation, Bayesian inference, copula dependence modeling, and numerical optimization within a single consistent API. Numerics provides this by implementing the full workflow from data retrieval and exploratory statistics through distribution fitting, model selection, uncertainty quantification, and joint probability analysis.

# Key Features

## Probability Distributions

The library implements 43 univariate distributions spanning continuous, discrete, and specialized families. Each distribution supports PDF, CDF, inverse CDF, random sampling, and summary statistics. Distributions commonly used in hydrology---including the Generalized Extreme Value, Log-Pearson Type III, Pearson Type III, and Generalized Pareto---implement L-moment estimation following @Hosking1990, which provides more robust parameter estimates than conventional moments for heavy-tailed distributions typical of flood data.

## Copulas

Seven bivariate copulas are available: Normal, Student-t, Clayton, Gumbel, Joe, Frank, and Ali-Mikhail-Haq. These enable modeling of dependence structures between random variables independently of their marginal distributions, following the framework of @Nelsen2006. The Student-t copula provides symmetric tail dependence, which is important for modeling joint extremes in dam safety applications. Parameter estimation is available via maximum pseudo-likelihood, inference functions for margins, and full maximum likelihood.

## MCMC Sampling

Eight MCMC samplers support Bayesian inference: Random Walk Metropolis-Hastings (RWMH), Adaptive RWMH [@Haario2001], Differential Evolution MCMC with snooker updates (DEMCz and DEMCzs) [@terBraak2008], Hamiltonian Monte Carlo (HMC) [@Neal2011], the No-U-Turn Sampler (NUTS) [@Hoffman2014], Gibbs sampling, and Sequential Nested Importance Sampling. All samplers support multi-chain execution, convergence diagnostics including split-$\hat{R}$ [@Vehtari2021], and configurable warmup, thinning, and adaptation. NUTS eliminates the manual step-count tuning required by HMC through an adaptive tree-building algorithm with dual averaging for step size selection.

## Optimization

Seventeen optimization algorithms include local methods (BFGS, Nelder-Mead, Powell, gradient descent, ADAM), global methods (Differential Evolution [@Storn1997], Simulated Annealing, Shuffled Complex Evolution [@Duan1994], Particle Swarm), and constrained optimization via augmented Lagrangian. These are used internally for maximum likelihood estimation and are available for user-defined objective functions.

# State of the Field

The most comparable libraries are SciPy's `scipy.stats` module in Python and MathNet.Numerics in .NET. SciPy provides broad statistical functionality but lacks L-moment estimation, USACE-specific distributions, and integrated MCMC. MathNet.Numerics offers core numerical methods for .NET but does not implement copulas, MCMC samplers, or hydrologic frequency analysis methods. R packages such as `lmom`, `copula`, and `rstan` provide individual components but require separate integration. Numerics combines all of these capabilities in a single library with a consistent API designed for the water resources engineering workflow.

# Acknowledgements

This work was funded by the U.S. Army Corps of Engineers Risk Management Center. The author thanks the USACE-RMC staff for testing and feedback. The views expressed herein are those of the author and do not necessarily reflect the views of the U.S. Army Corps of Engineers.

# References
