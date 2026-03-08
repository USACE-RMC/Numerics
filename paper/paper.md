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

<!-- General note, I do not use em dashes. Please use commas or parantheses instead -->

# Summary

Numerics is a free and open-source numerical computing library for .NET, developed by the U.S. Army Corps of Engineers Risk Management Center (USACE-RMC). The library provides over 40 univariate probability distributions with L-moment and maximum likelihood parameter estimation, six Markov Chain Monte Carlo (MCMC) samplers for Bayesian inference, bootstrap uncertainty quantification, bivariate copulas, optimization algorithms, and numerical methods for integration, differentiation, root finding, and linear algebra. Numerics targets engineers, scientists, and analysts working in infrastructure risk assessment, flood frequency analysis, and Monte Carlo simulation within .NET-based enterprise systems. The library is self-contained with zero external runtime dependencies, supports .NET 8.0 through 10.0 and .NET Framework 4.8.1, and is validated by over 1,000 unit tests against published reference values.

# Statement of Need

Infrastructure risk assessment for dams and levees requires the integration of extreme value statistics, uncertainty quantification, and Monte Carlo simulation. Practitioners face challenges including small sample sizes from historical records, regulatory requirements for specific statistical methods such as the Log-Pearson Type III distribution mandated by USGS Bulletin 17C [@england2019], robust outlier handling via the Multiple Grubbs-Beck Test [@cohn2013], and the need to incorporate expert judgment through Bayesian inference.

The .NET framework is the dominant platform for enterprise and government desktop applications in the United States, yet the ecosystem lacks comprehensive numerical libraries tailored for engineering risk analysis. Python offers mature tools such as NumPy [@harris2020] and SciPy [@virtanen2020], and R provides specialized hydrology packages <!-- add a reference to lmom, nsRFA, and SpatialExtremes packages -->, but neither ecosystem integrates with .NET-based engineering workflows without introducing language interop complexity, performance overhead, and deployment challenges for regulated government systems.

<!-- Need to stress that Numerics provides much more performant code than using R or Python. -->

Numerics fills this gap by providing domain-specific capabilities within .NET:

- **Hydrology-specific distributions**: Log-Pearson Type III, Generalized Extreme Value, Pearson Type III, Generalized Pareto, and Kappa-Four, with L-moment parameter estimation [@hosking1990; @hosking1997] that outperforms conventional moments for the small samples typical of flood records.
- **Bayesian inference**: Six MCMC samplers---Random Walk Metropolis-Hastings, Adaptive RWMH [@haario2001], Differential Evolution MCMC [@terbraak2008], Hamiltonian Monte Carlo [@neal2011], and Gibbs---with Gelman-Rubin convergence diagnostics [@vehtari2021]. <!-- Mention NUTS -->
- **Uncertainty quantification**: Bootstrap resampling methods [@efron1993] for confidence intervals on design estimates, a requirement in dam and levee safety risk assessments.
- **Global optimization**: Differential Evolution [@storn1997], Shuffled Complex Evolution [@duan1994], Particle Swarm Optimization [@kennedy1995], and Nelder-Mead [@nelder1965] for calibrating complex, multi-modal objective functions. <!-- Mention MLSL -->

# State of the Field

General-purpose .NET numerical libraries exist, most notably Math.NET Numerics [@mathnetnumerics], which provides linear algebra, probability distributions, and basic statistics. However, Math.NET Numerics does not include L-moment estimation, hydrology-specific distributions (Log-Pearson Type III, Kappa-Four), MCMC samplers <!-- they do provide MCMC but not the robust methods Numerics provides -->, bootstrap resampling, copulas, or global optimization algorithms---<!-- don't use em dashes --> capabilities that are essential for engineering risk analysis. Contributing these features to Math.NET Numerics was not pursued because the scope of domain-specific functionality, the distinct API design requirements for risk assessment workflows, and the need for long-term maintenance by a domain-expert organization warranted an independent library.

In other language ecosystems, Python's SciPy [@virtanen2020] provides broad numerical capabilities but lacks specialized hydrological distributions and L-moment estimation. R packages such as `lmom` [@hosking2019lmom] and `evd` [@stephenson2002evd] offer these features individually, but integrating R into .NET production applications introduces runtime dependencies and complicates deployment in regulated government environments.

Numerics consolidates these capabilities into a single, self-contained .NET library purpose-built for quantitative risk assessment in water resources engineering.

<!-- discuss that heavy computational methods like bootstrapping and MCMC are parallelized for performance. Numerics provides the ecosystem, long-term support, and performance needed for QRA in dam and levee saffety. -->

# Research Impact

Numerics serves as the computational engine for six USACE-RMC production applications used for infrastructure safety decisions:

<!-- improve bulleted text. Need to provide desciption to explain importance -->

- **[RMC-BestFit](https://github.com/USACE-RMC/RMC-BestFit)**: Bayesian estimation and fitting software that supports flood frequency analysis for flood risk management and dam and levee safety studies. 
- **[RMC-RFA](https://github.com/USACE-RMC/RMC-RFA)**: Reservoir frequency analysis software to facilitate flood hazard assessments for dam safety studies.
- **[RMC-TotalRisk](https://github.com/USACE-RMC/RMC-TotalRisk)**: Quantitative risk analysis software designed to support dam and levee safety investment decisions. 
- **[LifeSim](https://github.com/USACE-RMC/LifeSim)**: A consequence estimation engine used to simulate population redistribution during an evacuation.
- **[Levee Screening Tool](https://lst2.sec.usace.army.mil/)**: Portfolio-level levee risk screening across the national inventory. <!-- has been used to evaluate more than 7,500 levee segments -->
- **[Dam Screening Tool](https://dst.sec.usace.army.mil/)**: Portfolio-level dam risk screening across the national inventory. <!-- will be used to evaluate more than 90,000 dams-->

These applications support risk-informed decisions for thousands of dams and levees managed by USACE and other federal agencies <!-- and private dam and levee owners -->, directly affecting public safety for millions of Americans. By open-sourcing Numerics, USACE-RMC enables independent verification of the algorithms underpinning these critical assessments and invites contributions from the broader water resources engineering community.

<!-- our software is used globally as well in nations like Australian, the Netherlands, and Canada -->

# Software Design

Numerics is organized into namespaces that mirror its functional domains: `Numerics.Distributions` for probability distributions and copulas, `Numerics.Data.Statistics` for statistical analysis, `Numerics.Mathematics` for numerical methods (integration, optimization, root finding, linear algebra, ODE solvers), `Numerics.Sampling.MCMC` for Bayesian samplers, and `Numerics.MachineLearning` for supervised and unsupervised algorithms.

All univariate distributions inherit from a common base class that enforces a consistent API: probability density, cumulative distribution, inverse CDF (quantile), random variate generation, and parameter estimation via pluggable methods (L-moments, maximum likelihood, method of moments). This design allows analysts to switch between distributions without changing workflow code---a critical feature when comparing candidate models during frequency analysis.

The library has zero external runtime dependencies, relying only on .NET Base Class Libraries. This eliminates NuGet dependency conflicts in large government applications and simplifies deployment. Multi-targeting (.NET 8.0, 9.0, 10.0, and .NET Framework 4.8.1) ensures compatibility with both modern and legacy enterprise systems. Continuous integration via GitHub Actions runs the full test suite---over 1,000 tests validated against published references [@press2007; @rao2000; @england2019]---on every pull request.

<!-- 
need to reference the correct B17C reference. 
@techreport{USGS2018,
    title={Guidelines for Determining Flood Flow Frequency, Bulletin 17C},
    author={England, John F and Cohn, Timothy A and Faber, Beth A and Stedinger, Jery R and Thomas, W O and Veilleux, Andrea G and Kiang, Julie E and Mason, Robert R},
    institution = {U.S. Geological Survey},
    type = {Techniques and Methods},
    number = {4-B5},
    address = {Reston, VA, USA},
    year = {2018},
    doi = {10.3133/tm4B5}
} 
-->

<!-- Need to add example section. Maybe show running a bootstrap analysis with a LP3 distribution. Use mu = 3, sigma = 0.5, skew = 0.25 and N=80. Perform bootstrap with MOM and 10000 realizations etc. Show how to output BCa confidence intervals. I am open to ideas -->

# AI Usage Disclosure

Generative AI was used to assist with XML documentation comments, markdown documentation content, and code review during library development. All core library code, class architecture, numerical methods, and algorithms were designed and implemented by the authors. All AI-generated content was reviewed, edited, and validated by the authors, who made all design decisions and accept full responsibility for the work.

# Acknowledgements

The authors acknowledge the U.S. Army Corps of Engineers Risk Management Center for supporting the development and open-source release of this software.

<!-- Here some additional language for acknowledgements: The Numerics software would not exist without the support of Risk Management Center leadership, in~particular the RMC Director Bryant A. Robbins and RMC Lead Engineers David A. Margo (retired) and John F. England. The authors are grateful to all who have contributed to the development of Numerics and the content of this paper. -->

# References
