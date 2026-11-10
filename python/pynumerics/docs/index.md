# pynumerics

pynumerics brings the RMC.Numerics library to Python. The wheel ships the `Numerics.dll` assembly
of the published `RMC.Numerics` NuGet package, a loader that starts one .NET runtime per process,
and typed wrappers over every Numerics area: distributions, statistics, sampling and MCMC,
numerical methods, optimization, machine learning and data handling. Every value is computed by
Numerics; arrays cross between NumPy and .NET in bulk, bit for bit.

Start with the getting-started guide, then read the page for your area. The API reference lists
every public class and function with its full docstring.

```{toctree}
:maxdepth: 2
:caption: Guide

getting-started
```

```{toctree}
:maxdepth: 2
:caption: Areas

distributions
multivariate
regional
statistics
sampling
calculus
optimization
special-linalg-fourier
data
download
machine_learning
```

```{toctree}
:maxdepth: 1
:caption: Reference

api
```
