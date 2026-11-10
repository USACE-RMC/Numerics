"""RMC.Numerics for Python: the Numerics assembly, the single CoreCLR loader and typed wrappers.

Importing the package does not start the .NET runtime; the first wrapper call does, or
call :func:`initialize` explicitly. Submodules:

``runtime``
    The CoreCLR loader (one runtime and one ``Numerics`` assembly per process) and
    :func:`add_assembly` for dependent assemblies.
``marshal``
    Bit-exact bulk copies between NumPy and .NET arrays.
``identity``
    The version and SHA-256 of the loaded ``Numerics`` assembly.
``distributions``, ``multivariate``, ``copulas``
    Distributions; array inputs use the Numerics array overloads.
``regional``
    Hosking-Wallis regional frequency analysis: discordancy, heterogeneity, goodness of fit and growth curves.
``data``, ``statistics``, ``sampling``, ``mathematics``, ``bootstrap``
    Paired data and interpolation, statistics, sampling designs and the general bootstrap,
    the numerical utilities and Fourier transforms, and uncertainty analysis.
``mcmc``
    The MCMC samplers, their results and the convergence diagnostics.
``hypothesis_tests``, ``goodness_of_fit``, ``probability``
    Hypothesis tests, goodness-of-fit measures and probability combination rules.
``integration``, ``differentiation``, ``roots``, ``ode``
    Numerical integration, numerical derivatives, root finding and ODE solvers.
``optimization``, ``graphs``
    Local and global optimizers, knapsack and dynamic programming, and network shortest paths.
``special``, ``linalg``
    Special functions and the Numerics matrix, vector and decomposition types.
``functions``, ``timeseries``
    Univariate functions, parameter link functions and time series.
``download``
    Downloads of USGS, NOAA GHCN-Daily, Water Survey of Canada and BOM records, and offline parsers
    for saved responses.
``machine_learning``
    Regression, classification, clustering and dimension reduction.
``callbacks``
    Python callables and Roslyn-compiled C# source as the .NET delegates Numerics takes.
``enums``
    Python enums for the Numerics enums the wrappers accept.
``errors``
    The translation of Numerics exceptions into Python exceptions, applied by every wrapper.
"""

from pynumerics._version import (
    NUMERICS_ASSEMBLY_VERSION,
    NUMERICS_BUILD_COMMAND,
    NUMERICS_SHA256,
    NUMERICS_SOURCE_COMMIT,
    NUMERICS_VERSION,
    __version__,
)
from pynumerics.identity import (
    NumericsIdentity,
    assert_numerics_identity,
    numerics_identity,
)
from pynumerics.runtime import (
    RuntimeLoadError,
    add_assembly,
    get_runtime_info,
    initialize,
    is_initialized,
    loaded_assemblies,
    runtime_dir,
)

__all__ = [
    "NUMERICS_ASSEMBLY_VERSION", "NUMERICS_BUILD_COMMAND", "NUMERICS_SHA256", "NUMERICS_SOURCE_COMMIT",
    "NUMERICS_VERSION", "NumericsIdentity", "RuntimeLoadError",
    "__version__", "add_assembly", "assert_numerics_identity", "get_runtime_info", "initialize", "is_initialized",
    "loaded_assemblies", "numerics_identity", "runtime_dir",
]
