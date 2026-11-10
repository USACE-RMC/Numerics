# pynumerics

RMC.Numerics for Python. The wheel ships the `Numerics.dll` assembly from the published
`RMC.Numerics` NuGet package, byte for byte, with a loader that starts one .NET runtime and one
`Numerics` assembly per process, and typed wrappers over every Numerics area: univariate and
multivariate distributions, copulas, statistics, hypothesis tests, sampling, MCMC and the
bootstrap, integration, differentiation, root finding and ODEs, optimization, special functions,
linear algebra and Fourier transforms, data and time series, and machine learning.

The wrappers do no numerical work in Python. Arrays cross between NumPy and .NET in bulk, bit for
bit, every wrapper raises Python exceptions with the Numerics messages, and functions can be
Python callables or C# source compiled with the bundled Roslyn compiler. The distribution version
equals the RMC.Numerics version of the bundled assembly.

## Install

```bash
pip install pynumerics==2.3.0
```

The 2.3.0 wheel bundles `Numerics.dll` from the RMC.Numerics 2.3.0 NuGet package. Requirements:
CPython 3.10 to 3.14 on Linux or Windows, and the .NET 10 runtime (`Microsoft.NETCore.App` 10.0).
Set `DOTNET_ROOT` when .NET is not installed in the default location.

```python
from pynumerics import distributions

gev = distributions.GeneralizedExtremeValue(location=10.0, scale=2.0, shape=-0.1)
print(gev.inverse_cdf([0.5, 0.99, 0.999]))   # one call into InverseCDF(IList<double>)
```

## Documentation

The user guide and the API reference are in
[`docs/`](https://github.com/USACE-RMC/Numerics/tree/main/python/pynumerics/docs): start with
[getting started](https://github.com/USACE-RMC/Numerics/blob/main/python/pynumerics/docs/getting-started.md),
which covers installation, `PYNUMERICS_DLL_DIR` and runtime selection, arrays and errors, value
semantics, callbacks and migration from rmc-numerics. Build the HTML pages, including the API
reference generated from the docstrings, from `python/pynumerics`:

```bash
python -m pip install -e ".[docs]"
python -m sphinx -W --keep-going -b html docs docs/_build/html
```

Releases are published to PyPI by `.github/workflows/PythonPublish.yml` with trusted publishing:
a `pynumerics-vX.Y.Z` tag on a commit whose `_version.py` pins the published RMC.Numerics
package builds the wheel, tests it on Linux and Windows, and uploads it to TestPyPI and then PyPI.

## Build and test

From `python/pynumerics`:

```bash
python scripts/build_wheel.py --stage-only       # stage _bin/Numerics.dll and the Roslyn assemblies
python -m pip install -e ".[test]"
python -m pytest                                  # includes tests/test_docs.py: every docs code block
python scripts/build_wheel.py                    # build dist/pynumerics-<version>-py3-none-any.whl
```

The script takes `Numerics.dll` from `lib/net10.0` of the `RMC.Numerics` package in the
NuGet global packages folder (or `--nupkg PATH`, or `--download` from nuget.org) and refuses
any file whose SHA-256 differs from `NUMERICS_SHA256` in `src/pynumerics/_version.py`. It stages
`Microsoft.CodeAnalysis.dll` and `Microsoft.CodeAnalysis.CSharp.dll` from `lib/net9.0` of the
Microsoft.CodeAnalysis.Common and Microsoft.CodeAnalysis.CSharp 4.14.0 packages (NuGet cache, or
`--download`), each checked against its `ROSLYN_*_SHA256` pin; `--roslyn-only` stages only those.
To test against a locally built assembly, set `PYNUMERICS_DLL_DIR` to
`Numerics/bin/Release/net10.0`.
