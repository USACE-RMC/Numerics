# Getting started

## Installation

```bash
pip install pynumerics
```

The wheel is pure Python and runs on CPython 3.11 to 3.14 on Linux and Windows. It bundles
`Numerics.dll` and the Roslyn C# compiler assemblies, so the one system requirement is the
.NET 10 runtime (`Microsoft.NETCore.App` 10.0). Set `DOTNET_ROOT` when .NET is not installed in
its default location. The Python dependencies, pythonnet and NumPy, install with the wheel.

## The runtime and the assembly

Importing `pynumerics` does not start .NET. The first wrapper call does, or you can call
`pynumerics.initialize` yourself to choose the assembly directory or to surface loading problems
early. One process hosts one .NET runtime and one `Numerics` assembly, found in this order:

1. the directory passed to `pynumerics.initialize(dll_dir)`;
2. the `PYNUMERICS_DLL_DIR` environment variable;
3. the `_bin` directory packaged in the wheel.

`Numerics.dll` is referenced there by its full path, so a copy elsewhere on `sys.path` is never
loaded instead. To run against a locally built assembly, point the environment variable at the
build output before Python starts:

```bash
export PYNUMERICS_DLL_DIR=/path/to/Numerics/bin/Release/net10.0
```

```powershell
$env:PYNUMERICS_DLL_DIR = "C:\path\to\Numerics\bin\Release\net10.0"
```

After the first successful call, `initialize` returns the directory it loaded from and ignores
later arguments. `numerics_identity` reports what the process actually loaded:

```python
import pynumerics

directory = pynumerics.initialize()
identity = pynumerics.numerics_identity()
print(directory)
print(identity.version, identity.informational_version)
print(identity.sha256, identity.is_packaged)   # is_packaged: the file is the pinned assembly
print(pynumerics.get_runtime_info().kind)      # CoreCLR
```

Packages built on Numerics, such as Hydrologics, check the assembly with
`pynumerics.assert_numerics_identity(version, sha256)` and load their own engine into the same
runtime with `pynumerics.add_assembly(path)`. References from an added assembly to `Numerics`
bind to the one already loaded; `pynumerics.loaded_assemblies()` lists what was added.

Import `pynumerics` before anything else that uses pythonnet. A bare `import clr` before the
loader runs starts the legacy .NET Framework on Windows, which cannot be undone in that process,
and the loader refuses to continue in that state with `pynumerics.RuntimeLoadError`.

## A first fit and quantile

```python
from pynumerics import distributions

sample = distributions.GeneralizedExtremeValue(100.0, 20.0, -0.1).generate_random_values(80, seed=12345)

gev = distributions.GeneralizedExtremeValue()
gev.estimate(sample, "method_of_linear_moments")
print(gev.parameters)                       # location, scale, shape (Hosking's sign)
print(gev.inverse_cdf(0.99))                # the 100-year quantile
print(gev.inverse_cdf([0.5, 0.9, 0.99]))    # many quantiles in one call
print(gev.cdf(150.0), gev.ccdf(150.0))
```

Method names are the members of `pynumerics.enums.EstimationMethod`, given as the enum or as a
snake-case string; every pynumerics enum accepts both. `distributions.FAMILIES` maps each Numerics
type name to its typed class, and `create_distribution(name, parameters)` builds any family by
name. The distributions page lists every family and its parameters.

## Arrays and marshalling

Wrappers take Python scalars, lists or NumPy arrays. A scalar argument returns a float and an
array argument returns a NumPy array. Arrays cross into .NET once, through the Numerics array
overload where there is one, and results come back with one bulk copy. Each copy goes through
the buffer protocol of the .NET array, which pins the array only for the duration of the copy
and keeps every bit. `pynumerics.marshal` holds the copies the wrappers use, for your own direct calls into
Numerics:

```python
import numpy as np
from pynumerics import marshal

values = np.array([1.5, -0.0, np.inf, np.nan])
net = marshal.to_net_doubles(values)             # System.Double[]
back = marshal.to_numpy_doubles(net)
print(net.Length, back, np.signbit(back[1]))     # every bit survives, including -0.0
```

Each wrapper holds its .NET object in the `.net` attribute, so a Numerics member that has no
wrapper is one attribute away.

## Errors

Every public wrapper raises Python exceptions, never .NET ones, and keeps the Numerics message
text exactly. The .NET exception is chained as `__cause__`.

| .NET exception | Python exception |
| --- | --- |
| `ArgumentException`, `ArgumentOutOfRangeException` | `ValueError` |
| `OverflowException`, `DivideByZeroException`, other `ArithmeticException` | `OverflowError`, `ZeroDivisionError`, `ArithmeticError` |
| `NotImplementedException`, `NotSupportedException` | `NotImplementedError` |
| any other .NET exception | `RuntimeError` |

```python
from pynumerics.distributions import Normal

try:
    Normal().inverse_cdf(1.5)
except ValueError as error:
    print(error)                            # the Numerics message
    print(type(error.__cause__).__name__)   # the .NET exception
```

Code that calls .NET members directly can use the `pynumerics.errors.net_errors()` context manager
or the `pynumerics.errors.translated` decorator for the same translation. Members that need a newer
Numerics assembly than the one loaded raise `NotImplementedError`.

## Value semantics and frozen distributions

Distributions compare by value: two objects are equal when they wrap the same Numerics family
with bitwise-equal parameters. A mutable distribution is not hashable, because `estimate` and
`set_parameters` change its value. `frozen()` returns a read-only copy that is hashable by the
same value, so it can be a dictionary key or a field of a frozen record. `clone()` gives a
mutable copy again.

```python
from pynumerics.distributions import GeneralizedExtremeValue

fit = GeneralizedExtremeValue(100.0, 20.0, -0.1)
frozen = fit.frozen()
print(frozen == fit, {frozen: "annual maxima"}[frozen])
try:
    frozen.set_parameters([1.0, 2.0, 0.0])
except TypeError as error:
    print(error)
thawed = frozen.clone()
thawed.set_parameters([110.0, 20.0, -0.1])
print(thawed.parameters, frozen.parameters)
```

## Callbacks

Members that take a function, such as integrands, objectives, likelihoods and root-finding
functions, accept any of three forms:

- a Python callable, called from .NET on every evaluation. Measured on a desktop, a raw pythonnet
  delegate costs about 2 microseconds per call; pynumerics adds about 0.2 microseconds for scalar
  arguments and about 0.4 microseconds for an array of up to 1,024 values (staged through a
  per-thread buffer), and one copy at memory speed for longer arrays (about 30 microseconds for
  100,000 values). The callable always receives a new float64 NumPy array with the exact bits
  .NET passed;
- C# source text, compiled once with the bundled Roslyn compiler, which then runs entirely in .NET;
- a .NET delegate of the right type.

```python
import math
from pynumerics import callbacks, integration

result = integration.integrate_adaptive_gauss_kronrod(lambda x: math.exp(-x * x), -3.0, 3.0)
print(result.value, result.status)

f = callbacks.scalar_function(math.sqrt)   # an explicit Func<double, double> calling Python
print(f(2.0))
```

`pynumerics.callbacks.compile_csharp(source, delegate_type="Func<double,double>", *, usings=())`
compiles a lambda, an expression or static methods. The code sees `System`, `System.Linq`,
`System.Collections.Generic` and the main Numerics namespaces, and a compile error raises
`ValueError` with Roslyn's diagnostics, with line and column counted in your source. The Roslyn
assemblies ship in the wheel, pinned by SHA-256 and checked when they load; nothing is downloaded
at run time. Compiled delegates are cached by source: the first compilation in a process takes a
few seconds while Roslyn starts, later ones much less. Anywhere a callable is taken, the C# source
string can be passed directly:

```python
from pynumerics import callbacks, integration, roots

g = callbacks.compile_csharp("x => Math.Exp(-x * x)")
print(g(0.5))
print(integration.integrate_adaptive_gauss_kronrod("x => Math.Exp(-x * x)", -3.0, 3.0).value)
print(roots.brent_solve("x => x * x - 2.0", 0.0, 2.0))
```

The parallel rule: a Python callback holds the global interpreter lock, so a Numerics object that
would call it from several threads gains nothing and can contend. Whenever any function given to a
Numerics object is a Python callable, the wrapper sets that object's parallel switch
(`ParallelizeChains`, `UseParallel` and similar) to false. Compiled C# and .NET delegates leave the
switch as you set it. `callbacks.uses_python_callback` tells which case applies.

```python
from pynumerics import callbacks
from pynumerics.optimization import DifferentialEvolution

def objective(x):
    return float(x[0] ** 2 + x[1] ** 2)

de = DifferentialEvolution(objective, [-5.0, -5.0], [5.0, 5.0], use_parallel=True)
print(callbacks.uses_python_callback(objective), de.use_parallel)   # True False
```

Use C# for functions evaluated many times, such as MCMC likelihoods and global optimization
objectives: it removes the per-call crossing and keeps the parallel paths available.

## Migrating from rmc-numerics

pynumerics replaces the `rmc-numerics` distribution. The change is a hard break: there is no
`rmc_numerics` alias module and no fallback for the old environment variable.

| rmc-numerics | pynumerics |
| --- | --- |
| `pip install rmc-numerics` | `pip install pynumerics` |
| `import rmc_numerics` | `import pynumerics` |
| `from rmc_numerics import distributions` | `from pynumerics import distributions` |
| `RMC_NUMERICS_DLL_DIR` | `PYNUMERICS_DLL_DIR` |
| `rmc-numerics==...` in requirements | `pynumerics==2.3.0` |

The submodules keep their names, so most code needs only the import change. Set
`PYNUMERICS_DLL_DIR` wherever `RMC_NUMERICS_DLL_DIR` was set, in shells, CI and service
definitions; the old variable is ignored. Uninstall `rmc-numerics` from environments that move,
so one process never loads two copies of the loader.
