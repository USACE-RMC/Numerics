# Integration, derivatives, root finding and ODEs

pynumerics wraps the Numerics calculus routines in four modules:

| Module | Numerics types |
| --- | --- |
| `pynumerics.integration` | `Integration`, `TrapezoidalRule`, `SimpsonsRule`, `AdaptiveSimpsonsRule`, `AdaptiveSimpsonsRule2D`, `AdaptiveGaussKronrod`, `AdaptiveGaussKronrod2D`, `AdaptiveGaussLobatto`, `MonteCarloIntegration`, `Miser`, `Vegas`, `IntegrationStatus` |
| `pynumerics.differentiation` | `NumericalDerivative` |
| `pynumerics.roots` | `Bisection`, `Brent` (`brent_solve`, `brent_bracket`), `Secant`, `NewtonRaphson` |
| `pynumerics.ode` | `RungeKutta` |

Every function argument accepts any of three forms (see `pynumerics.callbacks`):

- a Python callable, which crosses from .NET into Python on each evaluation (about 2 microseconds);
- a string of C# source, compiled once with Roslyn, which runs entirely in .NET;
- a .NET delegate of the right type.

A Python callable and the equivalent C# source give the same bits, so code can start in Python
and move hot integrands to C# without changing results.

## Integration

The fixed rules return a float:

```python
import math
from pynumerics import integration

f = lambda x: math.exp(-x * x)
print(integration.gauss_legendre(f, -3.0, 3.0))        # 10-point Gauss-Legendre
print(integration.gauss_legendre_20(f, -3.0, 3.0))     # 20-point Gauss-Legendre
print(integration.trapezoidal_rule(f, -3.0, 3.0, 1000))
print(integration.simpsons_rule(f, -3.0, 3.0, 1000))
print(integration.midpoint_rule(f, -3.0, 3.0, 1000))
```

Each Numerics `Integrator` has an `integrate_*` function returning an `IntegrationResult` with
`value`, `standard_error` (NaN when the method has none), `iterations`, `function_evaluations` and
`status` (an `IntegrationStatus`). The `Integrator` settings are keyword arguments, listed in
`integration.INTEGRATOR_SETTINGS`; anything left out keeps the Numerics default.

```python
import math
from pynumerics import integration
from pynumerics.integration import IntegrationStatus

result = integration.integrate_adaptive_gauss_kronrod(
    lambda x: math.exp(-x * x), -3.0, 3.0, absolute_tolerance=1e-10, relative_tolerance=1e-10)
print(result.value, result.standard_error, result.status)
assert result.status is IntegrationStatus.SUCCESS

# The same integral with a compiled C# integrand: no Python call per evaluation.
compiled = integration.integrate_adaptive_gauss_kronrod("x => Math.Exp(-x * x)", -3.0, 3.0,
                                                         absolute_tolerance=1e-10, relative_tolerance=1e-10)
assert compiled == result

print(integration.integrate_adaptive_simpson(math.sin, 0.0, math.pi).value)
print(integration.integrate_adaptive_gauss_lobatto(math.sin, 0.0, math.pi).value)
print(integration.integrate_trapezoidal(math.sin, 0.0, math.pi).value)
print(integration.integrate_simpson(math.sin, 0.0, math.pi).value)
```

The Gauss-Kronrod recorder reports every node of the accepted intervals, and both adaptive 1-D
rules can integrate over stratification bins:

```python
import math
from pynumerics import integration, sampling

nodes = []
result = integration.integrate_adaptive_gauss_kronrod(
    lambda x: math.exp(-x * x), -3.0, 3.0, recorder=lambda x, w, fx: nodes.append((x, w, fx)))
print(len(nodes), sum(w * fx for _, w, fx in nodes), result.value)

bins = sampling.stratify_x_values(-3.0, 3.0, 4)
print(integration.integrate_adaptive_gauss_kronrod(lambda x: math.exp(-x * x), -3.0, 3.0, bins=bins).value)
```

`integrate_adaptive_gauss_kronrod` takes `method=`, an `AdaptiveGaussKronrodMethod` or its name. The
default is local bisection. `method="qags"` runs QUADPACK QAGS, the
algorithm behind `scipy.integrate.quad`, with the same evaluation counts and errors (results agree to
rounding, since the Kronrod sum order differs). Choose it for endpoint
singularities: `1 / sqrt(x)` on [0, 1] takes 231 evaluations and is accurate to about 1e-15, where local
bisection takes 1911 and leaves an error of 5.5e-9. Keep the default for cheap smooth integrands, where QAGS
is about 1.4 times slower, and for recorder sampling that must reproduce earlier nodes.

```python
import math
from pynumerics import integration

qags = integration.integrate_adaptive_gauss_kronrod(lambda x: 1.0 / math.sqrt(x), 0.0, 1.0, method="qags")
print(qags.value, qags.function_evaluations)  # 1.9999999999999993 231
```

Two-dimensional integrands take `(x, y)`:

```python
import math
from pynumerics import integration

f = lambda x, y: math.exp(-x * x - y * y)
print(integration.integrate_adaptive_simpson_2d(f, -2.0, 2.0, -2.0, 2.0).value)
print(integration.integrate_adaptive_gauss_kronrod_2d(f, -2.0, 2.0, -2.0, 2.0).value)
```

The sampling integrators take box bounds; the integrand receives a NumPy array. `seed` fixes the
Mersenne Twister. Vegas integrands also receive the sample weight, and Vegas returns a
`VegasResult` with the chi-squared and the adapted grid:

```python
import numpy as np
from pynumerics import integration

lower, upper = [0.0, 0.0, 0.0], [1.0, 1.0, 1.0]
square = lambda x: float(np.sum(x * x))   # the exact integral is 1
print(integration.integrate_monte_carlo(square, lower, upper, seed=7, max_iterations=40000).value)
print(integration.integrate_miser(square, lower, upper, max_function_evaluations=40000).value)

vegas = integration.integrate_vegas(lambda x, w: float(np.sum(x * x)), lower, upper)
print(vegas.value, vegas.standard_error, vegas.chi_squared, vegas.grid.shape)

# A rare-event probability: P(U1 + U2 > 1.98) = 0.0002.
tail = integration.integrate_vegas("(x, w) => x[0] + x[1] > 1.98 ? 1.0 : 0.0", [0.0, 0.0], [1.0, 1.0],
                                   rare_event_probability=2e-4)
print(tail.value)
```

None of these integrators has a parallel switch; every `integrate_*` function still applies
`callbacks.apply_parallel_policy`, so a Python callback never runs on several threads.

## Derivatives

```python
import math
import numpy as np
from pynumerics import differentiation

print(differentiation.derivative(math.exp, 1.0))             # central difference
print(differentiation.forward_difference(math.exp, 1.0))
print(differentiation.backward_difference(math.exp, 1.0))
print(differentiation.central_difference(math.exp, 1.0, step_size=1e-4))
print(differentiation.ridders(math.exp, 1.0))  # (derivative, error estimate); default step 0.1 * (1 + |x|)
# The default step 0.105 leaves the domain of the logarithm, where f returns NaN, so it is halved: about 20
print(differentiation.ridders(lambda x: math.log(x) if x > 0 else math.nan, 0.05))
print(differentiation.second_derivative(math.exp, 1.0))
print(differentiation.calculate_step_size(1.0))

f = lambda t: float(t[0] ** 2 * t[1])
print(differentiation.gradient(f, [1.0, 2.0]))
print(differentiation.hessian(f, [1.0, 2.0]))
print(differentiation.jacobian(lambda t: np.array([t[0] * t[1], t[0] ** 2]), [1.0, 2.0]))
print(differentiation.jacobian_at(lambda x, t: t[0] * math.exp(t[1] * x), [0.0, 1.0, 2.0], [1.5, 0.5]))
```

`gradient`, `hessian` and `jacobian` accept `lower_bounds` and `upper_bounds` and switch to
one-sided differences near a bound.

## Root finding

```python
import numpy as np
from pynumerics import roots

f = lambda x: x * x - 2.0
print(roots.bisection_solve(f, 1.0, 0.0, 4.0))
print(roots.secant_solve(f, 0.0, 4.0))
print(roots.newton_raphson_solve(f, lambda x: 2.0 * x, 1.0))
print(roots.newton_raphson_robust_solve(f, lambda x: 2.0 * x, 1.0, 0.0, 4.0))
print(roots.brent_solve(f, 0.0, 4.0))

# A system: x^2 + y^2 = 4 and x = y.
solution = roots.newton_raphson_solve_system(
    lambda x: np.array([x[0] ** 2 + x[1] ** 2 - 4.0, x[0] - x[1]]),
    lambda x: np.array([[2.0 * x[0], 2.0 * x[1]], [1.0, -1.0]]),
    [1.0, 1.0])
print(solution)
```

A solver that fails raises `ValueError` with the Numerics message; pass `report_failure=False` to
get the last iterate instead.

`brent_solve` is a literal port of SciPy's `brentq.c`: `brent_solve(f, a, b, xtol, maxiter,
relative_tolerance=rtol)` returns the same float as `scipy.optimize.brentq(f, a, b, xtol=xtol,
rtol=rtol, maxiter=maxiter)`. The default `relative_tolerance` is
`roots.BRENT_MINIMUM_RELATIVE_TOLERANCE` (4 * 2**-52), SciPy's default and smallest `rtol`.

```python
from pynumerics import roots

print(roots.brent_solve(lambda x: x * x - 2.0, 0.0, 2.0, 1e-12, relative_tolerance=1e-10))  # 1.4142135623731364
```

## ODE solvers

The right-hand side is `f(t, y)` for one state variable:

```python
import math
from pynumerics import ode

f = lambda t, y: -2.0 * t * y   # y = exp(-t^2)
print(ode.runge_kutta_fourth_order(f, 1.0, 0.0, 2.0, 11))     # values on 11 even times
print(ode.runge_kutta_second_order(f, 1.0, 0.0, 2.0, 11))
print(ode.runge_kutta_fourth_order_step(f, 1.0, 0.0, 0.1))     # one step
print(ode.runge_kutta_fehlberg(f, 1.0, 0.0, 2.0, 1e-6, 1e-10), math.exp(-4.0))
print(ode.runge_kutta_cash_karp("(t, y) => -2.0 * t * y", 1.0, 0.0, 2.0, 1e-6, 1e-10))
```
