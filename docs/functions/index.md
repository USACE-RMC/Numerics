# Univariate Functions

[← Previous: Multivariate Distributions](../distributions/multivariate.md) | [Back to Index](../index.md) | [Next: Machine Learning →](../machine-learning/machine-learning.md)

The `Numerics.Functions` namespace provides the uncertain-function toolkit: univariate
functional forms with optional stochastic residuals, sampled by confidence level, plus
serialization, a factory, and posterior-ensemble sampling. These are the shared math behind
risk input functions (rating curves, transforms, damage functions) in the consuming
engineering applications.

## The `IUnivariateFunction` contract

Every function implements `IUnivariateFunction`:

| Member | Meaning |
|--------|---------|
| `Function(x)` / `InverseFunction(y)` | Forward and inverse evaluation |
| `SetParameters(values)` / `ValidateParameters(...)` | The parameter-vector surface (fitted parameter sets apply directly) |
| `Minimum` / `Maximum` | The evaluation support (clamped); some types derive `Minimum` from a location parameter |
| `IsDeterministic` | Whether the stochastic residual is active |
| `ConfidenceLevel` | The runtime draw in [0, 1] driving the residual; values outside [0, 1] evaluate the mean/median form. Runtime state — never serialized |

## Available functions

| Function | Class | Form | Residual |
|----------|-------|------|----------|
| **Linear** | `LinearFunction` | $Y = \alpha + \beta X$ | Additive Gaussian, $\epsilon \sim N(0, \sigma)$ |
| **Power** | `PowerFunction` | $Y = \alpha (X - \xi)^\beta$, optional inverse form | Log-space Gaussian (multiplicative) |
| **Tabular** | `TabularFunction` | Interpolation over `UncertainOrderedPairedData` with axis transforms | Co-monotonic percentile of the ordinate distributions |
| **Segmented power** | `SegmentedPowerFunction` | BaRatin addition mode: $Q(h) = \sum_k 10^{\log_{10}\alpha_k} (h - h_k)^{\beta_k} \cdot \mathbb{1}\{h > h_k\}$ | Log₁₀-space Gaussian |
| **Composite** | `CompositeFunction` | Weighted average $\sum w_i f_i(x)$, or a mixture selected by one uniform | Inherited from the children |

### Segmented power (BaRatin addition mode)

`SegmentedPowerFunction` uses a parameter vector that can be applied directly to fitted draws. The
parameter vector is `[h₁, log₁₀α₁, β₁, h₂, log₁₀α₂, β₂, …, σ]` (length `3·segments + 1`), so
a fitted posterior `ParameterSet.Values` applies directly through `SetParameters`. Breakpoints
must be strictly ordered, exponents positive (the monotonicity constraint behind the
numeric Brent inverse), discharge is zero at and below the cease-to-flow stage `h₁`, and one
segment degenerates deterministically to the plain `PowerFunction` (the stochastic residual
spaces differ: log₁₀ here versus natural log in `PowerFunction`).

```cs
using Numerics.Functions;

// Two controls: main channel from stage 1, overbank activating at stage 3.
var rating = new SegmentedPowerFunction(new[] { 1.0, 1.5, 2.0, 3.0, 1.2, 1.5, 0.1 });
double q = rating.Function(5.0);            // deterministic (median) curve: ConfidenceLevel
                                            // defaults to -1, and any value outside [0, 1]
                                            // selects the deterministic curve
rating.ConfidenceLevel = 0.75;              // 75th-percentile curve via the log₁₀-space
                                            // residual: multiplies by 10^(z·σ)
double q75 = rating.Function(5.0);
```

### Composite (weighted average and mixture)

`CompositeFunction` combines child functions with non-negative weights summing to one. The
**weighted-average** mode evaluates $\sum w_i f_i(x)$ with the composite's confidence level
driving every child co-monotonically. The **mixture** mode composes with a single uniform: the
draw selects a child by cumulative weight and re-scales the remainder as the child's own draw
— deterministic composition sampling with no internal random source. Outside [0, 1] both modes
evaluate the weighted average of each child evaluated at its own configured confidence level.

## Serialization and the factory

Every concrete function serializes with an instance `ToXElement()` and a static
`FromXElement(XElement)`. `UnivariateFunctionFactory` dispatches on the
`UnivariateFunctionType` enum or on a serialized element's local name, and maps live instances
back to their enum type:

```cs
var element = rating.ToXElement();
IUnivariateFunction restored = UnivariateFunctionFactory.CreateFromXElement(element);
UnivariateFunctionType kind = UnivariateFunctionFactory.GetFunctionType(restored); // SegmentedPower
```

`ConfidenceLevel` is runtime sampling state and never serializes; `PowerFunction.Minimum`
derives from ξ and never serializes. The serialization surface lives on the concrete classes
(not the interface) so external `IUnivariateFunction` implementations remain source-compatible.

## Posterior ensembles

`EnsembleFunction` carries a template function plus an array of posterior `ParameterSet`
draws for use inside a simulation engine. Both sampling surfaces return a fresh
clone of the template configured with the selected draw, so concurrent realizations share no
mutable state.

```cs
var ensemble = new EnsembleFunction(rating, posteriorParameterSets);
IUnivariateFunction draw = ensemble.SampleAt(index);      // posterior draw by index
IUnivariateFunction pDraw = ensemble.Sample(0.37);        // min(⌊u·N⌋, N − 1)
```

## Link functions

The sibling `Link Functions` family (`ILinkFunction`: identity, log, logit, probit,
complementary log-log, Fisher-z, Yeo-Johnson) serves regression and machine-learning
transformations and has its own factory, `LinkFunctionFactory`.

---

[← Previous: Multivariate Distributions](../distributions/multivariate.md) | [Back to Index](../index.md) | [Next: Machine Learning →](../machine-learning/machine-learning.md)
