# Multivariate distributions and copulas

pynumerics wraps the multivariate distributions and the bivariate copulas of
`Numerics.Distributions` in two modules:

| Module | Numerics types |
| --- | --- |
| `pynumerics.multivariate` | `MultivariateDistribution`, `MultivariateNormal`, `MultivariateStudentT`, `Multinomial`, `BivariateEmpirical`, `Dirichlet`, `MultivariateDistributionType` |
| `pynumerics.copulas` | `BivariateCopula`, `ArchimedeanCopula`, the Normal, Student-t, Clayton, Gumbel, Frank, Joe, Ali-Mikhail-Haq and independence copulas, `BivariateCopulaEstimation`, `CopulaFactory`, `CopulaType`, `CopulaEstimationMethod` |

Every method calls one Numerics member; nothing is computed in Python. Numerics evaluates a
multivariate density or CDF at one point, so a two-dimensional array of points is evaluated row by
row. Samples come back as one `(n, d)` NumPy array in a single block copy. Copula members take
paired `u` and `v` arrays, which cross into Numerics once where Numerics has an array overload.

Arguments that take a random generator accept a `pynumerics.sampling.MersenneTwister` (or any
.NET `System.Random`); `generate_random_values` accepts either an integer seed or a generator.

## Multivariate normal

```python
import numpy as np
from pynumerics.multivariate import MultivariateNormal

mvn = MultivariateNormal([0.5, -1.0, 2.0], [[1.0, 0.5, 0.2], [0.5, 2.0, 0.3], [0.2, 0.3, 1.5]])
print(mvn.pdf([0.0, 0.0, 0.0]), mvn.log_pdf([0.0, 0.0, 0.0]))
print(mvn.pdf(np.array([[0.0, 0.0, 0.0], [1.0, -1.0, 3.0]])))   # one value per row
print(mvn.cdf([1.0, 0.0, 3.0]))                                  # Genz randomized lattice rule
print(mvn.interval([-1.0, -2.0, 1.0], [1.0, 0.0, 3.0]))          # P(lower <= X <= upper)
print(mvn.mahalanobis([1.0, 0.0, 3.0]))

conditional = mvn.conditional([1], [0.5])      # X0, X2 given X1 = 0.5
print(conditional.mean, conditional.covariance)
print(mvn.marginal(2, 0).covariance, mvn.independent_normal(1))

sample = mvn.generate_random_values(1000, seed=12)
print(sample.shape, mvn.latin_hypercube_random_values(10, 5).shape)
```

Above two dimensions the CDF is a randomized lattice rule whose shifts come from `mvnuni`. The
default generator is seeded, so results reproduce; assign a seed or generator to tie them to your
own seed, and tighten `absolute_error`, `relative_error` and `max_evaluations` for more accuracy.

```python
from pynumerics.multivariate import MultivariateNormal

mvn = MultivariateNormal([0.0, 0.0, 0.0])
mvn.mvnuni = 2024
mvn.absolute_error, mvn.max_evaluations = 1e-7, 100000
print(mvn.cdf([0.5, 0.2, 1.0]))
print(MultivariateNormal.bivariate_cdf(-1.0, 0.5, 0.4))   # Genz BVN upper orthant P(X > -1, Y > 0.5)
print(MultivariateNormal.bivariate(0.0, 1.0, 1.0, 2.0, 0.5).covariance)
```

## Multivariate Student-t

`MultivariateStudentT(degrees_of_freedom, location, scale_matrix)` takes the scale matrix, not the
covariance; the covariance is ν/(ν − 2)·Σ for ν > 2. The CDF is the Genz-Bretz method: in two
dimensions it is deterministic and exact to about 1e-14 (the Dunnett-Sobel closed form for integer ν,
adaptive quadrature over the χ(ν) mixing variable otherwise); above two dimensions it is a randomized
lattice rule controlled by `absolute_error` (default 1e-5), `relative_error` (default 0) and
`max_evaluations` (default 25,000), with shifts drawn from `mvnuni`.

```python
from pynumerics.multivariate import MultivariateStudentT

t = MultivariateStudentT(5.0, [1.0, -1.0], [[2.0, 0.6], [0.6, 1.0]])
print(t.pdf([0.5, 0.0]), t.cdf([0.5, 0.0]), t.ccdf([0.5, 0.0]))
print(t.covariance, t.standard_deviation)
print(t.inverse_cdf([0.5, 0.5, 0.9]))      # d normal probabilities, then the chi-square one
print(t.generate_random_values(5, seed=3))
```

## Multinomial and categorical draws

```python
import numpy as np
from pynumerics.multivariate import Multinomial
from pynumerics.sampling import MersenneTwister

m = Multinomial(6, [0.2, 0.3, 0.5])
print(m.pdf([1, 2, 3]), m.log_pmf([1, 2, 3]), m.mean, m.covariance(0, 1))
print(m.generate_random_values(4, seed=7))

cumulative = np.cumsum([0.1, 0.0, 0.4, 0.5])
uniforms = np.random.default_rng(1).random(10)
print(Multinomial.categorical_indices(cumulative, uniforms))   # clip(searchsorted(c, u, "right"), 0, n - 1)
print(Multinomial.categorical_indices(cumulative, [0.1, 0.5], inclusive=True))   # [0 2]: searchsorted "left"
print(Multinomial.sample([1.0, 2.0, 3.0], MersenneTwister(11)))
```

`categorical_indices` maps the whole array in one Numerics call. `inclusive=True` selects
min{k : c[k] ≥ u} instead of min{k : c[k] > u}, the inverse-transform count of a discrete CDF
table (for example a Poisson CDF from `scipy.stats.poisson.cdf`); the two differ only for a
variate equal to a table entry, and the default is unchanged.

## Bivariate empirical and Dirichlet

```python
from pynumerics.multivariate import BivariateEmpirical, Dirichlet

table = BivariateEmpirical([0.0, 1.0, 2.0], [0.0, 10.0], [[0.01, 0.1], [0.2, 0.5], [0.3, 0.99]])
print(table.cdf(0.5, 5.0), table.cdf([0.5, 1.5], [5.0, 2.0]))   # bilinear interpolation
table.probability_transform = "normal_z"            # interpolate in standard normal space
print(table.cdf(0.5, 5.0))

dirichlet = Dirichlet([2.0, 3.0, 5.0])
print(dirichlet.mean, dirichlet.covariance_matrix())
print(dirichlet.generate_random_values(3, seed=4))
```

## Copulas

Copulas take uniform margins. Attach marginal distributions to sample in real units and to use the
inference-from-margins and full likelihoods.

```python
import numpy as np
from pynumerics import copulas
from pynumerics.distributions import Gumbel, Normal

clayton = copulas.ClaytonCopula(2.0)
u, v = np.array([0.1, 0.5, 0.9]), np.array([0.2, 0.5, 0.7])
print(clayton.cdf(u, v), clayton.pdf(u, v), clayton.conditional_cdf(u, v))
print(clayton.and_joint_exceedance_probability(0.99, 0.99))
print(clayton.lower_tail_dependence, clayton.generator(0.3))

gumbel = copulas.GumbelCopula(1.8, Normal(2.0, 0.5), Gumbel(10.0, 3.0))
print(gumbel.generate_random_values(5, seed=3))       # marginal units
print(copulas.create_copula("student_t").copula_parameters)
```

### Estimation

`estimate` runs `BivariateCopulaEstimation.Estimate` in place. The pseudo likelihood takes
plotting positions; inference from margins and the full likelihood take the data and the attached
marginals (the full likelihood refits them too).

```python
import numpy as np
from pynumerics import copulas
from pynumerics.distributions import Normal

rng = np.random.default_rng(5)
z = rng.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], size=300)
x, y = 10.0 + 3.0 * z[:, 0], -2.0 + 0.5 * z[:, 1]
n = len(x)
u = (np.argsort(np.argsort(x)) + 1.0) / (n + 1.0)
v = (np.argsort(np.argsort(y)) + 1.0) / (n + 1.0)

normal = copulas.NormalCopula()
normal.estimate(u, v, "pseudo_likelihood")
print(normal.theta)

full = copulas.NormalCopula(0.0, Normal(), Normal())
full.estimate(x, y, copulas.CopulaEstimationMethod.FULL_LIKELIHOOD)
print(full.theta, full.marginal_distribution_x, full.marginal_distribution_y)

clayton = copulas.ClaytonCopula()
clayton.set_theta_from_tau(x, y)                  # invert Kendall's tau
print(clayton.theta)
```

### XML round trip

```python
from pynumerics import copulas
from pynumerics.distributions import Normal

copula = copulas.StudentTCopula(0.4, 6.0, Normal(0.0, 1.0), Normal(5.0, 2.0))
text = copula.to_xml(include_marginals=True)
restored = copulas.copula_from_xml(text, include_marginals=True)
print(text)
print(type(restored).__name__, restored.copula_parameters, restored.marginal_distribution_y)
```

Without `include_marginals` the element carries only the family and its parameters.
