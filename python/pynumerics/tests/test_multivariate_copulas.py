"""Multivariate distributions and bivariate copulas against SciPy and textbook values.

Each wrapper calls one Numerics member; the oracle computes the same quantity independently
(SciPy distributions, direct quadrature of the density, or closed forms from Nelsen (2006),
An Introduction to Copulas, 2nd ed., Springer), so the tolerances cover only the different
operation order or, where Numerics integrates numerically, its documented accuracy.
"""

from __future__ import annotations

import time

import numpy as np
import pytest

from pynumerics import copulas, distributions, multivariate, sampling
from pynumerics.marshal import to_net_doubles

stats = pytest.importorskip("scipy.stats")
integrate = pytest.importorskip("scipy.integrate")
optimize = pytest.importorskip("scipy.optimize")

SCALE = [[2.0, 0.6], [0.6, 1.0]]
LOCATION = [1.0, -1.0]


def _plotting_positions(values: np.ndarray) -> np.ndarray:
    return stats.rankdata(values) / (values.size + 1.0)


def _gaussian_pair(n: int, rho: float, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    z = rng.multivariate_normal([0.0, 0.0], [[1.0, rho], [rho, 1.0]], size=n)
    return z[:, 0], z[:, 1]


# ---------------------------------------------------------------------------------------------
# MultivariateStudentT


def test_student_t_density_moments_and_distance_equal_scipy() -> None:
    # Arrange
    t = multivariate.MultivariateStudentT(5.0, LOCATION, SCALE)
    oracle = stats.multivariate_t(LOCATION, SCALE, df=5.0)
    points = np.array([[0.5, 0.0], [1.0, -1.0], [3.0, 1.0], [-4.0, 2.5]])
    # Act / Assert
    np.testing.assert_allclose(t.pdf(points), oracle.pdf(points), rtol=1e-13)
    np.testing.assert_allclose(t.log_pdf(points), oracle.logpdf(points), rtol=1e-13)
    assert t.pdf(points[0]) == t.pdf(points)[0]
    np.testing.assert_allclose(t.covariance, 5.0 / 3.0 * np.array(SCALE), rtol=1e-14)
    np.testing.assert_allclose(t.variance, 5.0 / 3.0 * np.diag(SCALE), rtol=1e-14)
    np.testing.assert_allclose(t.standard_deviation, np.sqrt(5.0 / 3.0 * np.diag(SCALE)), rtol=1e-14)
    for vector in (t.mean, t.median, t.mode, t.location):
        np.testing.assert_array_equal(vector, LOCATION)
    np.testing.assert_array_equal(t.scale_matrix, SCALE)
    diff = np.array([0.0, 0.0]) - LOCATION
    assert t.mahalanobis([0.0, 0.0]) == pytest.approx(diff @ np.linalg.solve(SCALE, diff), rel=1e-14)
    assert t.dimension == 2 and t.degrees_of_freedom == 5.0 and t.is_positive_definite and t.parameters_valid
    assert t.type is multivariate.MultivariateDistributionType.MULTIVARIATE_STUDENT_T
    assert t.display_name and t.short_display_name


def _bivariate_t_cdf(point: list[float], df: float) -> float:
    """Exact bivariate t CDF: the conditional univariate t of the second coordinate integrated over the first."""
    s1, s2 = np.sqrt(SCALE[0][0]), np.sqrt(SCALE[1][1])
    rho = SCALE[0][1] / (s1 * s2)
    a, b = (point[0] - LOCATION[0]) / s1, (point[1] - LOCATION[1]) / s2

    def integrand(z: float) -> float:
        spread = np.sqrt((1.0 - rho**2) * (df + z * z) / (df + 1.0))
        return stats.t.pdf(z, df) * stats.t.cdf((b - rho * z) / spread, df + 1.0)

    return integrate.quad(integrand, -np.inf, a, epsabs=1e-15, epsrel=1e-13, limit=500)[0]


@pytest.mark.parametrize("df", [3.0, 5.5, 30.0])
def test_student_t_cdf_equals_exact_quadrature_and_scipy(df: float) -> None:
    # Arrange
    t = multivariate.MultivariateStudentT(df, LOCATION, SCALE)
    oracle = stats.multivariate_t(LOCATION, SCALE, df=df)
    for point in ([0.5, 0.0], [1.0, -1.0], [3.0, 1.0]):
        exact = _bivariate_t_cdf(point, df)
        # Act
        value = t.cdf(point)
        # Assert: the 2-D Genz-Bretz CDF is deterministic and exact to rounding (Dunnett-Sobel for
        # integer df, adaptive quadrature otherwise); scipy's randomized QMC is accurate to about 1e-6
        assert value == pytest.approx(exact, abs=1e-12)
        assert value == pytest.approx(oracle.cdf(point, maxpts=200000, random_state=1), abs=2e-6)
        assert t.ccdf(point) == pytest.approx(1.0 - value, abs=1e-15)
        assert t.log_cdf(point) == pytest.approx(np.log(value), rel=1e-14)
        assert t.log_ccdf(point) == pytest.approx(np.log1p(-value), rel=1e-13)


def test_student_t_three_dimensional_cdf_follows_its_seeded_generator() -> None:
    # Arrange
    def make() -> multivariate.MultivariateStudentT:
        t = multivariate.MultivariateStudentT(4.0, 3)
        t.mvnuni = 7
        return t

    oracle = stats.multivariate_t(np.zeros(3), np.eye(3), df=4.0)
    point = [0.5, 0.2, 1.0]
    # R 4.4.3 mvtnorm pmvt(upper, df = 4, sigma = diag(3), GenzBretz(maxpts = 5e7, abseps = 1e-11)), error 2.8e-9
    reference = 0.322384846867854
    tight = make()
    tight.absolute_error, tight.max_evaluations = 1e-7, 2_000_000
    # Act
    first, second = make().cdf(point), make().cdf(point)
    # Assert
    assert first == second
    assert first == pytest.approx(reference, abs=2e-5)
    assert first == pytest.approx(oracle.cdf(point, maxpts=200000, random_state=1), abs=2e-5)
    assert tight.cdf(point) == pytest.approx(reference, abs=1e-6)
    assert isinstance(make().mvnuni, sampling.MersenneTwister)


def test_student_t_cdf_error_control_members() -> None:
    # Arrange
    t = multivariate.MultivariateStudentT(4.0, 3)
    # Act
    defaults = (t.absolute_error, t.relative_error, t.max_evaluations)
    t.absolute_error, t.relative_error, t.max_evaluations = 1e-7, 1e-6, 123456
    clone = t.clone()
    # Assert
    assert defaults == (1e-5, 0.0, 25000)
    assert (clone.absolute_error, clone.relative_error, clone.max_evaluations) == (1e-7, 1e-6, 123456)
    with pytest.raises(ValueError):
        t.absolute_error = -1.0
    with pytest.raises(ValueError):
        t.max_evaluations = 0


def test_student_t_sampling_members() -> None:
    # Arrange
    t = multivariate.MultivariateStudentT(6.0, LOCATION, SCALE)
    # Act
    sample = t.generate_random_values(20000, seed=3)
    from_generator = t.generate_random_values(50, sampling.MersenneTwister(3))
    lhs = t.latin_hypercube_random_values(200, 5)
    bins = sampling.stratify_probabilities(0.01, 0.99, 25)
    stratified = t.stratified_random_values(bins, 9)
    point = t.inverse_cdf([0.5, 0.5, 0.5])
    # Assert
    np.testing.assert_array_equal(from_generator, sample[:50])
    np.testing.assert_allclose(sample.mean(axis=0), LOCATION, atol=0.05)
    np.testing.assert_allclose(np.cov(sample.T), 1.5 * np.array(SCALE), rtol=0.1)
    assert lhs.shape == (200, 2) and stratified.shape == (25, 2)
    np.testing.assert_allclose(point, LOCATION, atol=1e-12)
    with pytest.raises(ValueError, match="Dimension"):
        t.inverse_cdf([0.5, 0.5])


def test_student_t_parameters_and_construction() -> None:
    # Arrange
    t = multivariate.MultivariateStudentT(3.0, [0.0, 0.0])
    # Act
    t.set_parameters(7.0, LOCATION, SCALE)
    message = t.validate_parameters(-1.0, LOCATION, SCALE)
    clone = t.clone()
    # Assert
    assert t.degrees_of_freedom == 7.0 and t.validate_parameters(7.0, LOCATION, SCALE) is None
    assert message is not None
    assert isinstance(clone, multivariate.MultivariateStudentT) and clone.net is not t.net
    np.testing.assert_array_equal(clone.scale_matrix, SCALE)
    assert multivariate.MultivariateStudentT(2.0, 4).dimension == 4
    with pytest.raises(ValueError):
        multivariate.MultivariateStudentT(0.0, LOCATION)
    with pytest.raises(ValueError, match="location vector"):
        multivariate.MultivariateStudentT(2.0, 2, SCALE)


# ---------------------------------------------------------------------------------------------
# Multinomial


def test_multinomial_mass_and_moments_equal_scipy() -> None:
    # Arrange
    p = [0.2, 0.3, 0.5]
    m = multivariate.Multinomial(6, p)
    oracle = stats.multinomial(6, p)
    counts = np.array([[1, 2, 3], [6, 0, 0], [0, 0, 6], [2, 2, 2], [1, 1, 1]], dtype=float)
    # Act / Assert
    np.testing.assert_allclose(m.pdf(counts), oracle.pmf(counts), rtol=1e-13, atol=1e-300)
    np.testing.assert_allclose(m.log_pmf(counts[:4]), oracle.logpmf(counts[:4]), rtol=1e-13)
    np.testing.assert_allclose(m.log_pdf(counts[:4]), oracle.logpmf(counts[:4]), rtol=1e-13)
    assert m.pdf([1, 1, 1]) == 0.0 and m.log_pmf([1, 1, 1]) == -np.inf
    np.testing.assert_allclose(m.mean, oracle.mean(), rtol=1e-15)
    np.testing.assert_allclose(m.variance, np.diag(oracle.cov()), rtol=1e-15)
    cov = np.array([[m.covariance(i, j) for j in range(3)] for i in range(3)])
    np.testing.assert_allclose(cov, oracle.cov(), rtol=1e-14)
    assert m.number_of_trials == 6 and m.dimension == 3
    np.testing.assert_array_equal(m.probabilities, p)
    assert m.type is multivariate.MultivariateDistributionType.MULTINOMIAL
    with pytest.raises(NotImplementedError, match="closed-form"):
        m.cdf([1, 2, 3])


def test_multinomial_sampling_and_invalid_inputs() -> None:
    # Arrange
    m = multivariate.Multinomial(10, [0.1, 0.6, 0.3])
    # Act
    sample = m.generate_random_values(20000, seed=4)
    # Assert
    np.testing.assert_array_equal(sample.sum(axis=1), 10.0)
    np.testing.assert_allclose(sample.mean(axis=0), [1.0, 6.0, 3.0], atol=0.05)
    np.testing.assert_array_equal(sample, m.generate_random_values(20000, seed=4))
    for bad in ([0.5], [0.6, 0.6], [-0.1, 1.1]):
        with pytest.raises(ValueError):
            multivariate.Multinomial(3, bad)
    with pytest.raises(ValueError):
        multivariate.Multinomial(0, [0.5, 0.5])


def test_multinomial_sample_uses_the_supplied_generator() -> None:
    # Arrange
    weights = [1.0, 2.0, 3.0]
    draws = sampling.MersenneTwister(11).next_doubles(200) * 6.0
    expected = np.minimum(np.searchsorted(np.cumsum(weights), draws, side="left"), 2)
    generator = sampling.MersenneTwister(11)
    # Act
    indices = [multivariate.Multinomial.sample(weights, generator) for _ in range(200)]
    # Assert: a draw on a boundary selects the lower category
    np.testing.assert_array_equal(indices, expected)
    with pytest.raises(ValueError, match="positive"):
        multivariate.Multinomial.sample([0.0, 0.0], generator)


def test_categorical_indices_equal_clipped_right_searchsorted() -> None:
    # Arrange
    rng = np.random.default_rng(5)
    cumulative = np.cumsum([0.1, 0.0, 0.25, 0.4, 0.0, 0.25])
    uniforms = np.concatenate([rng.random(5000), cumulative, [-1.0, 0.0, 1.0, 2.0, np.inf, -np.inf]])
    # Act
    indices = multivariate.Multinomial.categorical_indices(cumulative, uniforms)
    grid = multivariate.Multinomial.categorical_indices(cumulative, uniforms[:12].reshape(3, 4))
    # Assert
    expected = np.clip(np.searchsorted(cumulative, uniforms, side="right"), 0, cumulative.size - 1)
    np.testing.assert_array_equal(indices, expected)
    assert indices.dtype == np.int32
    np.testing.assert_array_equal(grid, expected[:12].reshape(3, 4))
    assert multivariate.Multinomial.categorical_indices(cumulative, 0.5) == np.searchsorted(cumulative, 0.5, "right")
    assert multivariate.Multinomial.categorical_indices([2.0], [0.0, 5.0]).tolist() == [0, 0]


def test_categorical_indices_inclusive_equal_clipped_left_searchsorted() -> None:
    # Arrange: a Poisson(3.5) CDF table, whose exact entries are the ties of the two forms
    stats = pytest.importorskip("scipy.stats")
    cdf = stats.poisson(3.5).cdf(np.arange(40))
    rng = np.random.default_rng(8)
    uniforms = np.concatenate([rng.random(5000), cdf, [0.0, 1 - 2.0**-53, 1.0, -1.0, np.inf]])
    # Act
    inclusive = multivariate.Multinomial.categorical_indices(cdf, uniforms, inclusive=True)
    strict = multivariate.Multinomial.categorical_indices(cdf, uniforms, inclusive=False)
    # Assert
    np.testing.assert_array_equal(inclusive, np.clip(np.searchsorted(cdf, uniforms, side="left"), 0, cdf.size - 1))
    np.testing.assert_array_equal(strict, multivariate.Multinomial.categorical_indices(cdf, uniforms))
    np.testing.assert_array_equal(strict, np.clip(np.searchsorted(cdf, uniforms, side="right"), 0, cdf.size - 1))
    ties = multivariate.Multinomial.categorical_indices(cdf[:10], cdf[:9], inclusive=True)
    np.testing.assert_array_equal(ties, np.arange(9))
    assert multivariate.Multinomial.categorical_indices([0.2, 0.5, 1.0], 0.2, inclusive=True) == 0
    with pytest.raises(ValueError):
        multivariate.Multinomial.categorical_indices([0.5, 1.0], [np.nan], inclusive=True)


@pytest.mark.parametrize(("table", "uniforms", "kind"), [
    ([], [0.5], ValueError),
    ([0.5, 0.2], [0.5], ValueError),
    ([0.5, np.inf], [0.5], ValueError),
    ([0.5, 1.0], [np.nan], ValueError),
])
def test_categorical_indices_refuse_invalid_tables_and_variates(table, uniforms, kind) -> None:
    # Act / Assert
    with pytest.raises(kind):
        multivariate.Multinomial.categorical_indices(table, uniforms)


def test_categorical_indices_return_one_index_per_uniform() -> None:
    # Arrange
    cumulative = np.cumsum(np.full(50, 0.02))
    uniforms = np.random.default_rng(6).random(200000)
    # Act
    bulk = multivariate.Multinomial.categorical_indices(cumulative, uniforms)
    # Assert
    assert bulk.size == uniforms.size


@pytest.mark.benchmark
def test_categorical_indices_cross_once_and_beat_per_element_calls() -> None:
    # Arrange
    cumulative = np.cumsum(np.full(50, 0.02))
    uniforms = np.random.default_rng(6).random(200000)
    nd = multivariate._nd()
    net_table = to_net_doubles(cumulative)
    one = to_net_doubles([0.5])
    multivariate.Multinomial.categorical_indices(cumulative, uniforms[:10])
    # Act
    start = time.perf_counter()
    multivariate.Multinomial.categorical_indices(cumulative, uniforms)
    bulk_seconds = time.perf_counter() - start
    start = time.perf_counter()
    for _ in range(2000):
        nd.Multinomial.CategoricalIndices(net_table, one)
    per_call = (time.perf_counter() - start) / 2000
    # Assert: the whole array costs less than a twentieth of crossing per element
    assert bulk_seconds < per_call * uniforms.size / 20


# ---------------------------------------------------------------------------------------------
# BivariateEmpirical

X1 = [0.0, 1.0, 2.0]
X2 = [0.0, 10.0]
P = [[0.0, 0.1], [0.2, 0.5], [0.3, 1.0]]


def test_bivariate_empirical_cdf_equals_scipy_bilinear_interpolation() -> None:
    # Arrange
    empirical = multivariate.BivariateEmpirical(X1, X2, P)
    oracle = pytest.importorskip("scipy.interpolate").RegularGridInterpolator((X1, X2), np.array(P))
    rng = np.random.default_rng(7)
    a, b = rng.uniform(0.0, 2.0, 50), rng.uniform(0.0, 10.0, 50)
    # Act
    values = empirical.cdf(a, b)
    # Assert
    np.testing.assert_allclose(values, oracle(np.column_stack([a, b])), rtol=1e-13, atol=1e-15)
    assert empirical.cdf([a[0], b[0]]) == values[0] == empirical.cdf(a[0], b[0])
    np.testing.assert_allclose(empirical.cdf(np.column_stack([a, b])), values, rtol=0, atol=0)
    assert np.isnan(empirical.pdf(0.5, 5.0)) and np.isnan(empirical.pdf([0.5, 5.0]))
    assert empirical.dimension == 2
    assert empirical.type is multivariate.MultivariateDistributionType.BIVARIATE_EMPIRICAL_DISTRIBUTION


def test_bivariate_empirical_tables_transforms_and_validation() -> None:
    # Arrange
    empirical = multivariate.BivariateEmpirical(x1_transform="logarithmic")
    # Act
    empirical.set_parameters([1.0, 10.0, 100.0], X2, P)
    empirical.x2_transform = "none"
    empirical.probability_transform = "normal_z"
    clone = empirical.clone()
    # Assert
    np.testing.assert_array_equal(empirical.x1_values, [1.0, 10.0, 100.0])
    np.testing.assert_array_equal(empirical.x2_values, X2)
    np.testing.assert_array_equal(empirical.probability_values, P)
    assert empirical.x1_transform.name == "LOGARITHMIC" and empirical.x2_transform.name == "NONE"
    assert empirical.probability_transform.name == "NORMAL_Z"
    assert isinstance(clone, multivariate.BivariateEmpirical) and clone.x1_transform == empirical.x1_transform
    log_mid = empirical.cdf(np.sqrt(10.0), 0.0)
    assert 0.0 < log_mid < 0.2
    assert empirical.validate_parameters(X1, X2, P) is None
    assert "number of rows" in empirical.validate_parameters([0.0, 1.0], X2, P)


# ---------------------------------------------------------------------------------------------
# MultivariateNormal

COV3 = [[1.0, 0.5, 0.2], [0.5, 2.0, 0.3], [0.2, 0.3, 1.5]]
MEAN3 = [0.5, -1.0, 2.0]


def test_multivariate_normal_density_cdf_and_moments_equal_scipy() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal(MEAN3, COV3)
    oracle = stats.multivariate_normal(MEAN3, COV3)
    points = np.array([[0.0, 0.0, 0.0], [1.0, -1.0, 3.0], [-1.0, -3.0, 1.0]])
    # Act / Assert
    np.testing.assert_allclose(mvn.pdf(points), oracle.pdf(points), rtol=1e-13)
    np.testing.assert_allclose(mvn.log_pdf(points), oracle.logpdf(points), rtol=1e-13)
    mvn.absolute_error, mvn.relative_error, mvn.max_evaluations = 1e-7, 1e-7, 200000
    assert (mvn.absolute_error, mvn.relative_error, mvn.max_evaluations) == (1e-7, 1e-7, 200000)
    exact = stats.multivariate_normal.cdf(points, MEAN3, COV3, abseps=1e-8, releps=1e-8, maxpts=10**6)
    np.testing.assert_allclose(mvn.cdf(points), exact, atol=1e-6)
    np.testing.assert_allclose(mvn.variance, np.diag(COV3), rtol=1e-15)
    np.testing.assert_allclose(mvn.standard_deviation, np.sqrt(np.diag(COV3)), rtol=1e-15)
    for vector in (mvn.median, mvn.mode, mvn.mean):
        np.testing.assert_array_equal(vector, MEAN3)
    diff = points[1] - MEAN3
    assert mvn.mahalanobis(points[1]) == pytest.approx(diff @ np.linalg.solve(COV3, diff), rel=1e-13)
    assert mvn.is_positive_definite and mvn.is_density_valid and mvn.decomposition == "Cholesky"


def test_multivariate_normal_interval_mvndst_and_bivariate_functions_equal_scipy() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal(MEAN3, COV3)
    mvn.absolute_error, mvn.relative_error, mvn.max_evaluations = 1e-9, 1e-9, 10**7
    lower, upper = np.array([-1.0, -2.0, 1.0]), np.array([1.0, 0.0, 3.0])
    std = np.sqrt(np.diag(COV3))
    correlation = np.array(COV3) / np.outer(std, std)
    correl = [correlation[1, 0], correlation[2, 0], correlation[2, 1]]
    box = stats.multivariate_normal.cdf(upper, MEAN3, COV3, lower_limit=lower, abseps=1e-9, releps=1e-9,
                                        maxpts=10**7)
    bvn_cov = [[1.0, 0.4], [0.4, 1.0]]
    # Act
    interval = mvn.interval(lower, upper)
    error, value, inform = mvn.mvndst((lower - MEAN3) / std, (upper - MEAN3) / std, [2, 2, 2], correl,
                                      200000, 1e-8, 0.0)
    # Assert
    assert interval == pytest.approx(box, abs=1e-6)
    assert value == pytest.approx(box, abs=1e-6) and inform == 0 and 0.0 <= error < 1e-6
    for z1, z2 in [(0.0, 0.0), (1.0, -0.5), (-2.0, 1.5)]:
        lower_tail = stats.multivariate_normal.cdf([z1, z2], [0.0, 0.0], bvn_cov, abseps=1e-12, releps=1e-12)
        upper_tail = stats.multivariate_normal.cdf([-z1, -z2], [0.0, 0.0], bvn_cov, abseps=1e-12,
                                                   releps=1e-12)
        assert multivariate.MultivariateNormal.bivariate_cdf(-z1, -z2, 0.4) == pytest.approx(lower_tail, abs=1e-12)
        assert multivariate.MultivariateNormal.bivariate_cdf(z1, z2, 0.4) == pytest.approx(upper_tail, abs=1e-12)
        assert multivariate.MultivariateNormal.bvu(z1, z2, 0.4) == pytest.approx(upper_tail, abs=1e-12)
        assert multivariate.MultivariateNormal.bivnor(z1, z2, 0.4) == pytest.approx(upper_tail, abs=1e-9)
    for z in (-6.0, -1.0, 0.0, 2.5):
        assert multivariate.MultivariateNormal.mvnphi(z) == pytest.approx(stats.norm.cdf(z), rel=1e-14)
    for p in (1e-10, 0.025, 0.5, 0.9):
        assert multivariate.MultivariateNormal.phinvs(p) == pytest.approx(stats.norm.ppf(p), rel=1e-14)


def test_multivariate_normal_conditional_and_marginal_equal_textbook_formulas() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal(MEAN3, COV3)
    s = np.array(COV3)
    mu = np.array(MEAN3)
    observed, free = [1], [0, 2]
    gain = s[np.ix_(free, observed)] @ np.linalg.inv(s[np.ix_(observed, observed)])
    # Act
    conditional = mvn.conditional(observed, [0.5])
    marginal = mvn.marginal(2, 0)
    single = mvn.independent_normal(1)
    # Assert
    np.testing.assert_allclose(conditional.mean, mu[free] + gain @ (np.array([0.5]) - mu[observed]), rtol=1e-14)
    np.testing.assert_allclose(conditional.covariance,
                               s[np.ix_(free, free)] - gain @ s[np.ix_(observed, free)], rtol=1e-14)
    np.testing.assert_array_equal(marginal.mean, mu[[2, 0]])
    np.testing.assert_array_equal(marginal.covariance, s[np.ix_([2, 0], [2, 0])])
    assert isinstance(single, distributions.Normal)
    assert single.mean == -1.0 and single.standard_deviation == pytest.approx(np.sqrt(2.0), rel=1e-15)


def test_multivariate_normal_constructors_setters_and_generators() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal(3)
    # Act
    mvn.set_parameters(MEAN3, COV3)
    mvn.set_mean([0.0, 0.0, 0.0])
    accepted = mvn.try_set_covariance(np.eye(3))
    refused = mvn.try_set_parameters([0.0, 0.0, 0.0], [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    mvn.mvnuni = sampling.MersenneTwister(9)
    sample = mvn.generate_random_values(40, sampling.MersenneTwister(4))
    lhs = mvn.latin_hypercube_random_values(30, sampling.MersenneTwister(4))
    stratified = mvn.stratified_random_values(sampling.stratify_probabilities(0.05, 0.95, 12), 3)
    # Assert
    assert accepted and not refused
    np.testing.assert_array_equal(mvn.covariance, np.eye(3))
    np.testing.assert_array_equal(sample, mvn.generate_random_values(40, 4))
    np.testing.assert_array_equal(lhs, mvn.latin_hypercube_random_values(30, 4))
    assert stratified.shape == (12, 3)
    assert mvn.validate_parameters(MEAN3, COV3) is None
    assert "positive-definite" in mvn.validate_parameters(MEAN3, [[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    assert multivariate.MultivariateNormal.default_mvnuni_seed() == 12345
    bivariate = multivariate.MultivariateNormal.bivariate(0.0, 1.0, 1.0, 2.0, 0.5)
    np.testing.assert_allclose(bivariate.covariance, [[1.0, 1.0], [1.0, 4.0]], rtol=1e-15)
    univariate = multivariate.MultivariateNormal.univariate(2.0, 3.0)
    assert univariate.dimension == 1 and univariate.covariance[0, 0] == 9.0
    assert multivariate.MultivariateNormal([0.0, 0.0], np.eye(2), "Cholesky").decomposition == "Cholesky"
    with pytest.raises(ValueError, match="mean vector"):
        multivariate.MultivariateNormal(2, np.eye(2))


def test_wrap_returns_the_typed_class_and_shares_the_object() -> None:
    # Arrange
    nd = multivariate._nd()
    nets = [nd.MultivariateNormal(2), nd.Dirichlet(3, 2.0), nd.Multinomial(2, to_net_doubles([0.5, 0.5])),
            nd.MultivariateStudentT(2, 3.0), nd.BivariateEmpirical()]
    # Act
    wrapped = [multivariate.wrap(net) for net in nets]
    # Assert
    assert [type(w).__name__ for w in wrapped] == ["MultivariateNormal", "Dirichlet", "Multinomial",
                                                    "MultivariateStudentT", "BivariateEmpirical"]
    assert all(w.net is not None for w in wrapped)
    assert isinstance(multivariate.Dirichlet([1.0, 2.0]).clone(), multivariate.Dirichlet)
    assert multivariate.MultivariateDistributionType("multivariate_student_t").value == "MultivariateStudentT"


# ---------------------------------------------------------------------------------------------
# Copulas: densities, generators and tail dependence against Nelsen's closed forms


def _gaussian_density(u: np.ndarray, v: np.ndarray, rho: float) -> np.ndarray:
    a, b = stats.norm.ppf(u), stats.norm.ppf(v)
    return np.exp(-(rho**2 * (a**2 + b**2) - 2.0 * rho * a * b) / (2.0 * (1.0 - rho**2))) / np.sqrt(1.0 - rho**2)


def _clayton_density(u: np.ndarray, v: np.ndarray, theta: float) -> np.ndarray:
    return (1.0 + theta) * (u * v) ** (-theta - 1.0) * (u**-theta + v**-theta - 1.0) ** (-2.0 - 1.0 / theta)


def _student_density(u: np.ndarray, v: np.ndarray, rho: float, nu: float) -> np.ndarray:
    a, b = stats.t.ppf(u, nu), stats.t.ppf(v, nu)
    joint = stats.multivariate_t([0.0, 0.0], [[1.0, rho], [rho, 1.0]], df=nu).pdf(np.column_stack([a, b]))
    return joint / (stats.t.pdf(a, nu) * stats.t.pdf(b, nu))


UV = np.random.default_rng(8).uniform(0.02, 0.98, size=(2, 40))


def test_copula_densities_and_cdfs_equal_closed_forms() -> None:
    # Arrange
    u, v = UV
    theta = 1.7
    gumbel_cdf = np.exp(-((-np.log(u)) ** theta + (-np.log(v)) ** theta) ** (1.0 / theta))
    frank_cdf = -np.log1p(np.expm1(-theta * u) * np.expm1(-theta * v) / np.expm1(-theta)) / theta
    joe_cdf = 1.0 - ((1 - u) ** theta + (1 - v) ** theta - (1 - u) ** theta * (1 - v) ** theta) ** (1 / theta)
    amh_cdf = u * v / (1.0 - 0.4 * (1.0 - u) * (1.0 - v))
    gaussian = stats.multivariate_normal([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]])
    # Act / Assert
    np.testing.assert_allclose(copulas.NormalCopula(0.6).pdf(u, v), _gaussian_density(u, v, 0.6), rtol=1e-10)
    np.testing.assert_allclose(copulas.NormalCopula(0.6).cdf(u, v),
                               gaussian.cdf(np.column_stack([stats.norm.ppf(u), stats.norm.ppf(v)])), atol=1e-7)
    np.testing.assert_allclose(copulas.ClaytonCopula(theta).pdf(u, v), _clayton_density(u, v, theta), rtol=1e-10)
    np.testing.assert_allclose(copulas.ClaytonCopula(theta).cdf(u, v),
                               (u**-theta + v**-theta - 1.0) ** (-1.0 / theta), rtol=1e-13)
    np.testing.assert_allclose(copulas.StudentTCopula(0.5, 4.0).pdf(u, v), _student_density(u, v, 0.5, 4.0),
                               rtol=1e-7)
    np.testing.assert_allclose(copulas.GumbelCopula(theta).cdf(u, v), gumbel_cdf, rtol=1e-13)
    np.testing.assert_allclose(copulas.FrankCopula(theta).cdf(u, v), frank_cdf, rtol=1e-12)
    np.testing.assert_allclose(copulas.JoeCopula(theta).cdf(u, v), joe_cdf, rtol=1e-12)
    np.testing.assert_allclose(copulas.AMHCopula(0.4).cdf(u, v), amh_cdf, rtol=1e-13)
    clayton = copulas.ClaytonCopula(theta)
    np.testing.assert_allclose(clayton.log_pdf(u, v), np.log(_clayton_density(u, v, theta)), rtol=1e-10)
    np.testing.assert_allclose(clayton.or_joint_exceedance_probability(u, v), 1.0 - clayton.cdf(u, v), rtol=1e-15)
    np.testing.assert_allclose(clayton.and_joint_exceedance_probability(u, v), 1.0 - u - v + clayton.cdf(u, v),
                               rtol=1e-13)


def test_archimedean_generators_equal_nelsen_table_4_1() -> None:
    # Arrange: every generator, Clayton's included, is the one of Nelsen (2006, Table 4.1)
    t = np.array([0.1, 0.3, 0.7, 0.95])
    theta = 2.0
    cases = [
        (copulas.ClaytonCopula(theta), (t**-theta - 1.0) / theta),
        (copulas.FrankCopula(theta), -np.log(np.expm1(-theta * t) / np.expm1(-theta))),
        (copulas.JoeCopula(theta), -np.log1p(-((1.0 - t) ** theta))),
        (copulas.AMHCopula(0.5), np.log((1.0 - 0.5 * (1.0 - t)) / t)),
    ]
    for copula, expected in cases:
        # Act
        phi = copula.generator(t)
        # Assert
        np.testing.assert_allclose(phi, expected, rtol=1e-13)
        np.testing.assert_allclose(copula.generator_inverse(phi), t, rtol=1e-12)
        step = 1e-6
        numeric_prime = (copula.generator(t + step) - copula.generator(t - step)) / (2.0 * step)
        np.testing.assert_allclose(copula.generator_prime(t), numeric_prime, rtol=1e-6)
        numeric_prime2 = (copula.generator_prime(t + step) - copula.generator_prime(t - step)) / (2.0 * step)
        np.testing.assert_allclose(copula.generator_prime2(t), numeric_prime2, rtol=1e-5)
        np.testing.assert_allclose(copula.generator_prime_inverse(copula.generator_prime(t)), t, rtol=1e-8)
        assert isinstance(copula.generator(0.3), float)


def test_gumbel_generator_and_its_inverses_equal_nelsen_table_4_1() -> None:
    # Arrange
    t = np.array([0.1, 0.3, 0.7, 0.95])
    gumbel = copulas.GumbelCopula(2.0)
    # Act
    phi = gumbel.generator(t)
    # Assert
    np.testing.assert_allclose(phi, (-np.log(t)) ** 2.0, rtol=1e-13)
    np.testing.assert_allclose(gumbel.generator_inverse(phi), t, rtol=1e-12)
    np.testing.assert_allclose(gumbel.generator_prime_inverse(gumbel.generator_prime(t)), t, rtol=1e-8)


def test_gumbel_generator_derivatives_equal_finite_differences() -> None:
    # Arrange
    t = np.array([0.1, 0.3, 0.7, 0.95])
    gumbel = copulas.GumbelCopula(2.0)
    step = 1e-6
    # Act
    numeric_prime = (gumbel.generator(t + step) - gumbel.generator(t - step)) / (2.0 * step)
    numeric_prime2 = (gumbel.generator_prime(t + step) - gumbel.generator_prime(t - step)) / (2.0 * step)
    # Assert
    np.testing.assert_allclose(gumbel.generator_prime(t), numeric_prime, rtol=1e-6)
    np.testing.assert_allclose(gumbel.generator_prime(t), -2.0 * (-np.log(t)) / t, rtol=1e-14)
    np.testing.assert_allclose(gumbel.generator_prime2(t), numeric_prime2, rtol=1e-5)
    assert np.all(gumbel.generator_prime(t) < 0.0) and np.all(gumbel.generator_prime2(t) > 0.0)


def test_tail_dependence_equals_closed_forms() -> None:
    # Arrange
    theta, rho, nu = 2.5, 0.5, 4.0
    t_tail = 2.0 * stats.t.cdf(-np.sqrt((nu + 1.0) * (1.0 - rho) / (1.0 + rho)), nu + 1.0)
    # Act / Assert
    assert copulas.ClaytonCopula(theta).lower_tail_dependence == pytest.approx(2.0 ** (-1.0 / theta), rel=1e-14)
    assert copulas.GumbelCopula(theta).upper_tail_dependence == pytest.approx(2.0 - 2.0 ** (1.0 / theta), rel=1e-14)
    assert copulas.JoeCopula(theta).upper_tail_dependence == pytest.approx(2.0 - 2.0 ** (1.0 / theta), rel=1e-14)
    t_copula = copulas.StudentTCopula(rho, nu)
    assert t_copula.upper_tail_dependence == pytest.approx(t_tail, rel=1e-10)
    assert t_copula.lower_tail_dependence == t_copula.upper_tail_dependence
    for copula in (copulas.NormalCopula(0.9), copulas.FrankCopula(5.0), copulas.IndependenceCopula()):
        assert copula.upper_tail_dependence == 0.0 and copula.lower_tail_dependence == 0.0


@pytest.mark.parametrize(("factory", "theta_of_tau"), [
    (copulas.ClaytonCopula, lambda tau: 2.0 * tau / (1.0 - tau)),
    (copulas.GumbelCopula, lambda tau: 1.0 / (1.0 - tau)),
])
def test_set_theta_from_tau_inverts_kendalls_tau(factory, theta_of_tau) -> None:
    # Arrange
    x, y = _gaussian_pair(300, 0.6, 9)
    tau = stats.kendalltau(x, y).statistic
    copula = factory()
    # Act
    copula.set_theta_from_tau(x, y)
    # Assert
    assert copula.theta == pytest.approx(theta_of_tau(tau), rel=1e-12)


def test_frank_joe_and_amh_theta_from_tau_reproduce_the_sample_tau() -> None:
    # Arrange: Nelsen (2006) Example 5.4 gives tau for Frank through the Debye function
    x, y = _gaussian_pair(400, 0.3, 10)
    tau = stats.kendalltau(x, y).statistic

    def frank_tau(theta: float) -> float:
        debye = integrate.quad(lambda s: s / np.expm1(s), 0.0, theta)[0] / theta
        return 1.0 - 4.0 / theta * (1.0 - debye)

    def joe_tau(theta: float) -> float:
        return 1.0 + 4.0 / theta**2 * integrate.quad(
            lambda s: s * np.log(s) * (1.0 - s) ** (2.0 / theta - 2.0), 0.0, 1.0)[0]

    def amh_tau(theta: float) -> float:
        return 1.0 - 2.0 * (theta + (1.0 - theta) ** 2 * np.log1p(-theta)) / (3.0 * theta**2)

    # Act / Assert
    for copula, tau_of_theta in ((copulas.FrankCopula(), frank_tau), (copulas.JoeCopula(), joe_tau),
                                 (copulas.AMHCopula(), amh_tau)):
        copula.set_theta_from_tau(x, y)
        assert tau_of_theta(copula.theta) == pytest.approx(tau, abs=1e-6)


# ---------------------------------------------------------------------------------------------
# Copula estimation against independent maximizations of the same likelihoods


def _maximize(negative_log_likelihood, bounds: tuple[float, float]) -> float:
    return optimize.minimize_scalar(negative_log_likelihood, bounds=bounds, method="bounded",
                                    options={"xatol": 1e-10}).x


def test_pseudo_likelihood_fits_equal_scipy_maximizations() -> None:
    # Arrange
    x, y = _gaussian_pair(250, 0.55, 12)
    u, v = _plotting_positions(x), _plotting_positions(y)
    normal, clayton = copulas.NormalCopula(), copulas.ClaytonCopula()
    expected_rho = _maximize(lambda r: -np.sum(np.log(_gaussian_density(u, v, r))), (-0.99, 0.99))
    expected_theta = _maximize(lambda th: -np.sum(np.log(_clayton_density(u, v, th))), (0.01, 20.0))
    # Act
    normal.estimate(u, v)
    clayton.estimate(u, v, "pseudo_likelihood")
    # Assert
    assert normal.theta == pytest.approx(expected_rho, abs=1e-5)
    assert clayton.theta == pytest.approx(expected_theta, abs=1e-4)
    assert normal.pseudo_log_likelihood(u, v) == pytest.approx(np.sum(np.log(_gaussian_density(u, v, normal.theta))),
                                                               rel=1e-10)


def test_student_t_pseudo_likelihood_fit_equals_a_scipy_maximization() -> None:
    # Arrange
    rng = np.random.default_rng(13)
    z = stats.multivariate_t([0.0, 0.0], [[1.0, 0.6], [0.6, 1.0]], df=4.0).rvs(size=400, random_state=rng)
    u, v = _plotting_positions(z[:, 0]), _plotting_positions(z[:, 1])
    copula = copulas.StudentTCopula()

    def negative(p: np.ndarray) -> float:
        return -float(np.sum(np.log(_student_density(u, v, p[0], p[1]))))

    expected = optimize.minimize(negative, [0.5, 5.0], method="L-BFGS-B", bounds=[(-0.99, 0.99), (1.0, 60.0)]).x
    # Act
    copula.estimate(u, v, copulas.CopulaEstimationMethod.PSEUDO_LIKELIHOOD)
    # Assert
    assert copula.theta == pytest.approx(expected[0], abs=5e-3)
    assert copula.degrees_of_freedom == pytest.approx(expected[1], rel=0.05)
    assert -negative(copula.copula_parameters) >= -negative(expected) - 1e-3


def test_inference_from_margins_fit_equals_a_scipy_maximization() -> None:
    # Arrange
    x, y = _gaussian_pair(200, 0.5, 14)
    x, y = 3.0 + 2.0 * x, -1.0 + 0.5 * y
    mx, my = distributions.Normal(), distributions.Normal()
    mx.estimate(x)
    my.estimate(y)
    u = stats.norm.cdf(x, mx.mean, mx.standard_deviation)
    v = stats.norm.cdf(y, my.mean, my.standard_deviation)
    copula = copulas.ClaytonCopula(1.0, mx, my)
    expected = _maximize(lambda th: -np.sum(np.log(_clayton_density(u, v, th))), (0.01, 20.0))
    # Act
    copula.estimate(x, y, "inference_from_margins")
    # Assert
    assert copula.theta == pytest.approx(expected, abs=1e-4)
    assert copula.ifm_log_likelihood(x, y) == pytest.approx(np.sum(np.log(_clayton_density(u, v, copula.theta))),
                                                            rel=1e-9)


def test_full_likelihood_with_normal_margins_is_the_bivariate_normal_mle() -> None:
    # Arrange: a Gaussian copula with normal margins is the bivariate normal, whose MLE is closed form
    x, y = _gaussian_pair(300, 0.6, 15)
    x, y = 10.0 + 3.0 * x, -2.0 + 0.5 * y
    copula = copulas.NormalCopula(0.0, distributions.Normal(), distributions.Normal())
    rho = np.corrcoef(x, y)[0, 1]
    # Act
    copula.estimate(x, y, copulas.CopulaEstimationMethod.FULL_LIKELIHOOD)
    # Assert
    mx, my = copula.marginal_distribution_x, copula.marginal_distribution_y
    assert copula.theta == pytest.approx(rho, abs=2e-3)
    assert mx.mean == pytest.approx(x.mean(), abs=5e-3 * x.std()) and my.mean == pytest.approx(y.mean(), abs=5e-3)
    assert mx.standard_deviation == pytest.approx(x.std(), rel=5e-3)
    assert my.standard_deviation == pytest.approx(y.std(), rel=5e-3)
    sx, sy = mx.standard_deviation, my.standard_deviation
    joint = stats.multivariate_normal([mx.mean, my.mean], [[sx**2, copula.theta * sx * sy],
                                                           [copula.theta * sx * sy, sy**2]])
    assert copula.log_likelihood(x, y) == pytest.approx(joint.logpdf(np.column_stack([x, y])).sum(), rel=1e-8)


def test_estimation_refuses_missing_margins_and_ignores_the_independence_copula() -> None:
    # Arrange
    x, y = _gaussian_pair(50, 0.5, 16)
    independence = copulas.IndependenceCopula()
    # Act
    independence.estimate(x, y, "full_likelihood")
    # Assert
    assert independence.number_of_copula_parameters == 0
    with pytest.raises(ValueError, match="marginal"):
        copulas.ClaytonCopula(2.0).estimate(x, y, "full_likelihood")
    with pytest.raises(ValueError, match="marginal"):
        copulas.ClaytonCopula(2.0).ifm_log_likelihood(x, y)


def test_frozen_marginal_is_attached_as_a_clone() -> None:
    # Arrange
    frozen = distributions.Normal(0.0, 1.0).frozen()
    x, y = _gaussian_pair(100, 0.5, 17)
    copula = copulas.NormalCopula(0.2, frozen, distributions.Normal())
    # Act
    copula.estimate(x * 2.0 + 5.0, y, "full_likelihood")
    # Assert
    assert frozen.mean == 0.0 and copula.marginal_distribution_x.mean != 0.0


# ---------------------------------------------------------------------------------------------
# Factory, parameters and the XML round trip


def test_factory_creates_every_family_and_wrap_returns_the_typed_class() -> None:
    # Act / Assert
    for kind in copulas.CopulaType:
        copula = copulas.create_copula(kind)
        assert copula.type is kind
        assert type(copula).__name__ == type(copula.clone()).__name__
        assert copula.parameters_valid
    assert isinstance(copulas.create_copula("clayton"), copulas.ArchimedeanCopula)
    assert isinstance(copulas.wrap(copulas.FrankCopula(2.0).net), copulas.FrankCopula)


@pytest.mark.parametrize("copula", [
    copulas.NormalCopula(-0.35), copulas.StudentTCopula(0.4, 6.5), copulas.ClaytonCopula(1.0 / 3.0),
    copulas.GumbelCopula(2.25), copulas.FrankCopula(-4.0), copulas.JoeCopula(1.5), copulas.AMHCopula(0.3),
    copulas.IndependenceCopula(),
], ids=lambda c: c.type.name)
def test_xml_round_trip_restores_type_and_parameters_bitwise(copula) -> None:
    # Act
    restored = copulas.copula_from_xml(copula.to_xml())
    # Assert
    assert type(restored) is type(copula)
    np.testing.assert_array_equal(restored.copula_parameters, copula.copula_parameters)
    assert restored.marginal_distribution_x is None


def test_xml_round_trip_with_marginals() -> None:
    # Arrange
    copula = copulas.GumbelCopula(1.8, distributions.Normal(2.0, 0.5), distributions.Gumbel(10.0, 3.0))
    # Act
    text = copula.to_xml(include_marginals=True)
    restored = copulas.copula_from_xml(text, include_marginals=True)
    without = copulas.copula_from_xml(text)
    # Assert
    assert "MarginalDistributionX" in text and restored.theta == 1.8
    assert restored.marginal_distribution_x == distributions.Normal(2.0, 0.5)
    assert restored.marginal_distribution_y == distributions.Gumbel(10.0, 3.0)
    assert without.marginal_distribution_x is None
    np.testing.assert_array_equal(restored.generate_random_values(20, 3), copula.generate_random_values(20, 3))


@pytest.mark.parametrize("text", [
    '<Copula Type="NotACopula" Parameters="1" />',
    '<Copula Type="Clayton" Parameters="1|2" />',
    '<Copula Type="Clayton" Parameters="abc" />',
    '<Copula Type="Gumbel" Parameters="0.5" />',
])
def test_xml_refuses_malformed_elements(text: str) -> None:
    # Act / Assert
    with pytest.raises(ValueError):
        copulas.copula_from_xml(text)


def test_copula_parameter_members() -> None:
    # Arrange
    clayton = copulas.ClaytonCopula(2.0)
    student = copulas.StudentTCopula(0.3, 5.0)
    # Act
    clayton.set_copula_parameters([3.5])
    student.set_copula_parameters([0.6, 8.0])
    student.degrees_of_freedom = 9.0
    # Assert
    assert clayton.theta == 3.5 and clayton.number_of_copula_parameters == 1
    np.testing.assert_array_equal(student.copula_parameters, [0.6, 9.0])
    assert student.number_of_copula_parameters == 2 and student.degrees_of_freedom == 9.0
    assert clayton.validate_parameter(-2.0) is not None and clayton.validate_parameter(1.0) is None
    assert student.validate_parameters(0.5, 4.0) is None and student.validate_parameters(1.5, 4.0) is not None
    assert clayton.theta_minimum == -1.0 and clayton.theta_maximum == np.inf
    assert clayton.parameter_to_string[0][1] == "3.5" and clayton.parameter_name_short_form
    assert clayton.parameter_property_name == "Theta"
    assert clayton.display_name == "Clayton" and clayton.short_display_name
    bounds = clayton.parameter_constraints([0.1, 0.4, 0.7], [0.2, 0.5, 0.6])
    assert bounds.shape == (1, 2) and bounds[0, 0] < bounds[0, 1]
    clayton.theta = -5.0
    assert not clayton.parameters_valid
    with pytest.raises(TypeError, match="both"):
        copulas.StudentTCopula(0.5)


def test_copula_sampling_from_a_generator_reproduces_the_seed() -> None:
    # Arrange
    copula = copulas.FrankCopula(4.0)
    # Act
    seeded = copula.generate_random_values(100, 21)
    generated = copula.generate_random_values(100, sampling.MersenneTwister(21))
    # Assert
    np.testing.assert_array_equal(generated, seeded)
    assert stats.kendalltau(seeded[:, 0], seeded[:, 1]).statistic > 0.2
