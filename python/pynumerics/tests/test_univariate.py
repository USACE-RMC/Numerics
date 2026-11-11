"""Univariate distribution families against scipy.stats, closed forms and the Numerics members they wrap.

Each family is checked for its density or mass, CDF, quantile, survival quantile and moments against
an independent scipy.stats implementation (or a textbook closed form where scipy has none), so the
tolerances cover only the different operation order or the documented Numerics approximation.
"""

from __future__ import annotations

import copy
import math
import time

import numpy as np
import pytest

from pynumerics import bootstrap, distributions, enums, sampling

stats = pytest.importorskip("scipy.stats")
integrate = pytest.importorskip("scipy.integrate")
d = distributions


def _bits(values: object) -> list[int]:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64).tolist()


def _hosking_y(x: np.ndarray, location: float, scale: float, shape: float) -> np.ndarray:
    """Return Hosking's reduced variate of the generalized logistic and normal families."""
    z = (x - location) / scale
    return z if shape == 0 else -np.log1p(-shape * z) / shape


# family, scipy frozen oracle (or None), points, discrete
CONTINUOUS = [
    (lambda: d.Normal(3.0, 2.0), stats.norm(3.0, 2.0)),
    (lambda: d.GeneralizedExtremeValue(100.0, 15.0, -0.2), stats.genextreme(-0.2, 100.0, 15.0)),
    (lambda: d.GeneralizedPareto(10.0, 4.0, 0.15), stats.genpareto(-0.15, 10.0, 4.0)),
    (lambda: d.Gumbel(50.0, 8.0), stats.gumbel_r(50.0, 8.0)),
    (lambda: d.LogNormal(2.0, 0.3), stats.lognorm(0.3 * math.log(10.0), scale=100.0)),
    (lambda: d.LogNormal(1.5, 0.4, base=math.e), stats.lognorm(0.4, scale=math.exp(1.5))),
    (lambda: d.PearsonTypeIII(100.0, 20.0, 0.8), stats.pearson3(0.8, 100.0, 20.0)),
    (lambda: d.KappaFour(10.0, 2.0, 0.1, 0.2), stats.kappa4(0.2, 0.1, 10.0, 2.0)),
    (lambda: d.Uniform(2.0, 5.0), stats.uniform(2.0, 3.0)),
    (lambda: d.Beta(2.5, 4.0), stats.beta(2.5, 4.0)),
    (lambda: d.Cauchy(1.0, 2.0), stats.cauchy(1.0, 2.0)),
    (lambda: d.ChiSquared(7), stats.chi2(7)),
    (lambda: d.Exponential(3.0, 2.5), stats.expon(3.0, 2.5)),
    (lambda: d.GammaDistribution(2.0, 3.5), stats.gamma(3.5, scale=2.0)),
    (lambda: d.GeneralizedBeta(2.0, 3.0, 10.0, 30.0), stats.beta(2.0, 3.0, 10.0, 20.0)),
    (lambda: d.InverseChiSquared(9, 2.0), stats.invgamma(4.5, scale=9.0)),
    (lambda: d.InverseGamma(1.5, 6.0), stats.invgamma(6.0, scale=1.5)),
    (lambda: d.Logistic(5.0, 1.5), stats.logistic(5.0, 1.5)),
    (lambda: d.NoncentralT(8.0, 1.5), stats.nct(8.0, 1.5)),
    (lambda: d.Pareto(2.0, 6.0), stats.pareto(6.0, scale=2.0)),
    (lambda: d.Pert(10.0, 30.0, 100.0), stats.beta(1 + 4 * 20 / 90, 1 + 4 * 70 / 90, 10.0, 90.0)),
    (lambda: d.Rayleigh(3.0), stats.rayleigh(scale=3.0)),
    (lambda: d.StudentT(1.0, 2.0, 9.0), stats.t(9.0, 1.0, 2.0)),
    (lambda: d.Triangular(0.0, 3.0, 10.0), stats.triang(0.3, 0.0, 10.0)),
    (lambda: d.TruncatedNormal(1.0, 2.0, -1.0, 4.0), stats.truncnorm(-1.0, 1.5, 1.0, 2.0)),
    (lambda: d.TruncatedDistribution(d.Normal(1.0, 2.0), -1.0, 4.0), stats.truncnorm(-1.0, 1.5, 1.0, 2.0)),
    (lambda: d.Weibull(10.0, 2.5), stats.weibull_min(2.5, scale=10.0)),
]


@pytest.mark.parametrize(("make", "oracle"), CONTINUOUS, ids=lambda v: getattr(v, "__name__", None) or "")
def test_continuous_families_equal_scipy(make, oracle) -> None:
    # Arrange
    distribution = make()
    p = np.array([0.001, 0.05, 0.3, 0.5, 0.7, 0.95, 0.999])
    x = oracle.ppf(p)
    mean, variance, skew, excess = (float(v) for v in oracle.stats(moments="mvsk"))
    # Act / Assert
    np.testing.assert_allclose(distribution.pdf(x), oracle.pdf(x), rtol=1e-9)
    np.testing.assert_allclose(distribution.cdf(x), p, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(distribution.inverse_cdf(p), x, rtol=1e-8)
    np.testing.assert_allclose(distribution.inverse_ccdf(1.0 - p), oracle.isf(1.0 - p), rtol=1e-8)
    assert distribution.median == pytest.approx(float(oracle.median()), rel=1e-9)
    if math.isfinite(mean):
        assert distribution.mean == pytest.approx(mean, rel=1e-9, abs=1e-12)
    if math.isfinite(variance):
        assert distribution.variance == pytest.approx(variance, rel=1e-8)
        assert distribution.standard_deviation == pytest.approx(math.sqrt(variance), rel=1e-8)
    if math.isfinite(skew):
        assert distribution.skewness == pytest.approx(skew, rel=1e-6, abs=1e-9)
    if math.isfinite(excess):
        assert distribution.kurtosis == pytest.approx(excess + 3.0, rel=1e-6)  # Numerics: not excess


@pytest.mark.parametrize(("make", "oracle", "support"), [
    (lambda: d.Poisson(3.5), stats.poisson(3.5), np.arange(0, 12)),
    (lambda: d.Bernoulli(0.3), stats.bernoulli(0.3), np.array([0, 1])),
    (lambda: d.Binomial(0.3, 20), stats.binom(20, 0.3), np.arange(0, 21)),
    (lambda: d.Geometric(0.3), stats.geom(0.3, loc=-1), np.arange(0, 15)),
    (lambda: d.UniformDiscrete(2.0, 6.0), stats.randint(2, 7), np.arange(2, 7)),
])
def test_discrete_families_equal_scipy(make, oracle, support) -> None:
    # Arrange
    distribution = make()
    k = support.astype(float)
    # Act / Assert
    np.testing.assert_allclose(distribution.pdf(k), oracle.pmf(support), rtol=1e-12)
    np.testing.assert_allclose(distribution.cdf(k), oracle.cdf(support), rtol=1e-12)
    np.testing.assert_allclose(distribution.log_pdf(k), oracle.logpmf(support), rtol=1e-12)
    assert distribution.mean == pytest.approx(float(oracle.mean()), rel=1e-12)
    assert distribution.variance == pytest.approx(float(oracle.var()), rel=1e-12)
    for p in (0.1, 0.5, 0.9):
        assert distribution.inverse_cdf(p) == float(oracle.ppf(p))


def test_binomial_log_pdf_and_named_parameters() -> None:
    # Arrange
    binomial = d.Binomial(0.3, 20)
    # Act / Assert
    assert binomial.log_pdf(5.0) == pytest.approx(stats.binom(20, 0.3).logpmf(5), rel=1e-13)
    assert (binomial.probability, binomial.number_of_trials, binomial.complement) == (0.3, 20, 0.7)
    binomial.number_of_trials = 30
    assert binomial.parameters.tolist() == [0.3, 30.0]


@pytest.mark.parametrize(("make", "cdf"), [
    (lambda: d.GeneralizedLogistic(10.0, 2.0, 0.1), lambda x: 1.0 / (1.0 + np.exp(-_hosking_y(x, 10.0, 2.0, 0.1)))),
    (lambda: d.GeneralizedLogistic(10.0, 2.0, 0.0), lambda x: 1.0 / (1.0 + np.exp(-_hosking_y(x, 10.0, 2.0, 0.0)))),
    (lambda: d.GeneralizedNormal(10.0, 2.0, -0.2), lambda x: stats.norm.cdf(_hosking_y(x, 10.0, 2.0, -0.2))),
    (lambda: d.Deterministic(2.0), lambda x: (x >= 2.0).astype(float)),
])
def test_textbook_families_equal_closed_forms(make, cdf) -> None:
    # Arrange: Hosking and Wallis (1997), Appendix A.7 and A.8
    distribution = make()
    x = np.array([5.0, 8.0, 10.0, 12.5, 18.0]) if not isinstance(distribution, d.Deterministic) \
        else np.array([1.0, 2.0, 3.0])
    # Act / Assert
    np.testing.assert_allclose(distribution.cdf(x), cdf(x), rtol=1e-12, atol=1e-15)
    if not isinstance(distribution, d.Deterministic):
        np.testing.assert_allclose(distribution.inverse_cdf(distribution.cdf(x)), x, rtol=1e-10)


def test_von_mises_density_equals_scipy_and_cdf_accumulates_from_minus_pi() -> None:
    # Arrange
    distribution, oracle = d.VonMises(0.5, 2.0), stats.vonmises(2.0, loc=0.5)
    x = np.array([-2.0, 0.0, 1.0, 2.5])
    # Act
    expected = [integrate.quad(oracle.pdf, -math.pi, v, epsabs=1e-13)[0] for v in x]
    # Assert
    np.testing.assert_allclose(distribution.pdf(x), oracle.pdf(x), rtol=1e-12)
    np.testing.assert_allclose(distribution.cdf(x), expected, rtol=1e-6)


def test_percentile_pert_families_reproduce_their_percentiles() -> None:
    # Arrange
    pert = d.PertPercentile(10.0, 50.0, 100.0)
    pert_z = d.PertPercentileZ(0.1, 0.3, 0.7)
    # Act / Assert
    np.testing.assert_allclose(pert.inverse_cdf([0.05, 0.5, 0.95]), [10.0, 50.0, 100.0], rtol=1e-6)
    np.testing.assert_allclose(pert_z.inverse_cdf([0.05, 0.5, 0.95]), [0.1, 0.3, 0.7], rtol=1e-6)
    solved = pert.to_pert()
    assert isinstance(solved, d.Pert)
    np.testing.assert_allclose(solved.inverse_cdf([0.05, 0.5, 0.95]), [10.0, 50.0, 100.0], rtol=1e-6)
    pert.percentile_95th = 120.0
    pert.solve_parameters()
    assert pert.inverse_cdf(0.95) == pytest.approx(120.0, rel=1e-6)
    assert pert.min_allowable_value < pert.max_allowable_value
    assert isinstance(pert_z.to_pert(), d.Pert)


def test_competing_risks_minimum_equals_the_product_of_survivals() -> None:
    # Arrange
    weibull, exponential = d.Weibull(10.0, 2.0), d.Exponential(0.0, 8.0)
    risks = d.CompetingRisks([weibull, exponential])
    x = np.array([1.0, 3.0, 5.0, 8.0])
    # Act
    expected = 1.0 - stats.weibull_min(2.0, scale=10.0).sf(x) * stats.expon(0.0, 8.0).sf(x)
    # Assert
    np.testing.assert_allclose(risks.cdf(x), expected, rtol=1e-6)
    assert [type(c) for c in risks.distributions] == [d.Weibull, d.Exponential]
    assert risks.dependency == "Independent" and risks.minimum_of_random_variables
    incidence = risks.cumulative_incidence_functions()
    assert len(incidence) == 2 and all(isinstance(f, d.Empirical) for f in incidence)
    assert risks.generate_random_values_with_dependency(5, 3).tolist() == \
        risks.generate_random_values_with_dependency(5, 3).tolist()
    risks.set_parameters([d.Weibull(12.0, 2.0), d.Exponential(0.0, 8.0)])
    assert risks.distributions[0].scale == 12.0
    assert d.from_xml(risks.to_xml()) == risks


def test_truncated_distribution_keeps_its_base_and_compares_by_it() -> None:
    # Arrange
    truncated = d.TruncatedDistribution(d.Normal(0.0, 1.0), -1.0, 2.0)
    # Act / Assert
    assert (truncated.lower_bound, truncated.upper_bound) == (-1.0, 2.0)
    assert truncated.base_distribution == d.Normal(0.0, 1.0)
    assert truncated == d.TruncatedDistribution(d.Normal(0.0, 1.0), -1.0, 2.0)
    assert truncated != d.TruncatedDistribution(d.Logistic(0.0, 1.0), -1.0, 2.0)
    assert isinstance(d.wrap(truncated.net), d.TruncatedDistribution)
    assert "TruncatedDistribution" not in d.FAMILIES


def test_every_typed_class_defaults_to_the_numerics_default_instance() -> None:
    # Act / Assert
    for name, cls in d.FAMILIES.items():
        if name in {"Mixture", "CompetingRisks", "Empirical", "KernelDensity"}:
            continue
        default, factory = cls(), d.create_distribution(name)
        assert type(default) is type(factory) is cls, name
        assert _bits(default.parameters) == _bits(factory.parameters), name


@pytest.mark.parametrize(("make", "names"), [
    (lambda: d.GeneralizedExtremeValue(1.0, 2.0, -0.1), ("location", "scale", "shape")),
    (lambda: d.GeneralizedPareto(1.0, 2.0, 0.1), ("location", "scale", "shape")),
    (lambda: d.Gumbel(1.0, 2.0), ("location", "scale")),
    (lambda: d.KappaFour(1.0, 2.0, 0.1, 0.2), ("location", "scale", "shape", "shape2")),
    (lambda: d.Normal(1.0, 2.0), ("mean", "standard_deviation")),
    (lambda: d.LogNormal(1.0, 0.2), ("mean_of_log", "standard_deviation_of_log")),
    (lambda: d.PearsonTypeIII(10.0, 2.0, 0.5), ("mean", "standard_deviation", "skew")),
    (lambda: d.LogPearsonTypeIII(1.0, 0.2, 0.3), ("mean_of_log", "standard_deviation_of_log", "skew_of_log")),
    (lambda: d.Uniform(1.0, 2.0), ("minimum", "maximum")),
    (lambda: d.Beta(2.0, 3.0), ("alpha", "beta")),
    (lambda: d.Poisson(2.0), ("rate",)),
    (lambda: d.GammaDistribution(2.0, 3.0), ("scale", "shape")),
    (lambda: d.Weibull(2.0, 3.0), ("scale", "shape")),
    (lambda: d.StudentT(1.0, 2.0, 5.0), ("location", "scale", "degrees_of_freedom")),
    (lambda: d.TruncatedNormal(1.0, 2.0, 0.0, 3.0), ("mu", "sigma", "minimum", "maximum")),
    (lambda: d.Triangular(0.0, 1.0, 3.0), ("minimum", "most_likely", "maximum")),
    (lambda: d.VonMises(0.5, 2.0), ("mu", "kappa")),
])
def test_named_parameter_properties_read_and_set_the_parameters(make, names) -> None:
    # Arrange
    distribution = make()
    # Act / Assert: the names read GetParameters in order, and a setter writes the parameter
    assert [getattr(distribution, n) for n in names] == distribution.parameters.tolist()
    first = names[0]
    value = getattr(distribution, first) * 1.25
    setattr(distribution, first, value)
    assert distribution.parameters[0] == value
    frozen = distribution.frozen()
    with pytest.raises(AttributeError):
        setattr(frozen, first, value)


def test_log_families_take_a_base_that_joins_their_value() -> None:
    # Arrange
    base10, natural = d.LogNormal(1.0, 0.5), d.LogNormal(1.0, 0.5, base=math.e)
    lp3 = d.LogPearsonTypeIII(1.0, 0.3, 0.4, base=math.e)
    # Act / Assert
    assert base10.base == 10.0 and natural.base == math.e
    assert base10 != natural and hash(natural.frozen()) == hash(d.LogNormal(1.0, 0.5, base=math.e).frozen())
    assert natural.cdf(math.e) == pytest.approx(stats.lognorm(0.5, scale=math.e).cdf(math.e), rel=1e-12)
    expected = stats.pearson3(0.4, 1.0, 0.3).cdf(np.log([2.0, 3.0, 5.0]))
    np.testing.assert_allclose(lp3.cdf([2.0, 3.0, 5.0]), expected, rtol=1e-10)
    with pytest.raises(ValueError, match="greater than one"):
        d.LogNormal(1.0, 0.5, base=1.0)
    with pytest.raises(ValueError, match="greater than one"):
        natural.base = 0.5


def test_pearson_gamma_form_and_log_space() -> None:
    # Arrange: xi = mu - 2 sigma / gamma, beta = sigma gamma / 2, alpha = 4 / gamma^2
    p3 = d.PearsonTypeIII(100.0, 20.0, 0.8)
    lp3 = d.LogPearsonTypeIII(2.0, 0.25, -0.5)
    # Act / Assert
    assert p3.support_bound == pytest.approx(100.0 - 2 * 20.0 / 0.8, rel=1e-14)
    assert p3.scale == pytest.approx(20.0 * 0.8 / 2, rel=1e-14)
    assert p3.shape == pytest.approx(4 / 0.8 ** 2, rel=1e-14)
    assert lp3.support_bound == pytest.approx(2.0 + 2 * 0.25 / 0.5, rel=1e-14)
    log_space = lp3.log_space()
    assert isinstance(log_space, d.PearsonTypeIII) and log_space.parameters.tolist() == [2.0, 0.25, -0.5]
    assert lp3.cdf(100.0) == pytest.approx(log_space.cdf(2.0), rel=1e-12)


def _quad_linear_moments(ppf) -> np.ndarray:
    """Return [λ1, λ2, τ3, τ4] by integrating the quantile function against shifted Legendre polynomials."""
    weights = [lambda u: 1.0, lambda u: 2 * u - 1, lambda u: 6 * u * u - 6 * u + 1,
               lambda u: 20 * u ** 3 - 30 * u * u + 12 * u - 1]
    lam = [integrate.quad(lambda u, w=w: ppf(u) * w(u), 0, 1, epsabs=1e-12, epsrel=1e-12, limit=200)[0]
           for w in weights]
    return np.array([lam[0], lam[1], lam[2] / lam[1], lam[3] / lam[1]])


@pytest.mark.parametrize(("cls", "parameters", "oracle", "tolerance"), [
    (d.GeneralizedPareto, [10.0, 4.0, 0.15], stats.genpareto(-0.15, 10.0, 4.0), 1e-10),
    (d.GeneralizedExtremeValue, [100.0, 15.0, -0.2], stats.genextreme(-0.2, 100.0, 15.0), 2e-3),
    (d.KappaFour, [10.0, 2.0, 0.1, 0.2], stats.kappa4(0.2, 0.1, 10.0, 2.0), 1e-6),
])
def test_parameters_from_linear_moments_on_the_class_recover_the_parameters(cls, parameters, oracle,
                                                                            tolerance) -> None:
    # Arrange: L-moments of the scipy quantile function (GEV uses Hosking's 9e-4 approximation for |τ3| <= 0.5)
    moments = _quad_linear_moments(oracle.ppf)
    # Act
    fitted = cls.parameters_from_linear_moments(moments)
    back = cls.linear_moments_from_parameters(parameters)
    # Assert
    np.testing.assert_allclose(fitted, parameters, rtol=tolerance, atol=tolerance)
    np.testing.assert_allclose(back, moments, rtol=1e-8)
    np.testing.assert_array_equal(cls(*parameters).parameters_from_linear_moments(moments), fitted)


def test_parameters_from_moments_equal_scipy_moments() -> None:
    # Arrange
    oracle = stats.gumbel_r(50.0, 8.0)
    mean, variance, skew, excess = oracle.stats(moments="mvsk")
    # Act
    moments = d.Gumbel.moments_from_parameters([50.0, 8.0])
    parameters = d.Gumbel.parameters_from_moments([mean, math.sqrt(variance)])
    # Assert
    np.testing.assert_allclose(moments, [mean, math.sqrt(variance), skew, excess + 3], rtol=1e-12)
    np.testing.assert_allclose(parameters, [50.0, 8.0], rtol=1e-12)
    np.testing.assert_allclose(d.GeneralizedLogistic.moments_from_parameters([10.0, 2.0, 0.0])[:2],
                               [10.0, 2.0 * math.pi / math.sqrt(3)], rtol=1e-12)
    with pytest.raises(NotImplementedError, match="IMomentEstimation"):
        d.KappaFour.parameters_from_moments([1.0, 1.0, 0.0])
    with pytest.raises(TypeError):
        d.UnivariateDistribution.parameters_from_moments([1.0, 1.0])


def test_maximum_likelihood_and_constraints_equal_scipy_fits() -> None:
    # Arrange
    sample = stats.gumbel_r(50.0, 8.0).rvs(size=200, random_state=np.random.default_rng(5))
    gumbel = d.Gumbel()
    # Act
    fitted = gumbel.mle(sample)
    initial, lower, upper = gumbel.parameter_constraints(sample)
    # Assert
    reference = stats.gumbel_r.fit(sample)
    np.testing.assert_allclose(fitted, reference, rtol=2e-4)
    assert d.Gumbel(*fitted).log_likelihood(sample) >= d.Gumbel(*reference).log_likelihood(sample) - 1e-5
    assert gumbel.parameters.tolist() == [100.0, 10.0]  # unchanged
    assert initial.shape == lower.shape == upper.shape == (2,)
    assert np.all(lower <= initial) and np.all(initial <= upper)
    gumbel.set_parameters_from_mle(sample)
    np.testing.assert_allclose(gumbel.parameters, reference, rtol=1e-6)
    expon = d.Exponential().mle(sample)
    np.testing.assert_allclose(expon, stats.expon.fit(sample), rtol=5e-4)  # closed form: min and mean - min
    gpa = d.GeneralizedPareto()
    np.testing.assert_array_equal(gpa.mle(sample, location=20.0),
                                  d.GeneralizedPareto.fit_known_location(sample, 20.0).parameters)


def test_bootstrap_refits_a_reproducible_replicate() -> None:
    # Arrange
    gev = d.GeneralizedExtremeValue(100.0, 10.0, -0.1)
    # Act
    first = gev.bootstrap("maximum_likelihood", 60, seed=12)
    second = gev.bootstrap(enums.EstimationMethod.MAXIMUM_LIKELIHOOD, 60, seed=12)
    # Assert
    assert isinstance(first, d.GeneralizedExtremeValue) and first == second and first != gev
    with pytest.raises(NotImplementedError, match="IBootstrappable"):
        d.Bernoulli().bootstrap("method_of_moments", 10)


def test_standard_error_members_follow_the_delta_method() -> None:
    # Arrange
    gev = d.GeneralizedExtremeValue(100.0, 10.0, -0.1)
    p, n, step = 0.99, 80, 1e-6
    # Act
    gradient = gev.quantile_gradient(p)
    variance = gev.quantile_variance(p, n, "maximum_likelihood")
    jacobian, determinant = gev.quantile_jacobian([0.1, 0.5, 0.9])
    numeric = []
    for i in range(3):
        up, down = gev.parameters.copy(), gev.parameters.copy()
        up[i] += step
        down[i] -= step
        numeric.append((d.GeneralizedExtremeValue(*up).inverse_cdf(p) - d.GeneralizedExtremeValue(*down).inverse_cdf(p))
                       / (2 * step))
    # Assert
    np.testing.assert_allclose(gradient, numeric, rtol=1e-6)
    covariance = gev.parameter_covariance(n, "maximum_likelihood")
    assert variance == pytest.approx(gradient @ covariance @ gradient, rel=1e-10)
    np.testing.assert_allclose(jacobian[1], gev.quantile_gradient(0.5), rtol=1e-12)
    assert determinant == pytest.approx(np.linalg.det(jacobian), rel=1e-10)
    assert gev.log_abs_quantile_jacobian([0.1, 0.5, 0.9]) == pytest.approx(math.log(abs(determinant)), rel=1e-12)
    with pytest.raises(NotImplementedError, match="IStandardError"):
        d.Pareto().quantile_gradient(0.5)


GPA_SURVIVAL = [  # (u, sigma, xi hydrologic, q, reference bits, inverse_cdf(1 - q) bits): Stochastics vectors
    (3.681482536940442, 2.9923534363268085, -0.39815305233042264, 2.991269027366455e-06,
     "40264C96E312E3AB", "40264C96E312E3B4"),
    (3.064228242645925, 0.7161935809731612, 0.0, 6.485714300886042e-11, "4033DD83F2B0AFCE", "4033DD83F7035DBE"),
    (0.8234946457607617, 2.568463976798949, 0.30614590243611783, 4.218502413991999e-11,
     "40C89244C44F05B1", "40C89244A90DEE0A"),
    (2.692607598820238, 0.17804162514449712, 0.3667959761699061, 1.7511958387651114e-12,
     "40C376A75F5BDE44", "40C376B257423404"),
]


def _from_hex(text: str) -> float:
    return float(np.array([int(text, 16)], dtype=np.uint64).view(np.float64)[0])


@pytest.mark.parametrize(("u", "sigma", "xi", "q", "reference", "complement"), GPA_SURVIVAL)
def test_inverse_ccdf_keeps_survival_precision_for_the_generalized_pareto(u, sigma, xi, q, reference,
                                                                          complement) -> None:
    # Arrange: the hydrologic shape xi is Hosking's -kappa
    gpa = d.GeneralizedPareto(u, sigma, -xi)
    # Act
    quantile = gpa.inverse_ccdf(q)
    # Assert
    assert quantile == _from_hex(reference)
    assert gpa.inverse_cdf(1.0 - q) == _from_hex(complement)
    assert quantile == pytest.approx(float(stats.genpareto(xi, u, sigma).isf(q)), rel=1e-13)


def test_inverse_ccdf_array_equals_the_scalar_member_and_scipy() -> None:
    # Arrange
    gpa = d.GeneralizedPareto(1.0, 2.0, -0.2)
    q = np.geomspace(1e-14, 0.9, 20000)
    # Act
    array = gpa.inverse_ccdf(q)
    loop = [float(gpa.net.InverseCCDF(float(v))) for v in q]
    # Assert
    assert _bits(array) == _bits(loop)
    np.testing.assert_allclose(array, stats.genpareto(0.2, 1.0, 2.0).isf(q), rtol=1e-12)
    assert d.Normal().inverse_ccdf(0.025) == d.Normal().inverse_cdf(0.975)
    assert gpa.inverse_ccdf(np.array([[0.1, 0.2]])).shape == (1, 2)


@pytest.mark.benchmark
def test_inverse_ccdf_array_crosses_once_faster_than_scalar_calls() -> None:
    # Arrange
    gpa = d.GeneralizedPareto(1.0, 2.0, -0.2)
    q = np.geomspace(1e-14, 0.9, 20000)
    # Act
    start = time.perf_counter()
    gpa.inverse_ccdf(q)
    vector_time = time.perf_counter() - start
    scalar = gpa.net.InverseCCDF
    start = time.perf_counter()
    _ = [float(scalar(float(v))) for v in q]
    loop_time = time.perf_counter() - start
    # Assert
    assert vector_time * 5 < loop_time, (vector_time, loop_time)


POISSON_GPD = [  # (rate, u, sigma, xi hydrologic, location, scale): Stochastics reference vectors
    (10.0, 0.5, 0.25, 0.0, 1.0756462732485115, 0.25),
    (1.0, 2.0, 1.5, 0.2, 2.0, 1.5),
    (0.3, 2.0, 1.5, 0.2, 0.39502314197467103, 1.1790046283949343),
    (3.7, 10.0, 2.0, -0.15, 12.375909058299941, 1.6436136412550089),
    (250.0, 0.0, 0.8, 0.45, 19.550239167436285, 9.597607625346328),
    (1e-06, 1.0, 1.0, 0.1, -6.488113568490419, 0.251188643150958),
]


@pytest.mark.parametrize(("rate", "u", "sigma", "xi", "location", "scale"), POISSON_GPD)
def test_gev_from_poisson_generalized_pareto_equals_the_reference_vectors(rate, u, sigma, xi, location,
                                                                          scale) -> None:
    # Act
    gev = d.GeneralizedExtremeValue.from_poisson_generalized_pareto(rate, u, sigma, -xi)
    # Assert: the expm1 location agrees with the Pow form to rounding; the scale and shape are exact
    assert isinstance(gev, d.GeneralizedExtremeValue)
    assert gev.location == pytest.approx(location, rel=1e-14, abs=1e-15)
    assert gev.scale == pytest.approx(scale, rel=1e-15)
    assert gev.shape == -xi


def test_gev_from_poisson_generalized_pareto_is_the_annual_maximum_and_refuses_bad_rates() -> None:
    # Arrange: P(max <= x) = exp(-rate * S_gpd(x))
    rate, gpd = 3.7, d.GeneralizedPareto(10.0, 2.0, -0.15)
    gev = d.GeneralizedExtremeValue.from_poisson_generalized_pareto(rate, 10.0, 2.0, -0.15)
    x = np.array([12.0, 15.0, 25.0])
    # Act / Assert
    np.testing.assert_allclose(gev.cdf(x), np.exp(-rate * gpd.ccdf(x)), rtol=1e-12)
    with pytest.raises(ValueError):
        d.GeneralizedExtremeValue.from_poisson_generalized_pareto(0.0, 10.0, 2.0, 0.1)


POISSON_GEV = [  # (rate, location, scale, shape, GPD location, GPD scale): mpmath at 400 bits
    (3.7, 12.8911, 2.4337, -0.15, 9.999975062407715, 2.0000312593611573),
    (0.3, 2.0, 1.5, 0.2, 3.604976858025329, 1.1790046283949343),
    (250.0, 19.5, 9.6, -0.45, -0.05511241277857893, 0.8001994142496389),
    (1e-06, -6.5, 0.25, -0.1, 0.9526792638374316, 0.9952679263837432),
    (10.0, 1.0756, 0.25, 0.0, 0.49995372675148847, 0.25),
    (10.0, 1.0756, 0.25, 1e-12, 0.4999537267508257, 0.25000000000057565),
    (10.0, 1.0756, 0.25, -1e-12, 0.4999537267521512, 0.24999999999942435),
    (1e6, 3.5, 0.004, 0.4, 0.9981135684904192, 1.0047545726038323),
]


@pytest.mark.parametrize(("rate", "location", "scale", "shape", "gpd_location", "gpd_scale"), POISSON_GEV)
def test_gpd_from_poisson_generalized_extreme_value_equals_mpmath(rate, location, scale, shape, gpd_location,
                                                                  gpd_scale) -> None:
    # Act
    gpd = d.GeneralizedPareto.from_poisson_generalized_extreme_value(rate, location, scale, shape)
    # Assert: within 4 ulp of the terms combined; the shape is kept exactly
    unit = 2.0 ** -53
    tolerance = 4 * unit * (abs(location) + abs(location - gpd_location))
    assert isinstance(gpd, d.GeneralizedPareto)
    assert gpd.location == pytest.approx(gpd_location, rel=0, abs=tolerance)
    assert gpd.scale == pytest.approx(gpd_scale, rel=4 * unit)
    assert gpd.shape == shape


@pytest.mark.parametrize("rate", [1e-6, 0.3, 1.0, 3.7, 250.0, 1e6])
@pytest.mark.parametrize("shape", [-0.45, -1e-12, 0.0, 1e-12, 0.4])
def test_gpd_from_poisson_generalized_extreme_value_inverts_the_forward_mapping(rate, shape) -> None:
    # Act
    gev = d.GeneralizedExtremeValue.from_poisson_generalized_pareto(rate, 10.0, 2.0, shape)
    gpd = d.GeneralizedPareto.from_poisson_generalized_extreme_value(rate, gev.location, gev.scale, gev.shape)
    # Assert
    assert gpd.location == pytest.approx(10.0, rel=0, abs=1e-14 * (10.0 + 2 * abs(gev.location - 10.0)))
    assert gpd.scale == pytest.approx(2.0, rel=1e-14)
    assert gpd.shape == shape


def test_gpd_from_poisson_generalized_extreme_value_is_the_annual_maximum_and_refuses_bad_input() -> None:
    # Arrange: P(max <= x) = exp(-rate * S_gpd(x)) above the threshold
    rate, gev = 3.7, d.GeneralizedExtremeValue(12.0, 2.5, -0.15)
    gpd = d.GeneralizedPareto.from_poisson_generalized_extreme_value(rate, 12.0, 2.5, -0.15)
    x = gpd.inverse_cdf(np.array([0.1, 0.5, 0.99]))
    # Act / Assert
    np.testing.assert_allclose(gev.cdf(x), np.exp(-rate * gpd.ccdf(x)), rtol=1e-12)
    assert gev.cdf(gpd.location) == pytest.approx(math.exp(-rate), rel=1e-13)
    invalid = [(0.0, 1.0, 1.0, 0.1), (2.0, math.nan, 1.0, 0.1), (2.0, 1.0, -1.0, 0.1), (2.0, 1.0, 1.0, math.inf)]
    for arguments in invalid:
        with pytest.raises(ValueError):
            d.GeneralizedPareto.from_poisson_generalized_extreme_value(*arguments)


def test_parameter_set_quantiles_raise_the_row_error_when_validating() -> None:
    # Arrange
    gev = d.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    rows = np.array([[4.0, 1.0, -0.1], [5.0, 2.0, 0.1], [4.0, -1.0, 0.0]])
    # Act / Assert
    with pytest.raises(ValueError, match=r"^Parameter set 2: The scale parameter α \(alpha\) must be positive") as info:
        gev.parameter_set_quantiles(rows, 0.99, validate=True)
    assert isinstance(info.value.__cause__, ValueError)
    with pytest.raises(ValueError, match="^Parameter set 2: "):
        gev.parameter_set_quantile_table(rows, [0.5, 0.99], validate=True)
    quantiles, valid = gev.parameter_set_quantile_table(rows, [0.5, 0.99])
    assert valid.tolist() == [True, True, False] and np.isnan(quantiles[2]).all()
    np.testing.assert_array_equal(gev.parameter_set_quantiles(rows[:2], 0.99, validate=True),
                                  gev.parameter_set_quantiles(rows[:2], 0.99))


def test_base_members_equal_scipy_and_their_definitions() -> None:
    # Arrange
    gumbel, oracle = d.Gumbel(50.0, 8.0), stats.gumbel_r(50.0, 8.0)
    x = np.array([40.0, 55.0, 80.0])
    # Act / Assert
    np.testing.assert_allclose(gumbel.hf(x), oracle.pdf(x) / oracle.sf(x), rtol=1e-12)
    assert gumbel.mode == 50.0
    assert gumbel.coefficient_of_variation == pytest.approx(gumbel.standard_deviation / gumbel.mean, rel=1e-15)
    np.testing.assert_allclose(gumbel.central_moments(), [gumbel.mean, gumbel.standard_deviation, gumbel.skewness,
                                                          gumbel.kurtosis], rtol=1e-6)
    np.testing.assert_allclose(gumbel.central_moments(steps=2000)[:2], [gumbel.mean, gumbel.standard_deviation],
                               rtol=1e-3)
    mass = oracle.cdf(70.0) - oracle.cdf(45.0)
    assert gumbel.conditional_mean(45.0, 70.0) == pytest.approx(
        oracle.expect(lambda v: v, lb=45.0, ub=70.0) / mass, rel=1e-8)
    assert gumbel.conditional_variance(45.0, 70.0, mean=55.0) == pytest.approx(
        oracle.expect(lambda v: (v - 55.0) ** 2, lb=45.0, ub=70.0) / mass, rel=1e-7)
    assert gumbel.conditional_skewness(45.0, 70.0, mean=55.0) == pytest.approx(
        oracle.expect(lambda v: (v - 55.0) ** 3, lb=45.0, ub=70.0) / mass, rel=1e-6)
    assert gumbel.conditional_kurtosis(45.0, 70.0, mean=55.0) == pytest.approx(
        oracle.expect(lambda v: (v - 55.0) ** 4, lb=45.0, ub=70.0) / mass, rel=1e-7)
    assert gumbel.conditional_moments(45.0, 70.0)[0] == pytest.approx(gumbel.conditional_mean(45.0, 70.0), rel=1e-3)
    tail = oracle.expect(lambda v: v, lb=oracle.ppf(0.9)) / 0.1
    assert gumbel.conditional_expected_value(0.9) == pytest.approx(tail, rel=1e-7)
    assert gumbel.log_likelihood_left_censored(45.0, 3) == pytest.approx(3 * oracle.logcdf(45.0), rel=1e-12)
    assert gumbel.log_likelihood_right_censored(70.0, 2) == pytest.approx(2 * oracle.logsf(70.0), rel=1e-12)
    assert gumbel.log_likelihood_intervals(45.0, 70.0) == pytest.approx(math.log(mass), rel=1e-12)


def test_base_descriptive_members() -> None:
    # Arrange
    gev = d.GeneralizedExtremeValue(100.0, 10.0, -0.1)
    # Act / Assert
    assert gev.short_display_name == "GEV"
    assert gev.display_label.startswith("GEV")
    assert gev.parameter_display_names == ("Location (ξ)", "Scale (α)", "Shape (κ)")
    assert gev.parameter_names_short_form == ("ξ", "α", "κ")
    assert [row[0] for row in gev.parameters_to_string] == list(gev.parameter_display_names)
    assert gev.minimum_of_parameters.shape == gev.maximum_of_parameters.shape == (3,)
    assert np.all(gev.minimum_of_parameters <= gev.parameters) and np.all(gev.parameters <= gev.maximum_of_parameters)
    assert gev.compare_to(d.GeneralizedExtremeValue(200.0, 10.0, -0.1)) < 0 and gev.compare_to(None) > 0


def test_graph_tables_follow_the_distribution() -> None:
    # Arrange
    normal = d.Normal(0.0, 1.0)
    bins = sampling.stratify_x_values(-2.0, 2.0, 8)
    probabilities = sampling.stratify_probabilities(0.1, 0.9, 4)
    # Act
    pdf, cdf = normal.create_pdf_graph(), normal.create_cdf_graph(bins)
    quantiles = normal.create_inverse_cdf_graph([probabilities, probabilities])
    # Assert
    np.testing.assert_allclose(pdf[:, 1], stats.norm.pdf(pdf[:, 0]), rtol=1e-12)
    np.testing.assert_array_equal(cdf[:, 0], bins.midpoint)
    np.testing.assert_allclose(cdf[:, 1], stats.norm.cdf(bins.midpoint), rtol=1e-12)
    assert quantiles.shape == (8, 2)
    np.testing.assert_allclose(quantiles[:4, 0], stats.norm.ppf(quantiles[:4, 1]), rtol=1e-12)


def test_xml_round_trip_and_factory_members() -> None:
    # Arrange
    families = [d.GammaDistribution(2.0, 3.0), d.Pert(1.0, 2.0, 5.0), d.LogNormal(1.0, 0.2, base=math.e),
                d.Mixture([0.4, 0.6], [d.Normal(), d.Gumbel(1.0, 2.0)]), d.PertPercentile(10.0, 50.0, 100.0),
                d.Empirical([1.0, 2.0, 3.0], [0.1, 0.5, 0.9]), d.KernelDensity([1.0, 2.0, 4.0, 7.0])]
    # Act / Assert
    for distribution in families:
        assert d.from_xml(distribution.to_xml()) == distribution, type(distribution).__name__
    assert d.try_create_distribution("Mixture") is None
    assert d.try_create_distribution("weibull") == d.Weibull()
    with pytest.raises(ValueError):
        d.try_create_distribution("NotAFamily")
    with pytest.raises(ValueError):
        d.from_xml("<Distribution />")


def test_gamma_members_equal_the_pearson_frequency_factor() -> None:
    # Arrange
    gamma = d.GammaDistribution(2.0, 4.0)
    sample = stats.gamma(4.0, scale=2.0).rvs(size=300, random_state=np.random.default_rng(9))
    # Act / Assert
    assert gamma.rate == 0.5
    for skew, p in ((0.5, 0.99), (-0.3, 0.1)):
        expected = float(stats.pearson3(skew).ppf(p))  # the Cornish-Fisher approximation is within 1e-4
        assert d.GammaDistribution.frequency_factor_kp(skew, p) == pytest.approx(expected, rel=1e-4)
        step = 1e-5
        numeric = (stats.pearson3(skew + step).ppf(p) - stats.pearson3(skew - step).ppf(p)) / (2 * step)
        assert d.GammaDistribution.partial_kp(skew, p) == pytest.approx(numeric, rel=1e-2)
    assert gamma.wilson_hilferty_inverse_cdf(0.9) == pytest.approx(gamma.inverse_cdf(0.9), rel=1e-2)
    shape, _, scale = stats.gamma.fit(sample, floc=0.0)
    gamma.estimate_mle_newton_raphson(sample)
    assert gamma.shape == pytest.approx(shape, rel=1e-5)
    gamma.estimate_mle_bobee(sample)
    np.testing.assert_allclose(gamma.parameters, [scale, shape], rtol=1e-3)
    with pytest.raises(TypeError, match="frozen"):
        gamma.frozen().estimate_mle_bobee(sample)


def test_gamma_newton_raphson_scale_is_the_mean_over_the_fitted_shape() -> None:
    # Arrange
    sample = stats.gamma(4.0, scale=2.0).rvs(size=300, random_state=np.random.default_rng(9))
    gamma = d.GammaDistribution(2.0, 4.0)
    shape, _, scale = stats.gamma.fit(sample, floc=0.0, method="MLE")
    # Act
    gamma.estimate_mle_newton_raphson(sample)
    # Assert
    assert gamma.scale == pytest.approx(sample.mean() / gamma.shape, rel=1e-12)
    np.testing.assert_allclose([gamma.scale, gamma.shape], [scale, shape], rtol=1e-6)


def test_family_specific_members() -> None:
    # Arrange
    beta = d.GeneralizedBeta.pert(10.0, 30.0, 100.0)
    p3 = d.PearsonTypeIII(100.0, 20.0, 0.5)
    sample = stats.lognorm(0.4, scale=50.0).rvs(size=200, random_state=np.random.default_rng(2))
    # Act / Assert
    np.testing.assert_allclose(beta.parameters, [1 + 4 * 20 / 90, 1 + 4 * 70 / 90, 10.0, 100.0], rtol=1e-14)
    np.testing.assert_allclose(d.GeneralizedBeta.pert(10.0, 30.0, 100.0, 4.0).parameters, beta.parameters)
    beta.set_parameters_from_moments(0.5, 0.1, 0.0, 1.0)
    assert (beta.mean, beta.standard_deviation) == pytest.approx((0.5, 0.1), rel=1e-12)
    assert d.Pert(10.0, 30.0, 100.0).beta.parameters.tolist() == pytest.approx([1 + 80 / 90, 1 + 280 / 90, 10, 100])
    assert p3.wilson_hilferty_inverse_cdf(0.99) == pytest.approx(p3.inverse_cdf(0.99), rel=1e-2)
    assert p3.quantile_gradient_for_moments(0.99).shape == (3,)
    ln, log10 = np.log(sample), np.log10(sample)
    np.testing.assert_allclose(d.LnNormal().indirect_method_of_moments(sample)[:3],
                               [ln.mean(), ln.std(ddof=1), stats.skew(ln, bias=False)], rtol=1e-10)
    np.testing.assert_allclose(d.LogNormal().indirect_method_of_moments(sample)[:2],
                               [log10.mean(), log10.std(ddof=1)], rtol=1e-12)
    l1 = d.LogPearsonTypeIII().indirect_method_of_linear_moments(sample)
    assert l1[0] == pytest.approx(log10.mean(), rel=1e-12)
    np.testing.assert_allclose(d.LnNormal().indirect_method_of_linear_moments(sample)[:2],
                               [ln.mean(), l1[1] * math.log(10.0)], rtol=1e-12)
    mu, sigma = d.LogNormal().direct_method_of_moments(100.0, 30.0)
    assert d.LogNormal(mu, sigma).mean == pytest.approx(100.0, rel=1e-12)
    gev = d.GeneralizedExtremeValue(100.0, 10.0, -0.1)
    assert gev.solve_for_kappa(gev.skewness) == pytest.approx(-0.1, rel=1e-3)
    np.testing.assert_allclose(gev.direct_method_of_moments([gev.mean, gev.standard_deviation, gev.skewness]),
                               gev.parameters, rtol=1e-3)
    information = gev.expected_information_matrix(50)
    np.testing.assert_allclose(np.linalg.inv(information), gev.parameter_covariance(50), rtol=1e-8)
    gpa = d.GeneralizedPareto(0.0, 2.0, 0.1)
    assert gpa.peaks_per_block == 1.0 and gpa.modified_method_of_moments(sample).shape == (3,)
    assert d.Weibull().solve_mle(sample) == pytest.approx(stats.weibull_min.fit(sample, floc=0.0)[::-2], rel=1e-5)


def test_normal_interval_members_bracket_the_quantile() -> None:
    # Arrange
    normal = d.Normal(10.0, 2.0)
    # Act
    nci = normal.normal_confidence_intervals(30, [0.5, 0.99], [0.05, 0.95])
    tci = normal.noncentral_t_confidence_intervals(30, [0.5, 0.99], [0.05, 0.95])
    mci = normal.monte_carlo_confidence_intervals(30, 2000, [0.5, 0.99], [0.05, 0.95])
    # Assert: the median bound is mean ± z(0.95) σ / sqrt(n) for the normal interval
    half = stats.norm.ppf(0.95) * 2.0 / math.sqrt(30)
    np.testing.assert_allclose(nci[0], [10.0 - half, 10.0 + half], rtol=1e-10)
    k = stats.nct(29, stats.norm.ppf(0.99) * math.sqrt(30)).ppf([0.05, 0.95]) / math.sqrt(30)
    np.testing.assert_allclose(tci[1], 10.0 + 2.0 * k, rtol=1e-6)
    assert mci.shape == (2, 2) and mci[1, 0] < normal.inverse_cdf(0.99) < mci[1, 1]
    assert d.standard_log_pdf(np.array([0.0, 1.5])) == pytest.approx(stats.norm.logpdf([0.0, 1.5]), rel=1e-14)
    assert 0.99 > normal.expected_probability(30, 0.99) > 0.97


def test_table_family_members() -> None:
    # Arrange
    p = np.linspace(0.01, 0.99, 25)
    a = d.Empirical(stats.norm(5.0, 1.0).ppf(p), p)
    mixture = d.Mixture([0.3, 0.7], [d.Normal(), d.Normal(3.0, 1.0)])
    kde = d.KernelDensity([1.0, 2.0, 4.0])
    # Act
    total = d.Empirical.convolve([a, a])
    three = d.Empirical.convolve([a, a, a], 512)
    values, masses = d.Empirical.convolve_discrete([0.0, 1.0], [0.5, 0.5], [0.0, 1.0], [0.5, 0.5], 3)
    # Assert
    assert total.mean == pytest.approx(2 * a.mean, rel=1e-3) and three.mean == pytest.approx(3 * a.mean, rel=1e-3)
    assert d.Empirical.convolve([a, a], 256, True).mean == pytest.approx(2 * a.mean, rel=1e-2)
    with pytest.raises(ValueError, match="exactly two"):
        d.Empirical.convolve([a, a, a], log_spaced_output=True)
    assert masses.sum() == pytest.approx(1.0, rel=1e-12) and values.size == masses.size
    assert values @ masses == pytest.approx(1.0, rel=1e-12)
    assert a.x_value_order == "Ascending" and a.strict_probabilities
    assert isinstance(a.strict_x_values, bool)
    assert total.mean == pytest.approx(10.0, rel=1e-2)
    assert not mixture.is_zero_inflated and mixture.x_transform == "None"
    mixture.max_iterations, mixture.tolerance = 50, 1e-6
    assert (mixture.max_iterations, mixture.tolerance, mixture.iterations) == (50, 1e-6, 0)
    mixture.create_empirical_cdf()
    kde.set_sample_data([1.0, 2.0, 4.0, 8.0], [1.0, 1.0, 1.0, 2.0])
    assert kde.sample_size == 4
    assert d.KernelDensity.degenerate_relative_bandwidth == 1e-9 == kde.degenerate_absolute_bandwidth


def test_bootstrap_analysis_and_results_members() -> None:
    # Arrange
    gev = d.GeneralizedExtremeValue(100.0, 10.0, -0.1)
    analysis = bootstrap.BootstrapAnalysis(gev, "maximum_likelihood", 50, replications=200, seed=7)
    sample = gev.generate_random_values(50, 3)
    # Act
    results = analysis.estimate([0.5, 0.9, 0.99])
    values, fitness, _ = analysis.parameter_sets()
    # Assert
    assert analysis.distribution == gev and analysis.estimation_method == "MaximumLikelihood"
    assert (analysis.sample_size, analysis.replications, analysis.prng_seed) == (50, 200, 7)
    np.testing.assert_array_equal(values, analysis.parameters())
    assert np.isnan(fitness).all()
    assert analysis.product_moments().shape == analysis.linear_moments().shape == (200, 4)
    assert analysis.probabilities([120.0]).shape == (200, 1)
    low, high = analysis.compute_min_max_quantiles(0.01, 0.99)
    assert low < gev.inverse_cdf(0.01) < gev.inverse_cdf(0.99) < high
    for interval in (analysis.normal_quantile_ci([0.99]), analysis.bias_corrected_quantile_ci([0.99]),
                     analysis.bca_quantile_ci(sample, [0.99])):
        assert interval.shape == (1, 2) and interval[0, 0] < interval[0, 1]
    assert analysis.bootstrap_t_quantile_ci([0.99], 0.1).shape == (1, 2)
    assert results.parent_distribution == gev and results.parameter_sets[0].shape == (200, 3)
    results.aic, results.erl = 12.5, 40.0
    restored = bootstrap.UncertaintyAnalysisResults.from_bytes(results.to_bytes())
    np.testing.assert_array_equal(restored.mean_curve, results.mean_curve)
    assert restored.aic == 12.5 and restored.erl == 40.0
    from_xml = bootstrap.UncertaintyAnalysisResults.from_xml(results.to_xml())
    np.testing.assert_array_equal(from_xml.confidence_intervals, results.confidence_intervals)
    empty = bootstrap.UncertaintyAnalysisResults()
    assert empty.parameter_sets is None
    empty.process_parameter_sets([gev, None])
    assert np.isnan(empty.parameter_sets[0][1]).all()
    empty.parent_distribution = gev
    assert empty.parent_distribution == gev


def test_new_families_freeze_copy_and_wrap() -> None:
    # Arrange
    weibull = d.Weibull(3.0, 1.5)
    # Act
    frozen = weibull.frozen()
    # Assert
    assert isinstance(frozen, d.Weibull) and copy.deepcopy(frozen) is frozen
    assert d.FAMILIES["GammaDistribution"] is d.GammaDistribution
    assert isinstance(d.create_distribution("TruncatedNormal"), d.TruncatedNormal)
    assert d.GammaDistribution.parameter_names == ("Theta", "Kappa")
    with pytest.raises(TypeError, match="frozen"):
        frozen.set_parameters([1.0, 2.0])


def test_normal_lower_tail_matches_mpmath_to_rounding() -> None:
    # Arrange: mpmath ncdf at 40 digits
    z = np.array([-5.0, -10.0, -20.0, -37.0])
    expected = np.array([2.866515718791939e-07, 7.619853024160525e-24, 2.7536241186062337e-89, 5.725571222524577e-300])
    # Act / Assert
    np.testing.assert_allclose(d.standard_cdf(z), expected, rtol=1e-15, atol=0)
    assert d.standard_cdf(0.0) == 0.5
    assert d.Normal(1.0, 2.0).cdf(1.0) == 0.5


def test_gamma_tails_at_large_shape_match_mpmath() -> None:
    # Arrange: mpmath Q(21544.346900318822, 26930.433625398528) at 50 digits; scipy gammaincc agrees
    gamma = d.GammaDistribution(1.0, 21544.346900318822)
    x = 26930.433625398528
    expected = 5.637596774491294e-254
    # Act / Assert
    assert gamma.ccdf(x) == pytest.approx(expected, rel=2e-13)
    assert gamma.log_ccdf(x) == pytest.approx(math.log(expected), rel=1e-15)
    assert gamma.ccdf(x) == pytest.approx(float(stats.gamma(21544.346900318822).sf(x)), rel=1e-12)
    # Shapes far beyond the series range stay finite and match scipy near the median
    huge = d.GammaDistribution(1.0, 4e12)
    assert huge.cdf(3999994000000.0) == pytest.approx(float(stats.gamma(4e12).cdf(3999994000000.0)), rel=1e-9)
