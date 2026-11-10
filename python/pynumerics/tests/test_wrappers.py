"""Typed wrappers call Numerics and return the values Numerics computes."""

from __future__ import annotations

import math

import numpy as np
import pytest

from pynumerics import (
    bootstrap,
    copulas,
    data,
    distributions,
    functions,
    linalg,
    mathematics,
    multivariate,
    optimization,
    roots,
    sampling,
    special,
    statistics,
)


def test_standard_z_known_values_and_array_overload() -> None:
    # Arrange
    p = np.array([0.5, 0.975, 0.025, 0.999])
    # Act
    z = distributions.standard_z(p)
    # Assert
    np.testing.assert_allclose(z, [0.0, 1.959963984540054, -1.959963984540054, 3.090232306167813], atol=1e-12)
    assert z.shape == p.shape
    for value, expected in zip(p, z, strict=True):
        assert distributions.standard_z(float(value)) == expected  # scalar and array paths agree bitwise
    np.testing.assert_allclose(distributions.standard_cdf(z), p, atol=1e-14)


def test_stratify_value_and_probability_bins() -> None:
    # Act
    xbins = sampling.stratify_x_values(0.0, 10.0, 5)
    pbins = sampling.stratify_probabilities(0.001, 0.999, 4, importance="uniform", is_exhaustive=False)
    # Assert
    np.testing.assert_allclose(xbins.lower, [0, 2, 4, 6, 8])
    np.testing.assert_allclose(xbins.upper, [2, 4, 6, 8, 10])
    np.testing.assert_allclose(xbins.midpoint, [1, 3, 5, 7, 9])
    assert len(pbins) == 4
    np.testing.assert_allclose(pbins.lower[0], 0.001)
    np.testing.assert_allclose(pbins.upper[-1], 0.999)
    np.testing.assert_allclose(pbins.weight, pbins.upper - pbins.lower, atol=1e-15)


def test_latin_hypercube_one_point_per_stratum_and_reproducible() -> None:
    # Act
    design = sampling.latin_hypercube(10, dimension=3, seed=12345)
    again = sampling.latin_hypercube(10, dimension=3, seed=12345)
    # Assert
    assert design.shape == (10, 3)
    np.testing.assert_array_equal(design, again)
    for column in design.T:
        assert sorted(np.floor(column * 10).astype(int).tolist()) == list(range(10))


def test_gev_quantiles_through_array_overload_match_closed_form() -> None:
    # Arrange
    xi, alpha, kappa = 10.0, 2.0, -0.1
    gev = distributions.GeneralizedExtremeValue(xi, alpha, kappa)
    p = np.array([0.01, 0.5, 0.9, 0.99, 0.999])
    # Act
    q = gev.inverse_cdf(p)
    # Assert
    np.testing.assert_allclose(q, xi + alpha / kappa * (1.0 - (-np.log(p)) ** kappa), rtol=1e-12)
    for value, expected in zip(p, q, strict=True):
        assert gev.inverse_cdf(float(value)) == expected
    np.testing.assert_allclose(gev.cdf(q), p, rtol=1e-12)
    np.testing.assert_allclose(gev.ccdf(q), 1.0 - p, rtol=1e-10)
    sample = gev.generate_random_values(100, seed=7)
    np.testing.assert_array_equal(sample, gev.generate_random_values(100, seed=7))


def test_empirical_survival_orientation_against_hand_values() -> None:
    # Arrange: exceedance probabilities S(x) strictly descending; linear interpolation in both spaces
    dist = distributions.Empirical.survival([100.0, 200.0, 300.0], [0.5, 0.1, 0.01])
    dist.probability_transform = "none"
    # Act / Assert
    assert dist.probability_order == "Descending"
    assert dist.cdf(200.0) == 1.0 - 0.1
    assert dist.cdf(150.0) == pytest.approx(1.0 - 0.3, abs=1e-15)  # S(150) = (0.5 + 0.1) / 2
    assert dist.cdf(250.0) == pytest.approx(1.0 - 0.055, abs=1e-15)
    assert dist.inverse_cdf(0.7) == pytest.approx(150.0, abs=1e-10)  # reads S = 0.3 on the exceedance axis
    assert dist.inverse_cdf(0.9) == pytest.approx(200.0, abs=1e-10)
    np.testing.assert_allclose(dist.inverse_cdf(np.array([0.7, 0.945])), [150.0, 250.0], atol=1e-10)


def test_empirical_survival_clone_keeps_orientation() -> None:
    # Arrange
    dist = distributions.Empirical.survival([100.0, 200.0, 300.0], [0.5, 0.1, 0.01])
    clone = dist.clone()
    # Act / Assert
    assert clone.probability_order == "Descending"
    assert clone.cdf(250.0) == dist.cdf(250.0)


def test_empirical_from_sample_uses_plotting_positions() -> None:
    # Act
    dist = distributions.Empirical.from_sample([3.0, 1.0, 2.0, 4.0], "weibull")
    # Assert
    np.testing.assert_array_equal(dist.x_values, [1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(dist.probability_values, [0.2, 0.4, 0.6, 0.8])
    np.testing.assert_allclose(statistics.plotting_positions(4, "weibull"), [0.2, 0.4, 0.6, 0.8])
    np.testing.assert_allclose(statistics.plotting_positions(4, 0.5), [0.125, 0.375, 0.625, 0.875])


def test_uncertainty_analysis_results_on_a_tiny_ensemble() -> None:
    # Arrange: five normals with means 10..14; their medians are the means
    members = [distributions.Normal(10.0 + i, 1.0) for i in range(5)]
    uar = bootstrap.UncertaintyAnalysisResults()
    # Act
    ci = uar.process_confidence_intervals(members, [0.5], alpha=0.5)
    mean_curve = uar.process_mean_curve(members, [0.5])
    mode = uar.process_mode_curve(distributions.Normal(12.0, 1.0), [0.5])
    full = bootstrap.UncertaintyAnalysisResults.from_ensemble(distributions.Normal(12.0, 1.0), members, [0.1, 0.5])
    # Assert: Type-7 percentiles of {10, ..., 14} at 0.25 and 0.75
    np.testing.assert_allclose(ci, [[11.0, 13.0]], atol=1e-12)
    assert mean_curve[0] == pytest.approx(12.0, abs=1e-2)  # the symmetric mixture's median
    assert mode[0] == pytest.approx(12.0, abs=1e-12)
    assert full.confidence_intervals.shape == (2, 2)
    assert full.mean_curve.shape == (2,)


def test_bootstrap_analysis_smoke() -> None:
    # Arrange
    parent = distributions.Normal(100.0, 15.0)
    analysis = bootstrap.BootstrapAnalysis(parent, "method_of_moments", sample_size=50, replications=200, seed=12345)
    p = [0.1, 0.5, 0.9]
    # Act
    replicates = analysis.distributions()
    quantiles = analysis.quantiles(p)
    ci = analysis.percentile_quantile_ci(p, alpha=0.1)
    results = analysis.estimate(p, alpha=0.1)
    # Assert
    assert len(replicates) == 200
    assert isinstance(replicates[0], distributions.Normal)
    assert quantiles.shape == (200, 3)
    assert ci.shape == (3, 2) and np.all(ci[:, 0] < ci[:, 1])
    assert results.confidence_intervals.shape == (3, 2)
    assert analysis.failed_replications == 0


def test_statistics_mean_and_percentiles() -> None:
    # Arrange
    values = [4.0, 1.0, 3.0, 2.0, 5.0]
    # Act / Assert
    assert statistics.mean(values) == 3.0
    assert statistics.mean([1.0, 3.0], weights=[3.0, 1.0]) == 1.5
    assert statistics.percentile(values, 0.25) == 2.0
    assert statistics.median(values) == 3.0
    pair = [2.443080064681565, 10.0, 0.4453871940548014, -1.0]
    assert statistics.median(pair) == np.median(pair) == 1.4442336293681832
    assert statistics.median(sorted(pair), data_is_sorted=True) == np.median(pair)
    np.testing.assert_allclose(statistics.percentile(values, [0.1, 0.9]), [1.4, 4.6])
    np.testing.assert_allclose(statistics.percentile(values, [0.1, 0.9], weights=np.ones(5)), [1.4, 4.6])
    assert statistics.variance(values) == pytest.approx(2.5)
    running = statistics.RunningStatistics(values)
    running.push(6.0)
    assert running.count == 6
    assert running.mean == pytest.approx(3.5)
    assert running.variance == pytest.approx(np.var([1, 2, 3, 4, 5, 6], ddof=1))
    merged = statistics.RunningStatistics.combine(statistics.RunningStatistics([1.0, 2.0]),
                                                  statistics.RunningStatistics([3.0]))
    assert merged.count == 3 and merged.mean == pytest.approx(2.0)


def test_statistics_autocorrelation_kinds() -> None:
    # Arrange
    values = np.array([1.0, 3.0, 2.0, 5.0, 4.0, 6.0])
    d = values - values.mean()
    covariance = np.array([np.dot(d[k:], d[:values.size - k]) / values.size for k in range(4)])
    # Act
    lags, acvf = statistics.autocorrelation(values, 3, "covariance")
    _, acf = statistics.autocorrelation(values, 3)
    partial_lags, pacf = statistics.autocorrelation(values, 3, statistics.AutocorrelationType.PARTIAL)
    # Assert
    np.testing.assert_array_equal(lags, [0, 1, 2, 3])
    np.testing.assert_allclose(acvf, covariance, rtol=1e-14)
    np.testing.assert_allclose(acf, covariance / covariance[0], rtol=1e-14)
    np.testing.assert_array_equal(partial_lags, [1, 2, 3])
    assert pacf[0] == pytest.approx(acf[1], rel=1e-14)
    with pytest.raises(ValueError):
        statistics.autocorrelation(np.zeros(5), 2)


def test_paired_data_and_linear_interpolation() -> None:
    # Arrange
    opd = data.OrderedPairedData([0.0, 1.0, 2.0], [0.0, 10.0, 40.0])
    linear = data.Linear([1.0, 10.0, 100.0], [0.0, 1.0, 2.0], x_transform="logarithmic")
    # Act / Assert
    assert opd.is_valid and len(opd) == 3
    np.testing.assert_allclose(opd.get_y_from_x([0.5, 1.5, 5.0]), [5.0, 25.0, 40.0])
    assert opd.get_x_from_y(25.0) == pytest.approx(1.5)
    assert opd.get_y_from_x(3.0, extrapolation="both") == pytest.approx(70.0)
    np.testing.assert_allclose(linear([math.sqrt(10.0), 1000.0]), [0.5, 2.0], atol=1e-12)


def test_beta_poisson_mixture_and_elementary_functions() -> None:
    # Act / Assert
    assert distributions.Beta(2.0, 1.0).cdf(0.5) == pytest.approx(0.25, abs=1e-12)
    assert distributions.Beta(2.0, 3.0).inverse_cdf(distributions.Beta(2.0, 3.0).cdf(0.3)) == pytest.approx(0.3)
    assert distributions.Poisson(1.0).cdf(2.0) == pytest.approx(math.exp(-1.0) * 2.5, abs=1e-12)
    mixture = distributions.Mixture([0.5, 0.5], [distributions.Normal(0, 1), distributions.Normal(0, 1)])
    assert mixture.cdf(1.0) == pytest.approx(distributions.Normal(0, 1).cdf(1.0), abs=1e-12)
    assert mathematics.expm1(1e-10) == pytest.approx(1e-10 + 5e-21, rel=1e-15)
    np.testing.assert_allclose(mathematics.log1p(np.array([1e-12, 1.0])), [np.log1p(1e-12), np.log(2.0)], rtol=1e-15)


def test_other_families_construct_and_evaluate() -> None:
    # Arrange
    families = [
        distributions.GeneralizedPareto(0.0, 1.0, 0.1),
        distributions.Gumbel(10.0, 2.0),
        distributions.LogNormal(1.0, 0.2),
        distributions.LnNormal(10.0, 2.0),
        distributions.PearsonTypeIII(10.0, 2.0, 0.5),
        distributions.LogPearsonTypeIII(1.0, 0.2, 0.3),
        distributions.KappaFour(10.0, 2.0, 0.1, -0.2),
    ]
    p = np.array([0.1, 0.5, 0.9])
    for dist in families:
        # Act
        q = dist.inverse_cdf(p)
        # Assert
        assert dist.parameters_valid, dist.display_name
        np.testing.assert_allclose(dist.cdf(q), p, atol=1e-6, err_msg=dist.display_name)
    gumbel = distributions.Gumbel(10.0, 2.0)
    assert gumbel.inverse_cdf(0.5) == pytest.approx(10.0 - 2.0 * math.log(math.log(2.0)), rel=1e-12)
    # LogNormal is base 10: the median is 10**mean_of_log
    assert distributions.LogNormal(1.0, 0.2).inverse_cdf(0.5) == pytest.approx(10.0, rel=1e-12)


def test_estimate_fits_in_place() -> None:
    # Arrange
    sample = distributions.Normal(50.0, 5.0).generate_random_values(500, seed=3)
    fitted = distributions.Normal()
    # Act
    fitted.estimate(sample, "method_of_moments")
    # Assert
    assert fitted.mean == pytest.approx(float(np.mean(sample)), rel=1e-12)


def test_brent_root_and_bracket() -> None:
    # Act / Assert
    assert roots.brent_solve(lambda x: x * x - 2.0, 0.0, 2.0, tolerance=1e-12) == pytest.approx(math.sqrt(2.0))
    found, lo, hi, f1, f2 = roots.brent_bracket(lambda x: x - 5.0, 0.0, 1.0)
    assert found and lo <= 5.0 <= hi and f1 * f2 <= 0.0


def test_multivariate_normal() -> None:
    # Arrange
    mvn = multivariate.MultivariateNormal([1.0, 2.0], [[1.0, 0.5], [0.5, 2.0]])
    # Act / Assert
    assert mvn.dimension == 2
    np.testing.assert_allclose(mvn.covariance, [[1.0, 0.5], [0.5, 2.0]])
    np.testing.assert_allclose(mvn.inverse_cdf([0.5, 0.5]), [1.0, 2.0], atol=1e-12)
    sample = mvn.generate_random_values(2000, seed=11)
    assert sample.shape == (2000, 2)
    np.testing.assert_allclose(sample.mean(axis=0), [1.0, 2.0], atol=0.15)
    assert mvn.latin_hypercube_random_values(10, 5).shape == (10, 2)
    assert mvn.pdf([1.0, 2.0]) > 0.0


def test_copulas() -> None:
    # Act / Assert
    assert copulas.IndependenceCopula().cdf(0.3, 0.4) == pytest.approx(0.12, abs=1e-15)
    assert copulas.NormalCopula(0.0).cdf(0.3, 0.4) == pytest.approx(0.12, abs=1e-6)
    clayton = copulas.ClaytonCopula(2.0)
    np.testing.assert_allclose(clayton.cdf([0.3, 0.5], 0.4), [(0.3 ** -2 + 0.4 ** -2 - 1) ** -0.5,
                                                               (0.5 ** -2 + 0.4 ** -2 - 1) ** -0.5], rtol=1e-12)
    for copula in (copulas.GumbelCopula(1.5), copulas.FrankCopula(3.0), copulas.JoeCopula(1.5),
                   copulas.AMHCopula(0.3), copulas.StudentTCopula(0.5, 4.0)):
        sample = copula.generate_random_values(50, seed=5)
        assert sample.shape == (50, 2) and np.all((sample > 0) & (sample < 1))


def test_mersenne_twister_streams_agree() -> None:
    # Act
    bulk = sampling.MersenneTwister(12345).next_doubles(5)
    rng = sampling.MersenneTwister(12345)
    single = [rng.next_double() for _ in range(5)]
    # Assert
    np.testing.assert_array_equal(bulk, single)
    assert np.all((bulk >= 0.0) & (bulk < 1.0))


def test_scalar_members_are_bound_and_agree_with_the_array_overloads() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    p = np.array([0.01, 0.5, 0.99])
    # Act
    vector = gev.inverse_cdf(p)
    scalar = [gev.inverse_cdf(float(v)) for v in p]
    # Assert
    np.testing.assert_array_equal(vector, scalar)
    assert gev.log_cdf(5.0) == pytest.approx(math.log(gev.cdf(5.0)), rel=1e-14)
    assert gev.log_ccdf(5.0) == pytest.approx(math.log(gev.ccdf(5.0)), rel=1e-14)
    assert gev.log_pdf(5.0) == pytest.approx(math.log(gev.pdf(5.0)), rel=1e-14)


def test_array_members_bind_the_list_overload_once_and_keep_the_group_values() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    p = np.linspace(0.001, 0.999, 257)
    from pynumerics.marshal import to_net_doubles, to_numpy_doubles

    normal = distributions._nd().Normal
    # Act
    bound = gev._vector("InverseCDF")
    grouped = to_numpy_doubles(gev.net.InverseCDF(to_net_doubles(p)))
    z_grouped = to_numpy_doubles(normal.StandardZ(to_net_doubles(p)))
    linear = data.Linear([0.0, 1.0, 2.0], [0.0, 10.0, 40.0])
    x = np.array([0.25, 1.5, 3.0])
    # Assert
    assert gev._vector("InverseCDF") is bound  # bound once per instance
    assert distributions._static_list("StandardZ", normal) is distributions._static_list("StandardZ", normal)
    np.testing.assert_array_equal(gev.inverse_cdf(p).view(np.uint64), grouped.view(np.uint64))
    np.testing.assert_array_equal(distributions.standard_z(p).view(np.uint64), z_grouped.view(np.uint64))
    np.testing.assert_array_equal(gev.cdf(grouped), to_numpy_doubles(gev.net.CDF(to_net_doubles(grouped))))
    np.testing.assert_array_equal(gev.pdf(grouped), to_numpy_doubles(gev.net.PDF(to_net_doubles(grouped))))
    np.testing.assert_array_equal(distributions.standard_pdf(z_grouped),
                                  to_numpy_doubles(normal.StandardPDF(to_net_doubles(z_grouped))))
    np.testing.assert_array_equal(linear.interpolate(x), [2.5, 25.0, 40.0])


def test_parameter_names_validation_and_covariance() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    # Act / Assert
    assert gev.type_name == "GeneralizedExtremeValue"
    assert gev.parameter_names == ("Xi", "Alpha", "Kappa")
    assert gev.number_of_parameters == 3
    assert gev.validate_parameters([4.0, 1.0, -0.1]) is None
    assert "positive" in gev.validate_parameters([4.0, -1.0, 0.0])
    unit = gev.parameter_covariance(1)
    np.testing.assert_allclose(gev.parameter_covariance(10), unit / 10.0, rtol=1e-12)
    assert unit.shape == (3, 3) and np.allclose(unit, unit.T)
    with pytest.raises(NotImplementedError, match="exposes no parameter covariance"):
        distributions.Uniform(0.0, 1.0).parameter_covariance(1)


def test_parameter_set_members_equal_one_distribution_per_row() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    rows = np.array([[4.0, 1.0, -0.1], [5.0, 2.0, 0.1], [4.0, -1.0, 0.0]])
    # Act
    valid = gev.parameter_set_validity(rows)
    quantiles = gev.parameter_set_quantiles(rows[valid], 0.99)
    members = gev.parameter_set_distributions(rows, valid)
    # Assert
    np.testing.assert_array_equal(valid, [True, True, False])
    expected = [distributions.GeneralizedExtremeValue(*row).inverse_cdf(0.99) for row in rows[:2]]
    np.testing.assert_array_equal(quantiles, expected)
    assert members[2] is None
    assert [float(members[i].InverseCDF(0.99)) for i in range(2)] == expected


def test_parameter_set_quantile_table_equals_the_per_row_loop_bitwise() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    rng = np.random.default_rng(20261021)
    rows = np.column_stack([rng.normal(4.0, 0.3, 200), rng.lognormal(0.0, 0.2, 200), rng.normal(-0.1, 0.2, 200)])
    rows[7, 1] = -1.0  # an invalid scale
    probabilities = np.array([1e-8, 0.05, 0.5, 0.95, 0.99, 1.0 - 1e-9])
    # Act
    quantiles, valid = gev.parameter_set_quantile_table(rows, probabilities)
    # Assert
    np.testing.assert_array_equal(valid, gev.parameter_set_validity(rows))
    assert not valid[7] and np.isnan(quantiles[7]).all()
    assert quantiles.shape == (rows.shape[0], probabilities.size)
    for k, p in enumerate(probabilities):
        expected = gev.parameter_set_quantiles(rows[valid], p)
        np.testing.assert_array_equal(quantiles[valid, k].view(np.int64), expected.view(np.int64))
    with pytest.raises(ValueError, match="parameters per row"):
        gev.parameter_set_quantile_table(rows[:, :2], probabilities)
    with pytest.raises(ValueError, match="one-dimensional"):
        gev.parameter_set_quantile_table(rows, probabilities.reshape(2, 3))


def test_create_distribution_uniform_and_mixture_members() -> None:
    # Act
    gumbel = distributions.create_distribution("gumbel", [10.0, 2.0])
    uniform = distributions.Uniform(2.0, 4.0)
    mixture = distributions.Mixture([0.25, 0.75], [distributions.Normal(), distributions.Normal(1.0, 2.0)])
    # Assert
    assert isinstance(gumbel, distributions.Gumbel)
    np.testing.assert_array_equal(gumbel.parameters, [10.0, 2.0])
    assert uniform.inverse_cdf(0.5) == 3.0
    np.testing.assert_array_equal(mixture.weights, [0.25, 0.75])
    assert [type(c).__name__ for c in mixture.components] == ["Normal", "Normal"]
    assert gumbel.to_xml().startswith("<Distribution")


def test_linear_extrapolation_and_linear_moments() -> None:
    # Arrange
    linear = data.Linear([1.0, 2.0, 3.0], [1.0, 4.0, 9.0])
    # Act / Assert
    np.testing.assert_array_equal(linear.interpolate([0.0, 1.5, 4.0]), [1.0, 2.5, 9.0])
    np.testing.assert_array_equal(linear.extrapolate([0.0, 1.5, 4.0]), [-2.0, 2.5, 14.0])
    moments = statistics.linear_moments([1.0, 2.0, 3.0, 4.0])
    np.testing.assert_allclose(moments[:2], [2.5, 0.8333333333333334], rtol=1e-14)
    assert np.all(np.isnan(statistics.linear_moments([1.0, 2.0])))


def test_stratification_bin_records_and_value_bins_of_several_ranges() -> None:
    # Act
    bins = sampling.stratify_probabilities(0.01, 0.99, 4, "normal").bins
    many = sampling.stratify_x_values_many([(0.0, 10.0, 2), (10.0, 100.0, 3)])
    # Assert
    assert len(bins) == 4 and sum(b.weight for b in bins) == pytest.approx(1.0, abs=1e-12)
    assert bins[0].midpoint == pytest.approx((bins[0].lower_bound + bins[0].upper_bound) / 2.0)
    assert bins[1].contains(bins[1].lower_bound) and not bins[1].contains(bins[1].upper_bound)
    np.testing.assert_allclose(many.lower, [0.0, 5.0, 10.0, 40.0, 70.0])


def test_brent_search_special_functions_matrices_and_links() -> None:
    # Act
    result = optimization.BrentSearch(lambda x: (x - 1.3) ** 2, 0.0, 3.0, relative_tolerance=1e-12).minimize()
    maximum = optimization.BrentSearch(lambda x: -(x - 0.5) ** 2, 0.0, 2.0).maximize()
    links = functions.LinkController(["identity", "log", "identity"])
    # Assert
    assert result.parameters[0] == pytest.approx(1.3, abs=1e-9)
    assert result.status is optimization.OptimizationStatus.SUCCESS
    assert maximum.parameters[0] == pytest.approx(0.5, abs=1e-6)
    p = special.beta_incomplete(2.0, 3.0, 0.4)
    assert special.beta_incomplete_inverse(2.0, 3.0, p) == pytest.approx(0.4, abs=1e-12)
    np.testing.assert_array_equal(np.asarray(linalg.Matrix(np.eye(2) * 2.0).multiply([[1.0, 2.0], [3.0, 4.0]])),
                                  [[2.0, 4.0], [6.0, 8.0]])
    np.testing.assert_allclose(linalg.make_symmetric_positive_definite([[1.0, 0.5], [0.4, 1.0]]),
                               [[1.0, 0.45], [0.45, 1.0]])
    factor = linalg.CholeskyDecomposition([[4.0, 2.0], [2.0, 3.0]]).l
    np.testing.assert_allclose(factor @ factor.T, [[4.0, 2.0], [2.0, 3.0]], rtol=1e-14)
    np.testing.assert_allclose(links.link([4.0, 2.0, 0.1]), [4.0, math.log(2.0), 0.1])
    np.testing.assert_allclose(np.diag(links.link_jacobian([4.0, 2.0, 0.1])), [1.0, 0.5, 1.0])
    np.testing.assert_allclose(links.inverse_link([4.0, 0.0, 0.1]), [4.0, 1.0, 0.1])


def test_uncertainty_results_accept_a_net_member_array() -> None:
    # Arrange
    gev = distributions.GeneralizedExtremeValue(4.0, 1.0, -0.1)
    rows = np.array([[4.0 + 0.01 * i, 1.0, -0.1] for i in range(50)])
    members = gev.parameter_set_distributions(rows)
    wrapped = [distributions.GeneralizedExtremeValue(*row) for row in rows]
    # Act
    from_array = bootstrap.UncertaintyAnalysisResults.from_ensemble(gev, members, [0.5, 0.99])
    from_list = bootstrap.UncertaintyAnalysisResults.from_ensemble(gev, wrapped, [0.5, 0.99])
    # Assert
    np.testing.assert_array_equal(from_array.mean_curve, from_list.mean_curve)
    np.testing.assert_array_equal(from_array.confidence_intervals, from_list.confidence_intervals)


def test_net_errors_translate_the_numerics_exception_families() -> None:
    # Arrange
    from pynumerics.errors import net_errors

    # Act / Assert
    with pytest.raises(ValueError, match="between 0 and 1"):
        with net_errors():
            distributions.Normal().inverse_cdf(1.5)
    with pytest.raises(KeyError):
        with net_errors():
            raise KeyError("python errors pass through")


def test_bootstrap_members_and_parameters_stay_aligned() -> None:
    # Arrange
    analysis = bootstrap.BootstrapAnalysis(distributions.Gumbel(100.0, 10.0), "MethodOfMoments", 50, 200, 7)
    # Act
    members = analysis.members()
    rows = analysis.parameters()
    # Assert
    assert len(members) == 200 and rows.shape == (200, 2)
    assert [float(v) for v in members[3].GetParameters] == rows[3].tolist()


def test_ln_normal_log_parameter_construction() -> None:
    # Act
    ln = distributions.LnNormal.from_log_parameters(2.0, 0.1)
    mu, sigma = distributions.LnNormal.direct_method_of_moments(*ln.parameters)
    # Assert
    assert (float(ln.net.Mu), float(ln.net.Sigma)) == (2.0, 0.1)
    assert mu == pytest.approx(2.0, abs=1e-12) and sigma == pytest.approx(0.1, abs=1e-12)
    assert distributions.GeneralizedExtremeValue(4.0, 1.0, 0.0).has_parameter_covariance


HYDROLOGICS_MONTHLY_ALPHA = [3.0, 5.0, 9.0, 14.0, 11.0, 6.0, 2.5, 1.5, 4.0, 7.0, 3.5, 2.0]
"""Twelve monthly concentrations shaped like a flood-seasonality count table (sum 68.5)."""


def test_dirichlet_draws_are_bit_identical_to_the_net_call() -> None:
    # Arrange
    from pynumerics._util import namespace

    nd = namespace("Numerics.Distributions")
    System = namespace("System")
    alpha = System.Array.CreateInstance(System.Double, len(HYDROLOGICS_MONTHLY_ALPHA))
    for index, value in enumerate(HYDROLOGICS_MONTHLY_ALPHA):
        alpha[index] = value
    for count, seed in ((1, 1), (257, 12345), (1000, 2147483647)):
        reference = nd.Dirichlet(alpha).GenerateRandomValues(count, seed)
        expected = np.array([[float(reference[i, j]) for j in range(12)] for i in range(count)])
        # Act
        sample = multivariate.Dirichlet(HYDROLOGICS_MONTHLY_ALPHA).generate_random_values(count, seed)
        # Assert
        assert sample.shape == (count, 12) and sample.dtype == np.float64 and sample.flags.c_contiguous
        np.testing.assert_array_equal(sample.view(np.int64), expected.view(np.int64))
    wrapped = multivariate.Dirichlet.from_net(nd.Dirichlet(alpha))
    np.testing.assert_array_equal(wrapped.generate_random_values(50, 9).view(np.int64),
                                  multivariate.Dirichlet(HYDROLOGICS_MONTHLY_ALPHA)
                                  .generate_random_values(50, 9).view(np.int64))


def test_dirichlet_seeds_reproduce_and_rows_lie_on_the_simplex() -> None:
    # Arrange
    dirichlet = multivariate.Dirichlet([0.5, 2.0, 3.0])
    # Act
    first = dirichlet.generate_random_values(500, seed=42)
    again = dirichlet.generate_random_values(500, seed=42)
    other = dirichlet.generate_random_values(500, seed=43)
    # Assert
    np.testing.assert_array_equal(first, again)
    assert not np.array_equal(first, other)
    assert np.all(first >= 0.0) and np.all(first <= 1.0)
    np.testing.assert_allclose(first.sum(axis=1), 1.0, atol=1e-14)
    assert dirichlet.generate_random_values(0, seed=42).shape == (0, 3)


def test_dirichlet_moments_and_density_against_closed_forms() -> None:
    # Arrange
    alpha = np.array([2.0, 3.0, 5.0])
    a0 = alpha.sum()
    dirichlet = multivariate.Dirichlet(alpha)
    x = np.array([0.2, 0.3, 0.5])
    log_beta = sum(math.lgamma(a) for a in alpha) - math.lgamma(a0)
    # Act / Assert
    assert dirichlet.dimension == 3 and dirichlet.alpha_sum == 10.0 and dirichlet.parameters_valid
    np.testing.assert_array_equal(dirichlet.alpha, alpha)
    np.testing.assert_allclose(dirichlet.mean, alpha / a0, rtol=1e-15)
    np.testing.assert_allclose(dirichlet.variance, alpha * (a0 - alpha) / (a0 ** 2 * (a0 + 1.0)), rtol=1e-14)
    np.testing.assert_allclose(dirichlet.mode, (alpha - 1.0) / (a0 - 3.0), rtol=1e-14)
    expected_cov = (np.diag(alpha / a0) - np.outer(alpha, alpha) / a0 ** 2) / (a0 + 1.0)
    np.testing.assert_allclose(dirichlet.covariance_matrix(), expected_cov, rtol=1e-13, atol=1e-16)
    assert dirichlet.covariance(0, 1) == pytest.approx(expected_cov[0, 1], rel=1e-13)
    assert multivariate.Dirichlet.log_multivariate_beta(alpha) == pytest.approx(log_beta, rel=1e-12)
    expected_log_pdf = float(np.sum((alpha - 1.0) * np.log(x))) - log_beta
    assert dirichlet.log_pdf(x) == pytest.approx(expected_log_pdf, rel=1e-12)
    assert dirichlet.pdf(x) == pytest.approx(math.exp(expected_log_pdf), rel=1e-12)
    assert dirichlet.log_pdf([0.5, 0.5, 0.5]) == -math.inf
    symmetric = multivariate.Dirichlet.symmetric(4, 2.5)
    np.testing.assert_array_equal(symmetric.alpha, [2.5, 2.5, 2.5, 2.5])


def test_dirichlet_sample_mean_recovers_alpha_over_alpha_sum() -> None:
    # Arrange
    alpha = np.array(HYDROLOGICS_MONTHLY_ALPHA)
    dirichlet = multivariate.Dirichlet(alpha)
    count = 20000
    # Act
    sample = dirichlet.generate_random_values(count, seed=20261021)
    # Assert
    standard_error = np.sqrt(dirichlet.variance / count)
    assert np.all(np.abs(sample.mean(axis=0) - alpha / alpha.sum()) < 5.0 * standard_error)
    np.testing.assert_allclose(sample.var(axis=0, ddof=1), dirichlet.variance, rtol=0.05)


def test_dirichlet_argument_errors_map_to_python_exceptions() -> None:
    # Arrange
    from pynumerics._util import namespace
    from pynumerics.errors import net_errors

    System = namespace("System")

    # Act / Assert
    with pytest.raises(ValueError, match="must be positive") as raised:
        multivariate.Dirichlet([1.0, 0.0])
    assert isinstance(raised.value.__cause__, System.ArgumentOutOfRangeException)
    for alpha in ([1.0], [1.0, 0.0], [2.0, -1.0], [1.0, math.nan], [1.0, math.inf]):
        with pytest.raises(ValueError, match="alpha"):
            with net_errors():
                multivariate.Dirichlet(alpha)
    with pytest.raises(ValueError, match="dimension"):
        with net_errors():
            multivariate.Dirichlet.symmetric(1, 2.0)
    dirichlet = multivariate.Dirichlet([0.5, 2.0, 3.0])
    with pytest.raises(RuntimeError, match="mode"):
        with net_errors():
            _ = dirichlet.mode
    with pytest.raises(ValueError):
        with net_errors():
            dirichlet.covariance(0, 3)
    with pytest.raises(ArithmeticError):
        with net_errors():
            dirichlet.generate_random_values(-1, seed=1)

def test_poisson_inverse_cdf_matches_scipy_ppf_and_instance_bitwise() -> None:
    # Arrange: scipy 1.18.1 poisson.ppf, exact in integers
    stats = pytest.importorskip("scipy.stats")
    rates = np.array([0.1, 1.0, 3.7, 25.0, 150.0, 1000.0, 2e5])[:, None]
    p = np.array([1e-10, 0.01, 0.3, 0.5, 0.9, 0.999999])[None, :]
    # Act
    q = distributions.poisson_inverse_cdf(rates, p)
    # Assert
    assert q.shape == (7, 6)
    np.testing.assert_array_equal(q, stats.poisson.ppf(p, rates))
    for rate, row in zip(rates.ravel(), q, strict=True):
        np.testing.assert_array_equal(row, distributions.Poisson(rate).inverse_cdf(p.ravel()))
    assert distributions.poisson_inverse_cdf(2.0, 1.0) == math.inf
    assert isinstance(distributions.poisson_inverse_cdf(2.0, 0.5), float)


def test_poisson_log_pmf_and_cdf_match_scipy_and_instance_bitwise() -> None:
    # Arrange
    stats = pytest.importorskip("scipy.stats")
    rates = np.array([0.5, 3.7, 3.7, 25.0, 1000.0])
    k = np.array([0.0, 2.0, 7.0, 30.0, 980.0])
    # Act
    log_pmf = distributions.poisson_log_pmf(rates, k)
    cdf = distributions.poisson_cdf(rates, k)
    # Assert
    np.testing.assert_allclose(log_pmf, stats.poisson.logpmf(k, rates), rtol=1e-13)
    np.testing.assert_allclose(cdf, stats.poisson.cdf(k, rates), rtol=1e-13)
    for i in range(k.size):
        d = distributions.Poisson(rates[i])
        assert log_pmf[i] == d.log_pdf(k[i])
        assert cdf[i] == d.cdf(k[i])
    assert distributions.poisson_log_pmf(2.0, 1.5) == -math.inf


def test_binomial_log_pmf_and_cdf_match_scipy_and_instance_bitwise() -> None:
    # Arrange: the log mass forms ln C(n, k) from log factorials, so its error grows with ln n!
    stats = pytest.importorskip("scipy.stats")
    p = np.array([0.3, 0.5, 0.01, 0.9, 0.123, 0.6, 0.6])
    n = np.array([10, 1, 500, 40, 1000, 25, 25])
    k = np.array([3.0, 1.0, 2.0, 35.0, 100.0, 0.0, 25.0])
    # Act
    log_pmf = distributions.binomial_log_pmf(p, n, k)
    cdf = distributions.binomial_cdf(p, n, k)
    # Assert
    np.testing.assert_allclose(log_pmf, stats.binom.logpmf(k, n, p), rtol=0, atol=2e-12)
    np.testing.assert_allclose(cdf, stats.binom.cdf(k, n, p), rtol=1e-13)
    for i in range(k.size):
        d = distributions.Binomial(p[i], int(n[i]))
        assert log_pmf[i] == d.log_pdf(k[i])
        assert cdf[i] == d.cdf(k[i])
    broadcast = distributions.binomial_log_pmf(0.3, 10, np.arange(11.0))
    np.testing.assert_allclose(np.exp(broadcast).sum(), 1.0, rtol=1e-14)


def test_discrete_statics_reject_invalid_arguments() -> None:
    with pytest.raises(ValueError):
        distributions.poisson_inverse_cdf(-1.0, 0.5)
    with pytest.raises(ValueError):
        distributions.poisson_inverse_cdf(1.0, 1.5)
    with pytest.raises(ValueError):
        distributions.poisson_log_pmf([1.0, 2.0], [1.0, 2.0, 3.0])
    with pytest.raises(ValueError):
        distributions.binomial_log_pmf(1.5, 10, 2.0)
    with pytest.raises(ValueError):
        distributions.binomial_cdf(0.5, 0, 0.0)
    with pytest.raises(ValueError):
        distributions.binomial_cdf(0.5, 2.5, 0.0)
    with pytest.raises(ValueError):
        distributions.binomial_cdf(0.5, True, 0.0)
    with pytest.raises(ValueError):
        distributions.binomial_log_pmf(0.5, "10", 2.0)
    with pytest.raises(ValueError):
        distributions.binomial_log_pmf(0.5, 2**40, 2.0)
    with pytest.raises(ValueError):
        distributions.binomial_cdf(0.5, math.nan, 0.0)


@pytest.mark.parametrize("call", [
    "d.poisson_inverse_cdf([0.8, 2.5], [0.5, 0.9])",
    "d.poisson_log_pmf([0.8, 2.5], [1, 3])",
    "d.poisson_cdf([0.8, 2.5], [1, 3])",
    "d.binomial_log_pmf([0.3, 0.5], [10, 4], [3, 2])",
    "d.binomial_cdf([0.3, 0.5], [10, 4], [3, 2])",
])
def test_discrete_statics_work_as_first_call_in_fresh_process(call: str) -> None:
    # Arrange: the package source importable in a new interpreter that has not started the runtime
    import subprocess  # noqa: PLC0415
    import sys  # noqa: PLC0415
    from pathlib import Path  # noqa: PLC0415

    import pynumerics  # noqa: PLC0415

    source = str(Path(pynumerics.__file__).resolve().parents[1])
    code = f"import sys; sys.path.insert(0, {source!r}); from pynumerics import distributions as d; print({call})"
    # Act
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    # Assert
    assert result.returncode == 0, result.stderr
