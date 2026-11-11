"""Sampling designs, the Mersenne Twister and the bootstrap of ``pynumerics.sampling``.

Oracles: NumPy's legacy ``RandomState`` (the reference MT19937 with ``init_genrand`` and
``init_by_array``), SciPy's unscrambled Sobol sequence, SciPy's normal distribution for the
stratification transforms, and ``scipy.stats.bootstrap`` plus the textbook standard error of the
mean for the bootstrap intervals.
"""

from __future__ import annotations

import math
import threading
import time

import numpy as np
import pytest
from scipy import stats
from scipy.stats import qmc

from pynumerics import callbacks, distributions, sampling

needs_roslyn = pytest.mark.skipif(not callbacks.roslyn_available(), reason="Roslyn assemblies not staged")


# Mersenne Twister


def test_mersenne_twister_matches_the_reference_mt19937_streams() -> None:
    # Arrange
    reference = np.random.RandomState(5489)
    by_array = np.random.RandomState([1, 2, 3])
    twister = sampling.MersenneTwister(1)
    twister.initialize(5489)
    seeded = sampling.MersenneTwister([1, 2, 3])
    # Act
    raw = [twister.gen_rand_int32() for _ in range(5)]
    array_raw = [seeded.gen_rand_int32() for _ in range(5)]
    # Assert: init_genrand(5489) starts 3499211612, the published MT19937 value
    assert raw[0] == 3499211612
    assert raw == reference.randint(0, 2 ** 32, size=5, dtype=np.uint32).tolist()
    assert array_raw == by_array.randint(0, 2 ** 32, size=5, dtype=np.uint32).tolist()


def test_mersenne_twister_real_outputs_follow_the_reference_definitions() -> None:
    # Arrange
    twister, raw = sampling.MersenneTwister(42), sampling.MersenneTwister(42)
    res53 = np.random.RandomState(0)
    twister53 = sampling.MersenneTwister(0)
    twister53.initialize(0)
    # Act / Assert: genrand_int31 = int32 >> 1, real1 = int32 / (2^32 - 1), real2 = int32 / 2^32,
    # real3 = (int32 + 0.5) / 2^32, res53 = random_sample
    assert twister.gen_rand_int31() == raw.gen_rand_int32() >> 1
    assert twister.gen_rand_real1() == raw.gen_rand_int32() * (1.0 / 4294967295.0)
    assert twister.gen_rand_real2() == raw.gen_rand_int32() * (1.0 / 4294967296.0)
    assert twister.gen_rand_real3() == (raw.gen_rand_int32() + 0.5) * (1.0 / 4294967296.0)
    assert [twister53.gen_rand_res53() for _ in range(4)] == res53.random_sample(4).tolist()


def test_mersenne_twister_without_a_seed_draws_from_the_clock() -> None:
    # Act
    draw = sampling.MersenneTwister().next_double()
    # Assert
    assert 0.0 <= draw < 1.0


# Sobol and Latin hypercube


def test_unscrambled_sobol_sequence_matches_scipy_after_the_origin() -> None:
    # Arrange
    sequence = sampling.SobolSequence(3)
    expected = qmc.Sobol(3, scramble=False).random_base2(6)[1:]
    # Act
    points = sequence.next_doubles(63)
    # Assert
    assert sequence.dimension == 3 and sequence.seed is None
    assert np.array_equal(points, expected)
    assert np.array_equal(sampling.SobolSequence(3).skip_to(10), points[9])
    assert np.array_equal(sampling.SobolSequence(3).skip_to(0), points[0])


def test_scrambled_sobol_sequence_is_reproducible_and_in_the_unit_cube() -> None:
    # Act
    first = sampling.SobolSequence(2, seed=11).next_doubles(128)
    again = sampling.SobolSequence(2, seed=11).next_doubles(128)
    other = sampling.SobolSequence(2, seed=12).next_doubles(128)
    # Assert
    assert sampling.SobolSequence(2, seed=11).seed == 11
    assert np.array_equal(first, again) and not np.array_equal(first, other)
    assert np.all((first >= 0.0) & (first < 1.0))


def test_latin_hypercube_takes_a_generator_and_stratifies_every_column() -> None:
    # Act
    design = sampling.latin_hypercube(20, 3, sampling.MersenneTwister(5))
    again = sampling.latin_hypercube(20, 3, sampling.MersenneTwister(5))
    median = sampling.latin_hypercube(20, 2, sampling.MersenneTwister(5), median=True)
    # Assert: one point per 1/20 stratum in each column
    assert np.array_equal(design, again)
    for column in design.T:
        assert sorted(np.floor(column * 20).astype(int)) == list(range(20))
    assert np.allclose(np.sort(median[:, 0]), (np.arange(20) + 0.5) / 20)


# Stratification


def test_value_bins_map_to_probability_bins_through_a_cdf() -> None:
    # Arrange
    bins = sampling.stratify_x_values(0.0, 10.0, 5)
    normal = distributions.Normal(5.0, 2.0)
    # Act
    exhaustive = sampling.x_to_probability(bins, normal.cdf)
    partial = sampling.x_to_probability(bins, normal.cdf, is_exhaustive=False)
    exceedance = sampling.x_to_exceedance_probability(bins, lambda x: 1.0 - normal.cdf(x), is_exhaustive=False)
    # Assert
    cdf = stats.norm(5.0, 2.0).cdf(np.linspace(0.0, 10.0, 6))
    assert partial.weight == pytest.approx(np.diff(cdf), rel=1e-12)
    assert exhaustive.weight.sum() == pytest.approx(1.0, abs=1e-12)
    assert exceedance.weight == pytest.approx(np.diff(cdf), rel=1e-9)


def test_probability_bins_map_to_values_through_the_quantile_function() -> None:
    # Arrange
    probabilities = sampling.stratify_probabilities(0.01, 0.99, 4, is_exhaustive=False)
    exceedances = sampling.stratify_exceedance_probabilities(0.01, 0.99, 4, is_exhaustive=False)
    normal = stats.norm(0.0, 1.0)
    # Act
    x = sampling.probability_to_x(probabilities, distributions.Normal(0.0, 1.0).inverse_cdf)
    xe = sampling.exceedance_probability_to_x(exceedances, lambda q: float(normal.isf(q)))
    # Assert
    assert x.lower == pytest.approx(normal.ppf(probabilities.lower), rel=1e-12)
    assert x.upper == pytest.approx(normal.ppf(probabilities.upper), rel=1e-12)
    assert x.weight == pytest.approx(x.upper - x.lower, rel=1e-12)
    assert xe.lower == pytest.approx(normal.isf(exceedances.upper), rel=1e-12)
    assert xe.upper == pytest.approx(normal.isf(exceedances.lower), rel=1e-12)


def test_many_range_stratifications_concatenate_the_ranges() -> None:
    # Act
    many = sampling.stratify_probabilities_many([(0.0, 0.5, 2), sampling.StratificationOptions(0.5, 1.0, 3, True)])
    exceedance = sampling.stratify_exceedance_probabilities_many([(0.0, 0.5, 2), (0.5, 1.0, 2)])
    values = sampling.stratify_x_values_many([sampling.StratificationOptions(0.0, 1.0, 2), (1.0, 2.0, 2)])
    # Assert
    assert len(many) == 5 and len(exceedance) == 4 and len(values) == 4
    assert values.upper.tolist() == [0.5, 1.0, 1.5, 2.0]


def test_multivariate_probabilities_give_one_bin_list_per_dimension() -> None:
    # Act
    first = sampling.stratify_multivariate_probabilities(0.01, 0.99, 3, dimension=2, seed=3)
    again = sampling.stratify_multivariate_probabilities(0.01, 0.99, 3, dimension=2, seed=3)
    correlated = sampling.stratify_multivariate_probabilities(0.01, 0.99, 3, dimension=2, seed=3,
                                                              correlation=[[1.0, 0.8], [0.8, 1.0]])
    # Assert
    assert len(first) == 2 and all(len(b) == 3 for b in first)
    assert [b.lower.tolist() for b in first] == [b.lower.tolist() for b in again]
    assert len(correlated) == 2


def test_stratification_options_validate_clone_and_round_trip_through_xml() -> None:
    # Arrange
    options = sampling.StratificationOptions(0.0, 1.0, 4, is_probability=True)
    invalid = sampling.StratificationOptions(1.0, 0.0, 0)
    # Act
    copy = options.clone()
    restored = sampling.StratificationOptions.from_xml(options.to_xml())
    listed = sampling.StratificationOptions.list_from_xml(sampling.StratificationOptions.list_to_xml([options, copy]))
    # Assert
    assert options.is_valid and options.get_errors() == []
    assert not invalid.is_valid and len(invalid.get_errors()) >= 1
    for item in (copy, restored, *listed):
        assert (item.lower_bound, item.upper_bound, item.number_of_bins, item.is_probability) == (0.0, 1.0, 4, True)
    assert "number_of_bins" not in repr(options) and "4" in repr(options)


def test_bins_round_trip_through_xml_and_order_by_position() -> None:
    # Arrange
    bins = sampling.stratify_probabilities(0.001, 0.999, 4, "normal")
    # Act
    restored = sampling.StratificationBins.from_xml(bins.to_xml())
    records = bins.bins
    # Assert
    assert restored.lower.tolist() == bins.lower.tolist() and restored.weight.tolist() == bins.weight.tolist()
    assert records[0].compare_to(records[1]) == -1 and records[1].compare_to(records[0]) == 1
    assert records[0].to_net().Midpoint == records[0].midpoint
    overlapping = sampling.StratificationBin(0.0, 1.0, 1.0, 0.5)
    with pytest.raises(ValueError, match="overlapping"):
        overlapping.compare_to(sampling.StratificationBin(0.5, 1.5, 1.0, 1.0))


# Bootstrap


X = stats.norm(10.0, 2.0).ppf(sampling.MersenneTwister(5).next_doubles(40))
"""Forty observations from N(10, 2)."""


def _resample(data: np.ndarray, parameters: np.ndarray, prng: object) -> np.ndarray:
    """Draw a nonparametric resample with the replicate's Numerics generator."""
    return data[[prng.Next(data.size) for _ in range(data.size)]]


def _mean_bootstrap(replicates: int = 2000) -> sampling.Bootstrap:
    return sampling.Bootstrap(X, [X.mean()], resample=_resample, fit=lambda d: [d.mean()], statistic=lambda p: p,
                              jackknife=lambda d, i: np.delete(d, i), sample_size=lambda d: d.size,
                              replicates=replicates)


def test_bootstrap_of_the_mean_matches_scipy_and_the_textbook_standard_error() -> None:
    # Arrange
    boot = _mean_bootstrap()
    scipy_ci = {m: stats.bootstrap((X,), np.mean, n_resamples=20000, confidence_level=0.9, method=m,
                                   random_state=np.random.default_rng(1)).confidence_interval
                for m in ("percentile", "BCa")}
    se = X.std(ddof=0) / math.sqrt(X.size)
    # Act
    boot.run()
    percentile = boot.get_confidence_intervals("percentile", 0.1).parameter_results[0]
    bca_results = boot.get_confidence_intervals(sampling.BootstrapCIMethod.BCA, 0.1)
    bca = bca_results.statistic_results[0]
    normal = boot.get_confidence_intervals("normal", 0.1).statistic_results[0]
    bias = boot.get_confidence_intervals("bias_corrected", 0.1)
    # Assert: Monte Carlo agreement at a quarter of a standard error (about five Monte Carlo errors)
    assert boot.failed_replicates == 0 and boot.failed_attempts == 0 and boot.failed_jackknife_replicates == 0
    assert percentile.population_estimate == X.mean()
    assert percentile.standard_error == pytest.approx(se, rel=0.05)
    assert percentile.lower_ci == pytest.approx(scipy_ci["percentile"].low, abs=0.25 * se)
    assert percentile.upper_ci == pytest.approx(scipy_ci["percentile"].high, abs=0.25 * se)
    assert bca.lower_ci == pytest.approx(scipy_ci["BCa"].low, abs=0.25 * se)
    assert bca.upper_ci == pytest.approx(scipy_ci["BCa"].high, abs=0.25 * se)
    assert normal.lower_ci < X.mean() < normal.upper_ci
    assert bias.method is sampling.BootstrapCIMethod.BIAS_CORRECTED and bias.alpha == 0.1
    # Parameter intervals are percentile intervals whatever the method
    assert bca_results.parameter_results[0].lower_ci == percentile.lower_ci
    assert bca.lower_ci != percentile.lower_ci
    assert percentile.valid_count == percentile.total_count == 2000
    assert boot.bootstrap_parameter_sets.shape == (2000, 1) and boot.bootstrap_statistics.shape == (2000, 1)
    assert percentile.mean == pytest.approx(boot.bootstrap_statistics.mean())


def test_bootstrap_seed_reproduces_the_replicates() -> None:
    # Act
    first = _mean_bootstrap(300).run().bootstrap_statistics
    again = _mean_bootstrap(300).run().bootstrap_statistics
    other = _mean_bootstrap(300)
    other.seed = 1
    # Assert
    assert np.array_equal(first, again)
    assert not np.array_equal(first, other.run().bootstrap_statistics)


def test_python_bootstrap_functions_run_serially_with_the_parallel_bits() -> None:
    # Arrange: a Python callable turns UseParallel off, so the replicates run in order on this thread
    threads: list[int] = []

    def fit(d: np.ndarray) -> list[float]:
        threads.append(threading.get_ident())
        return [d.mean()]

    serial = sampling.Bootstrap(X, [X.mean()], resample=_resample, fit=fit, statistic=lambda p: p, replicates=300)
    parallel = _mean_bootstrap(300)
    parallel.use_parallel = True
    # Act
    serial.run()
    parallel.run()
    # Assert
    assert not serial.use_parallel and parallel.use_parallel
    assert len(threads) == 300 and set(threads) == {threading.get_ident()}
    np.testing.assert_array_equal(serial.bootstrap_statistics.view(np.uint64),
                                  parallel.bootstrap_statistics.view(np.uint64))


def test_python_bca_jackknife_runs_serially_with_the_parallel_bits() -> None:
    # Arrange: with use_parallel False the leave-one-out jackknife also runs in order on this thread
    calls: list[tuple[int, int]] = []

    def jackknife(d: np.ndarray, i: int) -> np.ndarray:
        calls.append((i, threading.get_ident()))
        return np.delete(d, i)

    serial = sampling.Bootstrap(X, [X.mean()], resample=_resample, fit=lambda d: [d.mean()], statistic=lambda p: p,
                                jackknife=jackknife, sample_size=lambda d: d.size, replicates=300)
    parallel = _mean_bootstrap(300)
    parallel.use_parallel = True
    # Act
    serial.run()
    parallel.run()
    serial_bca = serial.get_confidence_intervals("bca", 0.1).statistic_results[0]
    parallel_bca = parallel.get_confidence_intervals("bca", 0.1).statistic_results[0]
    # Assert
    assert not serial.use_parallel
    assert [i for i, _ in calls] == list(range(X.size))
    assert {t for _, t in calls} == {threading.get_ident()}
    assert np.float64(serial_bca.lower_ci).view(np.uint64) == np.float64(parallel_bca.lower_ci).view(np.uint64)
    assert np.float64(serial_bca.upper_ci).view(np.uint64) == np.float64(parallel_bca.upper_ci).view(np.uint64)


def _compiled_mean_bootstrap() -> sampling.Bootstrap:
    """Return the 2,000-replicate mean bootstrap of _mean_bootstrap with C# resample, fit and statistic."""
    return sampling.Bootstrap(
        X, [X.mean()],
        resample="(d, p, r) => { var o = new double[d.Length]; for (int i = 0; i < d.Length; i++) "
                 "o[i] = d[r.Next(d.Length)]; return o; }",
        fit="d => new Numerics.Mathematics.Optimization.ParameterSet(new[] { d.Average() }, double.NaN)",
        statistic="p => p.Values", replicates=2000)


@needs_roslyn
def test_csharp_bootstrap_functions_draw_the_same_replicates() -> None:
    # Arrange
    compiled = _compiled_mean_bootstrap()
    python = _mean_bootstrap()
    # Act
    python.run()
    compiled.run()
    # Assert: the same resamples; the means differ only in summation order
    assert compiled.bootstrap_statistics == pytest.approx(python.bootstrap_statistics, rel=1e-14)
    assert compiled.use_parallel and not python.use_parallel


@pytest.mark.benchmark
@needs_roslyn
def test_csharp_bootstrap_functions_draw_the_replicates_faster() -> None:
    # Arrange
    compiled = _compiled_mean_bootstrap()
    python = _mean_bootstrap()
    # Act
    start = time.perf_counter()
    python.run()
    python_seconds = time.perf_counter() - start
    start = time.perf_counter()
    compiled.run()
    compiled_seconds = time.perf_counter() - start
    # Assert
    assert compiled_seconds < python_seconds, (
        f"Bootstrap, 2000 replicates: Python {python_seconds:.3f} s, compiled {compiled_seconds:.3f} s")


@needs_roslyn
def test_studentized_and_double_bootstrap_runs() -> None:
    # Arrange: C# functions keep the inner replicates fast
    boot = sampling.Bootstrap(
        X, [X.mean()],
        resample="(d, p, r) => { var o = new double[d.Length]; for (int i = 0; i < d.Length; i++) "
                 "o[i] = d[r.Next(d.Length)]; return o; }",
        fit="d => new Numerics.Mathematics.Optimization.ParameterSet(new[] { d.Average() }, double.NaN)",
        statistic="p => p.Values", replicates=300)
    boot.inner_replicates = 50
    boot.transform, boot.inverse_transform = "x => x", "x => x"
    # Act
    t = boot.run_with_studentized_bootstrap().get_confidence_intervals("bootstrap_t", 0.1).parameter_results[0]
    double = boot.run_double_bootstrap(50).get_confidence_intervals("percentile").parameter_results[0]
    # Assert
    assert boot.inner_replicates == 50 and boot.transform(4.0) == 4.0 and boot.inverse_transform(4.0) == 4.0
    assert t.lower_ci < X.mean() < t.upper_ci
    assert double.lower_ci < X.mean() < double.upper_ci


def test_pivotal_bootstrap_with_covariance_fits() -> None:
    # Arrange
    def fit_with_covariance(data: np.ndarray) -> tuple[list[float], list[list[float]]]:
        return [data.mean()], [[data.var(ddof=1) / data.size]]

    parent = sampling.BootstrapFit([X.mean()], [[X.var(ddof=1) / X.size]])
    boot = sampling.Bootstrap(X, parent, resample=_resample, fit_with_covariance=fit_with_covariance,
                              statistic=lambda p: p, replicates=300)
    boot.pivotal_replicate_filter = lambda fit: fit.parameter_count == 1
    boot.pivotal_parameter_validator = lambda values: bool(np.all(np.isfinite(values)))
    boot.pivotal_link_factory = lambda context: [None] * context.parameter_count
    boot.pivotal_invalid_draw_policy = "use_parent"
    boot.pivotal_z_limit, boot.add_pivotal_jitter, boot.pivotal_jitter_scale = 6.0, True, 0.02
    boot.regularize_pivotal_covariances = True
    # Act
    boot.run_pivotal_bootstrap()
    pivotal = boot.get_confidence_intervals("percentile").parameter_results[0]
    raw = boot.get_raw_pivotal_confidence_intervals(0.1).parameter_results[0]
    diagnostics = boot.pivotal_diagnostics
    # Assert
    assert boot.original_covariance == pytest.approx(parent.covariance)
    assert boot.pivotal_invalid_draw_policy is sampling.PivotalBootstrapInvalidDrawPolicy.USE_PARENT
    assert (boot.pivotal_z_limit, boot.add_pivotal_jitter, boot.pivotal_jitter_scale) == (6.0, True, 0.02)
    assert pivotal.lower_ci < X.mean() < pivotal.upper_ci and raw.lower_ci < X.mean() < raw.upper_ci
    assert diagnostics.requested_replicates == 300 and diagnostics.accepted_raw_replicates == 300
    assert diagnostics.failed_raw_replicates == diagnostics.rejected_raw_replicates == 0
    assert diagnostics.retained_pivotal_replicates + diagnostics.invalid_pivotal_replicates >= 300
    assert diagnostics.resampling_time >= 0.0 and diagnostics.transformation_time >= 0.0
    assert len(boot.raw_bootstrap_fits) == 300 and boot.raw_bootstrap_parameter_sets.shape == (300, 1)
    assert boot.raw_bootstrap_statistics.shape == (300, 1) and boot.pivotal_links == [None]
    # The transformation alone, from the raw fits of the run
    boot.transform_pivotal_bootstrap(boot.raw_bootstrap_fits)
    assert boot.get_confidence_intervals("percentile").parameter_results[0].lower_ci < X.mean()


def test_bootstrap_support_records() -> None:
    # Arrange
    parent = sampling.BootstrapFit([1.0, 2.0], np.eye(2))
    context = sampling.PivotalBootstrapContext(parent, [sampling.BootstrapFit([1.5, 2.5], np.eye(2)),
                                                        sampling.BootstrapFit([0.5, 1.5], np.eye(2))])
    # Act / Assert
    assert parent.parameters.tolist() == [1.0, 2.0] and parent.covariance.tolist() == np.eye(2).tolist()
    assert context.parameter_count == 2 and context.parent_fit.parameters.tolist() == [1.0, 2.0]
    assert context.get_raw_parameter_values(1).tolist() == [2.5, 1.5] and len(context.raw_bootstrap_fits) == 2
    assert sampling.BootstrapResults().statistic_results == [] and sampling.BootstrapResults().failed_replicates == 0
    assert sampling.BootstrapStatisticResult().valid_count == 0
    assert sampling.PivotalBootstrapDiagnostics().requested_replicates == 0
    assert "parameters" in repr(parent)
    with pytest.raises(ValueError, match="covariance"):
        sampling.BootstrapFit([1.0])


def test_bootstrap_failures_and_settings() -> None:
    # Arrange: every fit raises, so every attempt fails and is retried
    boot = sampling.Bootstrap(X, [X.mean()], resample=_resample, fit=lambda d: 1 / 0, statistic=lambda p: p,
                              replicates=10)
    boot.max_retries = 2
    progress: list[int] = []
    boot.progress = progress.append
    # Act
    boot.run()
    # Assert
    assert boot.max_retries == 2 and boot.failed_replicates == 10 and boot.failed_attempts == 20
    assert np.isnan(boot.bootstrap_parameter_sets).all()
    assert sorted(progress) == list(range(1, 11))
    assert boot.progress is not None and boot.cancellation_token is not None
    boot.cancellation_token = boot.cancellation_token
    boot.original_covariance = None
    assert boot.original_covariance is None and boot.pivotal_z_limit is None
    for name in ("resample_function", "fit_function", "statistic_function", "fit_with_covariance_function",
                 "jackknife_function", "sample_size_function", "pivotal_link_factory", "pivotal_replicate_filter",
                 "pivotal_parameter_validator"):
        getattr(boot, name)
    boot.jackknife_function = None
    assert boot.jackknife_function is None
    with pytest.raises(RuntimeError, match="ResampleFunction must be set"):
        sampling.Bootstrap(X, [X.mean()]).run()


def test_a_raising_progress_callback_fails_the_run_instead_of_the_process() -> None:
    # Arrange
    boot = _mean_bootstrap(5)

    def progress(count: int) -> None:
        raise RuntimeError(f"stop at {count}")

    boot.progress = progress
    # Act / Assert
    with pytest.raises(Exception, match="stop at"):
        boot.run()
