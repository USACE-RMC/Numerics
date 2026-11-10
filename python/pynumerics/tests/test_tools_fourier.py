"""``Numerics.Tools`` and ``Numerics.Mathematics.Fourier`` wrappers against NumPy, SciPy and the standard library.

The FFT conventions are those of Numerical Recipes: the forward transform equals
``n * numpy.fft.ifft`` and the inverse equals ``numpy.fft.fft``. Arrays cross once each way.
"""

from __future__ import annotations

import math
import sys
import time
import zlib

import numpy as np
import pytest

from pynumerics import mathematics as m
from pynumerics._util import namespace
from pynumerics.marshal import to_net_doubles

sp = pytest.importorskip("scipy.special")
acf = pytest.importorskip("statsmodels.tsa.stattools").acf


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def test_log_sum_exp_matches_scipy() -> None:
    # Arrange
    rng = np.random.default_rng(1)
    values = rng.normal(scale=300.0, size=1000)
    # Act / Assert
    np.testing.assert_allclose(m.log_sum_exp(values), sp.logsumexp(values), rtol=1e-15)
    np.testing.assert_allclose(m.log_sum_exp(values[:50], values[50:100]),
                               np.logaddexp(values[:50], values[50:100]), rtol=1e-15)
    assert m.log_sum_exp([-math.inf, -math.inf]) == -math.inf
    assert m.log_sum_exp([]) == -math.inf
    assert m.log_sum_exp(1.0, 2.0) == float(namespace("Numerics").Tools.LogSumExp(1.0, 2.0))


def _scaled_ulps(actual: np.ndarray, expected: np.ndarray, largest: np.ndarray) -> np.ndarray:
    scale = np.maximum(np.abs(expected), np.abs(largest))
    return np.abs(actual - expected) / np.spacing(scale)


def test_log_sum_exp_special_values_follow_scipy() -> None:
    # SciPy 1.18: +inf gives +inf, NaN propagates, all -inf and an empty input give -inf.
    for values in ([1.0, math.inf], [-math.inf, math.inf], [-math.inf] * 3, [], [1.0, math.inf, math.nan]):
        expected = float(sp.logsumexp(values)) if values else -math.inf
        actual = m.log_sum_exp(values)
        assert (math.isnan(actual) and math.isnan(expected)) or actual == expected
    assert m.log_sum_exp(math.inf, 1.0) == math.inf
    assert math.isnan(m.log_sum_exp(math.inf, math.nan))


def test_log_sum_exp_near_zero_matches_scipy() -> None:
    # Arrange: results near zero, where the separated largest term and log1p keep relative accuracy.
    # Exp and Log come from the platform math library, so Linux results differ from Windows by about an
    # ulp; the bound is 2 scaled ulps on Windows and 3 elsewhere (measured worst on Ubuntu 24.04: 3).
    limit = 2.0 if sys.platform == "win32" else 3.0
    rng = np.random.default_rng(20261108)
    worst = 0.0
    for _ in range(300):
        values = rng.normal(scale=10 ** rng.uniform(-2, 2), size=int(rng.integers(2, 200)))
        values -= sp.logsumexp(values) - rng.uniform(-1e-3, 1e-3)
        # Act
        actual = m.log_sum_exp(values)
        # Assert
        worst = max(worst, float(_scaled_ulps(np.array(actual), np.array(sp.logsumexp(values)), values.max())))
    assert worst <= limit
    assert m.log_sum_exp([0.0, -30.0]) == sp.logsumexp([0.0, -30.0]) == 9.357622968839737e-14


def test_log_sum_exp_axis_matches_scipy_and_list() -> None:
    # Arrange
    rng = np.random.default_rng(7)
    matrix = rng.normal(scale=5.0, size=(40, 7))
    matrix[3, 2] = matrix[10, 2] = matrix[:, 2].max() + 1.0
    matrix[5, 4] = math.inf
    matrix[6, 5] = math.nan
    matrix[:, 6] = -math.inf
    for axis in (0, 1, -1, -2):
        # Act
        actual = m.log_sum_exp(matrix, axis=axis)
        # Assert: bitwise the list form of each column or row, and within 2 ulps of SciPy.
        lines = matrix.T if axis % 2 == 0 else matrix
        assert actual.shape == (lines.shape[0],)
        np.testing.assert_array_equal(_bits(actual), _bits([m.log_sum_exp(line) for line in lines]))
        expected = sp.logsumexp(matrix, axis=axis)
        finite = np.isfinite(expected)
        np.testing.assert_array_equal(np.isnan(actual), np.isnan(expected))
        np.testing.assert_array_equal(actual[~finite & ~np.isnan(expected)], expected[~finite & ~np.isnan(expected)])
        largest = np.nanmax(np.where(np.isfinite(lines), lines, -np.inf), axis=1)
        assert _scaled_ulps(actual[finite], expected[finite], largest[finite]).max() <= 2.0
    np.testing.assert_array_equal(m.log_sum_exp(np.empty((0, 3)), axis=0), [-math.inf] * 3)
    with pytest.raises(ValueError, match="2-D"):
        m.log_sum_exp([1.0, 2.0], axis=0)
    with pytest.raises(ValueError, match="axis must be"):
        m.log_sum_exp(matrix, axis=2)
    with pytest.raises(TypeError, match="axis must be an integer"):
        m.log_sum_exp(matrix, axis=True)
    with pytest.raises(TypeError, match="axis must be an integer"):
        m.log_sum_exp(matrix, axis=1.0)
    np.testing.assert_array_equal(_bits(m.log_sum_exp(matrix, axis=np.int64(1))), _bits(m.log_sum_exp(matrix, axis=1)))
    with pytest.raises(ValueError, match="not both"):
        m.log_sum_exp(matrix, 1.0, axis=0)


def test_log_sum_exp_reduces_in_one_crossing() -> None:
    # Arrange
    values = np.random.default_rng(2).normal(size=200_000)
    m.log_sum_exp(values[:10])
    # Act
    start = time.perf_counter()
    m.log_sum_exp(values)
    elapsed = time.perf_counter() - start
    # Assert: a crossing per element would take a tenth of a second or more
    assert elapsed < 0.05


def test_scalar_tools_match_python() -> None:
    # Act / Assert
    np.testing.assert_array_equal(m.copy_sign([3.0, -3.0, 2.0], [-1.0, 1.0, -0.0]), [-3.0, 3.0, 2.0])
    np.testing.assert_array_equal(m.sign([2.0, -2.0], [-5.0, 0.0]), [-2.0, 2.0])
    np.testing.assert_array_equal(m.clamp([-1.0, 0.5, 2.0], 0.0, 1.0), [0.0, 0.5, 1.0])
    assert m.sqr(3.0) == 9.0
    np.testing.assert_array_equal(m.power(2.0, [10, -2, 0]), [1024.0, 0.25, 1.0])
    assert m.power(0.0, -1) == math.inf
    np.testing.assert_array_equal(m.log10([1000.0, 0.0]), [3.0, -16.0])
    assert m.log(math.e) == 1.0 and m.log(0.0) == math.log(1e-16) and m.log(-0.0) == math.log(1e-16)
    assert m.is_finite(1.0) and not m.is_finite(math.nan) and not m.is_finite(math.inf)
    assert m.is_power_of_two(64) and not m.is_power_of_two(0) and not m.is_power_of_two(12)
    assert [m.next_power_of_two(n) for n in (0, 1, 5, 8, 1000)] == [1, 1, 8, 8, 1024]


def test_hypot_matches_numpy_without_overflow() -> None:
    # Arrange
    rng = np.random.default_rng(3)
    a = rng.standard_normal(500) * 10.0 ** rng.integers(-300, 300, 500)
    b = rng.standard_normal(500) * 10.0 ** rng.integers(-300, 300, 500)
    # Act
    h = m.hypot(a, b)
    # Assert
    # Tools.Hypot is within 2 ulps of the exact value and numpy.hypot within 1, so they differ by at most 3.
    np.testing.assert_array_max_ulp(h, np.hypot(a, b), maxulp=3)
    assert m.hypot(3.0, 4.0) == 5.0
    assert m.hypot(1e308, 1e308) == pytest.approx(math.sqrt(2) * 1e308)
    assert m.hypot(1e-320, 0.0) == 1e-320
    np.testing.assert_array_equal(m.hypot([3.0, -5.0], 4.0), [5.0, math.hypot(5.0, 4.0)])
    special = [math.inf, -math.inf, math.nan, 0.0, 2.0]
    for x in special:
        for y in special:
            np.testing.assert_array_equal(m.hypot(x, y), np.hypot(x, y))


def test_distances_match_numpy() -> None:
    # Arrange
    x, y = np.array([1.0, 2.0, 3.0]), np.array([4.0, -1.0, 0.5])
    lower, upper = np.zeros(3), np.array([10.0, 0.0, 5.0])
    w = np.array([1.0, 2.0, 0.5])
    # Act / Assert
    assert m.point_distance(0.0, 0.0, 3.0, 4.0) == 5.0
    np.testing.assert_allclose(m.distance(x, y), np.linalg.norm(x - y), rtol=1e-15)
    keep = upper > lower
    expected = np.linalg.norm(((x - lower) / np.where(keep, upper - lower, 1) - (y - lower)
                               / np.where(keep, upper - lower, 1))[keep])
    np.testing.assert_allclose(m.normalized_distance(x, y, lower, upper), expected, rtol=1e-15)
    np.testing.assert_allclose(m.weighted_distance(x, y, w), math.sqrt(np.sum(w * (x - y) ** 2)), rtol=1e-15)


def test_scaling_round_trips_match_numpy() -> None:
    # Arrange
    values = np.random.default_rng(3).normal(10.0, 3.0, size=50)
    # Act
    normalized = m.normalize(values)
    standardized = m.standardize(values)
    # Assert
    np.testing.assert_allclose(normalized, (values - values.min()) / np.ptp(values), rtol=1e-15, atol=1e-16)
    np.testing.assert_allclose(m.denormalize(normalized, values.min(), values.max()), values, rtol=1e-15)
    np.testing.assert_allclose(standardized, (values - values.mean()) / values.std(ddof=1), rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(m.destandardize(standardized, values.mean(), values.std(ddof=1)), values, rtol=1e-14)
    np.testing.assert_array_equal(m.normalize([2.0, 2.0]), [0.0, 0.0])


def test_reductions_match_numpy() -> None:
    # Arrange
    values = np.array([3.0, -1.0, 4.0, 1.5, -9.0])
    indicators = np.array([1, 0, 1, 1, 0])
    chosen, rest = values[indicators == 1], values[indicators == 0]
    # Act / Assert
    assert m.sum_values(values) == values.sum()
    assert m.sum_values(np.array([1, 2, 3])) == 6.0
    assert m.sum_values(values, indicators) == chosen.sum()
    assert m.sum_values(values, indicators, use_complement=True) == rest.sum()
    assert m.mean_values(values) == pytest.approx(values.mean(), rel=1e-15)
    assert m.mean_values(values, indicators) == pytest.approx(chosen.mean(), rel=1e-15)
    assert m.product(values) == np.prod(values)
    assert m.product(values, indicators, True) == np.prod(rest)
    assert m.min_value(values) == -9.0 and m.min_value(values, indicators) == chosen.min()
    assert m.max_value(values) == 4.0 and m.max_value(values, indicators, True) == rest.max()
    assert m.min_max(values) == (-9.0, 4.0)
    assert all(math.isnan(v) for v in m.min_max([1.0, math.nan]))
    assert m.arg_min(values) == 4 and m.arg_max([math.nan, 2.0, 7.0]) == 2 and m.arg_max([math.nan]) == -1
    assert m.sum_product(values, values) == pytest.approx(values @ values, rel=1e-15)


MIXED = [  # random.Random(20261107): 24 terms over twenty decades, both signs
    -4.415159360150449e-08, 7.317112953068489, 9.868669377084806e-07, -0.0006463016681697819, -7.502006205048129e-05,
    -0.01673876386175016, -0.0006856973900070548, -825842102090.3063, 75526526787.25735, 2.3726148037552075e-08,
    -9.761977681447285e-08, -51658895725.59849, 1.3059917476939964e-07, -65237681874.17293, 994405.7906876896,
    -2.0600436727650427e-07, -2.2111435313817664e-09, -349.5172509335271, -9.36662463129794e-09, -3.733797169087079e-09,
    0.024886648499566268, 45263.7648646135, 0.8450434057865094, 4.4580288798240146e-09,
]


@pytest.mark.parametrize("values", [
    [1e100, 1.0, -1e100],             # plain summation gives 0
    [1.0, 1e100, 1.0, -1e100],        # a term above the running sum: Kahan summation gives 0
    [1e16, 1.0, -1e16, 1e-16, 3.0, -2.0**-60],
    [0.1] * 10,
    MIXED,
    [5.0],
])
def test_compensated_sum_equals_math_fsum(values) -> None:
    # Act / Assert: bitwise the correctly rounded sum on these ill-conditioned inputs
    assert m.compensated_sum(values) == math.fsum(values)
    assert m.compensated_sum(np.asarray(values)) == math.fsum(values)


def test_compensated_sum_stays_within_the_neumaier_bound() -> None:
    # Arrange: cancelling terms over twelve decades
    rng = np.random.default_rng(11)
    values = rng.standard_normal(20000) * 10.0 ** rng.integers(-6, 6, 20000)
    values = np.concatenate([values, -values[:19000] * (1 + 1e-12)])
    exact = math.fsum(values)
    unit = 2.0 ** -53
    bound = 2 * unit * abs(exact) + 2 * values.size * unit * unit * np.abs(values).sum()
    plain = 0.0
    for value in values.tolist():
        plain += value
    # Act
    total = m.compensated_sum(values)
    # Assert
    assert abs(total - exact) <= bound
    assert m.sum_values(values) == plain
    assert abs(plain - exact) > abs(total - exact)


def test_compensated_sum_non_finite_and_empty_follow_ieee() -> None:
    # Act / Assert: the IEEE sum, as math.fsum for NaN and one infinity; math.fsum raises where IEEE gives NaN or inf
    assert math.isnan(m.compensated_sum([1.0, math.nan, 2.0])) and math.isnan(math.fsum([1.0, math.nan, 2.0]))
    assert m.compensated_sum([1.0, math.inf, -1e308]) == math.fsum([1.0, math.inf, -1e308]) == math.inf
    assert m.compensated_sum([-math.inf, 1.0]) == math.fsum([-math.inf, 1.0]) == -math.inf
    assert math.isnan(m.compensated_sum([math.inf, 1.0, -math.inf]))
    with pytest.raises(ValueError):
        math.fsum([math.inf, 1.0, -math.inf])
    with pytest.raises(OverflowError):
        math.fsum([1e308, 1e308, -1e308])
    assert m.compensated_sum([1e308, 1e308, -1e308]) == math.inf
    assert math.isnan(m.compensated_sum([])) and math.fsum([]) == 0.0


def test_compensated_sum_indicators_and_gather() -> None:
    # Arrange
    values = np.array([1e100, 1.0, -1e100, 7.0, 3.0])
    indicators = np.array([1, 1, 1, 0, 0])
    indices = np.array([4, 0, 1, 2, 1, 1])
    # Act / Assert
    assert m.compensated_sum(values, indicators) == math.fsum(values[indicators == 1]) == 1.0
    assert m.compensated_sum(values, indicators, use_complement=True) == math.fsum(values[indicators == 0])
    assert m.compensated_sum(values, [0, 0, 0, 0, 0]) == 0.0
    assert math.isnan(m.compensated_sum(values, [1, 1]))
    assert m.compensated_sum(values, indices=indices) == math.fsum(values[indices])
    assert m.compensated_sum(values, indices=indices, start=1, count=3) == math.fsum(values[indices[1:4]]) == 1.0
    assert m.compensated_sum(values, indices=indices, start=4) == 2.0
    assert math.isnan(m.compensated_sum(values, indices=indices, start=6))
    with pytest.raises(ValueError):
        m.compensated_sum(values, indices=indices, start=2, count=5)
    with pytest.raises(ValueError):
        m.compensated_sum(values, indices=[0, 5])
    with pytest.raises(ValueError):
        m.compensated_sum(values, indicators, indices=indices)


def test_compensated_cumulative_sum_equals_prefix_fsum() -> None:
    # Arrange
    values = np.array([1e100, 1.0, -1e100, 1e-100, 0.1, 0.2, *MIXED])
    # Act
    running = m.compensated_cumulative_sum(values)
    # Assert: each entry is bitwise compensated_sum of the prefix; these prefixes are also math.fsum
    assert running.dtype == np.float64 and running.shape == values.shape
    assert [m.compensated_sum(values[:i + 1]) for i in range(values.size)] == running.tolist()
    assert running[:6].tolist() == [math.fsum(values[:i + 1]) for i in range(6)]
    assert m.compensated_cumulative_sum([]).size == 0
    tail = m.compensated_cumulative_sum([1.0, math.inf, 2.0, -math.inf])
    assert tail[1] == math.inf and math.isnan(tail[3])


def test_compensated_sum_beats_math_fsum_on_a_large_array() -> None:
    # Arrange
    values = np.random.default_rng(12).standard_normal(1_000_000)
    m.compensated_sum(values[:10])
    # Act
    start = time.perf_counter()
    total = m.compensated_sum(values)
    numerics_seconds = time.perf_counter() - start
    start = time.perf_counter()
    exact = math.fsum(values)
    fsum_seconds = time.perf_counter() - start
    # Assert
    assert total == pytest.approx(exact, rel=1e-15, abs=1e-12)
    assert numerics_seconds < fsum_seconds, (numerics_seconds, fsum_seconds)


def test_sequence_and_compression() -> None:
    # Arrange
    data = b"the quick brown fox " * 50
    # Act
    packed = m.compress(data)
    # Assert
    np.testing.assert_array_equal(m.sequence(0, 10, 3), [0, 3, 6, 9])
    assert m.sequence(0, 10, 3).dtype == np.int32
    np.testing.assert_array_equal(m.sequence(5, 1, -2), [5, 3, 1])
    np.testing.assert_array_equal(m.sequence(0.0, 1.0, 0.25), [0.0, 0.25, 0.5, 0.75, 1.0])
    assert m.sequence(0, 3, 0).size == 0
    assert len(packed) < len(data)
    assert zlib.decompress(packed, -15) == data
    assert m.decompress(zlib.compress(data)[2:-4]) == data
    assert m.decompress(packed) == data


def test_constants_are_the_numerics_values() -> None:
    # Act / Assert
    c = m.constants
    assert c.double_machine_epsilon == pytest.approx(2.0**-53, rel=1e-14)
    # The single-precision unit roundoff 2**-24
    assert c.single_machine_epsilon == 2.0**-24
    assert c.euler == pytest.approx(np.euler_gamma, rel=1e-16)
    assert c.golden_ratio == pytest.approx((1 + math.sqrt(5)) / 2, rel=1e-16)
    assert c.log_sqrt_2pi == pytest.approx(0.5 * math.log(2 * math.pi), rel=1e-16)
    assert c.sqrt_2pi == pytest.approx(math.sqrt(2 * math.pi), rel=1e-16)
    assert c.sqrt2 == math.sqrt(2) and c.log2 == math.log(2)


@pytest.mark.parametrize("n", [2, 8, 1024])
def test_fft_matches_numpy(n: int) -> None:
    # Arrange
    rng = np.random.default_rng(n)
    x = rng.normal(size=n) + 1j * rng.normal(size=n)
    original = x.copy()
    # Act
    forward = m.fft(x)
    inverse = m.fft(forward, inverse=True)
    # Assert
    scale = np.abs(x).sum()
    np.testing.assert_allclose(forward, n * np.fft.ifft(x), rtol=0, atol=1e-14 * scale * max(1, math.log2(n)))
    np.testing.assert_allclose(m.fft(x, True), np.fft.fft(x), rtol=0, atol=1e-14 * scale * max(1, math.log2(n)))
    np.testing.assert_allclose(inverse / n, x, rtol=0, atol=1e-14 * scale)
    np.testing.assert_array_equal(x, original)


def test_fft_returns_the_bits_of_the_net_member() -> None:
    # Arrange
    x = np.random.default_rng(5).normal(size=64) + 1j * np.random.default_rng(6).normal(size=64)
    data = to_net_doubles(x.view(np.float64))
    # Act
    namespace("Numerics.Mathematics").Fourier.FFT(data, False)
    # Assert
    np.testing.assert_array_equal(_bits(m.fft(x).view(np.float64)), _bits([data[i] for i in range(128)]))


@pytest.mark.parametrize("n", [3, 6, 100, 1000])
def test_bluestein_fft_matches_numpy_for_any_length(n: int) -> None:
    # Arrange
    rng = np.random.default_rng(n)
    x = rng.normal(size=n) + 1j * rng.normal(size=n)
    # Act / Assert
    scale = np.abs(x).sum() * math.log2(n)
    np.testing.assert_allclose(m.bluestein_fft(x), n * np.fft.ifft(x), rtol=0, atol=1e-14 * scale)
    np.testing.assert_allclose(m.bluestein_fft(x, inverse=True), np.fft.fft(x), rtol=0, atol=1e-14 * scale)


def test_fft_refuses_a_length_that_is_not_a_power_of_two() -> None:
    # Act / Assert
    with pytest.raises(ValueError, match="power of 2"):
        m.fft(np.ones(6, dtype=complex))
    with pytest.raises(ValueError, match="1-D"):
        m.fft(np.ones((2, 2), dtype=complex))


def test_real_fft_packs_the_positive_frequencies() -> None:
    # Arrange
    n = 16
    r = np.random.default_rng(7).normal(size=n)
    full = n * np.fft.ifft(r)
    # Act
    packed = m.real_fft(r)
    back = m.real_fft(packed, inverse=True) * 2.0 / n
    # Assert
    atol = 1e-14 * np.abs(r).sum() * 4
    np.testing.assert_allclose(packed[0], full[0].real, atol=atol)
    np.testing.assert_allclose(packed[1], full[n // 2].real, atol=atol)
    np.testing.assert_allclose(packed[2::2] + 1j * packed[3::2], full[1:n // 2], atol=atol)
    np.testing.assert_allclose(back, r, atol=atol)


def test_convolution_and_correlation_match_numpy() -> None:
    # Arrange
    rng = np.random.default_rng(8)
    a, b = rng.normal(size=37), rng.normal(size=11)
    x, y = rng.normal(size=16), rng.normal(size=16)
    # Act
    conv = m.fft_convolve(a, b)
    corr = m.fft_correlation(x, y)
    # Assert
    direct = np.convolve(a, b)
    np.testing.assert_allclose(conv, direct, rtol=0, atol=1e-14 * np.abs(direct).max() * 8)
    expected = np.array([np.sum(x[(np.arange(16) + k) % 16] * y) for k in range(16)])
    np.testing.assert_allclose(corr, expected, rtol=0, atol=1e-14 * np.abs(expected).max() * 4)


def test_autocorrelation_matches_statsmodels() -> None:
    # Arrange
    series = np.cumsum(np.random.default_rng(9).normal(size=300))
    # Act
    table = m.fft_autocorrelation(series, 20)
    # Assert
    np.testing.assert_array_equal(table[:, 0], np.arange(21))
    np.testing.assert_allclose(table[:, 1], acf(series, nlags=20, fft=True), rtol=0, atol=1e-14)


def test_large_fft_crosses_once() -> None:
    # Arrange
    x = np.random.default_rng(10).normal(size=2**18).astype(complex)
    m.fft(x[:8])
    # Act
    start = time.perf_counter()
    m.fft(x)
    elapsed = time.perf_counter() - start
    # Assert: half a million doubles copied element by element would take seconds
    assert elapsed < 0.5


def _kahan(values: np.ndarray) -> tuple[float, list[float]]:
    t = c = 0.0
    running = []
    for v in values.tolist():
        a = v - c
        f = t + a
        c = (f - t) - a
        t = f
        running.append(t)
    return t, running


def test_kahan_compensated_sums_equal_the_classic_loop_bitwise() -> None:
    # Arrange
    rng = np.random.default_rng(21)
    # Act / Assert
    for n in (1, 2, 7, 50, 5000):
        for scale in (1.0, 1e-12, 1e8):
            w = rng.random(n) * scale
            w[rng.random(n) < 0.1] *= 1e6
            total, running = _kahan(w)
            assert m.compensated_sum(w, kahan=True) == total
            assert m.compensated_cumulative_sum(w, kahan=True).tolist() == running
            assert m.compensated_sum(w) == m.compensated_sum(w, kahan=False)
            picks = (rng.random(n) < 0.5).astype(np.int32)
            assert m.compensated_sum(w, picks, kahan=True) == _kahan(w[picks == 1])[0]
            idx = rng.integers(0, n, n)
            assert m.compensated_sum(w, indices=idx, kahan=True) == _kahan(w[idx])[0]
    assert math.isnan(m.compensated_sum([1.0, math.inf, 1.0], kahan=True))
