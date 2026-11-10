"""Hosking-Wallis regional frequency analysis against the North Cascades example of Hosking's LMOMENTS package.

Oracle: CASCADES.DAT and CASCADES.OUT of the LMOMENTS Fortran package, version 3.04
(https://lib.stat.cmu.edu/general/lmoments), PROGRAM XTEST with 500 simulations. Fixed values are checked to the
printed precision; H and Z within Monte Carlo tolerance, since the random number generators differ.
"""

from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pytest

from pynumerics import regional

N = [98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82]
MEAN = [19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932, 17.560, 69.518, 47.653,
        102.501, 52.413, 79.696, 44.643, 58.655, 39.024]
LCV = [0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179, 0.1308, 0.1119, 0.1018,
       0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046]
T3 = [0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492, 0.0940, -0.0429, 0.0435,
      0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128]
T4 = [0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900, 0.1273, 0.0927, 0.1446,
      0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583]
NAMES = [f"site{i + 1}" for i in range(19)]
PROBABILITIES = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.9, 0.95, 0.99, 0.999]


@pytest.fixture(scope="module")
def cascades() -> regional.RegionalData:
    return regional.RegionalData(NAMES, N, MEAN, LCV, T3, T4)


@pytest.fixture(scope="module")
def cascades_test(cascades: regional.RegionalData) -> regional.RegionalHomogeneityTest:
    return regional.RegionalHomogeneityTest(cascades, 500, 12345)


def test_discordancy_and_regional_average_match_cascades(cascades: regional.RegionalData) -> None:
    # Arrange
    expected = [0.60, 1.02, 0.38, 0.23, 0.93, 2.63, 2.12, 0.45, 0.11, 1.61, 2.08, 1.52, 0.31, 1.30, 1.58, 0.29,
                1.04, 0.43, 0.38]
    # Act
    d = cascades.discordancy()
    r = cascades.regional_average()
    # Assert
    np.testing.assert_allclose(d, expected, atol=0.005)
    np.testing.assert_allclose(r[:4], [1.0, 0.1103, 0.0279, 0.1366], atol=0.00005)
    assert np.isnan(r[4])
    assert regional.RegionalData.discordancy_critical_value(19) == 3.0
    assert regional.RegionalData.discordancy_critical_value(8, 0.05) == pytest.approx(2.2728)
    assert cascades.site_count == 19
    assert cascades.names[0] == "site1"
    np.testing.assert_array_equal(cascades.record_lengths, N)


def test_homogeneity_matches_cascades(cascades_test: regional.RegionalHomogeneityTest) -> None:
    # Assert: fixed values to the printed precision.
    np.testing.assert_allclose(cascades_test.kappa_parameters, [0.9542, 0.1533, 0.1236, -0.2955], atol=0.00005)
    np.testing.assert_allclose(cascades_test.observed_dispersion, [0.0104, 0.0339, 0.0405], atol=0.00005)
    np.testing.assert_allclose(cascades_test.distribution_tau4, [0.167, 0.111, 0.123, 0.123, 0.006], atol=0.0005)
    assert not cascades_test.is_generalized_logistic_parent
    # Assert: Monte Carlo values within about three standard errors.
    for observed, published in zip(cascades_test.heterogeneity, [0.62, -1.49, -2.37], strict=True):
        assert abs(observed - published) <= 0.3 + 0.1 * abs(published)
    for name, published in {"GLO": 3.46, "GEV": -2.94, "GNO": -1.51, "PE3": -1.60, "GPA": -14.75}.items():
        assert abs(cascades_test.z[name] - published) <= 0.3 + 0.1 * abs(published)
    assert cascades_test.get_z("GNO") == cascades_test.z["GNO"]
    assert set(cascades_test.acceptable_distributions()) == {"GeneralizedNormal", "PearsonTypeIII"}


def test_homogeneity_is_deterministic_across_thread_counts(cascades: regional.RegionalData) -> None:
    # Act
    serial = regional.RegionalHomogeneityTest(cascades, 100, 7, 1)
    parallel = regional.RegionalHomogeneityTest(cascades, 100, 7)
    # Assert
    np.testing.assert_array_equal(serial.heterogeneity, parallel.heterogeneity)
    assert serial.z == parallel.z


def test_growth_curve_matches_cascades(cascades: regional.RegionalData) -> None:
    # Act
    gno = regional.RegionalGrowthCurve("GNO", cascades)
    pe3 = regional.RegionalGrowthCurve("PearsonTypeIII", cascades)
    # Assert
    np.testing.assert_allclose(gno.parameters, [0.994, 0.195, -0.057], atol=0.0005)
    np.testing.assert_allclose(gno.growth_factors(PROBABILITIES),
                               [0.569, 0.616, 0.688, 0.753, 0.834, 0.994, 1.254, 1.331, 1.480, 1.654], atol=0.0005)
    np.testing.assert_allclose(pe3.growth_factors(PROBABILITIES),
                               [0.570, 0.616, 0.688, 0.753, 0.834, 0.994, 1.254, 1.331, 1.480, 1.653], atol=0.0005)
    assert gno.get_distribution().mean == pytest.approx(1.0, abs=1e-9)
    assert isinstance(gno.growth_factors(0.99), float)


def test_site_quantiles_scale_growth_factors(cascades: regional.RegionalData) -> None:
    # Act
    curve = regional.RegionalGrowthCurve.from_ratios("GEV", 0.1103, 0.0279)
    table = curve.site_quantiles(cascades, PROBABILITIES)
    one = curve.site_quantiles(MEAN[3], PROBABILITIES)
    # Assert
    assert table.shape == (19, 10)
    np.testing.assert_allclose(table[3], one, rtol=1e-15)
    np.testing.assert_allclose(one, MEAN[3] * curve.growth_factors(PROBABILITIES), rtol=1e-15)
    assert curve.data is None
    assert curve.distribution_type == "GeneralizedExtremeValue"


def test_from_samples_gives_exact_linear_moments() -> None:
    # Arrange: exact rational L-moments (Python fractions, mpmath 1.3.0 for printing).
    sample = [12.1, 15.3, 9.8, 22.4, 18.0, 13.7, 11.2, 25.9, 16.4, 14.8]
    # Act
    data = regional.RegionalData.from_samples([sample, [1.0, 2.0, np.nan, 3.5, 4.0]], names=["a", "b"])
    # Assert
    np.testing.assert_allclose(data.site_linear_moments(0),
                               [15.96, 2.9044444444444444, 0.23297628156082632, 0.18297081648267570,
                                0.045524100994644223], rtol=0, atol=1e-13)
    assert data.record_lengths.tolist() == [10, 4]
    assert np.isnan(data.tau5[1])


def test_invalid_inputs_raise(cascades: regional.RegionalData) -> None:
    # Act / Assert
    with pytest.raises(ValueError):
        regional.RegionalData(["a"], [10, 20], [1.0], [0.1], [0.0], [0.1])
    with pytest.raises(ValueError):
        regional.RegionalData.from_samples([[1.0, 2.0, 3.0]])
    with pytest.raises(ValueError):
        regional.RegionalHomogeneityTest(cascades, 1)
    with pytest.raises(ValueError):
        regional.RegionalGrowthCurve("Normal", cascades)
    with pytest.raises(ValueError):
        regional.RegionalData.discordancy_critical_value(10, 0.01)


COLD_STARTS = {
    "from_samples": "regional.RegionalData.from_samples([[1.0, 2.0, 3.5, 4.0, 6.0], [2.0, 3.0, 1.0, 5.0]])",
    "constructor": "regional.RegionalData(['a', 'b'], [10, 12], [1.0, 2.0], [0.1, 0.2], [0.0, 0.1], [0.1, 0.2])",
    "critical_value": "regional.RegionalData.discordancy_critical_value(10)",
    "from_ratios": "regional.RegionalGrowthCurve.from_ratios('GEV', 0.1, 0.05).growth_factors(0.99)",
    "homogeneity": ("regional.RegionalHomogeneityTest(regional.RegionalData(['a', 'b'], [10, 12], [1.0, 2.0],"
                    " [0.1, 0.2], [0.0, 0.1], [0.1, 0.2]), 10)"),
}
"""Each wrapper entry point, run as the first Numerics call of a fresh interpreter."""


@pytest.mark.parametrize("call", COLD_STARTS.values(), ids=COLD_STARTS.keys())
def test_entry_points_start_the_runtime(call: str) -> None:
    # Arrange
    code = f"from pynumerics import regional; print(type({call}).__name__)"
    # Act
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=os.environ.copy(),
                            timeout=120, check=False)
    # Assert
    assert result.returncode == 0, result.stderr


def test_results_report_nonpositive_means_and_validate_ratios() -> None:
    # Arrange: a large L-CV with records of four or five years gives simulated sites with a nonpositive mean.
    data = regional.RegionalData(list("abcde"), [4, 5, 4, 5, 4], [1.0] * 5, [0.70, 0.72, 0.68, 0.71, 0.69],
                                 [0.0, 0.01, -0.01, 0.02, -0.02], [0.12, 0.13, 0.11, 0.12, 0.12])
    # Act
    test = regional.RegionalHomogeneityTest(data, 200, 5)
    # Assert
    assert test.nonpositive_mean_simulations > 0
    with pytest.raises(ValueError):
        regional.RegionalData(["a"], [10], [1.0], [0.0], [0.0], [0.1])
    with pytest.raises(ValueError):
        regional.RegionalData(["a"], [10], [1.0], [0.1], [0.0], [0.1], [np.inf])
