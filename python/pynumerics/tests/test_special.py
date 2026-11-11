"""Special functions against scipy.special, and the array paths against the scalar Numerics members.

Each wrapper calls one Numerics member. The array path crosses into the ``IList<double>`` overload
once and must return the bits of the scalar member evaluated element by element. The scipy
tolerances are 1e-14 relative or better where the Numerics algorithm reaches double precision; the
oscillatory J and Y functions are compared absolutely next to their zeros above x = 8, where scipy
(Cephes) is accurate only in the absolute sense.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest

from pynumerics import special
from pynumerics._util import namespace

sp = pytest.importorskip("scipy.special")
integrate = pytest.importorskip("scipy.integrate")


def _bits(values: object) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.float64).view(np.uint64)


def _sf() -> object:
    return namespace("Numerics.Mathematics.SpecialFunctions")


X = np.linspace(0.05, 30.0, 200)


@pytest.mark.parametrize(("wrapper", "oracle", "x", "rtol"), [
    (special.gamma_function, sp.gamma, X, 1e-15),
    (special.gamma_function, sp.gamma, np.linspace(-9.7, -0.3, 50), 1e-15),
    (special.log_gamma, sp.gammaln, np.linspace(3.0, 300.0, 300), 1e-14),
    (special.digamma, sp.digamma, np.linspace(2.0, 30.0, 100), 1e-14),
    (special.trigamma, lambda v: sp.polygamma(1, v), X, 1e-15),
    (special.lanczos_gamma, sp.gamma, X, 1e-14),
    (special.erf, sp.erf, np.linspace(-5.0, 5.0, 201), 4e-15),
    (special.erfc, sp.erfc, np.linspace(-5.0, 5.0, 201), 4e-15),
    (special.inverse_erf, sp.erfinv, np.linspace(-0.999, 0.999, 199), 1e-15),
    (special.inverse_erfc, sp.erfcinv, np.linspace(0.001, 1.999, 199), 1e-15),
    (special.bessel_i0e, sp.i0e, np.linspace(0.1, 40.0, 150), 1e-15),
    (special.bessel_i1e, sp.i1e, np.linspace(0.1, 40.0, 150), 1e-15),
    (special.bessel_i0, sp.i0, np.linspace(0.1, 700.0, 300), 1e-14),
    (special.bessel_i1, sp.i1, np.linspace(0.1, 700.0, 300), 1e-14),
    (special.bessel_k0, sp.k0, np.linspace(0.01, 700.0, 300), 1e-14),
    (special.bessel_k1, sp.k1, np.linspace(0.01, 700.0, 300), 1e-14),
])
def test_double_precision_functions_match_scipy(wrapper, oracle, x, rtol) -> None:
    # Act
    values = wrapper(x)
    # Assert
    np.testing.assert_allclose(values, oracle(x), rtol=rtol, atol=0)


def test_log_gamma_and_digamma_near_their_zeros_match_scipy_absolutely() -> None:
    # Arrange
    x = np.linspace(0.05, 3.0, 120)
    # Act / Assert
    np.testing.assert_allclose(special.log_gamma(x), sp.gammaln(x), rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(special.digamma(x), sp.digamma(x), rtol=1e-14, atol=2e-15)


def test_log_gamma_matches_scipy_gammaln_across_the_double_range() -> None:
    # Arrange: scipy (Cephes) gammaln is accurate to about 4e-16 relative above 4, at the integers (with
    # ln Γ(1) = ln Γ(2) = 0 exactly) and at the half-integers, and only absolutely next to the zeros at 1
    # and 2, where it is compared absolutely.
    rng = np.random.default_rng(20261107)
    large = np.concatenate([10.0 ** rng.uniform(np.log10(4.0), 300.0, 2000), np.arange(1.0, 61.0),
                            np.arange(1.5, 200.0, 1.0)])
    near_zeros = np.concatenate([1.0 + s * 10.0 ** -np.arange(1.0, 16.0) for s in (1.0, -1.0)]
                                + [2.0 + s * 10.0 ** -np.arange(1.0, 15.0) for s in (1.0, -1.0)])
    tiny = 10.0 ** rng.uniform(-300.0, -3.0, 500)
    negative = -(10.0 ** rng.uniform(-3.0, 2.0, 500))
    negative = negative[np.abs(negative - np.round(negative)) > 1e-3]
    # Act / Assert
    np.testing.assert_allclose(special.log_gamma(large), sp.gammaln(large), rtol=1e-15, atol=0)
    np.testing.assert_allclose(special.log_gamma(tiny), sp.gammaln(tiny), rtol=1e-15, atol=0)
    np.testing.assert_allclose(special.log_gamma(near_zeros), sp.gammaln(near_zeros), rtol=0, atol=1e-15)
    # The reflection formula subtracts terms of order one, so negative arguments are compared relative to
    # the larger of 1 and the result
    reference = sp.gammaln(negative)
    error = np.abs(special.log_gamma(negative) - reference) / np.maximum(1.0, np.abs(reference))
    assert error.max() < 5e-15, error.max()


def test_log_gamma_keeps_its_relative_accuracy_next_to_the_zeros_at_1_and_2() -> None:
    # Arrange: ln Γ(1 + e) = -γe + (π²/12)e² + O(e³) and ln Γ(2 + e) = (1 - γ)e + (π²/12 - 1/2)e² + O(e³),
    # where the omitted terms are below 1e-30 relative for |e| = 1e-15 and 1e-14
    euler = 0.57721566490153286061
    e = np.array([1e-15, -1e-15, 1e-14, -1e-14])
    one = 1.0 + e
    two = 2.0 + e
    de_one = one - 1.0
    de_two = two - 2.0
    # Act / Assert: the arguments are rounded to doubles first, so the series is evaluated at the rounded offset
    np.testing.assert_allclose(special.log_gamma(one), -euler * de_one + (math.pi ** 2 / 12) * de_one ** 2,
                               rtol=1e-15, atol=0)
    np.testing.assert_allclose(special.log_gamma(two),
                               (1 - euler) * de_two + (math.pi ** 2 / 12 - 0.5) * de_two ** 2, rtol=1e-15, atol=0)


def test_log_gamma_returns_infinity_at_the_poles_and_nan_for_nan() -> None:
    # Act
    values = special.log_gamma(np.array([0.0, -0.0, -1.0, -2.0, -50.0, -1e20, np.inf, np.nan, -np.inf]))
    # Assert: scipy.special.gammaln also returns +inf at the poles and at +inf
    np.testing.assert_array_equal(values[:7], np.full(7, np.inf))
    np.testing.assert_array_equal(values[:7], sp.gammaln(np.array([0.0, -0.0, -1.0, -2.0, -50.0, -1e20, np.inf])))
    assert np.isnan(values[7]) and np.isnan(values[8])


@pytest.mark.parametrize("name", ["gamma_function", "log_gamma", "digamma", "trigamma", "stirling_gamma",
                                  "lanczos_gamma", "erf", "erfc", "debye1", "debye3", "bessel_i0", "bessel_i1",
                                  "bessel_i0e", "bessel_i1e", "bessel_k0", "bessel_k1", "bessel_j0", "bessel_j1",
                                  "bessel_y0", "bessel_y1"])
def test_array_path_returns_the_bits_of_the_scalar_member(name: str) -> None:
    # Arrange
    x = np.linspace(0.3, 7.0, 37).reshape(37, 1)
    wrapper = getattr(special, name)
    # Act
    values = wrapper(x)
    # Assert
    assert values.shape == x.shape
    np.testing.assert_array_equal(_bits(values.ravel()), _bits([wrapper(float(v)) for v in x.ravel()]))
    assert isinstance(wrapper(1.5), float)


def test_gamma_function_matches_the_net_member_bit_for_bit() -> None:
    # Arrange
    gamma = _sf().Gamma
    # Act / Assert
    np.testing.assert_array_equal(_bits(special.gamma_function(X)), _bits([gamma.Function(float(v)) for v in X]))


def test_stirling_is_an_asymptotic_approximation() -> None:
    # Arrange
    x = np.linspace(10.0, 100.0, 30)
    # Act / Assert
    np.testing.assert_allclose(special.stirling_gamma(x), sp.gamma(x), rtol=2e-11)


@pytest.mark.parametrize(("x", "expected"), [
    (0.5, 0.4795001221869535), (2.0, 0.004677734981047266), (3.0, 2.209049699858544e-05),
    (3.5, 7.430983723414128e-07), (6.0, 2.1519736712498913e-17), (10.0, 2.088487583762545e-45),
    (15.0, 7.212994172451207e-100), (26.0, 5.663192408856143e-296), (-2.0, 1.9953222650189528)])
def test_erfc_matches_mpmath_to_rounding(x: float, expected: float) -> None:
    # Arrange: mpmath erfc at 40 digits; scipy itself errs by up to about 6e-14 near x = 26
    # Act / Assert
    assert special.erfc(x) == pytest.approx(expected, rel=1e-15)


def test_regularized_incomplete_gamma_matches_scipy() -> None:
    # Arrange
    a = np.array([0.5, 1.5, 3.0, 10.0, 50.0])[:, None]
    x = np.linspace(0.01, 80.0, 60)[None, :]
    # Act
    lower = special.lower_incomplete_gamma(a, x)
    upper = special.upper_incomplete_gamma(a, x)
    # Assert
    assert lower.shape == (5, 60)
    # scipy itself errs by up to about 3.5e-14 at a = 50 (mpmath), and both err by about |ln P| units of
    # roundoff in the far tails, so the bulk is compared to 5e-14 and the tails to 1e-13
    for values, expected in ((lower, sp.gammainc(a, x)), (upper, sp.gammaincc(a, x))):
        bulk = expected > 1e-20
        np.testing.assert_allclose(values[bulk], expected[bulk], rtol=5e-14)
        np.testing.assert_allclose(values[~bulk], expected[~bulk], rtol=1e-13, atol=1e-300)
    assert special.lower_incomplete_gamma(2.0, -1.0) == 0.0
    assert special.upper_incomplete_gamma(2.0, -1.0) == 1.0


def test_inverse_incomplete_gamma_matches_scipy() -> None:
    # Arrange
    a = np.array([0.5, 1.5, 3.0, 10.0, 50.0])[:, None]
    y = np.linspace(0.001, 0.999, 50)[None, :]
    # Act / Assert
    np.testing.assert_allclose(special.inverse_lower_incomplete_gamma(a, y), sp.gammaincinv(a, y), rtol=5e-14)
    np.testing.assert_allclose(special.inverse_upper_incomplete_gamma(a, y), sp.gammainccinv(a, y), rtol=5e-14)
    # An upper probability near one is inverted through the exact complement; inverting Q directly errs by 2.6e-5.
    assert special.inverse_upper_incomplete_gamma(0.1, 1 - 1e-10) == pytest.approx(sp.gammainccinv(0.1, 1 - 1e-10),
                                                                                    rel=1e-13)
    assert special.inverse_lower_incomplete_gamma(2.0, 0.0) == 0.0
    assert special.inverse_lower_incomplete_gamma(2.0, 1.0) == math.inf


def test_as239_incomplete_gamma_takes_a_then_x_and_meets_its_tolerance() -> None:
    # Arrange
    a = np.array([0.5, 3.0, 10.0, 46.41588833612782])[:, None]
    x = np.linspace(0.5, 40.0, 30)[None, :]
    # Act
    values = special.incomplete_gamma(a, x)
    # Assert: the value of LowerIncomplete, nonzero below the exp(-180) floor of algorithm AS 239
    np.testing.assert_allclose(values, sp.gammainc(a, x), rtol=1e-13, atol=1e-300)
    np.testing.assert_array_equal(_bits(values), _bits(special.lower_incomplete_gamma(a, x)))
    assert special.incomplete_gamma(46.41588833612782, 0.004641588833612782) == pytest.approx(
        1.8173331985826257e-167, rel=1e-13)
    assert special.incomplete_gamma(3.0, 1.0) == _sf().Gamma.Incomplete(1.0, 3.0)
    with pytest.raises(ValueError, match="out of range"):
        special.incomplete_gamma(-1.0, 1.0)


def test_beta_functions_match_scipy() -> None:
    # Arrange
    a = np.array([0.3, 1.0, 2.5, 7.0])[:, None, None]
    b = np.array([0.4, 1.0, 3.0, 11.0])[None, :, None]
    x = np.linspace(0.01, 0.99, 25)[None, None, :]
    # Act
    complete = special.beta_function(a[..., 0], b[..., 0])
    ratio = special.beta_incomplete_ratio(a, b, x, sp.betaln(a, b))
    # Assert
    np.testing.assert_allclose(complete, sp.beta(a[..., 0], b[..., 0]), rtol=1e-13)
    np.testing.assert_allclose(ratio, sp.betainc(a, b, x), rtol=2e-14, atol=1e-300)
    np.testing.assert_allclose(special.beta_incomplete(2.0, 3.0, x.ravel()), sp.betainc(2.0, 3.0, x.ravel()),
                               rtol=4e-15)
    # The shapes broadcast like scipy.special.betainc, and an array crosses once with the scalar bits
    grid = special.beta_incomplete(a, b, x)
    np.testing.assert_allclose(grid, sp.betainc(a, b, x), rtol=2e-14, atol=1e-300)
    scalar = special.beta_incomplete(float(a.ravel()[0]), float(b.ravel()[0]), float(x.ravel()[3]))
    assert np.float64(scalar).view(np.uint64) == grid.reshape(-1)[3].view(np.uint64)
    # The inverse recovers x wherever I_x(a, b) is not flat against 0 or 1
    inverse = special.beta_incomplete_inverse(a, b, grid)
    resolved = (grid > 1e-6) & (grid < 1.0 - 1e-6)
    np.testing.assert_allclose(inverse[resolved], np.broadcast_to(x, grid.shape)[resolved], rtol=1e-10)


def test_log_beta_matches_scipy_betaln() -> None:
    # Arrange: shapes from 1e-3 to 1e7 on both sides of the branch points 1, 2 and 8 of algorithm 708
    rng = np.random.default_rng(708)
    a = 10.0 ** rng.uniform(-3.0, 7.0, 3000)
    b = 10.0 ** rng.uniform(-3.0, 7.0, 3000)
    # Act
    values = special.log_beta(a, b)
    reference = sp.betaln(a, b)
    # Assert: scipy (Cephes lbeta) is itself accurate only to a few 1e-9 relative to max(1, |ln B|) at the
    # largest shapes, so this checks agreement; the mpmath test below checks the accuracy
    error = np.abs(values - reference) / np.maximum(1.0, np.abs(reference))
    assert error.max() < 1e-8, error.max()
    small = (a < 8.0) & (b < 8.0)
    assert error[small].max() < 5e-15, error[small].max()
    np.testing.assert_array_equal(special.log_beta(b, a), values)


def test_log_beta_matches_scipy_betaln_for_small_first_shape() -> None:
    # Arrange: one shape below 1 and the other in (2, 8), the branch that reduces the larger shape to (1, 2]
    rng = np.random.default_rng(7085)
    a = rng.uniform(1e-3, 1.0, 3000)
    b = rng.uniform(2.0, 8.0, 3000)
    # Act
    values = special.log_beta(a, b)
    # Assert: scipy betaln is within about 2.3e-15 of mpmath here and log_beta within about 7e-16
    error = np.abs(values - sp.betaln(a, b)) / np.maximum(1.0, np.abs(sp.betaln(a, b)))
    assert error.max() < 3e-15, error.max()


@pytest.mark.parametrize(("a", "b", "expected"), [
    # mpmath 1.3.0 at 400 digits: loggamma(a) + loggamma(b) - loggamma(a + b)
    (0.5, 1e6, -6.335390211057437),
    (0.5, 1e15, -16.69702325453064),
    (3.7, 1e7, -58.208882081380345),
    (1.5, 1500.25, -11.091112700256792),
    (1000.0, 1000.0, -1388.4826016359023),
    (1e5, 1e5, -138633.92706134808),
    (8.0, 1e10, -175.68164608125824),
    (1e8, 0.25, -3.317147660352514),
    (2.5e7, 3e7, -37895515.41140297),
    (12.5, 1e300, -8615.959751215734),
    # mpmath 1.4.1 at 400 digits, one shape below 1 and the other in (2, 8)
    (0.9219, 7.5939, -1.8140586851054754),
    (0.7052453993501302, 7.681687846760826, -1.16997303546981),
    (0.8739311446570729, 7.6339892317246765, -1.6825461302652165),
    (0.01, 7.99, 4.5793301230650885),
])
def test_log_beta_matches_mpmath_at_large_and_asymmetric_shapes(a: float, b: float, expected: float) -> None:
    # Act
    value = special.log_beta(a, b)
    # Assert: the ln Gamma difference is in error by 1e-9 at (0.5, 1e6) and returns -16 at (0.5, 1e15)
    assert value == pytest.approx(expected, rel=0, abs=1e-15 * max(1.0, abs(expected)))


def test_log_beta_special_values_follow_scipy() -> None:
    # Act / Assert
    assert special.log_beta(-0.5, 2.0) == pytest.approx(sp.betaln(-0.5, 2.0), rel=1e-15)
    assert special.log_beta(np.inf, 2.0) == sp.betaln(np.inf, 2.0) == -np.inf
    assert np.isnan(special.log_beta(np.nan, 2.0))
    assert special.beta_function(40.0, 60.0) == math.exp(special.log_beta(40.0, 60.0))


def test_beta_series_members_reproduce_the_incomplete_beta_where_they_apply() -> None:
    # Arrange
    a, b, x = 2.0, 3.0, np.array([0.05, 0.1, 0.2])
    # Act
    series = special.beta_power_series(a, b, x)
    # Assert
    np.testing.assert_allclose(series, sp.betainc(a, b, x), rtol=1e-14)
    beta = _sf().Beta
    assert special.beta_continued_fraction1(a, b, 0.3) == beta.Incbcf(a, b, 0.3)
    assert special.beta_continued_fraction2(a, b, 0.3) == beta.Incbd(a, b, 0.3)


def test_inverse_erf_refuses_out_of_range_values() -> None:
    # Act / Assert
    with pytest.raises(ValueError):
        special.inverse_erf(1.5)
    assert special.inverse_erf(1.0) == math.inf


@pytest.mark.parametrize(("wrapper", "oracle"), [
    (special.bessel_j0, sp.j0), (special.bessel_j1, sp.j1), (special.bessel_y0, sp.y0), (special.bessel_y1, sp.y1),
])
def test_order_zero_and_one_bessel_functions_match_scipy(wrapper, oracle) -> None:
    # Arrange: below 8 relative (away from the zeros, where scipy is only absolutely accurate), above 8 absolute
    small = np.linspace(0.1, 8.0, 150)
    small = small[np.abs(oracle(small)) > 0.05]
    large = np.linspace(8.0, 400.0, 300)
    # Act / Assert
    np.testing.assert_allclose(wrapper(small), oracle(small), rtol=2e-14)
    np.testing.assert_allclose(wrapper(large), oracle(large), rtol=0, atol=2e-15)


@pytest.mark.parametrize(("wrapper", "zeros", "expected"), [
    (special.bessel_j0, [2.404825557695773, 5.520078110286311], [-5.191473894843614e-07, 3.402647757577026e-07]),
    (special.bessel_y0, [0.8935769662791675, 3.957678419314858, 7.086051060301773],
     [8.794203104203801e-07, -4.0254262101865376e-07, 3.000975938036521e-07]),
    (special.bessel_y1, [2.197141326031017, 5.429681040794135], [5.207862939856164e-07, -3.403180140906233e-07]),
])
def test_bessel_functions_keep_relative_accuracy_next_to_their_zeros(wrapper, zeros, expected) -> None:
    # Arrange: 1e-6 above a zero below 8; the expected values are mpmath at 40 digits. scipy (Cephes) is
    # accurate there only absolutely, with relative errors up to about 1e-10.
    x = np.array(zeros) + 1e-6
    # Act
    values = wrapper(x)
    # Assert
    np.testing.assert_allclose(values, expected, rtol=1e-15)


@pytest.mark.parametrize(("wrapper", "oracle"), [
    (special.bessel_in, sp.iv), (special.bessel_kn, sp.kn), (special.bessel_jn, sp.jv), (special.bessel_yn, sp.yn),
])
def test_integer_order_bessel_functions_broadcast_order_and_argument(wrapper, oracle) -> None:
    # Arrange
    n = np.array([2, 3, 5])[:, None]
    x = np.linspace(0.5, 40.0, 80)[None, :]
    # Act
    values = wrapper(n, x)
    # Assert
    assert values.shape == (3, 80)
    np.testing.assert_allclose(values, oracle(n, x), rtol=1e-13, atol=1e-15)
    assert values[1, 4] == wrapper(3, float(x[0, 4]))
    np.testing.assert_array_equal(_bits(values), _bits([[wrapper(int(o), float(v)) for v in x[0]] for o in n[:, 0]]))


def _debye(order: int, x: float) -> float:
    integral = integrate.quad(lambda t: t**order / math.expm1(t), 0.0, x, epsabs=0.0, epsrel=1e-13, limit=200)[0]
    return order / x**order * integral


def test_debye_functions_match_quadrature() -> None:
    # Arrange
    x = np.array([0.1, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
    # Act / Assert
    np.testing.assert_allclose(special.debye1(x), [_debye(1, v) for v in x], rtol=1e-14)
    np.testing.assert_allclose(special.debye3(x), [_debye(3, v) for v in x], rtol=1e-14)


def test_factorials_and_binomial_coefficients_match_exact_values() -> None:
    # Arrange: the oracle is the exact integer rounded to double. SciPy's floating factorial is a gamma
    # approximation whose accuracy depends on the SciPy version and platform (1.17 on Windows: 8e-14).
    n = np.arange(0, 171)
    # Act / Assert
    np.testing.assert_allclose(special.factorial(n), [float(math.factorial(k)) for k in n], rtol=2e-15)
    np.testing.assert_array_equal(special.factorial(np.arange(0, 23)), [math.factorial(k) for k in range(23)])
    assert special.factorial(171) == math.inf
    np.testing.assert_allclose(special.log_factorial(np.arange(2, 500)), sp.gammaln(np.arange(3, 501)), rtol=1e-15)
    nn, kk = np.meshgrid(np.arange(60), np.arange(60), indexing="ij")
    expected = np.array([[float(math.comb(a, b)) for b in range(60)] for a in range(60)])
    np.testing.assert_allclose(special.binomial_coefficient(nn, kk), expected, rtol=1e-15)
    assert special.binomial_coefficient(52, 5) == math.comb(52, 5)
    with pytest.raises(ValueError):
        special.factorial(-1)


def test_combinations_match_itertools() -> None:
    # Arrange
    import itertools  # noqa: PLC0415

    # Act
    found = special.find_combinations(2, 4)
    table = special.all_combinations(3)
    lazy = [c.tolist() for c in special.iter_all_combinations(3)]
    advanced, successor = special.next_combination(np.array([0, 3]), 4)
    last, _ = special.next_combination([2, 3], 4)
    # Assert
    np.testing.assert_array_equal(found, list(itertools.combinations(range(4), 2)))
    assert table.shape == (7, 3) and table.dtype == np.int32
    assert lazy == [list(np.flatnonzero(row)) for row in table]
    assert advanced and successor.tolist() == [1, 2]
    assert not last


def test_polynomials_match_numpy() -> None:
    # Arrange
    c = np.array([1.5, -2.0, 0.25, 3.0])
    x = np.linspace(-2.0, 2.0, 21)
    # Act / Assert
    np.testing.assert_allclose(special.polynomial(c, x), np.polynomial.polynomial.polyval(x, c), rtol=1e-15,
                               atol=1e-15)
    np.testing.assert_array_equal(special.polynomial_reversed(c, x), np.polyval(c, x))
    np.testing.assert_array_equal(special.polynomial_reversed(c, x, 2), np.polyval(c[:3], x))
    np.testing.assert_array_equal(special.polynomial_reversed_monic(c[1:], x), np.polyval(np.r_[1.0, c[1:]], x))
    with pytest.raises(ValueError):
        special.polynomial_reversed(c, 1.0, 9)


def test_chebyshev_matches_numpy_chebval() -> None:
    # Arrange: Evaluate.Chebyshev sums c0/2 + sum c_k T_k(x/2), numpy chebval sums c0 + sum c_k T_k(s)
    c = np.array([1.5, -2.0, 0.25, 3.0, -0.125])
    x = np.linspace(-2.0, 2.0, 41)
    halved = np.r_[c[0] / 2.0, c[1:]]
    # Act / Assert
    np.testing.assert_allclose(special.chebyshev(c, x), np.polynomial.chebyshev.chebval(x / 2.0, halved), rtol=1e-14,
                               atol=1e-14)
    np.testing.assert_allclose(special.chebyshev(c, x, 3), np.polynomial.chebyshev.chebval(x / 2.0, halved[:3]),
                               rtol=1e-14, atol=1e-14)
    assert special.chebyshev(c, 2.0) == float(np.sum(halved))
    with pytest.raises(ValueError):
        special.chebyshev(c, 1.0, 0)
    with pytest.raises(ValueError):
        special.chebyshev(c, 1.0, 6)


def test_gamma_max_and_poles() -> None:
    # Act / Assert
    assert special.gamma_max() == float(_sf().Gamma.GammaMax)
    assert math.isfinite(special.gamma_function(special.gamma_max() - 0.01))
    with pytest.raises(ArithmeticError):
        special.gamma_function(0.0)


def test_scalar_special_functions_are_bound_once() -> None:
    # Arrange
    x = np.linspace(0.5, 5.0, 2000)
    special.log_gamma(x)
    # Act
    start = time.perf_counter()
    special.log_gamma(x)
    elapsed = time.perf_counter() - start
    # Assert: a bound overload costs a few microseconds a call; resolving the group each call costs far more
    assert elapsed < 0.2


def test_array_arguments_cross_into_the_list_overloads_once() -> None:
    # Arrange: 20000 arguments; one crossing into the IList<double> overload, against one crossing per value
    x = np.linspace(0.5, 30.0, 20000)
    a = np.full_like(x, 3.5)
    special.bessel_j0(x)
    special.lower_incomplete_gamma(a, x)
    gamma = _sf().Gamma
    j0 = _sf().Bessel.J0.Overloads[float]
    # Act
    start = time.perf_counter()
    vector_j0 = special.bessel_j0(x)
    vector_p = special.lower_incomplete_gamma(a, x)
    vectorized = time.perf_counter() - start
    start = time.perf_counter()
    scalar_j0 = [j0(float(v)) for v in x]
    scalar_p = [gamma.LowerIncomplete(3.5, float(v)) for v in x]
    per_value = time.perf_counter() - start
    # Assert: the same bits, at least five times faster than a crossing per value
    np.testing.assert_array_equal(_bits(vector_j0), _bits(scalar_j0))
    np.testing.assert_array_equal(_bits(vector_p), _bits(scalar_p))
    assert vectorized * 5 < per_value


def test_integer_order_bessel_arrays_cross_once_per_order() -> None:
    # Arrange
    n = np.array([[2], [7], [2]])
    x = np.linspace(0.5, 20.0, 9)[None, :]
    # Act
    values = special.bessel_jn(n, x)
    # Assert: rows 0 and 2 share an order, and every value is the scalar member's
    assert values.shape == (3, 9)
    np.testing.assert_array_equal(_bits(values[0]), _bits(values[2]))
    np.testing.assert_array_equal(_bits(values), _bits([[special.bessel_jn(int(o), float(v)) for v in x[0]]
                                                         for o in n[:, 0]]))
