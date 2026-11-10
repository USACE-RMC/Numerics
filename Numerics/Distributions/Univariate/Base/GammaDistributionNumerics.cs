using System;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{
    /// <summary>Provides stable gamma-function tails, quantiles, derivatives, and asymptotic expansions for univariate distributions.</summary>
    internal static partial class DistributionNumerics
    {
        /// <summary>Log regularized lower gamma integral P(a,x).</summary>
        /// <param name="shape">The finite positive gamma shape.</param>
        /// <param name="x">The nonnegative unit-scale gamma variate.</param>
        /// <returns>The natural logarithm of the regularized lower incomplete gamma integral, including endpoint and not-a-number limits.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="shape"/> is not finite and positive.</exception>
        /// <remarks>The value is the logarithmic form of the incomplete gamma core of <see cref="Gamma.LowerIncomplete(double, double)"/>.</remarks>
        internal static double GammaLogCDF(double shape, double x) => GammaLogTailValue(shape, x, false);

        /// <summary>Log regularized upper gamma integral Q(a,x).</summary>
        /// <param name="shape">The finite positive gamma shape.</param>
        /// <param name="x">The nonnegative unit-scale gamma variate.</param>
        /// <returns>The natural logarithm of the regularized upper incomplete gamma integral, including endpoint and not-a-number limits.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="shape"/> is not finite and positive.</exception>
        /// <remarks>The value is the logarithmic form of the incomplete gamma core of <see cref="Gamma.UpperIncomplete(double, double)"/>.</remarks>
        internal static double GammaLogSurvival(double shape, double x) => GammaLogTailValue(shape, x, true);

        /// <summary>Evaluates a gamma log tail with the incomplete gamma core of <see cref="Gamma"/>.</summary>
        /// <param name="a">The finite positive gamma shape.</param>
        /// <param name="x">The unit-scale gamma variate.</param>
        /// <param name="upper"><see langword="true"/> to evaluate the upper tail; <see langword="false"/> to evaluate the lower tail.</param>
        /// <returns>The requested log tail, including endpoint and not-a-number limits.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="a"/> is not finite and positive.</exception>
        /// <exception cref="ArithmeticException">The selected series or continued fraction does not converge.</exception>
        /// <remarks>
        /// The regions, Temme expansion and factor are those of <see cref="Gamma.LowerIncomplete(double, double)"/>, with the
        /// factor x^a e^(-x) / Γ(a) carried as its logarithm, so tails below the double range keep finite logarithms.
        /// </remarks>
        private static double GammaLogTailValue(double a, double x, bool upper)
        {
            ValidateGammaShape(a);
            if (double.IsNaN(x)) return double.NaN;
            if (x <= 0) return upper ? 0 : double.NegativeInfinity;
            if (double.IsPositiveInfinity(x)) return upper ? double.NegativeInfinity : 0;
            return Gamma.LogIncompleteRatio(a, x, upper);
        }

        /// <summary>Log unit-scale gamma density, including one-sided endpoint limits.</summary>
        /// <param name="shape">The finite positive gamma shape.</param>
        /// <param name="x">The unit-scale gamma variate.</param>
        /// <returns>The log density, negative infinity outside support, the appropriate zero endpoint limit, or not-a-number when <paramref name="x"/> is not-a-number.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="shape"/> is not finite and positive.</exception>
        internal static double GammaLogDensity(double shape, double x)
        {
            ValidateGammaShape(shape);
            if (double.IsNaN(x)) return double.NaN;
            if (x < 0 || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            if (x == 0) return shape == 1 ? 0 : shape < 1 ? double.PositiveInfinity : double.NegativeInfinity;
            return shape == 1 ? -x : GammaLogKernel(shape, x) - Math.Log(x);
        }

        /// <summary>Implicit shape derivative of the actual unit-scale gamma quantile at its value.</summary>
        /// <param name="shape">The finite positive gamma shape.</param>
        /// <param name="unitQuantile">The nonnegative finite unit-scale gamma quantile.</param>
        /// <returns>The derivative of the unit-scale quantile with respect to <paramref name="shape"/>, including zero at the lower endpoint.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="shape"/> is invalid or <paramref name="unitQuantile"/> is negative or nonfinite.</exception>
        /// <exception cref="ArithmeticException">A required gamma series or continued fraction does not converge.</exception>
        /// <remarks>The fixed-x shape derivative of the smaller log tail comes from
        /// <see cref="Gamma.LogIncompleteRatio(double, double, bool, out double)"/>, which differentiates its series,
        /// continued fraction or Temme expansion together with the probability; dividing by the density, formed from the same factor x^a e^(-x) / Γ(a) as the tail so that its rounding cancels, gives the
        /// implicit quantile derivative. No frequency-factor approximation or perturbation across the shape boundary is used.</remarks>
        internal static double GammaQuantileShapeDerivative(double shape, double unitQuantile)
        {
            ValidateGammaShape(shape);
            if (unitQuantile == 0) return 0;
            if (!(unitQuantile > 0) || !Tools.IsFinite(unitQuantile))
                throw new ArgumentOutOfRangeException(nameof(unitQuantile));
            bool upper = unitQuantile >= shape;
            double log = Gamma.LogIncompleteRatio(shape, unitQuantile, upper, out double derivative);
            if (derivative == 0) return 0;
            return (upper ? Math.Sign(derivative) : -Math.Sign(derivative))
                * Math.Exp(Math.Log(Math.Abs(derivative)) + log - (Gamma.LogIncompleteFactor(shape, unitQuantile) - Math.Log(unitQuantile)));
        }

        /// <summary>Requires a positive finite gamma shape.</summary>
        /// <param name="shape">The shape to validate.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="shape"/> is not finite and positive.</exception>
        private static void ValidateGammaShape(double shape)
        {
            if (!(shape > 0) || !Tools.IsFinite(shape)) throw new ArgumentOutOfRangeException(nameof(shape), "Gamma shape must be positive and finite.");
        }

        /// <summary>Log(x^a exp(-x)/Gamma(a)) without subtracting large near-equal terms.</summary>
        /// <param name="a">The positive gamma shape.</param>
        /// <param name="x">The positive unit-scale gamma variate.</param>
        /// <returns><c>log(x^a*exp(-x)/Gamma(a))</c>.</returns>
        /// <remarks>
        /// This is the same quantity as <see cref="Gamma.LogIncompleteFactor(double, double)"/>, evaluated by different
        /// branches: the sum a ln x - x - ln Γ(a) for a &lt; 16 rather than the logarithm of the direct product up to
        /// a = 20, and above that the series for ln(1 + δ) - δ only for |δ| &lt; 0.25 rather than for every x &gt;= a / 2.
        /// The two therefore differ in the last bits, and where they differ the incomplete-factor form is usually the
        /// closer. They are kept separate because this kernel defines the gamma-family densities (gamma, Pearson type III,
        /// inverse gamma, inverse chi-squared, noncentral t and Poisson), whose values are reference results, while
        /// the incomplete factor is paired with the incomplete gamma ratios so that its rounding cancels in tail-to-density
        /// ratios.
        /// </remarks>
        private static double GammaLogKernel(double a, double x)
        {
            if (a < 16) return a * Math.Log(x) - x - Gamma.LogGamma(a);
            double delta = (x - a) / a;
            double deviance = Math.Abs(delta) < 0.25 ? a * Tools.Log1pMinusX(delta)
                : a * (Math.Log(x) - Math.Log(a)) - (x - a);
            return deviance + 0.5 * Math.Log(a) - Tools.LogSqrt2PI - Gamma.StirlingRemainder(a);
        }
    }
}
