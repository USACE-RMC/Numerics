using System;

namespace Numerics.Distributions
{
    /// <summary>Provides range-preserving covariance contractions and divided-exponential products for distribution uncertainty calculations.</summary>
    internal static partial class DistributionNumerics
    {
        /// <summary>Contracts an unchanged covariance and quantile gradient without forming avoidable overflowing products.</summary>
        /// <param name="covariance">The finite covariance in the supplied gradient coordinates.</param>
        /// <param name="gradient">The corresponding quantile gradient.</param>
        /// <param name="scale">An optional finite positive common quantile scale.</param>
        /// <returns>The nonnegative scalar variance, including true floating-point underflow or positive overflow.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Scale or dimensions are invalid.</exception>
        /// <exception cref="InvalidOperationException">The covariance, gradient or contracted variance is unresolved.</exception>
        /// <remarks>Forms each signed covariance-gradient product in logarithms, normalizes by the
        /// largest term, and uses compensated summation. No eigenvalue repair, clipping,
        /// diagonal inflation or negative-variance floor is applied.</remarks>
        internal static double ScaledQuantileVariance(double[,] covariance, double[] gradient, double scale = 1)
        {
            if (!(scale > 0) || !Tools.IsFinite(scale)) throw new ArgumentOutOfRangeException(nameof(scale));
            if (covariance.GetLength(0) != gradient.Length || covariance.GetLength(1) != gradient.Length)
                throw new ArgumentOutOfRangeException(nameof(covariance), "Covariance and gradient dimensions must agree.");
            var gradientLogs = new double[gradient.Length];
            for (int i = 0; i < gradient.Length; i++)
            {
                if (double.IsNaN(gradient[i])) throw new InvalidOperationException("The quantile gradient is undefined.");
                if (double.IsInfinity(gradient[i]))
                    throw new InvalidOperationException("An unrepresentable gradient requires finite physical coordinates before scalar variance contraction.");
                gradientLogs[i] = Math.Log(Math.Abs(gradient[i]));
                for (int j = 0; j < gradient.Length; j++)
                {
                    if (!Tools.IsFinite(covariance[i, j])) throw new InvalidOperationException("The parameter covariance is outside the finite floating-point range.");
                }
            }
            double largestLogTerm = double.NegativeInfinity;
            for (int i = 0; i < gradient.Length; i++)
            for (int j = 0; j < gradient.Length; j++)
            {
                if (covariance[i, j] == 0 || gradient[i] == 0 || gradient[j] == 0) continue;
                largestLogTerm = Math.Max(largestLogTerm, Math.Log(Math.Abs(covariance[i, j])) + gradientLogs[i] + gradientLogs[j]);
            }
            if (double.IsNegativeInfinity(largestLogTerm)) return 0;
            double quadratic = 0, correction = 0;
            for (int i = 0; i < gradient.Length; i++)
            for (int j = 0; j < gradient.Length; j++)
            {
                if (covariance[i, j] == 0 || gradient[i] == 0 || gradient[j] == 0) continue;
                double logTerm = Math.Log(Math.Abs(covariance[i, j])) + gradientLogs[i] + gradientLogs[j];
                double term = Math.Sign(covariance[i, j]) * Math.Sign(gradient[i]) * Math.Sign(gradient[j])
                    * Math.Exp(logTerm - largestLogTerm);
                double next = quadratic + term;
                correction += Math.Abs(quadratic) >= Math.Abs(term) ? (quadratic - next) + term : (term - next) + quadratic;
                quadratic = next;
            }
            quadratic += correction;
            if (!Tools.IsFinite(quadratic) || quadratic < 0)
                throw new InvalidOperationException("The quantile variance could not be resolved as a nonnegative quadratic form.");
            if (quadratic == 0) return 0;
            return Math.Exp(2 * Math.Log(scale) + largestLogTerm + Math.Log(quadratic));
        }

        /// <summary>Returns shape squared times the positive Gamma Fisher residual trigamma(shape)-1/shape.</summary>
        /// <param name="shape">The finite positive Gamma shape.</param>
        /// <returns>The scaled Fisher residual without subtracting rounded reciprocal leading terms.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Shape is not finite and positive.</exception>
        /// <remarks>Uses the trigamma recurrence below shape 32 and the same Bernoulli expansion as
        /// AccurateTrigamma after analytically removing 1/shape. Scaling before evaluation also
        /// avoids squaring a large shape or forming an underflowed unscaled residual.</remarks>
        internal static double GammaScaledFisherResidual(double shape)
        {
            if (!(shape > 0) || !Tools.IsFinite(shape)) throw new ArgumentOutOfRangeException(nameof(shape));
            double shifted = shape, recurrence = 0;
            while (shifted < 32)
            {
                double ratio = shape / shifted;
                recurrence += ratio * ratio;
                shifted++;
            }
            double r = 1 / shifted, s = r * r;
            double residual = .5 + r * (1.0 / 6 + s * (-1.0 / 30 + s * (1.0 / 42
                + s * (-1.0 / 30 + s * (5.0 / 66 - s * 691.0 / 2730)))));
            if (shifted == shape) return residual;
            double scaleRatio = shape / shifted;
            return recurrence - shape + shape * scaleRatio + scaleRatio * scaleRatio * residual;
        }

        /// <summary>Forms scale times value times exprel(argument), retaining a finite product after unit-scale overflow.</summary>
        /// <param name="scale">The signed outer multiplier.</param>
        /// <param name="value">The signed value multiplying the divided exponential.</param>
        /// <param name="argument">The argument of the exponential relative function.</param>
        /// <returns><paramref name="scale"/> times <paramref name="value"/> times <c>exprel(argument)</c>, evaluated with extended logarithmic range when necessary.</returns>
        internal static double ScaledExprelProduct(double scale, double value, double argument)
        {
            return ScaledExprelProductCore(scale, value, argument, false);
        }

        /// <summary>Forms scale times value squared times exprel'(argument) without squaring a tiny value prematurely.</summary>
        /// <param name="scale">The signed outer multiplier.</param>
        /// <param name="value">The value whose square multiplies the derivative.</param>
        /// <param name="argument">The argument of the exponential relative derivative.</param>
        /// <returns><paramref name="scale"/> times the square of <paramref name="value"/> times <c>exprel'(argument)</c>, evaluated with extended logarithmic range when necessary.</returns>
        internal static double ScaledExprelDerivativeProduct(double scale, double value, double argument)
        {
            return ScaledExprelProductCore(scale, value, argument, true);
        }

        /// <summary>Combines signed multipliers and divided-exponential logarithms when ordinary product arithmetic loses range.</summary>
        /// <param name="scale">The signed outer multiplier.</param>
        /// <param name="value">The signed value, or the value to square when <paramref name="derivative"/> is <see langword="true"/>.</param>
        /// <param name="argument">The divided-exponential argument.</param>
        /// <param name="derivative"><see langword="true"/> to use the derivative and square <paramref name="value"/>; otherwise, <see langword="false"/> to use the function itself.</param>
        /// <returns>The requested signed product, including its natural zero, infinity, or not-a-number limit.</returns>
        private static double ScaledExprelProductCore(double scale, double value, double argument, bool derivative)
        {
            if (value == 0) return 0;
            if (double.IsPositiveInfinity(argument)) return derivative ? double.PositiveInfinity : Math.Sign(value) * double.PositiveInfinity;
            if (double.IsNegativeInfinity(argument)) return 0;
            double divided = derivative ? ExprelDerivative(argument) : Exprel(argument);
            double factor = derivative ? value * value : value;
            double result = scale * (factor * divided);
            if (Tools.IsFinite(result) && result != 0) return result;
            double logDivided;
            if (argument > 50)
                logDivided = derivative ? argument + Math.Log(argument - 1) + Tools.Log1p(Math.Exp(-argument) / (argument - 1)) - 2 * Math.Log(argument)
                    : argument + Tools.Log1p(-Math.Exp(-argument)) - Math.Log(argument);
            else if (argument < -50)
                logDivided = derivative ? Tools.Log1p((argument - 1) * Math.Exp(argument)) - 2 * Math.Log(-argument)
                    : Math.Log(-Tools.Expm1(argument)) - Math.Log(-argument);
            else logDivided = Math.Log(divided);
            double logarithm = Math.Log(scale) + (derivative ? 2 : 1) * Math.Log(Math.Abs(value)) + logDivided;
            return (derivative ? 1 : Math.Sign(value)) * Math.Exp(logarithm);
        }
    }
}
