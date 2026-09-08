using System;

namespace Numerics.Distributions
{
    /// <summary>
    /// Evaluates the Kappa Four lower-support residual with compensated arithmetic when
    /// ordinary logarithms lose the distance from a positive-hondo boundary.
    /// </summary>
    internal static class KappaFourBoundary
    {
        private static readonly Pair LogTwo = new Pair(0.6931471805599453d, 2.3190468138462996E-17d);

        /// <summary>Computes a positive-hondo lower endpoint with compensated logarithm, exponential and affine arithmetic.</summary>
        /// <param name="xi">The finite location.</param>
        /// <param name="alpha">The positive finite scale.</param>
        /// <param name="k">The finite kappa shape.</param>
        /// <param name="h">The positive finite hondo shape.</param>
        /// <param name="endpoint">The accurately rounded finite endpoint when evaluation succeeds.</param>
        /// <returns>Whether the compensated finite-range evaluation succeeded.</returns>
        /// <remarks>
        /// A rounded log(h) multiplied by kappa can move the support by several doubles.
        /// Retaining the low parts prevents the declared support from excluding its first
        /// interior representable argument. No tolerance or support clipping is used.
        /// Extreme exponential or affine ranges remain the caller's responsibility.
        /// </remarks>
        internal static bool TryLowerEndpoint(double xi, double alpha, double k, double h, out double endpoint)
        {
            endpoint = double.NaN;
            if (!(h > 0) || !Tools.IsFinite(h)) return false;
            Pair logH = Log(new Pair(h));
            Pair exponent = Multiply(new Pair(-k), logH);
            if (!Tools.IsFinite(exponent.High) || Math.Abs(exponent.High) > 700) return false;
            Pair standard;
            if (Math.Abs(exponent.High) < .5)
            {
                Pair relative = new Pair(1), term = new Pair(1);
                for (int n = 1; n <= 32; n++)
                {
                    term = Divide(Multiply(term, exponent), new Pair(n + 1));
                    relative = Add(relative, term);
                }
                standard = Multiply(logH, relative);
            }
            else
            {
                int power = (int)Math.Round(exponent.High / LogTwo.High);
                Pair reduced = Subtract(exponent, Multiply(new Pair(power), LogTwo));
                Pair sum = new Pair(1), term = new Pair(1);
                for (int n = 1; n <= 32; n++)
                {
                    term = Divide(Multiply(term, reduced), new Pair(n));
                    sum = Add(sum, term);
                }
                standard = Divide(Subtract(new Pair(1), Scale(sum, power)), new Pair(k));
            }
            endpoint = Add(new Pair(xi), Multiply(new Pair(alpha), standard)).High;
            return Tools.IsFinite(endpoint);
        }

        /// <summary>Returns log F close to the lower support endpoint for positive hondo.</summary>
        /// <param name="x">The value being evaluated.</param>
        /// <param name="xi">The finite location parameter.</param>
        /// <param name="alpha">The positive finite scale parameter.</param>
        /// <param name="k">The finite kappa shape.</param>
        /// <param name="h">The positive finite hondo shape.</param>
        /// <returns>The log probability, including negative infinity at or below the endpoint.</returns>
        /// <remarks>
        /// Retains the low components of the standardized subtraction and logarithms. The
        /// logarithm uses an atanh series after binary range reduction; its truncation error
        /// is below the precision of a compensated pair. No endpoint clipping is applied.
        /// </remarks>
        internal static double LogProbability(double x, double xi, double alpha, double k, double h)
        {
            Pair y = StandardizedDifference(x, xi, alpha);
            if (y.High == 0d && h == 1d && x > xi)
            {
                // The standardized distance can underflow although its logarithm is finite.
                Pair difference = Subtract(new Pair(x), new Pair(xi));
                return Subtract(Log(difference), Log(new Pair(alpha))).High;
            }

            Pair logH = Log(new Pair(h));
            Pair z;
            if (k == 0d)
            {
                z = Subtract(logH, y);
            }
            else
            {
                Pair product = Multiply(new Pair(-k), y);
                if (Math.Abs(product.High) < 0.125d)
                {
                    // log(1-k*y)/k = -y*log1p(-k*y)/(-k*y), including k -> 0.
                    z = Subtract(logH, Multiply(y, LogOnePlusRelative(product)));
                }
                else
                {
                    Pair logBase;
                    if (double.IsPositiveInfinity(product.High))
                    {
                        // Here the reciprocal-product correction is below pair precision.
                        logBase = Add(Log(new Pair(Math.Abs(k))), Log(Abs(y)));
                    }
                    else
                    {
                        Pair basis = Add(new Pair(1d), product);
                        if (basis.High <= 0d) return 0d;
                        logBase = Log(basis);
                    }
                    z = Divide(Add(logBase, Multiply(new Pair(k), logH)), new Pair(k));
                }
            }

            if (z.High >= 0d) return double.NegativeInfinity;
            return Math.Log(-Tools.Expm1(z.High)) / h;
        }

        /// <summary>Evaluates the latent logarithm close to an endpoint where kappa times the standardized value approaches one.</summary>
        /// <param name="x">The value being evaluated.</param>
        /// <param name="xi">The finite location parameter.</param>
        /// <param name="alpha">The positive finite scale parameter.</param>
        /// <param name="k">The finite kappa shape.</param>
        /// <returns>Logarithm of the latent transform, with its continuous zero-kappa limit.</returns>
        /// <remarks>Retains subtraction and product corrections in the small positive base 1-k*(x-xi)/alpha.</remarks>
        internal static double LogT(double x, double xi, double alpha, double k)
        {
            Pair y = StandardizedDifference(x, xi, alpha);
            if (k == 0d) return -y.High;
            Pair basis = Subtract(new Pair(1d), Multiply(new Pair(k), y));
            if (basis.High <= 0d) return Math.Log(basis.High) / k;
            return Divide(Log(basis), new Pair(k)).High;
        }

        /// <summary>Forms the standardized difference in pairs, halving before a subtraction that would overflow.</summary>
        private static Pair StandardizedDifference(double x, double xi, double alpha)
        {
            Pair difference = Subtract(new Pair(x), new Pair(xi));
            if (double.IsInfinity(difference.High))
            {
                difference = Subtract(new Pair(x * 0.5d), new Pair(xi * 0.5d));
                return Multiply(Divide(difference, new Pair(alpha)), new Pair(2d));
            }
            return Divide(difference, new Pair(alpha));
        }

        /// <summary>Stores a leading double and its nonoverlapping rounding correction.</summary>
        private readonly struct Pair
        {
            internal readonly double High;
            internal readonly double Low;

            internal Pair(double high, double low = 0d)
            {
                High = high;
                Low = low;
            }
        }

        /// <summary>Renormalizes two components using an error-free finite sum.</summary>
        private static Pair Normalize(double high, double low)
        {
            double sum = high + low;
            if (double.IsInfinity(sum) || double.IsNaN(sum)) return new Pair(sum);
            double part = sum - high;
            return new Pair(sum, (high - (sum - part)) + (low - part));
        }

        /// <summary>Adds compensated pairs, retaining the leading addition error.</summary>
        private static Pair Add(Pair first, Pair second)
        {
            double sum = first.High + second.High;
            if (double.IsInfinity(sum) || double.IsNaN(sum)) return new Pair(sum);
            double part = sum - first.High;
            double error = (first.High - (sum - part)) + (second.High - part);
            return Normalize(sum, error + first.Low + second.Low);
        }

        /// <summary>Subtracts compensated pairs.</summary>
        private static Pair Subtract(Pair first, Pair second)
        {
            return Add(first, new Pair(-second.High, -second.Low));
        }

        /// <summary>Returns the absolute value of a compensated pair.</summary>
        private static Pair Abs(Pair value)
        {
            return value.High < 0d ? new Pair(-value.High, -value.Low) : value;
        }

        /// <summary>
        /// Multiplies compensated pairs using Dekker's product correction. Splitting the
        /// significand with a bit mask avoids overflowing a large splitter multiplication.
        /// </summary>
        private static Pair Multiply(Pair first, Pair second)
        {
            double product = first.High * second.High;
            if (double.IsInfinity(product) || double.IsNaN(product)) return new Pair(product);
            const long mask = ~((1L << 27) - 1L);
            double firstHigh = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(first.High) & mask);
            double secondHigh = BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(second.High) & mask);
            double firstLow = first.High - firstHigh;
            double secondLow = second.High - secondHigh;
            double error = ((firstHigh * secondHigh - product) + firstHigh * secondLow + firstLow * secondHigh)
                + firstLow * secondLow;
            error += first.High * second.Low + first.Low * second.High + first.Low * second.Low;
            return Normalize(product, error);
        }

        /// <summary>Divides compensated pairs using two residual corrections.</summary>
        private static Pair Divide(Pair numerator, Pair denominator)
        {
            double quotient = numerator.High / denominator.High;
            if (double.IsInfinity(quotient) || double.IsNaN(quotient)) return new Pair(quotient);
            Pair result = new Pair(quotient);
            Pair residual = Subtract(numerator, Multiply(denominator, result));
            double correction = residual.High / denominator.High;
            result = Add(result, new Pair(correction));
            residual = Subtract(residual, Multiply(denominator, new Pair(correction)));
            return Add(result, new Pair(residual.High / denominator.High));
        }

        /// <summary>Scales a pair by a power of two without overflowing a scaling factor.</summary>
        private static Pair Scale(Pair value, int exponent)
        {
            while (exponent > 512)
            {
                value = new Pair(value.High * 1.3407807929942597E154d, value.Low * 1.3407807929942597E154d);
                exponent -= 512;
            }
            while (exponent < -512)
            {
                value = new Pair(value.High * 7.458340731200207E-155d, value.Low * 7.458340731200207E-155d);
                exponent += 512;
            }
            double factor = Math.Pow(2d, exponent);
            return new Pair(value.High * factor, value.Low * factor);
        }

        /// <summary>Evaluates a positive pair's logarithm by binary reduction and an atanh series.</summary>
        private static Pair Log(Pair value)
        {
            if (!(value.High > 0d) || double.IsInfinity(value.High)) return new Pair(Math.Log(value.High));
            double leading = value.High;
            int correction = 0;
            if (leading < 2.2250738585072014E-308d)
            {
                leading *= 18014398509481984d;
                correction = -54;
            }
            int exponent = (int)((BitConverter.DoubleToInt64Bits(leading) >> 52) & 0x7ffL) - 1023 + correction;
            Pair reduced = Scale(value, -exponent);
            if (reduced.High > 1.4142135623730951d)
            {
                reduced = Scale(reduced, -1);
                exponent++;
            }
            Pair ratio = Divide(Subtract(reduced, new Pair(1d)), Add(reduced, new Pair(1d)));
            Pair square = Multiply(ratio, ratio);
            Pair term = ratio;
            Pair sum = ratio;
            for (int index = 1; index <= 24; index++)
            {
                term = Multiply(term, square);
                sum = Add(sum, Divide(term, new Pair(2d * index + 1d)));
            }
            return Add(Multiply(sum, new Pair(2d)), Multiply(LogTwo, new Pair(exponent)));
        }

        /// <summary>Evaluates log(1+v)/v continuously near zero without dividing by a tiny v.</summary>
        private static Pair LogOnePlusRelative(Pair value)
        {
            Pair sum = new Pair(1d);
            Pair term = new Pair(1d);
            Pair negative = new Pair(-value.High, -value.Low);
            for (int index = 1; index <= 40; index++)
            {
                term = Multiply(term, negative);
                sum = Add(sum, Divide(term, new Pair(index + 1d)));
            }
            return sum;
        }
    }
}
