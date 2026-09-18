using System;
using System.Collections.Generic;
using System.Linq;
using System.Xml.Linq;
using BigInteger = System.Numerics.BigInteger;

namespace Numerics.Distributions
{
    /// <summary>Numerical primitives shared by the reviewed univariate distributions.</summary>
    internal static partial class DistributionNumerics
    {

        /// <summary>Uses the established prior envelope whenever its initialization and bounds are usable.</summary>
        /// <param name="legacy">The original family-specific initialization and bounds.</param>
        /// <param name="fallback">The exceptional-input hardened initialization and bounds.</param>
        /// <returns>The legacy constraints if finite and within ordered bounds; otherwise the hardened constraints.</returns>
        /// <exception cref="ArgumentException">A selected constraint delegate rejects its input.</exception>
        /// <exception cref="InvalidOperationException">A selected constraint delegate cannot produce usable constraints.</exception>
        /// <exception cref="ArithmeticException">A selected constraint delegate encounters an arithmetic failure.</exception>
        /// <remarks>The fallback must not narrow previously valid prior envelopes or alter their rounding.</remarks>
        internal static Tuple<double[], double[], double[]> PreferLegacyConstraints(
            Func<Tuple<double[], double[], double[]>> legacy, Func<Tuple<double[], double[], double[]>> fallback)
        {
            Exception? legacyFailure = null;
            try
            {
                var constraints = legacy();
                bool usable = true;
                for (int i = 0; i < constraints.Item1.Length; i++)
                {
                    double initial = constraints.Item1[i], lower = constraints.Item2[i], upper = constraints.Item3[i];
                    if (!Tools.IsFinite(initial) || !Tools.IsFinite(lower) || !Tools.IsFinite(upper) || !(lower < upper && lower <= initial && initial <= upper))
                    {
                        usable = false;
                        break;
                    }
                }
                if (usable) return constraints;
            }
            catch (Exception exception) when (exception is ArgumentException || exception is InvalidOperationException || exception is ArithmeticException)
            {
                // Preserve initialization context if the exceptional-input path also fails.
                legacyFailure = exception;
            }
            try
            {
                return fallback();
            }
            catch (Exception exception) when (legacyFailure != null &&
                (exception is ArgumentException || exception is InvalidOperationException || exception is ArithmeticException))
            {
                exception.Data["LegacyParameterConstraintsFailure"] = legacyFailure;
                throw;
            }
        }

        /// <summary>A cache key including nested distribution settings omitted from flattened parameter vectors.</summary>
        /// <param name="distribution">The distribution whose complete mutable configuration is serialized into the key.</param>
        /// <returns>An invariant, deterministic representation of the distribution and any nested component configuration.</returns>
        internal static string ConfigurationState(UnivariateDistributionBase distribution)
        {
            // Composite serialization flattens scalar parameters and need not support an empirical
            // child. Cache identity therefore records configuration independently of that contract.
            if (distribution is Mixture mixture)
                return FormattableString.Invariant($"Mixture|{mixture.IsZeroInflated}|{mixture.ZeroWeight:R}|{mixture.XTransform}|{mixture.ProbabilityTransform}|")
                    + string.Join("|", mixture.Weights.Select(w => w.ToString("R", System.Globalization.CultureInfo.InvariantCulture)))
                    + "[" + string.Join("][", mixture.Distributions.Select(ConfigurationState)) + "]";
            if (distribution is CompetingRisks competing)
                return FormattableString.Invariant($"Competing|{competing.MinimumOfRandomVariables}|{competing.Dependency}|{competing.PRNGSeed}|{competing.XTransform}|{competing.ProbabilityTransform}|")
                    + (competing.CorrelationMatrix is null ? "null" : competing.CorrelationMatrix.GetLength(0) + "x" + competing.CorrelationMatrix.GetLength(1)
                        + ":" + string.Join("|", competing.CorrelationMatrix.Cast<double>().Select(v => v.ToString("R", System.Globalization.CultureInfo.InvariantCulture))))
                    + "[" + string.Join("][", competing.Distributions.Select(ConfigurationState)) + "]";
            return distribution.ToXElement().ToString(SaveOptions.DisableFormatting);
        }

        /// <summary>Forms an affine standardized value without overflowing a finite difference unnecessarily.</summary>
        /// <param name="x">The value in physical coordinates.</param>
        /// <param name="location">The location to subtract.</param>
        /// <param name="scale">The scale by which to divide the centered value.</param>
        /// <returns><c>(x-location)/scale</c>, using an algebraically equivalent form when subtraction of finite operands overflows.</returns>
        internal static double Standardize(double x, double location, double scale)
        {
            double difference = x - location;
            return double.IsInfinity(difference) && Tools.IsFinite(x) && Tools.IsFinite(location)
                ? x / scale - location / scale : difference / scale;
        }

        /// <summary>Computes log(1-exp(a)) for a nonpositive log probability, including its limits.</summary>
        /// <param name="a">The logarithm of a probability in the interval from negative infinity through zero.</param>
        /// <returns><c>log(1-exp(a))</c>, or not-a-number when <paramref name="a"/> is positive or not-a-number.</returns>
        internal static double Log1mExp(double a)
        {
            if (a > 0 || double.IsNaN(a)) return double.NaN;
            return a < -0.69314718055994530942 ? Tools.Log1p(-Math.Exp(a)) : Math.Log(-Tools.Expm1(a));
        }

        /// <summary>Computes log(exp(a)-exp(b)), with equal arguments representing zero mass.</summary>
        /// <param name="a">The logarithm of the minuend.</param>
        /// <param name="b">The logarithm of the subtrahend, which must not exceed <paramref name="a"/>.</param>
        /// <returns>The logarithmic difference, negative infinity for equal arguments, or not-a-number for invalid ordering or input.</returns>
        internal static double LogDifference(double a, double b)
        {
            if (double.IsNaN(a) || double.IsNaN(b) || b > a) return double.NaN;
            if (a == b) return double.NegativeInfinity;
            return a + Log1mExp(b - a);
        }

        /// <summary>Adds two nonnegative quantities represented by their logarithms.</summary>
        /// <param name="a">The logarithm of the first nonnegative quantity.</param>
        /// <param name="b">The logarithm of the second nonnegative quantity.</param>
        /// <returns>The logarithm of the sum, including the natural infinity and not-a-number limits.</returns>
        internal static double LogSum(double a, double b)
        {
            if (double.IsNaN(a) || double.IsNaN(b)) return double.NaN;
            if (double.IsPositiveInfinity(a) || double.IsPositiveInfinity(b)) return double.PositiveInfinity;
            if (double.IsNegativeInfinity(a)) return b;
            if (double.IsNegativeInfinity(b)) return a;
            double larger = Math.Max(a, b);
            return larger + Tools.Log1p(Math.Exp(Math.Min(a, b) - larger));
        }

        /// <summary>The normal log CDF, preserving the logarithm after the tail itself underflows.</summary>
        /// <param name="z">The standard Normal variate.</param>
        /// <returns>The natural logarithm of the standard Normal cumulative probability, including endpoint and not-a-number limits.</returns>
        /// <remarks>The far tail uses the convergent Laplace continued fraction for the Mills ratio.</remarks>
        internal static double NormalLogCDF(double z)
        {
            if (double.IsNaN(z)) return double.NaN;
            if (z > 0) return Log1mExp(NormalLogCDF(-z));
            if (z >= -10) return Math.Log(Normal.StandardCDF(z));
            if (double.IsNegativeInfinity(z)) return double.NegativeInfinity;
            double x = -z;
            // Q(x)/phi(x) = 1/(x+1/(x+2/(x+3/(...)))). At x>=10,
            // 64 backward levels are well beyond binary64 convergence.
            double fraction = 0;
            for (int i = 64; i >= 1; i--) fraction = i / (x + fraction);
            return -(0.5 * x) * x - Tools.LogSqrt2PI - Math.Log(x + fraction);
        }

        /// <summary>The normal log survival function, evaluated directly through reflection.</summary>
        /// <param name="z">The standard Normal variate.</param>
        /// <returns>The natural logarithm of the standard Normal probability above <paramref name="z"/>.</returns>
        internal static double NormalLogSurvival(double z) => NormalLogCDF(-z);

        /// <summary>The exponential divided difference expm1(x)/x, including x=0.</summary>
        /// <param name="x">The divided-exponential argument.</param>
        /// <returns><c>expm1(x)/x</c>, with the continuous value one at zero.</returns>
        internal static double Exprel(double x) => x == 0 ? 1 : Tools.Expm1(x) / x;

        /// <summary>The derivative of expm1(x)/x, with a convergent series at zero.</summary>
        /// <param name="x">The divided-exponential argument.</param>
        /// <returns>The derivative of <c>expm1(x)/x</c>, including its continuous value one-half at zero.</returns>
        internal static double ExprelDerivative(double x)
        {
            if (Math.Abs(x) >= 0.1) return ((x - 1) * Math.Exp(x) + 1) / x / x;
            double sum = 0.5, term = 0.5;
            for (int n = 1; n < 20; n++)
            {
                term *= x * (n + 1.0) / n / (n + 2.0);
                sum += term;
                if (Math.Abs(term) <= Math.Abs(sum) * 1E-17) break;
            }
            return sum;
        }

        /// <summary>Rejects nonfinite and endpoint probabilities for quantile uncertainty calculations.</summary>
        /// <param name="probability">The probability to validate.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="probability"/> is not finite and strictly between zero and one.</exception>
        internal static void ValidateProbability(double probability)
        {
            if (!(probability > 0 && probability < 1))
                throw new ArgumentOutOfRangeException(nameof(probability), "Quantile uncertainty requires a finite probability strictly between zero and one.");
        }

        /// <summary>Requires a positive sample size for asymptotic uncertainty.</summary>
        /// <param name="sampleSize">The number of independent observations.</param>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="sampleSize"/> is not positive.</exception>
        internal static void ValidateSampleSize(int sampleSize)
        {
            if (sampleSize <= 0) throw new ArgumentOutOfRangeException(nameof(sampleSize), "Sample size must be positive.");
        }

        /// <summary>Checks sample size and finite interior probabilities for interval approximations requiring n-1.</summary>
        /// <param name="sampleSize">The available number of observations.</param>
        /// <param name="quantiles">The quantile probabilities to validate.</param>
        /// <param name="percentiles">The confidence probabilities to validate.</param>
        /// <param name="minimumSampleSize">The smallest sample size accepted by the calling interval method.</param>
        /// <exception cref="ArgumentNullException"><paramref name="quantiles"/> or <paramref name="percentiles"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is too small, a probability list is empty, or a probability is not finite and strictly interior.</exception>
        internal static void ValidateConfidenceInputs(int sampleSize, IList<double> quantiles, IList<double> percentiles, int minimumSampleSize = 1)
        {
            if (sampleSize < minimumSampleSize) throw new ArgumentOutOfRangeException(nameof(sampleSize), "Insufficient observations for this confidence interval method.");
            if (quantiles == null) throw new ArgumentNullException(nameof(quantiles));
            if (percentiles == null) throw new ArgumentNullException(nameof(percentiles));
            if (quantiles.Count == 0 || percentiles.Count == 0) throw new ArgumentOutOfRangeException(nameof(quantiles), "Probability lists must not be empty.");
            foreach (double p in quantiles) ValidateProbability(p);
            foreach (double p in percentiles) ValidateProbability(p);
        }

        /// <summary>Checks that an initialization sample is finite, sufficiently long, and nonconstant.</summary>
        /// <param name="sample">The observations to validate.</param>
        /// <param name="minimumCount">The minimum accepted number of observations.</param>
        /// <param name="positive"><see langword="true"/> to require every observation to be strictly positive; otherwise, <see langword="false"/>.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is too short, contains an invalid observation, or is constant.</exception>
        internal static void ValidateSample(IList<double> sample, int minimumCount = 2, bool positive = false)
        {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < minimumCount) throw new ArgumentOutOfRangeException(nameof(sample), "Insufficient observations to initialize the distribution.");
            bool distinct = false;
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]) || (positive && sample[i] <= 0))
                    throw new ArgumentOutOfRangeException(nameof(sample), positive ? "Observations must be finite and strictly positive." : "Observations must be finite.");
                distinct |= sample[i] != sample[0];
            }
            if (!distinct) throw new ArgumentOutOfRangeException(nameof(sample), "A constant sample cannot initialize a positive scale.");
        }

        /// <summary>Assembles quantile gradients by observation row and obtains a determinant without artificial pivots.</summary>
        /// <param name="distribution">The distribution that supplies gradients in public parameter coordinates.</param>
        /// <param name="probabilities">One finite interior probability per matrix row and public parameter.</param>
        /// <param name="determinant">The signed determinant, including zero for an exactly singular Jacobian.</param>
        /// <returns>The square quantile-gradient matrix.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> or <paramref name="probabilities"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The probability count, an individual probability, or a returned gradient dimension is invalid.</exception>
        /// <exception cref="InvalidOperationException">A quantile derivative is nonfinite or the determinant cannot be evaluated from finite entries.</exception>
        internal static double[,] QuantileJacobian(IStandardError distribution, IList<double> probabilities, out double determinant)
        {
            var matrix = QuantileGradientMatrix(distribution, probabilities);
            double log = LogAbsDeterminant(matrix, out int sign);
            determinant = sign == 0 ? 0 : sign * Math.Exp(log);
            return matrix;
        }

        /// <summary>Builds the square Jacobian with one quantile per row in the public parameter coordinates.</summary>
        /// <param name="distribution">The distribution that supplies quantile gradients.</param>
        /// <param name="probabilities">One finite interior probability per matrix row and public parameter.</param>
        /// <returns>The square matrix whose rows are quantile gradients.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="distribution"/> or <paramref name="probabilities"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The probability count, an individual probability, or a returned gradient dimension is invalid.</exception>
        /// <exception cref="InvalidOperationException">A returned quantile derivative is nonfinite.</exception>
        internal static double[,] QuantileGradientMatrix(IStandardError distribution, IList<double> probabilities)
        {
            if (distribution == null) throw new ArgumentNullException(nameof(distribution));
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            int count = probabilities.Count;
            if (count == 0) throw new ArgumentOutOfRangeException(nameof(probabilities));
            if (distribution is IUnivariateDistribution univariate && count != univariate.NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(probabilities), "Provide one probability per public distribution parameter.");
            var matrix = new double[count, count];
            for (int i = 0; i < count; i++)
            {
                ValidateProbability(probabilities[i]);
                double[] gradient = distribution.QuantileGradient(probabilities[i]);
                if (gradient.Length != count) throw new ArgumentOutOfRangeException(nameof(probabilities), "The quantile Jacobian must be square.");
                for (int j = 0; j < count; j++)
                {
                    if (!Tools.IsFinite(gradient[j])) throw new InvalidOperationException("The quantile Jacobian contains a nonfinite derivative.");
                    matrix[i, j] = gradient[j];
                }
            }
            return matrix;
        }

        /// <summary>Log absolute determinant by row/column equilibration and partial pivoting.</summary>
        /// <param name="matrix">The finite square matrix to factor.</param>
        /// <param name="sign">The determinant sign, or zero when the matrix is exactly singular.</param>
        /// <returns>The natural logarithm of the absolute determinant, or negative infinity for exact singularity.</returns>
        /// <exception cref="ArgumentException"><paramref name="matrix"/> is not square.</exception>
        /// <exception cref="InvalidOperationException"><paramref name="matrix"/> contains a nonfinite entry.</exception>
        /// <remarks>Zero pivots retain exact singularity; no jitter or artificial tiny pivots are inserted.</remarks>
        internal static double LogAbsDeterminant(double[,] matrix, out int sign)
        {
            int n = matrix.GetLength(0);
            if (matrix.GetLength(1) != n) throw new ArgumentException("The matrix must be square.", nameof(matrix));
            var a = (double[,])matrix.Clone();
            double log = 0;
            sign = 1;
            for (int i = 0; i < n; i++)
            {
                double scale = 0;
                for (int j = 0; j < n; j++) scale = Math.Max(scale, Math.Abs(a[i, j]));
                if (!Tools.IsFinite(scale)) throw new InvalidOperationException("The determinant requires finite matrix entries.");
                if (scale == 0) { sign = 0; return double.NegativeInfinity; }
                log += Math.Log(scale);
                for (int j = 0; j < n; j++)
                {
                    a[i, j] /= scale;
                    if (a[i, j] == 0 && matrix[i, j] != 0) return ExactDyadicLogDeterminant(matrix, out sign);
                }
            }
            for (int j = 0; j < n; j++)
            {
                double scale = 0;
                for (int i = 0; i < n; i++) scale = Math.Max(scale, Math.Abs(a[i, j]));
                if (scale == 0) return ExactDyadicLogDeterminant(matrix, out sign);
                log += Math.Log(scale);
                for (int i = 0; i < n; i++) a[i, j] /= scale;
            }
            for (int j = 0; j < n; j++)
            {
                int pivot = j;
                for (int i = j + 1; i < n; i++) if (Math.Abs(a[i, j]) > Math.Abs(a[pivot, j])) pivot = i;
                // An exactly dependent row can leave a rounded residual after equilibration.
                // Resolve a small pivot using the original binary64 values as exact dyadic integers.
                if (Math.Abs(a[pivot, j]) < 1E-10) return ExactDyadicLogDeterminant(matrix, out sign);
                if (pivot != j)
                {
                    for (int k = j; k < n; k++) (a[j, k], a[pivot, k]) = (a[pivot, k], a[j, k]);
                    sign = -sign;
                }
                double diagonal = a[j, j];
                sign *= Math.Sign(diagonal);
                log += Math.Log(Math.Abs(diagonal));
                for (int i = j + 1; i < n; i++)
                {
                    double ratio = a[i, j] / diagonal;
                    for (int k = j + 1; k < n; k++) a[i, k] -= ratio * a[j, k];
                }
            }
            return log;
        }

        /// <summary>Fraction-free integer elimination for a determinant whose floating-point pivot is unresolved.</summary>
        /// <param name="matrix">The finite square binary64 matrix to evaluate exactly after dyadic scaling.</param>
        /// <param name="sign">The exact determinant sign, or zero when the matrix is singular.</param>
        /// <returns>The natural logarithm of the absolute determinant, or negative infinity for exact singularity.</returns>
        /// <remarks>Each finite binary64 row is scaled by an exact power of two to integers.
        /// Bareiss elimination then distinguishes exact dependence from a merely small determinant.
        /// The threshold selecting this path does not classify a matrix as singular.</remarks>
        private static double ExactDyadicLogDeterminant(double[,] matrix, out int sign)
        {
            int n = matrix.GetLength(0), binaryExponent = 0;
            var integers = new BigInteger[n, n];
            for (int i = 0; i < n; i++)
            {
                int minimumExponent = int.MaxValue;
                for (int j = 0; j < n; j++)
                {
                    if (matrix[i, j] == 0) continue;
                    long bits = BitConverter.DoubleToInt64Bits(matrix[i, j]);
                    int exponent = (int)((bits >> 52) & 0x7ff);
                    minimumExponent = Math.Min(minimumExponent, exponent == 0 ? -1074 : exponent - 1075);
                }
                if (minimumExponent == int.MaxValue) { sign = 0; return double.NegativeInfinity; }
                binaryExponent += minimumExponent;
                for (int j = 0; j < n; j++)
                {
                    long bits = BitConverter.DoubleToInt64Bits(matrix[i, j]);
                    int exponent = (int)((bits >> 52) & 0x7ff);
                    long mantissa = bits & 0xfffffffffffffL;
                    if (exponent != 0) mantissa |= 0x10000000000000L;
                    if (mantissa == 0) continue;
                    int power = exponent == 0 ? -1074 : exponent - 1075;
                    integers[i, j] = new BigInteger(bits < 0 ? -mantissa : mantissa) << (power - minimumExponent);
                }
            }
            sign = 1;
            BigInteger previous = BigInteger.One;
            for (int k = 0; k < n - 1; k++)
            {
                int pivot = k;
                while (pivot < n && integers[pivot, k].IsZero) pivot++;
                if (pivot == n) { sign = 0; return double.NegativeInfinity; }
                if (pivot != k)
                {
                    for (int j = k; j < n; j++) (integers[k, j], integers[pivot, j]) = (integers[pivot, j], integers[k, j]);
                    sign = -sign;
                }
                BigInteger diagonal = integers[k, k];
                for (int i = k + 1; i < n; i++)
                {
                    for (int j = k + 1; j < n; j++)
                        integers[i, j] = (diagonal * integers[i, j] - integers[i, k] * integers[k, j]) / previous;
                    integers[i, k] = BigInteger.Zero;
                }
                previous = diagonal;
            }
            BigInteger determinant = integers[n - 1, n - 1];
            sign *= determinant.Sign;
            return sign == 0 ? double.NegativeInfinity : BigInteger.Log(BigInteger.Abs(determinant)) + binaryExponent * Math.Log(2);
        }

        /// <summary>Whether a reviewed distribution has no probability atoms.</summary>
        /// <param name="distribution">The distribution or composite distribution to classify.</param>
        /// <returns><see langword="true"/> when the distribution and every nested component are recognized as continuous and the mixture has no hurdle atom; otherwise, <see langword="false"/>.</returns>
        private static bool IsContinuous(UnivariateDistributionBase distribution)
        {
            if (distribution is CompetingRisks competing) return competing.Distributions.All(IsContinuous);
            if (distribution is Mixture mixture) return !mixture.IsZeroInflated && mixture.Distributions.All(IsContinuous);
            return distribution is Normal || distribution is Logistic || distribution is LnNormal || distribution is LogNormal
                || distribution is PearsonTypeIII || distribution is LogPearsonTypeIII || distribution is Exponential
                || distribution is GammaDistribution || distribution is Weibull || distribution is Gumbel
                || distribution is GeneralizedExtremeValue || distribution is GeneralizedPareto
                || distribution is GeneralizedNormal || distribution is GeneralizedLogistic || distribution is KappaFour
                || distribution is Uniform;
        }

        /// <summary>Resolves a positive continuous interval when both pairs of log tails round to identical values.</summary>
        /// <param name="distribution">The continuous distribution whose interval probability is required.</param>
        /// <param name="lower">The open lower interval endpoint.</param>
        /// <param name="upper">The closed upper interval endpoint.</param>
        /// <returns>The logarithm of the integrated interval probability, or negative infinity when the distribution is not recognized as continuous or the clipped interval has no finite positive width.</returns>
        /// <remarks>Eight-point Gauss-Legendre integration is used only after tail subtraction collapses.
        /// The explicitly identified continuous families avoid treating an atom as a density contribution.</remarks>
        internal static double CollapsedContinuousLogInterval(UnivariateDistributionBase distribution, double lower, double upper)
        {
            if (!IsContinuous(distribution)) return double.NegativeInfinity;
            lower = Math.Max(lower, distribution.Minimum);
            upper = Math.Min(upper, distribution.Maximum);
            double width = upper - lower;
            if (!(width > 0) || !Tools.IsFinite(width)) return double.NegativeInfinity;
            double[] nodes = { .019855071751231884, .10166676129318663, .23723379504183551, .4082826787521751,
                .5917173212478249, .7627662049581645, .8983332387068134, .9801449282487681 };
            double[] weights = { .05061426814518813, .11119051722668724, .15685332293894365, .181341891689181,
                .181341891689181, .15685332293894365, .11119051722668724, .05061426814518813 };
            double sum = double.NegativeInfinity;
            for (int i = 0; i < nodes.Length; i++)
                sum = LogSum(sum, Math.Log(weights[i]) + distribution.LogPDF(lower + width * nodes[i]));
            return Math.Log(width) + sum;
        }
    }
}
