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

        /// <summary>
        /// Returns the initial values and bounds of a three-parameter (min, mode, max) distribution fitted to a sample.
        /// </summary>
        /// <param name="sample">The observations.</param>
        /// <param name="momentMode">The method-of-moments mode of the family, given the sample mean, minimum and maximum.</param>
        /// <returns>The initial values and the lower and upper bounds of the min, mode and max.</returns>
        /// <remarks>
        /// The established envelope bounds the min below by -10^(ceil(log10|min|) + 1) and the max above by
        /// 10^(ceil(log10|max|) + 1), and starts the mode at the moment estimate. When that envelope does not contain
        /// its initial values, for example when the sample minimum or maximum is exactly zero, the bound magnitude is
        /// S = max(|min|, |max|, max - min), or 1 when S is zero, the bounds are -10^(ceil(log10 S) + 1) and
        /// 10^(ceil(log10 S) + 1), which always contain the sample range, and the initial mode is clamped into the
        /// sample range. The min and max start one machine epsilon outside the sample range.
        /// </remarks>
        internal static Tuple<double[], double[], double[]> MinModeMaxConstraints(IList<double> sample, Func<double, double, double, double> momentMode)
        {
            return PreferLegacyConstraints(
                () => LegacyMinModeMaxConstraints(sample, momentMode), () => RobustMinModeMaxConstraints(sample, momentMode));
        }

        /// <summary>Returns the established (min, mode, max) initialization and prior envelope.</summary>
        /// <param name="sample">The observations.</param>
        /// <param name="momentMode">The method-of-moments mode, given the sample mean, minimum and maximum.</param>
        /// <returns>The legacy initial values and lower and upper bounds.</returns>
        private static Tuple<double[], double[], double[]> LegacyMinModeMaxConstraints(IList<double> sample, Func<double, double, double, double> momentMode)
        {
            var initialVals = new double[3];
            var lowerVals = new double[3];
            var upperVals = new double[3];
            // Estimate initial values using the method of moments (a.k.a product moments).
            var min = Tools.Min(sample);
            var max = Tools.Max(sample);
            var mean = Tools.Mean(sample);
            var mode = momentMode(mean, min, max);
            initialVals[0] = min - Tools.DoubleMachineEpsilon;
            initialVals[1] = mode;
            initialVals[2] = max + Tools.DoubleMachineEpsilon;
            // Get bounds of min
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(min)) + 1d));
            upperVals[0] = min;
            // Get bounds of mode
            lowerVals[1] = min;
            upperVals[1] = max;
            // Get bounds of max
            lowerVals[2] = max;
            upperVals[2] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(max)) + 1d));

            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Returns (min, mode, max) initial values inside bounds whose magnitude follows the sample scale.</summary>
        /// <param name="sample">The observations.</param>
        /// <param name="momentMode">The method-of-moments mode, given the sample mean, minimum and maximum.</param>
        /// <returns>Initial values inside bounds that contain the sample range.</returns>
        private static Tuple<double[], double[], double[]> RobustMinModeMaxConstraints(IList<double> sample, Func<double, double, double, double> momentMode)
        {
            var initialVals = new double[3];
            var lowerVals = new double[3];
            var upperVals = new double[3];
            var min = Tools.Min(sample);
            var max = Tools.Max(sample);
            var mean = Tools.Mean(sample);
            var mode = Tools.Clamp(momentMode(mean, min, max), min, max);
            double scale = Math.Max(Math.Max(Math.Abs(min), Math.Abs(max)), max - min);
            if (!(scale > 0d)) scale = 1d;
            double bound = Math.Pow(10d, Math.Ceiling(Math.Log10(scale) + 1d));
            initialVals[0] = min - Tools.DoubleMachineEpsilon;
            initialVals[1] = mode;
            initialVals[2] = max + Tools.DoubleMachineEpsilon;
            lowerVals[0] = -bound;
            upperVals[0] = min;
            lowerVals[1] = min;
            upperVals[1] = max;
            lowerVals[2] = max;
            upperVals[2] = bound;
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
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

        /// <summary>Writes the evaluation-affecting settings that lie outside a distribution's flattened parameter vector.</summary>
        /// <param name="distribution">The distribution whose settings are serialized.</param>
        /// <param name="element">The element receiving the optional attributes, or <see langword="null"/> to report only whether any would be written.</param>
        /// <returns><see langword="true"/> when at least one attribute applies; otherwise, <see langword="false"/>.</returns>
        /// <remarks>
        /// Each attribute is written only when its value differs from the family default, so the serialized form of a
        /// default-configured distribution is unchanged. <see cref="LogNormal"/> and <see cref="LogPearsonTypeIII"/>
        /// record a logarithm base other than 10. <see cref="LnNormal"/> records its natural-log coordinates
        /// <c>Mu</c> and <c>Sigma</c> when its stored physical-moment surface is inactive, because the flattened
        /// mean and standard deviation then do not determine those coordinates exactly.
        /// </remarks>
        internal static bool WriteConfigurationAttributes(UnivariateDistributionBase distribution, XElement? element)
        {
            var culture = System.Globalization.CultureInfo.InvariantCulture;
            switch (distribution)
            {
                case LogNormal logNormal when logNormal.Base != 10d:
                    element?.SetAttributeValue(nameof(LogNormal.Base), logNormal.Base.ToString("G17", culture));
                    return true;
                case LogPearsonTypeIII logPearson when logPearson.Base != 10d:
                    element?.SetAttributeValue(nameof(LogPearsonTypeIII.Base), logPearson.Base.ToString("G17", culture));
                    return true;
                case LnNormal lnNormal when !lnNormal.PhysicalMomentModeForSnapshot:
                    element?.SetAttributeValue(nameof(LnNormal.Mu), lnNormal.Mu.ToString("G17", culture));
                    element?.SetAttributeValue(nameof(LnNormal.Sigma), lnNormal.Sigma.ToString("G17", culture));
                    return true;
                default:
                    return false;
            }
        }

        /// <summary>Restores the optional settings written by <see cref="WriteConfigurationAttributes"/>.</summary>
        /// <param name="element">The serialized distribution element.</param>
        /// <param name="distribution">The distribution, already configured from its flattened parameters.</param>
        /// <exception cref="ArgumentException">An optional attribute is present but malformed or out of range.</exception>
        /// <remarks>Absent attributes leave the family defaults in place, so elements written without them restore as before.</remarks>
        internal static void ReadConfigurationAttributes(XElement element, UnivariateDistributionBase distribution)
        {
            switch (distribution)
            {
                case LogNormal logNormal:
                    if (TryReadAttribute(element, nameof(LogNormal.Base), out double logNormalBase))
                        logNormal.Base = logNormalBase > 1d ? logNormalBase : throw InvalidBase(element);
                    break;
                case LogPearsonTypeIII logPearson:
                    if (TryReadAttribute(element, nameof(LogPearsonTypeIII.Base), out double logPearsonBase))
                        logPearson.Base = logPearsonBase > 1d ? logPearsonBase : throw InvalidBase(element);
                    break;
                case LnNormal lnNormal:
                    bool hasMu = TryReadAttribute(element, nameof(LnNormal.Mu), out double mu);
                    bool hasSigma = TryReadAttribute(element, nameof(LnNormal.Sigma), out double sigma);
                    if (hasMu != hasSigma)
                        throw new ArgumentException("The serialized natural-log coordinates must include both Mu and Sigma.", nameof(element));
                    if (hasMu)
                    {
                        lnNormal.Mu = mu;
                        lnNormal.Sigma = sigma;
                    }
                    break;
            }
        }

        /// <summary>Parses an optional finite floating-point attribute.</summary>
        /// <param name="element">The element holding the attribute.</param>
        /// <param name="name">The attribute name.</param>
        /// <param name="value">The parsed value, or zero when the attribute is absent.</param>
        /// <returns><see langword="true"/> when the attribute is present; otherwise, <see langword="false"/>.</returns>
        /// <exception cref="ArgumentException">The attribute is present but is not a finite number.</exception>
        private static bool TryReadAttribute(XElement element, string name, out double value)
        {
            value = 0d;
            var attribute = element.Attribute(name);
            if (attribute == null) return false;
            if (!double.TryParse(attribute.Value, System.Globalization.NumberStyles.Any, System.Globalization.CultureInfo.InvariantCulture, out value)
                || !Tools.IsFinite(value))
                throw new ArgumentException("The serialized distribution setting '" + name + "' is invalid.", nameof(element));
            return true;
        }

        /// <summary>Creates the exception for a serialized logarithm base that is not greater than one.</summary>
        /// <param name="element">The serialized distribution element.</param>
        /// <returns>The argument exception to throw.</returns>
        private static ArgumentException InvalidBase(XElement element)
        {
            return new ArgumentException("The serialized logarithm base of '" + element.Attribute("Type")?.Value + "' must be greater than one.", nameof(element));
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

        /// <summary>Computes log(exp(a)-exp(b)), with equal arguments representing zero mass.</summary>
        /// <param name="a">The logarithm of the minuend.</param>
        /// <param name="b">The logarithm of the subtrahend, which must not exceed <paramref name="a"/>.</param>
        /// <returns>The logarithmic difference, negative infinity for equal arguments, or not-a-number for invalid ordering or input.</returns>
        internal static double LogDifference(double a, double b)
        {
            if (double.IsNaN(a) || double.IsNaN(b) || b > a) return double.NaN;
            if (a == b) return double.NegativeInfinity;
            return a + Tools.Log1mExp(b - a);
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
            if (z > 0) return Tools.Log1mExp(NormalLogCDF(-z));
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

        /// <summary>The element count from which <see cref="ForEachElementRange"/> runs its ranges in parallel.</summary>
        private const int ParallelElementThreshold = 4096;

        /// <summary>The number of elements in each parallel range of <see cref="ForEachElementRange"/>.</summary>
        private const int ParallelElementRange = 2048;

        /// <summary>
        /// Runs an element-wise kernel over [0, count), in parallel ranges when the count is large.
        /// </summary>
        /// <param name="count">The number of elements.</param>
        /// <param name="body">The kernel, called with the half-open range [start, end) it must fill.</param>
        /// <remarks>
        /// Each element is computed independently by the same kernel whatever range it falls in, so the results do
        /// not depend on the scheduling or the processor count. Counts below 4096 run serially in one range; larger
        /// counts are split into ranges of 2048 elements scheduled by <see cref="Parallel.For(int, int, Action{int})"/>.
        /// </remarks>
        internal static void ForEachElementRange(int count, Action<int, int> body)
        {
            if (count < ParallelElementThreshold)
            {
                body(0, count);
                return;
            }
            int ranges = (count + ParallelElementRange - 1) / ParallelElementRange;
            Parallel.For(0, ranges, r =>
            {
                int start = r * ParallelElementRange;
                body(start, Math.Min(count, start + ParallelElementRange));
            });
        }
    }
}
