using Numerics.Data.Statistics;
using System;
using System.Collections.Generic; 
using System.IO;
using System.IO.Compression;
using System.Runtime.CompilerServices;
using System.Threading;

namespace Numerics
{
    /// <summary>
    /// A class of public utility functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public static class Tools
    {

        /// <summary>
        /// Value used for comparing double precision equivalence: 1.11022302462516E-16, approximately 2^(-53).
        /// </summary>
        public const double DoubleMachineEpsilon = 1.11022302462516E-16;

        /// <summary>
        /// Value used for comparing single precision equivalence: the single-precision unit roundoff
        /// 2^(-24) = 5.9604644775390625E-08, exactly representable as a float.
        /// </summary>
        /// <remarks>
        /// This is the single-precision counterpart of <see cref="DoubleMachineEpsilon"/> (2^(-53)): half the
        /// spacing of floats just above one.
        /// </remarks>
        public const float SingleMachineEpsilon = 5.9604644775390625E-08f;

        /// <summary>
        /// Euler constant.
        /// </summary>
        public const double Euler = 0.5772156649015328606065120;

        /// <summary>
        /// The golden ratio. 
        /// </summary>
        public const double GoldenRatio = 1.61803398874989484820; // 0.5d * (1d + Math.Sqrt(5));

        /// <summary>
        /// Log of square root of 2*pi
        /// </summary>
        public const double LogSqrt2PI = 0.91893853320467274178032973640562;

        /// <summary>
        /// Square root of 2*pi
        /// </summary>
        public const double Sqrt2PI = 2.50662827463100050242E0;

        /// <summary>
        /// Square root of 2
        /// </summary>
        public const double Sqrt2 = 1.4142135623730950488016887;

        /// <summary>
        /// Log of 2.
        /// </summary>
        public const double Log2 = 0.69314718055994530941;


        /// <summary>
        /// Returns a value with the magnitude of x and the sign of y.
        /// </summary>
        public static double CopySign(double x, double y)
        {
            return (y >= 0.0) ? Math.Abs(x) : -Math.Abs(x);
        }

        /// <summary>
        /// Determines whether the specified value is finite
        /// </summary>
        /// <param name="x">The number to check.</param>
        public static bool IsFinite(double x) => !(double.IsNaN(x) || double.IsInfinity(x));

        /// <summary>
        /// Determine if integer is power of 2.
        /// </summary>
        /// <param name="n">Integer to test.</param>
        /// <returns>True if the integer is a power of 2.</returns>
        public static bool IsPowerOfTwo(int n)
        {
            return n > 0 && (n & (n - 1)) == 0;
        }

        /// <summary>
        /// Finds the next power of 2 greater than or equal to n.
        /// </summary>
        /// <param name="n">The input number.</param>
        /// <returns>The next power of 2.</returns>
        public static int NextPowerOfTwo(int n)
        {
            if (n < 1) return 1;

            // If already a power of 2, return it
            if ((n & (n - 1)) == 0)
                return n;

            // Find the next power of 2
            int power = 1;
            while (power < n)
            {
                power *= 2;
            }
            return power;
        }

        /// <summary>
        /// Returns value clamped to the inclusive range of min and max.
        /// </summary>
        /// <param name="x">The value to be clamped.</param>
        /// <param name="min">The lower bound of the result.</param>
        /// <param name="max">The upper bound of the result.</param>
        public static double Clamp(double x, double min, double max)
        {
            return x < min ? min : (x > max ? max : x);
        }

        /// <summary>
        /// Returns √(a² + b²) without destructive underflow or overflow.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns>The Euclidean norm of (a, b).</returns>
        /// <remarks>
        /// <para>
        /// The larger magnitude m is factored out, so the result is m·√(1 + r²) with r = min/max ≤ 1 (Press et al.,
        /// Numerical Recipes, 3rd ed., 2007, section 2.6). The ratio is squared by multiplication, which rounds
        /// correctly, and the result is within 2 ulps of the exact value.
        /// </para>
        /// <para>
        /// Special values follow IEEE 754 and C99 hypot: an infinite argument returns +∞, even when the other is NaN;
        /// otherwise a NaN argument returns NaN; and (0, 0) returns 0. The result is symmetric in the order and signs of
        /// the arguments.
        /// </para>
        /// </remarks>
        public static double Hypot(double a, double b)
        {
            double absa = Math.Abs(a), absb = Math.Abs(b);
            if (double.IsPositiveInfinity(absa) || double.IsPositiveInfinity(absb)) return double.PositiveInfinity;
            if (double.IsNaN(absa) || double.IsNaN(absb)) return double.NaN;
            if (absa > absb)
            {
                double r = absb / absa;
                return absa * Math.Sqrt(1d + r * r);
            }
            if (absb == 0d) return 0d;
            double s = absa / absb;
            return absb * Math.Sqrt(1d + s * s);
        }

        /// <summary>
        /// Returns a value with the same magnitude as a and the same sign as b.
        /// </summary>
        /// <param name="a">Value a.</param>
        /// <param name="b">Value b.</param>
        public static double Sign(double a, double b)
        {
            return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
        }

        /// <summary>
        /// Swap value a with value b.
        /// </summary>
        /// <param name="a">Value a.</param>
        /// <param name="b">Value b.</param>
        public static void Swap(ref double a, ref double b)
        {
            double temp = a;
            a = b;
            b = temp;
        }

        /// <summary>
        /// Returns the squared value of a.
        /// </summary>
        /// <param name="a">Value a.</param>
        public static double Sqr(double a)
        {
            return a * a;
        }

        /// <summary>
        /// Returns the value of a raised to the power of b.
        /// </summary>
        /// <param name="a">The value a.</param>
        /// <param name="b">The exponent b.</param>
        public static double Pow(double a, int b)
        {
            // Fast integer exponentiation with sensible edge handling.
            if (b == 0) return 1d;
            if (a == 1d) return 1d;
            if (a == -1d) return (b & 1) == 0 ? 1d : -1d;
            bool neg = b < 0;
            long n = b;
            if (neg) n = -n;
            double result = 1d;
            double baseVal = a;
            while (n > 0)
            {
                if ((n & 1L) != 0) result *= baseVal;
                baseVal *= baseVal;
                n >>= 1;
            }
            if (neg)
            {
                if (result == 0d) return double.PositiveInfinity; // 0^(-k)
                return 1d / result;
            }
            return result;
        }

        /// <summary>
        /// Returns the base 10 logarithm of a specified number. 
        /// </summary>
        /// <param name="x">The number whose logarithm is to be found.</param>
        /// <returns>The guarded base-10 logarithm of <paramref name="x"/>.</returns>
        /// <remarks>Values in [0, 1E-16) are raised to 1E-16, so <c>Log10(0)</c> returns -16 rather than negative infinity.</remarks>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        public static double Log10(double x)
        {
            if (x < 1E-16 && x >= 0d) x = 1E-16;
            return Math.Log10(x);
        }

        /// <summary>
        /// Returns the natural (base e) logarithm of a specified number. 
        /// </summary>
        /// <param name="x">The number whose logarithm is to be found.</param>
        /// <returns>The guarded natural logarithm of <paramref name="x"/>.</returns>
        /// <remarks>Values in [0, 1E-16), including negative zero, are raised to 1E-16, so <c>Log(0)</c> returns about -36.84 rather than negative infinity.</remarks>
        public static double Log(double x)
        {
            if (x < 1E-16 && Math.Sign(x) != -1) x = 1E-16;
            return Math.Log(x);
        }

        /// <summary>
        /// Computes <c>log(1 + x)</c> with improved numerical accuracy for small <paramref name="x"/>.
        /// </summary>
        /// <param name="x">
        /// The input value. The natural domain is <c>x &gt; -1</c>. 
        /// For <c>x = -1</c>, the result is <c>-∞</c>; for <c>x &lt; -1</c>, the result is <see cref="double.NaN"/>.
        /// </param>
        /// <returns>
        /// The natural logarithm of <c>1 + x</c>. For values of <paramref name="x"/> near zero,
        /// this method avoids catastrophic cancellation that occurs in <c>Math.Log(1 + x)</c>.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Applies the compensated evaluation <code>log(y) - ((y - 1) - x) / y</code> with <c>y = 1 + x</c> at every
        /// finite <c>x &gt; -1</c>: <c>(y - 1) - x</c> is exactly the rounding error of <c>1 + x</c>, and dividing it
        /// by <c>y</c> is the first-order correction of <c>log(y)</c> for that error, so the result is within one unit
        /// in the last place of the correctly rounded value wherever <c>log</c> is (it is not correctly rounded
        /// everywhere). The direct <c>log(1 + x)</c> loses up to about 1e-12 relative
        /// accuracy for <c>1e-4 &lt; |x| &lt; 1e-2</c>, where the rounding of <c>1 + x</c> is a large share of the
        /// result.
        /// </para>
        /// <para>
        /// Where <c>1 + x</c> is zero, negative, infinite or NaN the correction is not defined, and the value is
        /// <see cref="Math.Log(double)"/> of <c>1 + x</c>: −∞ at <c>x = −1</c>, NaN below it or for NaN, and +∞
        /// for +∞.
        /// </para>
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Log1p(double x)
        {
            double y = 1.0 + x;
            if (!(y > 0d) || double.IsInfinity(y))
                return Math.Log(y);
            return Math.Log(y) - ((y - 1.0) - x) / y;
        }

        /// <summary>
        /// Computes <c>exp(x) - 1</c> with improved numerical accuracy for small <paramref name="x"/>.
        /// </summary>
        /// <param name="x">The input value.</param>
        /// <returns>
        /// <c>exp(x) - 1</c>. For values of <paramref name="x"/> near zero, this method avoids the
        /// catastrophic cancellation that occurs in <c>Math.Exp(x) - 1</c>; deeply negative inputs
        /// return exactly -1 and large inputs overflow to positive infinity.
        /// </returns>
        /// <remarks>
        /// <para>
        /// For <c>|x| &gt; 1/2</c> the direct difference <c>exp(x) - 1</c> is returned. There the subtraction
        /// cancels at most a factor of about 2.6, so the result stays within about one unit in the last place
        /// (measured against mpmath: at most 1.06 units), which is closer than the compensated form below,
        /// whose division by <c>log(u)</c> adds its own rounding. The difference also gives exactly -1 once
        /// <c>exp(x)</c> falls below half a unit of 1, and positive infinity when <c>exp(x)</c> overflows.
        /// </para>
        /// <para>
        /// For <c>|x| &lt;= 1/2</c>, and for NaN, it uses the compensated evaluation
        /// <c>(u - 1) * x / log(u)</c> with <c>u = exp(x)</c>, which corrects the rounding of the exponential
        /// where the subtraction would cancel; when <c>u</c> rounds to one the input itself is returned. This is
        /// the companion of <see cref="Log1p(double)"/> for log-space probability arithmetic such as survival
        /// products of many small probabilities.
        /// </para>
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Expm1(double x)
        {
            if (Math.Abs(x) > 0.5) return Math.Exp(x) - 1.0;
            double u = Math.Exp(x);
            if (u == 1.0) return x;
            return (u - 1.0) * x / Math.Log(u);
        }

        /// <summary>Computes log(1-exp(a)) for a nonpositive log probability, including its limits.</summary>
        /// <param name="a">The logarithm of a probability in the interval from negative infinity through zero.</param>
        /// <returns><c>log(1-exp(a))</c>, or not-a-number when <paramref name="a"/> is positive or not-a-number.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static double Log1mExp(double a)
        {
            if (a > 0 || double.IsNaN(a)) return double.NaN;
            return a < -0.69314718055994530942 ? Tools.Log1p(-Math.Exp(a)) : Math.Log(-Tools.Expm1(a));
        }

        /// <summary>Cancellation-free exp(x)-1-x.</summary>
        /// <param name="x">The argument.</param>
        /// <returns><c>exp(x)-1-x</c>, evaluated by its Taylor series when the subtraction would cancel.</returns>
        /// <remarks>For |x| below 0.5 the series <c>Σₙ₌₂ xⁿ/n!</c> is summed to a relative tolerance of 1E-17; otherwise
        /// <c>expm1(x) - x</c> loses at most a factor of about 4.7 to cancellation.</remarks>
        internal static double Expm1MinusX(double x)
        {
            if (Math.Abs(x) >= 0.5) return Tools.Expm1(x) - x;
            double term = x * x / 2, sum = term;
            for (int n = 3; n < 60; n++)
            {
                term *= x / n;
                sum += term;
                if (Math.Abs(term) <= Math.Abs(sum) * 1E-17) break;
            }
            return sum;
        }

        /// <summary>Cancellation-free log(1+x)-x.</summary>
        /// <param name="x">The argument, greater than negative one in its probability-distribution uses.</param>
        /// <returns><c>log(1+x)-x</c>, evaluated by a local series when subtraction would cancel.</returns>
        internal static double Log1pMinusX(double x)
        {
            if (Math.Abs(x) >= 0.25) return Tools.Log1p(x) - x;
            double power = -x * x, sum = power / 2;
            for (int n = 3; n < 100; n++)
            {
                power *= -x;
                double term = power / n;
                sum += term;
                if (Math.Abs(term) <= Math.Abs(sum) * 1E-17) break;
            }
            return sum;
        }

        /// <summary>
        /// Computes <c>log(1 + x)</c> for each value of a list, with improved numerical accuracy for small values.
        /// </summary>
        /// <param name="values">The input values.</param>
        /// <returns>One result per value, each equal to <see cref="Log1p(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as the list forms of the distribution members
        /// (for example <see cref="Distributions.Normal.StandardZ(IList{double})"/>); a caller crossing a language
        /// boundary (for example from Python) pays the crossing once instead of once per value.
        /// </remarks>
        public static double[] Log1p(IList<double> values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = Log1p(values[i]);
            return result;
        }

        /// <summary>
        /// Computes <c>exp(x) - 1</c> for each value of a list, with improved numerical accuracy for small values.
        /// </summary>
        /// <param name="values">The input values.</param>
        /// <returns>One result per value, each equal to <see cref="Expm1(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Log1p(IList{double})"/>.
        /// </remarks>
        public static double[] Expm1(IList<double> values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = Expm1(values[i]);
            return result;
        }

        /// <summary>
        /// Returns the Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x1">X of point 1.</param>
        /// <param name="y1">Y of point 1.</param>
        /// <param name="x2">X of point 2.</param>
        /// <param name="y2">Y of point 2.</param>
        public static double Distance(double x1, double y1, double x2, double y2)
        {   
            double dx = x2 - x1;
            double dy = y2 - y1;
            return Math.Sqrt(dx * dx + dy * dy);
        }

        /// <summary>
        /// Returns the Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <returns>The Euclidean distance ||x - y||.</returns>
        /// <remarks>
        /// Arrays of compatible length are read directly, without interface dispatch, using the same
        /// operations in the same order.
        /// </remarks>
        public static double Distance(IList<double> x, IList<double> y)
        {
            double d = 0;
            if (x is double[] a && y is double[] b && b.Length >= a.Length)
            {
                for (int i = 0; i < a.Length; i++)
                {
                    double dx = a[i] - b[i];
                    d += dx * dx;
                }
                return Math.Sqrt(d);
            }
            for (int i = 0; i < x.Count; i++)
            {
                double dx = x[i] - y[i];
                d += dx * dx;
            }
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the Euclidean distance between two points after applying min-max normalization to each dimension.
        /// Normalization ensures each parameter dimension contributes equally to the distance metric,
        /// based on the specified lower and upper bounds.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <param name="lower">The lower bounds for each parameter dimension.</param>
        /// <param name="upper">The upper bounds for each parameter dimension.</param>
        /// <returns>The normalized Euclidean distance ||x - y||.</returns>
        public static double NormalizedDistance(IList<double> x, IList<double> y, IList<double> lower, IList<double> upper)
        {
            double d = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double range = upper[i] - lower[i];
                if (range <= 0d || double.IsNaN(range))
                {
                    // Degenerate dimension; contribute nothing (identical after normalization).
                    continue;
                }
                double xi = (x[i] - lower[i]) / range;
                double yi = (y[i] - lower[i]) / range;
                double dx = xi - yi;
                d += dx * dx;
            }
            return Math.Sqrt(d);
        }


        /// <summary>
        /// Returns the weighted Euclidean distance between two points ||x - y||.
        /// </summary>
        /// <param name="x">First point.</param>
        /// <param name="y">Second point.</param>
        /// <param name="weights">Weights.</param>
        /// <returns>The weighted Euclidean distance ||x - y||.</returns>
        public static double WeightedDistance(IList<double> x, IList<double> y, IList<double> weights)
        {
            double d = 0;
            for (int i = 0; i < x.Count; i++)
            {
                double dx = x[i] - y[i];
                d += (dx * dx) * weights[i];
            }
            return Math.Sqrt(d);
        }

        /// <summary>
        /// Returns the normalized values ranging between 0 and 1.
        /// <para> 
        /// <see href="https://www.codecademy.com/article/normalization"/>
        /// </para>
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double[] Normalize(IList<double> values)
        {
            var result = new double[values.Count];
            double min = 0, max = 0;
            MinMax(values, out min, out max);
            if (double.IsNaN(min) || double.IsNaN(max) || values.Count == 0)
                return result; // return zeros for empty/NaN input
            double range = max - min;
            if (range == 0d)
            {
                // All identical: define normalized values as 0.
                return result;
            }
            for (int i = 0; i < values.Count; i++)
                result[i] = (values[i] - min) / range;
            return result;
        }

        /// <summary>
        /// Returns the denormalized values. 
        /// </summary>
        /// <param name="values">The list of normalized values to denormalize.</param>
        /// <param name="min">The minimum value from the original data.</param>
        /// <param name="max">The maximum value from the original data.</param>
        public static double[] Denormalize(IList<double> values, double min, double max)
        {
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = values[i] * (max - min) + min;
            return result;
        }


        /// <summary>
        /// Returns the standardized values. 
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double[] Standardize(IList<double> values)
        {
            var result = new double[values.Count];
            var meanSD = Statistics.MeanStandardDeviation(values);
            double mu = meanSD.mean;
            double sd = meanSD.standardDeviation;
            if (sd <= 0d || double.IsNaN(sd))
            {
                // Degenerate variance 
                return result;
            }
            for (int i = 0; i < values.Count; i++)
                result[i] = (values[i] - mu) / sd;
            return result;
        }

        /// <summary>
        /// Returns the destandardized values. 
        /// </summary>
        /// <param name="values">The list of standardized values to destandardize.</param>
        /// <param name="mean">The mean of the original data.</param>
        /// <param name="standardDeviation">The standard deviation of the original data.</param>
        public static double[] Destandardize(IList<double> values, double mean, double standardDeviation)
        {
            var result = new double[values.Count];
            for (int i = 0; i < values.Count; i++)
                result[i] = values[i] * standardDeviation + mean;
            return result;
        }

        /// <summary>
        /// Estimates the sum of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Sum(IList<int> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum;
        }

        /// <summary>
        /// Estimates the sum of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Sum(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum;
        }

        /// <summary>
        /// Estimate the sum of a list of values. 
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param>
        public static double Sum(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                if (indicators[i] == (useComplement ? 0 : 1)) sum += values[i];      
            return sum;
        }

        /// <summary>
        /// Returns the sum product of two lists of values.
        /// </summary>
        /// <param name="values1">The first list of values.</param>
        /// <param name="values2">The second list of values.</param>
        public static double SumProduct(IList<double> values1, IList<double> values2)
        {
            if (values1.Count == 0) return double.NaN;
            if (values2.Count != values1.Count) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values1.Count; i++)
                sum += values1[i] * values2[i];
            return sum;
        }

        /// <summary>
        /// Returns the compensated (Neumaier) sum of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <returns>The compensated sum, or NaN when the list is empty.</returns>
        /// <remarks>
        /// <para>
        /// The values are added in index order with Neumaier's improvement of Kahan summation: a running
        /// correction collects the rounding error of each addition, taken exactly from whichever of the running
        /// sum and the new value is larger in magnitude, and is added to the running sum at the end. Unlike plain
        /// Kahan summation, the result stays accurate when a term is larger than the running sum. The error is
        /// bounded by 2u|s| + O(n u²) Σ|x_i| (u the unit roundoff, s the exact sum), so to first order it does not
        /// grow with n, whereas the bound of <see cref="Sum(IList{double})"/> grows as (n − 1) u Σ|x_i|.
        /// </para>
        /// <para>
        /// If the running sum is not finite, the IEEE sum is returned: NaN when any value is NaN or both infinities
        /// occur, otherwise the infinity. Python's <c>math.fsum</c> gives the same NaN or infinity, but raises where
        /// both infinities occur or finite values overflow, and returns 0 for an empty list, which returns NaN here as
        /// for <see cref="Sum(IList{double})"/>.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen.
        /// ZAMM - Zeitschrift für Angewandte Mathematik und Mechanik, 54(1), 39-51. https://doi.org/10.1002/zamm.19740540106
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public static double CompensatedSum(IList<double> values)
        {
            return CompensatedSum(values, false);
        }

        /// <summary>
        /// Returns the compensated sum of a list of values by Neumaier's method, or by classic Kahan summation.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="kahan">False for Neumaier's method, as <see cref="CompensatedSum(IList{double})"/>; true for
        /// classic Kahan summation.</param>
        /// <returns>The compensated sum, or NaN when the list is empty.</returns>
        /// <remarks>
        /// <para>
        /// With <paramref name="kahan"/> false the result is bitwise <see cref="CompensatedSum(IList{double})"/>.
        /// </para>
        /// <para>
        /// With <paramref name="kahan"/> true the values are added by classic Kahan summation instead: starting from
        /// t = c = 0, each value v gives a = v − c, f = t + a, c = (f − t) − a and t = f, and the result is t. Neumaier's
        /// method is the more accurate of the two, because Kahan's correction is lost when a term is larger than the
        /// running sum; the Kahan form is provided for reproducing results of the classic Kahan loop. When the running
        /// sum first becomes infinite the correction becomes NaN: a final infinite sum is returned as is, and any value
        /// added after it gives NaN, as in the classic loop. Neumaier's method returns the IEEE sum in both cases.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public static double CompensatedSum(IList<double> values, bool kahan)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count == 0) return double.NaN;
            double sum = 0d, correction = 0d;
            if (values is double[] array)
            {
                // The same additions without the interface call per value.
                if (kahan)
                    for (int i = 0; i < array.Length; i++) KahanAdd(ref sum, ref correction, array[i]);
                else
                    for (int i = 0; i < array.Length; i++) NeumaierAdd(ref sum, ref correction, array[i]);
            }
            else if (kahan)
            {
                for (int i = 0; i < values.Count; i++) KahanAdd(ref sum, ref correction, values[i]);
            }
            else
            {
                for (int i = 0; i < values.Count; i++) NeumaierAdd(ref sum, ref correction, values[i]);
            }
            return CompensatedResult(sum, correction, kahan);
        }

        /// <summary>
        /// Returns the compensated (Neumaier) sum of the values an indicator selects.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param>
        /// <returns>The compensated sum of the selected values, 0 when none is selected, or NaN when the list is empty
        /// or the indicators differ in length from the values.</returns>
        /// <remarks>
        /// The selection follows <see cref="Sum(IList{double}, IList{int}, bool)"/>, and the selected values are summed
        /// in index order as in <see cref="CompensatedSum(IList{double})"/>.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> or <paramref name="indicators"/> is null.</exception>
        public static double CompensatedSum(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            return CompensatedSum(values, indicators, useComplement, false);
        }

        /// <summary>
        /// Returns the compensated sum of the values an indicator selects, by Neumaier's method or by classic Kahan
        /// summation.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator.</param>
        /// <param name="kahan">False for Neumaier's method; true for classic Kahan summation.</param>
        /// <returns>The compensated sum of the selected values, 0 when none is selected, or NaN when the list is empty
        /// or the indicators differ in length from the values.</returns>
        /// <remarks>
        /// <para>
        /// The selected values are added in index order as in <see cref="CompensatedSum(IList{double}, bool)"/>; with
        /// <paramref name="kahan"/> false the result is bitwise <see cref="CompensatedSum(IList{double}, IList{int}, bool)"/>.
        /// </para>
        /// <para>
        /// With <paramref name="kahan"/> true the values are added by classic Kahan summation instead: starting from
        /// t = c = 0, each value v gives a = v − c, f = t + a, c = (f − t) − a and t = f, and the result is t. Neumaier's
        /// method is the more accurate of the two, because Kahan's correction is lost when a term is larger than the
        /// running sum; the Kahan form is provided for reproducing results of the classic Kahan loop. When the running
        /// sum first becomes infinite the correction becomes NaN: a final infinite sum is returned as is, and any value
        /// added after it gives NaN, as in the classic loop. Neumaier's method returns the IEEE sum in both cases.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> or <paramref name="indicators"/> is null.</exception>
        public static double CompensatedSum(IList<double> values, IList<int> indicators, bool useComplement, bool kahan)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (indicators == null) throw new ArgumentNullException(nameof(indicators));
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            int selected = useComplement ? 0 : 1;
            double sum = 0d, correction = 0d;
            for (int i = 0; i < values.Count; i++)
            {
                if (indicators[i] != selected) continue;
                if (kahan) KahanAdd(ref sum, ref correction, values[i]);
                else NeumaierAdd(ref sum, ref correction, values[i]);
            }
            return CompensatedResult(sum, correction, kahan);
        }

        /// <summary>
        /// Returns the compensated (Neumaier) sum of the values at a range of gathered indices,
        /// values[indices[i]] for i in [start, start + count).
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indices">The 0-based positions in <paramref name="values"/> to sum; a position may repeat.</param>
        /// <param name="start">The first entry of <paramref name="indices"/> to use.</param>
        /// <param name="count">The number of entries of <paramref name="indices"/> to use.</param>
        /// <returns>The compensated sum of the gathered values, or NaN when <paramref name="count"/> is zero.</returns>
        /// <remarks>
        /// The gathered values are summed in the order of <paramref name="indices"/> as in
        /// <see cref="CompensatedSum(IList{double})"/>, without copying them to a new list.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> or <paramref name="indices"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="start"/> or <paramref name="count"/>
        /// is negative, the range runs past the end of <paramref name="indices"/>, or a gathered index is outside
        /// <paramref name="values"/>.</exception>
        public static double CompensatedSum(IList<double> values, IList<int> indices, int start, int count)
        {
            return CompensatedSum(values, indices, start, count, false);
        }

        /// <summary>
        /// Returns the compensated sum of the values at a range of gathered indices, values[indices[i]] for i in
        /// [start, start + count), by Neumaier's method or by classic Kahan summation.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indices">The 0-based positions in <paramref name="values"/> to sum; a position may repeat.</param>
        /// <param name="start">The first entry of <paramref name="indices"/> to use.</param>
        /// <param name="count">The number of entries of <paramref name="indices"/> to use.</param>
        /// <param name="kahan">False for Neumaier's method; true for classic Kahan summation.</param>
        /// <returns>The compensated sum of the gathered values, or NaN when <paramref name="count"/> is zero.</returns>
        /// <remarks>
        /// <para>
        /// The gathered values are added in the order of <paramref name="indices"/> as in
        /// <see cref="CompensatedSum(IList{double}, bool)"/>; with <paramref name="kahan"/> false the result is bitwise
        /// <see cref="CompensatedSum(IList{double}, IList{int}, int, int)"/>.
        /// </para>
        /// <para>
        /// With <paramref name="kahan"/> true the values are added by classic Kahan summation instead: starting from
        /// t = c = 0, each value v gives a = v − c, f = t + a, c = (f − t) − a and t = f, and the result is t. Neumaier's
        /// method is the more accurate of the two, because Kahan's correction is lost when a term is larger than the
        /// running sum; the Kahan form is provided for reproducing results of the classic Kahan loop. When the running
        /// sum first becomes infinite the correction becomes NaN: a final infinite sum is returned as is, and any value
        /// added after it gives NaN, as in the classic loop. Neumaier's method returns the IEEE sum in both cases.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> or <paramref name="indices"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="start"/> or <paramref name="count"/>
        /// is negative, the range runs past the end of <paramref name="indices"/>, or a gathered index is outside
        /// <paramref name="values"/>.</exception>
        public static double CompensatedSum(IList<double> values, IList<int> indices, int start, int count, bool kahan)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (indices == null) throw new ArgumentNullException(nameof(indices));
            if (start < 0) throw new ArgumentOutOfRangeException(nameof(start), "The start must not be negative.");
            if (count < 0 || count > indices.Count - start)
                throw new ArgumentOutOfRangeException(nameof(count), "The range must lie within the indices.");
            if (count == 0) return double.NaN;
            int n = values.Count;
            double sum = 0d, correction = 0d;
            for (int i = start; i < start + count; i++)
            {
                int j = indices[i];
                if ((uint)j >= (uint)n)
                    throw new ArgumentOutOfRangeException(nameof(indices), "Each gathered index must lie within the values.");
                if (kahan) KahanAdd(ref sum, ref correction, values[j]);
                else NeumaierAdd(ref sum, ref correction, values[j]);
            }
            return CompensatedResult(sum, correction, kahan);
        }

        /// <summary>
        /// Returns the compensated (Neumaier) cumulative sums of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <returns>An array whose entry i is the compensated sum of values[0] through values[i]; empty when the list is empty.</returns>
        /// <remarks>
        /// One pass carries the running sum and correction of <see cref="CompensatedSum(IList{double})"/>, so entry i
        /// equals <see cref="CompensatedSum(IList{double})"/> of the first i + 1 values, bit for bit, including the
        /// IEEE result once the running sum is not finite.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public static double[] CompensatedCumulativeSum(IList<double> values)
        {
            return CompensatedCumulativeSum(values, false);
        }

        /// <summary>
        /// Returns the compensated cumulative sums of a list of values, by Neumaier's method or by classic Kahan summation.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="kahan">False for Neumaier's method, as <see cref="CompensatedCumulativeSum(IList{double})"/>; true
        /// for classic Kahan summation.</param>
        /// <returns>An array whose entry i is the compensated sum of values[0] through values[i]; empty when the list is empty.</returns>
        /// <remarks>
        /// <para>
        /// Entry i equals <see cref="CompensatedSum(IList{double}, bool)"/> of the first i + 1 values with the same
        /// <paramref name="kahan"/>, bit for bit; for Kahan it is the running t after value i.
        /// </para>
        /// <para>
        /// With <paramref name="kahan"/> true the values are added by classic Kahan summation instead: starting from
        /// t = c = 0, each value v gives a = v − c, f = t + a, c = (f − t) − a and t = f, and the result is t. Neumaier's
        /// method is the more accurate of the two, because Kahan's correction is lost when a term is larger than the
        /// running sum; the Kahan form is provided for reproducing results of the classic Kahan loop. When the running
        /// sum first becomes infinite the correction becomes NaN: a final infinite sum is returned as is, and any value
        /// added after it gives NaN, as in the classic loop. Neumaier's method returns the IEEE sum in both cases.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public static double[] CompensatedCumulativeSum(IList<double> values, bool kahan)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            var result = new double[values.Count];
            double sum = 0d, correction = 0d;
            if (values is double[] array)
            {
                // The same additions without the interface call per value.
                if (kahan)
                {
                    for (int i = 0; i < array.Length; i++)
                    {
                        KahanAdd(ref sum, ref correction, array[i]);
                        result[i] = sum;
                    }
                }
                else
                {
                    for (int i = 0; i < array.Length; i++)
                    {
                        NeumaierAdd(ref sum, ref correction, array[i]);
                        result[i] = NeumaierResult(sum, correction);
                    }
                }
                return result;
            }
            for (int i = 0; i < values.Count; i++)
            {
                if (kahan) KahanAdd(ref sum, ref correction, values[i]);
                else NeumaierAdd(ref sum, ref correction, values[i]);
                result[i] = CompensatedResult(sum, correction, kahan);
            }
            return result;
        }

        /// <summary>
        /// Adds one value to a Neumaier running sum and its compensation.
        /// </summary>
        /// <param name="sum">The running (uncompensated) sum, updated to sum + value.</param>
        /// <param name="compensation">The accumulated rounding error of the running sum, updated with the error of
        /// this addition.</param>
        /// <param name="value">The value to add.</param>
        /// <remarks>
        /// <para>
        /// This is the step of <see cref="CompensatedSum(IList{double})"/>. Start from sum = compensation = 0, call it
        /// once per value in order, and read the result with <see cref="NeumaierResult(double, double)"/>; the result is
        /// bitwise the one <see cref="CompensatedSum(IList{double})"/> returns for the same values in the same order.
        /// </para>
        /// <para>
        /// With t = sum + value, the rounding error of the addition is recovered exactly by the Fast2Sum step with the
        /// operands ordered by magnitude: (sum − t) + value when |sum| ≥ |value|, otherwise (value − t) + sum. Unlike
        /// <see cref="KahanAdd(ref double, ref double, double)"/>, the error is not lost when a value is larger than the
        /// running sum.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Neumaier, A. (1974). Rundungsfehleranalyse einiger Verfahren zur Summation endlicher Summen.
        /// ZAMM - Zeitschrift für Angewandte Mathematik und Mechanik, 54(1), 39-51. https://doi.org/10.1002/zamm.19740540106
        /// </para>
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void NeumaierAdd(ref double sum, ref double compensation, double value)
        {
            double t = sum + value;
            if (Math.Abs(sum) >= Math.Abs(value)) compensation += (sum - t) + value;
            else compensation += (value - t) + sum;
            sum = t;
        }

        /// <summary>
        /// Returns the result of a Neumaier running sum formed by <see cref="NeumaierAdd(ref double, ref double, double)"/>.
        /// </summary>
        /// <param name="sum">The running sum.</param>
        /// <param name="compensation">The running compensation.</param>
        /// <returns>The sum plus its compensation, or the IEEE sum as it stands once it is not finite.</returns>
        /// <remarks>
        /// Once the running sum is infinite or NaN the compensation carries no information (it is NaN or infinite
        /// itself), so the IEEE sum is returned: NaN when any value was NaN or both infinities occurred, otherwise the
        /// infinity. This is the final step of <see cref="CompensatedSum(IList{double})"/>.
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double NeumaierResult(double sum, double compensation)
        {
            return IsFinite(sum) ? sum + compensation : sum;
        }

        /// <summary>
        /// Adds one value to a classic Kahan running sum and its compensation.
        /// </summary>
        /// <param name="sum">The running compensated sum t, which is also the result.</param>
        /// <param name="compensation">The compensation c, the negated low-order part lost by the last addition.</param>
        /// <param name="value">The value v to add.</param>
        /// <remarks>
        /// <para>
        /// The step is a = v − c, f = t + a, c = (f − t) − a, t = f (Kahan, 1965). Start from sum = compensation = 0
        /// and call it once per value in order; <paramref name="sum"/> is then bitwise the result of
        /// <see cref="CompensatedSum(IList{double}, bool)"/> with kahan true for the same values in the same order.
        /// </para>
        /// <para>
        /// Neumaier's step <see cref="NeumaierAdd(ref double, ref double, double)"/> is the more accurate of the two,
        /// because Kahan's compensation is lost when a value is larger than the running sum; this step is provided for
        /// reproducing results of the classic Kahan loop. Once the sum becomes infinite the compensation becomes NaN, and
        /// any value added after it gives NaN.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Kahan, W. (1965). Further remarks on reducing truncation errors. Communications of the ACM, 8(1), 40.
        /// https://doi.org/10.1145/363707.363723
        /// </para>
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void KahanAdd(ref double sum, ref double compensation, double value)
        {
            double a = value - compensation;
            double f = sum + a;
            compensation = (f - sum) - a;
            sum = f;
        }

        /// <summary>
        /// Returns the result of a compensated running sum.
        /// </summary>
        /// <param name="sum">The running sum.</param>
        /// <param name="correction">The running correction.</param>
        /// <param name="kahan">True when the sum was formed by <see cref="KahanAdd(ref double, ref double, double)"/>.</param>
        /// <returns>The Kahan running sum as it stands; for Neumaier, <see cref="NeumaierResult(double, double)"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double CompensatedResult(double sum, double correction, bool kahan)
        {
            return kahan ? sum : NeumaierResult(sum, correction);
        }


        /// <summary>
        /// Estimates the mean of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Mean(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double sum = 0d;
            for (int i = 0; i < values.Count; i++)
                sum += values[i];
            return sum / values.Count;
        }

        /// <summary>
        /// Estimates the mean of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param>
        public static double Mean(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double sum = 0d;
            int count = 0;
            for (int i = 0; i < values.Count; i++)
            {
                if (indicators[i] == (useComplement ? 0 : 1))
                {
                    sum += values[i];
                    count += 1;
                }
            }
            return count == 0 ? double.NaN : sum / count;
        }

        /// <summary>
        /// Estimates the product of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Product(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double product = values[0];
            if (product == 0) return 0d;
            for (int i = 1; i < values.Count; i++)
            {
                product *= values[i];
                if (product == 0) return 0d;
            }
            return product;
        }

        /// <summary>
        /// Estimates the product of a list of values.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Product(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double product = 1;
            int flag = useComplement ? 0 : 1;
            for (int i = 0; i < values.Count; i++)
            {
                if (indicators[i] == flag)
                {
                    product *= values[i];
                    if (product == 0) return 0d;
                }
            }
            return product;
        }

        /// <summary>
        /// Returns the minimum and maximum values from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="min">Output. Minimum value.</param>
        /// <param name="max">Output. Maximum value.</param>
        public static void MinMax(IList<double> values, out double min, out double max)
        {
            min = double.PositiveInfinity;
            max = double.NegativeInfinity;
            if (values.Count == 0) { min = double.NaN; max = double.NaN; return; };
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (double.IsNaN(v)) { min = double.NaN; max = double.NaN; return; }
                if (v < min) min = v;
                if (v > max) max = v;
            }

        }

        /// <summary>
        /// Returns the index of the minimum value. NaN entries are ignored.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <returns>The zero-based index of the smallest value, or -1 if the list is empty or contains only NaN.</returns>
        public static int ArgMin(IList<double> values)
        {
            double min = double.PositiveInfinity;
            int index = -1;
            for (int i = 0; i < values.Count; i++)
            {
                if (!double.IsNaN(values[i]) && (index == -1 || values[i] < min))
                {
                    min = values[i];
                    index = i;
                }
            }
            return index;
        }

        /// <summary>
        /// Returns the index of the maximum value. NaN entries are ignored.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <returns>The zero-based index of the largest value, or -1 if the list is empty or contains only NaN.</returns>
        public static int ArgMax(IList<double> values)
        {
            double max = double.NegativeInfinity;
            int index = -1;
            for (int i = 0; i < values.Count; i++)
            {
                if (!double.IsNaN(values[i]) && (index == -1 || values[i] > max))
                {
                    max = values[i];
                    index = i;
                }
            }
            return index;
        }

        /// <summary>
        /// Returns the smallest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Min(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double min = double.PositiveInfinity;
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (double.IsNaN(v)) { return double.NaN; }
                if (v < min) min = v;        
            }
            return min;
        }

        /// <summary>
        /// Returns the smallest indicated value from a list of values.
        /// Returns NaN if the list is empty, the indicator list length differs, no entry is
        /// indicated, or any indicated entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Min(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double min = double.PositiveInfinity;
            bool any = false;
            int flag = useComplement ? 0 : 1;
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (indicators[i] == flag)
                {
                    if (double.IsNaN(v)) return double.NaN;
                    if (v < min) min = v;
                    any = true;
                } 
            }
            return any? min : double.NaN;
        }

        /// <summary>
        /// Returns the largest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        public static double Max(IList<double> values)
        {
            if (values.Count == 0) return double.NaN;
            double max = double.NegativeInfinity;
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (double.IsNaN(v)) { return double.NaN; }
                if (v > max) max = v;
            }
            return max;
        }

        /// <summary>
        /// Returns the largest value from a list of values.
        /// Returns NaN if the list is empty or any entry is NaN.
        /// </summary>
        /// <param name="values">The list of values.</param>
        /// <param name="indicators">The list of indicators (0's or 1's).</param>
        /// <param name="useComplement">If false, 1's are indicator, if true, all 0's are indicator. Default = false.</param> 
        public static double Max(IList<double> values, IList<int> indicators, bool useComplement = false)
        {
            if (values.Count == 0) return double.NaN;
            if (indicators.Count != values.Count) return double.NaN;
            double max = double.NegativeInfinity;
            bool any = false;
            int flag = useComplement ? 0 : 1;
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (indicators[i] == flag)
                {
                    if (double.IsNaN(v)) return double.NaN;
                    if (v > max) max = v;
                    any = true;
                }
            }
            return any ? max : double.NaN;
        }

        /// <summary>
        /// Supporting function used to add doubles from an interlocked parallel loop.
        /// </summary>
        /// <param name="valueToAddTo">The accumulator to update atomically.</param>
        /// <param name="valueToAdd">The value to add.</param>
        /// <returns>The accumulator value produced by the atomic addition.</returns>
        /// <remarks>
        /// All concurrent writes to <paramref name="valueToAddTo"/> must use this method or another
        /// interlocked operation. The compare-and-swap operation can retry under contention. Concurrent
        /// updates are not lost, but their order depends on thread scheduling, so floating-point reductions
        /// using this method are not guaranteed to be bit-reproducible.
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double ParallelAdd(ref double valueToAddTo, double valueToAdd)
        {
            double currentValue = valueToAddTo;
            double newValue = currentValue + valueToAdd;
            double observedValue = Interlocked.CompareExchange(ref valueToAddTo, newValue, currentValue);
            if (observedValue == currentValue &&
                (currentValue != 0d ||
                 BitConverter.DoubleToInt64Bits(observedValue) == BitConverter.DoubleToInt64Bits(currentValue)))
            {
                return newValue;
            }
            if (double.IsNaN(observedValue) && double.IsNaN(currentValue)) return newValue;
            return ParallelAddRetry(ref valueToAddTo, valueToAdd, observedValue);
        }

        /// <summary>
        /// Retries an atomic double addition after a failed compare-and-swap operation.
        /// </summary>
        /// <param name="valueToAddTo">The accumulator to update atomically.</param>
        /// <param name="valueToAdd">The value to add.</param>
        /// <param name="observedValue">The accumulator value observed by the failed operation.</param>
        /// <returns>The accumulator value produced by the atomic addition.</returns>
        /// <remarks>
        /// The observed value seeds the lock-free retry loop so that another non-volatile read is unnecessary.
        /// Signed zeros require a representation comparison because numerical equality does not distinguish
        /// positive zero from negative zero, while compare-and-swap does.
        /// </remarks>
        private static double ParallelAddRetry(ref double valueToAddTo, double valueToAdd, double observedValue)
        {
            double newCurrentValue = observedValue;
            while (true)
            {
                double currentValue = newCurrentValue;
                double newValue = currentValue + valueToAdd;
                newCurrentValue = Interlocked.CompareExchange(ref valueToAddTo, newValue, currentValue);
                if (newCurrentValue == currentValue &&
                    (currentValue != 0d ||
                     BitConverter.DoubleToInt64Bits(newCurrentValue) == BitConverter.DoubleToInt64Bits(currentValue)))
                {
                    return newValue;
                }
                if (double.IsNaN(newCurrentValue) && double.IsNaN(currentValue)) return newValue;
            }
        }

        /// <summary>
        /// Returns log(exp(u) + exp(v)) without overflow or underflow.
        /// </summary>
        /// <param name="u">The first log value.</param>
        /// <param name="v">The second log value.</param>
        /// <returns>log(exp(u) + exp(v)); NaN when either value is NaN, +∞ when either is +∞ and neither is NaN,
        /// and −∞ when both are −∞.</returns>
        /// <remarks>
        /// <para>
        /// With m the larger and n the smaller value, the result is log1p(exp(n − m)) + m, or log1p(0) + log(2) + m when
        /// the two are equal; it is bitwise <see cref="LogSumExp(IList{double})"/> of the two values. Separating the
        /// largest term and using <see cref="Log1p(double)"/> keeps the result accurate when it is near zero, where
        /// m + log(1 + exp(n − m)) loses relative accuracy (Blanchard, Higham and Higham, 2021).
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Blanchard, P., Higham, D. J., &amp; Higham, N. J. (2021). Accurately computing the log-sum-exp and softmax
        /// functions. IMA Journal of Numerical Analysis, 41(4), 2311-2330. https://doi.org/10.1093/imanum/draa038
        /// </para>
        /// </remarks>
        public static double LogSumExp(double u, double v)
        {
            if (double.IsNaN(u) || double.IsNaN(v)) return double.NaN;
            double max = Math.Max(u, v);
            if (double.IsInfinity(max)) return max;
            if (u == v) return LogSumExpResult(max, 2, 0d);
            return LogSumExpResult(max, 1, Math.Exp(Math.Min(u, v) - max));
        }

        /// <summary>
        /// Returns log(Σ exp(x_i)) of a list of log values without overflow or underflow.
        /// </summary>
        /// <param name="values">The log values.</param>
        /// <returns>The log of the sum of the exponentials; −∞ for an empty list or when every value is −∞, NaN when any
        /// value is NaN, and +∞ when any value is +∞ and none is NaN.</returns>
        /// <remarks>
        /// <para>
        /// With m the largest value and k the number of values equal to it, the result is log1p(s / k) + log(k) + m,
        /// where s is the compensated sum (<see cref="NeumaierAdd(ref double, ref double, double)"/>) of exp(x_i − m)
        /// over the other values, added in index order. This is the evaluation of <c>scipy.special.logsumexp</c>:
        /// separating the largest terms and using <see cref="Log1p(double)"/> keeps the result accurate when it is near
        /// zero, where m + log(Σ exp(x_i − m)) loses relative accuracy (Blanchard, Higham and Higham, 2021).
        /// </para>
        /// <para>
        /// The special values also follow SciPy: an empty list returns −∞, the log of an empty sum, as does a list of
        /// −∞ values; +∞ is returned when any value is +∞, since the sum is then infinite; and NaN propagates.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Blanchard, P., Higham, D. J., &amp; Higham, N. J. (2021). Accurately computing the log-sum-exp and softmax
        /// functions. IMA Journal of Numerical Analysis, 41(4), 2311-2330. https://doi.org/10.1093/imanum/draa038
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        public static double LogSumExp(IList<double> values)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (values.Count == 0) return double.NegativeInfinity;
            double max = Max(values);
            if (double.IsNaN(max) || double.IsInfinity(max)) return max;
            int count = 0;
            double sum = 0d, compensation = 0d;
            for (int i = 0; i < values.Count; i++)
            {
                double v = values[i];
                if (v == max) count++;
                else NeumaierAdd(ref sum, ref compensation, Math.Exp(v - max));
            }
            return LogSumExpResult(max, count, NeumaierResult(sum, compensation));
        }

        /// <summary>
        /// Returns the log-sum-exponential of each column or each row of a matrix of log values.
        /// </summary>
        /// <param name="values">The log values.</param>
        /// <param name="axis">0 to reduce over the rows, giving one value per column; 1 to reduce over the columns,
        /// giving one value per row, as the <c>axis</c> argument of <c>scipy.special.logsumexp</c>.</param>
        /// <returns>One value per column (axis 0) or per row (axis 1), each bitwise
        /// <see cref="LogSumExp(IList{double})"/> of that column or row; −∞ for every entry when the reduced dimension
        /// is empty.</returns>
        /// <remarks>
        /// The matrix is read in row-major order for both axes. For axis 0 the column maxima, counts and compensated
        /// sums are carried across each row, so every column is reduced in the same order and with the same arithmetic as
        /// <see cref="LogSumExp(IList{double})"/> without being copied out of the matrix.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="axis"/> is not 0 or 1.</exception>
        public static double[] LogSumExp(double[,] values, int axis)
        {
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (axis != 0 && axis != 1) throw new ArgumentOutOfRangeException(nameof(axis), "The axis must be 0 or 1.");
            int rows = values.GetLength(0), cols = values.GetLength(1);
            if (axis == 1)
            {
                var rowResult = new double[rows];
                for (int r = 0; r < rows; r++)
                {
                    double max = double.NegativeInfinity;
                    for (int c = 0; c < cols; c++)
                    {
                        double v = values[r, c];
                        if (double.IsNaN(v)) { max = double.NaN; break; }
                        if (v > max) max = v;
                    }
                    if (double.IsNaN(max) || double.IsInfinity(max))
                    {
                        rowResult[r] = max;
                        continue;
                    }
                    int count = 0;
                    double sum = 0d, compensation = 0d;
                    for (int c = 0; c < cols; c++)
                    {
                        double v = values[r, c];
                        if (v == max) count++;
                        else NeumaierAdd(ref sum, ref compensation, Math.Exp(v - max));
                    }
                    rowResult[r] = LogSumExpResult(max, count, NeumaierResult(sum, compensation));
                }
                return rowResult;
            }

            var maxima = new double[cols];
            for (int c = 0; c < cols; c++) maxima[c] = double.NegativeInfinity;
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    if (double.IsNaN(maxima[c])) continue;
                    double v = values[r, c];
                    if (double.IsNaN(v)) maxima[c] = double.NaN;
                    else if (v > maxima[c]) maxima[c] = v;
                }
            }
            var counts = new int[cols];
            var sums = new double[cols];
            var compensations = new double[cols];
            for (int r = 0; r < rows; r++)
            {
                for (int c = 0; c < cols; c++)
                {
                    double max = maxima[c];
                    if (double.IsNaN(max) || double.IsInfinity(max)) continue;
                    double v = values[r, c];
                    if (v == max) counts[c]++;
                    else NeumaierAdd(ref sums[c], ref compensations[c], Math.Exp(v - max));
                }
            }
            var result = new double[cols];
            for (int c = 0; c < cols; c++)
            {
                double max = maxima[c];
                result[c] = double.IsNaN(max) || double.IsInfinity(max)
                    ? max
                    : LogSumExpResult(max, counts[c], NeumaierResult(sums[c], compensations[c]));
            }
            return result;
        }

        /// <summary>
        /// Returns log1p(s / k) + log(k) + max, the log-sum-exponential with its k largest terms separated.
        /// </summary>
        /// <param name="max">The finite largest log value.</param>
        /// <param name="count">The number k of values equal to <paramref name="max"/>, at least one.</param>
        /// <param name="rest">The sum s of exp(x − max) over the other values.</param>
        /// <returns>log(k + s) + max, evaluated in the order of <c>scipy.special.logsumexp</c>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double LogSumExpResult(double max, int count, double rest)
        {
            if (count > 1) rest /= count;
            return Log1p(rest) + Math.Log(count) + max;
        }

        /// <summary>
        /// Returns a sequence of integers.
        /// </summary>
        /// <param name="start">The starting value.</param>
        /// <param name="end">The ending value.</param>
        /// <param name="step">The step size.</param>
        public static int[] Sequence(int start, int end, int step = 1)
        {
            var result = new List<int>();
            if (step == 0) return result.ToArray();
            if (step > 0)
            {
                for (int v = start; v <= end; v += step) result.Add(v);
            }
            else
            {
                for (int v = start; v >= end; v += step) result.Add(v);
            }
            return result.ToArray();
        }

        /// <summary>
        /// Returns a sequence of doubles.
        /// </summary>
        /// <param name="start">The starting value.</param>
        /// <param name="end">The ending value.</param>
        /// <param name="step">The step size.</param>
        public static double[] Sequence(double start, double end, double step = 1)
        {
            var result = new List<double>();
            if (step == 0) return result.ToArray();
            if (step > 0)
            {
                for (double v = start; v <= end; v += step) result.Add(v);
            }
            else
            {
                for (double v = start; v >= end; v += step) result.Add(v);
            }
            return result.ToArray();
        }


        /// <summary>
        /// Returns a compressed a byte array.
        /// </summary>
        /// <param name="data">An array of bytes.</param>
        public static byte[]? Compress(byte[]? data)
        {
            if (data is null) return null!;
            var output = new MemoryStream();
            using (var dstream = new DeflateStream(output, CompressionLevel.Optimal))
            {
                dstream.Write(data, 0, data.Length);
            }
            return output.ToArray();
        }

        /// <summary>
        /// Returns a decompressed byte array.
        /// </summary>
        /// <param name="data">An array of bytes.</param>
        public static byte[]? Decompress(byte[]? data)
        {
            if (data is null) return null!;
            var input = new MemoryStream(data);
            var output = new MemoryStream();
            using (var dstream = new DeflateStream(input, CompressionMode.Decompress))
            {
                dstream.CopyTo(output);
            }
            return output.ToArray();
        }

    }
}
