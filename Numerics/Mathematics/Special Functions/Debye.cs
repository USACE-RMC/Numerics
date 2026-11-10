using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.SpecialFunctions
{
    /// <summary>
    /// The Debye function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// In mathematics, the family of Debye functions is given by the equation:
    /// </para>
    /// <code>
    ///                     x
    ///     D_n(x) = n/x^n ∫  t^n / (e^t - 1) dt
    ///                     0
    /// </code>
    /// <para>
    /// The order n is fixed by each method on this class rather than passed as an argument.
    /// <see cref="Function(double)"/> is the <b>order-3</b> Debye function D₃, and
    /// <see cref="FunctionOrderOne(double)"/> is the <b>order-1</b> Debye function D₁. The two are
    /// different functions and are not interchangeable: D₃(1) = 0.6744156 while D₁(1) = 0.7775046.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Abramowitz, M. and Stegun, I.A. (1972). "Handbook of Mathematical Functions." National Bureau of Standards,
    /// Applied Mathematics Series 55. Section 27.1.
    /// </description></item>
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Debye_function" />
    /// </description></item>
    /// </list>
    /// </remarks>
    public class Debye
    {
        /// <summary>
        /// The Maclaurin coefficients of the order-3 Debye function for the even powers x², x⁴, ... x⁵².
        /// </summary>
        /// <remarks>
        /// The k-th entry is 3 B(2k) / ((2k + 3) (2k)!), where B(2k) is the 2k-th Bernoulli number. The entries
        /// were evaluated with mpmath at 50 digits and rounded once to double.
        /// </remarks>
        private static readonly double[] Debye3SeriesCoefficients =
        {
            0.05, -0.0005952380952380953, 1.1022927689594357E-05,
            -2.2546897546897547E-07, 4.8177131510464845E-09, -1.0568380277374986E-10,
            2.3616240936502374E-12, -5.352126783667236E-14, 1.2265802937539779E-15,
            -2.8367852589887764E-17, 6.610803394032275E-19, -1.5504960762013913E-20,
            3.656593489271863E-22, -8.664694284229884E-24, 2.061774956670621E-25,
            -4.924106287604745E-27, 1.1798695748228635E-28, -2.8353807235887003E-30,
            6.831756773484179E-32, -1.6500156388609047E-33, 3.99378117498809E-35,
            -9.68588224128696E-37, 2.3533213012984824E-38, -5.727266513321583E-40,
            1.3959888929437818E-41, -3.4074964046294964E-43
        };

        /// <summary>
        /// The limit of the order-3 Debye integral, ∫[0 to ∞] t³ / (e^t - 1) dt = π⁴/15.
        /// </summary>
        private const double Debye3Limit = 6.493939402266829;

        /// <summary>
        /// Computes the order-3 Debye function D₃(x).
        /// </summary>
        /// <param name="x">The point to evaluate. Must be non-negative.</param>
        /// <returns>
        /// The order-3 Debye function evaluated at the given x, with D₃(0) = 1 and D₃(+∞) = 0.
        /// </returns>
        /// <remarks>
        /// <para>
        /// D₃(x) = (3/x³) ∫[0 to x] t³ / (e^t - 1) dt. Two branches are used, the expansions of Abramowitz and
        /// Stegun (1972), 27.1.1 and 27.1.2. For 0 &lt; x ≤ 3 the Maclaurin series
        /// D₃(x) = 1 - 3x/8 + Σ[k ≥ 1] 3 B(2k) x^(2k) / ((2k + 3) (2k)!) is summed by Horner's rule through x⁵²,
        /// where the next term is below 1E-18 of the value; the series converges for x &lt; 2π. For x &gt; 3 the
        /// integral is written against its limit π⁴/15, giving
        /// D₃(x) = (3/x³) (π⁴/15 - Σ[k ≥ 1] e^(-kx) (x³/k + 3x²/k² + 6x/k³ + 6/k⁴)), summed until a term falls
        /// below 1E-18 of π⁴/15; from x = 60 the sum is below double precision and is skipped. The divisions by x
        /// are applied one at a time, so x³ never overflows before the result underflows.
        /// </para>
        /// <para>
        /// <b> Accuracy. </b> Measured against mpmath at 50 digits over x in [1E-6, 1000], the worst relative
        /// error is about 1E-15 (see docs/mathematics/oracles/special-functions).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when x &lt; 0.</exception>
        public static double Function(double x)
        {
            if (x < 0.0) { throw new ArgumentOutOfRangeException(nameof(x), "X must be positive."); }
            if (double.IsNaN(x)) return double.NaN;
            if (x == 0.0) return 1.0;

            if (x <= 3.0)
            {
                double z = x * x;
                double sum = 0.0;
                for (int k = Debye3SeriesCoefficients.Length - 1; k >= 0; k--)
                    sum = sum * z + Debye3SeriesCoefficients[k];
                return 1.0 - 0.375 * x + sum * z;
            }

            double remainder = 0.0;
            if (x < 60.0)
            {
                double e = Math.Exp(-x);
                double ek = 1.0;
                for (int k = 1; k <= DebyeMaximumTerms; k++)
                {
                    ek *= e;
                    double kd = k;
                    double term = ek * (((x / kd + 3.0 / (kd * kd)) * x + 6.0 / (kd * kd * kd)) * x
                        + 6.0 / (kd * kd * kd * kd));
                    remainder += term;
                    if (term < 1E-18 * Debye3Limit) break;
                }
            }
            return 3.0 * (Debye3Limit - remainder) / x / x / x;
        }

        /// <summary>
        /// The Maclaurin coefficients of the order-1 Debye function for the even powers x², x⁴, ... x²⁴.
        /// </summary>
        /// <remarks>
        /// The k-th entry is B(2k) / ((2k + 1) * (2k)!), where B(2k) is the 2k-th Bernoulli number. The entries
        /// were evaluated as exact rationals and rounded once to double.
        /// </remarks>
        private static readonly double[] DebyeSeriesCoefficients =
        {
            2.7777777777777776E-02, -2.7777777777777778E-04, 4.7241118669690098E-06,
            -9.1857730746619641E-08, 1.8978869988971000E-09, -4.0647616451442256E-11,
            8.9216910204564523E-13, -1.9939295860721074E-14, 4.5189800296199183E-16,
            -1.0356517612181247E-17, 2.3952186210261870E-19, -5.5817858743250090E-21
        };

        /// <summary>
        /// The largest number of exponential terms summed in the large-argument branch of the order-1 Debye function.
        /// </summary>
        private const int DebyeMaximumTerms = 1000;

        /// <summary>
        /// The size at which an exponential term is small enough to end the large-argument Debye summation.
        /// </summary>
        private const double DebyeTermTolerance = 1E-20;

        /// <summary>
        /// Computes the order-1 Debye function D₁(x) for any real argument.
        /// </summary>
        /// <param name="x">The point to evaluate. Every real value is admissible.</param>
        /// <returns>The order-1 Debye function evaluated at the given x.</returns>
        /// <remarks>
        /// <para>
        /// The order-1 Debye function is D₁(x) = (1/x) ∫[0 to x] t / (e^t - 1) dt. The integrand tends to 1 as
        /// t tends to 0, so the removable limit D₁(0) = 1 is returned exactly. This is a different function from
        /// <see cref="Function(double)"/>, which is the order-3 Debye function; the two are not interchangeable.
        /// </para>
        /// <para>
        /// Three branches are used. A negative argument is reduced by the reflection D₁(-y) = D₁(y) + y/2 for
        /// y &gt; 0, so the function is finite and smooth on the whole real line. For 0 &lt; x ≤ 1 the Maclaurin
        /// series D₁(x) = 1 - x/4 + Σ[k ≥ 1] B(2k) x^(2k) / ((2k + 1) (2k)!) is used, truncated after x²⁴, where
        /// the next term is below 1E-22. For x &gt; 1 the integral is written against its limit π²/6, giving
        /// D₁(x) = π²/(6x) - Σ[k ≥ 1] e^(-k x) (1/k + 1/(k² x)), which converges geometrically and is summed
        /// until the term falls below 1E-20.
        /// </para>
        /// <para>
        /// <b> Accuracy. </b> Measured against mpmath at 60 decimal digits, the worst relative error over
        /// x in [1E-8, 100] and its negative mirror is 3 ulp, or 3.4E-16 absolute, at x near 1.75. At the
        /// endpoints D₁(100) = 0.016449340668482266 matches the asymptote π²/600 = 0.016449340668482264,
        /// and the reflection D₁(-1) - D₁(1) = 0.5 holds exactly.
        /// </para>
        /// </remarks>
        public static double FunctionOrderOne(double x)
        {
            if (x == 0d) return 1d;

            // D₁(-y) = D₁(y) + y/2 for y > 0.
            if (x < 0d) return FunctionOrderOne(-x) - 0.5d * x;

            if (x <= 1d)
            {
                double squared = x * x;
                double power = 1d;
                double sum = 0d;
                for (int k = 0; k < DebyeSeriesCoefficients.Length; k++)
                {
                    power *= squared;
                    sum += DebyeSeriesCoefficients[k] * power;
                }
                return 1d - 0.25d * x + sum;
            }

            double remainder = 0d;
            for (int k = 1; k <= DebyeMaximumTerms; k++)
            {
                double term = Math.Exp(-k * x) * (1d / k + 1d / ((double)k * k * x));
                remainder += term;
                if (term < DebyeTermTolerance) break;
            }
            return Math.PI * Math.PI / (6d * x) - remainder;
        }

        #region List Overloads

        /// <summary>
        /// Computes the order-3 Debye function for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Function(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Function(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Function(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Function(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the order-1 Debye function for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="FunctionOrderOne(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="FunctionOrderOne(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] FunctionOrderOne(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = FunctionOrderOne(x[i]);
            return result;
        }

        #endregion

    }
}
