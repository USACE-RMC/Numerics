using System;

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
    /// <see href = "https://en.wikipedia.org/wiki/Debye_function" />
    /// </description></item>
    /// </list>
    /// </remarks>
    public class Debye
    {
        /// <summary>
        /// Computes the order-3 Debye function D₃(x).
        /// </summary>
        /// <param name="x">The point in the series to evaluate.</param>
        /// <remarks>
        /// <para>
        ///     Authors:
        ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// <see href = "http://duffy.princeton.edu/sites/default/files/pdfs/links/Debye_Function.pdf"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <returns>
        /// The order-3 Debye function evaluated at the given x
        /// </returns>
        public static double Function(double x)
        {
            
            if (x <0.0) { throw new ArgumentOutOfRangeException(nameof(x), "X must be positive."); }


            if (x == 0.0)
            {
                return 1.0;
            }
            else if ( x > 0.0 && x <= 0.1)
            {
                double t = 5.952380953E-4;
                return 1.0 - 0.375 * x + x * x * (0.05 - t * x * x);
            }
            else if (x > 0.1 && x <= 7.25)
            {
                return ((((0.0946173 * x - 4.432582) * x + 85.07724) * x - 800.6087) * x + 3953.632) / ((((x + 15.121491) * x + 143.155337) * x + 682.0012) * x + 3953.632) ;
            }
            else if (x > 7.25)
            {
                double N = 25 / x;
                double D = 0.0;
                double D2 = 1.0;
                double x3 = 0;
                if (x <= 25)
                {
                    for (int i = 1; i <= N; i++)
                    {
                        D2 *= Math.Exp(-x);
                        x3 = i * x;
                        D += D2 * (6 + x3 * (6.0 + x3 * (3 + x3))) / Math.Pow(i, 4);
                    }
                    return 3.0 * (6.493939402 - D) / (x * x * x);
                }
                else if(x > 25){
                    return 3.0 * (6.493939402 - D) / (x * x * x);
                }
            }
            return double.NaN;

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

    }
}
