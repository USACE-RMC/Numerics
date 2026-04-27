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
    /// In mathematics, the Debye function is given by the equation:
    /// </para>
    /// <code>
    ///                  x
    ///     D(x) = x/x^n ∫  t^n / (e^t - 1) dt
    ///                  0
    /// </code>
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
        /// Computes the Debye function.
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
        /// The Debye function evaluated at the given x
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

    }
}
