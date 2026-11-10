using System;
using System.Collections.Generic;
using Numerics.Distributions;

namespace Numerics.Mathematics.SpecialFunctions
{
    /// <summary>
    /// The error function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The error function, also known as the Gauss error function, is defined 
    /// as:
    /// </para>
    /// <code>
    ///                        X
    ///     erf(X) = 2/sqrt(𝜋) ∫  e^(-t^2) dt
    ///                        0           
    /// </code>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// <see href = "https://en.wikipedia.org/wiki/Error_function#cite_note-2"/>
    ///  </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public partial class Erf
    {
        /// <summary>
        /// 2/sqrt(π), the slope of erf at the origin.
        /// </summary>
        private const double TwoOverSqrtPi = 1.1283791670955126d;

        /// <summary>
        /// Computes the error function
        /// </summary>
        /// <remarks>
        /// For |X| &lt; 1E-8 the leading Taylor term 2X/sqrt(π) is returned (the next term is below
        /// 4E-17 relative), so tiny arguments do not underflow to zero through X².
        /// </remarks>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The error function evaluated with the given upper bound
        /// </returns>
        public static double Function(double X)
        {
            if (Math.Abs(X) < 1E-8)
            {
                return TwoOverSqrtPi * X;
            }

            if (X < 0d)
            {
                return -Gamma.LowerIncomplete(0.5d, X * X);
            }
            else
            {
                return Gamma.LowerIncomplete(0.5d, X * X);
            }
        }

        /// <summary>
        /// Computes the complement of the error function
        /// </summary>
        /// <remarks>
        /// <para>
        /// For x &gt;= 0, erfc(x) = exp(-x²) erfcx(x), where the scaled complement erfcx(x) = exp(x²) erfc(x) is
        /// summed from the Chebyshev expansion of Schonfelder (1978) in t = (8x - 30) / (4x + 15) (see
        /// <see cref="ScaledErfc(double)"/>), and exp(-x²) is formed by <see cref="ExpNegativeScaledSquare(double, double)"/>
        /// without the rounding of x². For x &lt; 0, erfc(x) = 2 - erfc(-x). The argument is used as given, and the
        /// relative error measured against mpmath is at most about 5E-16 (root mean square 2E-16) until the result
        /// becomes subnormal near x = 26.5; erfc(x) underflows to zero above about 27.3, erfc(-∞) = 2 and NaN returns NaN. Below
        /// |x| = 2^-27 the two-term Taylor series 1 - 2x / sqrt(π) is returned, so erfc(0) = 1 exactly.
        /// </para>
        /// <para>
        /// References: Schonfelder, J. L. (1978). Chebyshev expansions for the error and related functions.
        /// Mathematics of Computation 32(144), 1232-1240. Cody, W. J. (1969). Rational Chebyshev approximations for
        /// the error function. Mathematics of Computation 23(107), 631-637 (the splitting of x² in the exponential).
        /// </para>
        /// </remarks>
        /// <param name="X"> The upper bound </param>>
        /// <returns>
        /// The complement of the error function evaluated with the given upper bound
        /// </returns>
        public static double Erfc(double X)
        {
            if (double.IsNaN(X))
            {
                return double.NaN;
            }

            double ax = Math.Abs(X);
            double tail;
            if (ax < TaylorLimit)
                tail = 1d - TwoOverSqrtPi * ax;
            else if (ax >= ErfcUnderflow)
                tail = 0d;
            else
                tail = ExpNegativeScaledSquare(ax, 1d) * ScaledErfc(ax);
            return X < 0d ? 2d - tail : tail;
        }

        /// <summary>
        /// 2^-27, below which erfc(x) = 1 - 2x / sqrt(π) to double precision (the next term, 2x³ / (3 sqrt(π)), is below 1E-24).
        /// </summary>
        internal const double TaylorLimit = 7.450580596923828E-09;

        /// <summary>
        /// The argument above which erfc underflows to zero: exp(-x²) is below the smallest subnormal double.
        /// </summary>
        private const double ErfcUnderflow = 27.3d;

        /// <summary>
        /// Returns the scaled complementary error function erfcx(x) = exp(x²) erfc(x) for finite x &gt;= 0.
        /// </summary>
        /// <param name="x">The finite non-negative argument.</param>
        /// <returns>erfcx(x), which decreases from 1 at x = 0 like 1 / (x sqrt(π)).</returns>
        /// <remarks>
        /// For x &lt; 3 the Chebyshev series of Schonfelder (1978) in t = (8x - 30) / (4x + 15), which maps [0, ∞)
        /// onto [-2, 2], is summed by the Clenshaw recurrence over its first 27 coefficients. The series sums to a
        /// small value from coefficients near one as t approaches 2, so for x &gt;= 3 the value is instead
        /// g(1 / x²) / (x sqrt(π)) with g from <see cref="ScaledErfcAsymptotic(double)"/>. Measured against mpmath,
        /// the two pieces are within 5E-16 and 7E-17 relative in double precision
        /// (docs/mathematics/oracles/special-functions/erfc_coefficients.py).
        /// </remarks>
        internal static double ScaledErfc(double x)
        {
            if (x >= ScaledErfcSplit)
            {
                return ScaledErfcAsymptotic(1d / (x * x)) / (x * SqrtPi);
            }

            double t = (8d * x - 30d) / (4d * x + 15d);
            return Evaluate.Chebyshev(ScaledErfcCoefficients, t, ScaledErfcTerms);
        }

        /// <summary>
        /// Returns g(s) = sqrt(π) x erfcx(x) at s = 1 / x², for x &gt;= 3 (0 &lt;= s &lt;= 1/9).
        /// </summary>
        /// <param name="s">The inverse square 1 / x², at most about 1/9.</param>
        /// <returns>g(s), which is 1 - s/2 + 3s²/4 - ... and lies in [0.95, 1].</returns>
        /// <remarks>
        /// The polynomial of degree 21 is the near-minimax Chebyshev fit of mpmath.chebyfit on [0, 1/9],
        /// within 2E-20 of g (docs/mathematics/oracles/special-functions/erfc_coefficients.py). Its leading
        /// coefficients are those of the asymptotic series of erfc, and Horner's rule evaluates it to about
        /// 6E-17 relative because the terms decrease rapidly.
        /// </remarks>
        internal static double ScaledErfcAsymptotic(double s)
        {
            double[] c = ScaledErfcAsymptoticCoefficients;
            double v = c[0];
            for (int i = 1; i < c.Length; i++)
            {
                v = v * s + c[i];
            }
            return v;
        }

        /// <summary>
        /// The argument from which <see cref="ScaledErfc(double)"/> uses the polynomial in 1 / x².
        /// </summary>
        internal const double ScaledErfcSplit = 3d;

        /// <summary>
        /// sqrt(π).
        /// </summary>
        internal const double SqrtPi = 1.7724538509055160273d;

        /// <summary>
        /// Coefficients of g(s) = sqrt(π) x erfcx(x) in powers of s = 1 / x², highest degree first.
        /// </summary>
        private static readonly double[] ScaledErfcAsymptoticCoefficients =
        {
            -12407552698613.125d,
            16137455100736.809d,
            -9905325518072.896d,
            3825478120604.5835d,
            -1046965049244.0111d,
            217275927511.00812d,
            -35847526752.09136d,
            4893068876.67497d,
            -575171732.2155399d,
            60945642.98519396d,
            -6139944.585513929d,
            622665.1497328572d,
            -66909.49745838967d,
            7910.721267393985d,
            -1055.6366634952337d,
            162.42074314166672d,
            -29.531241309060057d,
            6.562499954955915d,
            -1.874999999855718d,
            0.7499999999997532d,
            -0.49999999999999983d,
            1.0d
        };

        /// <summary>
        /// Returns exp(-s x²) for 0 &lt;= x &lt; 64 and s = 1 or 0.5, without the rounding error of x².
        /// </summary>
        /// <param name="x">The non-negative argument, below 64.</param>
        /// <param name="scale">The exact factor s, 1 or 0.5.</param>
        /// <returns>exp(-s x²).</returns>
        /// <remarks>
        /// A rounded x² carries a relative error of up to 2^-53, which the exponential multiplies by s x², so
        /// exp(-x²) would lose about x² units in the last place. As in Cody (1969), x is split as x = h + l, with
        /// h = x truncated to 20 fraction bits; below 64, h² then has at most 52 significant bits and is exact, and
        /// exp(-s x²) = exp(-s h²) exp(-s (x - h)(x + h)), where the second factor, within 1.3E-4 of one, is
        /// summed from its Taylor series.
        /// </remarks>
        internal static double ExpNegativeScaledSquare(double x, double scale)
        {
            const double Split = 1048576d;
            double h = Math.Truncate(x * Split) / Split;
            // 0 <= x - h < 2^-20 and x + h < 128, so |v| < 1.3E-4 and the cubic Taylor polynomial of exp(v) is within v⁴/24 < 1E-17
            double v = -scale * ((x - h) * (x + h));
            return Math.Exp(-scale * (h * h)) * (1d + v * (1d + v * (0.5d + v / 6d)));
        }

        /// <summary>
        /// The number of Chebyshev coefficients of <see cref="ScaledErfcCoefficients"/> summed by <see cref="ScaledErfc(double)"/>.
        /// </summary>
        private const int ScaledErfcTerms = 27;

        /// <summary>
        /// Chebyshev coefficients of Schonfelder (1978) for erfcx(x) = exp(x²) erfc(x) = (b_0 - b_2) / 2, the Clenshaw
        /// sums in t = (8x - 30) / (4x + 15).
        /// </summary>
        /// <remarks>
        /// The 44 coefficients are those of Schonfelder (1978) as used by Genz's MVNPHI; the multivariate normal
        /// integrands sum the first 25. docs/mathematics/oracles/special-functions/erfc_coefficients.py
        /// recomputes them with mpmath (within 5E-31) and checks the truncated sums against erfcx.
        /// </remarks>
        internal static readonly double[] ScaledErfcCoefficients =
        {
            6.10143081923200417926465815756E-1,
            -4.34841272712577471828182820888E-1,
            1.76351193643605501125840298123E-1,
            -6.0710795609249414860051215825E-2,
            1.7712068995694114486147141191E-2,
            -4.321119385567293818599864968E-3,
            8.54216676887098678819832055E-4,
            -1.27155090609162742628893940E-4,
            1.1248167243671189468847072E-5,
            3.13063885421820972630152E-7,
            -2.70988068537762022009086E-7,
            3.0737622701407688440959E-8,
            2.515620384817622937314E-9,
            -1.028929921320319127590E-9,
            2.9944052119949939363E-11,
            2.6051789687266936290E-11,
            -2.634839924171969386E-12,
            -6.43404509890636443E-13,
            1.12457401801663447E-13,
            1.7281533389986098E-14,
            -4.264101694942375E-15,
            -5.45371977880191E-16,
            1.58697607761671E-16,
            2.0899837844334E-17,
            -5.900526869409E-18,
            -9.41893387554E-19,
            2.14977356470E-19,
            4.6660985008E-20,
            -7.243011862E-21,
            -2.387966824E-21,
            1.91177535E-22,
            1.20482568E-22,
            -6.72377E-25,
            -5.747997E-24,
            -4.28493E-25,
            2.44856E-25,
            4.3793E-26,
            -8.151E-27,
            -3.089E-27,
            9.3E-29,
            1.74E-28,
            1.6E-29,
            -8.0E-30,
            -2.0E-30
        };

        /// <summary>
        /// Computes the inverse error function
        /// </summary>
        /// <remarks>
        /// The function is evaluated as an odd function of y. For |y| &lt; 0.5 a normal-quantile start
        /// is refined with two Halley steps on erf(x) - |y|, which keeps full relative precision near
        /// zero. For |y| &gt;= 0.5 the result is <see cref="InverseErfc(double)"/>(1 - |y|), where 1 - |y| is exact.
        /// y = ±1 returns ±∞.
        /// </remarks>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse error function evaluated at the given y
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when |y| &gt; 1.</exception>
        public static double InverseErf(double y)
        {
            if (double.IsNaN(y))
            {
                return double.NaN;
            }

            double ay = Math.Abs(y);
            double r;
            if (ay == 1.0d)
            {
                r = double.PositiveInfinity;
            }
            else if (ay >= 0.5d)
            {
                r = InverseErfc(1.0d - ay);
            }
            else
            {
                r = Normal.StandardZ(0.5d * ay + 0.5d) * Tools.Sqrt2 / 2.0d;
                for (int i = 0; i < 2; i++)
                {
                    // Halley step: f = erf(r) - |y|, f' = 2/sqrt(π) exp(-r²), f''/f' = -2r
                    double u = (Function(r) - ay) / (TwoOverSqrtPi * Math.Exp(-r * r));
                    r -= u / (1.0d + r * u);
                }
            }

            return y < 0d ? -r : r;
        }

        /// <summary>
        /// Computes the inverse of the complement of the error function
        /// </summary>
        ///<param name="y"> The value to be evaluated (such that y = erf(erf^-1(y)) ) </param>
        /// <returns>
        /// The inverse of the complement of the error function evaluated at the given y
        /// </returns>
        public static double InverseErfc(double y)
        {
            // The identity erfc(x) = 2*Phi(-x*sqrt(2)) inverts through the lower tail, so a small
            // argument maps to StandardZ of y/2 rather than to a probability that rounds to one
            return -Normal.StandardZ(0.5d * y) / Tools.Sqrt2;
        }

        #region List Overloads

        /// <summary>
        /// Computes the error function for each value of a list.
        /// </summary>
        /// <param name="X">The values of X.</param>
        /// <returns>One result per value, each equal to <see cref="Function(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Function(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Function(IList<double> X)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            var result = new double[X.Count];
            for (int i = 0; i < X.Count; i++)
                result[i] = Function(X[i]);
            return result;
        }

        /// <summary>
        /// Computes the complementary error function for each value of a list.
        /// </summary>
        /// <param name="X">The values of X.</param>
        /// <returns>One result per value, each equal to <see cref="Erfc(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Erfc(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Erfc(IList<double> X)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            var result = new double[X.Count];
            for (int i = 0; i < X.Count; i++)
                result[i] = Erfc(X[i]);
            return result;
        }

        /// <summary>
        /// Computes the inverse error function for each value of a list.
        /// </summary>
        /// <param name="y">The values of y.</param>
        /// <returns>One result per value, each equal to <see cref="InverseErf(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="InverseErf(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] InverseErf(IList<double> y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            var result = new double[y.Count];
            for (int i = 0; i < y.Count; i++)
                result[i] = InverseErf(y[i]);
            return result;
        }

        /// <summary>
        /// Computes the inverse complementary error function for each value of a list.
        /// </summary>
        /// <param name="y">The values of y.</param>
        /// <returns>One result per value, each equal to <see cref="InverseErfc(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="InverseErfc(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] InverseErfc(IList<double> y)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            var result = new double[y.Count];
            for (int i = 0; i < y.Count; i++)
                result[i] = InverseErfc(y[i]);
            return result;
        }

        #endregion

    }
}