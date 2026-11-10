// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
// 
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

using System;
using System.Collections.Generic;

using Numerics.Distributions;

namespace Numerics.Mathematics.SpecialFunctions
{

    /// <summary>
    /// Gamma Γ(x) functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description </b>
    ///     In mathematics, the gamma function (represented by the capital Greek
    ///     letter Γ) is an extension of the factorial function, with its argument
    ///     shifted down by 1, to real and complex numbers. That is, if <c>n</c> is
    ///     a positive integer:
    /// </para>
    /// <code>
    /// 
    ///     Γ(n) = (n-1)!
    /// 
    /// </code>
    /// <para>
    ///     The gamma function is defined for all complex numbers except the negative
    ///     integers and zero. For complex numbers with a positive real part, it is
    ///     defined via an improper integral that converges:
    /// </para>
    /// <code>
    ///            ∞
    ///     Γ(z) = ∫  t^(z-1)e^(-t) dt
    ///            0
    /// </code>
    /// <para>
    ///     This integral function is extended by analytic continuation to all
    ///     complex numbers except the non-positive integers (where the function
    ///     has simple poles), yielding the meromorphic function called the gamma
    ///     function.
    ///     The gamma function is a component in various probability-distribution
    ///     functions, and as such it is applicable in the fields of probability
    ///     and statistics, as well as combinatorics.
    /// </para>
    /// <para>
    /// <b> References: </b>
    ///     This code was copied and modified from two primary sources: 1) THe LMOMENTS FORTAN package from J. R. M. Hosking;
    ///     and 2) the Accord Math Library.
    /// <list type="bullet">
    /// <item><description>
    ///     LMOMENTS package. Available at: <see href = "http://ftp.uni-bayreuth.de/math/statlib/general/lmoments"/>
    /// </description></item>
    /// <item><description>
    ///     Accord Math Library, <see href="http://accord-framework.net"/>
    /// </description></item>
    /// <item><description>
    ///     Wikipedia contributors, "Gamma function,". Wikipedia, The Free
    ///     Encyclopedia. Available at: <see href = "http://en.wikipedia.org/wiki/Gamma_function"/>
    /// </description></item>
    /// <item><description>
    ///     Cephes Math Library, <see href = "http://www.netlib.org/cephes/"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public sealed class Gamma
    {
        private static readonly double[] _p = [0.99999999999980993d, 676.5203681218851d, -1259.1392167224028d, 771.32342877765313d, -176.61502916214059d, 12.507343278686905d, -0.13857109526572012d, 0.0000099843695780195716d, 0.00000015056327351493116d];
        
        /// <summary>
        /// The maximum Gamma value. 
        /// </summary>
        public const double GammaMax = 171.624376956302725;

        private static readonly double[] gamma_P =
        {
            1.60119522476751861407E-4,
            1.19135147006586384913E-3,
            1.04213797561761569935E-2,
            4.76367800457137231464E-2,
            2.07448227648435975150E-1,
            4.94214826801497100753E-1,
            9.99999999999999996796E-1
        };

        private static readonly double[] gamma_Q =
        {
            -2.31581873324120129819E-5,
            5.39605580493303397842E-4,
            -4.45641913851797240494E-3,
            1.18139785222060435552E-2,
            3.58236398605498653373E-2,
            -2.34591795718243348568E-1,
            7.14304917030273074085E-2,
            1.00000000000000000320E0
        };

        private static readonly double[] STIR =
        {
            7.87311395793093628397E-4,
            -2.29549961613378126380E-4,
            -2.68132617805781232825E-3,
            3.47222221605458667310E-3,
            8.33333333333482257126E-2,
        };


        // Rational minimax approximations to ln Γ(x) of Cody and Hillstrom (1967), with the
        // coefficients of the SPECFUN routine DLGAMA (Cody, 1993). Numerator coefficients are listed
        // from the highest power down; the denominators of P1/Q1 and P2/Q2 have an implicit leading x^8.

        // ln Γ(x) = (x - 1)(D1 + (x - 1) P1(x - 1) / Q1(x - 1)) on (0.5, 1.5), where D1 = -γ.
        private const double LogGammaD1 = -5.772156649015328605195174E-1;

        private static readonly double[] LogGammaP1 =
        {
            4.945235359296727046734888E0, 2.018112620856775083915565E2,
            2.290838373831346393026739E3, 1.131967205903380828685045E4,
            2.855724635671635335736389E4, 3.848496228443793359990269E4,
            2.637748787624195437963534E4, 7.225813979700288197698961E3
        };

        private static readonly double[] LogGammaQ1 =
        {
            6.748212550303777196073036E1, 1.113332393857199323513008E3,
            7.738757056935398733233834E3, 2.763987074403340708898585E4,
            5.499310206226157329794414E4, 6.161122180066002127833352E4,
            3.635127591501940507276287E4, 8.785536302431013170870835E3
        };

        // ln Γ(x) = (x - 2)(D2 + (x - 2) P2(x - 2) / Q2(x - 2)) on (1.5, 4), where D2 = 1 - γ.
        private const double LogGammaD2 = 4.227843350984671393993777E-1;

        private static readonly double[] LogGammaP2 =
        {
            4.974607845568932035012064E0, 5.424138599891070494101986E2,
            1.550693864978364947665077E4, 1.847932904445632425417223E5,
            1.088204769468828767498470E6, 3.338152967987029735917223E6,
            5.106661678927352456275255E6, 3.074109054850539556250927E6
        };

        private static readonly double[] LogGammaQ2 =
        {
            1.830328399370592604055942E2, 7.765049321445005871323047E3,
            1.331903827966074194402448E5, 1.136705821321969608938755E6,
            5.267964117437946917577538E6, 1.346701454311101692290052E7,
            1.782736530353274213975932E7, 9.533095591844353613395747E6
        };

        // ln Γ(x) = D4 + (x - 4) P4(x - 4) / Q4(x - 4) on (4, 12), where D4 = ln 6. The leading
        // denominator term is -x^8, stored as the first coefficient of Q4.
        private const double LogGammaD4 = 1.791759469228055000094023E0;

        private static readonly double[] LogGammaP4 =
        {
            1.474502166059939948905062E4, 2.426813369486704502836312E6,
            1.214755574045093227939592E8, 2.663432449630976949898078E9,
            2.940378956634553899906876E10, 1.702665737765398868392998E11,
            4.926125793377430887588120E11, 5.606251856223951465078242E11
        };

        private static readonly double[] LogGammaQ4 =
        {
            -1.0,
            2.690530175870899333379843E3, 6.393885654300092398984238E5,
            4.135599930241388052042842E7, 1.120872109616147941376570E9,
            1.488613728678813811542398E10, 1.016803586272438228077304E11,
            3.417476345507377132798597E11, 4.463158187419713286462081E11
        };

        // Correction to Stirling's formula on (12, ∞): a minimax polynomial in 1/x², listed from the
        // highest power down, whose value divided by x is added to (x - 1/2) ln x - x + ln √(2π).
        private static readonly double[] LogGammaStirling =
        {
            5.7083835261E-03, -1.910444077728E-03, 8.4171387781295E-04,
            -5.952379913043012E-04, 7.93650793500350248E-04,
            -2.777777777777681622553E-03, 8.333333333333333331554247E-02
        };

        /// <summary>
        /// The Stirling approximation
        /// </summary>
        /// <remarks>
        /// <para>
        /// In mathematics, Stirling's approximation (or Stirling's formula) is an approximation for factorials.
        /// It is a good approximation, leading to accurate results even for small values of n.
        /// </para>
        /// References: 
        /// <list type="bullet">
        /// <item><description>
        /// <see href = "https://en.wikipedia.org/wiki/Stirling%27s_approximation"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="x">  The value to be evaluated </param>
        /// <returns>
        /// The Stirling approximation for the given x
        /// </returns>
        public static double Stirling(double x)
        {
            double MAXSTIR = 143.01608;

            double w = 1.0 / x;
            double y = Math.Exp(x);

            w = 1.0 + w * Evaluate.PolynomialRev(STIR, w, 4);

            if (x > MAXSTIR)
            {
                double v = Math.Pow(x, 0.5 * x - 0.25);
                if (double.IsPositiveInfinity(v) && double.IsPositiveInfinity(y))
                {
                    y = double.PositiveInfinity;
                }
                else
                {
                    y = v * (v / y);
                }
            }
            else
            {
                y = Math.Pow(x, x - 0.5) / y;
            }

            y = Tools.Sqrt2PI * y * w;
            return y;
        }

        /// <summary>
        /// The Lancoz approximation
        /// </summary>
        /// <remarks>
        /// <para>
        /// In mathematics, the Lanczos approximation is a method for computing the gamma function numerically,
        /// published by Cornelius Lanczos in 1964. It is a practical alternative to the more popular
        /// Stirling's approximation for calculating the gamma function with fixed precision.
        /// </para>
        ///  References: 
        /// <list type="bullet">
        /// <item><description>
        /// <see herf = "https://en.wikipedia.org/wiki/Lanczos_approximation"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="x"> The value to be evaluated </param>
        /// <returns>
        /// The Lancoz approximation of the Gamma function evaluated at the given x
        /// </returns>
        public static double Lanczos(double x) 
        {
            int g = 7;
            if (x < 0.5d)
            {
                return Math.PI / (Math.Sin(Math.PI * x) * Function(1d - x));
            }
            else
            {
                x -= 1d;
                double y = _p[0];
                for (int i = 1; i < g + 2; i++)
                    y += _p[i] / (x + i);
                double t = x + g + 0.5d;
                return Tools.Sqrt2PI * Math.Pow(t, x + 0.5d) * Math.Exp(-t) * y;
            }
        }

        /// <summary>
        /// The Gamma function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Large positive arguments return positive infinity rather than throwing.
        /// </para>
        /// </remarks>
        /// <param name="x"> The value to be evaluated </param>
        /// <returns>
        /// The Gamma function evaluated at the given x
        /// </returns>
        /// <exception cref="ArithmeticException">Thrown at zero and at negative integers of magnitude 33 or less (poles).</exception>
        /// <exception cref="OverflowException">Thrown at negative integers of magnitude above 33 (poles).</exception>
        public static double Function(double x)
        {
            double p, z;
            double q = Math.Abs(x);

            if (q > 33.0)
            {
                if (x < 0.0)
                {
                    p = Math.Floor(q);

                    if (p == q)
                        throw new OverflowException();

                    z = q - p;
                    if (z > 0.5)
                    {
                        p += 1.0;
                        z = q - p;
                    }
                    z = q * Math.Sin(Math.PI * z);

                    if (z == 0.0)
                        throw new OverflowException();

                    z = Math.Abs(z);
                    z = Math.PI / (z * Stirling(q));

                    return -z;
                }
                else
                {
                    return Stirling(x);
                }
            }

            z = 1.0;
            while (x >= 3.0)
            {
                x -= 1.0;
                z *= x;
            }

            while (x < 0.0)
            {
                if (x == 0.0)
                {
                    throw new ArithmeticException();
                }
                else if (x > -1.0E-9)
                {
                    return (z / ((1.0 + Tools.Euler * x) * x));
                }
                z /= x;
                x += 1.0;
            }

            while (x < 2.0)
            {
                if (x == 0.0)
                {
                    throw new ArithmeticException();
                }
                else if (x < 1.0E-9)
                {
                    return (z / ((1.0 + Tools.Euler * x) * x));
                }

                z /= x;
                x += 1.0;
            }

            if (x == 2.0)
                return z;

            x -= 2.0;
            p = Evaluate.PolynomialRev(gamma_P, x, 6);
            q = Evaluate.PolynomialRev(gamma_Q, x, 7);
            return z * p / q;
        }

        /// <summary>
        /// The digamma function.
        /// </summary>
        /// <remarks>
        /// The digamma function is also know as Euler's PSI function.
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Based algorithm AS103, Applied Statistics, 1976, Vol. 25, No. 3.
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="X"> The value to be evaluated </param>
        /// <returns>
        /// The first derivative of LogGamma given X.
        /// </returns>
        public static double Digamma(double X)
        {
            double SMALL = 0.000000001d;
            double CRIT = 13d;
            double C1 = 0.083333333333333329d;
            double C2 = -0.0083333333333333332d;
            double C3 = 0.003968253968253968d;
            double C4 = -0.0041666666666666666d;
            double C5 = 0.007575757575757576d;
            double C6 = -0.021092796092796094d;
            double C7 = 0.083333333333333329d;
            double D1 = -0.57721566490153287d;
            double _digamma = 0.0d;
            if (X <= 0.0d)
            {
                throw new ArgumentOutOfRangeException("X", "X must be greater than zero.");
            }

            if (X > SMALL)
            {
                double Y = X;
                while (Y < CRIT)
                {
                    _digamma = _digamma - 1.0d / Y;
                    Y = Y + 1.0d;
                }

                _digamma = _digamma + Math.Log(Y) - 0.5d / Y;
                Y = 1.0d / (Y * Y);
                double SUM = ((((((C7 * Y + C6) * Y + C5) * Y + C4) * Y + C3) * Y + C2) * Y + C1) * Y;
                _digamma = _digamma - SUM;
            }
            else
            {
                _digamma = D1 - 1.0d / X;
            }

            return _digamma;
        }

        /// <summary>
        /// The Trigamma function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The argument is shifted with ψ'(x) = ψ'(x + 1) + 1/x² until it reaches 20, and the
        /// asymptotic series ψ'(z) ~ 1/z + 1/(2z²) + Σ B₂ₖ / z^(2k+1) is then summed with the exact
        /// Bernoulli numbers B₂ through B₁₄. The shift terms are added from smallest to largest, so
        /// the result is accurate to a few ulps for every x &gt; 0, including tiny x where ψ'(x) ≈ 1/x² + π²/6.
        /// </para>
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Based algorithm AS121, Applied Statistics, 1978, Vol. 27, No. 1.
        /// </description></item>
        /// <item><description>
        /// This code has been adapted from the FORTRAN77 and subsequent
        /// C code by B. E. Schneider and John Burkardt. The code had been
        /// made public under the GNU LGPL license.
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="X"> The value to be evaluated </param>
        /// <returns>
        /// Calculates Trigamma(x) = d^2 logGamma(x) / dx^2
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when X &lt;= 0.</exception>
        public static double Trigamma(double X)
        {
            const double B = 20.0d;
            // Bernoulli numbers B2 through B14
            const double B2 = 1.0d / 6.0d;
            const double B4 = -1.0d / 30.0d;
            const double B6 = 1.0d / 42.0d;
            const double B8 = -1.0d / 30.0d;
            const double B10 = 5.0d / 66.0d;
            const double B12 = -691.0d / 2730.0d;
            const double B14 = 7.0d / 6.0d;
            if (X <= 0d)
            {
                throw new ArgumentOutOfRangeException("X", "X must be greater than zero.");
            }

            // Number of unit shifts needed to bring the argument to B or above
            int n = X < B ? (int)Math.Ceiling(B - X) : 0;
            double Z = X + n;

            // Asymptotic series at Z >= B
            double Y = 1.0d / (Z * Z);
            double S = Y * (B2 + Y * (B4 + Y * (B6 + Y * (B8 + Y * (B10 + Y * (B12 + Y * B14))))));
            double _trigamma = (1.0d + 0.5d / Z + S) / Z;

            // Recurrence terms, smallest first
            for (int k = n - 1; k >= 0; k--)
            {
                double t = X + k;
                _trigamma += 1.0d / (t * t);
            }

            return _trigamma;
        }

        /// <summary>
        /// Computes the logarithm of the Gamma function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// For z &gt; 0 the value comes from the rational minimax approximations of Cody and Hillstrom (1967),
        /// with the coefficients and interval scheme of the SPECFUN routine DLGAMA (Cody, 1993):
        /// </para>
        /// <list type="bullet">
        /// <item><description>
        /// z ≤ 2.22E-16: ln Γ(z) = -ln z.
        /// </description></item>
        /// <item><description>
        /// (2.22E-16, 1.5]: ln Γ(z) = (z - 1)(D1 + (z - 1) R1(z - 1)), or ln Γ(1 + z) - ln z evaluated with the same
        /// approximation (or the one about 2) when z &lt; 0.6796875. There ln z carries its rounding residual from one
        /// Newton step, ln z = L + (z - e^L) / e^L with L = fl(ln z), so the sum is rounded once.
        /// </description></item>
        /// <item><description>
        /// (1.5, 4]: ln Γ(z) = (z - 2)(D2 + (z - 2) R2(z - 2)).
        /// </description></item>
        /// <item><description>
        /// (4, 12]: ln Γ(z) = ln 6 + (z - 4) R4(z - 4).
        /// </description></item>
        /// <item><description>
        /// z &gt; 12: Stirling's formula (z - 1/2) ln z - z + ln √(2π) plus a minimax correction in 1/z².
        /// </description></item>
        /// </list>
        /// <para>
        /// The approximations about 1 and 2 carry the factor (z - 1) or (z - 2) explicitly, so the relative
        /// accuracy holds at the zeros ln Γ(1) = ln Γ(2) = 0. Against 40-digit references the relative error
        /// is below 1E-15 over (0, 1E300), including integer and half-integer arguments.
        /// </para>
        /// <para>
        /// For z &lt; 0 the reflection formula ln|Γ(z)| = ln π - ln|sin(πz)| - ln Γ(1 - z) is used, with sin(πz)
        /// evaluated from the distance of z to the nearest integer, which is exact in double precision. For
        /// negative non-integer z the result is ln|Γ(z)|, matching scipy.special.gammaln.
        /// </para>
        /// <para>
        /// Zero and the negative integers are poles and return positive infinity, as does z = +∞; arguments
        /// above about 2.55E305 overflow to positive infinity. NaN and z = -∞ return NaN.
        /// </para>
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Cody, W. J., and Hillstrom, K. E. (1967). Chebyshev approximations for the natural logarithm of the
        /// Gamma function. Mathematics of Computation, 21(98), 198-203.
        /// </description></item>
        /// <item><description>
        /// Cody, W. J. (1993). Algorithm 715: SPECFUN, a portable FORTRAN package of special function routines
        /// and test drivers. ACM Transactions on Mathematical Software, 19(1), 22-30.
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="z">The argument of the gamma function.</param>
        /// <returns>The natural logarithm of the absolute value of Γ(z).</returns>
        public static double LogGamma(double z)
        {
            const double LnPi = 1.1447298858494001741434273513530587116472948129153d;

            if (double.IsNaN(z)) return double.NaN;
            if (double.IsPositiveInfinity(z)) return double.PositiveInfinity;

            if (z <= 0.0)
            {
                // Distance to the nearest integer; the subtraction is exact in double precision.
                double d = Math.Abs(z - Math.Round(z));
                if (d == 0.0) return double.PositiveInfinity;
                return LnPi - Math.Log(Math.Sin(Math.PI * d)) - LogGamma(1.0 - z);
            }

            if (z <= 2.0 * Tools.DoubleMachineEpsilon)
            {
                return -Math.Log(z);
            }
            else if (z <= 1.5)
            {
                // Below 0.6796875, ln Γ(z) = ln Γ(1 + z) - ln z. The logarithm carries its rounding residual:
                // with L = fl(ln z), ln z = L + ln(1 + r) ≈ L + r, where r = (z - e^L) / e^L.
                double corr = 0.0;
                double corrResidual = 0.0;
                double xm1 = z - 1.0;
                if (z < 0.6796875)
                {
                    double logZ = Math.Log(z);
                    double expLogZ = Math.Exp(logZ);
                    corr = -logZ;
                    corrResidual = -(z - expLogZ) / expLogZ;
                    xm1 = z;
                }

                if (z <= 0.5 || z >= 0.6796875)
                {
                    return corr + (corrResidual + xm1 * (LogGammaD1 + xm1 * (Evaluate.PolynomialRev(LogGammaP1, xm1) / Evaluate.PolynomialRev_1(LogGammaQ1, xm1))));
                }

                // 0.5 < z < 0.6796875: 1 + z lies in (1.5, 1.68), where the approximation about 2 is used.
                double xm2 = z - 1.0;
                return corr + (corrResidual + xm2 * (LogGammaD2 + xm2 * (Evaluate.PolynomialRev(LogGammaP2, xm2) / Evaluate.PolynomialRev_1(LogGammaQ2, xm2))));
            }
            else if (z <= 4.0)
            {
                double xm2 = z - 2.0;
                return xm2 * (LogGammaD2 + xm2 * (Evaluate.PolynomialRev(LogGammaP2, xm2) / Evaluate.PolynomialRev_1(LogGammaQ2, xm2)));
            }
            else if (z <= 12.0)
            {
                double xm4 = z - 4.0;
                return LogGammaD4 + xm4 * (Evaluate.PolynomialRev(LogGammaP4, xm4) / Evaluate.PolynomialRev(LogGammaQ4, xm4));
            }
            else
            {
                // The correction is below the rounding of the result beyond about 2.25E76.
                double res = 0.0;
                if (z <= 2.25E76)
                {
                    res = Evaluate.PolynomialRev(LogGammaStirling, 1.0 / (z * z));
                }

                double lnz = Math.Log(z);
                res = res / z + Tools.LogSqrt2PI - 0.5 * lnz;
                return res + z * (lnz - 1.0);
            }
        }

        /// <summary>
        /// The incomplete gamma integral.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Returns the regularized lower incomplete gamma ratio P(alpha, X) = γ(alpha, X) / Γ(alpha), the quantity
        /// of algorithm AS 239 (Shea, 1988), with the argument order and validation of that algorithm. The value is
        /// the one of <see cref="LowerIncomplete(double, double)"/>, whose remarks describe the method (series,
        /// continued fraction and Temme's uniform expansion). The value is accurate to about 1E-14 relative and
        /// underflows only below the double range.
        /// </para>
        /// <para>
        /// A NaN argument returns NaN, and X = +∞ returns 1.
        /// </para>
        /// <para>
        /// References: Shea, B. L. (1988). Algorithm AS 239: Chi-squared and incomplete gamma integral. Applied
        /// Statistics 37(3), 466-473.
        /// </para>
        /// </remarks>
        /// <param name="X">Argument of function (upper limit of integration).</param>
        /// <param name="alpha">The shape parameter.</param>
        /// <returns>
        /// The incomplete gamma integral evaluated with the given upper limit x
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when alpha &lt;= 0 or X &lt; 0.</exception>
        public static double Incomplete(double X, double alpha)
        {
            if (double.IsNaN(X) || double.IsNaN(alpha))
            {
                return double.NaN;
            }

            if (alpha <= 0.0d)
            {
                throw new ArgumentOutOfRangeException("alpha", "Shape parameter (alpha) out of range.");
            }

            if (X < 0.0d)
            {
                throw new ArgumentOutOfRangeException("X", "Argument of function (X) out of range.");
            }

            return LowerIncomplete(alpha, X);
        }

        /// <summary>
        /// Upper incomplete regularized Gamma function Q
        /// (a.k.a the incomplete complemented Gamma function)
        /// </summary>
        /// <remarks>
        /// This function is equivalent to Q(a, x) = Γ(a, x) / Γ(a).
        /// <para>
        /// Edge values: Q(a, x) = 1 for a &gt; 0 and x &lt;= 0 (the CDF convention below the support),
        /// Q(0, x) = 0 for x &gt; 0 as in scipy.special.gammaincc, and a &lt; 0, or a = 0 with
        /// x &lt;= 0, return NaN. An infinite shape gives Q = 1 for finite x, and NaN with x = +∞.
        /// </para>
        /// <para>
        /// The method is that of <see cref="LowerIncomplete(double, double)"/>, with the smaller of P and Q always
        /// evaluated directly, so a small upper tail keeps its relative precision.
        /// </para>
        /// </remarks>
        /// <param name="x"> The value to be evaluated (lower limit of integration) </param>
        /// <param name="a">The shape parameter, a &gt;= 0.</param>
        /// <returns>
        /// The upper incomplete Gamma function evaluated with the lower limit x
        /// </returns>
        public static double UpperIncomplete(double a, double x)
        {
            if (a < 0d)
            {
                return double.NaN;
            }

            if (a == 0d)
            {
                return x > 0d ? 0.0d : double.NaN;
            }

            if (x <= 0d)
            {
                return 1.0d;
            }

            if (double.IsNaN(a) || double.IsNaN(x))
            {
                return double.NaN;
            }

            if (double.IsPositiveInfinity(a))
            {
                return double.IsPositiveInfinity(x) ? double.NaN : 1d;
            }

            if (double.IsPositiveInfinity(x))
            {
                return 0d;
            }

            return IncompleteRatio(a, x, true);
        }

        /// <summary>Immutable integer zeta values used by the GEV series.</summary>
        private static readonly double[] IntegerZetaValues = CreateIntegerZetaValues();

        /// <summary>Precomputes the exact integer zeta range consumed by the bounded series.</summary>
        /// <returns>The zeta values for integer arguments 2 through 59, indexed by argument minus two.</returns>
        /// <remarks>The existing ordered power sum for arguments 17 through 59 is evaluated once during static initialization.</remarks>
        private static double[] CreateIntegerZetaValues()
        {
            double[] values = { 1.6449340668482264365, 1.2020569031595942854, 1.0823232337111381915,
                1.0369277551433699263, 1.0173430619844491397, 1.0083492773819228268,
                1.0040773561979443394, 1.0020083928260822144, 1.0009945751278180853,
                1.0004941886041194646, 1.0002460865533080483, 1.0001227133475784891,
                1.0000612481350587048, 1.0000305882363070205, 1.0000152822594086519 };
            Array.Resize(ref values, 58);
            for (int n = 17; n < 60; n++)
            {
                double sum = 1;
                for (int k = 2; k <= 32; k++) sum += Math.Pow(k, -n);
                values[n - 2] = sum;
            }
            return values;
        }

        /// <summary>
        /// Returns the gamma function ratio Γ(x) / Γ(x + 1/2) for x &gt; 0.
        /// </summary>
        /// <param name="x">The positive argument.</param>
        /// <returns>Γ(x) / Γ(x + 1/2).</returns>
        /// <remarks>
        /// For x &lt; 1 the ratio is formed directly from <see cref="Function(double)"/>. For x &gt;= 1 it is
        /// exp(ln B(x, 1/2) - ln Γ(1/2)), since Γ(x) / Γ(x + 1/2) = B(x, 1/2) / Γ(1/2), with ln B from
        /// <see cref="Beta.LogFunction(double, double)"/>, which does not subtract large ln Γ values and does not overflow.
        /// Against mpmath the relative error is below 3 ulps for x &lt; 1, below 2.5 ulps for 1 &lt;= x &lt; 1000 and below 8
        /// ulps up to 1E8, where the exponential carries the rounding of ln B; the direct quotient reaches 300 ulps near
        /// x = 150.
        /// </remarks>
        internal static double RatioHalf(double x)
        {
            if (x < 1.0) return Function(x) / Function(x + 0.5);
            return Math.Exp(Beta.LogFunction(x, 0.5) - LogGamma(0.5));
        }

        /// <summary>
        /// Returns ln Γ(x + c) - ln Γ(x) for x &gt; 0 and x + c &gt; 0 without subtracting large log-gamma values.
        /// </summary>
        /// <param name="x">The positive base argument.</param>
        /// <param name="c">The increment, of either sign, with x + c &gt; 0.</param>
        /// <returns>ln Γ(x + c) - ln Γ(x), or NaN when x or x + c is not a positive finite value.</returns>
        /// <remarks>
        /// The smaller of x and x + c is first raised to at least 8 with the recurrence ln Γ(z + 1) = ln Γ(z) + ln z,
        /// each step contributing ln(1 + c/z), which has the sign of c, so the steps do not cancel. The remaining
        /// difference is ln(Γ(b) / Γ(a + b)) with b ≥ 8 and a = |c| (routine ALGDIV of DiDonato and Morris, 1992,
        /// in <see cref="Beta"/>), which keeps the absolute accuracy of the result when ln Γ(x) and ln Γ(x + c) are
        /// large and nearly equal, for example for a small increment c and a large x.
        /// </remarks>
        internal static double LogGammaDifference(double x, double c)
        {
            double low = Math.Min(x, x + c);
            if (!(low > 0d) || !Tools.IsFinite(x) || !Tools.IsFinite(c)) return double.NaN;
            if (c == 0d) return 0d;
            double sum = 0d;
            while (low < 8d)
            {
                double ratio = c / x;
                sum += Math.Abs(ratio) < 0.5d ? Tools.Log1p(ratio) : Math.Log((x + c) / x);
                x += 1d;
                low += 1d;
            }
            return c > 0d ? -sum - Beta.LogGammaRatio(c, x) : -sum + Beta.LogGammaRatio(-c, x + c);
        }

        /// <summary>
        /// Returns ln Γ(1 + a) without rounding 1 + a, for a &gt; -1.
        /// </summary>
        /// <param name="a">The increment from one.</param>
        /// <returns>ln Γ(1 + a).</returns>
        /// <remarks>
        /// For -0.3203125 &lt;= a &lt;= 0.5 the value is a (D1 + a P1(a) / Q1(a)), the rational approximation of
        /// <see cref="LogGamma(double)"/> about 1 (Cody and Hillstrom, 1967; SPECFUN DLGAMA) evaluated at a itself, so
        /// it keeps full relative accuracy as a approaches zero. Elsewhere it is <see cref="LogGamma(double)"/> at 1 + a.
        /// </remarks>
        internal static double LogGammaOnePlus(double a)
        {
            if (a > 0.5 || a < -0.3203125) return LogGamma(a + 1.0);
            return a * (LogGammaD1 + a * (Evaluate.PolynomialRev(LogGammaP1, a) / Evaluate.PolynomialRev_1(LogGammaQ1, a)));
        }

        /// <summary>Integer zeta constants for the small-shape GEV moments.</summary>
        /// <param name="n">The integer zeta argument, which must be at least two.</param>
        /// <returns>The Riemann zeta value at <paramref name="n"/>.</returns>
        internal static double ZetaInteger(int n)
        {
            if (n < 60) return IntegerZetaValues[n - 2];
            double sum = 1;
            for (int k = 2; k <= 32; k++) sum += Math.Pow(k, -n);
            return sum;
        }

        /// <summary>
        /// Lower incomplete regularized gamma function P
        /// (a.k.a. the incomplete Gamma function).
        /// </summary>
        /// <remarks>
        /// This function is equivalent to P(a, x) = γ(a, x) / Γ(a).
        /// <para>
        /// Edge values: P(a, x) = 0 for a &gt; 0 and x &lt;= 0 (the CDF convention below the support),
        /// P(0, x) = 1 for x &gt; 0 as in scipy.special.gammainc, and a &lt; 0, or a = 0 with
        /// x &lt;= 0, return NaN. An infinite shape gives P = 0 for finite x, and NaN with x = +∞.
        /// </para>
        /// <para>
        /// <b> Method. </b> The organization of Cephes igam and igamc (Moshier, 1989, as revised in SciPy) is
        /// followed, and the smaller of P and Q is always evaluated directly. For a &gt; 20 with x close to a
        /// (|x - a| / a &lt; 0.3 for a &lt; 200, &lt; 4.5 / sqrt(a) above, and also &lt; 0.1 for a &gt;= 1E4, where the
        /// series would need of order sqrt(a) terms), Temme's uniform asymptotic
        /// expansion Q = erfc(η sqrt(a/2)) / 2 + exp(-aη²/2) / sqrt(2πa) Σ c_k(η) a^(-k) is summed (Temme, 1979;
        /// DiDonato and Morris, 1986; NIST DLMF 8.12), with 20 rows of the Maclaurin coefficients of c_k computed
        /// with mpmath by docs/mathematics/oracles/special-functions/temme_coefficients.py. Elsewhere the power
        /// series of P (x &lt; a or x &lt;= 1), the Legendre continued fraction of Q (DLMF 8.9.2) by the modified
        /// Lentz method, or, for small x, the series of Q in powers of x are used. The factor x^a e^(-x) / Γ(a)
        /// is formed directly from <see cref="Function(double)"/> for a &lt;= 20 and, above, from the Stirling
        /// form a (ln(1 + δ) - δ) + ln(a / (2π)) / 2 minus the Stirling series remainder, with δ = (x - a) / a,
        /// so it does not cancel and does not underflow before the result does.
        /// </para>
        /// <para>
        /// <b> Accuracy. </b> Measured against mpmath at 50 digits over a in [1E-3, 1E6] and x from 1E-4 a to
        /// 30 a, the relative error is about 1E-14 or below where the result exceeds about 1E-30; deeper in the
        /// tails it grows in proportion to |ln P| (or |ln Q|), the conditioning of the exponential factor
        /// (see docs/mathematics/oracles/special-functions).
        /// </para>
        /// <para>
        /// References: Temme, N. M. (1979). The asymptotic expansion of the incomplete gamma functions. SIAM
        /// Journal on Mathematical Analysis 10(4), 757-766. DiDonato, A. R. and Morris, A. H. (1986). Computation
        /// of the incomplete gamma function ratios and their inverse. ACM Transactions on Mathematical Software
        /// 12(4), 377-393. Moshier, S. L. (1989). Methods and Programs for Mathematical Functions. Ellis Horwood.
        /// </para>
        /// </remarks>
        /// <param name="a">The shape parameter, a &gt;= 0.</param>
        /// <param name="x"> The value to be evaluated (upper limit of integration) </param>
        /// <returns>
        /// The lower incomplete Gamma function evaluated with the upper limit x
        /// </returns>
        public static double LowerIncomplete(double a, double x)
        {
            if (a < 0d)
            {
                return double.NaN;
            }

            if (a == 0d)
            {
                return x > 0d ? 1.0d : double.NaN;
            }

            if (x <= 0d)
            {
                return 0.0d;
            }

            if (double.IsNaN(a) || double.IsNaN(x))
            {
                return double.NaN;
            }

            if (double.IsPositiveInfinity(a))
            {
                return double.IsPositiveInfinity(x) ? double.NaN : 0d;
            }

            if (double.IsPositiveInfinity(x))
            {
                return 1d;
            }

            return IncompleteRatio(a, x, false);
        }

        /// <summary>
        /// The largest shape for which the incomplete gamma factor x^a e^(-x) / Γ(a) is formed from Γ(a) directly.
        /// </summary>
        private const double IncompleteDirectShape = 20.0;

        /// <summary>
        /// The shape from which Temme's expansion also covers |x - a| / a &lt; 0.1, where the series and continued fraction
        /// need of order sqrt(a) terms.
        /// </summary>
        private const double TemmeWideShape = 1E4;

        /// <summary>
        /// Returns true where <see cref="IncompleteTemme(double, double, bool)"/> evaluates the ratio.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <returns>True for a &gt; 20 with |x - a| / a below 0.3 (a &lt; 200), below 4.5 / sqrt(a) (a &gt;= 200) or below 0.1 (a &gt;= 1E4).</returns>
        /// <remarks>
        /// The first two bands are those of Cephes igam. The series of P needs about 37 a / |x - a| terms, so
        /// outside the 4.5 / sqrt(a) band it would need of order 8 sqrt(a) terms and fail to converge for a above
        /// about 1E10; for a &gt;= 1E4 the expansion is therefore also used for |x - a| / a &lt; 0.1, where η is below
        /// 0.11 and its Maclaurin coefficients converge to double precision.
        /// </remarks>
        private static bool UseTemme(double a, double x)
        {
            if (!(a > IncompleteDirectShape)) return false;
            double delta = Math.Abs(x - a) / a;
            return (a < 200.0 && delta < 0.3) || (a >= 200.0 && delta < 4.5 / Math.Sqrt(a))
                || (a >= TemmeWideShape && delta < 0.1);
        }

        /// <summary>
        /// Evaluates the regularized incomplete gamma ratio P(a, x) or Q(a, x) for finite positive a and x.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <param name="upper">True for Q(a, x); false for P(a, x).</param>
        /// <returns>The requested ratio.</returns>
        /// <remarks>
        /// The region selection is that of Cephes igam and igamc as revised in SciPy: the asymptotic band
        /// for a &gt; 20 (<see cref="UseTemme(double, double)"/>, widened for a &gt;= 1E4), then the series of P where x &lt; a or x &lt;= 1, the continued fraction of Q where
        /// x &gt; 1.1 and x &gt;= a, and the small-x series of Q where it converges faster than that of P. The
        /// complement 1 - P or 1 - Q is taken only for the larger of the two ratios.
        /// </remarks>
        /// <exception cref="ArithmeticException">Thrown when the series or continued fraction does not converge.</exception>
        private static double IncompleteRatio(double a, double x, bool upper)
        {
            if (UseTemme(a, x))
                return IncompleteTemme(a, x, upper);

            if (x > 1.1)
            {
                if (x < a)
                {
                    double p = IncompleteSeries(a, x);
                    return upper ? 1.0 - p : p;
                }
                double q = IncompleteContinuedFraction(a, x);
                return upper ? q : 1.0 - q;
            }

            // Small x: the series of Q in powers of x converges faster than the series of P when a is small.
            bool qSeries = x <= 0.5 ? -0.4 / Math.Log(x) >= a : x * 1.1 >= a;
            if (qSeries)
            {
                double q = Math.Exp(LogSmallUpperIncomplete(a, x, out _));
                return upper ? q : 1.0 - q;
            }
            double ps = IncompleteSeries(a, x);
            return upper ? 1.0 - ps : ps;
        }

        /// <summary>
        /// Returns the natural logarithm of P(a, x) or Q(a, x) for finite positive a and x, including tails below the
        /// double range.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <param name="upper">True for ln Q(a, x); false for ln P(a, x).</param>
        /// <returns>The requested logarithm, at most zero.</returns>
        /// <remarks>
        /// The regions and sums are those of <see cref="IncompleteRatio(double, double, bool)"/>; the factor
        /// x^a e^(-x) / Γ(a) is carried as its logarithm (<see cref="LogIncompleteFactor(double, double)"/>), so a ratio
        /// that underflows keeps a finite logarithm, and the larger ratio is ln(1 - e^ln(smaller)) by log1p.
        /// </remarks>
        /// <exception cref="ArithmeticException">Thrown when the series or continued fraction does not converge.</exception>
        internal static double LogIncompleteRatio(double a, double x, bool upper)
        {
            if (UseTemme(a, x))
                return LogIncompleteTemme(a, x, upper);

            double logP, logQ;
            if (x > 1.1)
            {
                if (x < a)
                {
                    logP = LogIncompleteFactor(a, x) + Math.Log(IncompleteSeriesSum(a, x)) - Math.Log(a);
                    return upper ? Tools.Log1mExp(logP) : logP;
                }
                logQ = LogIncompleteFactor(a, x) + Math.Log(IncompleteContinuedFractionSum(a, x));
                return upper ? logQ : Tools.Log1mExp(logQ);
            }

            bool qSeries = x <= 0.5 ? -0.4 / Math.Log(x) >= a : x * 1.1 >= a;
            if (qSeries)
            {
                logQ = LogSmallUpperIncomplete(a, x, out _);
                return upper ? logQ : Tools.Log1mExp(logQ);
            }
            logP = LogIncompleteFactor(a, x) + Math.Log(IncompleteSeriesSum(a, x)) - Math.Log(a);
            return upper ? Tools.Log1mExp(logP) : logP;
        }

        /// <summary>
        /// Returns the natural logarithm of P(a, x) or Q(a, x) and its derivative with respect to the shape a at fixed x, for
        /// finite positive a and x.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <param name="upper">True for ln Q(a, x); false for ln P(a, x).</param>
        /// <param name="shapeDerivative">The derivative of the returned logarithm with respect to a at fixed x.</param>
        /// <returns>The requested logarithm, equal to <see cref="LogIncompleteRatio(double, double, bool)"/>.</returns>
        /// <remarks>
        /// The regions and sums are those of <see cref="LogIncompleteRatio(double, double, bool)"/>, and each sum carries its
        /// derivative with respect to a alongside: the power series of P term by term, the continued fraction of Q through
        /// the derivative of the modified Lentz recurrence, the small-x series of Q as in
        /// <see cref="LogSmallUpperIncomplete(double, double, out double)"/>, and Temme's expansion through dη/da = -δ / (η a)
        /// with δ = (x - a) / a. Each sum returns the value at the step where the value alone converges, so the logarithm is
        /// the same as without the derivative, and continues until the derivative has also converged. The derivative of the
        /// factor ln(x^a e^(-x) / Γ(a)) is ln x - ψ(a), formed for a &gt;= 16 as ln(x / a) + ln a - ψ(a) with the
        /// asymptotic series of ln a - ψ(a), so it does not cancel when x is close to a. The derivative of the larger ratio
        /// follows from that of the smaller one, d ln(1 - e^L) = -e^L / (1 - e^L) dL.
        /// </remarks>
        /// <exception cref="ArithmeticException">Thrown when the series or continued fraction does not converge.</exception>
        internal static double LogIncompleteRatio(double a, double x, bool upper, out double shapeDerivative)
        {
            double log, derivative;
            bool lowerComputed;
            if (UseTemme(a, x))
            {
                return LogIncompleteTemme(a, x, upper, out shapeDerivative);
            }

            if (x > 1.1 && x >= a)
            {
                double fraction = IncompleteContinuedFractionSum(a, x, out double fractionDerivative);
                log = LogIncompleteFactor(a, x) + Math.Log(fraction);
                derivative = LogIncompleteFactorShapeDerivative(a, x) + fractionDerivative / fraction;
                lowerComputed = false;
            }
            else if (x <= 1.1 && (x <= 0.5 ? -0.4 / Math.Log(x) >= a : x * 1.1 >= a))
            {
                log = LogSmallUpperIncomplete(a, x, out derivative);
                lowerComputed = false;
            }
            else
            {
                double sum = IncompleteSeriesSum(a, x, out double sumDerivative);
                log = LogIncompleteFactor(a, x) + Math.Log(sum) - Math.Log(a);
                derivative = LogIncompleteFactorShapeDerivative(a, x) + sumDerivative / sum - 1.0 / a;
                lowerComputed = true;
            }

            if (upper != lowerComputed)
            {
                shapeDerivative = derivative;
                return log;
            }
            double complement = Tools.Log1mExp(log);
            shapeDerivative = -derivative * Math.Exp(log - complement);
            return complement;
        }

        /// <summary>
        /// Returns the derivative of ln(x^a e^(-x) / Γ(a)) with respect to a at fixed x, ln x - ψ(a).
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <returns>The shape derivative of the logarithm of the incomplete gamma factor.</returns>
        /// <remarks>
        /// For a &lt; 16 it is ln x - <see cref="Digamma(double)"/>. Above, it is ln(x / a) + ln a - ψ(a), with
        /// ln a - ψ(a) = 1/(2a) + 1/(12a²) - 1/(120a⁴) + 1/(252a⁶) - 1/(240a⁸) (DLMF 5.11.2) and ln(x / a) = ln(1 + δ) by
        /// <see cref="Tools.Log1p(double)"/> for |δ| &lt; 0.25, δ = (x - a) / a, so the difference does not cancel.
        /// </remarks>
        private static double LogIncompleteFactorShapeDerivative(double a, double x)
        {
            if (a < 16) return Math.Log(x) - Digamma(a);
            double inverse = 1 / a, square = inverse * inverse;
            double residual = 0.5 * inverse + square * (1.0 / 12 - square * (1.0 / 120 - square * (1.0 / 252 - square / 240)));
            double delta = (x - a) / a;
            return (Math.Abs(delta) < 0.25 ? Tools.Log1p(delta) : Math.Log(x) - Math.Log(a)) + residual;
        }

        /// <summary>Returns ln Q(a, x) from the series of Q in powers of x, for small a and x, avoiding the subtraction 1 - P near one.</summary>
        /// <param name="a">The positive gamma shape in the small-shape branch.</param>
        /// <param name="x">The positive unit-scale variate in the small-argument branch.</param>
        /// <param name="derivative">The fixed-<paramref name="x"/> derivative of the returned log upper tail with respect to <paramref name="a"/>.</param>
        /// <returns>The natural logarithm of the regularized upper incomplete gamma integral.</returns>
        internal static double LogSmallUpperIncomplete(double a, double x, out double derivative)
        {
            double sum = 0, dsum = 0, power = 1;
            for (int n = 1; n <= 1000; n++)
            {
                power *= -x / n;
                double term = power / (a + n);
                sum += term;
                dsum -= term / (a + n);
                if (Math.Abs(term) < Math.Abs(sum) * 2E-16) break;
            }
            double logx = Math.Log(x), u = a * logx - LogGammaOnePlus(a);
            double leading = Math.Exp(u), factor = a * leading;
            double q = -Tools.Expm1(u) - factor * sum;
            double du = logx - Digamma(a + 1);
            derivative = (-leading * du - factor * (du * sum + dsum) - leading * sum) / q;
            return Math.Log(q);
        }

        /// <summary>Returns the Stirling series remainder ln Γ(a) - ((a - 1/2) ln a - a + ln √(2π)), through the a^(-11) term, for a of at least sixteen.</summary>
        /// <param name="a">The gamma argument, expected to be at least sixteen.</param>
        /// <returns>The retained Stirling-series correction to the leading log-gamma terms.</returns>
        internal static double StirlingRemainder(double a)
        {
            double r = 1 / a, r2 = r * r;
            return r * (1.0 / 12 + r2 * (-1.0 / 360 + r2 * (1.0 / 1260 + r2 * (-1.0 / 1680 + r2 * (1.0 / 1188 - r2 * 691.0 / 360360)))));
        }

        /// <summary>
        /// Returns ln(x^a e^(-x) / Γ(a)) for finite positive a and x, the logarithm of <see cref="IncompleteFactor(double, double)"/>.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <returns>The logarithm of the incomplete gamma factor.</returns>
        /// <remarks>
        /// For a &lt;= 20 it is the logarithm of the direct product x^a e^(-x) / Γ(a) where that product neither overflows
        /// nor underflows, which avoids the cancellation of a ln x - x - ln Γ(a), and that sum otherwise; above 20, the
        /// Stirling form of <see cref="IncompleteFactor(double, double)"/>.
        /// </remarks>
        internal static double LogIncompleteFactor(double a, double x)
        {
            if (a <= IncompleteDirectShape)
            {
                double gamma = Function(a);
                double logFactor = a * Math.Log(x) - x - Math.Log(gamma);
                if (x < 700.0 && logFactor > -700.0 && Math.Abs(a * Math.Log(x)) < 700.0)
                    return Math.Log(Math.Pow(x, a) * Math.Exp(-x) / gamma);
                return logFactor;
            }
            double delta = (x - a) / a;
            double deviance = x < 0.5 * a
                ? a * Math.Log(x / a) + (a - x)
                : a * Tools.Log1pMinusX(delta);
            return deviance + 0.5 * Math.Log(a) - Tools.LogSqrt2PI - StirlingRemainder(a);
        }

        /// <summary>
        /// Returns x^a e^(-x) / Γ(a) for finite positive a and x without cancellation or premature underflow.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument.</param>
        /// <returns>The incomplete gamma factor, zero only when it underflows the double range.</returns>
        /// <remarks>
        /// For a &lt;= 20 the factor is x^a e^(-x) / Γ(a) from <see cref="Function(double)"/>, multiplied
        /// directly where no part overflows or underflows and through its logarithm otherwise. Above 20 its
        /// logarithm is a (ln(1 + δ) - δ) + ln(a / (2π)) / 2 - r(a) with δ = (x - a) / a and r the remainder of the
        /// Stirling series of ln Γ(a) (Abramowitz and Stegun 6.1.41), where ln(1 + δ) - δ is summed as a series
        /// for |δ| &lt; 0.25 and, for x &lt; a / 2, taken as ln(x / a) + (a - x) / a so that 1 + δ is never rounded.
        /// </remarks>
        private static double IncompleteFactor(double a, double x)
        {
            if (a <= IncompleteDirectShape)
            {
                double gamma = Function(a);
                double logFactor = a * Math.Log(x) - x - Math.Log(gamma);
                if (x < 700.0 && logFactor > -700.0 && Math.Abs(a * Math.Log(x)) < 700.0)
                    return Math.Pow(x, a) * Math.Exp(-x) / gamma;
                return Math.Exp(logFactor);
            }

            double delta = (x - a) / a;
            double deviance = x < 0.5 * a
                ? a * Math.Log(x / a) + (a - x)
                : a * Tools.Log1pMinusX(delta);
            return Math.Exp(deviance + 0.5 * Math.Log(a) - Tools.LogSqrt2PI
                - StirlingRemainder(a));
        }

        /// <summary>
        /// The power series P(a, x) = x^a e^(-x) / Γ(a + 1) Σ x^n / ((a + 1) ... (a + n)) (DLMF 8.7.1).
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument, below a or at most 1.1.</param>
        /// <returns>P(a, x).</returns>
        /// <exception cref="ArithmeticException">Thrown when the series does not converge.</exception>
        private static double IncompleteSeries(double a, double x)
        {
            double factor = IncompleteFactor(a, x);
            if (factor == 0.0) return 0.0;
            return IncompleteSeriesSum(a, x) * factor / a;
        }

        /// <summary>
        /// The sum Σ x^n / ((a + 1) ... (a + n)), n &gt;= 0, of the power series of P(a, x) (DLMF 8.7.1).
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument, below a or at most 1.1.</param>
        /// <returns>The series sum, at least one.</returns>
        /// <exception cref="ArithmeticException">Thrown when the series does not converge.</exception>
        private static double IncompleteSeriesSum(double a, double x)
        {
            double r = a;
            double c = 1.0;
            double sum = 1.0;
            for (int n = 0; n < 1_000_000; n++)
            {
                r += 1.0;
                c *= x / r;
                sum += c;
                if (c <= sum * Tools.DoubleMachineEpsilon)
                    return sum;
            }
            throw new ArithmeticException("The incomplete gamma series did not converge.");
        }

        /// <summary>
        /// The sum of <see cref="IncompleteSeriesSum(double, double)"/> with its derivative with respect to a.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite positive argument, below a or at most 1.1.</param>
        /// <param name="derivative">The derivative of the sum with respect to a.</param>
        /// <returns>The series sum, the same value as <see cref="IncompleteSeriesSum(double, double)"/>.</returns>
        /// <remarks>
        /// Each term c_n = c_(n-1) x / (a + n) has the derivative c_n' = c_(n-1)' x / (a + n) - c_n / (a + n). The sum is
        /// taken at the step where the value converges, and the derivative terms are added until they are also below the
        /// double precision of their running total.
        /// </remarks>
        /// <exception cref="ArithmeticException">Thrown when the series does not converge.</exception>
        private static double IncompleteSeriesSum(double a, double x, out double derivative)
        {
            double r = a;
            double c = 1.0, dc = 0.0;
            double sum = 1.0, dsum = 0.0;
            double value = double.NaN;
            for (int n = 0; n < 1_000_000; n++)
            {
                r += 1.0;
                double ratio = x / r;
                dc = dc * ratio - c * ratio / r;
                c *= ratio;
                dsum += dc;
                if (double.IsNaN(value))
                {
                    sum += c;
                    if (c <= sum * Tools.DoubleMachineEpsilon) value = sum;
                }
                if (!double.IsNaN(value) && Math.Abs(dc) <= Math.Abs(dsum) * Tools.DoubleMachineEpsilon)
                {
                    derivative = dsum;
                    return value;
                }
            }
            throw new ArithmeticException("The incomplete gamma series did not converge.");
        }

        /// <summary>
        /// The Legendre continued fraction for Q(a, x) (DLMF 8.9.2), summed by the modified Lentz method.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite argument, greater than 1.1 and at least a.</param>
        /// <returns>Q(a, x).</returns>
        /// <exception cref="ArithmeticException">Thrown when the continued fraction does not converge.</exception>
        private static double IncompleteContinuedFraction(double a, double x)
        {
            double factor = IncompleteFactor(a, x);
            if (factor == 0.0) return 0.0;
            return factor * IncompleteContinuedFractionSum(a, x);
        }

        /// <summary>
        /// The Legendre continued fraction for Q(a, x) divided by the factor x^a e^(-x) / Γ(a) (DLMF 8.9.2), summed by the
        /// modified Lentz method.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite argument, greater than 1.1 and at least a.</param>
        /// <returns>The continued fraction.</returns>
        /// <exception cref="ArithmeticException">Thrown when the continued fraction does not converge.</exception>
        private static double IncompleteContinuedFractionSum(double a, double x)
        {
            const double Tiny = 1E-300;
            double b = x + 1.0 - a;
            double c = 1.0 / Tiny;
            double d = 1.0 / b;
            double h = d;
            for (int i = 1; i < 1_000_000; i++)
            {
                double an = -i * (i - a);
                b += 2.0;
                d = an * d + b;
                if (Math.Abs(d) < Tiny) d = Tiny;
                c = b + an / c;
                if (Math.Abs(c) < Tiny) c = Tiny;
                d = 1.0 / d;
                double delta = d * c;
                h *= delta;
                if (Math.Abs(delta - 1.0) <= Tools.DoubleMachineEpsilon)
                    return h;
            }
            throw new ArithmeticException("The incomplete gamma continued fraction did not converge.");
        }

        /// <summary>
        /// The continued fraction of <see cref="IncompleteContinuedFractionSum(double, double)"/> with its derivative with
        /// respect to a.
        /// </summary>
        /// <param name="a">The finite positive shape.</param>
        /// <param name="x">The finite argument, greater than 1.1 and at least a.</param>
        /// <param name="derivative">The derivative of the continued fraction with respect to a.</param>
        /// <returns>The continued fraction, the same value as <see cref="IncompleteContinuedFractionSum(double, double)"/>.</returns>
        /// <remarks>
        /// The Lentz quantities d, c and h are differentiated step by step, with da_i/da = i and db_i/da = -1 for the partial
        /// numerators a_i = -i (i - a) and denominators b_i = x + 1 - a + 2i. The fraction is taken at the step where the value
        /// converges, and the recurrence continues until the derivative factor is also below double precision.
        /// </remarks>
        /// <exception cref="ArithmeticException">Thrown when the continued fraction does not converge.</exception>
        private static double IncompleteContinuedFractionSum(double a, double x, out double derivative)
        {
            const double Tiny = 1E-300;
            double b = x + 1.0 - a;
            double c = 1.0 / Tiny, dc = 0.0;
            double d = 1.0 / b, dd = d * d;
            double h = d, dh = dd;
            double value = double.NaN;
            for (int i = 1; i < 1_000_000; i++)
            {
                double an = -i * (i - a);
                b += 2.0;
                double ddenominator = i * d + an * dd - 1.0;
                d = an * d + b;
                if (Math.Abs(d) < Tiny) d = Tiny;
                dc = -1.0 + i / c - an * dc / (c * c);
                c = b + an / c;
                if (Math.Abs(c) < Tiny) c = Tiny;
                d = 1.0 / d;
                dd = -ddenominator * d * d;
                double delta = d * c;
                double ddelta = dd * c + d * dc;
                dh = dh * delta + h * ddelta;
                h *= delta;
                if (double.IsNaN(value) && Math.Abs(delta - 1.0) <= Tools.DoubleMachineEpsilon) value = h;
                if (!double.IsNaN(value) && Math.Abs(ddelta) <= Tools.DoubleMachineEpsilon)
                {
                    derivative = dh;
                    return value;
                }
            }
            throw new ArithmeticException("The incomplete gamma continued fraction did not converge.");
        }

        /// <summary>
        /// Temme's uniform asymptotic expansion of P(a, x) or Q(a, x) for large a with x near a.
        /// </summary>
        /// <param name="a">The shape, greater than 20.</param>
        /// <param name="x">The argument, close to a.</param>
        /// <param name="upper">True for Q(a, x); false for P(a, x).</param>
        /// <returns>The requested ratio.</returns>
        /// <remarks>
        /// With λ = x/a and η = sign(λ - 1) sqrt(-2 (ln λ - (λ - 1))), Q = erfc(η sqrt(a/2)) / 2 + R and
        /// P = erfc(-η sqrt(a/2)) / 2 - R, where R = exp(-aη²/2) / sqrt(2πa) Σ c_k(η) a^(-k) (NIST DLMF 8.12.3,
        /// 8.12.4 and 8.12.8). Each c_k(η) is summed from its Maclaurin coefficients, and the sums stop once a
        /// term is below the double precision of the running total, as in Cephes igam_asymptotic_series.
        /// </remarks>
        private static double IncompleteTemme(double a, double x, bool upper)
        {
            double eta = TemmeEta(a, x);
            double sum = TemmeSum(a, eta);
            double sign = upper ? 1.0 : -1.0;
            double r = Math.Exp(-0.5 * a * eta * eta) * sum / Math.Sqrt(2.0 * Math.PI * a);
            return 0.5 * Erf.Erfc(sign * eta * Math.Sqrt(0.5 * a)) + sign * r;
        }

        /// <summary>
        /// Temme's expansion of <see cref="IncompleteTemme(double, double, bool)"/> in logarithmic form.
        /// </summary>
        /// <param name="a">The shape, greater than 20.</param>
        /// <param name="x">The argument, close to a.</param>
        /// <param name="upper">True for ln Q(a, x); false for ln P(a, x).</param>
        /// <returns>The logarithm of the requested ratio, finite where the ratio underflows.</returns>
        /// <remarks>
        /// With z = η sqrt(a), ln Q = ln Φ(-z) + ln(1 + φ(z) / (sqrt(a) Φ(-z)) Σ c_k(η) a^(-k)) and
        /// ln P = ln Φ(z) + ln(1 - φ(z) / (sqrt(a) Φ(z)) Σ c_k(η) a^(-k)), where φ is the standard normal density.
        /// </remarks>
        private static double LogIncompleteTemme(double a, double x, bool upper)
        {
            double eta = TemmeEta(a, x);
            double sum = TemmeSum(a, eta);
            double z = eta * Math.Sqrt(a);
            double normal = DistributionNumerics.NormalLogCDF(upper ? -z : z);
            double logDensity = -(0.5 * z) * z - Tools.LogSqrt2PI;
            double correction = (upper ? 1.0 : -1.0) * sum * Math.Exp(logDensity - 0.5 * Math.Log(a) - normal);
            return normal + Tools.Log1p(correction);
        }

        /// <summary>
        /// Temme's expansion of <see cref="LogIncompleteTemme(double, double, bool)"/> with the derivative of the logarithm
        /// with respect to a at fixed x.
        /// </summary>
        /// <param name="a">The shape, greater than 20.</param>
        /// <param name="x">The argument, close to a.</param>
        /// <param name="upper">True for ln Q(a, x); false for ln P(a, x).</param>
        /// <param name="shapeDerivative">The derivative of the returned logarithm with respect to a.</param>
        /// <returns>The same logarithm as <see cref="LogIncompleteTemme(double, double, bool)"/>.</returns>
        /// <remarks>
        /// With δ = (x - a) / a, dη/da = -δ / (η a) (-1/a at δ = 0) and dz/da = η / (2 sqrt(a)) + sqrt(a) dη/da. The sum
        /// S = Σ c_k(η) a^(-k) is differentiated in η through the derivatives of the Maclaurin polynomials and in a through
        /// -k c_k(η) a^(-k-1). Its value is that of <see cref="TemmeSum(double, double)"/>, and the derivative terms are added
        /// until they are below the double precision of their running total.
        /// </remarks>
        private static double LogIncompleteTemme(double a, double x, bool upper, out double shapeDerivative)
        {
            double eta = TemmeEta(a, x);
            double sum = 0.0, sumEta = 0.0, sumA = 0.0;
            double inversePower = 1.0;
            bool converged = false;
            for (int k = 0; k < TemmeCoefficients.Length; k++)
            {
                double[] row = TemmeCoefficients[k];
                double ck = 0.0, dck = 0.0;
                for (int n = row.Length - 1; n >= 0; n--)
                {
                    dck = dck * eta + ck;
                    ck = ck * eta + row[n];
                }
                double term = ck * inversePower;
                double termEta = dck * inversePower;
                double termA = -k * term / a;
                if (!converged)
                {
                    sum += term;
                    if (Math.Abs(term) <= Tools.DoubleMachineEpsilon * Math.Abs(sum)) converged = true;
                }
                sumEta += termEta;
                sumA += termA;
                if (converged && Math.Abs(termEta) <= Tools.DoubleMachineEpsilon * Math.Abs(sumEta)
                    && Math.Abs(termA) <= Tools.DoubleMachineEpsilon * Math.Abs(sumA)) break;
                inversePower /= a;
            }

            double root = Math.Sqrt(a);
            double z = eta * root;
            double normal = DistributionNumerics.NormalLogCDF(upper ? -z : z);
            double logDensity = -(0.5 * z) * z - Tools.LogSqrt2PI;
            double correction = (upper ? 1.0 : -1.0) * sum * Math.Exp(logDensity - 0.5 * Math.Log(a) - normal);
            double log = normal + Tools.Log1p(correction);

            double delta = (x - a) / a;
            double etaDerivative = delta == 0.0 ? -1.0 / a : -delta / eta / a;
            double zDerivative = eta / (2.0 * root) + root * etaDerivative;
            double correctionDerivative = (sumEta * etaDerivative + sumA - (z * zDerivative + 0.5 / a) * sum) / root;
            shapeDerivative = Math.Exp(logDensity - log) * (upper ? -zDerivative + correctionDerivative : zDerivative - correctionDerivative);
            return log;
        }

        /// <summary>
        /// Returns η = sign(λ - 1) sqrt(-2 (ln λ - (λ - 1))) with λ = x / a, the variable of Temme's expansion.
        /// </summary>
        /// <param name="a">The positive shape.</param>
        /// <param name="x">The positive argument.</param>
        /// <returns>η, with the sign of x - a.</returns>
        private static double TemmeEta(double a, double x)
        {
            double sigma = (x - a) / a;
            return sigma == 0.0 ? 0.0
                : Math.Sign(sigma) * Math.Sqrt(-2.0 * Tools.Log1pMinusX(sigma));
        }

        /// <summary>
        /// Returns Σ c_k(η) a^(-k), each c_k summed from its Maclaurin coefficients, stopping once a term is below the
        /// double precision of the running total, as in Cephes igam_asymptotic_series.
        /// </summary>
        /// <param name="a">The shape, greater than 20.</param>
        /// <param name="eta">Temme's variable η.</param>
        /// <returns>The sum of the expansion.</returns>
        private static double TemmeSum(double a, double eta)
        {
            double sum = 0.0;
            double inversePower = 1.0;
            for (int k = 0; k < TemmeCoefficients.Length; k++)
            {
                double[] row = TemmeCoefficients[k];
                double ck = 0.0;
                for (int n = row.Length - 1; n >= 0; n--)
                    ck = ck * eta + row[n];
                double term = ck * inversePower;
                sum += term;
                if (Math.Abs(term) <= Tools.DoubleMachineEpsilon * Math.Abs(sum)) break;
                inversePower /= a;
            }
            return sum;
        }

        /// <summary>
        /// Maclaurin coefficients d[k][n] of Temme's c_k(η) = Σ d[k][n] η^n for k = 0 to 19 and n = 0 to 24.
        /// </summary>
        /// <remarks>
        /// c_0(η) = 1/(λ - 1) - 1/η and c_k(η) = (1/η) c_(k-1)'(η) + (-1)^k g_k / (λ - 1), with g_k the Stirling
        /// coefficients of Γ*(a) (NIST DLMF 8.12.9 to 8.12.11 and 5.11.3). The table was computed with mpmath power
        /// series at 80 digits by docs/mathematics/oracles/special-functions/temme_coefficients.py, which also checks
        /// the expansion against mpmath to 2E-28 relative for a from 25 to 2000; the leading rows agree with the
        /// table of Cephes igam.
        /// </remarks>
        internal static readonly double[][] TemmeCoefficients =
        {
            new double[]
            {
                -0.3333333333333333, 0.08333333333333333, -0.014814814814814815, 0.0011574074074074073,
                0.0003527336860670194, -0.0001787551440329218, 3.919263178522438E-05, -2.185448510679992E-06,
                -1.85406221071516E-06, 8.296711340953087E-07, -1.7665952736826078E-07, 6.707853543401498E-09,
                1.0261809784240309E-08, -4.382036018453353E-09, 9.14769958223679E-10, -2.5514193994946248E-11,
                -5.830772132550426E-11, 2.4361948020667415E-11, -5.0276692801141755E-12, 1.1004392031956135E-13,
                3.371763262400985E-13, -1.392388722418162E-13, 2.8534893807047445E-14, -5.139111834242572E-16,
                -1.9752288294349442E-15
            },
            new double[]
            {
                -0.001851851851851852, -0.003472222222222222, 0.0026455026455026454, -0.0009902263374485596,
                0.00020576131687242798, -4.018775720164609E-07, -1.8098550334489977E-05, 7.64916091608111E-06,
                -1.6120900894563446E-06, 4.647127802807434E-09, 1.378633446915721E-07, -5.752545603517705E-08,
                1.1951628599778148E-08, -1.7543241719747647E-11, -1.0091543710600413E-09, 4.162792991842583E-10,
                -8.56390702649298E-11, 6.067215101604758E-14, 7.1624989648114856E-12, -2.933186643771437E-12,
                5.996696365683689E-13, -2.1671786527323313E-16, -4.978339972369262E-14, 2.0291628823713425E-14,
                -4.13125571381061E-15
            },
            new double[]
            {
                0.004133597883597883, -0.0026813271604938273, 0.0007716049382716049, 2.0093878600823047E-06,
                -0.0001073665322636516, 5.2923448829120125E-05, -1.2760635188618728E-05, 3.423578734096138E-08,
                1.3721957309062934E-06, -6.298992138380055E-07, 1.4280614206064242E-07, -2.0477098421990866E-10,
                -1.409252991086752E-08, 6.228974084922022E-09, -1.3670488396617114E-09, 9.428356159014678E-13,
                1.2872252400089318E-10, -5.5645956134363323E-11, 1.197593554636698E-11, -4.1689782251838634E-15,
                -1.0940640427884595E-12, 4.662239946390136E-13, -9.905105763906907E-14, 1.8931876768373515E-17,
                8.859221872591127E-15
            },
            new double[]
            {
                0.0006494341563786008, 0.00022947209362139917, -0.0004691894943952557, 0.00026772063206283885,
                -7.561801671883977E-05, -2.396505113867297E-07, 1.1082654115347302E-05, -5.6749528269915965E-06,
                1.4230900732435883E-06, -2.7861080291528143E-11, -1.6958404091930278E-07, 8.099464905388083E-08,
                -1.9111168485973655E-08, 2.3928620439808118E-12, 2.0620131815488797E-09, -9.460496661855133E-10,
                2.1541049775774907E-10, -1.388823336813903E-14, -2.1894761681963938E-11, 9.790998951171684E-12,
                -2.178219188018096E-12, 6.208819573407901E-17, 2.126978363279737E-13, -9.344688791517433E-14,
                2.045367122678285E-14
            },
            new double[]
            {
                -0.0008618882909167117, 0.0007840392217200666, -0.0002990724803031902, -1.4638452578843418E-06,
                6.641498215465122E-05, -3.968365047179435E-05, 1.1375726970678419E-05, 2.507497226237533E-10,
                -1.6954149536558305E-06, 8.907507532205309E-07, -2.292934834000805E-07, 2.956794137544049E-11,
                2.8865829742708783E-08, -1.4189739437803219E-08, 3.4463580499464896E-09, -2.3024517174528067E-13,
                -3.9409233028046403E-10, 1.86023389685045E-10, -4.356323005056618E-11, 1.278600101629623E-15,
                4.67927502665792E-12, -2.149246470613483E-12, 4.908815614809652E-13, -6.33859148489156E-18,
                -5.045332069080094E-14
            },
            new double[]
            {
                -0.00033679855336635813, -6.972813758365857E-05, 0.0002772753244959392, -0.00019932570516188847,
                6.797780477937208E-05, 1.419062920643967E-07, -1.3594048189768693E-05, 8.018470256334202E-06,
                -2.291481176508095E-06, -3.252473551298454E-10, 3.4652846491085265E-07, -1.8447187191171344E-07,
                4.8240967037894184E-08, -1.7989466721743514E-14, -6.306194500013523E-09, 3.162417628774568E-09,
                -7.840924253697429E-10, 5.192679165254041E-15, 9.358944242306784E-11, -4.513426216163278E-11,
                1.0799129993116828E-11, -3.661886712685252E-17, -1.210902069055155E-12, 5.680743584990564E-13,
                -1.3249659916340829E-13
            },
            new double[]
            {
                0.0005313079364639922, -0.0005921664373536939, 0.0002708782096718045, 7.902353232660328E-07,
                -8.153969367561969E-05, 5.61168275310625E-05, -1.8329116582843375E-05, -3.0796134506033047E-09,
                3.465155368803609E-06, -2.0291327396058603E-06, 5.788792863149004E-07, 2.338630673826657E-13,
                -8.828600746330484E-08, 4.7435958880408125E-08, -1.2545415020710383E-08, 8.649648858010293E-14,
                1.6846058979264062E-09, -8.575492823577594E-10, 2.1598224929232125E-10, -7.613230520476153E-16,
                -2.6639822008536144E-11, 1.3065700536611057E-11, -3.1799163902367977E-12, 4.710976121367431E-18,
                3.6902800842763465E-13
            },
            new double[]
            {
                0.00034436760689237765, 5.171790908260592E-05, -0.00033493161081142234, 0.0002812695154763237,
                -0.00010976582244684731, -1.2741009095484485E-07, 2.7744451511563645E-05, -1.8263488805711332E-05,
                5.7876949497350525E-06, 4.93875893393627E-10, -1.0595367014026043E-06, 6.166714376110408E-07,
                -1.7562973359060463E-07, -1.297447328701544E-12, 2.695423606288966E-08, -1.4578352908731272E-08,
                3.887645959386175E-09, -3.881002251019412E-17, -5.327994173877286E-10, 2.7437977643314844E-10,
                -6.995796092070568E-11, 2.589986387486848E-17, 8.856689099669639E-12, -4.403168815871311E-12,
                1.0865561947091654E-12
            },
            new double[]
            {
                -0.0006526239185953094, 0.0008394987206720873, -0.000438297098541721, -6.969091458420552E-07,
                0.00016644846642067547, -0.00012783517679769218, 4.629953263691304E-05, 4.557909867922708E-09,
                -1.0595271125805195E-05, 6.783342904865167E-06, -2.1075476666258803E-06, -1.7213731432817144E-11,
                3.773587741611098E-07, -2.1867506700122867E-07, 6.220228804018927E-08, 6.597703826733E-16,
                -9.590386497425686E-09, 5.213214492280807E-09, -1.3991589583935709E-09, 5.382058999060575E-16,
                1.9484714275467745E-10, -1.0127287556389682E-10, 2.6077347197254926E-11, -5.090418699993299E-18,
                -3.3721464474854593E-12
            },
            new double[]
            {
                -0.0005967612901927463, -7.204895416020011E-05, 0.0006782308837667328, -0.0006401475260262758,
                0.00027750107634328704, 1.819700838046515E-07, -8.479507117068503E-05, 6.105192082501531E-05,
                -2.1073920183404862E-05, -8.858589014125599E-10, 4.5284535953805374E-06, -2.8427815022504407E-06,
                8.708234177864641E-07, 3.6886101871706966E-12, -1.534469519070206E-07, 8.862466778790695E-08,
                -2.5184812301826817E-08, -1.0225912098215092E-14, 3.896947075815478E-09, -2.1267304792235634E-09,
                5.737013552805138E-10, -1.8877498501697116E-19, -8.093153869465787E-11, 4.23827232834492E-11,
                -1.1002224534207725E-11
            },
            new double[]
            {
                0.0013324454494800656, -0.0019144384985654776, 0.0011089369134596636, 9.9324041226423E-07,
                -0.0005087450129309319, 0.00042735056665392886, -0.00016858853767910798, -8.1301893922785E-09,
                4.5284402370562144E-05, -3.127053674781734E-05, 1.044986828530338E-05, 4.8435226265680926E-11,
                -2.148256587345626E-06, 1.329369701097492E-06, -4.029569309210103E-07, -1.756787766632329E-13,
                7.014504316366825E-08, -4.040787734999483E-08, 1.1474026743371964E-08, 3.964274685356394E-18,
                -1.7804938269892715E-09, 9.748026254873165E-10, -2.6405338676507616E-10, 5.79487516340376E-18,
                3.764774955354384E-11
            },
            new double[]
            {
                0.001579727660730835, 0.00016251626278391583, -0.0020633421035543276, 0.00213896861856891,
                -0.0010108559391263003, -3.99127055299192E-07, 0.0003623502508476469, -0.00028143901463712157,
                0.00010449513336495887, 2.12114184918303E-09, -2.5779417251947842E-05, 1.7281818956040464E-05,
                -5.641377387290428E-06, -1.1024320105776174E-11, 1.1223224418895174E-06, -6.869339637952674E-07,
                2.0653236975414888E-07, 4.6714772409838506E-14, -3.5609886164949055E-08, 2.0470855345905963E-08,
                -5.809173863328336E-09, -1.3328212875828647E-16, 9.035460439133513E-10, -4.959878251733084E-10,
                1.3481607129399748E-10
            },
            new double[]
            {
                -0.004072512119514016, 0.00640336283380807, -0.004041016108167662, -2.1837328028662328E-06,
                0.002174044180125464, -0.001970044051841889, 0.0008359546974796246, 1.9445447567109655E-08,
                -0.000257793871204217, 0.00019009987368139304, -6.769649993743896E-05, -1.4440629666426571E-10,
                1.5712512518742267E-05, -1.0304008744776894E-05, 3.304517767401387E-06, 7.982976024232571E-13,
                -6.4097794149313E-07, 3.8894624761300054E-07, -1.161834764494887E-07, -2.8168086305964423E-15,
                1.9878012911297094E-08, -1.1407719956357511E-08, 3.2355857064185554E-09, 4.1759462466484876E-20,
                -5.042311271810582E-10
            },
            new double[]
            {
                -0.0059475779383993, -0.0005401647678926045, 0.00879104135507679, -0.009857631558785612,
                0.005013469503102154, 1.2807521786221875E-06, -0.0020626019342754685, 0.0017109128573523059,
                -0.000676953127141338, -6.901154567656214E-09, 0.00018855128143995903, -0.0001339521566349197,
                4.626318303352804E-05, 4.003423061332135E-11, -1.0255652921494033E-05, 6.612086372797651E-06,
                -2.0913022027253007E-06, -2.095177564960382E-13, 3.975602904199325E-07, -2.395621197881589E-07,
                7.118288338214586E-08, 8.925574871713252E-16, -1.2101547235064677E-08, 6.935061824833439E-09,
                -1.966146445385609E-09
            },
            new double[]
            {
                0.01740202778752271, -0.02952788094569912, 0.020045875571402798, 7.0289515966903405E-06,
                -0.012375421071343148, 0.011976293444235255, -0.0054156038466518525, -6.329089339641862E-08,
                0.0018855118129005065, -0.001473473274825001, 0.0005551581009770838, 5.240683441255066E-10,
                -0.00014357913535784835, 9.91812932249433E-05, -3.346083474947831E-05, -3.5755837291098967E-12,
                7.1560851960630075E-06, -4.551680262815553E-06, 1.4236576649271474E-06, 1.8803149079275236E-14,
                -2.662340389892921E-07, 1.5950642189595716E-07, -4.718751467384107E-08, -6.510781264821694E-17,
                7.979509102674677E-09
            },
            new double[]
            {
                0.03024912416090589, 0.0024817436002649977, -0.049939134373457025, 0.05991564300930787,
                -0.03248320760162339, -5.721296865210344E-06, 0.015085251778569354, -0.013261324005088445,
                0.0055515262632426145, 3.026318225703001E-08, -0.0017229548406756724, 0.0012893570099929638,
                -0.00046845138348319875, -1.8302599378930445E-10, 0.00011449739014822654, -7.737856522124447E-05,
                2.5625836246985202E-05, 1.0766165332658074E-12, -5.324680928242262E-06, 3.3496348630644643E-06,
                -1.038125312868401E-06, -5.608908533478749E-15, 1.9150821930676722E-07, -1.1418365800203775E-07,
                3.365442520915233E-08
            },
            new double[]
            {
                -0.09905102088015905, 0.17954011706123485, -0.12989606383463778, -3.1478872752284355E-05,
                0.09051063527684813, -0.0928288244111844, 0.04441211283987781, 2.7779236316835886E-07,
                -0.017229543805449696, 0.014182925050891573, -0.005621416163374734, -2.3959850918638095E-09,
                0.0016029634366079909, -0.0011606784674435774, 0.00041001337768153875, 1.8365800753181603E-11,
                -9.58442565636559E-05, 6.364306233776471E-05, -2.0762506244890635E-05, -1.1806017999805486E-13,
                4.213180823912094E-06, -2.6262241337013133E-06, 8.077062049488396E-07, 5.996409690563338E-16,
                -1.47297373744462E-07
            },
            new double[]
            {
                -0.19994542198219728, -0.015056113040026424, 0.3647023946934849, -0.4643519231173355,
                0.26640934719197895, 3.403826602714719E-05, -0.13784338709329624, 0.1276467178337056,
                -0.056213828755200985, -1.7531508854830108E-07, 0.019235592956768112, -0.015088821281095316,
                0.005740185445135012, 1.0622382710173866E-09, -0.0015335082692563998, 0.0010819320643228215,
                -0.0003737251019394563, -6.617090419433389E-12, 8.42636173809102E-05, -5.5150706827484874E-05,
                1.7769536448337793E-05, 3.879070571006568E-14, -3.5351369749902462E-06, 2.1865832127706725E-06,
                -6.681284949240542E-07
            },
            new double[]
            {
                0.7243860850402943, -1.3918010932653375, 1.0654143352413967, 0.0001876173868950258,
                -0.827055011761527, 0.8935243334782841, -0.44971003995291337, -1.6107401567546651E-06,
                0.1923559016527109, -0.1659770216004261, 0.06888222268181433, 1.391009172443142E-08,
                -0.021469115615086628, 0.016228980898865892, -0.005979601617258422, -1.1287468171928069E-10,
                0.001516745111978496, -0.0010478634293554165, 0.00035539072889105875, 8.162616536150783E-13,
                -7.777301344470886E-05, 5.029141389162907E-05, -1.6035083877747675E-05, -5.0260094840152135E-15,
                3.1369106037108428E-06
            },
            new double[]
            {
                1.6668949727276812, 0.1165462765994632, -3.3288393225018904, 4.469232548286404,
                -2.6977693045875806, -0.0002600667859891061, 1.5389017615694538, -1.4937962361134611,
                0.6888196463323315, 1.3077482004532885E-06, -0.2576296332559629, 0.2109767610212545,
                -0.08371440835921982, -7.792042747000674E-09, 0.024267923064833764, -0.017813678334552763,
                0.006397033038886307, 4.9415957307924725E-11, -0.0015554602758916928, 0.0010561196918773697,
                -0.00035277184484116647, -3.0939103274354235E-13, 7.528585452877513E-05, -4.8186515801973E-05,
                1.5227272261135095E-05
            }
        };

        /// <summary>
        /// Inverse of the <see cref="LowerIncomplete(double, double)">
        /// incomplete Gamma integral (LowerIncomplete, P)</see>.
        /// </summary>
        /// <remarks>
        /// Returns x with P(a, x) = y from <see cref="InverseIncompleteRatio(double, double, bool, double)"/>: for y &gt; 1/2
        /// the root of Q(a, x) = 1 - y is found, where 1 - y is exact, so both tails keep their relative precision.
        /// y = 0 returns 0 and y = 1 returns +∞; a quantile below the smallest positive double returns 0. NaN, a &lt;= 0 or
        /// y outside [0, 1] return NaN.
        /// </remarks>
        /// <param name="a">The shape, a &gt; 0.</param>
        /// <param name="y">The lower-tail probability.</param>
        /// <returns>
        /// The inverse of the lower incomplete Gamma function evaluated at y
        /// </returns>
        public static double InverseLowerIncomplete(double a, double y)
        {
            return InverseIncompleteRatio(a, y, false);
        }

        /// <summary>
        /// Inverse of the <see cref="UpperIncomplete(double, double)">complemented
        /// incomplete Gamma integral (UpperIncomplete, Q)</see>.
        /// </summary>
        /// <remarks>
        /// Returns x with Q(a, x) = y from <see cref="InverseIncompleteRatio(double, double, bool, double)"/>: for y &gt; 1/2
        /// the root of P(a, x) = 1 - y is found, where 1 - y is exact, so an upper probability near one (a small lower
        /// tail) keeps its relative precision. y = 0 returns +∞ and y = 1 returns 0; NaN, a &lt;= 0 or y outside [0, 1]
        /// return NaN.
        /// </remarks>
        /// <param name="a">The shape, a &gt; 0.</param>
        /// <param name="y">The upper-tail probability.</param>
        /// <returns>
        /// The inverse of the upper incomplete Gamma function evaluated at y
        /// </returns>
        public static double InverseUpperIncomplete(double a, double y)
        {
            return InverseIncompleteRatio(a, y, true);
        }

        /// <summary>
        /// Returns x with P(a, x) = p, or Q(a, x) = p when <paramref name="upper"/> is true, for a &gt; 0, computing ln Γ(a).
        /// </summary>
        /// <param name="a">The shape, a &gt; 0.</param>
        /// <param name="probability">The tail probability in [0, 1].</param>
        /// <param name="upper">True to invert Q(a, x); false to invert P(a, x).</param>
        /// <returns>The quantile of <see cref="InverseIncompleteRatio(double, double, bool, double)"/>.</returns>
        internal static double InverseIncompleteRatio(double a, double probability, bool upper)
        {
            return InverseIncompleteRatio(a, probability, upper, LogGamma(a));
        }

        /// <summary>
        /// Returns x with P(a, x) = p, or Q(a, x) = p when <paramref name="upper"/> is true, for a &gt; 0.
        /// </summary>
        /// <param name="a">The shape, a &gt; 0.</param>
        /// <param name="probability">The tail probability in [0, 1].</param>
        /// <param name="upper">True to invert Q(a, x); false to invert P(a, x).</param>
        /// <param name="logGammaA">ln Γ(a) from <see cref="LogGamma(double)"/>, so a caller that inverts many
        /// probabilities at one shape computes it once.</param>
        /// <returns>
        /// The quantile x: 0 or +∞ at the ends of [0, 1], 0 where the quantile is below the smallest positive double, and
        /// NaN for NaN, a &lt;= 0 or a probability outside [0, 1].
        /// </returns>
        /// <remarks>
        /// <para>
        /// The smaller tail is inverted (p &gt; 1/2 is replaced by the exact complement 1 - p and the other tail), so both
        /// tails keep their relative precision. a = 1 is the closed form x = -ln q or -ln(1 - p).
        /// </para>
        /// <para>
        /// The root of g(x) = ln T(a, x) - ln p, with T the tail from <see cref="LogIncompleteRatio(double, double, bool)"/>,
        /// is found by Halley's method as in Press et al. (2007, section 6.2.1): with the hazard h = f(x) / T, g' = ±h and
        /// g''/g' = (a - 1)/x - 1 ∓ h, and each step is x ← x - t / (1 - r / 2) with t = g / g' and r = t g''/g' bounded
        /// to [-1, 1]. ln p is carried as k ln 2 + ln m with p = m 2^k, the leading part exact, so g keeps its absolute
        /// accuracy for probabilities as small as 1E-300. The iterates stay inside a bracket of the root, and a step that
        /// leaves it is replaced by bisection (geometric once both ends are positive, and by squaring an upper bound below
        /// 1/2 while the lower end is 0), so the method always returns a value.
        /// It stops once a step is below 1E-8 x, which leaves an error of order 1E-24 x after the cubically convergent last
        /// step, or after 100 steps.
        /// </para>
        /// <para>
        /// Starting values: for the lower tail, the leading series term x = (p Γ(a + 1))^(1/a), formed in logarithms, where
        /// it is below a / 10; when its logarithm is below -708 that term is the result, since the next term changes ln x by
        /// about x / (a + 1). Otherwise, for a &gt; 1 the Wilson and Hilferty (1931) approximation
        /// a (1 - 1/(9a) + z/(3 sqrt(a)))³ with z the standard normal quantile of the lower tail, and for a &lt;= 1 the start
        /// of Press et al. (2007), x = (u / s)^(1/a) for u &lt; s = 1 - a (0.253 + 0.12 a) and x = 1 - ln((1 - u) / (1 - s))
        /// above, with u the lower-tail probability. For an upper tail far beyond the shape (ln q &lt; -10 and
        /// a &lt; -ln q / 10), the asymptotic relation x = -ln q - ln Γ(a) + (a - 1) ln x (DiDonato and Morris, 1986) is
        /// iterated twice instead.
        /// </para>
        /// <para>
        /// References: DiDonato, A. R. and Morris, A. H. (1986). Computation of the incomplete gamma function ratios and their
        /// inverse. ACM Transactions on Mathematical Software 12(4), 377-393. Press, W. H., Teukolsky, S. A., Vetterling,
        /// W. T. and Flannery, B. P. (2007). Numerical Recipes, 3rd ed., section 6.2.1. Cambridge University Press. Wilson,
        /// E. B. and Hilferty, M. M. (1931). The distribution of chi-square. Proceedings of the National Academy of Sciences
        /// 17(12), 684-688.
        /// </para>
        /// </remarks>
        internal static double InverseIncompleteRatio(double a, double probability, bool upper, double logGammaA)
        {
            if (double.IsNaN(a) || double.IsNaN(probability) || !(a > 0.0) || probability < 0.0 || probability > 1.0)
                return double.NaN;
            if (probability == 0.0) return upper ? double.PositiveInfinity : 0.0;
            if (probability == 1.0) return upper ? 0.0 : double.PositiveInfinity;
            if (double.IsPositiveInfinity(a)) return double.PositiveInfinity;

            // Invert the smaller tail; 1 - p is exact for p > 1/2.
            if (probability > 0.5)
            {
                probability = 1.0 - probability;
                upper = !upper;
            }
            if (a == 1.0) return upper ? -Math.Log(probability) : -Tools.Log1p(-probability);

            // ln p = k ln 2 + ln m with p = m 2^k and m in (1/sqrt(2), sqrt(2)]; k Ln2High is exact.
            const double Ln2High = 6.93147180369123816490E-01;
            const double Ln2Low = 1.90821492927058770002E-10;
            const double MinimumNormal = 2.2250738585072014E-308;
            int k = 0;
            double scaled = probability;
            if (scaled < MinimumNormal) { scaled *= 18014398509481984.0; k = -54; } // 2^54, exact
            long bits = BitConverter.DoubleToInt64Bits(scaled);
            k += (int)((bits >> 52) & 0x7FF) - 1023;
            double m = BitConverter.Int64BitsToDouble((bits & 0x000FFFFFFFFFFFFFL) | 0x3FF0000000000000L);
            if (m > Tools.Sqrt2) { m *= 0.5; k++; }
            double logTargetHigh = k * Ln2High;
            double logTargetLow = k * Ln2Low + Math.Log(m);
            double logTarget = logTargetHigh + logTargetLow;

            // Starting value
            double x;
            if (!upper)
            {
                double logSmall = (logTarget + LogGammaOnePlus(a)) / a;
                if (logSmall < -708.0) return Math.Exp(logSmall);
                double small = Math.Exp(logSmall);
                if (small < 0.1 * a) x = small;
                else if (a > 1.0) x = WilsonHilferty(a, Normal.StandardZ(probability), small);
                else x = SmallShapeStart(a, probability, 1.0 - probability);
            }
            else if (logTarget < -10.0 && a < -0.1 * logTarget)
            {
                x = Math.Max(1.0, -logTarget - logGammaA);
                for (int i = 0; i < 2; i++)
                {
                    double next = -logTarget - logGammaA + (a - 1.0) * Math.Log(x);
                    if (next > 0.0) x = next;
                }
            }
            else if (a > 1.0) x = WilsonHilferty(a, -Normal.StandardZ(probability), a);
            else x = SmallShapeStart(a, 1.0 - probability, probability);

            // Safeguarded Halley iteration on g = ln T(a, x) - ln p
            double lo = 0.0, hi = double.PositiveInfinity;
            for (int iteration = 0; iteration < 100; iteration++)
            {
                // g from the ratio T / p while both are normal doubles, so it carries the relative precision of T
                // rather than the absolute precision of ln T; from the logarithms below that range.
                double tail = IncompleteRatio(a, x, upper);
                double logTail, g;
                if (tail >= MinimumNormal && probability >= MinimumNormal)
                {
                    logTail = Math.Log(tail);
                    g = Math.Log(tail / probability);
                }
                else
                {
                    logTail = LogIncompleteRatio(a, x, upper);
                    g = (logTail - logTargetHigh) - logTargetLow;
                }
                if (g == 0.0) return x;
                if ((g < 0.0) != upper) lo = x; else hi = x;

                double logDensity = a <= IncompleteDirectShape
                    ? (a - 1.0) * Math.Log(x) - x - logGammaA
                    : LogIncompleteFactor(a, x) - Math.Log(x);
                double hazard = Math.Exp(logDensity - logTail);
                double next = double.NaN;
                if (hazard > 0.0 && !double.IsInfinity(hazard))
                {
                    double t = upper ? -g / hazard : g / hazard;
                    double curvature = (a - 1.0) / x - 1.0 + (upper ? hazard : -hazard);
                    double r = Math.Max(-1.0, Math.Min(1.0, t * curvature));
                    double step = t / (1.0 - 0.5 * r);
                    next = x - step;
                    // Halley's method is cubic, so a step below 1E-8 x leaves an error of order 1E-24 x
                    if (Math.Abs(step) < 1E-8 * x) return next > lo && next < hi ? next : x;
                }

                if (!(next > lo && next < hi))
                {
                    // Bisect the bracket, geometrically once both ends are positive. While the lower end is 0, a
                    // bound below one is squared, which doubles ln x and reaches the smallest double in 11 steps.
                    if (double.IsPositiveInfinity(hi)) next = 2.0 * x;
                    else if (lo > 0.0) next = Math.Sqrt(lo) * Math.Sqrt(hi);
                    else if (hi < 0.5) next = Math.Max(hi * hi, double.Epsilon);
                    else next = 0.5 * hi;
                    // The root lies below the smallest positive double.
                    if (lo == 0.0 && hi == double.Epsilon) return 0.0;
                    if (next == lo || next == hi) return x;
                }
                x = next;
            }
            return x;
        }

        /// <summary>
        /// Returns the Wilson and Hilferty (1931) start a (1 - 1/(9a) + z/(3 sqrt(a)))³ for the inverse of the incomplete
        /// gamma ratio, or <paramref name="fallback"/> where the cube root is not positive.
        /// </summary>
        /// <param name="a">The shape, a &gt; 1.</param>
        /// <param name="z">The standard normal quantile of the lower-tail probability.</param>
        /// <param name="fallback">The start used where the approximation is not positive.</param>
        /// <returns>The starting value.</returns>
        private static double WilsonHilferty(double a, double z, double fallback)
        {
            double w = 1.0 - 1.0 / (9.0 * a) + z / (3.0 * Math.Sqrt(a));
            return w > 0.0 ? a * w * w * w : fallback;
        }

        /// <summary>
        /// Returns the start of Press et al. (2007, Gamma::invgammp) for a &lt;= 1: (u / s)^(1/a) for u &lt; s and
        /// 1 - ln((1 - u) / (1 - s)) above, with s = 1 - a (0.253 + 0.12 a).
        /// </summary>
        /// <param name="a">The shape, 0 &lt; a &lt;= 1.</param>
        /// <param name="lower">The lower-tail probability u.</param>
        /// <param name="complement">The upper-tail probability 1 - u, exact when it is the tail being inverted.</param>
        /// <returns>The starting value.</returns>
        private static double SmallShapeStart(double a, double lower, double complement)
        {
            double s = 1.0 - a * (0.253 + a * 0.12);
            return lower < s ? Math.Pow(lower / s, 1.0 / a) : 1.0 - Math.Log(complement / (1.0 - s));
        }

        #region List Overloads

        /// <summary>
        /// Computes Stirling's approximation of the gamma function for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Stirling(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Stirling(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Stirling(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Stirling(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the Lanczos approximation of the gamma function for each value of a list.
        /// </summary>
        /// <param name="x">The values of x.</param>
        /// <returns>One result per value, each equal to <see cref="Lanczos(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Lanczos(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Lanczos(IList<double> x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = Lanczos(x[i]);
            return result;
        }

        /// <summary>
        /// Computes the gamma function for each value of a list.
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
        /// Computes the digamma function for each value of a list.
        /// </summary>
        /// <param name="X">The values of X.</param>
        /// <returns>One result per value, each equal to <see cref="Digamma(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Digamma(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Digamma(IList<double> X)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            var result = new double[X.Count];
            for (int i = 0; i < X.Count; i++)
                result[i] = Digamma(X[i]);
            return result;
        }

        /// <summary>
        /// Computes the trigamma function for each value of a list.
        /// </summary>
        /// <param name="X">The values of X.</param>
        /// <returns>One result per value, each equal to <see cref="Trigamma(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Trigamma(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Trigamma(IList<double> X)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            var result = new double[X.Count];
            for (int i = 0; i < X.Count; i++)
                result[i] = Trigamma(X[i]);
            return result;
        }

        /// <summary>
        /// Computes the logarithm of the gamma function for each value of a list.
        /// </summary>
        /// <param name="z">The values of z.</param>
        /// <returns>One result per value, each equal to <see cref="LogGamma(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="LogGamma(double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] LogGamma(IList<double> z)
        {
            if (z == null) throw new ArgumentNullException(nameof(z));
            var result = new double[z.Count];
            for (int i = 0; i < z.Count; i++)
                result[i] = LogGamma(z[i]);
            return result;
        }

        /// <summary>
        /// Computes the incomplete gamma integral (AS 239 argument order) for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="X">The values of X, one per result.</param>
        /// <param name="alpha">The values of alpha, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="Incomplete(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Incomplete(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Incomplete(IList<double> X, IList<double> alpha)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            if (alpha == null) throw new ArgumentNullException(nameof(alpha));
            if (alpha.Count != X.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(alpha));
            var result = new double[X.Count];
            for (int i = 0; i < X.Count; i++)
                result[i] = Incomplete(X[i], alpha[i]);
            return result;
        }

        /// <summary>
        /// Computes the regularized upper incomplete gamma function Q(a, x) for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="UpperIncomplete(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="UpperIncomplete(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] UpperIncomplete(IList<double> a, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = UpperIncomplete(a[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes the regularized lower incomplete gamma function P(a, x) for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="LowerIncomplete(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="LowerIncomplete(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] LowerIncomplete(IList<double> a, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = LowerIncomplete(a[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes the inverse of the regularized lower incomplete gamma function for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="y">The values of y, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="InverseLowerIncomplete(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="InverseLowerIncomplete(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] InverseLowerIncomplete(IList<double> a, IList<double> y)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (y.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(y));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = InverseLowerIncomplete(a[i], y[i]);
            return result;
        }

        /// <summary>
        /// Computes the inverse of the regularized upper incomplete gamma function for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="y">The values of y, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="InverseUpperIncomplete(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="InverseUpperIncomplete(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] InverseUpperIncomplete(IList<double> a, IList<double> y)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (y.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(y));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = InverseUpperIncomplete(a[i], y[i]);
            return result;
        }

        #endregion

    }
}
