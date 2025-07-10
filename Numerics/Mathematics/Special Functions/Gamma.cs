/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
    ///     has simple poles), yielding the meromorphic function we call the gamma
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
        /// </para>
        /// </remarks>
        /// <param name="x"> The value to be evaluated </param>
        /// <returns>
        /// The Gamma function evaluated at the given x
        /// </returns>
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

            if ((x == 2.0) || (x == 3.0))
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
        public static double Trigamma(double X)
        {
            double A = 0.0001d;
            double B = 5.0d;
            double B2 = 0.1666666667d;
            double B4 = -0.03333333333d;
            double B6 = 0.02380952381d;
            double B8 = -0.03333333333d;
            double Y;
            double Z;
            double _trigamma = 0.0d;
            if (X <= 0d)
            {
                throw new ArgumentOutOfRangeException("X", "X must be greater than zero.");
            }

            Z = X;
            if (X <= A)
            {
                // Use small value approximation if X <= A
                _trigamma = 1.0d / X / X;
            }
            else
            {
                // Increase argument to ( X + I ) >= B
                _trigamma = 0.0d;
                while (Z < B)
                {
                    _trigamma = _trigamma + 1.0d / Z / Z;
                    Z += 1.0d;
                }
                // Apply asymptotic formula if argument is B or greater
                Y = 1.0d / Z / Z;
                _trigamma = _trigamma + 0.5d * Y + (1.0d + Y * (B2 + Y * (B4 + Y * (B6 + Y * B8)))) / Z;
            }

            return _trigamma;
        }

        /// <summary>
        /// Natural logarithm of the Gamma function.
        /// </summary>
        /// <remarks>
        /// References: 
        /// <list type="bullet">
        /// <item><description>
        /// Based algorithm ACM291, Commun. Assoc. Comput. Mach. (1966)
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="X"> The value to be evaluated </param>
        /// <returns>
        /// The natural logarithm of the Gamma function evaluated at the given x
        /// </returns>
        public static double LogGamma(double X)
        {
            double SMALL = 0.0000001d;
            double CRIT = 13d;
            double BIG = 1000000000.0d;
            double TOOBIG = 2.0E+36d;
            double C0 = 0.5d * Math.Log(2d * Math.PI);
            double C1 = 0.083333333333333329d;
            double C2 = -0.0027777777777777779d;
            double C3 = 0.00079365079365079365d;
            double C4 = -0.00059523809523809529d;
            double C5 = 0.00084175084175084171d;
            double C6 = -0.0019175269175269176d;
            double C7 = 0.00641025641025641d;
            double S1 = -0.57721566490153287d;
            double S2 = 0.8224670334241132d;
            double _logGamma = 0.0d;
            double XX;
            double Y;
            double Z;
            double SUM1;
            double SUM2;
            if (X <= 0.0d)
            {
                throw new ArgumentOutOfRangeException("X", "X must be greater than zero.");
            }

            if (X > TOOBIG)
            {
                throw new ArgumentOutOfRangeException("X", "X is too big. It must be less than 2E36.");
            }

            if (Math.Abs(X - 2.0d) > SMALL)
            {
                if (Math.Abs(X - 1.0d) > SMALL)
                {
                    if (X > SMALL)
                    {
                        // reduce to LogGamma(X+N) where X+N >= CRIT
                        SUM1 = 0.0d;
                        Y = X;
                        if (Y < CRIT)
                        {
                            Z = 1.0d;
                            while (Y < CRIT)
                            {
                                Z = Z * Y;
                                Y = Y + 1.0d;
                            }

                            SUM1 = SUM1 - Math.Log(Z);
                        }

                        // use asymptotic expansion if Y >= CRIT
                        SUM1 = SUM1 + (Y - 0.5d) * Math.Log(Y) - Y + C0;
                        SUM2 = 0.0d;
                        if (Y >= BIG)
                        {
                            _logGamma = SUM1 + SUM2;
                        }
                        else
                        {
                            Z = 1.0d / (Y * Y);
                            SUM2 = ((((((C7 * Z + C6) * Z + C5) * Z + C4) * Z + C3) * Z + C2) * Z + C1) / Y;
                            _logGamma = SUM1 + SUM2;
                        }
                    }
                    else
                    {
                        _logGamma = -Math.Log(X) + S1 * X;
                    }
                }
                else
                {
                    XX = X - 1.0d;
                    _logGamma = XX * (S1 + XX * S2);
                }
            }
            else
            {
                XX = X - 2.0d;
                _logGamma = Math.Log(X - 1.0d) + XX * (S1 + XX * S2);
            }

            return _logGamma;
        }

        /// <summary>
        /// The incomplete gamma integral.
        /// </summary>
        /// <remarks>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Based on algorithm AS239, Applied Statistics, 1988, Vol. 37, No. 3.
        /// </description></item>
        /// <item><description>
        /// N.L. Johnson And S. Kotz, "Continuous Univariate Distributions 1", P.180
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="X">Argument of function (upper limit of integration).</param>
        /// <param name="alpha">The shape parameter.</param>
        /// <returns>
        /// The incomplete gamma integral evaluated with the given upper limit x
        /// </returns>
        public static double Incomplete(double X, double alpha)
        {
            double X13 = 13.0d;
            double X36 = 36.0d;
            double X42 = 42.0d;
            double X119 = 119.0d;
            double X1620 = 1620.0d;
            double X38880 = 38880.0d;
            double RTHALF = 0.70710678118654757d;
            double EPS = 0.000000000001d;
            int MAXIT = 100000;
            double OFL = 1.0E+30d;
            double UFL = -180.0d;
            double AHILL = 10000.0d;
            double _gammaind = 0.0d;
            if (alpha <= 0.0d)
            {
                throw new ArgumentOutOfRangeException("alpha", "Shape parameter (alpha) out of range.");
            }

            if (X < 0.0d)
            {
                throw new ArgumentOutOfRangeException("X", "Argument of function (X) out of range.");
            }

            if (X == 0d)
            {
                return _gammaind;
            }

            if (alpha > AHILL)
            {

                // alpha is large, so use Hill's approximation 
                // (N.L. Johnson And S. Kotz, "Continuous Univariate Distributions 1", P.180)
                // The for-Loop calculates 2*(X-ALPHA-ALPHA*DLOG(X/ALPHA)),
                // using power-series expansion to avoid rounding error

                // variables
                double SUM;
                double TERM;
                double R;
                double Z;
                double WW;
                double W;
                double H1;
                double H2;
                double H3;
                double H4;
                R = 1.0d / Math.Sqrt(alpha);
                Z = (X - alpha) * R;
                TERM = Z * Z;
                SUM = 0.5d * TERM;
                for (short i = 1; i <= 12; i++)
                {
                    TERM = -TERM * Z * R;
                    SUM = SUM + TERM / (i + 2.0d);
                    if (TERM <= EPS)
                    {
                        break;
                    }
                }

                WW = 2.0d * SUM;
                W = Math.Sqrt(WW);
                if (X < alpha)
                {
                    W = -W;
                }

                H1 = 1.0d / 3.0d;
                H2 = -W / X36;
                H3 = (-WW + X13) / X1620;
                H4 = (X42 * WW + X119) * W / X38880;
                Z = (((H4 * R + H3) * R + H2) * R + H1) * R + W;
                _gammaind = 0.5d + 0.5d * Erf.Function(Z * RTHALF);
            }
            else if (X > 1.0d && X >= alpha)
            {

                // continued-fraction expansion

                // variables
                double TERM;
                double ARG;
                double A;
                double B;
                double PN1;
                double PN2;
                double PN3;
                double PN4;
                double RATIO;
                double AN;
                double PN5;
                double PN6;
                double RN;
                double DIFF;
                A = 1.0d - alpha;
                B = A + X + 1.0d;
                TERM = 0.0d;
                PN1 = 1.0d;
                PN2 = X;
                PN3 = X + 1.0d;
                PN4 = X * B;
                RATIO = PN3 / PN4;
                for (int i = 1; i <= MAXIT; i++)
                {
                    A += 1.0d;
                    B += 2.0d;
                    TERM += 1.0d;
                    AN = A * TERM;
                    PN5 = B * PN3 - AN * PN1;
                    PN6 = B * PN4 - AN * PN2;
                    if (PN6 != 0.0d)
                    {
                        RN = PN5 / PN6;
                        DIFF = Math.Abs(RATIO - RN);
                        if (DIFF <= EPS && DIFF <= EPS * RN)
                        {
                            break;
                        }
                        else
                        {
                            RATIO = RN;
                        }
                    }

                    PN1 = PN3;
                    PN2 = PN4;
                    PN3 = PN5;
                    PN4 = PN6;
                    if (Math.Abs(PN5) < OFL)
                    {
                        // ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE.'
                        // Consider adding a message box or something
                        break;
                    }
                    else
                    {
                        PN1 = PN1 / OFL;
                        PN2 = PN2 / OFL;
                        PN3 = PN3 / OFL;
                        PN4 = PN4 / OFL;
                    }

                    if (i == MAXIT)
                    {
                        // ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE.'
                        // Consider adding a message box or something
                    }
                }

                ARG = alpha * Math.Log(X) - X - LogGamma(alpha) + Math.Log(RATIO);
                _gammaind = 1.0d;
                if (ARG >= UFL)
                {
                    _gammaind = 1.0d - Math.Exp(ARG);
                }
            }
            else
            {

                // series expansion

                // variables
                double SUM;
                double TERM;
                double ARG;
                double A;
                SUM = 1.0d;
                TERM = 1.0d;
                A = alpha;
                for (int i = 1; i <= MAXIT; i++)
                {
                    A += 1.0d;
                    TERM = TERM * X / A;
                    SUM = SUM + TERM;
                    if (TERM <= EPS)
                    {
                        break;
                    }

                    if (i == MAXIT)
                    {
                        // ITERATION HAS NOT CONVERGED. RESULT MAY BE UNRELIABLE.'
                        // Consider adding a message box or something
                    }
                }

                ARG = alpha * Math.Log(X) - X - LogGamma(alpha) + Math.Log(SUM / alpha);
                _gammaind = 0.0d;
                if (ARG >= UFL)
                {
                    _gammaind = Math.Exp(ARG);
                }
            }

            return _gammaind;
        }

        /// <summary>
        /// Upper incomplete regularized Gamma function Q
        /// (a.k.a the incomplete complemented Gamma function)
        /// </summary>
        /// <remarks>
        /// This function is equivalent to Q(x) = Γ(a, x) / Γ(s).
        /// </remarks>
        /// <param name="x"> The value to be evaluated (lower limit of integration) </param>
        /// <param name="a">  complex parameter, such that the real part of a is positive </param>
        /// <returns>
        /// The upper incomplete Gamma function evaluated with the lower limit x
        /// </returns>
        public static double UpperIncomplete(double a, double x)
        {
            const double LogMax = 709.782712893384d;
            const double DoubleEpsilon = 0.000000000000000111022302462516d;
            const double big = 4.5035996273705E+15d;
            const double biginv = 0.000000000000000222044604925031d;
            double ans;
            double ax;
            double c;
            double yc;
            double r;
            double t;
            double y;
            double z;
            double pk;
            double pkm1;
            double pkm2;
            double qk;
            double qkm1;
            double qkm2;
            if (x <= 0d || a <= 0d)
            {
                return 1.0d;
            }

            if (x < 1.0d || x < a)
            {
                return 1.0d - LowerIncomplete(a, x);
            }

            if (double.IsPositiveInfinity(x))
            {
                return 0d;
            }

            ax = a * Math.Log(x) - x - LogGamma(a);
            if (ax < -LogMax)
            {
                return 0.0d;
            }

            ax = Math.Exp(ax);

            // continued fraction
            y = 1.0d - a;
            z = x + y + 1.0d;
            c = 0.0d;
            pkm2 = 1.0d;
            qkm2 = x;
            pkm1 = x + 1.0d;
            qkm1 = z * x;
            ans = pkm1 / qkm1;
            do
            {
                c += 1.0d;
                y += 1.0d;
                z += 2.0d;
                yc = y * c;
                pk = pkm1 * z - pkm2 * yc;
                qk = qkm1 * z - qkm2 * yc;
                if (qk != 0d)
                {
                    r = pk / qk;
                    t = Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                {
                    t = 1.0d;
                }

                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (Math.Abs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }
            }
            while (t > DoubleEpsilon);
            return ans * ax;
        }

        /// <summary>
        /// Lower incomplete regularized gamma function P
        /// (a.k.a. the incomplete Gamma function).
        /// </summary>
        /// <remarks>
        /// This function is equivalent to P(x) = γ(a, x) / Γ(s).
        /// </remarks>
        /// <param name="a">  complex parameter, such that the real part of a is positive </param>
        /// <param name="x"> The value to be evaluated (upper limit of integration) </param>
        /// <returns>
        /// The lower incomplete Gamma function evaluated with the upper limit x
        /// </returns>
        public static double LowerIncomplete(double a, double x)
        {
            const double LogMax = 709.782712893384d;
            const double DoubleEpsilon = 0.000000000000000111022302462516d;
            if (a <= 0d)
            {
                return 1.0d;
            }

            if (x <= 0d)
            {
                return 0.0d;
            }

            if (x > 1.0d && x > a)
            {
                return 1.0d - UpperIncomplete(a, x);
            }

            double ax = a * Math.Log(x) - x - LogGamma(a);
            if (ax < -LogMax)
            {
                return 0.0d;
            }

            ax = Math.Exp(ax);
            double r = a;
            double c = 1.0d;
            double ans = 1.0d;
            do
            {
                r += 1.0d;
                c *= x / r;
                ans += c;
            }
            while (c / ans > DoubleEpsilon);
            return ans * ax / a;
        }

        /// <summary>
        /// Inverse of the <see cref="LowerIncomplete">
        /// incomplete Gamma integral (LowerIncomplete, P)</see>.
        /// </summary>
        /// <param name="a">  complex parameter, such that the real part of a is positive </param>
        /// <param name="y">  the value to evaluate the inverse at </param>
        /// <returns>
        /// The inverse of the lower incomplete Gamma function evaluated at y
        /// </returns>
        public static double InverseLowerIncomplete(double a, double y)
        {
            return Inverse(a, 1d - y);
        }

        /// <summary>
        /// Inverse of the <see cref="UpperIncomplete">complemented
        /// incomplete Gamma integral (UpperIncomplete, Q)</see>.
        /// </summary>
        /// /// <param name="a">  complex parameter, such that the real part of a is positive </param>
        /// <param name="y"> the value to evaluate the inverse at </param>
        /// <returns>
        /// The inverse of the upper incomplete Gamma function evaluated at y
        /// </returns>
        public static double InverseUpperIncomplete(double a, double y)
        {
            return Inverse(a, y);
        }

        /// <summary>
        /// Inverse of the <see cref="Function">complemented
        /// Gamma function (Function, P)</see>.
        /// </summary>
        /// <param name="a">  complex parameter, such that the real part of a is positive </param>
        /// <param name="y"> the value to evaluate the inverse at </param>
        /// <returns>
        /// The inverse of the Gamma function evaluated at y
        /// </returns>
        private static double Inverse(double a, double y)
        {
            const double LogMax = 709.782712893384d;
            const double DoubleEpsilon = 0.000000000000000111022302462516d;

            // bound the solution
            double x0 = double.MaxValue;
            double yl = 0d;
            double x1 = 0d;
            double yh = 1.0d;
            double dithresh = 5.0d * DoubleEpsilon;

            // approximation to inverse function
            double d = 1.0d / (9.0d * a);
            double yy = 1.0d - d - Normal.StandardZ(y) * Math.Sqrt(d);
            double x = a * yy * yy * yy;
            double lgm = LogGamma(a);
            for (int i = 0; i <= 9; i++)
            {
                if (x > x0 || x < x1)
                {
                    goto ihalve;
                }

                yy = UpperIncomplete(a, x);
                if (yy < yl || yy > yh)
                {
                    goto ihalve;
                }

                if (yy < y)
                {
                    x0 = x;
                    yl = yy;
                }
                else
                {
                    x1 = x;
                    yh = yy;
                }

                // compute the derivative of the function at this point
                d = (a - 1.0d) * Math.Log(x) - x - lgm;
                if (d < -LogMax)
                {
                    goto ihalve;
                }

                d = -Math.Exp(d);

                // compute the step to the next approximation of x
                d = (yy - y) / d;
                if (Math.Abs(d / x) < DoubleEpsilon)
                {
                    return x;
                }

                x = x - d;
            }

        ihalve:          
            // Resort to interval halving if Newton iteration did not converge. 

            d = 0.0625d;
            if (x0 == double.MaxValue)
            {
                if (x <= 0.0d)
                {
                    x = 1.0d;
                }

                while (x0 == double.MaxValue && !double.IsNaN(x) && !double.IsInfinity(x))
                {
                    x = (1.0d + d) * x;
                    yy = UpperIncomplete(a, x);
                    if (yy < y)
                    {
                        x0 = x;
                        yl = yy;
                        break;
                    }

                    d = d + d;
                }
            }

            d = 0.5d;
            double dir = 0d;
            for (int i = 0; i <= 399; i++)
            {
                double t = x1 + d * (x0 - x1);
                if (double.IsNaN(t))
                {
                    break;
                }

                x = t;
                yy = UpperIncomplete(a, x);
                lgm = (x0 - x1) / (x1 + x0);
                if (Math.Abs(lgm) < dithresh)
                {
                    break;
                }

                lgm = (yy - y) / y;
                if (Math.Abs(lgm) < dithresh)
                {
                    break;
                }

                if (x <= 0.0d)
                {
                    break;
                }

                if (yy >= y)
                {
                    x1 = x;
                    yh = yy;
                    if (dir < 0d)
                    {
                        dir = 0d;
                        d = 0.5d;
                    }
                    else if (dir > 1d)
                    {
                        d = 0.5d * d + 0.5d;
                    }
                    else
                    {
                        d = (y - yl) / (yh - yl);
                    }

                    dir += 1d;
                }
                else
                {
                    x0 = x;
                    yl = yy;
                    if (dir > 0d)
                    {
                        dir = 0d;
                        d = 0.5d;
                    }
                    else if (dir < -1)
                    {
                        d = 0.5d * d;
                    }
                    else
                    {
                        d = (y - yl) / (yh - yl);
                    }

                    dir -= 1d;
                }
            }

            if (x == 0.0d || double.IsNaN(x))
            {
                throw new ArithmeticException();
            }

            return x;
        }
    }
}