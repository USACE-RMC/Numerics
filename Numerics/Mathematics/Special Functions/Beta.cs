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
    /// A class for Beta functions.
    /// </summary>
    /// <remarks>  
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The beta function is also know as Euler's integral. Denoted with the Greek
    /// capital beta, it has the following equation:
    /// <code>
    ///              1
    ///     B(a,b) = ∫ t^(a-1)*(t-1)^(b-1) dt
    ///              0
    /// </code>
    /// <para>
    /// It can also be written in terms of the Gamma function as:
    /// </para>
    /// <code>
    ///     B(a,b) = Γ(a)*Γ(b) / Γ(a+b)
    ///
    /// </code>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott,
    /// Cambridge University Press, 1991.
    /// https://en.wikipedia.org/wiki/Beta_function
    /// </description></item>
    /// <item><description>
    /// This class uses beta functions from the Accord Math Library
    /// <see href = "http://accord-framework.net" />
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public sealed class Beta
    {
        private const double LogMin = -745.13321910194122d;
        private const double LogMax = 709.782712893384d;
        private const double GammaMax = 171.62437695630271d;

        /// <summary>
        /// The Beta function
        /// </summary>
        /// <remarks>
        /// References:
        /// <list type="bullet">
        /// <item><description> 
        /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott,
        /// Cambridge University Press, 1991.
        ///  </description></item>
        /// </list>
        /// </remarks>
        /// <param name="a">The lower limit </param>
        /// <param name="b">The upper limit </param>
        /// <returns>
        /// The beta function evaluated at the given upper and lower limit
        /// </returns>
        public static double Function(double a, double b)
        {
            return Math.Exp(LogFunction(a, b));
        }

        /// <summary>
        /// Computes the natural logarithm of the beta function, ln B(a, b) = ln Γ(a) + ln Γ(b) - ln Γ(a + b).
        /// </summary>
        /// <remarks>
        /// <para>
        /// For a, b &gt; 0 the value comes from the routines BETALN, ALGDIV and BCORR of algorithm 708 (DiDonato and
        /// Morris, 1992). With a = min(a, b) and b = max(a, b):
        /// </para>
        /// <list type="bullet">
        /// <item><description>
        /// a ≥ 8: the Stirling terms are combined analytically, -(a - 1/2) ln(a / (a + b)) - b ln(1 + a / b) - ln(b)/2 + ln √(2π),
        /// and the difference of the Stirling corrections δ(a) + δ(b) - δ(a + b) is evaluated directly, where
        /// ln Γ(x) = (x - 1/2) ln x - x + ln √(2π) + δ(x).
        /// </description></item>
        /// <item><description>
        /// b ≥ 8 &gt; a: ln B(a, b) = ln Γ(a) + ln(Γ(b) / Γ(a + b)), with the ratio from the same Stirling-correction difference.
        /// </description></item>
        /// <item><description>
        /// Otherwise each argument above 2 is reduced to (1, 2] by the recurrence Γ(x + 1) = xΓ(x), with products accumulated
        /// before the logarithm is taken, and the remaining ln Γ values are evaluated with <see cref="Gamma.LogGamma(double)"/>.
        /// For a &lt; 1 only b is reduced, so ln Γ(b) - ln Γ(a + b) is never taken between two values near ln Γ(8).
        /// </description></item>
        /// </list>
        /// <para>
        /// No ln Γ of a large argument is subtracted, so the result keeps its absolute accuracy when ln Γ(a) and ln Γ(a + b)
        /// are large and nearly equal, for example at a = b = 1000 or at a = 1/2, b = 1E6.
        /// </para>
        /// <para>
        /// When either argument is zero or negative the result is ln|Γ(a)| + ln|Γ(b)| - ln|Γ(a + b)| evaluated with
        /// <see cref="Gamma.LogGamma(double)"/>, as for scipy.special.betaln. A NaN argument returns NaN, and a positive
        /// infinite argument with the other argument positive returns negative infinity.
        /// </para>
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// DiDonato, A. R., and Morris, A. H. (1992). Algorithm 708: Significant digit computation of the incomplete beta
        /// function ratios. ACM Transactions on Mathematical Software, 18(3), 360-373.
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="a">The first shape argument.</param>
        /// <param name="b">The second shape argument.</param>
        /// <returns>The natural logarithm of the beta function evaluated at a and b.</returns>
        public static double LogFunction(double a, double b)
        {
            if (double.IsNaN(a) || double.IsNaN(b)) return double.NaN;
            if (a <= 0d || b <= 0d) return Gamma.LogGamma(a) + Gamma.LogGamma(b) - Gamma.LogGamma(a + b);
            if (double.IsPositiveInfinity(a) || double.IsPositiveInfinity(b)) return double.NegativeInfinity;

            double p = Math.Min(a, b);
            double q = Math.Max(a, b);

            if (p >= 8d)
            {
                double w = LogBetaStirlingCorrection(p, q);
                double h = p / q;
                double u = -(p - 0.5d) * Math.Log(h / (h + 1d));
                double v = q * Tools.Log1p(h);
                return u > v ? -0.5d * Math.Log(q) + Tools.LogSqrt2PI + w - v - u
                             : -0.5d * Math.Log(q) + Tools.LogSqrt2PI + w - u - v;
            }

            double logProduct = 0d;
            if (p < 1d)
            {
                if (q >= 8d) return Gamma.LogGamma(p) + LogGammaRatio(p, q);
                if (q <= 2d) return Gamma.LogGamma(p) + (Gamma.LogGamma(q) - Gamma.LogGamma(p + q));
            }
            else if (p < 2d)
            {
                if (q <= 2d) return Gamma.LogGamma(p) + Gamma.LogGamma(q) - Gamma.LogGamma(p + q);
            }
            else
            {
                int n = (int)(p - 1d);
                double product = 1d;
                if (q > 1000d)
                {
                    // Reduction of p when q > 1000: B(p, q) = (p - 1)/(p - 1 + q) B(p - 1, q)
                    for (int i = 1; i <= n; i++)
                    {
                        p -= 1d;
                        product *= p / (p / q + 1d);
                    }
                    return Math.Log(product) - n * Math.Log(q) + (Gamma.LogGamma(p) + LogGammaRatio(p, q));
                }

                // Reduction of p when q <= 1000
                for (int i = 1; i <= n; i++)
                {
                    p -= 1d;
                    double h = p / q;
                    product *= h / (h + 1d);
                }
                logProduct = Math.Log(product);
            }

            if (q >= 8d) return logProduct + Gamma.LogGamma(p) + LogGammaRatio(p, q);

            // Reduction of q < 8: B(p, q) = (q - 1)/(p + q - 1) B(p, q - 1)
            int m = (int)(q - 1d);
            double z = 1d;
            for (int i = 1; i <= m; i++)
            {
                q -= 1d;
                z *= q / (p + q);
            }
            return logProduct + Math.Log(z) + (Gamma.LogGamma(p) + (Gamma.LogGamma(q) - Gamma.LogGamma(p + q)));
        }

        /// <summary>
        /// Coefficients of the Stirling correction δ(x) = Σ c_k / x^(2k+1), k = 0..5, of algorithm 708, listed from c_5 down to c_0.
        /// </summary>
        private static readonly double[] StirlingCorrection =
        {
            -.00165322962780713, 8.37308034031215E-4, -5.9520293135187E-4,
            7.9365066682539E-4, -.00277777777760991, .0833333333333333
        };

        /// <summary>
        /// Computes δ(b) - δ(a + b) for b ≥ 8 and a positive increment a, where δ is the Stirling correction of ln Γ
        /// (routines ALGDIV and BCORR of algorithm 708).
        /// </summary>
        /// <param name="b">The base argument, at least 8.</param>
        /// <param name="c">The factor a / (a + b), or b / (a + b) when a &gt; b.</param>
        /// <param name="x">The ratio b / (a + b), or a / (a + b) when a &gt; b.</param>
        /// <returns>The difference of the Stirling corrections.</returns>
        /// <remarks>
        /// δ(b) - δ(a + b) = (c / b) Σ c_k s_(2k+1) / b^(2k), where s_n = (1 - x^n) / (1 - x), so the difference is formed
        /// without cancellation.
        /// </remarks>
        private static double StirlingCorrectionDifference(double b, double c, double x)
        {
            double x2 = x * x;
            double s3 = x + x2 + 1d;
            double s5 = x + x2 * s3 + 1d;
            double s7 = x + x2 * s5 + 1d;
            double s9 = x + x2 * s7 + 1d;
            double s11 = x + x2 * s9 + 1d;
            double t = 1d / (b * b);
            double[] k = StirlingCorrection;
            double w = ((((k[0] * s11 * t + k[1] * s9) * t + k[2] * s7) * t + k[3] * s5) * t + k[4] * s3) * t + k[5];
            return w * c / b;
        }

        /// <summary>
        /// Computes ln(Γ(b) / Γ(a + b)) for b ≥ 8 (routine ALGDIV of algorithm 708).
        /// </summary>
        /// <param name="a">The positive increment.</param>
        /// <param name="b">The base argument, at least 8.</param>
        /// <returns>The logarithm of the gamma-function ratio.</returns>
        internal static double LogGammaRatio(double a, double b)
        {
            double h, c, x, d;
            if (a > b)
            {
                h = b / a;
                c = 1d / (h + 1d);
                x = h / (h + 1d);
                d = a + (b - 0.5d);
            }
            else
            {
                h = a / b;
                c = h / (h + 1d);
                x = 1d / (h + 1d);
                d = b + (a - 0.5d);
            }

            double w = StirlingCorrectionDifference(b, c, x);
            double u = d * Tools.Log1p(a / b);
            double v = a * (Math.Log(b) - 1d);
            return u > v ? w - v - u : w - u - v;
        }

        /// <summary>
        /// Computes δ(a) + δ(b) - δ(a + b) for 8 ≤ a ≤ b, where δ is the Stirling correction of ln Γ (routine BCORR of
        /// algorithm 708).
        /// </summary>
        /// <param name="a">The smaller argument, at least 8.</param>
        /// <param name="b">The larger argument.</param>
        /// <returns>The combined Stirling correction.</returns>
        private static double LogBetaStirlingCorrection(double a, double b)
        {
            double h = a / b;
            double w = StirlingCorrectionDifference(b, h / (h + 1d), 1d / (h + 1d));
            double t = 1d / (a * a);
            return Evaluate.PolynomialRev(StirlingCorrection, t) / a + w;
        }

        /// <summary>
        /// The Incomplete (regularized) Beta function Ix(a, b)
        /// </summary>
        /// <remarks>
        /// <para>
        /// This code was copied from the Accord Math Library (the Cephes routine incbet). When one shape is at least 15
        /// and the other at most 1, and the variable of the large shape is at least 0.7, the asymptotic expansion of
        /// algorithm 708 is used instead (<see cref="AsymptoticIncomplete"/>), because the continued fractions converge
        /// too slowly near x = a / (a + b) there. On a 92-point grid (one shape from 30 to 10^4, the
        /// other from 0.01 to 1, and 1 - x from 10^-6 to 0.33), the largest error against mpmath is 2.6E-14 relative. A
        /// small shape below 0.1 with the other shape between 1 and 15 has the same slow region close to x = 1; there the
        /// large shape is first raised to at least 15 by the upward recurrence (routine BUP of algorithm 708).
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="a">The first shape parameter, a &gt; 0.</param>
        /// <param name="b">The second shape parameter, b &gt; 0.</param>
        /// <param name="x">The value to be evaluated, in [0, 1].</param>
        /// <returns>The regularized incomplete beta function Ix(a, b).</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a &lt;= 0, b &lt;= 0, or x is outside [0, 1].</exception>
        public static double Incomplete(double a, double b, double x)
        {
            double aa, bb, t, xx, xc, w, y;
            bool flag;

            // validate inputs
            if (a <= 0.0)
            {
                throw new ArgumentOutOfRangeException("a", "Lower limit a must be greater than zero.");
            }

            if (b <= 0.0)
            {
                throw new ArgumentOutOfRangeException("b", "Upper limit b must be greater than zero.");
            }

            if (x <= 0.0 || x >= 1.0)
            {
                if (x == 0.0)
                {
                    return 0.0;
                }

                if (x == 1.0)
                {
                    return 1.0;
                }

                throw new ArgumentOutOfRangeException("x", "The x value must be between 0 and 1.");
            }

            // One shape at least 15 and the other at most 1, where the variable of the large shape is at least 0.7:
            // the continued fractions converge slowly near x = a / (a + b) there, so the asymptotic expansion is used.
            // A small shape below 0.1 with the other shape between 1 and 15 has the same slow region close to x = 1,
            // and there the large shape is first raised to 15 or more (ShiftedAsymptoticIncomplete).
            double small = Math.Min(a, b), large = Math.Max(a, b);
            if (small <= 1.0 && (large >= AsymptoticMinimumShape || (large > 1.0 && small < 0.1)))
            {
                if (a >= b)
                {
                    double logX = Math.Log(x);
                    if (logX >= AsymptoticMinimumLogVariable) return ShiftedAsymptoticIncomplete(a, b, x, logX);
                }
                else
                {
                    // the expansion gives the complement I(1-x)(b, a); one minus it keeps its accuracy when it is at most 1/2
                    double logComplement = Tools.Log1p(-x);
                    if (logComplement >= AsymptoticMinimumLogVariable)
                    {
                        double complement = ShiftedAsymptoticIncomplete(b, a, 1.0 - x, logComplement);
                        if (complement <= 0.5) return 1.0 - complement;
                    }
                }
            }

            flag = false;
            if (b * x <= 1.0 && x <= 0.95)
            {
                t = PowerSeries(a, b, x);
                return t;
            }

            w = 1.0 - x;
            if (x > a / (a + b))
            {
                flag = true;
                aa = b;
                bb = a;
                xc = x;
                xx = w;
            }
            else
            {
                aa = a;
                bb = b;
                xc = w;
                xx = x;
            }

            if (flag && bb * xx <= 1.0 && xx <= 0.95)
            {
                t = PowerSeries(aa, bb, xx);
                if (t <= Tools.DoubleMachineEpsilon)
                {
                    t = 1.0 - Tools.DoubleMachineEpsilon;
                }
                else
                {
                    t = 1.0 - t;
                }

                return t;
            }

            y = xx * (aa + bb - 2.0) - (aa - 1.0);
            if (y < 0.0)
            {
                w = Incbcf(aa, bb, xx);
            }
            else
            {
                w = Incbd(aa, bb, xx) / xc;
            }

            y = aa * Math.Log(xx);
            t = bb * Math.Log(xc);
            if (aa + bb < GammaMax && Math.Abs(y) < LogMax && Math.Abs(t) < LogMax)
            {
                t = Math.Pow(xc, bb);
                t *= Math.Pow(xx, aa);
                t /= aa;
                t *= w;
                t *= Gamma.Function(aa + bb) / (Gamma.Function(aa) * Gamma.Function(bb));
                if (flag)
                {
                    if (t <= Tools.DoubleMachineEpsilon)
                    {
                        t = 1.0 - Tools.DoubleMachineEpsilon;
                    }
                    else
                    {
                        t = 1.0 - t;
                    }
                }

                return t;
            }

            y += t - LogFunction(aa, bb);
            y += Math.Log(w / aa);
            if (y < LogMin)
            {
                t = 0.0;
            }
            else
            {
                t = Math.Exp(y);
            }

            if (flag)
            {
                if (t <= Tools.DoubleMachineEpsilon)
                {
                    t = 1.0 - Tools.DoubleMachineEpsilon;
                }
                else
                {
                    t = 1.0 - t;
                }
            }

            return t;
        }

        /// <summary>
        /// Returns Ix(a, b) for a large shape a &gt; 1 and a small shape b ≤ 1 with x ≥ 0.7, by the asymptotic expansion,
        /// first raising a below 15 by the upward recurrence.
        /// </summary>
        /// <param name="a">The large shape, above 1.</param>
        /// <param name="b">The small shape, in (0, 1].</param>
        /// <param name="x">The variable, at least 0.7.</param>
        /// <param name="logX">ln x.</param>
        /// <returns>Ix(a, b).</returns>
        /// <remarks>
        /// This is routine BUP of algorithm 708: with n the smallest integer such that a + n ≥ 15,
        /// Ix(a, b) = Ix(a + n, b) + Σ tᵢ, i = 0, ..., n - 1, where t₀ = x^a (1 - x)^b / (a B(a, b)) and
        /// tᵢ₊₁ = tᵢ x (a + b + i) / (a + i + 1). Every term is positive, so nothing cancels, and Ix(a + n, b) comes from
        /// <see cref="AsymptoticIncomplete"/>. For a ≥ 15 the expansion is used directly.
        /// </remarks>
        private static double ShiftedAsymptoticIncomplete(double a, double b, double x, double logX)
        {
            if (a >= AsymptoticMinimumShape) return AsymptoticIncomplete(a, b, logX, 0.0, false);
            int n = (int)Math.Ceiling(AsymptoticMinimumShape - a);
            double term = Math.Exp(a * logX + b * Tools.Log1p(-x) - Math.Log(a) - LogFunction(a, b));
            double sum = 0.0;
            for (int i = 0; i < n; i++)
            {
                sum += term;
                term *= x * (a + b + i) / (a + i + 1.0);
            }
            return sum + AsymptoticIncomplete(a + n, b, logX, 0.0, false);
        }

        /// <summary>
        /// The smallest large shape for which <see cref="AsymptoticIncomplete"/> is used, as in algorithm 708.
        /// </summary>
        internal const double AsymptoticMinimumShape = 15.0;

        /// <summary>
        /// ln 0.7: <see cref="AsymptoticIncomplete"/> is used for x of at least 0.7, as in algorithm 708.
        /// </summary>
        internal const double AsymptoticMinimumLogVariable = -0.35667494393873238;

        /// <summary>
        /// Returns Ix(a, b), or its logarithm, by the asymptotic expansion for a large shape a and a small shape b.
        /// </summary>
        /// <param name="a">The large shape, at least 15.</param>
        /// <param name="b">The small shape, in (0, 1].</param>
        /// <param name="logX">ln x, at least ln 0.7: the leading part of a double-double value.</param>
        /// <param name="logXLow">The trailing part of ln x, or zero.</param>
        /// <param name="logarithm">True to return ln Ix(a, b), which stays finite after Ix(a, b) underflows.</param>
        /// <returns>Ix(a, b), or ln Ix(a, b) when <paramref name="logarithm"/> is true.</returns>
        /// <remarks>
        /// <para>
        /// This is routine BGRAT of algorithm 708. With ν = a + (b - 1)/2 and z = -ν ln x,
        /// Ix(a, b) = Γ(a + b) / (Γ(a) ν^b) [Q(b, z) + e^(-z) z^b / Γ(b) Σ dₙ Jₙ], a series in powers of (ln x / (2ν))²
        /// whose terms Jₙ follow from J₀ = e^z z^(-b) Γ(b) Q(b, z) by a recurrence; up to 30 terms are summed to a
        /// relative 1E-17. The gamma ratio is formed from the Stirling corrections and logarithms of numbers near one,
        /// without subtracting large log-gamma values. For b = 1/2, Q(1/2, z) = erfc(√z) = e^(-z) erfcx(√z), so the
        /// factor e^(-z) is taken out exactly and the trailing part of z, carried in double-double arithmetic, corrects
        /// it; other b use <see cref="Gamma.UpperIncomplete(double, double)"/>.
        /// </para>
        /// <para>
        /// References: A. R. DiDonato and A. H. Morris Jr., "Algorithm 708: Significant digit computation of the
        /// incomplete beta function ratios," ACM Transactions on Mathematical Software, vol. 18, no. 3, pp. 360-373, 1992.
        /// </para>
        /// </remarks>
        internal static double AsymptoticIncomplete(double a, double b, double logX, double logXLow, bool logarithm)
        {
            double bm1 = (b - 0.5) - 0.5;
            double nu = a + 0.5 * bm1;
            // ν is not exact when a exceeds 2^51; its residual, and z = -ν ln x, are carried as double-double values
            double nuLow = (a - nu) + 0.5 * bm1;
            double product = nu * logX;
            double zHigh = -product;
            double zLow = -(ProductError(nu, logX, product) + nu * logXLow + nuLow * logX);
            if (!(zHigh > 0.0)) return logarithm ? 0.0 : 1.0;

            // ln(Γ(a + b) / (Γ(a) ν^b)) = a (ln(1 + h) - h) + (b - 1/2) ln(1 + h) - b ln(1 + (b - 1)/(2a))
            //                              - (δ(a) - δ(a + b)), with h = b/a and δ the Stirling correction
            double h = b / a;
            double logRatio = a * Tools.Log1pMinusX(h) + (b - 0.5) * Tools.Log1p(h)
                              - b * Tools.Log1p(bm1 / (2.0 * a)) - StirlingCorrectionDifference(a, h / (h + 1.0), 1.0 / (h + 1.0));

            bool half = b == 0.5;
            double w = Math.Sqrt(zHigh);
            double q = 0.0, r = 0.0, j0;
            if (half)
            {
                j0 = Erf.SqrtPi * Erf.ScaledErfc(w) / w;
            }
            else
            {
                q = Gamma.UpperIncomplete(b, zHigh);
                r = Math.Exp(-zHigh + b * Math.Log(zHigh) - Gamma.LogGamma(b));
                if (r == 0.0)
                {
                    double value = Math.Exp(logRatio) * q;
                    return logarithm ? Math.Log(value) : value;
                }
                j0 = q / r;
            }

            double v = 0.25 / (nu * nu);
            double t2 = 0.25 * logX * logX;
            Span<double> c = stackalloc double[31];
            Span<double> d = stackalloc double[31];
            double jn = j0, sum = 0.0, t = 1.0, cn = 1.0, n2 = 0.0;
            for (int n = 1; n <= 30; n++)
            {
                double bp2n = b + n2;
                jn = (bp2n * (bp2n + 1.0) * jn + (zHigh + bp2n + 1.0) * t) * v;
                n2 += 2.0;
                t *= t2;
                cn /= n2 * (n2 + 1.0);
                c[n] = cn;
                double s = 0.0;
                double coef = b - n;
                for (int i = 1; i < n; i++)
                {
                    s += coef * c[i] * d[n - i];
                    coef += b;
                }
                d[n] = bm1 * cn + s / n;
                double dj = d[n] * jn;
                sum += dj;
                if (Math.Abs(dj) <= 1E-17 * (j0 + sum)) break;
            }

            if (half)
            {
                // Ix(a, 1/2) = G e^(-z) [erfcx(√z) + √z Σ dₙ Jₙ / √π]
                double e = Erf.ScaledErfc(w) + w * sum / Erf.SqrtPi;
                if (logarithm) return logRatio - zHigh - zLow + Math.Log(e);
                // e^(-z) underflows long before the trailing part of z exceeds one, where 1 - zLow would change sign
                double scale = Math.Exp(-zHigh);
                return scale == 0.0 ? 0.0 : scale * Math.Exp(logRatio) * (1.0 - zLow) * e;
            }

            double result = Math.Exp(logRatio) * (q + r * sum);
            return logarithm ? Math.Log(result) : result;
        }

        /// <summary>
        /// Returns the rounding error a b - p of a floating-point product p = a b, exactly.
        /// </summary>
        /// <param name="a">The first factor.</param>
        /// <param name="b">The second factor.</param>
        /// <param name="product">The rounded product a * b.</param>
        /// <returns>The exact residual a b - <paramref name="product"/> while the product and the partial products neither
        /// overflow nor underflow.</returns>
        /// <remarks>
        /// Dekker's two-product: each factor is split by Veltkamp's method into halves of 26 bits, whose products are
        /// exact. It gives the same value as a fused multiply-add on every target framework. The split multiplies a factor
        /// by 2^27 + 1, so a factor above 2^996 is first scaled by 2^-28, together with the product, and the residual is
        /// scaled back; both scalings are exact.
        /// </remarks>
        internal static double ProductError(double a, double b, double product)
        {
            const double Splitter = 134217729.0; // 2^27 + 1
            const double Large = 6.696928794914171E+299; // 2^996
            const double Down = 3.725290298461914E-09; // 2^-28
            const double Up = 268435456.0; // 2^28
            if (Math.Abs(a) > Large && !double.IsInfinity(a)) return Up * ProductError(a * Down, b, product * Down);
            if (Math.Abs(b) > Large && !double.IsInfinity(b)) return Up * ProductError(a, b * Down, product * Down);
            double ca = Splitter * a, aHigh = ca - (ca - a), aLow = a - aHigh;
            double cb = Splitter * b, bHigh = cb - (cb - b), bLow = b - bHigh;
            return ((aHigh * bHigh - product) + aHigh * bLow + aLow * bHigh) + aLow * bLow;
        }

        /// <summary>
        /// Continued fraction expansion #1 for incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="a">The lower limit </param>
        /// <param name="b">The upper limit </param>
        /// <param name="x">The value to be evaluated </param> 
        /// <returns>
        /// The continue fraction expansion #1 for the incomplete beta integral evaluated at the given x, given a lower and upper limit
        /// </returns>
        public static double Incbcf(double a, double b, double x)
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, thresh;
            double big = 4503599627370496.0;
            double biginv = 0.00000000000000022204460492503131d;
            k1 = a;
            k2 = a + b;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = b - 1.0;
            k7 = k4;
            k8 = a + 2.0;
            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            ans = 1.0;
            r = 1.0;
            thresh = 3.0 * Tools.DoubleMachineEpsilon;
            for (int i = 1; i <= 300; i++)
            {
                xk = -(x * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                xk = x * k5 * k6 / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (qk != 0)
                    r = pk / qk;
                if (r != 0)
                {
                    t = Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                {
                    t = 1.0;
                }

                if (t < thresh)
                    return ans;
                k1 += 1.0;
                k2 += 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 -= 1.0;
                k7 += 2.0;
                k8 += 2.0;
                if (Math.Abs(qk) + Math.Abs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }

                if (Math.Abs(qk) < biginv || Math.Abs(pk) < biginv)
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            }

            return ans;
        }

        /// <summary>
        /// Continued fraction expansion #2 for incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="a">The lower limit </param>
        /// <param name="b">The upper limit </param>
        /// <param name="x">The value to be evaluated </param> 
        /// <returns>
        /// The continue fraction expansion #2 for the incomplete beta integral evaluated at the given x, given a lower and upper limit
        /// </returns>
        public static double Incbd(double a, double b, double x)
        {
            double xk, pk, pkm1, pkm2, qk, qkm1, qkm2;
            double k1, k2, k3, k4, k5, k6, k7, k8;
            double r, t, ans, z, thresh;
            double big = 4503599627370496.0;
            double biginv = 0.00000000000000022204460492503131d;
            k1 = a;
            k2 = b - 1.0;
            k3 = a;
            k4 = a + 1.0;
            k5 = 1.0;
            k6 = a + b;
            k7 = a + 1.0;
            k8 = a + 2.0;
            pkm2 = 0.0;
            qkm2 = 1.0;
            pkm1 = 1.0;
            qkm1 = 1.0;
            z = x / (1.0 - x);
            ans = 1.0;
            r = 1.0;
            thresh = 3.0 * Tools.DoubleMachineEpsilon;
            for (int i = 1; i <= 300; i++)
            {
                xk = -(z * k1 * k2) / (k3 * k4);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                xk = z * k5 * k6 / (k7 * k8);
                pk = pkm1 + pkm2 * xk;
                qk = qkm1 + qkm2 * xk;
                pkm2 = pkm1;
                pkm1 = pk;
                qkm2 = qkm1;
                qkm1 = qk;
                if (qk != 0.0)
                    r = pk / qk;
                if (r != 0.0)
                {
                    t = Math.Abs((ans - r) / r);
                    ans = r;
                }
                else
                {
                    t = 1.0;
                }

                if (t < thresh)
                    return ans;
                k1 += 1.0;
                k2 -= 1.0;
                k3 += 2.0;
                k4 += 2.0;
                k5 += 1.0;
                k6 += 1.0;
                k7 += 2.0;
                k8 += 2.0;
                if (Math.Abs(qk) + Math.Abs(pk) > big)
                {
                    pkm2 *= biginv;
                    pkm1 *= biginv;
                    qkm2 *= biginv;
                    qkm1 *= biginv;
                }

                if (Math.Abs(qk) < biginv || Math.Abs(pk) < biginv)
                {
                    pkm2 *= big;
                    pkm1 *= big;
                    qkm2 *= big;
                    qkm1 *= big;
                }
            }

            return ans;
        }

        /// <summary>
        /// Inverse of incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library. Where the search ends on a target above 1/2 (the complement
        /// of a small tail, either as given or after its bisection reflects the problem past x = 0.75), the root is
        /// refined by Newton's method on the logarithm of the smaller tail, so a small tail probability keeps its relative
        /// accuracy; the bisection does not reflect a target below half an ulp of one, whose complement rounds to one.
        /// A root below e^-708 is the leading series term (y a B(a, b))^(1/a), a subnormal number or 0.
        /// <para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <param name="aa"> must be positive </param>
        /// <param name="bb"> must positive </param>
        /// <param name="yy0"> the value to be evaluated by the inverse (such that y = Ix(aa, bb)) </param>
        /// <returns>
        /// The inverse of the incomplete beta integral at the given y
        /// </returns>
        /// <exception cref="ArithmeticException">The search underflows and the root is not below e^-708.</exception>
        public static double IncompleteInverse(double aa, double bb, double yy0)
        {
            const double LogMin = -745.13321910194122d;
            double a, b, y0, d, y, x, x0, x1, lgm, yp, di, dithresh, yl, yh;
            int i, dir;
            bool nflg;
            bool rflg;

            if (yy0 <= 0)
                return 0.0;
            if (yy0 >= 1.0)
                return 1.0;

            if (aa <= 1.0 || bb <= 1.0)
            {
                nflg = true;
                dithresh = 4.0 * Tools.DoubleMachineEpsilon;
                rflg = false;
                a = aa;
                b = bb;
                y0 = yy0;
                x = a / (a + b);
                y = Incomplete(a, b, x);
                goto ihalve;
            }
            else
            {
                nflg = false;
                dithresh = 0.0001;
            }

            yp = -Normal.StandardZ(yy0);
            if (yy0 > 0.5d)
            {
                rflg = true;
                a = bb;
                b = aa;
                y0 = 1.0 - yy0;
                yp = -yp;
            }
            else
            {
                rflg = false;
                a = aa;
                b = bb;
                y0 = yy0;
            }

            lgm = (yp * yp - 3.0) / 6.0;
            x0 = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
            y = yp * Math.Sqrt(x0 + lgm) / x0 - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (lgm + 5.0 / 6.0 - 2.0 / (3.0 * x0));
            y = 2.0 * y;
            if (y < LogMin)
            {
                return UnderflowRoot(aa, bb, yy0);
            }

            x = a / (a + b * Math.Exp(y));
            y = Incomplete(a, b, x);
            yp = (y - y0) / y0;
            if (Math.Abs(yp) < 0.01d)
            {
                goto newt;
            }

        ihalve:

            x0 = 0.0;
            yl = 0.0;
            x1 = 1.0;
            yh = 1.0;
            di = 0.5;
            dir = 0;
            for (i = 0; i < 400; i++)
            {
                if (i != 0)
                {
                    x = x0 + di * (x1 - x0);
                    if (x == 1.0)
                        x = 1.0 - Tools.DoubleMachineEpsilon;
                    y = Incomplete(a, b, x);
                    yp = (x1 - x0) / (x1 + x0);
                    if (Math.Abs(yp) < dithresh)
                    {
                        x0 = x;
                        goto newt;
                    }
                }

                if (y < y0)
                {
                    x0 = x;
                    yl = y;
                    if (dir < 0)
                    {
                        dir = 0;
                        di = 0.5;
                    }
                    else if (dir > 1)
                    {
                        di = 0.5 * di + 0.5;
                    }
                    else
                    {
                        di = (y0 - y) / (yh - yl);
                    }

                    dir += 1;
                    // The reflected problem has target 1 - y0, which rounds to one, and the bisection then
                    // cycles between the two orientations, when y0 is below half an ulp of one.
                    if (x0 > 0.75 && 1.0 - y0 < 1.0)
                    {
                        if (rflg)
                        {
                            rflg = false;
                            a = aa;
                            b = bb;
                            y0 = yy0;
                        }
                        else
                        {
                            rflg = true;
                            a = bb;
                            b = aa;
                            y0 = 1.0 - yy0;
                        }

                        x = 1.0 - x;
                        y = Incomplete(a, b, x);
                        goto ihalve;
                    }
                }
                else
                {
                    x1 = x;
                    if (rflg && x1 < Tools.DoubleMachineEpsilon)
                    {
                        x0 = 0.0;
                        goto done;
                    }

                    yh = y;
                    if (dir > 0)
                    {
                        dir = 0;
                        di = 0.5;
                    }
                    else if (dir < -1)
                    {
                        di = 0.5 * di;
                    }
                    else
                    {
                        di = (y - y0) / (yh - yl);
                    }

                    dir -= 1;
                }
            }

            if (x0 >= 1.0)
            {
                x0 = 1.0 - Tools.DoubleMachineEpsilon;
                goto done;
            }

            if (x == 0.0)
                return UnderflowRoot(aa, bb, yy0);

        newt:           
            if (nflg)
            {
                goto done;
            }

            x0 = x;
            lgm = -LogFunction(a, b);
            for (i = 0; i < 10; i++)
            {
                if (i != 0)
                    y = Incomplete(a, b, x0);

                d = (a - 1.0) * Math.Log(x0) + (b - 1.0) * Math.Log(1.0 - x0) + lgm;
                if (d < LogMin)
                    return UnderflowRoot(aa, bb, yy0);
                d = Math.Exp(d);
                d = (y - y0) / d;
                x = x0;
                x0 -= d;
                if (x0 <= 0.0)
                    return UnderflowRoot(aa, bb, yy0);
                if (x0 >= 1.0)
                {
                    x0 = 1.0 - Tools.DoubleMachineEpsilon;
                    goto done;
                }

                if (Math.Abs(d / x0) < 64.0 * Tools.DoubleMachineEpsilon)
                {
                    goto done;
                }
            }

        done:
            if (rflg)
            {
                if (x0 <= double.Epsilon)
                {
                    x0 = 1.0 - double.Epsilon;
                }
                else
                {
                    x0 = 1.0 - x0;
                }
            }

            // A target above 1/2 (the complement of a small tail, after a reflection inside the bisection or as given)
            // is resolved only to its absolute accuracy, so the root is refined against the smaller tail.
            return y0 > 0.5 ? RefineIncompleteInverse(aa, bb, yy0, x0) : x0;
        }

        /// <summary>
        /// Returns the root of Ix(a, b) = y where it is below e^-708, from the leading term of the series.
        /// </summary>
        /// <param name="a">The first shape, positive.</param>
        /// <param name="b">The second shape, positive.</param>
        /// <param name="y">The target probability, in (0, 1).</param>
        /// <returns>x = (y a B(a, b))^(1/a), a subnormal number or 0 where it is below the smallest positive double.</returns>
        /// <exception cref="ArithmeticException">The root is not below e^-708, or y &gt; 1/2.</exception>
        /// <remarks>
        /// For x → 0, Ix(a, b) = x^a / (a B(a, b)) (1 + O(x)), so below e^-708 the leading term gives the root to double
        /// precision. The power is formed in logarithms, so it underflows only where the root does.
        /// </remarks>
        private static double UnderflowRoot(double a, double b, double y)
        {
            if (y <= 0.5)
            {
                double logRoot = (Math.Log(y) + Math.Log(a) + LogFunction(a, b)) / a;
                if (logRoot < -708.0) return Math.Exp(logRoot);
            }
            throw new ArithmeticException("Underflow");
        }

        /// <summary>
        /// Refines a root of Ix(a, b) = y by Newton's method on the logarithm of the smaller tail, where its first Newton
        /// correction is at least 64 ε relative, the tolerance of the Newton stage of the inverse, and the tail is more
        /// than 1E-13 relative from the target.
        /// </summary>
        /// <param name="a">The first shape, positive.</param>
        /// <param name="b">The second shape, positive.</param>
        /// <param name="y">The target probability, in (0, 1).</param>
        /// <param name="x">The approximate root, in [0, 1].</param>
        /// <returns><paramref name="x"/> where it already meets either tolerance; otherwise the refined root.</returns>
        /// <remarks>
        /// Once its bisection passes x = 0.75, <see cref="IncompleteInverse(double, double, double)"/> solves the reflected
        /// problem with target 1 - y, so a small tail probability keeps only its absolute accuracy. With T the smaller tail (Ix(a, b) for y ≤ 1/2, its
        /// complement above) and f the beta density, g(x) = ln T(x) - ln T* has g' = ±f/T, and each step
        /// x ← x ∓ g T/f carries the relative accuracy of T. The upper tail is evaluated as I(1-x)(b, a) for x ≥ 1/2,
        /// where 1 - x is exact, or where f &lt; 1, so the rounding of 1 - x moves the tail by less than the rounding of
        /// 1 - Ix(a, b); otherwise as 1 - Ix(a, b). Up to 20 steps are taken, each halved while it would leave
        /// (0, 1); the iteration stops once either tolerance is met, applying the last step when it is below 64 ε x.
        /// </remarks>
        private static double RefineIncompleteInverse(double a, double b, double y, double x)
        {
            if (!(x > 0.0 && x < 1.0)) return x;
            bool lower = y <= 0.5;
            double logTarget = lower ? Math.Log(y) : Tools.Log1p(-y);
            double logBeta = LogFunction(a, b);
            double current = x;
            for (int i = 0; i < 20; i++)
            {
                double logDensity = (a - 1.0) * Math.Log(current) + (b - 1.0) * Tools.Log1p(-current) - logBeta;
                double tail = lower ? Incomplete(a, b, current)
                                    : current >= 0.5 || logDensity < 0.0 ? Incomplete(b, a, 1.0 - current) : 1.0 - Incomplete(a, b, current);
                if (!(tail > 0.0)) return i == 0 ? x : current;
                double g = Math.Log(tail) - logTarget;
                double step = g * Math.Exp(Math.Log(tail) - logDensity);
                if (!lower) step = -step;
                if (double.IsNaN(step) || double.IsInfinity(step)) return i == 0 ? x : current;
                bool smallStep = Math.Abs(step) < 64.0 * Tools.DoubleMachineEpsilon * current;
                if (smallStep || Math.Abs(g) <= 1E-13)
                {
                    if (i == 0) return x;
                    return smallStep ? current - step : current;
                }
                double next = current - step;
                while (!(next > 0.0 && next < 1.0))
                {
                    step *= 0.5;
                    next = current - step;
                }
                current = next;
            }
            return current;
        }

        /// <summary>
        /// Power series for incomplete beta integral
        /// </summary>
        /// <remarks>
        /// Use when b*x is small and x not too close to 1.
        /// <para>
        /// This code was copied from the Accord Math Library.
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="a">The lower limit </param>
        /// <param name="b">The upper limit </param>
        /// <param name="x">The value to be evaluated </param> 
        /// <returns>
        /// The power series for the incomplete beta integral evaluated at the given x, given an upper and lower limit
        /// </returns>
        public static double PowerSeries(double a, double b, double x)
        {
            double s, t, u, v, n, t1, z, ai;
            ai = 1.0 / a;
            u = (1.0 - b) * x;
            v = u / (a + 1.0);
            t1 = v;
            t = u;
            n = 2.0;
            s = 0.0;
            z = Tools.DoubleMachineEpsilon * ai;
            while (Math.Abs(v) > z)
            {
                u = (n - b) * x / n;
                t *= u;
                v = t / (a + n);
                s += v;
                n += 1.0;
            }

            s += t1;
            s += ai;
            u = a * Math.Log(x);
            if (a + b < GammaMax && Math.Abs(u) < LogMax)
            {
                t = Gamma.Function(a + b) / (Gamma.Function(a) * Gamma.Function(b));
                s = s * t * Math.Pow(x, a);
            }
            else
            {
                t = -LogFunction(a, b) + u + Math.Log(s);
                if (t < LogMin)
                {
                    s = 0.0;
                }
                else
                {
                    s = Math.Exp(t);
                }
            }

            return s;
        }

        /// <summary>
        /// Computes incomplete beta function ratio
        /// </summary>
        /// <remarks>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Algorithm as 63 applied statistics (1973), vol.22, no.3.
        /// <see href = "http://people.sc.fsu.edu/~jburkardt/f77_src/asa243/asa243.f" />
        /// </description></item>
        /// </list>
        /// </remarks>
        /// <param name="x">The argument, between 0 and 1 </param>
        /// <param name="p"> must be positive </param>
        /// <param name="q"> must be positive </param>
        /// <param name="beta">The logarithm of the complete beta function </param> 
        /// <returns>
        /// The incomplete beta function ratio (incomplete / complete)
        /// </returns>
        public static double IncompleteRatio(double x, double p, double q, double beta)
        {

            // Define accuracy and initialize parameters
            double acu = 0.000000000000001d;
            bool indx;
            double psq;
            double cx;
            double xx;
            double pp;
            double qq;
            double term;
            double ai;
            int ns;
            double rx;
            double temp;
            double b;

            // Check the input arguments.
            if (p <= 0.0 | q <= 0.0)
            {
                throw new Exception("Invalid inputs! P and Q must be positive.");
            }

            if (x < 0.0 | x > 1.0)
            {
                throw new Exception("Invalid inputs! X must be between 0 and 1.");
            }

            // Special cases.
            if (x == 0.0 | x == 1.0)
            {
                return x;
            }

            // Change tail if necessary and determine S.
            psq = p + q;
            cx = 1.0 - x;
            if (p < psq * x)
            {
                xx = cx;
                cx = x;
                pp = q;
                qq = p;
                indx = true;
            }
            else
            {
                xx = x;
                pp = p;
                qq = q;
                indx = false;
            }

            term = 1.0;
            ai = 1.0;
            b = 1.0;
            ns = (int)(qq + cx * psq);

            // Use the Soper reduction formula.
            rx = xx / cx;
            temp = qq - ai;
            if (ns == 0) rx = xx;
            rx = xx / cx;
        three:
            ;
            temp = qq - ai;
            if (ns == 0) rx = xx;
        four:
            ;
            term = term * temp * rx / (pp + ai);
            b += term;
            temp = Math.Abs(term);
            if (temp <= acu & temp <= acu * b)
                goto five;
            ai += 1;
            ns -= 1;
            if (ns >= 0)
                goto three;
            temp = psq;
            psq += 1d;
            goto four;

        // calculate result
        five:
            ;
            b = b * Math.Exp(pp * Math.Log(xx) + (qq - 1.0) * Math.Log(cx) - beta) / pp;
            if (indx)
                b = 1.0 - b;
            return b;
        }

        #region List Overloads

        /// <summary>
        /// Computes the beta function for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="Function(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Function(double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Function(IList<double> a, IList<double> b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = Function(a[i], b[i]);
            return result;
        }

        /// <summary>
        /// Computes the natural logarithm of the beta function for each set of corresponding values of 2 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="LogFunction(double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="LogFunction(double, double)"/>.
        /// </remarks>
        public static double[] LogFunction(IList<double> a, IList<double> b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = LogFunction(a[i], b[i]);
            return result;
        }

        /// <summary>
        /// Computes the regularized incomplete beta integral for each set of corresponding values of 3 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="Incomplete(double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Incomplete(double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Incomplete(IList<double> a, IList<double> b, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = Incomplete(a[i], b[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes continued fraction expansion 1 of the incomplete beta integral for each set of corresponding values of 3 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="Incbcf(double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Incbcf(double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Incbcf(IList<double> a, IList<double> b, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = Incbcf(a[i], b[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes continued fraction expansion 2 of the incomplete beta integral for each set of corresponding values of 3 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="Incbd(double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="Incbd(double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] Incbd(IList<double> a, IList<double> b, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = Incbd(a[i], b[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes the inverse of the regularized incomplete beta integral for each set of corresponding values of 3 lists.
        /// </summary>
        /// <param name="aa">The values of aa, one per result.</param>
        /// <param name="bb">The values of bb, one per result.</param>
        /// <param name="yy0">The values of yy0, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="IncompleteInverse(double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="IncompleteInverse(double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] IncompleteInverse(IList<double> aa, IList<double> bb, IList<double> yy0)
        {
            if (aa == null) throw new ArgumentNullException(nameof(aa));
            if (bb == null) throw new ArgumentNullException(nameof(bb));
            if (yy0 == null) throw new ArgumentNullException(nameof(yy0));
            if (bb.Count != aa.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(bb));
            if (yy0.Count != aa.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(yy0));
            var result = new double[aa.Count];
            for (int i = 0; i < aa.Count; i++)
                result[i] = IncompleteInverse(aa[i], bb[i], yy0[i]);
            return result;
        }

        /// <summary>
        /// Computes the power series of the regularized incomplete beta integral for each set of corresponding values of 3 lists.
        /// </summary>
        /// <param name="a">The values of a, one per result.</param>
        /// <param name="b">The values of b, one per result.</param>
        /// <param name="x">The values of x, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="PowerSeries(double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="PowerSeries(double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] PowerSeries(IList<double> a, IList<double> b, IList<double> x)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (b.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(b));
            if (x.Count != a.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(x));
            var result = new double[a.Count];
            for (int i = 0; i < a.Count; i++)
                result[i] = PowerSeries(a[i], b[i], x[i]);
            return result;
        }

        /// <summary>
        /// Computes the incomplete beta function ratio (AS 63) for each set of corresponding values of 4 lists.
        /// </summary>
        /// <param name="x">The values of x, one per result.</param>
        /// <param name="p">The values of p, one per result.</param>
        /// <param name="q">The values of q, one per result.</param>
        /// <param name="beta">The values of beta, one per result.</param>
        /// <returns>One result per value, each equal to <see cref="IncompleteRatio(double, double, double, double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when a list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the lists do not have the same number of values.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="Tools.Expm1(IList{double})"/>; a caller
        /// crossing a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// The scalar member is called for each value, so the results are bit for bit those of <see cref="IncompleteRatio(double, double, double, double)"/>,
        /// and an exception it throws for a value propagates unchanged.
        /// </remarks>
        public static double[] IncompleteRatio(IList<double> x, IList<double> p, IList<double> q, IList<double> beta)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (p == null) throw new ArgumentNullException(nameof(p));
            if (q == null) throw new ArgumentNullException(nameof(q));
            if (beta == null) throw new ArgumentNullException(nameof(beta));
            if (p.Count != x.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(p));
            if (q.Count != x.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(q));
            if (beta.Count != x.Count)
                throw new ArgumentException("The lists must have the same number of values.", nameof(beta));
            var result = new double[x.Count];
            for (int i = 0; i < x.Count; i++)
                result[i] = IncompleteRatio(x[i], p[i], q[i], beta[i]);
            return result;
        }

        #endregion

    }
}