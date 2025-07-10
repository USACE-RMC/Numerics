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
            return Math.Exp(Gamma.LogGamma(a) + Gamma.LogGamma(b) - Gamma.LogGamma(a + b));
        }

        /// <summary>
        /// The Incomplete (regularized) Beta function Ix(a, b)
        /// </summary>
        /// <remarks>
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
        /// The incomplete beta function evaluated at the given x, given an upper and lower limit
        /// </returns>
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

            y += t + Gamma.LogGamma(aa + bb) - Gamma.LogGamma(aa) - Gamma.LogGamma(bb);
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
        /// <param name="aa"> must be positive </param>
        /// <param name="bb"> must positive </param>
        /// <param name="yy0"> the value to be evaluated by the inverse (such that y = Ix(aa, bb)) </param>
        /// <returns>
        /// The inverse of the incomplete beta integral at the given y
        /// </returns>
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
                x0 = 1.0;
                throw new ArithmeticException("Underflow");
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
                    if (x0 > 0.75)
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
                throw new ArithmeticException("Underflow");

        newt:           
            if (nflg)
            {
                goto done;
            }

            x0 = x;
            lgm = Gamma.LogGamma(a + b) - Gamma.LogGamma(a) - Gamma.LogGamma(b);
            for (i = 0; i < 10; i++)
            {
                if (i != 0)
                    y = Incomplete(a, b, x0);

                d = (a - 1.0) * Math.Log(x0) + (b - 1.0) * Math.Log(1.0 - x0) + lgm;
                if (d < LogMin)
                    throw new ArithmeticException("Underflow");
                d = Math.Exp(d);
                d = (y - y0) / d;
                x = x0;
                x0 -= d;
                if (x0 <= 0.0)
                    throw new ArithmeticException("underflow");
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

            return x0;
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
                t = Gamma.LogGamma(a + b) - Gamma.LogGamma(a) - Gamma.LogGamma(b) + u + Math.Log(s);
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
    }
}