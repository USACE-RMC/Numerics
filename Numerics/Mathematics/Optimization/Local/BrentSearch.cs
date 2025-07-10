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

using System;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// The Brent optimization algorithm. The function need not be differentiable, and no derivatives are taken.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains methods for finding the minimum or maximum of a function using Brent's method.
    /// Brent's method is a hybrid root-finding algorithm that combines the bisection method, secant method, and inverse
    /// quadratic interpolation.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Brent%27s_method"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class BrentSearch : Optimizer
    {
        /// <summary>
        /// Construct a new Brent optimization method. 
        /// </summary>
        /// <param name="objectiveFunction">The objective function to evaluate.</param>
        /// <param name="lowerBound">The lower bound (inclusive) of the interval containing the optimal point.</param>
        /// <param name="upperBound">The upper bound (inclusive) of the interval containing the optimal point.</param>
        public BrentSearch(Func<double, double> objectiveFunction, double lowerBound, double upperBound) : base((x) => objectiveFunction(x[0]), 1)
        {
            // validate inputs
            if (upperBound < lowerBound)
            {
                throw new ArgumentOutOfRangeException("upperBound", "The upper bound cannot be less than the lower bound.");
            }
            LowerBound = lowerBound;
            UpperBound = upperBound;
        }

        /// <summary>
        /// The lower bound (inclusive) of the interval containing the optimal point. 
        /// </summary>
        public double LowerBound { get; private set; }

        /// <summary>
        /// The upper bound (inclusive) of the interval containing the optimal point.
        /// </summary>
        public double UpperBound { get; private set; }

        /// <inheritdoc/>
        protected override void Optimize()
        { 
            // Define variables
            bool cancel = false;
            double ax = LowerBound, bx = 0.5 * (UpperBound + LowerBound), cx = UpperBound;
            // Golder ratio and a small number which protects against trying to achieve 
            // fractional accuracy for a minimum that happens to be exactly zero.
            double CGOLD = 0.381966d;
            double ZEPS = Tools.DoubleMachineEpsilon * 1.0e-3;
            double a, b, d = 0.0, etemp, fu, fv, fw, fx;
            double p, q, r, tol1, tol2, u, v, w, x, xm;
            double e = 0.0;

            a = (ax < cx ? ax : cx);
            b = (ax > cx ? ax : cx);
            x = w = v = bx;
            fw = fv = fx = Evaluate(new double[] { x }, ref cancel);
            for (int i = 1; i <= MaxIterations; i++)
            {
                xm = 0.5 * (a + b);
                tol2 = 2.0 * (tol1 = RelativeTolerance * Math.Abs(x) + ZEPS);
                // Test for done here
                if (Math.Abs(x - xm) <= tol2 - 0.5d * (b - a))
                {
                    UpdateStatus(OptimizationStatus.Success);
                    return;
                }
                // Construct a trial parabolic fit
                if (Math.Abs(e) > tol1)
                {
                    r = (x - w) * (fx - fv);
                    q = (x - v) * (fx - fw);
                    p = (x - v) * q - (x - w) * r;
                    q = 2.0 * (q - r);
                    if (q > 0.0) p = -p;
                    q = Math.Abs(q);
                    etemp = e;
                    e = d;
                    if (Math.Abs(p) >= Math.Abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                        d = CGOLD * (e = (x >= xm ? a - x : b - x));
                    // The above conditions determine the acceptability of the parabolic fit. 
                    // here we take the golden section step into the larger of the two segments. 
                    else
                    {
                        d = p / q;
                        u = x + d;
                        if (u - a < tol2 || b - u < tol2)
                            d = Tools.Sign(tol1, xm - x);
                    }
                }
                else
                {
                    d = CGOLD * (e = (x >= xm ? a - x : b - x));
                }

                u = (Math.Abs(d) >= tol1 ? x + d : x + Tools.Sign(tol1, d));
                fu = Evaluate(new double[] { u }, ref cancel);
                if (cancel == true) return;
                // This is the one function evaluation per iteration.
                if (fu <= fx)
                {
                    if (u >= x) a = x; else b = x;
                    v = w;
                    w = x;
                    x = u;
                    fv = fw;
                    fw = fx;
                    fx = fu;
                }
                else
                {
                    if (u < x) a = u; else b = u;
                    if (fu <= fw || w == x)
                    {
                        v = w;
                        w = u;
                        fv = fw;
                        fw = fu;
                    }
                    else if (fu <= fv || v == x || v == w)
                    {
                        v = u;
                        fv = fu;
                    }
                }
                // Done with housekeeping. Back for another iteration
            }

            // If we made it to here, the maximum iterations were reached.
            UpdateStatus(OptimizationStatus.MaximumIterationsReached);
        }

        /// <summary>
        /// Bracket the objective function minimum.
        /// </summary>
        /// <param name="s">Starting step size. Default = 1E-2.</param>
        /// <param name="k">Expansion factor. Default = 2.</param>
        public void Bracket(double s = 1E-2, double k = 2d)
        {
            double a = LowerBound, b = a + s;
            double fa = ObjectiveFunction(new double[] { a });
            double fb = ObjectiveFunction(new double[] { b });
            double c, fc, temp;
            if (fb > fa)
            {
                temp = a;
                a = b;
                b = temp;
                temp = fa;
                fa = fb;
                fb = temp;
                s *= -1;
            }
            while (true)
            {
                c = b + s;
                fc = ObjectiveFunction(new double[] { c });
                if (fc > fb) break;
                a = b;          
                b = c;
                fa = fb;
                fb = fc;
            }
            if (a < c)
            {
                LowerBound = a;
                UpperBound = c;
            }
            else
            {
                LowerBound = c;
                UpperBound = a;
            }

        }

    }
}