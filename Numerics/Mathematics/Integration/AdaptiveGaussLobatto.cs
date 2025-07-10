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

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class that performs adaptive integration by the Gauss-Lobatto method with a Kronrod extension. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description>
    /// "Numerical Recipes: The art of Scientific Computing, Third Edition. Press et al. 2017.
    /// </description></item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class AdaptiveGaussLobatto : Integrator
    {

        /// <summary>
        /// Constructs a new adaptive Gauss-Lobatto method.
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="min">The minimum value under which the integral must be computed.</param>
        /// <param name="max">The maximum value under which the integral must be computed.</param>
        public AdaptiveGaussLobatto(Func<double, double> function, double min, double max)
        {
            if (max <= min) throw new ArgumentNullException(nameof(max), "The maximum value cannot be less than or equal to the minimum value.");
            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            a = min;
            b = max;
        }

        private double a;
        private double b;

        /// <summary>
        /// The unidimensional function to integrate.
        /// </summary>
        public Func<double, double> Function { get; }

        /// <summary>
        /// The minimum value under which the integral must be computed.
        /// </summary>
        public double Min => a;

        /// <summary>
        /// The maximum value under which the integral must be computed. 
        /// </summary>
        public double Max => b;

        // Abscissas for Gauss-Lobatto-Kronrod quadrature
        private static double alpha = Math.Sqrt(2.0 / 3.0);
        private static double beta = 1.0 / Math.Sqrt(5.0);
        private static double x1 = 0.942882415695480;
        private static double x2 = 0.641853342345781;
        private static double x3 = 0.236383199662150;
        private static double[] x = new double[] { 0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1 };

        // convergence variables
        private static bool terminate = true;
        private static bool outOfTolerance = false;
        private double toler;

        /// <inheritdoc/>
        public override void Integrate()
        {
            ClearResults();
            Validate();

            try
            {
                double m, h, fa, fb, i1, i2, iS, erri1, erri2, r;
                var y = new double[13];
                m = 0.5 * (a + b);
                h = 0.5 * (b - a);
                fa = y[0] = Function(a);
                fb = y[12] = Function(b);
                FunctionEvaluations += 2;
                for (int i = 1; i < 12; i++)
                {
                    y[i] = Function(m + x[i] * h);
                    FunctionEvaluations++;
                }
                    
                // 4-point Gauss-Lobatto formula
                i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));
                // 7-point Kronrod extension
                i1 = (h / 1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]);
                // 13- point Kronrod extension
                iS= h * (0.0158271919734802 * (y[0] + y[12]) + 0.0942738402188500 * (y[1] + y[11]) + 0.155071987336585 * (y[2] + y[10]) + 0.188821573960182 * (y[3] + y[9]) + 0.199773405226859 * (y[4] + y[8]) + 0.224926465333340 * (y[5] + y[7]) + 0.242611071901408 * y[6]);
                erri1 = Math.Abs(i1 -iS);
                erri2 = Math.Abs(i2 -iS);
                r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;
                toler = (r > 0.0 && r < 1.0) ? RelativeTolerance / r : RelativeTolerance;
                if (iS == 0.0) iS = b - a;
                iS= Math.Abs(iS);

                Result = adaptlob(Function, a, b, fa, fb, iS);

            }
            catch (Exception ex)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw ex;
            }

        }

        /// <summary>
        /// Helper function for recursion.
        /// </summary>
        private double adaptlob(Func<double, double> function, double a, double b, double fa, double fb, double iS)
        {
            double m, h, mll, ml, mr, mrr, fmll, fml, fm, fmrr, fmr, i1, i2;
            m = 0.5 * (a + b);
            h = 0.5 * (b - a);
            mll = m - alpha * h;
            ml = m - beta * h;
            mr = m + beta * h;
            mrr = m + alpha * h;
            fmll = function(mll);
            fml = function(ml);
            fm = function(m);
            fmr = function(mr);
            fmrr = function(mrr);

            Iterations++;
            FunctionEvaluations += 5;

            // 4-point Gauss-Lobatto formula
            i2 = h / 6.0 * (fa + fb + 5.0 * (fml + fmr)); 
            // 7-point Kronrod extension
            i1 = h / 1470.0 * (77.0 * (fa + fb) + 432.0 * (fmll + fmrr) + 625.0 * (fml + fmr) + 672.0 * fm);

            // Check iterations
            if (Iterations >= MaxIterations)
            {
                // Terminate recursion
                UpdateStatus(IntegrationStatus.MaximumIterationsReached);
                return i1;
            }
            // Check function evaluations
            if (FunctionEvaluations >= MaxFunctionEvaluations)
            {
                // Terminate recursion
                UpdateStatus(IntegrationStatus.MaximumFunctionEvaluationsReached);           
                return i1;
            }
            // Check tolerance
            if (Math.Abs(i1 - i2) <= toler * iS || mll <= a || b <= mrr)
            {
                if ((mll <= a || b <= mrr) && terminate)
                {
                    // Interval contains no more machine numbers
                    outOfTolerance = true;
                    terminate = false;
                }
                // Terminate recursion
                UpdateStatus(IntegrationStatus.Success);
                return i1;
            }
            else
                // Subdivide interval
                return adaptlob(function, a, mll, fa, fmll,iS) +
                       adaptlob(function, mll, ml, fmll, fml,iS) +
                       adaptlob(function, ml, m, fml, fm,iS) +
                       adaptlob(function, m, mr, fm, fmr,iS) +
                       adaptlob(function, mr, mrr, fmr, fmrr,iS) +
                       adaptlob(function, mrr, b, fmrr, fb,iS);
        }

    }
}
