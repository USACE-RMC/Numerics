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
using Numerics.Sampling;

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
    /// <para>
    /// <b> Description: </b>
    /// Adaptive Gauss-Lobatto-Kronrod quadrature uses nested quadrature rules of 4, 7, and 13 points.
    /// The Gauss-Lobatto rules include the endpoints, making them particularly suitable for functions
    /// with endpoint singularities. The method compares estimates from different order rules to 
    /// estimate error and adaptively subdivides intervals until the error tolerance is met.
    /// The 13-point Kronrod extension provides high accuracy while the comparison between 4-point and 
    /// 7-point estimates helps determine convergence reliability.
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

        private double a, b, _squaredError;

        // Abscissas for Gauss-Lobatto-Kronrod quadrature
        private static readonly double alpha = Math.Sqrt(2.0 / 3.0);
        private static readonly double beta = 1.0 / Math.Sqrt(5.0);
        private static readonly double x1 = 0.942882415695480;
        private static readonly double x2 = 0.641853342345781;
        private static readonly double x3 = 0.236383199662150;
        private static readonly double[] x = new double[] { 0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1 };

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

        /// <summary>
        /// The minimum recursion depth. Default = 0.
        /// </summary>
        public int MinDepth { get; set; } = 0;

        /// <summary>
        /// The maximum recursion depth. Default = 100.
        /// </summary>
        public int MaxDepth { get; set; } = 100;

        /// <summary>
        /// Returns an approximate measure of the standard error of the integration. 
        /// </summary>
        public double StandardError { get; private set; }

        /// <inheritdoc/>
        public override void Integrate()
        {
            _squaredError = 0;
            StandardError = 0;
            ClearResults();
            Validate();

            try
            {
                double m = 0.5 * (a + b);
                double h = 0.5 * (b - a);
                var y = new double[13];

                // Evaluate at all 13 points
                y[0] = Function(a);
                y[12] = Function(b);
                FunctionEvaluations += 2;

                for (int i = 1; i < 12; i++)
                {
                    y[i] = Function(m + x[i] * h);
                    FunctionEvaluations++;
                }

                // 4-point Gauss-Lobatto formula
                double i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));

                // 7-point Kronrod extension
                double i1 = (h / 1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]);

                // 13-point Kronrod extension
                double iS = h * (0.0158271919734802 * (y[0] + y[12]) +
                                 0.0942738402188500 * (y[1] + y[11]) +
                                 0.155071987336585 * (y[2] + y[10]) +
                                 0.188821573960182 * (y[3] + y[9]) +
                                 0.199773405226859 * (y[4] + y[8]) +
                                 0.224926465333340 * (y[5] + y[7]) +
                                 0.242611071901408 * y[6]);

                double erri1 = Math.Abs(i1 - iS);
                double erri2 = Math.Abs(i2 - iS);
                double r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;
                double toler = (r > 0.0 && r < 1.0) ? RelativeTolerance / r : RelativeTolerance;

                double iSabs = (iS == 0.0) ? (b - a) : Math.Abs(iS);

                // Recursively sub-divide
                Result = AdaptiveLobatto(Function, a, b, y[0], y[12], MaxDepth, iSabs, toler);

                // Standard error calculated after recursion completes
                StandardError = Math.Sqrt(_squaredError);

                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else
                {
                    Status = IntegrationStatus.Success;
                }
            }
            catch (Exception ex)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw ex;
            }
        }

        /// <summary>
        /// Evaluates the integral over user-provided stratification bins.
        /// </summary>
        /// <param name="bins">The stratification bins to integrate over.</param>
        public void Integrate(List<StratificationBin> bins)
        {
            _squaredError = 0;
            StandardError = 0;
            ClearResults();
            Validate();

            try
            {
                double mu = 0;
                for (int i = 0; i < bins.Count; i++)
                {
                    double binA = bins[i].LowerBound;
                    double binB = bins[i].UpperBound;
                    double m = 0.5 * (binA + binB);
                    double h = 0.5 * (binB - binA);
                    var y = new double[13];

                    // Evaluate at all 13 points
                    y[0] = Function(binA);
                    y[12] = Function(binB);
                    FunctionEvaluations += 2;

                    for (int j = 1; j < 12; j++)
                    {
                        y[j] = Function(m + x[j] * h);
                        FunctionEvaluations++;
                    }

                    // 4-point Gauss-Lobatto formula
                    double i2 = (h / 6.0) * (y[0] + y[12] + 5.0 * (y[4] + y[8]));

                    // 7-point Kronrod extension
                    double i1 = (h / 1470.0) * (77.0 * (y[0] + y[12]) + 432.0 * (y[2] + y[10]) + 625.0 * (y[4] + y[8]) + 672.0 * y[6]);

                    // 13-point Kronrod extension
                    double iS = h * (0.0158271919734802 * (y[0] + y[12]) +
                                     0.0942738402188500 * (y[1] + y[11]) +
                                     0.155071987336585 * (y[2] + y[10]) +
                                     0.188821573960182 * (y[3] + y[9]) +
                                     0.199773405226859 * (y[4] + y[8]) +
                                     0.224926465333340 * (y[5] + y[7]) +
                                     0.242611071901408 * y[6]);

                    double erri1 = Math.Abs(i1 - iS);
                    double erri2 = Math.Abs(i2 - iS);
                    double r = (erri2 != 0.0) ? erri1 / erri2 : 1.0;
                    double toler = (r > 0.0 && r < 1.0) ? RelativeTolerance / r : RelativeTolerance;

                    double iSabs = (iS == 0.0) ? (binB - binA) : Math.Abs(iS);

                    // Recursively sub-divide
                    mu += AdaptiveLobatto(Function, binA, binB, y[0], y[12], MaxDepth, iSabs, toler);
                }

                // Final result and standard error
                Result = mu;
                StandardError = Math.Sqrt(_squaredError);

                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else
                {
                    Status = IntegrationStatus.Success;
                }
            }
            catch (Exception ex)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw ex;
            }
        }

        /// <summary>
        /// A helper function for adaptive Gauss-Lobatto integration.
        /// Recursively subdivides the integration domain and applies nested Gauss-Lobatto-Kronrod rules.
        /// </summary>
        /// <param name="f">The unidimensional function to integrate.</param>
        /// <param name="a">The minimum value under which the integral must be computed.</param>
        /// <param name="b">The maximum value under which the integral must be computed.</param>
        /// <param name="fa">The function evaluated at a.</param>
        /// <param name="fb">The function evaluated at b.</param>
        /// <param name="depth">The current recursion depth remaining.</param>
        /// <param name="iSabs">The absolute value of the integral estimate (for scaling tolerance).</param>
        /// <param name="toler">The adaptive tolerance for this interval.</param>
        /// <returns>
        /// An evaluation of the integral using adaptive Gauss-Lobatto-Kronrod with error less than the specified tolerance. 
        /// This is accomplished by subdividing the interval into 6 subintervals until the error between different order 
        /// estimates is sufficiently small.
        /// </returns>
        private double AdaptiveLobatto(Func<double, double> f, double a, double b, double fa, double fb,
            int depth, double iSabs, double toler)
        {
            double m = 0.5 * (a + b);
            double h = 0.5 * (b - a);
            double mll = m - alpha * h;
            double ml = m - beta * h;
            double mr = m + beta * h;
            double mrr = m + alpha * h;

            double fmll = f(mll);
            double fml = f(ml);
            double fm = f(m);
            double fmr = f(mr);
            double fmrr = f(mrr);
            FunctionEvaluations += 5;

            // 4-point Gauss-Lobatto formula on current interval
            double i2 = h / 6.0 * (fa + fb + 5.0 * (fml + fmr));

            // 7-point Kronrod extension on current interval
            double i1 = h / 1470.0 * (77.0 * (fa + fb) + 432.0 * (fmll + fmrr) + 625.0 * (fml + fmr) + 672.0 * fm);

            // Error estimate: difference between 4-point and 7-point estimates
            double error = Math.Abs(i1 - i2);

            // Check if convergence criteria are met
            // The original algorithm uses toler * iSabs as the tolerance threshold
            if (depth <= 0 || Math.Abs(a - b) <= Tools.DoubleMachineEpsilon ||
                mll <= a || b <= mrr || FunctionEvaluations >= MaxFunctionEvaluations ||
                error <= toler * iSabs)
            {
                // Convergence is reached
                _squaredError += error * error; // Accumulate squared errors
                return i1; // Return the more accurate 7-point estimate
            }
            else
            {
                // Subdivide interval into 6 subintervals
                // Pass the same iSabs and toler through all recursive calls
                var r1 = AdaptiveLobatto(f, a, mll, fa, fmll, depth - 1, iSabs, toler);
                var r2 = AdaptiveLobatto(f, mll, ml, fmll, fml, depth - 1, iSabs, toler);
                var r3 = AdaptiveLobatto(f, ml, m, fml, fm, depth - 1, iSabs, toler);
                var r4 = AdaptiveLobatto(f, m, mr, fm, fmr, depth - 1, iSabs, toler);
                var r5 = AdaptiveLobatto(f, mr, mrr, fmr, fmrr, depth - 1, iSabs, toler);
                var r6 = AdaptiveLobatto(f, mrr, b, fmrr, fb, depth - 1, iSabs, toler);

                return r1 + r2 + r3 + r4 + r5 + r6;
            }
        }
    }
}
