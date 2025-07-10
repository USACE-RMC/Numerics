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

using Numerics.Sampling;
using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Integration
{

    /// <summary>
    /// A class that performs adaptive Simpson's integration. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Adaptive Simpson's rule uses an estimate of the error that comes from calculating a definite integral with Simpson's rule. If the error between the previous
    /// evaluation of the rule and the current evaluation of the rule exceeds a certain specified tolerance, the rule calls for subdividing the interval. Adaptive Simpson's
    /// rule is applied to each subinterval in a recursive manner until the error qualifications are met.
    /// <code>
    ///             | S(a, b) - S(a, m) + S(m, b) | &lt; epsilon
    /// </code>
    /// where a and b are the bounds of integration and m is the midpoint between them.
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Adaptive_Simpson%27s_method"/>
    /// </remarks>
    [Serializable]
    public class AdaptiveSimpsonsRule : Integrator
    {

        /// <summary>
        /// Constructs a new adaptive Simpson's rule.
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="min">The minimum value under which the integral must be computed.</param>
        /// <param name="max">The maximum value under which the integral must be computed.</param>
        public AdaptiveSimpsonsRule(Func<double, double> function, double min, double max)
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
            StandardError = 0;
            ClearResults();
            Validate();

            try
            {
                double fa = Function(a);
                double fb = Function(b);
                FunctionEvaluations += 2;
                double m = 0, fm = 0, whole = Simpsons(Function, a, fa, b, fb, ref m, ref fm);       
                Result = AdaptiveSimpsons(Function, a, fa, b, fb, RelativeTolerance, MaxDepth, whole, m, fm);

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
        /// Evaluates the integral.
        /// </summary>
        /// <param name="bins">The stratification bins to integrate over.</param>
        public void Integrate(List<StratificationBin> bins)
        {
            StandardError = 0;
            ClearResults();
            Validate();

            try
            {
                double mu = 0;
                double sigma = 0;
                for (int i = 0; i < bins.Count; i++)
                {
                    double a = bins[i].LowerBound;
                    double b = bins[i].UpperBound;
                    double fa = Function(a);
                    double fb = Function(b);
                    FunctionEvaluations += 2;
                    double m = 0, fm = 0, whole = Simpsons(Function, a, fa, b, fb, ref m, ref fm);
                    mu += AdaptiveSimpsons(Function, a, fa, b, fb, RelativeTolerance, MaxDepth, whole, m, fm);
                    sigma += StandardError;
                }

                Result = mu;
                StandardError = sigma;
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
        /// A helper function to the Integrate() function
        /// </summary>
        /// <param name="f"> The unidimensional function to integrate </param>
        /// <param name="a"> The minimum value under which the integral must be computed </param>
        /// <param name="fa"> The function evaluated at a </param>
        /// <param name="b"> The maximum value under which the integral must be computed </param>
        /// <param name="fb"> The function evaluated at b </param>
        /// <param name="m"> The midpoint between a and b </param>
        /// <param name="fm"> The function evaluated at m </param>
        /// <returns>
        /// A three point Simpson's Rule evaluation on [a,b]
        /// </returns>
        private double Simpsons(Func<double, double> f, double a, double fa, double b, double fb, ref double m, ref double fm)
        {
            m = (a + b) / 2d;
            fm = f(m);
            FunctionEvaluations++;           
            return Math.Abs(b - a) / 6d * (fa + 4d * fm + fb);
        }

        /// <summary>
        /// A helper function to the Integrate() function
        /// </summary>
        /// <param name="f"> The unidimensional function to integrate </param>
        /// <param name="a"> The minimum value under which the integral must be computed </param>
        /// <param name="fa"> The function evaluated at a </param>
        /// <param name="b"> The maximum value under which the integral must be computed </param>
        /// <param name="fb"> The function evaluated at b </param>
        /// <param name="epsilon"> Machine epsilon </param>
        /// <param name="depth"> Less than or equal to 0 (max recursions have been reached) </param>
        /// <param name="whole"> The original whole three point Simpson's Rule evaluation on [a,b] </param>
        /// <param name="m"> The midpoint between a and b </param>
        /// <param name="fm"> The function evaluated at m </param>
        /// <returns>
        /// An evaluation of Simpson's Rule with the error less than a certain tolerance. This is accomplished but subdividing the interval the rule is
        /// evaluated on until the error between the last evaluation and the current evaluation is sufficiently small.
        /// </returns>
        private double AdaptiveSimpsons(Func<double, double> f, double a, double fa, double b, double fb, double epsilon, int depth, double whole, double m, double fm)
        {

            double lm = 0, flm = 0, left = Simpsons(f, a, fa, m, fm, ref lm, ref flm);
            double rm = 0, frm = 0, right = Simpsons(f, m, fm, b, fb, ref rm, ref frm);
            double delta = (left + right - whole) / 15d;

            // Check tolerance 
            // - Depth is less than or equal to 0 (max recursions have been reached)
            // - Abs(a-b) is smaller than machine epsilon
            // - Maximum number of function evaluations
            // - Minimum number of function evaluations, the minimum depth, and also delta less than tolerance
            if (depth <= 0 ||  Math.Abs(a - b) <= Tools.DoubleMachineEpsilon || FunctionEvaluations >= MaxFunctionEvaluations ||  
                (FunctionEvaluations >= MinFunctionEvaluations && depth <= MaxDepth - MinDepth && Math.Abs(delta) <= epsilon + epsilon * Math.Abs(whole)))
            {
                // convergence is reached
                // Terminate recursion
                StandardError += Math.Abs(delta);
                return left + right + delta; 
            }
            else
            {
                // Subdivide interval
                var l = AdaptiveSimpsons(f, a, fa, m, fm, epsilon, depth - 1, left, lm, flm);
                var r = AdaptiveSimpsons(f, m, fm, b, fb, epsilon, depth - 1, right, rm, frm);
                return l + r;
            }           
        }


    }
}
