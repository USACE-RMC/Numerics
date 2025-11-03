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

        private double a, b, _squaredError;

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
                // Fist evaluation of Simpson's Rule on the whole interval
                double m = (a + b) / 2d;
                double fa = Function(a);
                double fb = Function(b);
                double fm = Function(m);
                FunctionEvaluations +=3; // Count the three evaluations: fa, fb, fm
                double whole = Math.Abs(b - a) / 6d * (fa + 4d * fm + fb);

                // Recursively sub-divide
                Result = AdaptiveSimpsons(Function, a, fa, b, fb, MaxDepth, whole, m, fm, a, b);

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
            catch (Exception)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw;
            }

        }

        /// <summary>
        /// Evaluates the integral.
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
                double sigmaSquared = 0;
                for (int i = 0; i < bins.Count; i++)
                {
                    // Fist evaluation of Simpson's Rule on the whole interval
                    double a = bins[i].LowerBound;
                    double b = bins[i].UpperBound;
                    double m = (a + b) / 2d;
                    double fa = Function(a);
                    double fb = Function(b);
                    double fm = Function(m);
                    FunctionEvaluations += 3; // Count the three evaluations: fa, fb, fm
                    double whole = Math.Abs(b - a) / 6d * (fa + 4d * fm + fb);

                    // Recursively sub-divide
                    mu += AdaptiveSimpsons(Function, a, fa, b, fb, MaxDepth, whole, m, fm, a, b);

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
            catch (Exception)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw;
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
        /// <param name="epsilon"> Machine epsilon </param>
        /// <param name="depth"> Less than or equal to 0 (max recursions have been reached) </param>
        /// <param name="whole"> The original whole three point Simpson's Rule evaluation on [a,b] </param>
        /// <param name="m"> The midpoint between a and b </param>
        /// <param name="fm"> The function evaluated at m </param>
        /// <param name="a0"> The original lower bound of the integral </param>
        /// <param name="b0"> The original upper bound of the integral </param>
        /// <returns>
        /// An evaluation of Simpson's Rule with the error less than a certain tolerance. This is accomplished by subdividing the interval 
        /// until the error between the last evaluation and the current evaluation is sufficiently small.
        /// </returns>
        private double AdaptiveSimpsons(Func<double, double> f, double a, double fa, double b, double fb, int depth, double whole, double m, double fm, double a0, double b0)
        {
            double h = (b - a) * 0.5;
            double lm = a + h * 0.5;  // left mid
            double rm = a + h + h * 0.5;  // right mid
            double flm = f(lm), frm = f(rm);
            FunctionEvaluations += 2;

            double left = h / 6 * (fa + 4 * flm + fm);
            double right = h / 6 * (fm + 4 * frm + fb);

            // Calculate error for current interval
            double error = (left + right - whole);
            double delta = error / 15d; // Richardson correction

            // Richardson-based convergence tolerance
            double toleranceScaled = 15d * (RelativeTolerance * Math.Abs(b - a) / Math.Abs(b0 - a0));

            // Absolute and Relative tolerance checks
            bool absoluteToleranceReached = Math.Abs(error) <= AbsoluteTolerance;
            bool relativeToleranceReached = Math.Abs(error) <= toleranceScaled;

            // Check if convergence criteria are met
            if (depth <= 0 || Math.Abs(a - b) <= Tools.DoubleMachineEpsilon || FunctionEvaluations >= MaxFunctionEvaluations ||
                (FunctionEvaluations >= MinFunctionEvaluations && depth <= MaxDepth - MinDepth && (absoluteToleranceReached || relativeToleranceReached)))
            {
                // Convergence is reached
                _squaredError += delta * delta; // Accumulate squared errors
                return left + right + delta;
            }
            else
            {
                // Recursively subdivide the intervals and accumulate results
                var leftResult = AdaptiveSimpsons(f, a, fa, m, fm, depth - 1, left, lm, flm, a0, b0);
                var rightResult = AdaptiveSimpsons(f, m, fm, b, fb, depth - 1, right, rm, frm, a0, b0);
                return leftResult + rightResult;
            }
        }

    }
}
