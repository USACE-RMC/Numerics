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
    /// A class for Trapezoidal rule integration. Integration steps are refined until convergence.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// In its most basic form, the Trapezoidal rule approximates the region under the graph of f(x) as a trapezoid and calculating its area as:
    /// <code>
    ///         b
    ///         ∫ f(x) dx ~ T(a, b) = (b - a) * [ f(a) + f(b) ] / 2
    ///         a
    /// </code>
    /// This method can also be seen as averaging the left and right Riemann sums. Moreover, the interval can be subdivided to improve the accuracy of the approximation, with
    /// the trapezoidal rule applied to every sub-interval. Here the integration steps are refined, i.e. the interval is subdivided, until convergence.
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Trapezoidal_rule"/>
    /// </remarks>
    [Serializable]
    public class TrapezoidalRule : Integrator
    {
        /// <summary>
        /// Construct a new trapezoidal rule class. 
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="min">The minimum value under which the integral must be computed.</param>
        /// <param name="max">The maximum value under which the integral must be computed.</param>
        public TrapezoidalRule(Func<double, double> function, double min, double max)
        {
            if (function == null) throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            if (max <= min) throw new ArgumentNullException(nameof(max), "The maximum value cannot be less than or equal to the minimum value.");

            Function = function;
            a = min;
            b = max;
        }

        private double a;
        private double b;
        private double _s = 0;

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

        /// <inheritdoc/>
        public override void Integrate()
        {
            ClearResults();
            Validate();

            try
            {
                // Integrate
                double s0 = 0, s1 = 0;
                for (int i = 0; i < MaxIterations; i++)
                {
                    s0 = Next();

                    // Check function evaluations
                    if (FunctionEvaluations >= MaxFunctionEvaluations)
                    {
                        UpdateStatus(IntegrationStatus.MaximumFunctionEvaluationsReached);
                        return;
                    }
                       
                    // Check convergence
                    if (i > 2)
                    {
                        if (EvaluateConvergence(s1, s0) || (s0 == 0.0 && s1 == 0.0))
                        {
                            Result = s0;
                            UpdateStatus(IntegrationStatus.Success);
                            return;
                        }
                    }
                    s1 = s0;
                }
                // If we get to here, then the maximum number of steps were reached before converging. 
                UpdateStatus(IntegrationStatus.MaximumIterationsReached);
            }
            catch (Exception ex)
            {
                UpdateStatus(IntegrationStatus.Failure, ex);
            }       
        }

        /// <summary>
        /// Returns the value of the integral at the nth step of refinement.
        /// </summary>
        private double Next()
        {
            double x, tnm, sum, del;
            int it, j;
            Iterations++;
            if (Iterations == 1)
            {
                _s = 0.5 * (b - a) * (Function(a) + Function(b));
                FunctionEvaluations += 2;
            }
            else
            {
                for (it = 1, j = 1; j < Iterations - 1; j++) it <<= 1;
                tnm = it;
                del = (b - a) / tnm;
                x = a + 0.5 * del;
                for (sum = 0.0, j = 0; j < it; j++, x += del)
                {
                    sum += Function(x);
                    FunctionEvaluations += 1;
                }
                _s = 0.5 * (_s + (b - a) * sum / tnm);
            }
            return _s;
        }

    }
}
