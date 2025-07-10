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
    /// A class for Simpson's rule integration. Integration steps are refined until convergence.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Simpson's rule approximates the integral by splitting the interval and evaluating the function at those split points. The most basic implementation of this is the
    /// following equation for the integral of f(x) over [a, b] with 2 subdivisions:
    /// <code>
    ///         b
    ///         ∫ f(x) dx ~ S(a, b) = (b - a) / 6 * [ f(a) + 4*f(m) + f(b) ]
    ///         a
    /// </code>
    /// where m is the midpoint [ (a+b) / 2 ] between a and b. The more subdivisions there are, i.e. the more refined the steps are, the more accurate the approximation becomes.
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Simpson%27s_rule"/>
    /// </remarks>
    [Serializable]
    public class SimpsonsRule : Integrator
    {
        /// <summary>
        /// Construct a new Simpson's rule class. 
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="min">Start point for integration.</param>
        /// <param name="max">End point for integration.</param>
        public SimpsonsRule(Func<double, double> function, double min, double max)
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
                double s = 0, os = 0;
                double st = 0, ost = 0;
                for (int i = 0; i < MaxIterations; i++)
                {
                    st = Next();
                    s = (4d * st - ost) / 3d;

                    // Check function evaluations
                    if (FunctionEvaluations >= MaxFunctionEvaluations)
                    {
                        UpdateStatus(IntegrationStatus.MaximumFunctionEvaluationsReached);
                        return;
                    }

                    // Check convergence
                    if (i > 2)
                    {
                        if (EvaluateConvergence(os, s) || (s == 0.0 && os == 0.0))
                        {
                            Result = s;
                            UpdateStatus(IntegrationStatus.Success);
                            return;
                        }
                    }
                    os = s;
                    ost = st;
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
