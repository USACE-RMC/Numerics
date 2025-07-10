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
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class for Monte Carlo integration for multidimensional integration.
    /// </summary>
    /// <remarks>
    /// <para>
    ///    <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This method numerically computes a definite integral by randomly choosing points at which the integrand is evaluated (rather than evaluating
    /// the integral at regular grid like other methods).
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Monte_Carlo_integration"/>
    /// </remarks>
    [Serializable]
    public class MonteCarloIntegration : Integrator
    {
        /// <summary>
        /// Creates a new Monte Carlo Integration class.
        /// </summary>
        /// <param name="function">The multidimensional function to integrate.</param>
        /// <param name="dimensions">The number of dimensions in the function to evaluate.</param>
        /// <param name="min">The minimum values under which the integral must be computed.</param>
        /// <param name="max">The maximum values under which the integral must be computed.</param>
        public MonteCarloIntegration(Func<double[], double> function, int dimensions, IList<double> min, IList<double> max)
        {
            if (dimensions < 1 ) throw new ArgumentOutOfRangeException(nameof(dimensions), "There must be at least 1 dimension to evaluate.");

            // Check if the length of the min and max values equal the number of dimensions
            if (min.Count != dimensions || max.Count != dimensions)
            {
                throw new ArgumentOutOfRangeException(nameof(min), "The minimum and maximum values must be the same length as the number of dimensions.");
            }

            // Check if the minimum values are less than the maximum values
            for (int i = 0; i < min.Count; i++)
            {
                if (max[i] <= min[i])
                {
                    throw new ArgumentOutOfRangeException(nameof(max), "The maximum values cannot be less than or equal to the minimum values.");
                }
            }

            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            Dimensions = dimensions;
            Min = min.ToArray();
            Max = max.ToArray();
            Random = new Random();
        }

        /// <summary>
        /// The multidimensional function to integrate.
        /// </summary>
        public Func<double[], double> Function { get; }

        /// <summary>
        /// The number of dimensions in the function to evaluate./>.
        /// </summary>
        public int Dimensions { get; }

        /// <summary>
        /// The minimum values under which the integral must be computed.
        /// </summary>
        public double[] Min { get; }

        /// <summary>
        /// The maximum values under which the integral must be computed. 
        /// </summary>
        public double[] Max { get; }

        /// <summary>
        /// The random number generator to be used within the Monte Carlo integration.
        /// </summary>
        public Random Random { get; set; }

        /// <summary>
        /// The integration error. 
        /// </summary>
        public double StandardError { get; protected set; }

        /// <summary>
        /// Determines whether to use a Sobol sequence or a pseudo-Random number generator. 
        /// </summary>
        public bool UseSobolSequence { get; set; } = true;

        /// <inheritdoc/>
        public override void Integrate()
        {
            ClearResults();
            Validate();

            try
            {
                var sample = new double[Dimensions];

                double sum = 0, sum2 = 0;
                double avg = 0, avg2 = 0;

                // Get dx
                double dx = 1;
                for (int i = 0; i < Dimensions; i++)
                    dx *= (Max[i] - Min[i]);


                for (int i = 1; i <= MaxIterations; i++)
                {
                    for (int j = 0; j < sample.Length; j++)
                        sample[j] = Min[j] + Random.NextDouble() * (Max[j] - Min[j]);

                    double f = Function(sample);

                    Iterations++;
                    FunctionEvaluations++;
                    sum += f;
                    sum2 += f * f;

                    avg = sum / Iterations;
                    avg2 = sum2 / Iterations;

                    Result = avg * dx;
                    StandardError = Math.Sqrt((avg2 - avg * avg) / Iterations) * dx;
                    
                    // Check tolerance
                    if (Iterations > MinIterations && (Math.Abs(StandardError / Result) < RelativeTolerance))
                    {
                        break;
                    }
                }

                // Update status
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
                UpdateStatus(IntegrationStatus.Failure, ex);
            }
        
        }
    }
}
