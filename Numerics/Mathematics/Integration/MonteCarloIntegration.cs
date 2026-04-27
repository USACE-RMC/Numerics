using Numerics.Sampling;
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
            Random = new MersenneTwister();
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
