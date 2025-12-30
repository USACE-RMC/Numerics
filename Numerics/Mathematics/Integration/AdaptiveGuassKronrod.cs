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
    /// A class that performs adaptive Gauss-Kronrod integration. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Adaptive Gauss-Kronrod quadrature uses a pair of nested quadrature rules: an n-point Gauss rule 
    /// and a (2n+1)-point Kronrod extension. The Kronrod rule reuses all the Gauss points and adds 
    /// n+1 optimally-placed points to provide both a higher-order estimate and an accurate error estimate. 
    /// This implementation uses the 10-point Gauss, 21-point Kronrod rule (G10K21), which is 21st order 
    /// accurate for smooth functions. The difference between the two estimates provides a reliable error 
    /// indicator, and intervals are adaptively subdivided until the error tolerance is met.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><see href="https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula"/></item>
    /// <item>Piessens, R., et al. (1983). QUADPACK: A Subroutine Package for Automatic Integration. Springer-Verlag.</item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class AdaptiveGaussKronrod : Integrator
    {
        /// <summary>
        /// Constructs a new adaptive Gauss-Kronrod rule.
        /// </summary>
        /// <param name="function">The function to integrate.</param>
        /// <param name="min">The minimum value under which the integral must be computed.</param>
        /// <param name="max">The maximum value under which the integral must be computed.</param>
        public AdaptiveGaussKronrod(Func<double, double> function, double min, double max)
        {
            if (max <= min) throw new ArgumentNullException(nameof(max), "The maximum value cannot be less than or equal to the minimum value.");
            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            a = min;
            b = max;
        }

        private double a, b, _squaredError;

        // Gauss-Kronrod G10K21 nodes and weights (for interval [-1, 1])
        // 10-point Gauss nodes (symmetric, so we only store positive values)
        private static readonly double[] xGauss = new double[]
        {
            0.973906528517171720077964012084452,
            0.865063366688984510732096688423493,
            0.679409568299024406234327365114874,
            0.433395394129247190799265943165784,
            0.148874338981631210884826001129720
        };

        // 10-point Gauss weights
        private static readonly double[] wGauss = new double[]
        {
            0.066671344308688137593568809893332,
            0.149451349150580593145776339657697,
            0.219086362515982043995534934228163,
            0.269266719309996355091226921569469,
            0.295524224714752870173892994651338
        };

        // 21-point Kronrod nodes (symmetric, so we only store positive values)
        // These include all 10 Gauss nodes plus 11 additional nodes
        private static readonly double[] xKronrod = new double[]
        {
            0.995657163025808080735527280689003,
            0.973906528517171720077964012084452,
            0.930157491355708226001207180059508,
            0.865063366688984510732096688423493,
            0.780817726586416897063717578345042,
            0.679409568299024406234327365114874,
            0.562757134668604683339000099272694,
            0.433395394129247190799265943165784,
            0.294392862701460198131126603103866,
            0.148874338981631210884826001129720,
            0.000000000000000000000000000000000
        };

        // 21-point Kronrod weights
        private static readonly double[] wKronrod = new double[]
        {
            0.011694638867371874278064396062192,
            0.032558162307964727478818972459390,
            0.054755896574351996031381300244580,
            0.075039674810919952767043140916190,
            0.093125454583697605535065465083366,
            0.109387158802297641899210590325805,
            0.123491976262065851077958109831074,
            0.134709217311473325928054001771707,
            0.142775938577060080797094273138717,
            0.147739104901338491374841515972068,
            0.149445554002916905664936468389821
        };

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
                // Initial evaluation using Gauss-Kronrod rule on the whole interval
                var (kronrodResult, gaussResult) = EvaluateGaussKronrod(a, b);

                // Recursively sub-divide
                Result = AdaptiveGK(Function, a, b, MaxDepth, kronrodResult, gaussResult, a, b);

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
                    // Initial evaluation using Gauss-Kronrod rule on the bin interval
                    double binA = bins[i].LowerBound;
                    double binB = bins[i].UpperBound;
                    var (kronrodResult, gaussResult) = EvaluateGaussKronrod(binA, binB);

                    // Recursively sub-divide
                    mu += AdaptiveGK(Function, binA, binB, MaxDepth, kronrodResult, gaussResult, binA, binB);
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
        /// Evaluates the Gauss-Kronrod G10K21 rule over the interval [a, b].
        /// </summary>
        /// <param name="a">The lower bound of integration.</param>
        /// <param name="b">The upper bound of integration.</param>
        /// <returns>A tuple containing (Kronrod estimate, Gauss estimate).</returns>
        private (double kronrod, double gauss) EvaluateGaussKronrod(double a, double b)
        {
            double center = 0.5 * (a + b);
            double halfLength = 0.5 * (b - a);

            double resultGauss = 0.0;
            double resultKronrod = 0.0;

            // Evaluate at center point (x = 0)
            double f0 = Function(center);
            resultKronrod += wKronrod[10] * f0;
            FunctionEvaluations++;

            // Evaluate at symmetric pairs of points
            for (int i = 0; i < 10; i++)
            {
                double abscissa = halfLength * xKronrod[i];
                double fval1 = Function(center - abscissa);
                double fval2 = Function(center + abscissa);
                double fsum = fval1 + fval2;

                resultKronrod += wKronrod[i] * fsum;

                // Add to Gauss result if this is a Gauss node
                // Gauss nodes are at indices 1, 3, 5, 7, 9 of the Kronrod array
                if (i % 2 == 1)
                {
                    int gaussIndex = i / 2;
                    resultGauss += wGauss[gaussIndex] * fsum;
                }

                FunctionEvaluations += 2;
            }

            resultGauss *= halfLength;
            resultKronrod *= halfLength;

            return (resultKronrod, resultGauss);
        }

        /// <summary>
        /// A helper function for adaptive Gauss-Kronrod integration.
        /// </summary>
        /// <param name="f">The unidimensional function to integrate.</param>
        /// <param name="a">The minimum value under which the integral must be computed.</param>
        /// <param name="b">The maximum value under which the integral must be computed.</param>
        /// <param name="depth">The current recursion depth remaining.</param>
        /// <param name="kronrodWhole">The Kronrod evaluation on [a,b].</param>
        /// <param name="gaussWhole">The Gauss evaluation on [a,b].</param>
        /// <param name="a0">The original lower bound of the integral.</param>
        /// <param name="b0">The original upper bound of the integral.</param>
        /// <returns>
        /// An evaluation of the integral using adaptive Gauss-Kronrod with error less than the specified tolerance. 
        /// This is accomplished by subdividing the interval until the error between the Gauss and Kronrod estimates 
        /// is sufficiently small.
        /// </returns>
        private double AdaptiveGK(Func<double, double> f, double a, double b, int depth,
            double kronrodWhole, double gaussWhole, double a0, double b0)
        {
            // Error estimate: difference between Kronrod and Gauss results
            double error = Math.Abs(kronrodWhole - gaussWhole);

            // Scaled tolerance based on interval length relative to original domain
            double toleranceScaled = RelativeTolerance * Math.Abs(b - a) / Math.Abs(b0 - a0);

            // Absolute and Relative tolerance checks
            bool absoluteToleranceReached = error <= AbsoluteTolerance;
            bool relativeToleranceReached = error <= toleranceScaled * Math.Abs(kronrodWhole);

            // Check if convergence criteria are met
            if (depth <= 0 || Math.Abs(a - b) <= Tools.DoubleMachineEpsilon || FunctionEvaluations >= MaxFunctionEvaluations ||
                (FunctionEvaluations >= MinFunctionEvaluations && depth <= MaxDepth - MinDepth && (absoluteToleranceReached || relativeToleranceReached)))
            {
                // Convergence is reached
                _squaredError += error * error; // Accumulate squared errors
                return kronrodWhole; // Return the more accurate Kronrod estimate
            }
            else
            {
                // Subdivide the interval at the midpoint
                double m = (a + b) / 2.0;

                // Evaluate Gauss-Kronrod on left half
                var (kronrodLeft, gaussLeft) = EvaluateGaussKronrod(a, m);

                // Evaluate Gauss-Kronrod on right half
                var (kronrodRight, gaussRight) = EvaluateGaussKronrod(m, b);

                // Recursively subdivide the intervals and accumulate results
                var leftResult = AdaptiveGK(f, a, m, depth - 1, kronrodLeft, gaussLeft, a0, b0);
                var rightResult = AdaptiveGK(f, m, b, depth - 1, kronrodRight, gaussRight, a0, b0);

                return leftResult + rightResult;
            }
        }
    }
}
