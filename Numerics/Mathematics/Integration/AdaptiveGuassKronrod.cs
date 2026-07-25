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
            if (max <= min) throw new ArgumentOutOfRangeException(nameof(max), "The maximum value cannot be less than or equal to the minimum value.");
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

        // The unscaled Kronrod weight per capture slot: slot 0 is the center node (wKronrod[10]),
        // then the symmetric pair for each abscissa index i occupies slots 2i+1 and 2i+2.
        private static readonly double[] wCapture = BuildCaptureWeights();

        /// <summary>
        /// Builds the unscaled Kronrod weight layout matching the node-capture order.
        /// </summary>
        /// <returns>The 21 unscaled weights in capture order.</returns>
        private static double[] BuildCaptureWeights()
        {
            var weights = new double[21];
            weights[0] = wKronrod[10];
            for (int i = 0; i < 10; i++)
            {
                weights[2 * i + 1] = wKronrod[i];
                weights[2 * i + 2] = wKronrod[i];
            }
            return weights;
        }

        /// <summary>
        /// The unidimensional function to integrate.
        /// </summary>
        public Func<double, double> Function { get; }

        /// <summary>
        /// An optional recorder invoked as (x, weight, f(x)) for every node of every ACCEPTED
        /// interval — the intervals of the final composite rule. Null (the default) records
        /// nothing and leaves the integration path unchanged with zero overhead.
        /// </summary>
        /// <remarks>
        /// A naive per-evaluation weight hand-off double-counts under adaptivity: a rejected
        /// interval's 21 evaluations are superseded by its children's, so their weights must
        /// never carry measure. The recorder therefore fires only when an interval is accepted
        /// (tolerance reached, or the depth/evaluation caps force acceptance), with each node's
        /// Kronrod weight scaled by the interval half-length. Two identities follow: the weights
        /// of one accepted interval sum to exactly its width, so the recorded weights sum to the
        /// integration domain's width (any double count would overshoot it); and Σ weight·f(x)
        /// over all recorded nodes reproduces <see cref="Integrator.Result"/> to floating-point
        /// reassociation (the composite rule IS that sum). Recorded abscissas are strictly
        /// interior to their interval (the G10K21 property), so adjacent intervals never repeat
        /// a node.
        /// </remarks>
        public Action<double, double, double>? Recorder { get; set; }

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
                var (kronrodResult, gaussResult, nodes, values) = EvaluateGaussKronrod(a, b);

                // Recursively sub-divide
                Result = AdaptiveGK(Function, a, b, MaxDepth, kronrodResult, gaussResult, a, b, nodes, values);

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
                    var (kronrodResult, gaussResult, nodes, values) = EvaluateGaussKronrod(binA, binB);

                    // Recursively sub-divide
                    mu += AdaptiveGK(Function, binA, binB, MaxDepth, kronrodResult, gaussResult, binA, binB, nodes, values);
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
        /// Evaluates the Gauss-Kronrod G10K21 rule over the interval [a, b]. When a
        /// <see cref="Recorder"/> is attached, the 21 node abscissas and function values are
        /// also captured (in the fixed capture order matching the static weight layout) so the
        /// interval can flush them if it is later accepted; with no recorder the capture is
        /// skipped entirely and the evaluation path is unchanged.
        /// </summary>
        /// <param name="a">The lower bound of integration.</param>
        /// <param name="b">The upper bound of integration.</param>
        /// <returns>A tuple containing (Kronrod estimate, Gauss estimate, captured nodes, captured values).</returns>
        private (double kronrod, double gauss, double[]? nodes, double[]? values) EvaluateGaussKronrod(double a, double b)
        {
            double center = 0.5 * (a + b);
            double halfLength = 0.5 * (b - a);

            double resultGauss = 0.0;
            double resultKronrod = 0.0;
            double[]? nodes = Recorder != null ? new double[21] : null;
            double[]? values = Recorder != null ? new double[21] : null;

            // Evaluate at center point (x = 0)
            double f0 = Function(center);
            resultKronrod += wKronrod[10] * f0;
            FunctionEvaluations++;
            if (nodes != null)
            {
                nodes[0] = center;
                values![0] = f0;
            }

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
                if (nodes != null)
                {
                    nodes[2 * i + 1] = center - abscissa;
                    nodes[2 * i + 2] = center + abscissa;
                    values![2 * i + 1] = fval1;
                    values[2 * i + 2] = fval2;
                }
            }

            resultGauss *= halfLength;
            resultKronrod *= halfLength;

            return (resultKronrod, resultGauss, nodes, values);
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
        /// <param name="nodes">The interval's captured node abscissas (null when no recorder is attached).</param>
        /// <param name="values">The interval's captured function values (null when no recorder is attached).</param>
        /// <returns>
        /// An evaluation of the integral using adaptive Gauss-Kronrod with error less than the specified tolerance.
        /// This is accomplished by subdividing the interval until the error between the Gauss and Kronrod estimates
        /// is sufficiently small.
        /// </returns>
        private double AdaptiveGK(Func<double, double> f, double a, double b, int depth,
            double kronrodWhole, double gaussWhole, double a0, double b0, double[]? nodes, double[]? values)
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

                // The interval is ACCEPTED: it is part of the final composite rule, so its
                // nodes carry quadrature measure — flush them to the recorder with the
                // half-length-scaled Kronrod weights. Rejected (subdivided) intervals never
                // reach this point, so their superseded evaluations never carry weight.
                if (Recorder != null && nodes != null && values != null)
                {
                    double halfLength = 0.5 * (b - a);
                    for (int j = 0; j < 21; j++)
                    {
                        Recorder(nodes[j], wCapture[j] * halfLength, values[j]);
                    }
                }
                return kronrodWhole; // Return the more accurate Kronrod estimate
            }
            else
            {
                // Subdivide the interval at the midpoint
                double m = (a + b) / 2.0;

                // Evaluate Gauss-Kronrod on left half
                var (kronrodLeft, gaussLeft, nodesLeft, valuesLeft) = EvaluateGaussKronrod(a, m);

                // Evaluate Gauss-Kronrod on right half
                var (kronrodRight, gaussRight, nodesRight, valuesRight) = EvaluateGaussKronrod(m, b);

                // Recursively subdivide the intervals and accumulate results
                var leftResult = AdaptiveGK(f, a, m, depth - 1, kronrodLeft, gaussLeft, a0, b0, nodesLeft, valuesLeft);
                var rightResult = AdaptiveGK(f, m, b, depth - 1, kronrodRight, gaussRight, a0, b0, nodesRight, valuesRight);

                return leftResult + rightResult;
            }
        }
    }
}
