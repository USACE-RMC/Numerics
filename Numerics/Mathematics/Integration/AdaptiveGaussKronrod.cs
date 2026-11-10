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
    /// accurate for smooth functions. Two adaptive strategies are available through <see cref="Method"/>.
    /// </para>
    /// <para>
    /// <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/> (the default) bisects recursively and accepts an
    /// interval when the difference between its Gauss and Kronrod estimates meets the tolerance scaled by the
    /// interval's share of the range.
    /// </para>
    /// <para>
    /// <see cref="AdaptiveGaussKronrodMethod.Qags"/> is the QUADPACK QAGS algorithm. Each interval carries the
    /// QUADPACK error estimate resasc * min(1, (200 |K - G| / resasc)^1.5), floored at 50 eps times the integral
    /// of |f|; the intervals are kept in a priority queue by error estimate and the one with the largest error is
    /// always bisected next. When that interval becomes small, the sequence of composite sums is accelerated with
    /// Wynn's epsilon algorithm (QUADPACK QELG), which recovers the limit for integrands with endpoint or interior
    /// algebraic and logarithmic singularities. The roundoff, extrapolation and divergence tests are those of
    /// QUADPACK DQAGSE, and integration stops when the summed error estimate is at most
    /// max(<see cref="Integrator.AbsoluteTolerance"/>, <see cref="Integrator.RelativeTolerance"/> |result|).
    /// Departures from DQAGSE keep the class contract: the subinterval limit is derived from
    /// <see cref="Integrator.MaxFunctionEvaluations"/>; <see cref="MinDepth"/> starts the integration from a uniform
    /// partition; the tolerance test, and the roundoff test, wait until <see cref="Integrator.MinFunctionEvaluations"/>
    /// evaluations are made; an interval whose bisection level has reached <see cref="MaxDepth"/>, or whose width
    /// is at the floating-point resolution test of DQAGSE, is retired instead of ending the integration; and an
    /// interval whose estimate is NaN or infinite is bisected depth first until its pieces are finite. A
    /// non-integrable singularity that lies exactly on a bisection point can then be integrated as a principal
    /// value and reported as <see cref="IntegrationStatus.Success"/>: for 1/(x - 0.5) on [0, 1] the result is
    /// about 1E-13. The rule sums are accumulated in the order of QUADPACK DQK21.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><see href="https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula"/></item>
    /// <item>Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., and Kahaner, D. K. (1983). QUADPACK: A Subroutine Package for Automatic Integration. Springer-Verlag. <see href="https://doi.org/10.1007/978-3-642-61786-7"/></item>
    /// <item>Wynn, P. (1956). On a device for computing the e_m(S_n) transformation. Mathematical Tables and Other Aids to Computation, 10(54), 91-96. <see href="https://doi.org/10.2307/2002183"/></item>
    /// <item>QUADPACK source: DQAGS, DQAGSE, DQELG and DQK21, <see href="https://www.netlib.org/quadpack/"/></item>
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
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a limit is not finite or <paramref name="max"/> is not greater than <paramref name="min"/>.</exception>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="function"/> is null.</exception>
        public AdaptiveGaussKronrod(Func<double, double> function, double min, double max)
        {
            if (!Tools.IsFinite(min)) throw new ArgumentOutOfRangeException(nameof(min), "The minimum value must be finite.");
            if (!Tools.IsFinite(max)) throw new ArgumentOutOfRangeException(nameof(max), "The maximum value must be finite.");
            if (max <= min) throw new ArgumentOutOfRangeException(nameof(max), "The maximum value cannot be less than or equal to the minimum value.");
            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            a = min;
            b = max;
        }

        private double a, b, _squaredError;

        // Number of intervals accepted at the depth or width limit without meeting the tolerance.
        private int _forcedAcceptances;


        // QUADPACK machine constants: d1mach(4) = 2^-52, d1mach(1) and d1mach(2).
        private const double Epmach = 2.220446049250313E-16;
        private const double Uflow = 2.2250738585072014E-308;
        private const double Oflow = double.MaxValue;

        // Maximum number of elements of the epsilon table (QUADPACK limexp).
        private const int EpsilonLimit = 50;

        // Gauss-Kronrod G10K21 nodes and weights (for interval [-1, 1])
        // 10-point Gauss weights (symmetric, so only the positive values are stored)
        private static readonly double[] wGauss = new double[]
        {
            0.066671344308688137593568809893332,
            0.149451349150580593145776339657697,
            0.219086362515982043995534934228163,
            0.269266719309996355091226921569469,
            0.295524224714752870173892994651338
        };

        // 21-point Kronrod nodes (symmetric, so only the positive values are stored)
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

        // The unscaled Kronrod weight per node slot: slot 0 is the center node (wKronrod[10]),
        // then the symmetric pair for each abscissa index i occupies slots 2i+1 and 2i+2.
        private static readonly double[] wCapture = BuildCaptureWeights();

        /// <summary>
        /// The QK21 result for one interval.
        /// </summary>
        [Serializable]
        private struct Estimate
        {
            /// <summary>The Kronrod estimate of the integral.</summary>
            public double Area;
            /// <summary>The QUADPACK error estimate.</summary>
            public double Error;
            /// <summary>The Kronrod estimate of the integral of |f| (QUADPACK resabs).</summary>
            public double AbsArea;
            /// <summary>The Kronrod estimate of the integral of |f - mean| (QUADPACK resasc).</summary>
            public double AscArea;

            /// <summary>True when the estimate and its error are finite.</summary>
            public readonly bool IsFinite => Tools.IsFinite(Area) && Tools.IsFinite(Error);
        }

        /// <summary>
        /// The unscaled DQK21 sums accumulated while the nodes are evaluated.
        /// </summary>
        [Serializable]
        private struct RuleSums
        {
            /// <summary>The Gauss sum (QUADPACK resg).</summary>
            public double Gauss;
            /// <summary>The Kronrod sum (QUADPACK resk).</summary>
            public double Kronrod;
            /// <summary>The Kronrod sum of |f| (QUADPACK resabs).</summary>
            public double Abs;
        }

        /// <summary>
        /// An entry of the error queue.
        /// </summary>
        [Serializable]
        private struct QueueEntry
        {
            /// <summary>The error estimate of the interval.</summary>
            public double Error;
            /// <summary>The sequence number of the interval.</summary>
            public int Sequence;
            /// <summary>The interval slot.</summary>
            public int Slot;
        }

        /// <summary>
        /// The outcome of one QAGS run over a single range.
        /// </summary>
        [Serializable]
        private struct Outcome
        {
            /// <summary>The integral estimate.</summary>
            public double Result;
            /// <summary>The error estimate.</summary>
            public double Error;
            /// <summary>The number of intervals in the final partition.</summary>
            public int Intervals;
            /// <summary>True when the interval limit derived from the evaluation budget stopped the run.</summary>
            public bool LimitReached;
            /// <summary>True when the run stopped without meeting the tolerance for any other reason.</summary>
            public bool NotConverged;
        }

        // Node values of the two halves evaluated last, slots 0-20 and 21-41.
        private double[] _scratch = Array.Empty<double>();

        // The interval list, interleaved four doubles per slot: the QUADPACK alist, blist, rlist and elist
        // entries (lower limit, upper limit, estimate, error estimate).
        private double[] _intervals = Array.Empty<double>();

        // Two ints per slot: the bisection level and the sequence number (a later store wins ties in the queue).
        private int[] _slots = Array.Empty<int>();

        // The binary max-heap of intervals ordered by error estimate, and the intervals set aside while the
        // extrapolation step passes over small intervals. Each entry carries its error and sequence number.
        private QueueEntry[] _queue = Array.Empty<QueueEntry>();
        private QueueEntry[] _skipped = Array.Empty<QueueEntry>();

        // The node values kept for the recorder, 21 per slot.
        private double[] _values = Array.Empty<double>();
        private int _capacity, _count, _sequenceCounter, _heapCount, _skippedCount;

        // The epsilon table (QUADPACK rlist2) and the last three extrapolated results (res3la).
        private double[] _epsilonTable = Array.Empty<double>();
        private double[] _lastResults = Array.Empty<double>();

        // Set when the subinterval limit stops the bisection of an interval with a non-finite estimate.
        private bool _budgetExhausted;

        /// <summary>Gets the lower limit of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref double Lower(int s) => ref _intervals[4 * s];

        /// <summary>Gets the upper limit of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref double Upper(int s) => ref _intervals[4 * s + 1];

        /// <summary>Gets the estimate of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref double AreaOf(int s) => ref _intervals[4 * s + 2];

        /// <summary>Gets the error estimate of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref double ErrorOf(int s) => ref _intervals[4 * s + 3];

        /// <summary>Gets the bisection level of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref int Level(int s) => ref _slots[2 * s];

        /// <summary>Gets the sequence number of the interval in a slot.</summary>
        /// <param name="s">The slot.</param>
        /// <returns>A reference to the value.</returns>
        private ref int Sequence(int s) => ref _slots[2 * s + 1];

        /// <summary>Gets the slot at the top of the queue.</summary>
        private int Top => _queue[0].Slot;

        /// <summary>
        /// Per-slot capture buffers for the recorder's node abscissas, grown on demand. Only
        /// allocated when <see cref="Recorder"/> is set.
        /// </summary>
        private double[]?[] _nodePool = Array.Empty<double[]?>();

        /// <summary>
        /// Per-slot capture buffers for the recorder's node values.
        /// </summary>
        private double[]?[] _valuePool = Array.Empty<double[]?>();

        /// <summary>
        /// Returns the capture buffer for a recursion slot, growing the pool as needed.
        /// </summary>
        /// <param name="pool">The pool to rent from.</param>
        /// <param name="slot">The slot index; distinct for every interval alive at once.</param>
        /// <returns>The buffer.</returns>
        private static double[] RentCapture(ref double[]?[] pool, int slot)
        {
            if (slot >= pool.Length)
            {
                var grown = new double[]?[Math.Max(slot + 1, pool.Length * 2 + 8)];
                Array.Copy(pool, grown, pool.Length);
                pool = grown;
            }
            return pool[slot] ??= new double[21];
        }

        /// <summary>
        /// Builds the unscaled Kronrod weight layout matching the node order.
        /// </summary>
        /// <returns>The 21 unscaled weights in node order.</returns>
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
        /// Gets or sets the adaptive strategy. Default = <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/>.
        /// </summary>
        /// <remarks>
        /// The recorder, depth limits, stratification bins, evaluation budget and status are honoured by both
        /// methods; their precise meaning under each is given with the member that uses them. The two methods
        /// accept different intervals, so their results and recorded nodes differ.
        /// </remarks>
        public AdaptiveGaussKronrodMethod Method { get; set; } = AdaptiveGaussKronrodMethod.LocalBisection;

        /// <summary>
        /// Gets or sets a callback invoked as (x, weight, f(x)) for the quadrature nodes of the final partition.
        /// A null callback disables recording.
        /// </summary>
        /// <remarks>
        /// Evaluations from bisected intervals are omitted because their halves replace them in the composite
        /// rule. Under <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/> the nodes are reported as each
        /// interval is accepted; under <see cref="AdaptiveGaussKronrodMethod.Qags"/> they are reported after
        /// integration, interval by interval from left to right. Each
        /// reported Kronrod weight includes the interval half-length, so the weights sum to the
        /// integration-domain width and their weighted function values reproduce the composite Kronrod sum up
        /// to floating-point reassociation. That sum is <see cref="Integrator.Result"/> unless the QAGS
        /// extrapolated value was returned; the composite sum is then the less accurate of the two, because extrapolation
        /// stops refinement before the composite rule converges (for x^-0.5 on [0, 1] the result is 2 and the
        /// composite sum 1.994). Recording does not change the result or the evaluation count. The callback is
        /// snapshotted when integration begins.
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
        /// <remarks>
        /// Under <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/> no interval is accepted above this depth;
        /// under <see cref="AdaptiveGaussKronrodMethod.Qags"/> integration starts from a uniform partition into
        /// 2^MinDepth intervals, limited by <see cref="MaxDepth"/> and by the interval count the evaluation budget
        /// allows.
        /// </remarks>
        public int MinDepth { get; set; } = 0;

        /// <summary>
        /// The maximum recursion depth. Default = 100.
        /// </summary>
        /// <remarks>
        /// An interval produced by this many successive bisections of the original range is not bisected further.
        /// </remarks>
        public int MaxDepth { get; set; } = 100;

        /// <summary>
        /// Returns an approximate measure of the standard error of the integration.
        /// </summary>
        /// <remarks>
        /// Under <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/> it is the square root of the sum of the
        /// squared Kronrod-Gauss differences of the accepted intervals. Under
        /// <see cref="AdaptiveGaussKronrodMethod.Qags"/> it is the QUADPACK absolute error estimate of the returned
        /// result, and for stratified integration the square root of the sum of the squared estimates of the bins.
        /// </remarks>
        public double StandardError { get; private set; }

        /// <inheritdoc/>
        /// <remarks>
        /// <see cref="AdaptiveGaussKronrodMethod.LocalBisection"/> applies the base rule, both tolerances in
        /// [1E-15, 1]. <see cref="AdaptiveGaussKronrodMethod.Qags"/> applies the QUADPACK rule: both tolerances in
        /// [0, 1], and when <see cref="Integrator.AbsoluteTolerance"/> is zero the relative tolerance must be at
        /// least max(50 eps, 5E-29) with eps = 2^-52, so a positive absolute tolerance allows a zero relative
        /// tolerance.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when an iteration or evaluation limit is less than one,
        /// a tolerance is NaN, negative or greater than one, or both tolerances are below the QUADPACK floor.</exception>
        protected override void Validate()
        {
            if (Method != AdaptiveGaussKronrodMethod.Qags)
            {
                base.Validate();
                return;
            }
            if (MinIterations < 1) throw new ArgumentOutOfRangeException(nameof(MinIterations), "The minimum number of iterations must be greater than or equal to 1.");
            if (MinFunctionEvaluations < 1) throw new ArgumentOutOfRangeException(nameof(MinFunctionEvaluations), "The minimum number of function evaluations must be greater than or equal to 1.");
            if (MaxIterations < 1) throw new ArgumentOutOfRangeException(nameof(MaxIterations), "The maximum number of iterations must be greater than 1.");
            if (MaxFunctionEvaluations < 1) throw new ArgumentOutOfRangeException(nameof(MaxFunctionEvaluations), "The maximum number of function evaluations must be greater than 1.");
            if (!(RelativeTolerance >= 0 && RelativeTolerance <= 1)) throw new ArgumentOutOfRangeException(nameof(RelativeTolerance), "The relative tolerance must be between 0 and 1.");
            if (!(AbsoluteTolerance >= 0 && AbsoluteTolerance <= 1)) throw new ArgumentOutOfRangeException(nameof(AbsoluteTolerance), "The absolute tolerance must be between 0 and 1.");
            if (AbsoluteTolerance == 0 && RelativeTolerance < Math.Max(50 * Epmach, 5E-29))
                throw new ArgumentOutOfRangeException(nameof(RelativeTolerance), "When the absolute tolerance is zero, the relative tolerance must be at least max(50 eps, 5E-29).");
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// <b> LocalBisection: </b>
        /// Refinement is bounded by <see cref="Integrator.MaxFunctionEvaluations"/> and <see cref="MaxDepth"/>; the
        /// inherited <see cref="Integrator.MinIterations"/> and <see cref="Integrator.MaxIterations"/> are not
        /// consulted. <see cref="Integrator.Status"/> is <see cref="IntegrationStatus.MaximumFunctionEvaluationsReached"/>
        /// when the evaluation budget is exhausted, <see cref="IntegrationStatus.MaximumIterationsReached"/> when any
        /// interval was accepted at the depth or width limit without meeting the tolerance, and
        /// <see cref="IntegrationStatus.Success"/> otherwise. If an interval with a NaN or infinite estimate would be
        /// accepted, the result cannot be finite, so integration stops with <see cref="IntegrationStatus.Failure"/>
        /// and throws an <see cref="ArithmeticException"/> when <see cref="Integrator.ReportFailure"/> is true. A
        /// non-finite estimate on an interval that is subdivided further does not stop integration, because its
        /// children may be finite.
        /// </para>
        /// <para>
        /// <b> Qags: </b>
        /// The evaluation budget <see cref="Integrator.MaxFunctionEvaluations"/> sets the QUADPACK subinterval
        /// limit, floor((floor(MaxFunctionEvaluations / 21) + 1) / 2), so the evaluation count does not exceed
        /// the budget once the initial partition is evaluated. <see cref="Integrator.Iterations"/> is the number
        /// of bisections; the inherited <see cref="Integrator.MinIterations"/> and <see cref="Integrator.MaxIterations"/>
        /// are not consulted. Integration does not stop on the tolerance before
        /// <see cref="Integrator.MinFunctionEvaluations"/> evaluations.
        /// </para>
        /// <para>
        /// <see cref="Integrator.Status"/> is <see cref="IntegrationStatus.MaximumFunctionEvaluationsReached"/>
        /// when the evaluation budget is exhausted, <see cref="IntegrationStatus.MaximumIterationsReached"/>
        /// when the tolerance was not met for another reason (an interval retired at the depth or width limit,
        /// roundoff, extrapolation failure or a divergent integral, QUADPACK codes 2 to 5), and
        /// <see cref="IntegrationStatus.Success"/> otherwise. If an interval with a NaN or infinite estimate cannot
        /// be bisected further, the result cannot be finite, so integration stops with
        /// <see cref="IntegrationStatus.Failure"/> and throws an <see cref="ArithmeticException"/> when
        /// <see cref="Integrator.ReportFailure"/> is true.
        /// </para>
        /// </remarks>
        public override void Integrate()
        {
            if (Method != AdaptiveGaussKronrodMethod.Qags)
            {
                IntegrateLocal();
                return;
            }
            StandardError = 0;
            ClearResults();
            Validate();
            Action<double, double, double>? recorder = Recorder;

            try
            {
                var outcome = Qags(a, b, IntervalLimit(), recorder);
                Result = outcome.Result;
                StandardError = outcome.Error;
                Status = StatusOf(outcome.LimitReached, outcome.NotConverged);
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
        /// <remarks>
        /// Each bin is integrated separately to the tolerance, in order, and the evaluation budget is shared.
        /// Under <see cref="AdaptiveGaussKronrodMethod.Qags"/> a bin may use the subintervals the budget allows
        /// less one for each bin still to come. The status rules match <see cref="Integrate()"/>.
        /// </remarks>
        public void Integrate(List<StratificationBin> bins)
        {
            if (Method != AdaptiveGaussKronrodMethod.Qags)
            {
                IntegrateLocal(bins);
                return;
            }
            StandardError = 0;
            ClearResults();
            Validate();
            Action<double, double, double>? recorder = Recorder;

            try
            {
                int limit = IntervalLimit();
                int used = 0;
                bool limitReached = false, notConverged = false;
                double mu = 0, squaredError = 0;
                for (int i = 0; i < bins.Count; i++)
                {
                    int binLimit = Math.Max(1, limit - used - (bins.Count - 1 - i));
                    var outcome = Qags(bins[i].LowerBound, bins[i].UpperBound, binLimit, recorder);
                    mu += outcome.Result;
                    squaredError += outcome.Error * outcome.Error;
                    used += outcome.Intervals;
                    limitReached |= outcome.LimitReached;
                    notConverged |= outcome.NotConverged;
                }

                Result = mu;
                StandardError = Math.Sqrt(squaredError);
                Status = StatusOf(limitReached, notConverged);
            }
            catch (Exception)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw;
            }
        }

        /// <summary>
        /// Integrates [a, b] with recursive local bisection.
        /// </summary>
        private void IntegrateLocal()
        {
            _squaredError = 0;
            _forcedAcceptances = 0;
            StandardError = 0;
            ClearResults();
            Validate();
            Action<double, double, double>? recorder = Recorder;

            try
            {
                // Initial evaluation using Gauss-Kronrod rule on the whole interval
                var (kronrodResult, gaussResult, nodes, values) = EvaluateGaussKronrod(a, b, 0, recorder);

                // Recursively sub-divide
                Result = AdaptiveGK(Function, a, b, MaxDepth, kronrodResult, gaussResult, a, b, nodes, values, 0, recorder);

                // Standard error calculated after recursion completes
                StandardError = Math.Sqrt(_squaredError);

                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else if (_forcedAcceptances > 0)
                {
                    Status = IntegrationStatus.MaximumIterationsReached;
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
        /// Integrates the stratification bins with recursive local bisection.
        /// </summary>
        /// <param name="bins">The stratification bins to integrate over.</param>
        private void IntegrateLocal(List<StratificationBin> bins)
        {
            _squaredError = 0;
            _forcedAcceptances = 0;
            StandardError = 0;
            ClearResults();
            Validate();
            Action<double, double, double>? recorder = Recorder;

            try
            {
                double mu = 0;
                for (int i = 0; i < bins.Count; i++)
                {
                    // Initial evaluation using Gauss-Kronrod rule on the bin interval
                    double binA = bins[i].LowerBound;
                    double binB = bins[i].UpperBound;
                    var (kronrodResult, gaussResult, nodes, values) = EvaluateGaussKronrod(binA, binB, 0, recorder);

                    // Recursively sub-divide
                    mu += AdaptiveGK(Function, binA, binB, MaxDepth, kronrodResult, gaussResult, binA, binB, nodes, values, 0, recorder);
                }

                // Final result and standard error
                Result = mu;
                StandardError = Math.Sqrt(_squaredError);

                if (FunctionEvaluations >= MaxFunctionEvaluations)
                {
                    Status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                }
                else if (_forcedAcceptances > 0)
                {
                    Status = IntegrationStatus.MaximumIterationsReached;
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
        /// interval can report them if it is accepted; with no recorder the capture is
        /// skipped entirely and the evaluation path is unchanged.
        /// </summary>
        /// <param name="a">The lower bound of integration.</param>
        /// <param name="b">The upper bound of integration.</param>
        /// <param name="slot">The capture-buffer slot.</param>
        /// <param name="recorder">The recorder snapshot for this integration.</param>
        /// <returns>A tuple containing (Kronrod estimate, Gauss estimate, captured nodes, captured values).</returns>
        private (double kronrod, double gauss, double[]? nodes, double[]? values) EvaluateGaussKronrod(double a, double b, int slot, Action<double, double, double>? recorder)
        {
            double center = 0.5 * (a + b);
            double halfLength = 0.5 * (b - a);

            double resultGauss = 0.0;
            double resultKronrod = 0.0;
            double[]? nodes = recorder != null ? RentCapture(ref _nodePool, slot) : null;
            double[]? values = recorder != null ? RentCapture(ref _valuePool, slot) : null;

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
        /// <param name="level">The recursion level used to select capture slots.</param>
        /// <param name="recorder">The recorder snapshot for this integration.</param>
        /// <returns>
        /// An evaluation of the integral using adaptive Gauss-Kronrod with error less than the specified tolerance.
        /// This is accomplished by subdividing the interval until the error between the Gauss and Kronrod estimates
        /// is sufficiently small.
        /// </returns>
        private double AdaptiveGK(Func<double, double> f, double a, double b, int depth,
            double kronrodWhole, double gaussWhole, double a0, double b0, double[]? nodes, double[]? values, int level,
            Action<double, double, double>? recorder)
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
                // Once a non-finite estimate enters the composite rule the total cannot be finite, so stop
                // instead of spending the remaining budget on other intervals.
                if (!Tools.IsFinite(kronrodWhole))
                    throw new ArithmeticException("The integrand produced a NaN or infinite estimate on an accepted interval, so the integral cannot be finite.");

                // Count intervals accepted only because a depth, width or budget limit was reached.
                if (!absoluteToleranceReached && !relativeToleranceReached) _forcedAcceptances++;

                // Convergence is reached
                _squaredError += error * error; // Accumulate squared errors

                // Report only nodes that belong to the final composite rule, using weights
                // scaled by the accepted interval's half-length.
                if (recorder != null && nodes != null && values != null)
                {
                    double halfLength = 0.5 * (b - a);
                    for (int j = 0; j < 21; j++)
                    {
                        recorder(nodes[j], wCapture[j] * halfLength, values[j]);
                    }
                }
                return kronrodWhole; // Return the more accurate Kronrod estimate
            }
            else
            {
                // Subdivide the interval at the midpoint
                double m = (a + b) / 2.0;

                // Evaluate Gauss-Kronrod on left half
                // Both halves are evaluated before either recurses, so they need distinct capture
                // slots; their own children take the slots of the next level, which are free by
                // the time each subtree runs.
                var (kronrodLeft, gaussLeft, nodesLeft, valuesLeft) = EvaluateGaussKronrod(a, m, 2 * level + 1, recorder);

                // Evaluate Gauss-Kronrod on right half
                var (kronrodRight, gaussRight, nodesRight, valuesRight) = EvaluateGaussKronrod(m, b, 2 * level + 2, recorder);

                // Recursively subdivide the intervals and accumulate results
                var leftResult = AdaptiveGK(f, a, m, depth - 1, kronrodLeft, gaussLeft, a0, b0, nodesLeft, valuesLeft, level + 1, recorder);
                var rightResult = AdaptiveGK(f, m, b, depth - 1, kronrodRight, gaussRight, a0, b0, nodesRight, valuesRight, level + 1, recorder);

                return leftResult + rightResult;
            }
        }
        /// <summary>
        /// Returns the QUADPACK subinterval limit implied by the evaluation budget.
        /// </summary>
        /// <returns>The maximum number of subintervals, at least one.</returns>
        private int IntervalLimit()
        {
            return Math.Max(1, (MaxFunctionEvaluations / 21 + 1) / 2);
        }

        /// <summary>
        /// Maps the run flags to the integration status.
        /// </summary>
        /// <param name="limitReached">True when the subinterval limit stopped a run.</param>
        /// <param name="notConverged">True when a run stopped without meeting the tolerance for another reason.</param>
        /// <returns>The status.</returns>
        private IntegrationStatus StatusOf(bool limitReached, bool notConverged)
        {
            if (limitReached || FunctionEvaluations >= MaxFunctionEvaluations) return IntegrationStatus.MaximumFunctionEvaluationsReached;
            if (notConverged) return IntegrationStatus.MaximumIterationsReached;
            return IntegrationStatus.Success;
        }

        /// <summary>
        /// Integrates over [lower, upper] with the QUADPACK DQAGSE algorithm.
        /// </summary>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="limit">The maximum number of subintervals.</param>
        /// <param name="recorder">The recorder snapshot for this integration.</param>
        /// <returns>The result, error estimate and termination flags.</returns>
        /// <exception cref="ArithmeticException">Thrown when an interval with a non-finite estimate cannot be bisected.</exception>
        /// <remarks>
        /// The control flow follows DQAGSE; the labels of the Fortran source are noted where the flow jumps.
        /// The descending error list of DQPSRT is replaced by a binary heap, and the intervals passed over by
        /// the extrapolation step (the positions before nrmax in DQAGSE) are held aside until the next
        /// extrapolation.
        /// </remarks>
        private Outcome Qags(double lower, double upper, int limit, Action<double, double, double>? recorder)
        {
            _count = 0;
            _heapCount = 0;
            _skippedCount = 0;
            _sequenceCounter = 0;
            _budgetExhausted = false;
            bool record = recorder != null;
            if (_scratch.Length == 0) _scratch = new double[42];
            double epsabs = AbsoluteTolerance, epsrel = RelativeTolerance;

            // First approximation: the whole range, or the uniform MinDepth partition of it.
            int split = Math.Max(0, Math.Min(Math.Min(MinDepth, MaxDepth), 30));
            while (split > 0 && (1 << split) > limit) split--;

            double result, abserr, defabs, resabs;
            double sumArea = 0, sumError = 0, sumAbs = 0, sumAsc = 0;
            bool rootTested = false;
            if (split == 0)
            {
                var root = Kronrod21(lower, upper, 0);
                if (root.IsFinite)
                {
                    result = root.Area;
                    abserr = root.Error;
                    defabs = root.AbsArea;
                    resabs = root.AscArea;
                    double rootBound = Math.Max(epsabs, epsrel * Math.Abs(result));
                    bool rootConverged = FunctionEvaluations >= MinFunctionEvaluations &&
                        ((abserr <= rootBound && abserr != resabs) || abserr == 0);
                    bool rootRoundoff = abserr <= 100 * Epmach * defabs && abserr > rootBound;
                    if (rootConverged || rootRoundoff || limit == 1)
                    {
                        // DQAGSE label 140: the single-interval rule is the answer.
                        if (record) ReportInterval(recorder!, lower, upper, _scratch, 0);
                        return new Outcome
                        {
                            Result = result,
                            Error = abserr,
                            Intervals = 1,
                            LimitReached = !rootConverged && !rootRoundoff,
                            NotConverged = !rootConverged && rootRoundoff
                        };
                    }
                    PlaceFinite(Allocate(), lower, upper, 0, root, 0, record, ref sumArea, ref sumError);
                    rootTested = true;
                }
                else
                {
                    int slot = Allocate();
                    Store(slot, lower, upper, 0, root, 0, record);
                    Resolve(slot, lower, upper, 0, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
                    result = sumArea;
                    abserr = sumError;
                    defabs = sumAbs;
                    resabs = sumAsc;
                }
            }
            else
            {
                Partition(lower, upper, 0, split, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
                result = sumArea;
                abserr = sumError;
                defabs = sumAbs;
                resabs = sumAsc;
            }

            // The subinterval limit stopped the bisection of a non-finite interval.
            if (_budgetExhausted) return Composite(sumError, 1, false, recorder);

            // Test on accuracy for a first approximation made of a partition.
            double dres = Math.Abs(result);
            double errbnd = Math.Max(epsabs, epsrel * dres);
            int ier = 0;
            if (!rootTested)
            {
                bool converged = FunctionEvaluations >= MinFunctionEvaluations &&
                    ((abserr <= errbnd && abserr != resabs) || abserr == 0);
                if (converged) return Composite(abserr, 0, false, recorder);
                if (abserr <= 100 * Epmach * defabs && abserr > errbnd) ier = 2;
                if (_count >= limit) ier = 1;
                if (ier != 0) return Composite(abserr, ier, false, recorder);
            }

            // Initialization.
            double firstSum = result, secondSum = 0;
            bool tableReady = false;
            double area = result;
            double errsum = abserr;
            abserr = Oflow;
            int nres = 0;
            int numrl2 = 2;
            int ktmin = 0;
            bool extrap = false, noext = false;
            int iroff1 = 0, iroff2 = 0, iroff3 = 0, ierro = 0;
            int ksgn = dres >= (1 - 50 * Epmach) * defabs ? 1 : -1;
            double small = 0, erlarg = 0, ertest = 0, correc = 0;
            double retiredError = 0;
            bool retired = false, forced = false;
            int bisections = 0;
            double smallStart = 0.375 * Math.Abs(upper - lower);
            if (!rootTested)
            {
                // A partition stands where DQAGSE is after its first bisection (label 80): small is 0.75 of
                // the widest piece, between the widths of that level and the next, and the composite sum is
                // the first element of the epsilon table.
                double widest = 0;
                for (int k = 0; k < _count; k++) widest = Math.Max(widest, Width(k));
                small = 0.75 * widest;
                erlarg = errsum;
                ertest = errbnd;
                numrl2 = 1;
                bisections = 1;
            }

            while (true)
            {
                // Retire intervals that cannot be bisected; stop when the error that can still be
                // reduced meets the tolerance.
                while (_heapCount > 0 && !CanBisect(Top))
                {
                    int s = PopTop();
                    retiredError += ErrorOf(s);
                    retired = true;
                }
                if (_heapCount == 0)
                {
                    if (_skippedCount > 0)
                    {
                        RestoreSkipped();
                        continue;
                    }
                    forced = true;
                    break;
                }
                if (retired && errsum - retiredError <= errbnd && FunctionEvaluations >= MinFunctionEvaluations)
                {
                    forced = true;
                    break;
                }

                // Bisect the subinterval with the largest error estimate.
                int maxerr = PopTop();
                double errmax = ErrorOf(maxerr);
                double a1 = Lower(maxerr);
                double b2 = Upper(maxerr);
                double b1 = 0.5 * (a1 + b2);
                double a2 = b1;
                double erlast = errmax;
                int level = Level(maxerr) + 1;
                Bisect(a1, b2, out var left, out var right);
                bisections++;
                Iterations++;

                double area12, erro12;
                if (left.IsFinite && right.IsFinite)
                {
                    // Improve previous approximations to integral and error and test for accuracy.
                    area12 = left.Area + right.Area;
                    erro12 = left.Error + right.Error;
                    errsum = errsum + erro12 - errmax;
                    area = area + area12 - AreaOf(maxerr);
                    // Refinement that MinFunctionEvaluations forces past convergence is not evidence of roundoff,
                    // so the roundoff counters start once the minimum is met.
                    if (left.AscArea != left.Error && right.AscArea != right.Error && FunctionEvaluations >= MinFunctionEvaluations)
                    {
                        if (!(Math.Abs(AreaOf(maxerr) - area12) > 1E-5 * Math.Abs(area12) || erro12 < 0.99 * errmax))
                        {
                            if (extrap) iroff2++;
                            else iroff1++;
                        }
                        if (_count + 1 > 10 && erro12 > errmax) iroff3++;
                    }

                    // The half with the larger error takes the parent's slot, as in DQAGSE; it also wins ties
                    // in the queue, matching the insertion order of DQPSRT.
                    int last = Allocate();
                    if (right.Error > left.Error)
                    {
                        Store(last, a1, b1, level, left, 0, record);
                        Store(maxerr, a2, b2, level, right, 21, record);
                    }
                    else
                    {
                        Store(last, a2, b2, level, right, 21, record);
                        Store(maxerr, a1, b1, level, left, 0, record);
                    }
                    Push(last);
                    Push(maxerr);
                }
                else
                {
                    // A half with a non-finite estimate is bisected until its pieces are finite; finite halves
                    // are stored first because resolving reuses the node-value scratch.
                    double oldArea = AreaOf(maxerr);
                    double halfArea = 0, halfError = 0, halfAbs = 0, halfAsc = 0;
                    int last = Allocate();
                    if (left.IsFinite) PlaceFinite(maxerr, a1, b1, level, left, 0, record, ref halfArea, ref halfError);
                    else Store(maxerr, a1, b1, level, left, 0, record);
                    if (right.IsFinite) PlaceFinite(last, a2, b2, level, right, 21, record, ref halfArea, ref halfError);
                    else Store(last, a2, b2, level, right, 21, record);
                    if (!left.IsFinite) Resolve(maxerr, a1, b1, level, limit, record, ref halfArea, ref halfError, ref halfAbs, ref halfAsc);
                    if (!right.IsFinite) Resolve(last, a2, b2, level, limit, record, ref halfArea, ref halfError, ref halfAbs, ref halfAsc);
                    area12 = halfArea;
                    erro12 = halfError;
                    errsum = errsum + erro12 - errmax;
                    area = area + area12 - oldArea;
                    if (_budgetExhausted) return Composite(errsum, 1, false, recorder);
                }

                errbnd = Math.Max(epsabs, epsrel * Math.Abs(area));

                // Test for roundoff error and eventually set error flag.
                if (iroff1 + iroff2 >= 10 || iroff3 >= 20) ier = 2;
                if (iroff2 >= 5) ierro = 3;

                // Set error flag in the case that the number of subintervals equals limit.
                if (_count >= limit) ier = 1;

                bool minimumMet = FunctionEvaluations >= MinFunctionEvaluations;
                if (errsum <= errbnd && minimumMet) return Composite(errsum, 0, false, recorder);
                if (ier != 0) break;
                if (bisections == 1)
                {
                    // DQAGSE label 80: small is 0.375 of the range, between the widths of the first and second levels.
                    small = smallStart;
                    erlarg = errsum;
                    ertest = errbnd;
                    secondSum = area;
                    continue;
                }
                if (noext) continue;
                erlarg -= erlast;
                if (Math.Abs(b1 - a1) > small) erlarg += erro12;
                if (!extrap)
                {
                    // Test whether the interval to be bisected next is the smallest interval.
                    if (_heapCount == 0 || Width(Top) > small) continue;
                    extrap = true;
                    SkipTop();
                }

                // DQAGSE label 40.
                bool extrapolate = true;
                if (ierro != 3 && erlarg > ertest)
                {
                    // The smallest interval has the largest error. Before bisecting, decrease the sum of the
                    // errors over the larger intervals (erlarg) and perform extrapolation.
                    while (_heapCount > 0)
                    {
                        if (Width(Top) > small)
                        {
                            extrapolate = false;
                            break;
                        }
                        SkipTop();
                    }
                }
                if (!extrapolate) continue;

                // Perform extrapolation (DQAGSE label 60).
                if (!tableReady)
                {
                    // The epsilon table is allocated when extrapolation first runs.
                    if (_epsilonTable.Length == 0)
                    {
                        _epsilonTable = new double[EpsilonLimit + 2];
                        _lastResults = new double[3];
                    }
                    _epsilonTable[0] = firstSum;
                    _epsilonTable[1] = secondSum;
                    Array.Clear(_lastResults, 0, _lastResults.Length);
                    tableReady = true;
                }
                numrl2++;
                _epsilonTable[numrl2 - 1] = area;

                // DQAGSE always reaches this point with at least three elements in the table. A first
                // approximation made of a partition starts the table with one element, so its first composite
                // sum is only stored, as DQAGSE stores its second at label 80.
                if (numrl2 >= 3)
                {
                    Extrapolate(ref numrl2, out double reseps, out double abseps, ref nres);
                    ktmin++;
                    if (ktmin > 5 && abserr < 1E-3 * errsum) ier = 5;
                    if (abseps < abserr)
                    {
                        ktmin = 0;
                        abserr = abseps;
                        result = reseps;
                        correc = erlarg;
                        ertest = Math.Max(epsabs, epsrel * Math.Abs(reseps));
                        if (abserr <= ertest && minimumMet) break;
                    }
                }

                // Prepare bisection of the smallest interval (DQAGSE label 70).
                if (numrl2 == 1) noext = true;
                if (ier == 5) break;
                RestoreSkipped();
                extrap = false;
                small *= 0.5;
                erlarg = errsum;
            }

            // Set final result and error estimate (DQAGSE label 100).
            if (abserr == Oflow) return Composite(errsum, ier, forced, recorder);
            if (ier + ierro != 0)
            {
                if (ierro == 3) abserr += correc;
                if (ier == 0) ier = 3;
                if (result != 0 && area != 0)
                {
                    if (abserr / Math.Abs(result) > errsum / Math.Abs(area)) return Composite(errsum, ier, forced, recorder);
                }
                else if (abserr > errsum)
                {
                    return Composite(errsum, ier, forced, recorder);
                }
                else if (area == 0)
                {
                    return Finish(result, abserr, ier, forced, recorder);
                }
            }

            // Test on divergence (DQAGSE label 110).
            if (!(ksgn == -1 && Math.Max(Math.Abs(result), Math.Abs(area)) <= defabs * 0.01))
            {
                if (0.01 > result / area || result / area > 100 || errsum > Math.Abs(area)) ier = 6;
            }
            return Finish(result, abserr, ier, forced, recorder);
        }

        /// <summary>
        /// Returns the composite Kronrod sum over the current partition as the result (DQAGSE label 115).
        /// </summary>
        /// <param name="error">The error estimate.</param>
        /// <param name="ier">The QUADPACK code before the final decrement.</param>
        /// <param name="forced">True when a retired interval ended the run.</param>
        /// <param name="recorder">The recorder snapshot.</param>
        /// <returns>The outcome.</returns>
        private Outcome Composite(double error, int ier, bool forced, Action<double, double, double>? recorder)
        {
            double sum = 0;
            for (int k = 0; k < _count; k++) sum += AreaOf(k);
            return Finish(sum, error, ier, forced, recorder);
        }

        /// <summary>
        /// Reports the final partition to the recorder and builds the outcome (DQAGSE label 130).
        /// </summary>
        /// <param name="result">The integral estimate.</param>
        /// <param name="error">The error estimate.</param>
        /// <param name="ier">The QUADPACK code before the final decrement.</param>
        /// <param name="forced">True when a retired interval ended the run.</param>
        /// <param name="recorder">The recorder snapshot.</param>
        /// <returns>The outcome.</returns>
        private Outcome Finish(double result, double error, int ier, bool forced, Action<double, double, double>? recorder)
        {
            if (ier > 2) ier--;
            if (recorder != null) ReportPartition(recorder);
            return new Outcome
            {
                Result = result,
                Error = error,
                Intervals = _count,
                LimitReached = ier == 1,
                NotConverged = forced || ier > 1
            };
        }

        /// <summary>
        /// Evaluates the 21-point Gauss-Kronrod rule on [lower, upper] (QUADPACK DQK21).
        /// </summary>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="offset">The scratch offset receiving the 21 node values in node order.</param>
        /// <returns>The Kronrod estimate, error estimate, and the integrals of |f| and |f - mean|.</returns>
        /// <remarks>
        /// The node values are kept in the scratch for the pass over |f - mean| and for the recorder.
        /// </remarks>
        private Estimate Kronrod21(double lower, double upper, int offset)
        {
            double halfLength = 0.5 * (upper - lower);
            var sums = EvaluateNodes(0.5 * (lower + upper), halfLength, offset);
            return RuleEstimate(halfLength, offset, sums);
        }

        /// <summary>
        /// Evaluates the 21-point Gauss-Kronrod rule on both halves of [lower, upper], the left half first.
        /// </summary>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="left">The estimate on [lower, mid], with its node values in scratch slots 0 to 20.</param>
        /// <param name="right">The estimate on [mid, upper], with its node values in scratch slots 21 to 41.</param>
        /// <remarks>
        /// Both halves are evaluated before either pass over |f - mean|, so the two passes and error estimates
        /// are independent and the processor overlaps them. The estimates equal those of two calls of
        /// <see cref="Kronrod21"/>.
        /// </remarks>
        private void Bisect(double lower, double upper, out Estimate left, out Estimate right)
        {
            double mid = 0.5 * (lower + upper);
            double leftHalf = 0.5 * (mid - lower), rightHalf = 0.5 * (upper - mid);
            var leftSums = EvaluateNodes(0.5 * (lower + mid), leftHalf, 0);
            var rightSums = EvaluateNodes(0.5 * (mid + upper), rightHalf, 21);
            left = RuleEstimate(leftHalf, 0, leftSums);
            right = RuleEstimate(rightHalf, 21, rightSums);
        }

        /// <summary>
        /// Evaluates the integrand at the 21 Kronrod nodes in the order of DQK21, stores the values in the
        /// scratch in node order, and accumulates the Gauss and Kronrod sums and the Kronrod sum of |f|.
        /// </summary>
        /// <param name="center">The interval midpoint.</param>
        /// <param name="halfLength">The interval half-length.</param>
        /// <param name="offset">The scratch offset receiving the 21 node values.</param>
        /// <returns>The unscaled Gauss, Kronrod and |f| sums.</returns>
        /// <remarks>
        /// DQK21 evaluates the center, then the pairs at the five Gauss nodes from the outermost inward, then the
        /// pairs at the five Kronrod-only nodes from the outermost inward, and accumulates each sum in that order.
        /// </remarks>
        private RuleSums EvaluateNodes(double center, double halfLength, int offset)
        {
            double[] fv = _scratch;
            double[] xk = xKronrod, wk = wKronrod, wg = wGauss;
            Func<double, double> f = Function;
            double fc = f(center);
            FunctionEvaluations++;
            fv[offset] = fc;
            double resultGauss = 0;
            double resultKronrod = wk[10] * fc;
            double resultAbs = Math.Abs(resultKronrod);
            for (int i = 1; i < 10; i += 2)
            {
                double abscissa = halfLength * xk[i];
                double fval1 = f(center - abscissa);
                double fval2 = f(center + abscissa);
                FunctionEvaluations += 2;
                fv[offset + 2 * i + 1] = fval1;
                fv[offset + 2 * i + 2] = fval2;
                double fsum = fval1 + fval2;
                resultGauss += wg[i >> 1] * fsum;
                resultKronrod += wk[i] * fsum;
                resultAbs += wk[i] * (Math.Abs(fval1) + Math.Abs(fval2));
            }
            for (int i = 0; i < 10; i += 2)
            {
                double abscissa = halfLength * xk[i];
                double fval1 = f(center - abscissa);
                double fval2 = f(center + abscissa);
                FunctionEvaluations += 2;
                fv[offset + 2 * i + 1] = fval1;
                fv[offset + 2 * i + 2] = fval2;
                resultKronrod += wk[i] * (fval1 + fval2);
                resultAbs += wk[i] * (Math.Abs(fval1) + Math.Abs(fval2));
            }
            return new RuleSums { Gauss = resultGauss, Kronrod = resultKronrod, Abs = resultAbs };
        }

        /// <summary>
        /// Completes DQK21 from the sums of <see cref="EvaluateNodes"/>: the integral of |f - mean| and the
        /// error estimate.
        /// </summary>
        /// <param name="halfLength">The interval half-length.</param>
        /// <param name="offset">The scratch offset holding the node values.</param>
        /// <param name="sums">The unscaled Gauss, Kronrod and |f| sums.</param>
        /// <returns>The Kronrod estimate, error estimate, and the integrals of |f| and |f - mean|.</returns>
        /// <remarks>
        /// The sum of |f - mean| runs over the center and then the ten pairs from the outermost inward, as in
        /// DQK21. (200 r)^1.5 is evaluated as q sqrt(q) with q = 200 r.
        /// </remarks>
        private Estimate RuleEstimate(double halfLength, int offset, RuleSums sums)
        {
            double[] fv = _scratch;
            double[] wk = wKronrod;
            double resultKronrod = sums.Kronrod;
            double mean = resultKronrod * 0.5;
            double resultAsc = wk[10] * Math.Abs(fv[offset] - mean);
            for (int i = 0; i < 10; i++)
                resultAsc += wk[i] * (Math.Abs(fv[offset + 2 * i + 1] - mean) + Math.Abs(fv[offset + 2 * i + 2] - mean));

            double resultGauss = sums.Gauss, resultAbs = sums.Abs;

            double absHalfLength = Math.Abs(halfLength);
            var estimate = new Estimate
            {
                Area = resultKronrod * halfLength,
                AbsArea = resultAbs * absHalfLength,
                AscArea = resultAsc * absHalfLength
            };
            double error = Math.Abs((resultKronrod - resultGauss) * halfLength);
            if (estimate.AscArea != 0 && error != 0)
            {
                double q = 200 * error / estimate.AscArea;
                error = estimate.AscArea * (q < 1 ? q * Math.Sqrt(q) : 1);
            }
            if (estimate.AbsArea > Uflow / (50 * Epmach)) error = Math.Max(50 * Epmach * estimate.AbsArea, error);
            estimate.Error = error;
            return estimate;
        }

        /// <summary>
        /// Evaluates the uniform initial partition of [lower, upper] to the given depth, left to right.
        /// </summary>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="level">The bisection level of [lower, upper].</param>
        /// <param name="depth">The partition depth.</param>
        /// <param name="limit">The maximum number of subintervals.</param>
        /// <param name="record">True when node values are kept for the recorder.</param>
        /// <param name="sumArea">Accumulates the estimates.</param>
        /// <param name="sumError">Accumulates the error estimates.</param>
        /// <param name="sumAbs">Accumulates the integrals of |f|.</param>
        /// <param name="sumAsc">Accumulates the integrals of |f - mean|.</param>
        private void Partition(double lower, double upper, int level, int depth, int limit, bool record,
            ref double sumArea, ref double sumError, ref double sumAbs, ref double sumAsc)
        {
            if (level < depth && CanBisect(lower, upper, level))
            {
                double mid = 0.5 * (lower + upper);
                Partition(lower, mid, level + 1, depth, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
                Partition(mid, upper, level + 1, depth, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
                return;
            }
            var estimate = Kronrod21(lower, upper, 0);
            int slot = Allocate();
            if (estimate.IsFinite)
            {
                PlaceFinite(slot, lower, upper, level, estimate, 0, record, ref sumArea, ref sumError);
                sumAbs += estimate.AbsArea;
                sumAsc += estimate.AscArea;
            }
            else
            {
                Store(slot, lower, upper, level, estimate, 0, record);
                Resolve(slot, lower, upper, level, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
            }
        }

        /// <summary>
        /// Bisects an interval with a non-finite estimate depth first until every piece is finite, storing the
        /// pieces in the queue. Both halves are evaluated and stored before either is resolved further. When the
        /// subinterval limit is reached first, the non-finite interval stays in the partition outside the queue,
        /// its estimate enters the sums, and the run stops with the limit status.
        /// </summary>
        /// <param name="slot">The slot of the interval, already stored; it receives the left half.</param>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="level">The bisection level of [lower, upper].</param>
        /// <param name="limit">The maximum number of subintervals.</param>
        /// <param name="record">True when node values are kept for the recorder.</param>
        /// <param name="sumArea">Accumulates the estimates of the finite pieces.</param>
        /// <param name="sumError">Accumulates their error estimates.</param>
        /// <param name="sumAbs">Accumulates their integrals of |f|.</param>
        /// <param name="sumAsc">Accumulates their integrals of |f - mean|.</param>
        /// <exception cref="ArithmeticException">Thrown when the interval is at the depth or width limit.</exception>
        private void Resolve(int slot, double lower, double upper, int level, int limit, bool record,
            ref double sumArea, ref double sumError, ref double sumAbs, ref double sumAsc)
        {
            if (_budgetExhausted) return;
            if (!CanBisect(lower, upper, level))
                throw new ArithmeticException("The integrand produced a NaN or infinite estimate on an interval that cannot be bisected further, so the integral cannot be finite.");

            if (_count >= limit)
            {
                _budgetExhausted = true;
                sumArea += AreaOf(slot);
                sumError += ErrorOf(slot);
                return;
            }

            double mid = 0.5 * (lower + upper);
            Bisect(lower, upper, out var left, out var right);
            Iterations++;
            int other = Allocate();
            if (left.IsFinite)
            {
                PlaceFinite(slot, lower, mid, level + 1, left, 0, record, ref sumArea, ref sumError);
                sumAbs += left.AbsArea;
                sumAsc += left.AscArea;
            }
            if (right.IsFinite)
            {
                PlaceFinite(other, mid, upper, level + 1, right, 21, record, ref sumArea, ref sumError);
                sumAbs += right.AbsArea;
                sumAsc += right.AscArea;
            }
            if (!left.IsFinite) Store(slot, lower, mid, level + 1, left, 0, record);
            if (!right.IsFinite) Store(other, mid, upper, level + 1, right, 21, record);
            if (!left.IsFinite) Resolve(slot, lower, mid, level + 1, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
            if (!right.IsFinite) Resolve(other, mid, upper, level + 1, limit, record, ref sumArea, ref sumError, ref sumAbs, ref sumAsc);
        }

        /// <summary>
        /// Stores a finite interval, queues it and adds its estimate and error to the running sums.
        /// </summary>
        /// <param name="slot">The interval slot.</param>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="level">The bisection level.</param>
        /// <param name="estimate">The interval estimate.</param>
        /// <param name="offset">The scratch offset holding its node values.</param>
        /// <param name="record">True when node values are kept for the recorder.</param>
        /// <param name="sumArea">Accumulates the estimate.</param>
        /// <param name="sumError">Accumulates the error estimate.</param>
        private void PlaceFinite(int slot, double lower, double upper, int level, Estimate estimate, int offset, bool record,
            ref double sumArea, ref double sumError)
        {
            Store(slot, lower, upper, level, estimate, offset, record);
            Push(slot);
            sumArea += estimate.Area;
            sumError += estimate.Error;
        }

        /// <summary>
        /// Returns true when the interval in a slot may be bisected.
        /// </summary>
        /// <param name="slot">The interval slot.</param>
        /// <returns>True when its level is below <see cref="MaxDepth"/> and its width passes the DQAGSE resolution test.</returns>
        private bool CanBisect(int slot)
        {
            return CanBisect(Lower(slot), Upper(slot), Level(slot));
        }

        /// <summary>
        /// Returns true when [lower, upper] may be bisected: its level is below <see cref="MaxDepth"/> and it
        /// passes the DQAGSE resolution test max(|a1|, |b2|) &gt; (1 + 100 eps)(|a2| + 1000 uflow) on its halves.
        /// </summary>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="level">The bisection level.</param>
        /// <returns>True when the interval may be bisected.</returns>
        private bool CanBisect(double lower, double upper, int level)
        {
            if (level >= MaxDepth) return false;
            double mid = 0.5 * (lower + upper);
            return Math.Max(Math.Abs(lower), Math.Abs(upper)) > (1 + 100 * Epmach) * (Math.Abs(mid) + 1000 * Uflow);
        }

        /// <summary>
        /// Returns the width of the interval in a slot.
        /// </summary>
        /// <param name="slot">The interval slot.</param>
        /// <returns>|upper - lower|.</returns>
        private double Width(int slot)
        {
            return Math.Abs(Upper(slot) - Lower(slot));
        }

        /// <summary>
        /// Reserves the next interval slot, growing the interval storage as needed.
        /// </summary>
        /// <returns>The slot index.</returns>
        private int Allocate()
        {
            if (_count == _capacity)
            {
                _capacity = Math.Max(16, 2 * _capacity);
                Array.Resize(ref _intervals, 4 * _capacity);
                Array.Resize(ref _slots, 2 * _capacity);
                Array.Resize(ref _queue, _capacity);
                Array.Resize(ref _skipped, _capacity);
            }
            return _count++;
        }

        /// <summary>
        /// Writes an interval into a slot and, when recording, copies its node values from the scratch.
        /// </summary>
        /// <param name="slot">The interval slot.</param>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="level">The bisection level.</param>
        /// <param name="estimate">The interval estimate.</param>
        /// <param name="offset">The scratch offset holding its node values.</param>
        /// <param name="record">True when node values are kept for the recorder.</param>
        private void Store(int slot, double lower, double upper, int level, Estimate estimate, int offset, bool record)
        {
            Lower(slot) = lower;
            Upper(slot) = upper;
            AreaOf(slot) = estimate.Area;
            ErrorOf(slot) = estimate.Error;
            Level(slot) = level;
            Sequence(slot) = ++_sequenceCounter;
            if (record)
            {
                if (_values.Length < 21 * _capacity) Array.Resize(ref _values, 21 * _capacity);
                Array.Copy(_scratch, offset, _values, 21 * slot, 21);
            }
        }

        /// <summary>
        /// Returns true when entry a precedes entry b in the queue: a larger error, or an equal error and a
        /// later store.
        /// </summary>
        /// <param name="a">The first entry.</param>
        /// <param name="b">The second entry.</param>
        /// <returns>True when entry a has priority.</returns>
        private static bool Precedes(in QueueEntry a, in QueueEntry b)
        {
            return a.Error > b.Error || (a.Error == b.Error && a.Sequence > b.Sequence);
        }

        /// <summary>
        /// Adds a slot to the queue.
        /// </summary>
        /// <param name="slot">The interval slot.</param>
        private void Push(int slot)
        {
            Push(new QueueEntry { Error = ErrorOf(slot), Sequence = Sequence(slot), Slot = slot });
        }

        /// <summary>
        /// Adds an entry to the queue.
        /// </summary>
        /// <param name="entry">The entry.</param>
        private void Push(QueueEntry entry)
        {
            var q = _queue;
            int k = _heapCount++;
            while (k > 0)
            {
                int parent = (k - 1) >> 1;
                if (!Precedes(entry, q[parent])) break;
                q[k] = q[parent];
                k = parent;
            }
            q[k] = entry;
        }

        /// <summary>
        /// Removes and returns the entry with the largest error estimate.
        /// </summary>
        /// <returns>The entry.</returns>
        private QueueEntry PopEntry()
        {
            var q = _queue;
            var top = q[0];
            var moved = q[--_heapCount];
            int k = 0;
            while (true)
            {
                int child = 2 * k + 1;
                if (child >= _heapCount) break;
                if (child + 1 < _heapCount && Precedes(q[child + 1], q[child])) child++;
                if (!Precedes(q[child], moved)) break;
                q[k] = q[child];
                k = child;
            }
            if (_heapCount > 0) q[k] = moved;
            return top;
        }

        /// <summary>
        /// Removes and returns the slot with the largest error estimate.
        /// </summary>
        /// <returns>The interval slot.</returns>
        private int PopTop()
        {
            return PopEntry().Slot;
        }

        /// <summary>
        /// Sets the top entry aside while the extrapolation step passes over small intervals.
        /// </summary>
        private void SkipTop()
        {
            var entry = PopEntry();
            _skipped[_skippedCount++] = entry;
        }

        /// <summary>
        /// Returns the intervals set aside by the extrapolation step to the queue.
        /// </summary>
        private void RestoreSkipped()
        {
            for (int i = 0; i < _skippedCount; i++) Push(_skipped[i]);
            _skippedCount = 0;
        }

        /// <summary>
        /// Applies Wynn's epsilon algorithm to the sequence of composite sums (QUADPACK DQELG).
        /// </summary>
        /// <param name="n">The number of elements in the epsilon table, at least three; adjusted on return.</param>
        /// <param name="result">The extrapolated value.</param>
        /// <param name="abserr">The error estimate of the extrapolated value.</param>
        /// <param name="nres">The number of calls so far; incremented.</param>
        /// <remarks>
        /// The table holds the newest diagonal of the epsilon scheme with at most 50 elements. The error
        /// estimate compares the result with the last three results once four are available, and is never
        /// below 5 eps |result|.
        /// </remarks>
        private void Extrapolate(ref int n, out double result, out double abserr, ref int nres)
        {
            double[] epstab = _epsilonTable;
            double[] res3la = _lastResults;
            nres++;
            abserr = Oflow;
            result = epstab[n - 1];

            epstab[n + 1] = epstab[n - 1];
            int newelm = (n - 1) / 2;
            epstab[n - 1] = Oflow;
            int num = n;
            int k1 = n;
            for (int i = 1; i <= newelm; i++)
            {
                int k2 = k1 - 1;
                int k3 = k1 - 2;
                double res = epstab[k1 + 1];
                double e0 = epstab[k3 - 1];
                double e1 = epstab[k2 - 1];
                double e2 = res;
                double e1abs = Math.Abs(e1);
                double delta2 = e2 - e1;
                double err2 = Math.Abs(delta2);
                double tol2 = Math.Max(Math.Abs(e2), e1abs) * Epmach;
                double delta3 = e1 - e0;
                double err3 = Math.Abs(delta3);
                double tol3 = Math.Max(e1abs, Math.Abs(e0)) * Epmach;
                if (err2 <= tol2 && err3 <= tol3)
                {
                    // e0, e1 and e2 are equal to within machine accuracy; convergence is assumed.
                    result = res;
                    abserr = Math.Max(err2 + err3, 5 * Epmach * Math.Abs(result));
                    return;
                }

                double e3 = epstab[k1 - 1];
                epstab[k1 - 1] = e1;
                double delta1 = e1 - e3;
                double err1 = Math.Abs(delta1);
                double tol1 = Math.Max(e1abs, Math.Abs(e3)) * Epmach;

                // If two elements are very close to each other, or the table behaves irregularly, omit a part
                // of the table by adjusting the value of n.
                if (err1 <= tol1 || err2 <= tol2 || err3 <= tol3)
                {
                    n = i + i - 1;
                    break;
                }
                double ss = 1 / delta1 + 1 / delta2 - 1 / delta3;
                double epsinf = Math.Abs(ss * e1);
                if (epsinf <= 1E-4)
                {
                    n = i + i - 1;
                    break;
                }

                // Compute a new element and eventually adjust the value of result.
                res = e1 + 1 / ss;
                epstab[k1 - 1] = res;
                k1 -= 2;
                double error = err2 + Math.Abs(res - e2) + err3;
                if (error <= abserr)
                {
                    abserr = error;
                    result = res;
                }
            }

            // Shift the table.
            if (n == EpsilonLimit) n = 2 * (EpsilonLimit / 2) - 1;
            int ib = num % 2 == 0 ? 2 : 1;
            int ie = newelm + 1;
            for (int i = 1; i <= ie; i++)
            {
                int ib2 = ib + 2;
                epstab[ib - 1] = epstab[ib2 - 1];
                ib = ib2;
            }
            if (num != n)
            {
                int indx = num - n + 1;
                for (int i = 1; i <= n; i++)
                {
                    epstab[i - 1] = epstab[indx - 1];
                    indx++;
                }
            }

            if (nres < 4)
            {
                res3la[nres - 1] = result;
                abserr = Oflow;
            }
            else
            {
                abserr = Math.Abs(result - res3la[2]) + Math.Abs(result - res3la[1]) + Math.Abs(result - res3la[0]);
                res3la[0] = res3la[1];
                res3la[1] = res3la[2];
                res3la[2] = result;
            }
            abserr = Math.Max(abserr, 5 * Epmach * Math.Abs(result));
        }

        /// <summary>
        /// Reports the nodes of every interval of the final partition, intervals ordered by lower limit.
        /// </summary>
        /// <param name="recorder">The recorder snapshot.</param>
        private void ReportPartition(Action<double, double, double> recorder)
        {
            var order = new int[_count];
            var keys = new double[_count];
            for (int k = 0; k < _count; k++)
            {
                order[k] = k;
                keys[k] = Lower(k);
            }
            Array.Sort(keys, order);
            for (int k = 0; k < _count; k++)
            {
                int slot = order[k];
                ReportInterval(recorder, Lower(slot), Upper(slot), _values, 21 * slot);
            }
        }

        /// <summary>
        /// Reports the 21 nodes of one interval with weights scaled by its half-length.
        /// </summary>
        /// <param name="recorder">The recorder snapshot.</param>
        /// <param name="lower">The lower limit.</param>
        /// <param name="upper">The upper limit.</param>
        /// <param name="values">The node values in node order.</param>
        /// <param name="offset">The offset of the first value.</param>
        /// <remarks>
        /// The abscissas are recomputed with the same arithmetic as <see cref="Kronrod21"/>, so they match the
        /// evaluated points exactly.
        /// </remarks>
        private static void ReportInterval(Action<double, double, double> recorder, double lower, double upper, double[] values, int offset)
        {
            double center = 0.5 * (lower + upper);
            double halfLength = 0.5 * (upper - lower);
            recorder(center, wCapture[0] * halfLength, values[offset]);
            for (int i = 0; i < 10; i++)
            {
                double abscissa = halfLength * xKronrod[i];
                recorder(center - abscissa, wCapture[2 * i + 1] * halfLength, values[offset + 2 * i + 1]);
                recorder(center + abscissa, wCapture[2 * i + 2] * halfLength, values[offset + 2 * i + 2]);
            }
        }
    }
}
