namespace Numerics.Mathematics.Integration
{
    /// <summary>
    /// A class that performs globally adaptive Gauss-Kronrod integration over a two-dimensional rectangle.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Each rectangular region is evaluated with the tensor product of the 10-point Gauss, 21-point
    /// Kronrod rule (G10K21) on both axes: the 21 x 21 Kronrod product supplies the region estimate,
    /// and the embedded 10 x 10 Gauss product - read from the same 441 function values - supplies the
    /// error estimate as the absolute difference between the two. Refinement is global rather than
    /// recursive: live regions are kept on a max-heap keyed by error, the worst region is bisected at
    /// its midpoint along the axis with the larger one-axis error indicator (the mixed Gauss-Kronrod
    /// products, also read from the same function values), and integration stops when the sum of
    /// region errors meets the absolute or relative tolerance, or a budget is exhausted. All product
    /// Kronrod weights are positive, so the accepted regions form a composite rule whose node weights
    /// partition the domain area exactly; the optional <see cref="Recorder"/> exposes that composite
    /// rule for callers that account for integration mass.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><see href="https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula"/></item>
    /// <item>Piessens, R., et al. (1983). QUADPACK: A Subroutine Package for Automatic Integration. Springer-Verlag.
    /// (The global error-directed refinement strategy follows the QAG pattern.)</item>
    /// </list>
    /// </remarks>
    [Serializable]
    public class AdaptiveGaussKronrod2D : Integrator
    {
        /// <summary>
        /// Constructs a new globally adaptive two-dimensional Gauss-Kronrod rule.
        /// </summary>
        /// <param name="function">The two-dimensional function to integrate.</param>
        /// <param name="minX">The minimum x-value under which the integral must be computed.</param>
        /// <param name="maxX">The maximum x-value under which the integral must be computed.</param>
        /// <param name="minY">The minimum y-value under which the integral must be computed.</param>
        /// <param name="maxY">The maximum y-value under which the integral must be computed.</param>
        /// <exception cref="ArgumentNullException">Thrown when the function is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when a bound is not finite, or when a maximum bound is less than or equal to its minimum.
        /// </exception>
        public AdaptiveGaussKronrod2D(Func<double, double, double> function, double minX, double maxX, double minY, double maxY)
        {
            Function = function ?? throw new ArgumentNullException(nameof(function), "The function cannot be null.");
            if (!Tools.IsFinite(minX)) throw new ArgumentOutOfRangeException(nameof(minX), "The minimum x-value must be finite.");
            if (!Tools.IsFinite(maxX)) throw new ArgumentOutOfRangeException(nameof(maxX), "The maximum x-value must be finite.");
            if (!Tools.IsFinite(minY)) throw new ArgumentOutOfRangeException(nameof(minY), "The minimum y-value must be finite.");
            if (!Tools.IsFinite(maxY)) throw new ArgumentOutOfRangeException(nameof(maxY), "The maximum y-value must be finite.");
            if (maxX <= minX) throw new ArgumentOutOfRangeException(nameof(maxX), "The maximum x-value cannot be less than or equal to the minimum x-value.");
            if (maxY <= minY) throw new ArgumentOutOfRangeException(nameof(maxY), "The maximum y-value cannot be less than or equal to the minimum y-value.");
            ax = minX;
            bx = maxX;
            ay = minY;
            by = maxY;
        }

        private double ax, bx, ay, by;

        // Gauss-Kronrod G10K21 positive abscissas and weights on [-1, 1], identical constants to the
        // one-dimensional AdaptiveGaussKronrod rule. Kronrod indices 1, 3, 5, 7, 9 are the Gauss nodes.
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

        private static readonly double[] wGauss = new double[]
        {
            0.066671344308688137593568809893332,
            0.149451349150580593145776339657697,
            0.219086362515982043995534934228163,
            0.269266719309996355091226921569469,
            0.295524224714752870173892994651338
        };

        /// <summary>
        /// The 21 signed node abscissas on [-1, 1] in ascending order.
        /// </summary>
        private static readonly double[] Nodes = BuildNodes();

        /// <summary>
        /// The Kronrod weights aligned with <see cref="Nodes"/>.
        /// </summary>
        private static readonly double[] WeightsKronrod = BuildWeights(false);

        /// <summary>
        /// The Gauss weights aligned with <see cref="Nodes"/>; zero at non-Gauss nodes.
        /// </summary>
        private static readonly double[] WeightsGauss = BuildWeights(true);

        /// <summary>
        /// Builds the ascending signed node layout from the positive-abscissa table.
        /// </summary>
        /// <returns>The 21 signed abscissas in ascending order.</returns>
        private static double[] BuildNodes()
        {
            var nodes = new double[21];
            for (int i = 0; i < 10; i++)
            {
                nodes[i] = -xKronrod[i];
                nodes[20 - i] = xKronrod[i];
            }
            nodes[10] = 0d;
            return nodes;
        }

        /// <summary>
        /// Builds a weight vector aligned with the ascending node layout.
        /// </summary>
        /// <param name="gauss">True for the embedded Gauss weights (zero off the Gauss nodes); false for the Kronrod weights.</param>
        /// <returns>The 21 aligned weights.</returns>
        private static double[] BuildWeights(bool gauss)
        {
            var weights = new double[21];
            for (int p = 0; p < 21; p++)
            {
                int i = p <= 10 ? p : 20 - p;
                if (gauss)
                {
                    // Gauss membership: odd indices of the positive-abscissa Kronrod table.
                    weights[p] = i % 2 == 1 ? wGauss[i / 2] : 0d;
                }
                else
                {
                    weights[p] = wKronrod[i];
                }
            }
            return weights;
        }

        /// <summary>
        /// A live or frozen rectangular region of the composite rule.
        /// </summary>
        [Serializable]
        private sealed class Region
        {
            /// <summary>The region bounds.</summary>
            public double Ax, Bx, Ay, By;

            /// <summary>The tensor Kronrod estimate over the region.</summary>
            public double Kronrod;

            /// <summary>The region error estimate |K - G|.</summary>
            public double Error;

            /// <summary>The x-axis error indicator |K - (Gx tensor Ky)|.</summary>
            public double ErrorX;

            /// <summary>The y-axis error indicator |K - (Kx tensor Gy)|.</summary>
            public double ErrorY;

            /// <summary>The bisection depth of the region (the root is zero).</summary>
            public int Depth;

            /// <summary>True once the region has been replaced by its children.</summary>
            public bool Split;

            /// <summary>True once the region can no longer be refined (depth or width exhausted).</summary>
            public bool Frozen;

            /// <summary>The 21 x-node abscissas; allocated only when a recorder is attached.</summary>
            public double[]? XNodes;

            /// <summary>The 21 y-node abscissas; allocated only when a recorder is attached.</summary>
            public double[]? YNodes;

            /// <summary>The 441 function values in row-major (x, y) order; allocated only when a recorder is attached.</summary>
            public double[]? Values;
        }

        /// <summary>
        /// The two-dimensional function to integrate.
        /// </summary>
        public Func<double, double, double> Function { get; }

        /// <summary>
        /// The minimum x-value under which the integral must be computed.
        /// </summary>
        public double MinX => ax;

        /// <summary>
        /// The maximum x-value under which the integral must be computed.
        /// </summary>
        public double MaxX => bx;

        /// <summary>
        /// The minimum y-value under which the integral must be computed.
        /// </summary>
        public double MinY => ay;

        /// <summary>
        /// The maximum y-value under which the integral must be computed.
        /// </summary>
        public double MaxY => by;

        /// <summary>
        /// The minimum bisection depth. Every region is refined to at least this depth before
        /// convergence may be declared. Default = 0.
        /// </summary>
        public int MinDepth { get; set; } = 0;

        /// <summary>
        /// The maximum bisection depth of any region. A region at this depth is accepted with its
        /// current error estimate. Default = 100.
        /// </summary>
        public int MaxDepth { get; set; } = 100;

        /// <summary>
        /// Returns the global error bound of the integration: the sum of the absolute Kronrod-Gauss
        /// differences over the final regions.
        /// </summary>
        /// <remarks>
        /// This is the quantity tested against the tolerances, following the QUADPACK global
        /// strategy. It is a conservative bound rather than the root-sum-square the one-dimensional
        /// rule reports, so equal names do not imply equal scales across the two classes.
        /// </remarks>
        public double StandardError { get; private set; }

        /// <summary>
        /// Gets or sets a callback invoked as (x, y, weight, f(x, y)) for the quadrature nodes of the
        /// final composite rule. A null callback disables recording.
        /// </summary>
        /// <remarks>
        /// The callback is snapshotted when integration begins and flushed once at completion over the
        /// regions that make up the final composite rule; nodes of subdivided parent regions are never
        /// reported. Each weight carries the tensor Kronrod weight scaled by the region half-lengths,
        /// so the reported weights sum to the domain area and the weighted function values reproduce
        /// <see cref="Integrator.Result"/> up to floating-point reassociation. Nothing is reported
        /// when integration fails. Node capture allocates three arrays per region, so leave the
        /// recorder null when the composite rule is not needed.
        /// </remarks>
        public Action<double, double, double, double>? Recorder { get; set; }

        /// <inheritdoc/>
        /// <remarks>
        /// <see cref="Integrator.Status"/> reports <see cref="IntegrationStatus.Success"/> when the
        /// error bound meets a tolerance, <see cref="IntegrationStatus.MaximumFunctionEvaluationsReached"/>
        /// when the evaluation budget stops refinement first, and
        /// <see cref="IntegrationStatus.MaximumIterationsReached"/> when every region has reached
        /// <see cref="MaxDepth"/> (or a width bisection can no longer reduce) without meeting a
        /// tolerance - the depth-exhausted result is returned rather than silently reported as
        /// converged. Refinement is bounded by <see cref="Integrator.MaxFunctionEvaluations"/> and
        /// <see cref="MaxDepth"/>; the inherited <see cref="Integrator.MinIterations"/> and
        /// <see cref="Integrator.MaxIterations"/> are not consulted, and
        /// <see cref="Integrator.Iterations"/> reports the number of region splits. The computation
        /// is sequential and deterministic: identical inputs produce bit-identical results.
        /// </remarks>
        public override void Integrate()
        {
            StandardError = 0;
            ClearResults();
            Validate();
            if (MinDepth < 0) throw new ArgumentOutOfRangeException(nameof(MinDepth), "The minimum depth cannot be negative.");
            if (MaxDepth < MinDepth) throw new ArgumentOutOfRangeException(nameof(MaxDepth), "The maximum depth cannot be less than the minimum depth.");
            Action<double, double, double, double>? recorder = Recorder;
            bool capture = recorder != null;

            try
            {
                // All regions ever created, in creation order. Split parents remain in the list but are
                // excluded from the final sums; final sums iterate this list so the floating-point
                // association order is independent of the heap's internal layout.
                var regions = new List<Region> { Evaluate(ax, bx, ay, by, 0, capture) };

                // Max-heap of refinable region indices keyed by error estimate.
                var heap = new List<(int Index, double Error)>();
                Push(heap, (0, regions[0].Error));

                // Running sums drive the convergence test; the reported result and error bound are
                // re-summed over the final regions below.
                double resultSum = regions[0].Kronrod;
                double errorSum = regions[0].Error;

                // Refine every region to the minimum depth before convergence may be declared.
                for (int i = 0; i < regions.Count; i++)
                {
                    var region = regions[i];
                    if (region.Split || region.Frozen || region.Depth >= MinDepth) continue;
                    if (FunctionEvaluations >= MaxFunctionEvaluations) break;
                    SplitRegion(regions, heap, region, capture, ref resultSum, ref errorSum);
                }

                var status = IntegrationStatus.Success;
                while (true)
                {
                    double tolerance = Math.Max(AbsoluteTolerance, RelativeTolerance * Math.Abs(resultSum));
                    if (FunctionEvaluations >= MinFunctionEvaluations && errorSum <= tolerance)
                    {
                        status = IntegrationStatus.Success;
                        break;
                    }
                    if (FunctionEvaluations >= MaxFunctionEvaluations)
                    {
                        status = IntegrationStatus.MaximumFunctionEvaluationsReached;
                        break;
                    }
                    if (heap.Count == 0)
                    {
                        // Every region is depth- or width-exhausted and the tolerance is unmet.
                        status = IntegrationStatus.MaximumIterationsReached;
                        break;
                    }

                    var (index, _) = Pop(heap);
                    var worst = regions[index];
                    // Lazy deletion: entries whose region was split during the minimum-depth pass (or
                    // frozen) are stale - the children carry their own entries.
                    if (worst.Split || worst.Frozen) continue;
                    if (worst.Depth >= MaxDepth)
                    {
                        worst.Frozen = true;
                        continue;
                    }
                    SplitRegion(regions, heap, worst, capture, ref resultSum, ref errorSum);
                }

                // Re-sum the final composite rule in creation order.
                double result = 0, error = 0;
                for (int i = 0; i < regions.Count; i++)
                {
                    var region = regions[i];
                    if (region.Split) continue;
                    result += region.Kronrod;
                    error += region.Error;
                }
                Result = result;
                StandardError = error;
                Status = status;

                if (recorder != null)
                {
                    FlushRecorder(regions, recorder);
                }
            }
            catch (Exception)
            {
                Status = IntegrationStatus.Failure;
                if (ReportFailure) throw;
            }
        }

        /// <summary>
        /// Bisects a region along its dominant-error axis, replacing it with two evaluated children.
        /// A region is frozen instead when both axes are at machine-epsilon width, or when the chosen
        /// axis can no longer be bisected in floating point because its midpoint rounds onto an
        /// endpoint.
        /// </summary>
        /// <param name="regions">The region list; the children are appended.</param>
        /// <param name="heap">The refinable-region heap; the children are pushed.</param>
        /// <param name="region">The region to bisect.</param>
        /// <param name="capture">True to capture nodes and values for the recorder.</param>
        /// <param name="resultSum">The running result sum, updated in place.</param>
        /// <param name="errorSum">The running error sum, updated in place.</param>
        private void SplitRegion(List<Region> regions, List<(int Index, double Error)> heap, Region region,
            bool capture, ref double resultSum, ref double errorSum)
        {
            // Prefer the axis with the larger one-axis error indicator; fall back to the other when
            // the preferred axis has collapsed to machine width.
            bool splitX = region.ErrorX >= region.ErrorY;
            bool xTooNarrow = Math.Abs(region.Bx - region.Ax) <= Tools.DoubleMachineEpsilon;
            bool yTooNarrow = Math.Abs(region.By - region.Ay) <= Tools.DoubleMachineEpsilon;
            if (splitX && xTooNarrow) splitX = false;
            if (!splitX && yTooNarrow)
            {
                if (xTooNarrow)
                {
                    region.Frozen = true;
                    return;
                }
                splitX = true;
            }

            // Far from the origin an axis can be wider than the absolute floor yet sit at one unit
            // in the last place, where its midpoint rounds onto an endpoint and bisection would
            // reproduce the region bit for bit - re-splitting an identical child once per depth
            // level, at a full tensor evaluation per wasted split, until the requested depth is
            // reached. The chosen axis carries the dominant error indicator, so when it can no
            // longer be bisected the region's error is irreducible at floating-point resolution and
            // the region is frozen with its current estimate; splitting the other axis instead could
            // not reduce the dominant component and would multiply regions without converging.
            Region left, right;
            if (splitX)
            {
                double mx = 0.5 * (region.Ax + region.Bx);
                if (mx <= region.Ax || mx >= region.Bx)
                {
                    region.Frozen = true;
                    return;
                }
                left = Evaluate(region.Ax, mx, region.Ay, region.By, region.Depth + 1, capture);
                right = Evaluate(mx, region.Bx, region.Ay, region.By, region.Depth + 1, capture);
            }
            else
            {
                double my = 0.5 * (region.Ay + region.By);
                if (my <= region.Ay || my >= region.By)
                {
                    region.Frozen = true;
                    return;
                }
                left = Evaluate(region.Ax, region.Bx, region.Ay, my, region.Depth + 1, capture);
                right = Evaluate(region.Ax, region.Bx, my, region.By, region.Depth + 1, capture);
            }

            region.Split = true;
            region.XNodes = null;
            region.YNodes = null;
            region.Values = null;
            resultSum += left.Kronrod + right.Kronrod - region.Kronrod;
            errorSum += left.Error + right.Error - region.Error;
            Iterations++;

            regions.Add(left);
            Push(heap, (regions.Count - 1, left.Error));
            regions.Add(right);
            Push(heap, (regions.Count - 1, right.Error));
        }

        /// <summary>
        /// Evaluates the G10K21 tensor product over a rectangle, returning the region record with its
        /// Kronrod estimate, error estimate, and the two one-axis error indicators.
        /// </summary>
        /// <param name="rax">The lower x-bound of the region.</param>
        /// <param name="rbx">The upper x-bound of the region.</param>
        /// <param name="ray">The lower y-bound of the region.</param>
        /// <param name="rby">The upper y-bound of the region.</param>
        /// <param name="depth">The bisection depth of the region.</param>
        /// <param name="capture">True to retain nodes and values for the recorder.</param>
        /// <returns>The evaluated region.</returns>
        private Region Evaluate(double rax, double rbx, double ray, double rby, int depth, bool capture)
        {
            double cx = 0.5 * (rax + rbx);
            double hx = 0.5 * (rbx - rax);
            double cy = 0.5 * (ray + rby);
            double hy = 0.5 * (rby - ray);

            var xs = new double[21];
            var ys = new double[21];
            for (int i = 0; i < 21; i++)
            {
                xs[i] = cx + hx * Nodes[i];
                ys[i] = cy + hy * Nodes[i];
            }
            double[]? values = capture ? new double[441] : null;

            // Accumulate the four weight products in one sweep: Kronrod-Kronrod (the estimate),
            // Gauss-Gauss (the error pair), and the two mixed products (the axis indicators).
            double sumKK = 0, sumGG = 0, sumGxKy = 0, sumKxGy = 0;
            for (int i = 0; i < 21; i++)
            {
                double rowK = 0, rowG = 0;
                double x = xs[i];
                for (int j = 0; j < 21; j++)
                {
                    double f = Function(x, ys[j]);
                    if (values != null) values[i * 21 + j] = f;
                    rowK += WeightsKronrod[j] * f;
                    rowG += WeightsGauss[j] * f;
                }
                sumKK += WeightsKronrod[i] * rowK;
                sumKxGy += WeightsKronrod[i] * rowG;
                sumGxKy += WeightsGauss[i] * rowK;
                sumGG += WeightsGauss[i] * rowG;
            }
            FunctionEvaluations += 441;

            double scale = hx * hy;
            double kronrod = sumKK * scale;
            var region = new Region
            {
                Ax = rax,
                Bx = rbx,
                Ay = ray,
                By = rby,
                Kronrod = kronrod,
                Error = Math.Abs(kronrod - sumGG * scale),
                ErrorX = Math.Abs(kronrod - sumGxKy * scale),
                ErrorY = Math.Abs(kronrod - sumKxGy * scale),
                Depth = depth,
                XNodes = capture ? xs : null,
                YNodes = capture ? ys : null,
                Values = values
            };
            return region;
        }

        /// <summary>
        /// Reports the final composite rule to the recorder in region-creation order.
        /// </summary>
        /// <param name="regions">The region list.</param>
        /// <param name="recorder">The recorder snapshot.</param>
        private static void FlushRecorder(List<Region> regions, Action<double, double, double, double> recorder)
        {
            for (int r = 0; r < regions.Count; r++)
            {
                var region = regions[r];
                if (region.Split || region.XNodes == null || region.YNodes == null || region.Values == null) continue;
                double scale = 0.25 * (region.Bx - region.Ax) * (region.By - region.Ay);
                for (int i = 0; i < 21; i++)
                {
                    double wx = WeightsKronrod[i];
                    for (int j = 0; j < 21; j++)
                    {
                        recorder(region.XNodes[i], region.YNodes[j], wx * WeightsKronrod[j] * scale, region.Values[i * 21 + j]);
                    }
                }
            }
        }

        /// <summary>
        /// Pushes an entry onto the max-heap.
        /// </summary>
        /// <param name="heap">The heap storage.</param>
        /// <param name="entry">The region index and its error key.</param>
        private static void Push(List<(int Index, double Error)> heap, (int Index, double Error) entry)
        {
            heap.Add(entry);
            int child = heap.Count - 1;
            while (child > 0)
            {
                int parent = (child - 1) / 2;
                if (heap[parent].Error >= heap[child].Error) break;
                (heap[parent], heap[child]) = (heap[child], heap[parent]);
                child = parent;
            }
        }

        /// <summary>
        /// Pops the maximum-error entry from the heap.
        /// </summary>
        /// <param name="heap">The heap storage.</param>
        /// <returns>The region index and its error key.</returns>
        private static (int Index, double Error) Pop(List<(int Index, double Error)> heap)
        {
            var top = heap[0];
            int last = heap.Count - 1;
            heap[0] = heap[last];
            heap.RemoveAt(last);
            int parent = 0;
            while (true)
            {
                int left = 2 * parent + 1;
                if (left >= heap.Count) break;
                int right = left + 1;
                int larger = right < heap.Count && heap[right].Error > heap[left].Error ? right : left;
                if (heap[parent].Error >= heap[larger].Error) break;
                (heap[parent], heap[larger]) = (heap[larger], heap[parent]);
                parent = larger;
            }
            return top;
        }
    }
}
