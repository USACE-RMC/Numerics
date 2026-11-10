using System;
using System.Collections.Generic;

namespace Numerics.Mathematics.Optimization
{

    /// <summary>
    /// Pareto-front utilities for multi-objective optimization: dominance, non-dominated sorting,
    /// crowding distance, hypervolume, and the generational distance indicators.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// All objectives are minimized. A point is given as one row of objective values, either as a
    /// jagged array <c>double[n][m]</c> or a rectangular array <c>double[n, m]</c> with n points and
    /// m objectives. To maximize an objective, negate it before calling these methods.
    /// </para>
    /// <para>
    /// <b> Input policy: </b>
    /// Arrays must not be null, every row must have the same positive number of objectives, and NaN is
    /// rejected with an <see cref="ArgumentException"/>. Dominance and sorting accept positive and
    /// negative infinity (for example, a penalized infeasible point). Crowding distance, hypervolume,
    /// and the distance indicators require finite values.
    /// </para>
    /// <para>
    ///     <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Deb, K., Pratap, A., Agarwal, S., and Meyarivan, T. (2002). A fast and elitist multiobjective genetic
    /// algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation, 6(2), 182-197.
    /// </description></item>
    /// <item><description>
    /// Zitzler, E., and Thiele, L. (1999). Multiobjective evolutionary algorithms: A comparative case study and
    /// the strength Pareto approach. IEEE Transactions on Evolutionary Computation, 3(4), 257-271.
    /// </description></item>
    /// <item><description>
    /// Beume, N., Fonseca, C. M., Lopez-Ibanez, M., Paquete, L., and Vahrenhold, J. (2009). On the complexity of
    /// computing the hypervolume indicator. IEEE Transactions on Evolutionary Computation, 13(5), 1075-1082.
    /// </description></item>
    /// <item><description>
    /// While, L., Bradstreet, L., and Barone, L. (2012). A fast way of calculating exact hypervolumes.
    /// IEEE Transactions on Evolutionary Computation, 16(1), 86-95.
    /// </description></item>
    /// <item><description>
    /// Van Veldhuizen, D. A., and Lamont, G. B. (2000). On measuring multiobjective evolutionary algorithm
    /// performance. Proceedings of the 2000 Congress on Evolutionary Computation, 204-211.
    /// </description></item>
    /// <item><description>
    /// Coello Coello, C. A., and Reyes Sierra, M. (2004). A study of the parallelization of a coevolutionary
    /// multi-objective evolutionary algorithm. MICAI 2004, Lecture Notes in Computer Science, 688-697.
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public static class ParetoFront
    {

        #region Dominance

        /// <summary>
        /// Determines whether point <paramref name="a"/> Pareto dominates point <paramref name="b"/> (minimization).
        /// </summary>
        /// <param name="a">The objective values of the first point.</param>
        /// <param name="b">The objective values of the second point.</param>
        /// <returns>True if <paramref name="a"/> is no worse than <paramref name="b"/> in every objective and strictly
        /// better in at least one; otherwise, false. Equal points do not dominate each other.</returns>
        /// <exception cref="ArgumentNullException">Thrown when either point is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the points are empty, have different lengths, or contain NaN.</exception>
        public static bool Dominates(double[] a, double[] b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));
            if (a.Length == 0) throw new ArgumentException("A point must have at least one objective.", nameof(a));
            if (a.Length != b.Length) throw new ArgumentException("The points must have the same number of objectives.", nameof(b));
            bool strict = false;
            for (int k = 0; k < a.Length; k++)
            {
                if (double.IsNaN(a[k])) throw new ArgumentException("Objective values cannot be NaN.", nameof(a));
                if (double.IsNaN(b[k])) throw new ArgumentException("Objective values cannot be NaN.", nameof(b));
            }
            for (int k = 0; k < a.Length; k++)
            {
                if (a[k] > b[k]) return false;
                if (a[k] < b[k]) strict = true;
            }
            return strict;
        }

        #endregion

        #region Non-Dominated Sorting

        /// <summary>
        /// Sorts points into non-dominated fronts (minimization).
        /// </summary>
        /// <param name="points">The points, one row of objective values per point.</param>
        /// <returns>The fronts, front 0 (the non-dominated set) first. Each front lists point indices in ascending order.
        /// An empty input returns an empty list.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="points"/> or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when rows are empty or ragged, or a value is NaN.</exception>
        /// <remarks>
        /// <para>
        /// The fronts are those of the fast non-dominated sort of Deb et al. (2002): front k + 1 holds the points
        /// dominated only by points of fronts 0 to k. Equivalently, the front of a point is one plus the largest front
        /// among the points that dominate it. Duplicate points share a front.
        /// </para>
        /// <para>
        /// The points are first ordered lexicographically, so a point can only be dominated by points earlier in that
        /// order, and the front of a point is found by a binary search over the fronts built so far, because "front k
        /// contains a dominator of p" is monotone in k. For two objectives each test looks at the last point added to
        /// the front, so the sort is O(N log N). For three objectives each test is a binary search in a staircase of the
        /// front's second and third objectives: O(N log^2 N) comparisons, with at most O(N^2) element moves. For four or
        /// more objectives each point is compared with every earlier point, in O(M N^2) time and O(M N) memory, without
        /// storing the domination lists of the original algorithm. Every path gives exactly the fronts of the definition.
        /// </para>
        /// </remarks>
        public static List<List<int>> NonDominatedSort(double[][] points)
        {
            var data = Flatten(points, nameof(points), false, out int n, out int m);
            return SortFlat(data, n, m);
        }

        /// <inheritdoc cref="NonDominatedSort(double[][])"/>
        public static List<List<int>> NonDominatedSort(double[,] points)
        {
            var data = Flatten(points, nameof(points), false, out int n, out int m);
            return SortFlat(data, n, m);
        }

        /// <summary>
        /// Returns the indices of the non-dominated points (front 0), in ascending order (minimization).
        /// </summary>
        /// <param name="points">The points, one row of objective values per point.</param>
        /// <returns>The ascending indices of the points that no other point dominates. Duplicates of a non-dominated
        /// point are all returned.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="points"/> or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when rows are empty or ragged, or a value is NaN.</exception>
        public static int[] NonDominatedIndices(double[][] points)
        {
            var data = Flatten(points, nameof(points), false, out int n, out int m);
            return FirstFrontFlat(data, n, m);
        }

        /// <inheritdoc cref="NonDominatedIndices(double[][])"/>
        public static int[] NonDominatedIndices(double[,] points)
        {
            var data = Flatten(points, nameof(points), false, out int n, out int m);
            return FirstFrontFlat(data, n, m);
        }

        /// <summary>
        /// Computes the front number of every point from the flattened data.
        /// </summary>
        private static int[] RanksFlat(double[] data, int n, int m, out int frontCount)
        {
            var rank = new int[n];
            frontCount = 0;
            if (n == 0) return rank;
            var order = LexicographicOrder(data, n, m);
            if (m == 1)
            {
                // Equal values share a front; each larger distinct value starts the next front.
                int r = 0;
                rank[order[0]] = 0;
                for (int i = 1; i < n; i++)
                {
                    if (data[order[i]] > data[order[i - 1]]) r++;
                    rank[order[i]] = r;
                }
                frontCount = r + 1;
            }
            else if (m == 2)
            {
                // Last point added to each front, in lexicographic order. That point has the
                // smallest second objective in its front, and "front k dominates p" is monotone in k.
                var lastX = new double[n];
                var lastY = new double[n];
                for (int i = 0; i < n; i++)
                {
                    int p = order[i];
                    double x = data[2 * p], y = data[2 * p + 1];
                    int lo = 0, hi = frontCount;
                    while (lo < hi)
                    {
                        int mid = (lo + hi) >> 1;
                        bool dominated = lastY[mid] < y || (lastY[mid] == y && lastX[mid] < x);
                        if (dominated) lo = mid + 1; else hi = mid;
                    }
                    rank[p] = lo;
                    lastX[lo] = x;
                    lastY[lo] = y;
                    if (lo == frontCount) frontCount++;
                }
            }
            else if (m == 3)
            {
                frontCount = RanksThree(data, n, order, rank);
            }
            else
            {
                // Earlier points (lexicographically) are the only possible dominators.
                var sorted = new double[n * m];
                for (int i = 0; i < n; i++) Array.Copy(data, order[i] * m, sorted, i * m, m);
                var sortedRank = new int[n];
                for (int i = 0; i < n; i++)
                {
                    int r = 0;
                    int bi = i * m;
                    for (int j = 0; j < i; j++)
                    {
                        if (sortedRank[j] < r) continue;
                        int bj = j * m;
                        bool weak = true, strict = false;
                        for (int k = 0; k < m; k++)
                        {
                            double aj = sorted[bj + k], ai = sorted[bi + k];
                            if (aj > ai) { weak = false; break; }
                            if (aj < ai) strict = true;
                        }
                        if (weak && strict) r = sortedRank[j] + 1;
                    }
                    sortedRank[i] = r;
                    rank[order[i]] = r;
                    if (r + 1 > frontCount) frontCount = r + 1;
                }
            }
            return rank;
        }

        /// <summary>
        /// Front numbers for three objectives, from the points in lexicographic order.
        /// </summary>
        /// <remarks>
        /// Every point already placed precedes the current point p lexicographically, so its first objective is no
        /// larger. A front then contains a dominator of p exactly when its two-dimensional staircase of
        /// (second, third) objectives contains an entry that is no worse than p in both, and strictly better in one of
        /// them or in the first objective. Each front keeps that staircase (second objective strictly ascending, third
        /// strictly descending, with the smallest first objective of each entry), and the front of p is found by a
        /// binary search over fronts, because "front k contains a dominator of p" is monotone in k. Each query is
        /// O(log N); an insertion may shift a staircase, so the worst case is O(N^2) element moves.
        /// </remarks>
        private static int RanksThree(double[] data, int n, int[] order, int[] rank)
        {
            var stairs = new List<Staircase>();
            for (int i = 0; i < n; i++)
            {
                int p = order[i];
                double x = data[3 * p], y = data[3 * p + 1], z = data[3 * p + 2];
                int lo = 0, hi = stairs.Count;
                while (lo < hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (stairs[mid].Dominates(x, y, z)) lo = mid + 1; else hi = mid;
                }
                if (lo == stairs.Count) stairs.Add(new Staircase());
                stairs[lo].Insert(x, y, z);
                rank[p] = lo;
            }
            return stairs.Count;
        }

        /// <summary>
        /// The two-dimensional staircase of (second, third) objectives of one front, used by <see cref="RanksThree"/>.
        /// </summary>
        private sealed class Staircase
        {
            private double[] _x = new double[4], _y = new double[4], _z = new double[4];
            private int _count;

            /// <summary>
            /// Whether a stored point dominates (x, y, z), given that every stored point has a first objective no larger than x.
            /// </summary>
            public bool Dominates(double x, double y, double z)
            {
                // Last entry with second objective no larger than y: it has the smallest third objective among them.
                int lo = 0, hi = _count;
                while (lo < hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (_y[mid] <= y) lo = mid + 1; else hi = mid;
                }
                int e = lo - 1;
                if (e < 0 || _z[e] > z) return false;
                if (_z[e] < z || _y[e] < y) return true;
                return _x[e] < x;
            }

            /// <summary>
            /// Adds a point, keeping only the two-dimensional minima.
            /// </summary>
            public void Insert(double x, double y, double z)
            {
                int lo = 0, hi = _count;
                while (lo < hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (_y[mid] < y) lo = mid + 1; else hi = mid;
                }
                int pos = lo;
                if (pos > 0 && _z[pos - 1] <= z) return;
                if (pos < _count && _y[pos] == y && _z[pos] <= z) return;
                int end = pos;
                while (end < _count && _z[end] >= z) end++;
                int shift = end - pos - 1;
                if (shift < 0 && _count == _y.Length)
                {
                    int cap = 2 * _count;
                    Array.Resize(ref _x, cap);
                    Array.Resize(ref _y, cap);
                    Array.Resize(ref _z, cap);
                }
                if (shift != 0)
                {
                    Array.Copy(_x, end, _x, pos + 1, _count - end);
                    Array.Copy(_y, end, _y, pos + 1, _count - end);
                    Array.Copy(_z, end, _z, pos + 1, _count - end);
                    _count -= shift;
                }
                _x[pos] = x;
                _y[pos] = y;
                _z[pos] = z;
            }
        }

        /// <summary>
        /// Builds the fronts from the flattened data.
        /// </summary>
        private static List<List<int>> SortFlat(double[] data, int n, int m)
        {
            var rank = RanksFlat(data, n, m, out int frontCount);
            var fronts = new List<List<int>>(frontCount);
            for (int f = 0; f < frontCount; f++) fronts.Add(new List<int>());
            for (int i = 0; i < n; i++) fronts[rank[i]].Add(i);
            return fronts;
        }

        /// <summary>
        /// Returns the first front from the flattened data.
        /// </summary>
        private static int[] FirstFrontFlat(double[] data, int n, int m)
        {
            if (n == 0) return new int[0];
            var order = LexicographicOrder(data, n, m);
            // A point is non-dominated if no earlier point (lexicographically) dominates it,
            // and it suffices to check the non-dominated points kept so far.
            var kept = new List<int>();
            var isFirst = new bool[n];
            for (int i = 0; i < n; i++)
            {
                int p = order[i];
                int bp = p * m;
                bool dominated = false;
                for (int t = 0; t < kept.Count && !dominated; t++)
                {
                    int bq = kept[t] * m;
                    bool weak = true, strict = false;
                    for (int k = 0; k < m; k++)
                    {
                        double aq = data[bq + k], ap = data[bp + k];
                        if (aq > ap) { weak = false; break; }
                        if (aq < ap) strict = true;
                    }
                    dominated = weak && strict;
                }
                if (!dominated)
                {
                    kept.Add(p);
                    isFirst[p] = true;
                }
            }
            var result = new int[kept.Count];
            int c = 0;
            for (int i = 0; i < n; i++) if (isFirst[i]) result[c++] = i;
            return result;
        }

        #endregion

        #region Crowding Distance

        /// <summary>
        /// Computes the crowding distance of every point of one front (Deb et al. 2002).
        /// </summary>
        /// <param name="front">The points of the front, one row of finite objective values per point.</param>
        /// <returns>The crowding distance of each point, aligned with the rows of <paramref name="front"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="front"/> or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when rows are empty or ragged, or a value is NaN or infinite.</exception>
        /// <remarks>
        /// <para>
        /// For each objective m, the points are ordered by that objective (ties by index). The first and last points
        /// get positive infinity, and every interior point i adds (f_m[next] - f_m[previous]) / (max f_m - min f_m),
        /// so each objective is normalized by its range over the front. An objective whose range is zero contributes
        /// nothing to any point, boundary points included. Every point of a front with one or two points gets positive
        /// infinity.
        /// </para>
        /// <para>
        /// This is the sum over objectives of Deb et al. (2002). pymoo 0.6.2 (<c>calc_crowding_distance</c>) uses the
        /// same terms but divides the sum by the number of objectives. By default (<c>filter_out_duplicates=True</c>)
        /// pymoo also removes duplicate points before computing the distances and gives each removed duplicate a
        /// distance of zero; this method keeps every point, so the values differ from pymoo's on a front with duplicates.
        /// </para>
        /// </remarks>
        public static double[] CrowdingDistance(double[][] front)
        {
            var data = Flatten(front, nameof(front), true, out int n, out int m);
            return CrowdingFlat(data, n, m);
        }

        /// <inheritdoc cref="CrowdingDistance(double[][])"/>
        public static double[] CrowdingDistance(double[,] front)
        {
            var data = Flatten(front, nameof(front), true, out int n, out int m);
            return CrowdingFlat(data, n, m);
        }

        /// <summary>
        /// Computes the crowding distance of the points of one front given by indices into a population.
        /// </summary>
        /// <param name="points">The population, one row of objective values per point.</param>
        /// <param name="frontIndices">The indices of the front members, for example one front returned by
        /// <see cref="NonDominatedSort(double[][])"/>.</param>
        /// <returns>The crowding distance of each front member, aligned with <paramref name="frontIndices"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when an argument or a referenced row is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when an index is outside the population.</exception>
        /// <exception cref="ArgumentException">Thrown when referenced rows are empty or ragged, or a value is NaN or infinite.</exception>
        /// <remarks>
        /// Only the referenced rows are read and validated. The definition is that of <see cref="CrowdingDistance(double[][])"/>.
        /// </remarks>
        public static double[] CrowdingDistance(double[][] points, IList<int> frontIndices)
        {
            if (points == null) throw new ArgumentNullException(nameof(points));
            if (frontIndices == null) throw new ArgumentNullException(nameof(frontIndices));
            var rows = new double[frontIndices.Count][];
            for (int i = 0; i < rows.Length; i++)
            {
                int idx = frontIndices[i];
                if (idx < 0 || idx >= points.Length) throw new ArgumentOutOfRangeException(nameof(frontIndices), "A front index is outside the population.");
                rows[i] = points[idx];
            }
            var data = Flatten(rows, nameof(points), true, out int n, out int m);
            return CrowdingFlat(data, n, m);
        }

        /// <summary>
        /// Crowding distance from the flattened data.
        /// </summary>
        private static double[] CrowdingFlat(double[] data, int n, int m)
        {
            var cd = new double[n];
            if (n <= 2)
            {
                for (int i = 0; i < n; i++) cd[i] = double.PositiveInfinity;
                return cd;
            }
            var idx = new int[n];
            var keys = new double[n];
            for (int k = 0; k < m; k++)
            {
                for (int i = 0; i < n; i++)
                {
                    idx[i] = i;
                    keys[i] = data[i * m + k];
                }
                StableSortByKey(keys, idx);
                double range = keys[n - 1] - keys[0];
                if (range == 0) continue;
                cd[idx[0]] = double.PositiveInfinity;
                cd[idx[n - 1]] = double.PositiveInfinity;
                for (int i = 1; i < n - 1; i++)
                    cd[idx[i]] += (keys[i + 1] - keys[i - 1]) / range;
            }
            return cd;
        }

        #endregion

        #region Hypervolume

        /// <summary>
        /// Computes the exact hypervolume dominated by a set of points and bounded by a reference point (minimization).
        /// </summary>
        /// <param name="points">The points, one row of finite objective values per point.</param>
        /// <param name="referencePoint">The finite reference point, with one value per objective.</param>
        /// <returns>The Lebesgue measure of the union of the boxes [p, r] over the points p that strictly dominate the
        /// reference point r in every objective. Points that do not are ignored, and an empty set gives zero.</returns>
        /// <exception cref="ArgumentNullException">Thrown when an argument or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when rows are empty or ragged, the reference point length differs
        /// from the number of objectives, or a value is NaN or infinite.</exception>
        /// <remarks>
        /// <para>
        /// The hypervolume indicator is that of Zitzler and Thiele (1999). Dominated and duplicate points do not change it.
        /// </para>
        /// <para>
        /// One objective: r - min f. Two objectives: sort and sweep, O(N log N). Three objectives: a sweep in the third
        /// objective that keeps the two-dimensional staircase of the points seen so far and updates its area by the
        /// exclusive area of each insertion (Beume et al. 2009, with an array in place of a balanced tree): O(N log N)
        /// comparisons and at most O(N^2) element moves. Four or more objectives: the WFG algorithm of While et al. (2012),
        /// slicing on the last objective, with exclusive hypervolumes computed from non-dominated limit sets and the
        /// three-objective sweep as the base case. The worst case of WFG is exponential in the number of objectives.
        /// </para>
        /// </remarks>
        public static double Hypervolume(double[][] points, double[] referencePoint)
        {
            var data = Flatten(points, nameof(points), true, out int n, out int m);
            return HypervolumeFlat(data, n, m, referencePoint);
        }

        /// <inheritdoc cref="Hypervolume(double[][], double[])"/>
        public static double Hypervolume(double[,] points, double[] referencePoint)
        {
            var data = Flatten(points, nameof(points), true, out int n, out int m);
            return HypervolumeFlat(data, n, m, referencePoint);
        }

        /// <summary>
        /// Validates the reference point, filters the points that strictly dominate it, and dispatches by dimension.
        /// </summary>
        private static double HypervolumeFlat(double[] data, int n, int m, double[] referencePoint)
        {
            if (referencePoint == null) throw new ArgumentNullException(nameof(referencePoint));
            if (n > 0 && referencePoint.Length != m)
                throw new ArgumentException("The reference point must have one value per objective.", nameof(referencePoint));
            for (int k = 0; k < referencePoint.Length; k++)
            {
                if (!Tools.IsFinite(referencePoint[k]))
                    throw new ArgumentException("The reference point must be finite.", nameof(referencePoint));
            }
            if (n == 0) return 0d;
            // Keep the points strictly better than the reference point in every objective.
            var kept = new double[n * m];
            int c = 0;
            for (int i = 0; i < n; i++)
            {
                bool inside = true;
                for (int k = 0; k < m; k++)
                {
                    if (!(data[i * m + k] < referencePoint[k])) { inside = false; break; }
                }
                if (inside)
                {
                    Array.Copy(data, i * m, kept, c * m, m);
                    c++;
                }
            }
            if (c == 0) return 0d;
            return HvRecursive(kept, c, m, referencePoint);
        }

        /// <summary>
        /// Hypervolume of <paramref name="n"/> points stored row-major with stride <paramref name="d"/>, all strictly
        /// inside the reference box, using the first <paramref name="d"/> reference coordinates.
        /// </summary>
        /// <remarks>
        /// For d of four or more this is the WFG recursion. With the points sorted by the last objective ascending,
        /// HV = sum_i (r_d - z_i) * [inclhv(p_i) - HV_(d-1)(nds(limit({p_0..p_(i-1)}, p_i)))], where limit replaces each
        /// earlier point q by max(q, p_i) in the first d - 1 objectives and nds keeps the non-dominated distinct points.
        /// The input array is not modified.
        /// </remarks>
        private static double HvRecursive(double[] pts, int n, int d, double[] r)
        {
            if (n == 0) return 0d;
            if (n == 1)
            {
                double v = 1d;
                for (int k = 0; k < d; k++) v *= r[k] - pts[k];
                return v;
            }
            if (d == 1)
            {
                double min = pts[0];
                for (int i = 1; i < n; i++) if (pts[i] < min) min = pts[i];
                return r[0] - min;
            }
            if (d == 2) return Hv2(pts, n, r);
            if (d == 3) return Hv3(pts, n, r);

            // Sort by the last objective ascending; any order of ties gives the same volume.
            int last = d - 1;
            var order = new int[n];
            var keys = new double[n];
            for (int i = 0; i < n; i++) { order[i] = i; keys[i] = pts[i * d + last]; }
            Array.Sort(keys, order);
            var sorted = new double[n * d];
            for (int i = 0; i < n; i++) Array.Copy(pts, order[i] * d, sorted, i * d, d);

            int dm = d - 1;
            var limit = new double[n * dm];
            var removed = new bool[n];
            double volume = 0d;
            for (int i = 0; i < n; i++)
            {
                int bi = i * d;
                double height = r[last] - sorted[bi + last];
                double incl = 1d;
                for (int k = 0; k < dm; k++) incl *= r[k] - sorted[bi + k];

                // Limit set of the earlier points, projected to d - 1 objectives.
                for (int j = 0; j < i; j++)
                {
                    int bj = j * d, lj = j * dm;
                    for (int k = 0; k < dm; k++)
                    {
                        double q = sorted[bj + k], p = sorted[bi + k];
                        limit[lj + k] = q > p ? q : p;
                    }
                }
                int count = NonDominatedCompact(limit, i, dm, removed);
                double excl = incl - HvRecursive(limit, count, dm, r);
                volume += height * excl;
            }
            return volume;
        }

        /// <summary>
        /// Removes dominated and duplicate rows in place and returns the number of rows kept at the front of the array.
        /// </summary>
        private static int NonDominatedCompact(double[] pts, int n, int d, bool[] removed)
        {
            for (int i = 0; i < n; i++) removed[i] = false;
            for (int i = 0; i < n; i++)
            {
                if (removed[i]) continue;
                int bi = i * d;
                for (int j = i + 1; j < n; j++)
                {
                    if (removed[j]) continue;
                    int bj = j * d;
                    bool iLeq = true, jLeq = true;
                    for (int k = 0; k < d; k++)
                    {
                        double a = pts[bi + k], b = pts[bj + k];
                        if (a > b) iLeq = false;
                        else if (a < b) jLeq = false;
                        if (!iLeq && !jLeq) break;
                    }
                    if (iLeq) removed[j] = true;        // i weakly dominates j (or equals it)
                    else if (jLeq) { removed[i] = true; break; }
                }
            }
            int c = 0;
            for (int i = 0; i < n; i++)
            {
                if (removed[i]) continue;
                if (c != i) Array.Copy(pts, i * d, pts, c * d, d);
                c++;
            }
            return c;
        }

        /// <summary>
        /// Two-objective hypervolume by sort and sweep. Points are stored with stride 2 and lie strictly inside the box.
        /// </summary>
        private static double Hv2(double[] pts, int n, double[] r)
        {
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++) { x[i] = pts[2 * i]; y[i] = pts[2 * i + 1]; }
            Array.Sort(x, y);
            double area = 0d, prevY = r[1];
            int s = 0;
            while (s < n)
            {
                // Among equal first objectives only the smallest second objective matters.
                double xs = x[s], ys = y[s];
                int e = s + 1;
                while (e < n && x[e] == xs) { if (y[e] < ys) ys = y[e]; e++; }
                if (ys < prevY)
                {
                    area += (r[0] - xs) * (prevY - ys);
                    prevY = ys;
                }
                s = e;
            }
            return area;
        }

        /// <summary>
        /// Three-objective hypervolume by a sweep in the third objective over an incrementally maintained
        /// two-dimensional staircase. Points are stored with stride 3 and lie strictly inside the box.
        /// </summary>
        private static double Hv3(double[] pts, int n, double[] r)
        {
            var order = new int[n];
            var z = new double[n];
            for (int i = 0; i < n; i++) { order[i] = i; z[i] = pts[3 * i + 2]; }
            Array.Sort(z, order);

            // Staircase sorted by x ascending and y strictly descending.
            var sx = new double[n];
            var sy = new double[n];
            int c = 0;
            double area = 0d, volume = 0d;
            for (int t = 0; t < n; t++)
            {
                int p = order[t];
                double x = pts[3 * p], y = pts[3 * p + 1];

                // pos = first staircase index with sx >= x.
                int lo = 0, hi = c;
                while (lo < hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (sx[mid] < x) lo = mid + 1; else hi = mid;
                }
                int pos = lo;
                bool dominated = (pos > 0 && sy[pos - 1] <= y) || (pos < c && sx[pos] == x && sy[pos] <= y);
                if (!dominated)
                {
                    // Points from pos with sy >= y are dominated by the new point.
                    int end = pos;
                    while (end < c && sy[end] >= y) end++;
                    double curX = x, curH = pos > 0 ? sy[pos - 1] : r[1];
                    double delta = 0d;
                    for (int k = pos; k < end; k++)
                    {
                        delta += (sx[k] - curX) * (curH - y);
                        curX = sx[k];
                        curH = sy[k];
                    }
                    double endX = end < c ? sx[end] : r[0];
                    delta += (endX - curX) * (curH - y);
                    area += delta;

                    // Replace [pos, end) by the new point.
                    int shift = end - pos - 1;
                    if (shift != 0)
                    {
                        Array.Copy(sx, end, sx, pos + 1, c - end);
                        Array.Copy(sy, end, sy, pos + 1, c - end);
                        c -= shift;
                    }
                    sx[pos] = x;
                    sy[pos] = y;
                }
                double zNext = t + 1 < n ? z[t + 1] : r[2];
                volume += area * (zNext - z[t]);
            }
            return volume;
        }

        #endregion

        #region Distance Indicators

        /// <summary>
        /// Computes the generational distance (GD) of an approximation set to a reference front.
        /// </summary>
        /// <param name="approximation">The approximation set, one row of finite objective values per point.</param>
        /// <param name="referenceFront">The reference (true) Pareto front, one row per point.</param>
        /// <returns>The mean, over the approximation points, of the Euclidean distance to the nearest reference point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when an argument or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a set is empty, rows are empty or ragged, the sets have different
        /// numbers of objectives, or a value is NaN or infinite.</exception>
        /// <remarks>
        /// GD = (1/n) sum_i min_j ||a_i - z_j||, computed in O(n R M) for R reference points. This matches pymoo 0.6.2
        /// <c>GD</c> without normalization. Van Veldhuizen and Lamont (2000) define GD as (sum_i d_i^2)^(1/2) / n, which
        /// differs from this mean distance; objectives are not normalized here, so scale them first if their units differ.
        /// </remarks>
        public static double GenerationalDistance(double[][] approximation, double[][] referenceFront)
        {
            var a = Flatten(approximation, nameof(approximation), true, out int na, out int ma);
            var z = Flatten(referenceFront, nameof(referenceFront), true, out int nz, out int mz);
            CheckIndicatorSets(na, ma, nz, mz);
            return MeanNearestDistance(a, na, z, nz, ma);
        }

        /// <inheritdoc cref="GenerationalDistance(double[][], double[][])"/>
        public static double GenerationalDistance(double[,] approximation, double[,] referenceFront)
        {
            var a = Flatten(approximation, nameof(approximation), true, out int na, out int ma);
            var z = Flatten(referenceFront, nameof(referenceFront), true, out int nz, out int mz);
            CheckIndicatorSets(na, ma, nz, mz);
            return MeanNearestDistance(a, na, z, nz, ma);
        }

        /// <summary>
        /// Computes the inverted generational distance (IGD) of an approximation set to a reference front.
        /// </summary>
        /// <param name="approximation">The approximation set, one row of finite objective values per point.</param>
        /// <param name="referenceFront">The reference (true) Pareto front, one row per point.</param>
        /// <returns>The mean, over the reference points, of the Euclidean distance to the nearest approximation point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when an argument or a row is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a set is empty, rows are empty or ragged, the sets have different
        /// numbers of objectives, or a value is NaN or infinite.</exception>
        /// <remarks>
        /// IGD = (1/R) sum_j min_i ||z_j - a_i|| (Coello Coello and Reyes Sierra 2004), computed in O(n R M). This matches
        /// pymoo 0.6.2 <c>IGD</c> without normalization. Objectives are not normalized here.
        /// </remarks>
        public static double InvertedGenerationalDistance(double[][] approximation, double[][] referenceFront)
        {
            var a = Flatten(approximation, nameof(approximation), true, out int na, out int ma);
            var z = Flatten(referenceFront, nameof(referenceFront), true, out int nz, out int mz);
            CheckIndicatorSets(na, ma, nz, mz);
            return MeanNearestDistance(z, nz, a, na, ma);
        }

        /// <inheritdoc cref="InvertedGenerationalDistance(double[][], double[][])"/>
        public static double InvertedGenerationalDistance(double[,] approximation, double[,] referenceFront)
        {
            var a = Flatten(approximation, nameof(approximation), true, out int na, out int ma);
            var z = Flatten(referenceFront, nameof(referenceFront), true, out int nz, out int mz);
            CheckIndicatorSets(na, ma, nz, mz);
            return MeanNearestDistance(z, nz, a, na, ma);
        }

        /// <summary>
        /// Validates that both indicator sets are non-empty and have the same number of objectives.
        /// </summary>
        private static void CheckIndicatorSets(int na, int ma, int nz, int mz)
        {
            if (na == 0) throw new ArgumentException("The approximation set must contain at least one point.", "approximation");
            if (nz == 0) throw new ArgumentException("The reference front must contain at least one point.", "referenceFront");
            if (ma != mz) throw new ArgumentException("The approximation set and reference front must have the same number of objectives.", "referenceFront");
        }

        /// <summary>
        /// Mean over the rows of <paramref name="from"/> of the Euclidean distance to the nearest row of <paramref name="to"/>.
        /// </summary>
        private static double MeanNearestDistance(double[] from, int nFrom, double[] to, int nTo, int m)
        {
            double sum = 0d;
            for (int i = 0; i < nFrom; i++)
            {
                int bi = i * m;
                double best = double.PositiveInfinity;
                for (int j = 0; j < nTo; j++)
                {
                    int bj = j * m;
                    double s = 0d;
                    for (int k = 0; k < m && s < best; k++)
                    {
                        double diff = from[bi + k] - to[bj + k];
                        s += diff * diff;
                    }
                    if (s < best) best = s;
                }
                sum += Math.Sqrt(best);
            }
            return sum / nFrom;
        }

        #endregion

        #region Helpers

        /// <summary>
        /// Validates a jagged point array and copies it to a row-major buffer.
        /// </summary>
        private static double[] Flatten(double[][] points, string paramName, bool requireFinite, out int n, out int m)
        {
            if (points == null) throw new ArgumentNullException(paramName);
            n = points.Length;
            m = 0;
            if (n == 0) return new double[0];
            if (points[0] == null) throw new ArgumentNullException(paramName, "A point cannot be null.");
            m = points[0].Length;
            if (m == 0) throw new ArgumentException("Each point must have at least one objective.", paramName);
            var data = new double[n * m];
            for (int i = 0; i < n; i++)
            {
                var row = points[i];
                if (row == null) throw new ArgumentNullException(paramName, "A point cannot be null.");
                if (row.Length != m) throw new ArgumentException("Every point must have the same number of objectives.", paramName);
                for (int k = 0; k < m; k++)
                {
                    double v = row[k];
                    CheckValue(v, paramName, requireFinite);
                    data[i * m + k] = v;
                }
            }
            return data;
        }

        /// <summary>
        /// Validates a rectangular point array and copies it to a row-major buffer.
        /// </summary>
        private static double[] Flatten(double[,] points, string paramName, bool requireFinite, out int n, out int m)
        {
            if (points == null) throw new ArgumentNullException(paramName);
            n = points.GetLength(0);
            m = points.GetLength(1);
            if (n == 0) { m = 0; return new double[0]; }
            if (m == 0) throw new ArgumentException("Each point must have at least one objective.", paramName);
            var data = new double[n * m];
            for (int i = 0; i < n; i++)
            {
                for (int k = 0; k < m; k++)
                {
                    double v = points[i, k];
                    CheckValue(v, paramName, requireFinite);
                    data[i * m + k] = v;
                }
            }
            return data;
        }

        /// <summary>
        /// Rejects NaN, and infinity when finite values are required.
        /// </summary>
        private static void CheckValue(double v, string paramName, bool requireFinite)
        {
            if (double.IsNaN(v)) throw new ArgumentException("Objective values cannot be NaN.", paramName);
            if (requireFinite && double.IsInfinity(v)) throw new ArgumentException("Objective values must be finite.", paramName);
        }

        /// <summary>
        /// Returns the point indices in lexicographic order of their objective rows, ties by index.
        /// </summary>
        private static int[] LexicographicOrder(double[] data, int n, int m)
        {
            var order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            if (m == 1)
            {
                var keys = new double[n];
                Array.Copy(data, keys, n);
                StableSortByKey(keys, order);
                return order;
            }
            Array.Sort(order, (a, b) =>
            {
                int ba = a * m, bb = b * m;
                for (int k = 0; k < m; k++)
                {
                    double x = data[ba + k], y = data[bb + k];
                    if (x < y) return -1;
                    if (x > y) return 1;
                }
                return a.CompareTo(b);
            });
            return order;
        }

        /// <summary>
        /// Sorts <paramref name="keys"/> ascending and applies the same permutation to <paramref name="items"/>,
        /// breaking ties by the item value so the order is stable for index items that start in ascending order.
        /// </summary>
        private static void StableSortByKey(double[] keys, int[] items)
        {
            int n = keys.Length;
            var pairs = new KeyValuePair<double, int>[n];
            for (int i = 0; i < n; i++) pairs[i] = new KeyValuePair<double, int>(keys[i], items[i]);
            Array.Sort(pairs, (p, q) =>
            {
                int c = p.Key.CompareTo(q.Key);
                return c != 0 ? c : p.Value.CompareTo(q.Value);
            });
            for (int i = 0; i < n; i++) { keys[i] = pairs[i].Key; items[i] = pairs[i].Value; }
        }

        #endregion

    }
}
