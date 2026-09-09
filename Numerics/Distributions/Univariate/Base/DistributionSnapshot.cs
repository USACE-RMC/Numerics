using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Numerics.Distributions
{
    /// <summary>
    /// An immutable bitwise snapshot of a distribution's mutable configuration whose equality
    /// implies an identical canonical configuration string.
    /// </summary>
    /// <remarks>
    /// Composite owners publish a snapshot after refreshing their caches and compare it against
    /// live state on later evaluations, replacing per-call canonical-string serialization with a
    /// scalar walk that allocates nothing. Capture succeeds only for exact built-in types whose
    /// canonical form is fully determined by the captured scalars; derived types and
    /// table-backed families return <see langword="null"/> so their owners retain the generic
    /// canonical-string path. A bitwise match implies the canonical configuration string is
    /// unchanged; a mismatch falls back to the string comparison, which remains the deciding
    /// authority, so incomplete capture can never invalidate correctly cached state.
    /// </remarks>
    internal sealed class DistributionSnapshot
    {
        /// <summary>The exact runtime type of the captured distribution.</summary>
        private readonly Type _type;

        /// <summary>The jump-table family code for the captured type, from <see cref="FamilyOf"/>.</summary>
        private readonly byte _family;

        /// <summary>The captured scalar state, as raw bits, in the fixed per-family order.</summary>
        private readonly long[] _scalarBits;

        /// <summary>Captured child snapshots for composite nodes; otherwise <see langword="null"/>.</summary>
        private readonly DistributionSnapshot[]? _children;

        /// <summary>Initializes an immutable snapshot node.</summary>
        /// <param name="type">The exact runtime type of the captured distribution.</param>
        /// <param name="family">The jump-table family code for the captured type.</param>
        /// <param name="scalarBits">The captured scalar state in the fixed per-family order.</param>
        /// <param name="children">Child snapshots for composite nodes, or <see langword="null"/>.</param>
        private DistributionSnapshot(Type type, byte family, long[] scalarBits, DistributionSnapshot[]? children)
        {
            _type = type;
            _family = family;
            _scalarBits = scalarBits;
            _children = children;
        }

        /// <summary>Captures the mutable state that determines a distribution's canonical configuration.</summary>
        /// <param name="distribution">The distribution to capture.</param>
        /// <returns>An immutable snapshot, or <see langword="null"/> when any node in the tree is not an exact supported built-in.</returns>
        internal static DistributionSnapshot? TryCapture(UnivariateDistributionBase? distribution)
        {
            if (distribution is null || !IsSupportedTree(distribution)) return null;
            return Capture(distribution);
        }

        /// <summary>Compares live state bitwise without allocating wrappers, parameter arrays, or strings.</summary>
        /// <param name="distribution">The distribution whose live state is compared with this snapshot.</param>
        /// <returns><see langword="true"/> when every captured scalar and child is bitwise unchanged; otherwise, <see langword="false"/>.</returns>
        /// <remarks>The compare arms are written inline against the stored bits rather than through
        /// the capture cursor: the consuming libraries link the Debug build of this assembly, whose
        /// minimal-optimization jitting keeps every helper call, so the hot path minimizes call
        /// count. A capture/compare round-trip test per supported family guards the two switches
        /// against drifting apart.</remarks>
        internal bool Matches(UnivariateDistributionBase? distribution)
        {
            if (distribution is null || distribution.GetType() != _type) return false;
            var bits = _scalarBits;
            int index = 0;
            if (!MatchScalars(distribution, _family, bits, ref index) || index != bits.Length) return false;
            if (_children is null) return true;
            var children = ChildrenOf(distribution);
            if (children is null || children.Length != _children.Length) return false;
            for (int i = 0; i < children.Length; i++)
                if (!_children[i].Matches(children[i])) return false;
            return true;
        }

        /// <summary>Compares one node's live scalars against stored bits with inline arithmetic.</summary>
        /// <param name="distribution">The node whose scalars are compared; its exact type must already be verified.</param>
        /// <param name="family">The node's family code, stored at capture for jump-table dispatch.</param>
        /// <param name="bits">The stored scalar bits.</param>
        /// <param name="index">The read position within <paramref name="bits"/>, advanced past this node's scalars.</param>
        /// <returns><see langword="true"/> when every scalar is bitwise unchanged; otherwise, <see langword="false"/>.</returns>
        private static bool MatchScalars(UnivariateDistributionBase distribution, byte family, long[] bits, ref int index)
        {
            int k = index;
            switch (family)
            {
                case 1:
                {
                    var x = (Normal)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Mu)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Sigma)) return false;
                    index = k + 2;
                    return true;
                }
                case 2:
                {
                    var x = (LogNormal)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Mu)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Sigma)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Base)) return false;
                    index = k + 3;
                    return true;
                }
                case 3:
                {
                    var x = (LnNormal)distribution;
                    if (k + 5 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Mu)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Sigma)
                        || bits[k + 2] != (x.PhysicalMomentModeForSnapshot ? 1L : 0L)
                        || bits[k + 3] != BitConverter.DoubleToInt64Bits(x.PhysicalMeanForSnapshot)
                        || bits[k + 4] != BitConverter.DoubleToInt64Bits(x.PhysicalStandardDeviationForSnapshot)) return false;
                    index = k + 5;
                    return true;
                }
                case 4:
                {
                    var x = (LogPearsonTypeIII)distribution;
                    if (k + 4 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Mu)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Sigma)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Gamma)
                        || bits[k + 3] != BitConverter.DoubleToInt64Bits(x.Base)) return false;
                    index = k + 4;
                    return true;
                }
                case 5:
                {
                    var x = (PearsonTypeIII)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Mu)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Sigma)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Gamma)) return false;
                    index = k + 3;
                    return true;
                }
                case 6:
                {
                    var x = (GammaDistribution)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Theta)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 2;
                    return true;
                }
                case 7:
                {
                    var x = (Weibull)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Lambda)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 2;
                    return true;
                }
                case 8:
                {
                    var x = (Gumbel)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)) return false;
                    index = k + 2;
                    return true;
                }
                case 9:
                {
                    var x = (GeneralizedExtremeValue)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 3;
                    return true;
                }
                case 10:
                {
                    var x = (GeneralizedLogistic)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 3;
                    return true;
                }
                case 11:
                {
                    var x = (GeneralizedNormal)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 3;
                    return true;
                }
                case 12:
                {
                    var x = (GeneralizedPareto)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Kappa)) return false;
                    index = k + 3;
                    return true;
                }
                case 13:
                {
                    var x = (Exponential)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)) return false;
                    index = k + 2;
                    return true;
                }
                case 14:
                {
                    var x = (KappaFour)distribution;
                    if (k + 4 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Kappa)
                        || bits[k + 3] != BitConverter.DoubleToInt64Bits(x.Hondo)) return false;
                    index = k + 4;
                    return true;
                }
                case 15:
                {
                    var x = (Uniform)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Min)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Max)) return false;
                    index = k + 2;
                    return true;
                }
                case 16:
                {
                    var x = (Triangular)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Min)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.MostLikely)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Max)) return false;
                    index = k + 3;
                    return true;
                }
                case 17:
                {
                    var x = (Pert)distribution;
                    if (k + 3 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Min)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.MostLikely)
                        || bits[k + 2] != BitConverter.DoubleToInt64Bits(x.Max)) return false;
                    index = k + 3;
                    return true;
                }
                case 18:
                {
                    var x = (Deterministic)distribution;
                    if (k + 1 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Value)) return false;
                    index = k + 1;
                    return true;
                }
                case 19:
                {
                    var x = (Logistic)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.Xi)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Alpha)) return false;
                    index = k + 2;
                    return true;
                }
                case 20:
                {
                    var x = (Cauchy)distribution;
                    if (k + 2 > bits.Length || bits[k] != BitConverter.DoubleToInt64Bits(x.X0)
                        || bits[k + 1] != BitConverter.DoubleToInt64Bits(x.Gamma)) return false;
                    index = k + 2;
                    return true;
                }
                case 21:
                {
                    var x = (Mixture)distribution;
                    var weights = x.Weights;
                    if (weights is null || k + 5 + weights.Length > bits.Length
                        || bits[k] != (x.IsZeroInflated ? 1L : 0L)
                        || bits[k + 1] != (long)x.XTransform
                        || bits[k + 2] != (long)x.ProbabilityTransform
                        || bits[k + 3] != BitConverter.DoubleToInt64Bits(x.ZeroWeight)
                        || bits[k + 4] != weights.Length) return false;
                    k += 5;
                    for (int i = 0; i < weights.Length; i++, k++)
                        if (bits[k] != BitConverter.DoubleToInt64Bits(weights[i])) return false;
                    index = k;
                    return true;
                }
                case 22:
                {
                    var x = (CompetingRisks)distribution;
                    if (k + 6 > bits.Length
                        || bits[k] != (x.MinimumOfRandomVariables ? 1L : 0L)
                        || bits[k + 1] != (long)x.Dependency
                        || bits[k + 2] != x.PRNGSeed
                        || bits[k + 3] != (long)x.XTransform
                        || bits[k + 4] != (long)x.ProbabilityTransform) return false;
                    k += 5;
                    var matrix = x.CorrelationMatrixArray;
                    if (matrix is null)
                    {
                        if (bits[k] != 0L) return false;
                        index = k + 1;
                        return true;
                    }
                    int rows = matrix.GetLength(0), columns = matrix.GetLength(1);
                    if (k + 3 + rows * columns > bits.Length || bits[k] != 1L
                        || bits[k + 1] != rows || bits[k + 2] != columns) return false;
                    k += 3;
                    int rowStart = matrix.GetLowerBound(0), columnStart = matrix.GetLowerBound(1);
                    for (int row = 0; row < rows; row++)
                        for (int column = 0; column < columns; column++, k++)
                            if (bits[k] != BitConverter.DoubleToInt64Bits(matrix[rowStart + row, columnStart + column])) return false;
                    index = k;
                    return true;
                }
            }
            return false;
        }

        /// <summary>Walks the tree checking that every node is an exact supported built-in, without allocating.</summary>
        /// <param name="distribution">The root of the tree to check.</param>
        /// <returns><see langword="true"/> when every node is capturable; otherwise, <see langword="false"/>.</returns>
        private static bool IsSupportedTree(UnivariateDistributionBase? distribution)
        {
            if (distribution is null) return false;
            Type type = distribution.GetType();
            if (type == typeof(Mixture) || type == typeof(CompetingRisks))
            {
                var children = ChildrenOf(distribution);
                if (children is null || children.Length == 0) return false;
                for (int i = 0; i < children.Length; i++)
                    if (!IsSupportedTree(children[i])) return false;
                return true;
            }
            return FamilyOf(type) != 0;
        }

        /// <summary>Maps an exact runtime type to its jump-table family code, or zero when unsupported.</summary>
        /// <param name="type">The exact runtime type to classify.</param>
        /// <returns>The family code for a supported type; otherwise, zero.</returns>
        private static byte FamilyOf(Type type)
        {
            if (type == typeof(Normal)) return 1;
            if (type == typeof(LogNormal)) return 2;
            if (type == typeof(LnNormal)) return 3;
            if (type == typeof(LogPearsonTypeIII)) return 4;
            if (type == typeof(PearsonTypeIII)) return 5;
            if (type == typeof(GammaDistribution)) return 6;
            if (type == typeof(Weibull)) return 7;
            if (type == typeof(Gumbel)) return 8;
            if (type == typeof(GeneralizedExtremeValue)) return 9;
            if (type == typeof(GeneralizedLogistic)) return 10;
            if (type == typeof(GeneralizedNormal)) return 11;
            if (type == typeof(GeneralizedPareto)) return 12;
            if (type == typeof(Exponential)) return 13;
            if (type == typeof(KappaFour)) return 14;
            if (type == typeof(Uniform)) return 15;
            if (type == typeof(Triangular)) return 16;
            if (type == typeof(Pert)) return 17;
            if (type == typeof(Deterministic)) return 18;
            if (type == typeof(Logistic)) return 19;
            if (type == typeof(Cauchy)) return 20;
            if (type == typeof(Mixture)) return 21;
            if (type == typeof(CompetingRisks)) return 22;
            return 0;
        }

        /// <summary>Returns a composite node's live child array without collection wrappers.</summary>
        /// <param name="distribution">The composite distribution.</param>
        /// <returns>The live child array, or <see langword="null"/> for a leaf or an unpopulated composite.</returns>
        private static UnivariateDistributionBase[]? ChildrenOf(UnivariateDistributionBase distribution)
        {
            if (distribution is Mixture mixture && mixture.GetType() == typeof(Mixture)) return mixture.Distributions;
            if (distribution is CompetingRisks competing && competing.GetType() == typeof(CompetingRisks)) return competing.ComponentArray;
            return null;
        }

        /// <summary>Builds a snapshot node for a tree already verified by <see cref="IsSupportedTree"/>.</summary>
        /// <param name="distribution">The distribution to capture.</param>
        /// <returns>The captured node.</returns>
        private static DistributionSnapshot Capture(UnivariateDistributionBase distribution)
        {
            byte family = FamilyOf(distribution.GetType());
            var sink = new List<long>();
            var cursor = new ScalarCursor(sink);
            VisitScalars(distribution, family, ref cursor);
            var children = ChildrenOf(distribution);
            DistributionSnapshot[]? captured = null;
            if (children is not null)
            {
                captured = new DistributionSnapshot[children.Length];
                for (int i = 0; i < children.Length; i++)
                    captured[i] = Capture(children[i]);
            }
            return new DistributionSnapshot(distribution.GetType(), family, sink.ToArray(), captured);
        }

        /// <summary>A dual-mode scalar walker that either records bits or compares them against a stored array.</summary>
        private struct ScalarCursor
        {
            /// <summary>The capture sink; <see langword="null"/> in compare mode.</summary>
            private readonly List<long>? _sink;

            /// <summary>The stored bits to compare against; <see langword="null"/> in capture mode.</summary>
            private readonly long[]? _stored;

            /// <summary>The next compare index into the stored bits.</summary>
            internal int Index;

            /// <summary>Initializes a capture-mode cursor.</summary>
            /// <param name="sink">The list receiving captured bits.</param>
            internal ScalarCursor(List<long> sink) { _sink = sink; _stored = null; Index = 0; }

            /// <summary>Initializes a compare-mode cursor.</summary>
            /// <param name="stored">The stored bits to compare against.</param>
            internal ScalarCursor(long[] stored) { _sink = null; _stored = stored; Index = 0; }

            /// <summary>Records or compares one double as raw bits.</summary>
            /// <param name="value">The live value.</param>
            /// <returns><see langword="true"/> to continue the walk; <see langword="false"/> on a compare mismatch.</returns>
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            internal bool Visit(double value) => Visit(BitConverter.DoubleToInt64Bits(value));

            /// <summary>Records or compares one integral value.</summary>
            /// <param name="value">The live value.</param>
            /// <returns><see langword="true"/> to continue the walk; <see langword="false"/> on a compare mismatch.</returns>
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            internal bool Visit(long value)
            {
                if (_sink is not null) { _sink.Add(value); return true; }
                var stored = _stored!;
                if (Index >= stored.Length || stored[Index] != value) { Index = int.MinValue; return false; }
                Index++;
                return true;
            }
        }

        /// <summary>Visits every scalar that determines a node's canonical configuration, in a fixed per-family order.</summary>
        /// <param name="distribution">The node whose scalars are visited; its exact type must already be verified.</param>
        /// <param name="family">The node's family code, stored at capture for jump-table dispatch.</param>
        /// <param name="cursor">The dual-mode cursor receiving or comparing the scalars.</param>
        /// <returns><see langword="true"/> when the walk completes; <see langword="false"/> on a compare mismatch or an unsupported type.</returns>
        /// <remarks>Leaf lists mirror each family's <c>GetParameters</c> order, extended by the
        /// evaluation-affecting settings that sit outside the flattened parameters (the logarithm
        /// base on the log families, the physical-moment surface on <see cref="LnNormal"/>, and
        /// the composite flags, weights, seed, and correlation entries). Capturing more than the
        /// canonical string records is deliberate: extra strictness can only force a string
        /// re-comparison, never a wrong match.</remarks>
        private static bool VisitScalars(UnivariateDistributionBase distribution, byte family, ref ScalarCursor cursor)
        {
            switch (family)
            {
            case 1:
            {
                var x = (Normal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma);
            }
            case 2:
            {
                var x = (LogNormal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Base);
            }
            case 3:
            {
                // The physical-moment fields are captured raw so the reported moment surface is
                // pinned without evaluating it: two distinct physical pairs can round to identical
                // log coordinates, so the log coordinates alone would under-determine the
                // canonical form.
                var x = (LnNormal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma)
                    && cursor.Visit(x.PhysicalMomentModeForSnapshot ? 1L : 0L)
                    && cursor.Visit(x.PhysicalMeanForSnapshot) && cursor.Visit(x.PhysicalStandardDeviationForSnapshot);
            }
            case 4:
            {
                var x = (LogPearsonTypeIII)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Gamma) && cursor.Visit(x.Base);
            }
            case 5:
            {
                var x = (PearsonTypeIII)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Gamma);
            }
            case 6:
            {
                var x = (GammaDistribution)distribution;
                return cursor.Visit(x.Theta) && cursor.Visit(x.Kappa);
            }
            case 7:
            {
                var x = (Weibull)distribution;
                return cursor.Visit(x.Lambda) && cursor.Visit(x.Kappa);
            }
            case 8:
            {
                var x = (Gumbel)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            case 9:
            {
                var x = (GeneralizedExtremeValue)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            case 10:
            {
                var x = (GeneralizedLogistic)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            case 11:
            {
                var x = (GeneralizedNormal)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            case 12:
            {
                var x = (GeneralizedPareto)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            case 13:
            {
                var x = (Exponential)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            case 14:
            {
                var x = (KappaFour)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa) && cursor.Visit(x.Hondo);
            }
            case 15:
            {
                var x = (Uniform)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.Max);
            }
            case 16:
            {
                var x = (Triangular)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.MostLikely) && cursor.Visit(x.Max);
            }
            case 17:
            {
                var x = (Pert)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.MostLikely) && cursor.Visit(x.Max);
            }
            case 18:
            {
                var x = (Deterministic)distribution;
                return cursor.Visit(x.Value);
            }
            case 19:
            {
                var x = (Logistic)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            case 20:
            {
                var x = (Cauchy)distribution;
                return cursor.Visit(x.X0) && cursor.Visit(x.Gamma);
            }
            case 21:
            {
                var x = (Mixture)distribution;
                var weights = x.Weights;
                if (weights is null || !cursor.Visit(x.IsZeroInflated ? 1L : 0L) || !cursor.Visit((long)x.XTransform)
                    || !cursor.Visit((long)x.ProbabilityTransform) || !cursor.Visit(x.ZeroWeight)
                    || !cursor.Visit(weights.Length)) return false;
                for (int i = 0; i < weights.Length; i++)
                    if (!cursor.Visit(weights[i])) return false;
                return true;
            }
            case 22:
            {
                var x = (CompetingRisks)distribution;
                if (!cursor.Visit(x.MinimumOfRandomVariables ? 1L : 0L) || !cursor.Visit((long)x.Dependency)
                    || !cursor.Visit(x.PRNGSeed) || !cursor.Visit((long)x.XTransform)
                    || !cursor.Visit((long)x.ProbabilityTransform)) return false;
                var matrix = x.CorrelationMatrixArray;
                if (matrix is null) return cursor.Visit(0L);
                int rows = matrix.GetLength(0), columns = matrix.GetLength(1);
                if (!cursor.Visit(1L) || !cursor.Visit(rows) || !cursor.Visit(columns)) return false;
                int rowStart = matrix.GetLowerBound(0), columnStart = matrix.GetLowerBound(1);
                for (int row = 0; row < rows; row++)
                    for (int column = 0; column < columns; column++)
                        if (!cursor.Visit(matrix[rowStart + row, columnStart + column])) return false;
                return true;
            }
            }
            return false;
        }
    }
}
