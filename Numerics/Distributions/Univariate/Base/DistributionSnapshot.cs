using System;
using System.Collections.Generic;

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

        /// <summary>The captured scalar state, as raw bits, in the fixed per-family order.</summary>
        private readonly long[] _scalarBits;

        /// <summary>Captured child snapshots for composite nodes; otherwise <see langword="null"/>.</summary>
        private readonly DistributionSnapshot[]? _children;

        /// <summary>Initializes an immutable snapshot node.</summary>
        /// <param name="type">The exact runtime type of the captured distribution.</param>
        /// <param name="scalarBits">The captured scalar state in the fixed per-family order.</param>
        /// <param name="children">Child snapshots for composite nodes, or <see langword="null"/>.</param>
        private DistributionSnapshot(Type type, long[] scalarBits, DistributionSnapshot[]? children)
        {
            _type = type;
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
        internal bool Matches(UnivariateDistributionBase? distribution)
        {
            if (distribution is null || distribution.GetType() != _type) return false;
            var cursor = new ScalarCursor(_scalarBits);
            if (!VisitScalars(distribution, ref cursor) || cursor.Index != _scalarBits.Length) return false;
            if (_children is null) return true;
            var children = ChildrenOf(distribution);
            if (children is null || children.Length != _children.Length) return false;
            for (int i = 0; i < children.Length; i++)
                if (!_children[i].Matches(children[i])) return false;
            return true;
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
            return IsSupportedLeaf(type);
        }

        /// <summary>Whether a leaf type is in the exact built-in set the scalar visitor understands.</summary>
        /// <param name="type">The exact runtime type to check.</param>
        /// <returns><see langword="true"/> for a supported leaf family; otherwise, <see langword="false"/>.</returns>
        private static bool IsSupportedLeaf(Type type)
        {
            return type == typeof(Normal) || type == typeof(LogNormal) || type == typeof(LnNormal)
                || type == typeof(LogPearsonTypeIII) || type == typeof(PearsonTypeIII)
                || type == typeof(GammaDistribution) || type == typeof(Weibull) || type == typeof(Gumbel)
                || type == typeof(GeneralizedExtremeValue) || type == typeof(GeneralizedLogistic)
                || type == typeof(GeneralizedNormal) || type == typeof(GeneralizedPareto)
                || type == typeof(Exponential) || type == typeof(KappaFour) || type == typeof(Uniform)
                || type == typeof(Triangular) || type == typeof(Pert) || type == typeof(Deterministic)
                || type == typeof(Logistic) || type == typeof(Cauchy);
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
            var sink = new List<long>();
            var cursor = new ScalarCursor(sink);
            VisitScalars(distribution, ref cursor);
            var children = ChildrenOf(distribution);
            DistributionSnapshot[]? captured = null;
            if (children is not null)
            {
                captured = new DistributionSnapshot[children.Length];
                for (int i = 0; i < children.Length; i++)
                    captured[i] = Capture(children[i]);
            }
            return new DistributionSnapshot(distribution.GetType(), sink.ToArray(), captured);
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
            internal bool Visit(double value) => Visit(BitConverter.DoubleToInt64Bits(value));

            /// <summary>Records or compares one integral value.</summary>
            /// <param name="value">The live value.</param>
            /// <returns><see langword="true"/> to continue the walk; <see langword="false"/> on a compare mismatch.</returns>
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
        /// <param name="cursor">The dual-mode cursor receiving or comparing the scalars.</param>
        /// <returns><see langword="true"/> when the walk completes; <see langword="false"/> on a compare mismatch or an unsupported type.</returns>
        /// <remarks>Leaf lists mirror each family's <c>GetParameters</c> order, extended by the
        /// evaluation-affecting settings that sit outside the flattened parameters (the logarithm
        /// base on the log families, the physical-moment surface on <see cref="LnNormal"/>, and
        /// the composite flags, weights, seed, and correlation entries). Capturing more than the
        /// canonical string records is deliberate: extra strictness can only force a string
        /// re-comparison, never a wrong match.</remarks>
        private static bool VisitScalars(UnivariateDistributionBase distribution, ref ScalarCursor cursor)
        {
            Type type = distribution.GetType();
            if (type == typeof(Normal))
            {
                var x = (Normal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma);
            }
            if (type == typeof(LogNormal))
            {
                var x = (LogNormal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Base);
            }
            if (type == typeof(LnNormal))
            {
                var x = (LnNormal)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Mean) && cursor.Visit(x.StandardDeviation);
            }
            if (type == typeof(LogPearsonTypeIII))
            {
                var x = (LogPearsonTypeIII)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Gamma) && cursor.Visit(x.Base);
            }
            if (type == typeof(PearsonTypeIII))
            {
                var x = (PearsonTypeIII)distribution;
                return cursor.Visit(x.Mu) && cursor.Visit(x.Sigma) && cursor.Visit(x.Gamma);
            }
            if (type == typeof(GammaDistribution))
            {
                var x = (GammaDistribution)distribution;
                return cursor.Visit(x.Theta) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(Weibull))
            {
                var x = (Weibull)distribution;
                return cursor.Visit(x.Lambda) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(Gumbel))
            {
                var x = (Gumbel)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            if (type == typeof(GeneralizedExtremeValue))
            {
                var x = (GeneralizedExtremeValue)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(GeneralizedLogistic))
            {
                var x = (GeneralizedLogistic)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(GeneralizedNormal))
            {
                var x = (GeneralizedNormal)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(GeneralizedPareto))
            {
                var x = (GeneralizedPareto)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa);
            }
            if (type == typeof(Exponential))
            {
                var x = (Exponential)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            if (type == typeof(KappaFour))
            {
                var x = (KappaFour)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha) && cursor.Visit(x.Kappa) && cursor.Visit(x.Hondo);
            }
            if (type == typeof(Uniform))
            {
                var x = (Uniform)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.Max);
            }
            if (type == typeof(Triangular))
            {
                var x = (Triangular)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.MostLikely) && cursor.Visit(x.Max);
            }
            if (type == typeof(Pert))
            {
                var x = (Pert)distribution;
                return cursor.Visit(x.Min) && cursor.Visit(x.MostLikely) && cursor.Visit(x.Max);
            }
            if (type == typeof(Deterministic))
            {
                var x = (Deterministic)distribution;
                return cursor.Visit(x.Value);
            }
            if (type == typeof(Logistic))
            {
                var x = (Logistic)distribution;
                return cursor.Visit(x.Xi) && cursor.Visit(x.Alpha);
            }
            if (type == typeof(Cauchy))
            {
                var x = (Cauchy)distribution;
                return cursor.Visit(x.X0) && cursor.Visit(x.Gamma);
            }
            if (type == typeof(Mixture))
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
            if (type == typeof(CompetingRisks))
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
            return false;
        }
    }
}
