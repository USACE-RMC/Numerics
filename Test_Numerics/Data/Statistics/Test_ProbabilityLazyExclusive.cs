using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Mathematics.SpecialFunctions;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for lazy exclusive-probability enumeration, including dense-order agreement,
    /// convergence, residual mass, row reuse, and dimensions beyond dense enumeration.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_ProbabilityLazyExclusive
    {
        /// <summary>Builds the binomial combination counts for n events.</summary>
        private static int[] BinomialCounts(int n)
        {
            var counts = new int[n];
            for (int i = 1; i <= n; i++) counts[i - 1] = (int)Math.Round(Factorial.BinomialCoefficient(n, i));
            return counts;
        }

        /// <summary>
        /// Verifies the lazy generator reproduces the dense combination matrix row for row.
        /// </summary>
        [TestMethod]
        public void Test_AllCombinationsLazy_MatchesDenseRowOrder()
        {
            for (int n = 1; n <= 12; n++)
            {
                var dense = Factorial.AllCombinations(n);
                int row = 0;
                foreach (var combination in Factorial.AllCombinationsLazy(n))
                {
                    var expanded = new int[n];
                    for (int t = 0; t < combination.Length; t++) expanded[combination[t]] = 1;
                    for (int column = 0; column < n; column++)
                    {
                        Assert.AreEqual(dense[row, column], expanded[column], $"n={n}, row {row}, column {column}");
                    }
                    row++;
                }
                Assert.AreEqual(dense.GetLength(0), row, $"n={n} row count");
            }
        }

        /// <summary>
        /// Verifies the lazy expansion emits the same rows and probabilities as the dense overload
        /// when the expansion runs to completion.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_MatchesDense_WhenComplete()
        {
            // Probabilities small enough that the bracket cannot converge before the last size.
            var probabilities = new double[] { 0.4d, 0.35d, 0.3d, 0.25d, 0.2d };
            int n = probabilities.Length;

            var denseProbabilities = new List<double>();
            var denseIndicators = new List<int[]>();
            bool truncated = Probability.IndependentExclusive(probabilities, BinomialCounts(n),
                Factorial.AllCombinations(n), denseProbabilities, denseIndicators, 0d, 0d);
            Assert.IsFalse(truncated);

            var lazyProbabilities = new List<double>();
            var lazyIndicators = new List<int[]>();
            var status = Probability.IndependentExclusiveLazy(probabilities, lazyProbabilities, lazyIndicators,
                includeNoEventRow: false, maxEmittedCombinations: 0, absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Complete, status);
            Assert.HasCount(denseProbabilities.Count, lazyProbabilities);
            for (int i = 0; i < denseProbabilities.Count; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseProbabilities[i]),
                    BitConverter.DoubleToInt64Bits(lazyProbabilities[i]), $"probability at row {i}");
                CollectionAssert.AreEqual(denseIndicators[i], lazyIndicators[i], $"indicators at row {i}");
            }
        }

        /// <summary>
        /// Verifies the exclusive probabilities sum to one over a complete expansion, with the
        /// no-event row included.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_WithNoEventRow_SumsToOne()
        {
            var probabilities = new double[] { 0.4d, 0.35d, 0.3d, 0.25d };
            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();

            var status = Probability.IndependentExclusiveLazy(probabilities, eventProbabilities, eventIndicators,
                includeNoEventRow: true, maxEmittedCombinations: 0, absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Complete, status);
            Assert.HasCount(1 << probabilities.Length, eventProbabilities);
            double total = 0d;
            for (int i = 0; i < eventProbabilities.Count; i++) total += eventProbabilities[i];
            Assert.AreEqual(1d, total, 1e-12);
        }

        /// <summary>
        /// Verifies the cap stops the expansion and closes with the exact unenumerated mass.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_Capped_ClosesWithResidualMass()
        {
            var probabilities = new double[] { 0.4d, 0.35d, 0.3d, 0.25d, 0.2d, 0.15d };
            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();

            var status = Probability.IndependentExclusiveLazy(probabilities, eventProbabilities, eventIndicators,
                includeNoEventRow: true, maxEmittedCombinations: 10, absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Capped, status);

            // Ten enumerated combinations, plus the no-event row, plus the closing row.
            Assert.HasCount(12, eventProbabilities);
            double total = 0d;
            for (int i = 0; i < eventProbabilities.Count; i++) total += eventProbabilities[i];
            Assert.AreEqual(1d, total, 1e-12, "The closing row must carry exactly the unenumerated mass.");

            // The residual is attributed to the all-ones combination.
            var closing = eventIndicators[eventIndicators.Count - 1];
            for (int i = 0; i < closing.Length; i++) Assert.AreEqual(1, closing[i]);
        }

        /// <summary>
        /// Verifies the expansion runs beyond the dimension at which the dense matrix cannot be
        /// built, converging on realistic failure probabilities rather than enumerating 2^n rows.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_BeyondDenseDimensionLimit()
        {
            // The dense form throws above thirty events; nothing here allocates per row.
            const int n = 40;
            Assert.Throws<ArgumentOutOfRangeException>(() => Factorial.AllCombinations(n));

            // Failure probabilities of the order a risk model carries, where the bracket closes
            // after the triples.
            var probabilities = new double[n];
            for (int i = 0; i < n; i++) probabilities[i] = 5e-4d + i * 1e-5d;

            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();
            var status = Probability.IndependentExclusiveLazy(probabilities, eventProbabilities, eventIndicators,
                includeNoEventRow: true, maxEmittedCombinations: 1_000_000);

            // Converged on the bracket rather than stopping at the cap, having enumerated a
            // vanishing fraction of the 2^40 combinations the dense form would have needed.
            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, status);
            double enumeratedFraction = eventProbabilities.Count / Math.Pow(2d, n);
            Assert.IsLessThanOrEqualTo(1e-5, enumeratedFraction,
                $"Emitted {eventProbabilities.Count} of 2^{n} combinations.");
        }

        /// <summary>
        /// Verifies indicator row arrays are reused across calls, so a repeated caller allocates
        /// no output rows after the first call.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_ReusesIndicatorRows()
        {
            var probabilities = new double[] { 0.3d, 0.25d, 0.2d };
            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();

            Probability.IndependentExclusiveLazy(probabilities, eventProbabilities, eventIndicators,
                absoluteTolerance: 0d, relativeTolerance: 0d);
            var first = new int[eventIndicators.Count][];
            eventIndicators.CopyTo(first);

            Probability.IndependentExclusiveLazy(probabilities, eventProbabilities, eventIndicators,
                absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.HasCount(first.Length, eventIndicators);
            for (int i = 0; i < first.Length; i++)
            {
                Assert.AreSame(first[i], eventIndicators[i], $"Row {i} was reallocated.");
            }
        }

        /// <summary>Builds a positive-definite equicorrelation matrix.</summary>
        private static double[,] CorrelationMatrix(int n, double correlation)
        {
            var matrix = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) matrix[i, j] = i == j ? 1d : correlation;
            }

            return matrix;
        }

        /// <summary>
        /// Verifies complete lazy PCM union and exclusive enumeration are bit-identical to the
        /// established dense calculations, including indicator order.
        /// </summary>
        [TestMethod]
        public void Test_LazyPCM_MatchesDense_WhenComplete()
        {
            var probabilities = new double[] { 0.32d, 0.27d, 0.19d, 0.11d };
            int n = probabilities.Length;
            int[] counts = BinomialCounts(n);
            int[,] indicators = Factorial.AllCombinations(n);
            double[,] correlation = CorrelationMatrix(n, 0.25d);

            double denseUnion = Probability.UnionPCM(
                probabilities, counts, indicators, correlation, 0d, 0d);
            double lazyUnion = Probability.UnionPCMLazy(
                probabilities, correlation, out var unionStatus, 0d, 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Complete, unionStatus);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseUnion),
                BitConverter.DoubleToInt64Bits(lazyUnion));

            Probability.ExclusivePCM(probabilities, counts, indicators, correlation,
                out List<double> denseProbabilities, out List<int[]> denseIndicators, 0d, 0d);
            var lazyProbabilities = new List<double>();
            var lazyIndicators = new List<int[]>();
            var exclusiveStatus = Probability.ExclusivePCMLazy(
                probabilities, correlation, lazyProbabilities, lazyIndicators, 0d, 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Complete, exclusiveStatus);
            Assert.HasCount(denseProbabilities.Count, lazyProbabilities);
            for (int i = 0; i < denseProbabilities.Count; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseProbabilities[i]),
                    BitConverter.DoubleToInt64Bits(lazyProbabilities[i]), $"probability at row {i}");
                CollectionAssert.AreEqual(denseIndicators[i], lazyIndicators[i], $"indicators at row {i}");
            }
        }

        /// <summary>
        /// Verifies early-converged lazy PCM preserves the dense row order, half-gap row,
        /// exclusive values, union value, and completion status.
        /// </summary>
        [TestMethod]
        public void Test_LazyPCM_MatchesDense_WhenConverged()
        {
            var probabilities = new double[] { 0.005d, 0.004d, 0.003d, 0.002d, 0.001d, 0.0005d };
            int n = probabilities.Length;
            int[] counts = BinomialCounts(n);
            int[,] indicators = Factorial.AllCombinations(n);
            double[,] correlation = CorrelationMatrix(n, 0.1d);

            double denseUnion = Probability.UnionPCM(
                probabilities, counts, indicators, correlation, 1E-4, 1E-4);
            double lazyUnion = Probability.UnionPCMLazy(
                probabilities, correlation, out var unionStatus);

            Probability.ExclusivePCM(probabilities, counts, indicators, correlation,
                out List<double> denseProbabilities, out List<int[]> denseIndicators);
            var lazyProbabilities = new List<double>();
            var lazyIndicators = new List<int[]>();
            var exclusiveStatus = Probability.ExclusivePCMLazy(
                probabilities, correlation, lazyProbabilities, lazyIndicators);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, unionStatus);
            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, exclusiveStatus);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseUnion),
                BitConverter.DoubleToInt64Bits(lazyUnion));
            Assert.HasCount(denseProbabilities.Count, lazyProbabilities);
            for (int i = 0; i < denseProbabilities.Count; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseProbabilities[i]),
                    BitConverter.DoubleToInt64Bits(lazyProbabilities[i]), $"probability at row {i}");
                CollectionAssert.AreEqual(denseIndicators[i], lazyIndicators[i], $"indicators at row {i}");
            }
        }

        /// <summary>
        /// Verifies lazy perfectly-positive enumeration preserves the dense algorithm and reuses
        /// caller-owned indicator rows.
        /// </summary>
        [TestMethod]
        public void Test_PositiveLazy_MatchesDenseAndReusesRows()
        {
            var probabilities = new double[] { 0.4d, 0.3d, 0.2d, 0.1d };
            int n = probabilities.Length;
            Probability.PositivelyDependentExclusive(probabilities, BinomialCounts(n),
                Factorial.AllCombinations(n), out List<double> denseProbabilities,
                out List<int[]> denseIndicators, 0d, 0d);

            var lazyProbabilities = new List<double>();
            var lazyIndicators = new List<int[]>();
            var status = Probability.PositivelyDependentExclusiveLazy(
                probabilities, lazyProbabilities, lazyIndicators, 0d, 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Complete, status);
            Assert.HasCount(denseProbabilities.Count, lazyProbabilities);
            for (int i = 0; i < denseProbabilities.Count; i++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(denseProbabilities[i]),
                    BitConverter.DoubleToInt64Bits(lazyProbabilities[i]), $"probability at row {i}");
                CollectionAssert.AreEqual(denseIndicators[i], lazyIndicators[i], $"indicators at row {i}");
            }

            int[][] firstRows = lazyIndicators.ToArray();
            Probability.PositivelyDependentExclusiveLazy(
                probabilities, lazyProbabilities, lazyIndicators, 0d, 0d);
            for (int i = 0; i < firstRows.Length; i++)
            {
                Assert.AreSame(firstRows[i], lazyIndicators[i], $"Row {i} was reallocated.");
            }
        }

        /// <summary>
        /// Verifies dependent lazy enumeration supports more than twenty events, while the
        /// independent lazy path separately operates above the dense matrix limit.
        /// </summary>
        [TestMethod]
        public void Test_DependentLazyEnumeration_ExceedsTwentyEvents()
        {
            const int dimension = 24;
            var probabilities = new double[dimension];
            probabilities[0] = 7E-4d;
            probabilities[1] = 6E-4d;
            probabilities[2] = 5E-4d;

            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();
            var status = Probability.PositivelyDependentExclusiveLazy(
                probabilities, eventProbabilities, eventIndicators);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, status);
            Assert.IsGreaterThan(0, eventProbabilities.Count);
            Assert.IsLessThan(Math.Pow(2d, dimension), eventProbabilities.Count);
        }

        /// <summary>
        /// Verifies dense convenience methods reject dimensions whose combination count exceeds
        /// the signed 32-bit array limit before attempting exponential allocation or recursion.
        /// </summary>
        [TestMethod]
        public void Test_DenseExclusiveConvenienceMethods_RejectInt32OverflowDimension()
        {
            var probabilities = new double[31];

            Assert.Throws<ArgumentOutOfRangeException>(
                () => Probability.IndependentExclusive(probabilities));
            Assert.Throws<ArgumentOutOfRangeException>(
                () => Probability.PositivelyDependentExclusive(probabilities));
        }

        /// <summary>
        /// Verifies PCM enumeration remains convergence-driven beyond the dense matrix limit and
        /// returns only finite probabilities without imposing a combination cap.
        /// </summary>
        [TestMethod]
        public void Test_LazyPCM_BeyondDenseDimensionLimit_ConvergesToProbabilities()
        {
            const int dimension = 32;
            var probabilities = new double[dimension];
            for (int i = 0; i < dimension; i++) probabilities[i] = 5E-6d + i * 1E-8d;
            double[,] correlation = CorrelationMatrix(dimension, 0.1d);

            Assert.Throws<ArgumentOutOfRangeException>(
                () => Factorial.AllCombinations(dimension));

            double union = Probability.UnionPCMLazy(
                probabilities, correlation, out var unionStatus);
            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, unionStatus);
            Assert.IsTrue(!double.IsNaN(union) && !double.IsInfinity(union) &&
                union >= 0d && union <= 1d, $"UnionPCM returned {union:R}.");

            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();
            var exclusiveStatus = Probability.ExclusivePCMLazy(
                probabilities, correlation, eventProbabilities, eventIndicators);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Converged, exclusiveStatus);
            Assert.HasCount(eventProbabilities.Count, eventIndicators);
            Assert.IsLessThan((double)int.MaxValue, eventProbabilities.Count);
            foreach (double probability in eventProbabilities)
            {
                Assert.IsTrue(!double.IsNaN(probability) && !double.IsInfinity(probability) &&
                    probability >= 0d && probability <= 1d,
                    $"ExclusivePCM returned {probability:R}.");
            }
        }

        /// <summary>
        /// Pins every new lazy API's default convergence tolerances at 1E-4.
        /// </summary>
        [TestMethod]
        public void Test_LazyProbability_DefaultTolerancesRemainOneEminusFour()
        {
            string[] methodNames =
            {
                nameof(Probability.UnionPCMLazy),
                nameof(Probability.ExclusivePCMLazy),
                nameof(Probability.PositivelyDependentExclusiveLazy),
            };

            foreach (string methodName in methodNames)
            {
                var method = typeof(Probability).GetMethod(methodName);
                Assert.IsNotNull(method, methodName);
                var parameters = method.GetParameters();
                Assert.AreEqual(1E-4, parameters[parameters.Length - 2].DefaultValue, methodName);
                Assert.AreEqual(1E-4, parameters[parameters.Length - 1].DefaultValue, methodName);
            }
        }

        /// <summary>
        /// Verifies the lazy PCM output rows are reused across calls.
        /// </summary>
        [TestMethod]
        public void Test_LazyPCM_ReusesIndicatorRows()
        {
            var probabilities = new double[] { 0.2d, 0.15d, 0.1d, 0.05d };
            var eventProbabilities = new List<double>();
            var eventIndicators = new List<int[]>();
            double[,] correlation = CorrelationMatrix(probabilities.Length, 0.2d);

            Probability.ExclusivePCMLazy(
                probabilities, correlation, eventProbabilities, eventIndicators, 0d, 0d);
            int[][] firstRows = eventIndicators.ToArray();

            Probability.ExclusivePCMLazy(
                probabilities, correlation, eventProbabilities, eventIndicators, 0d, 0d);

            Assert.HasCount(firstRows.Length, eventIndicators);
            for (int i = 0; i < firstRows.Length; i++)
            {
                Assert.AreSame(firstRows[i], eventIndicators[i], $"Row {i} was reallocated.");
            }
        }

        /// <summary>
        /// Verifies finite scalar and array probability outputs are clipped while NaN signaling is
        /// retained, and proves clipping is cellwise rather than proportional normalization.
        /// </summary>
        [TestMethod]
        public void Test_ProbabilityOutputs_AreClippedWithoutNormalization()
        {
            static void AssertProbability(double value, string label)
            {
                Assert.IsTrue(double.IsNaN(value) || value >= 0d && value <= 1d,
                    $"{label} returned {value:R}.");
            }

            AssertProbability(Probability.AAndB(1d, 2d), nameof(Probability.AAndB));
            AssertProbability(Probability.IndependentJointProbability(new[] { 2d, 2d }),
                nameof(Probability.IndependentJointProbability));
            AssertProbability(Probability.PositiveJointProbability(new[] { 2d, 3d }),
                nameof(Probability.PositiveJointProbability));
            AssertProbability(Probability.NegativeJointProbability(new[] { 2d, 3d }),
                nameof(Probability.NegativeJointProbability));
            AssertProbability(Probability.IndependentUnion(new[] { -1d, -1d }),
                nameof(Probability.IndependentUnion));
            AssertProbability(Probability.PositivelyDependentUnion(new[] { 2d, 3d }),
                nameof(Probability.PositivelyDependentUnion));
            AssertProbability(Probability.NegativelyDependentUnion(new[] { 2d, 3d }),
                nameof(Probability.NegativelyDependentUnion));

            var adversarial = new[] { 2d, 0.5d };
            int[,] indicators = Factorial.AllCombinations(adversarial.Length);
            double[] independent = Probability.IndependentExclusive(adversarial, indicators);
            double[] positive = Probability.PositivelyDependentExclusive(adversarial, indicators);
            foreach (double probability in independent) AssertProbability(probability, "IndependentExclusive");
            foreach (double probability in positive) AssertProbability(probability, "PositivelyDependentExclusive");

            Assert.AreEqual(1d, independent[0], 0d);
            Assert.AreEqual(0d, independent[1], 0d);
            Assert.AreEqual(1d, independent[2], 0d);
            Assert.AreEqual(2d, independent[0] + independent[1] + independent[2], 0d,
                "Cellwise clipping must not proportionally renormalize an exclusive partition.");

            double nan = Probability.IndependentExclusive(
                new[] { double.NaN, 0.5d }, new[] { 1, 0 });
            Assert.IsTrue(double.IsNaN(nan), "Tools.Clamp must preserve the established NaN signal.");
        }

        /// <summary>
        /// Test that a capped enumeration without the no-event row still closes on the exact
        /// union mass: the emitted rows plus the closing pseudo-row sum to one minus the
        /// no-event probability.
        /// </summary>
        [TestMethod]
        public void Test_LazyExclusive_CappedWithoutNoEventRow_PreservesUnionMass()
        {
            double[] probabilities = { 0.4d, 0.35d, 0.3d, 0.25d, 0.2d, 0.15d };
            var output = new List<double>();
            var indicators = new List<int[]>();

            var status = Probability.IndependentExclusiveLazy(probabilities, output, indicators,
                includeNoEventRow: false, maxEmittedCombinations: 10,
                absoluteTolerance: 0d, relativeTolerance: 0d);

            Assert.AreEqual(Probability.ExclusiveEnumerationStatus.Capped, status);
            double noEventMass = probabilities.Aggregate(1d, (mass, probability) => mass * (1d - probability));
            Assert.AreEqual(1d - noEventMass, output.Sum(), 1E-12);
            Assert.HasCount(11, output);
        }
    }
}
