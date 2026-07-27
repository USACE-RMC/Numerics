using System;
using System.Collections.Generic;
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
    }
}
