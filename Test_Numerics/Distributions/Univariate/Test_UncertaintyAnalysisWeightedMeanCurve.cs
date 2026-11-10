using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for the weighted <see cref="UncertaintyAnalysisResults.ProcessMeanCurve(UnivariateDistributionBase[], double[], double[], double, double)"/>:
    /// uniform weights reproduce the unweighted mean curve bit for bit, integer weights equal the replicated
    /// ensemble, and the curve is the quantile function of the weighted mixture.
    /// </summary>
    [TestClass]
    public class Test_UncertaintyAnalysisWeightedMeanCurve
    {
        /// <summary>The non-exceedance probabilities evaluated in every case.</summary>
        private static readonly double[] Probabilities = { 0.01, 0.1, 0.5, 0.9, 0.99, 0.999 };

        /// <summary>A small GEV ensemble with distinct parameters, standing in for logic-tree branches.</summary>
        private static UnivariateDistributionBase[] Branches()
        {
            return new UnivariateDistributionBase[]
            {
                new GeneralizedExtremeValue(100d, 20d, -0.1),
                new GeneralizedExtremeValue(110d, 25d, 0d),
                new GeneralizedExtremeValue(95d, 15d, 0.1),
                new GeneralizedExtremeValue(105d, 30d, -0.2),
            };
        }

        /// <summary>Asserts two curves are bitwise identical.</summary>
        /// <param name="expected">The reference curve.</param>
        /// <param name="actual">The curve under test.</param>
        private static void AssertBitwise(double[] expected, double[] actual)
        {
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(actual[i]), $"Index {i}.");
        }

        /// <summary>Equal weights, whatever their common value, give the unweighted mean curve bit for bit.</summary>
        [TestMethod]
        public void UniformWeightsReproduceTheUnweightedMeanCurve()
        {
            var branches = Branches();
            var unweighted = new UncertaintyAnalysisResults();
            unweighted.ProcessMeanCurve(branches, Probabilities);

            foreach (double common in new[] { 1d, 0.25, 1d / 3d, 7.5 })
            {
                var weighted = new UncertaintyAnalysisResults();
                weighted.ProcessMeanCurve(branches, Probabilities, new[] { common, common, common, common });
                AssertBitwise(unweighted.MeanCurve!, weighted.MeanCurve!);
            }

            // A null member is skipped by both, whatever weight it carries.
            var withNull = new UnivariateDistributionBase[] { branches[0], null!, branches[1], branches[2] };
            var unweightedNull = new UncertaintyAnalysisResults();
            unweightedNull.ProcessMeanCurve(withNull, Probabilities);
            var weightedNull = new UncertaintyAnalysisResults();
            weightedNull.ProcessMeanCurve(withNull, Probabilities, new[] { 2d, 9d, 2d, 2d });
            AssertBitwise(unweightedNull.MeanCurve!, weightedNull.MeanCurve!);
        }

        /// <summary>
        /// Integer weights equal the unweighted curve of the ensemble with each member replicated in place:
        /// w F is exactly F + ... + F for small integers, and the chunk sums, range and grid agree.
        /// </summary>
        [TestMethod]
        public void IntegerWeightsEqualTheReplicatedEnsemble()
        {
            var a = new GeneralizedExtremeValue(100d, 20d, -0.1);
            var b = new GeneralizedExtremeValue(120d, 10d, 0.05);

            var weighted = new UncertaintyAnalysisResults();
            weighted.ProcessMeanCurve(new UnivariateDistributionBase[] { a, b }, Probabilities, new[] { 2d, 1d });
            var replicated = new UncertaintyAnalysisResults();
            replicated.ProcessMeanCurve(new UnivariateDistributionBase[] { a, a, b }, Probabilities);

            // B = 2 gives the chunks {2 F_a} and {F_b}; B = 3 gives {F_a}, {F_a} and {F_b}. The totals
            // (F_a + F_a) + F_b and 2 F_a + F_b agree exactly, and so do the divisors 3 and 2 + 1.
            AssertBitwise(replicated.MeanCurve!, weighted.MeanCurve!);
        }

        /// <summary>The weighted mean curve is the quantile function of the mixture with the normalized weights.</summary>
        [TestMethod]
        public void WeightedMeanCurveIsTheMixtureQuantile()
        {
            var branches = Branches();
            double[] weights = { 0.5, 0.2, 0.2, 0.1 };
            var results = new UncertaintyAnalysisResults();
            results.ProcessMeanCurve(branches, Probabilities, weights);
            var mixture = new Mixture(weights, branches);

            for (int i = 0; i < Probabilities.Length; i++)
            {
                double exact = mixture.InverseCDF(Probabilities[i]);
                Assert.AreEqual(exact, results.MeanCurve![i], 1E-3 * Math.Abs(exact), $"p = {Probabilities[i]}.");
            }

            // Scaling every weight leaves the mixture unchanged.
            var scaled = new UncertaintyAnalysisResults();
            scaled.ProcessMeanCurve(branches, Probabilities, new[] { 5d, 2d, 2d, 1d });
            for (int i = 0; i < Probabilities.Length; i++)
                Assert.AreEqual(results.MeanCurve![i], scaled.MeanCurve![i], 1E-12 * Math.Abs(results.MeanCurve[i]));
        }

        /// <summary>A zero weight removes the member from the sums and from the quantile range.</summary>
        [TestMethod]
        public void ZeroWeightDropsTheMember()
        {
            var branches = Branches();
            var outlier = new GeneralizedExtremeValue(1000d, 5d, 0d);
            var withOutlier = new UnivariateDistributionBase[] { branches[0], outlier, branches[1], branches[2] };

            var weighted = new UncertaintyAnalysisResults();
            weighted.ProcessMeanCurve(withOutlier, Probabilities, new[] { 0.5, 0d, 0.3, 0.2 });
            var dropped = new UncertaintyAnalysisResults();
            dropped.ProcessMeanCurve(new UnivariateDistributionBase[] { branches[0], branches[1], branches[2] }, Probabilities, new[] { 0.5, 0.3, 0.2 });

            for (int i = 0; i < Probabilities.Length; i++)
                Assert.AreEqual(dropped.MeanCurve![i], weighted.MeanCurve![i], 1E-12 * Math.Abs(dropped.MeanCurve[i]), $"p = {Probabilities[i]}.");
        }

        /// <summary>The weighted overload leaves the other results untouched and refuses invalid weights.</summary>
        [TestMethod]
        public void InvalidWeightsAreRefused()
        {
            var branches = Branches();
            var results = new UncertaintyAnalysisResults();

            Assert.ThrowsExactly<ArgumentNullException>(() => results.ProcessMeanCurve(branches, Probabilities, (double[])null!));
            Assert.ThrowsExactly<ArgumentException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { 1d, 1d }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { 1d, -1d, 1d, 1d }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { 1d, double.NaN, 1d, 1d }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { 1d, double.PositiveInfinity, 1d, 1d }));
            Assert.ThrowsExactly<ArgumentException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { 0d, 0d, 0d, 0d }));
            Assert.ThrowsExactly<ArgumentException>(() => results.ProcessMeanCurve(branches, Array.Empty<double>(), new[] { 1d, 1d, 1d, 1d }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => results.ProcessMeanCurve(branches, Probabilities, new[] { double.MaxValue, double.MaxValue, 1d, 1d }));
            Assert.ThrowsExactly<ArgumentException>(() => results.ProcessMeanCurve(null!, Probabilities, new[] { 1d }));
            Assert.ThrowsExactly<InvalidOperationException>(() => results.ProcessMeanCurve(new UnivariateDistributionBase[] { null!, null! }, Probabilities, new[] { 1d, 2d }));
            // A positive weight on a null member does not count.
            Assert.ThrowsExactly<ArgumentException>(() => results.ProcessMeanCurve(new UnivariateDistributionBase[] { branches[0], null! }, Probabilities, new[] { 0d, 2d }));
            Assert.IsNull(results.MeanCurve);
        }

        /// <summary>
        /// A single positive unit weight among zeros gives bitwise the unweighted curve of that member alone:
        /// the zero-weight members leave the range and the sums, and 1 × F / 1 is F / 1.
        /// </summary>
        [TestMethod]
        public void SinglePositiveUnitWeightEqualsTheMemberAlone()
        {
            var branches = Branches();
            var weighted = new UncertaintyAnalysisResults();
            weighted.ProcessMeanCurve(branches, Probabilities, new[] { 0d, 1d, 0d, 0d });
            var alone = new UncertaintyAnalysisResults();
            alone.ProcessMeanCurve(new[] { branches[1] }, Probabilities);
            AssertBitwise(alone.MeanCurve!, weighted.MeanCurve!);
            for (int i = 0; i < Probabilities.Length; i++)
                Assert.AreEqual(branches[1].InverseCDF(Probabilities[i]), weighted.MeanCurve![i], 1E-3 * Math.Abs(weighted.MeanCurve[i]));
        }

        /// <summary>
        /// A large non-uniform ensemble (more members than reduction chunks, with null members) is bitwise
        /// reproducible from call to call, so the result does not depend on the thread schedule; the range
        /// arguments pass through, and uniform weights with non-default ranges still equal the unweighted curve.
        /// </summary>
        [TestMethod]
        public void LargeEnsembleIsReproducibleAndRangesPassThrough()
        {
            var random = new Random(20261023);
            int count = 301;
            var members = new UnivariateDistributionBase[count];
            var weights = new double[count];
            for (int j = 0; j < count; j++)
            {
                members[j] = j % 37 == 5 ? null! : new GeneralizedExtremeValue(100 + 10 * random.NextDouble(), 15 + 10 * random.NextDouble(), 0.3 * (random.NextDouble() - 0.5));
                weights[j] = random.NextDouble();
            }
            var first = new UncertaintyAnalysisResults();
            first.ProcessMeanCurve(members, Probabilities, weights);
            for (int repeat = 0; repeat < 5; repeat++)
            {
                var again = new UncertaintyAnalysisResults();
                again.ProcessMeanCurve(members, Probabilities, weights);
                AssertBitwise(first.MeanCurve!, again.MeanCurve!);
            }
            var mixture = new Mixture(NonNullWeights(members, weights, out var nonNull), nonNull);
            for (int i = 0; i < Probabilities.Length; i++)
                Assert.AreEqual(mixture.InverseCDF(Probabilities[i]), first.MeanCurve![i], 2E-3 * Math.Abs(first.MeanCurve[i]), $"p = {Probabilities[i]}.");

            var uniform = new double[count];
            for (int j = 0; j < count; j++) uniform[j] = 0.5;
            var plain = new UncertaintyAnalysisResults();
            plain.ProcessMeanCurve(members, Probabilities, 0.01, 1 - 1E-6);
            var equal = new UncertaintyAnalysisResults();
            equal.ProcessMeanCurve(members, Probabilities, uniform, 0.01, 1 - 1E-6);
            AssertBitwise(plain.MeanCurve!, equal.MeanCurve!);
        }

        /// <summary>The non-null members and their normalized weights, for the mixture reference.</summary>
        /// <param name="members">The ensemble, with null members.</param>
        /// <param name="weights">The aligned weights.</param>
        /// <param name="nonNull">On return, the non-null members.</param>
        /// <returns>The normalized weights of the non-null members.</returns>
        private static double[] NonNullWeights(UnivariateDistributionBase[] members, double[] weights, out UnivariateDistributionBase[] nonNull)
        {
            var kept = new System.Collections.Generic.List<UnivariateDistributionBase>();
            var w = new System.Collections.Generic.List<double>();
            double total = 0;
            for (int j = 0; j < members.Length; j++)
            {
                if (members[j] is null) continue;
                kept.Add(members[j]);
                w.Add(weights[j]);
                total += weights[j];
            }
            nonNull = kept.ToArray();
            var normalized = w.ToArray();
            for (int j = 0; j < normalized.Length; j++) normalized[j] /= total;
            return normalized;
        }
    }
}
