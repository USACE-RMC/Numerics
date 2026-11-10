using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for <see cref="UnivariateDistributionBase.ParameterSetQuantiles(double[,], IList{double}, out bool[])"/>:
    /// the quantile matrix must be bitwise identical to the per-row loop it replaces.
    /// </summary>
    [TestClass]
    public class Test_ParameterSetQuantiles
    {
        /// <summary>The probabilities evaluated in every case, including both ends of the support.</summary>
        private static readonly double[] Probabilities = { 0.0, 1E-8, 0.001, 0.05, 0.5, 0.95, 0.99, 1.0 - 1E-9, 1.0 };

        /// <summary>
        /// Compares the batched member against the per-row loop ValidateParameters → SetParameters →
        /// ParametersValid → InverseCDF on a clone, bit for bit, and checks the instance is unchanged.
        /// </summary>
        /// <param name="distribution">The distribution whose family is evaluated.</param>
        /// <param name="parameterSets">The parameter rows, valid and invalid.</param>
        /// <param name="expectedInvalid">The number of rows expected to be invalid.</param>
        private static void AssertMatchesPerRowLoop(UnivariateDistributionBase distribution, double[,] parameterSets, int expectedInvalid)
        {
            var before = distribution.GetParameters;

            var quantiles = distribution.ParameterSetQuantiles(parameterSets, Probabilities, out bool[] valid);

            int rows = parameterSets.GetLength(0), arity = parameterSets.GetLength(1);
            Assert.AreEqual(rows, quantiles.GetLength(0));
            Assert.AreEqual(Probabilities.Length, quantiles.GetLength(1));
            Assert.HasCount(rows, valid);
            var probe = distribution.Clone();
            int invalid = 0;
            for (int i = 0; i < rows; i++)
            {
                var row = new double[arity];
                for (int j = 0; j < arity; j++) row[j] = parameterSets[i, j];
                bool expectedValid = probe.ValidateParameters(row, false) is null;
                Assert.AreEqual(expectedValid, valid[i], $"row {i}");
                probe.SetParameters(row);
                Assert.AreEqual(probe.ParametersValid, valid[i], $"row {i}: SetParameters validity");
                if (!expectedValid)
                {
                    invalid++;
                    for (int k = 0; k < Probabilities.Length; k++)
                        Assert.IsTrue(double.IsNaN(quantiles[i, k]), $"row {i}, p {Probabilities[k]}");
                    continue;
                }
                for (int k = 0; k < Probabilities.Length; k++)
                {
                    double expected = probe.InverseCDF(Probabilities[k]);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(quantiles[i, k]),
                        $"row {i}, p {Probabilities[k]}: expected {expected:R}, actual {quantiles[i, k]:R}");
                }
            }
            Assert.AreEqual(expectedInvalid, invalid);
            CollectionAssert.AreEqual(before, distribution.GetParameters);
        }

        /// <summary>GEV rows on both sides of κ = 0, at κ = 0, and invalid scales.</summary>
        [TestMethod]
        public void GeneralizedExtremeValueMatchesPerRowLoop()
        {
            var parameterSets = new double[,]
            {
                { 4.8294, 1.2409, -0.1136 },
                { 4.5, 1.0, 0.25 },
                { 5.1, 1.4, 0.0 },
                { 4.9, 1.3, -0.45 },
                { 4.9, -1.3, -0.1 },
                { double.NaN, 1.3, -0.1 },
            };
            AssertMatchesPerRowLoop(new GeneralizedExtremeValue(4.8294, 1.2409, -0.1136), parameterSets, 2);
        }

        /// <summary>
        /// The GEV override (the probability transform computed once for all rows) stays bitwise identical to
        /// the per-row loop over many random rows, extreme and tiny shapes, and probabilities at and next to the
        /// ends of the support.
        /// </summary>
        [TestMethod]
        public void GeneralizedExtremeValueOverrideMatchesPerRowLoopOverManyRows()
        {
            var random = new Random(20261023);
            int rows = 2000;
            var parameterSets = new double[rows, 3];
            for (int i = 0; i < rows; i++)
            {
                parameterSets[i, 0] = 100 + 50 * (random.NextDouble() - 0.5);
                parameterSets[i, 1] = 20 * Math.Exp(random.NextDouble() - 0.5);
                parameterSets[i, 2] = 0.8 * (random.NextDouble() - 0.5);
            }
            double[] shapes = { 0.0, -0.0, 1E-300, -1E-300, 1E-17, -1E-17, 5.0, -5.0, 50.0, -50.0 };
            for (int s = 0; s < shapes.Length; s++) parameterSets[s, 2] = shapes[s];
            parameterSets[20, 0] = 1E300;
            parameterSets[21, 1] = 1E300;
            parameterSets[22, 1] = 0.0;
            parameterSets[23, 2] = double.PositiveInfinity;

            var probabilities = new[] { 0.0, double.Epsilon, 1E-300, 1E-12, 0.01, 0.5, 0.99, 1.0 - 1E-12, 1.0 - 1E-16, 1.0 };
            var distribution = new GeneralizedExtremeValue(100, 20, -0.1);
            var quantiles = distribution.ParameterSetQuantiles(parameterSets, probabilities, out bool[] valid);

            var probe = new GeneralizedExtremeValue();
            for (int i = 0; i < rows; i++)
            {
                var row = new[] { parameterSets[i, 0], parameterSets[i, 1], parameterSets[i, 2] };
                Assert.AreEqual(probe.ValidateParameters(row, false) is null, valid[i], $"row {i}");
                if (!valid[i])
                {
                    for (int k = 0; k < probabilities.Length; k++) Assert.IsTrue(double.IsNaN(quantiles[i, k]));
                    continue;
                }
                probe.SetParameters(row);
                for (int k = 0; k < probabilities.Length; k++)
                {
                    double expected = probe.InverseCDF(probabilities[k]);
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(quantiles[i, k]),
                        $"row {i}, p {probabilities[k]:R}: expected {expected:R}, actual {quantiles[i, k]:R}");
                }
            }
            Assert.AreEqual(2, valid.Count(v => !v));

            // An invalid probability throws as the per-row InverseCDF does, and only when a row is valid.
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.ParameterSetQuantiles(parameterSets, new[] { 0.5, 1.5 }, out _));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.ParameterSetQuantiles(parameterSets, new[] { double.NaN }, out _));
            var none = distribution.ParameterSetQuantiles(new double[,] { { 1.0, -1.0, 0.0 } }, new[] { 1.5 }, out bool[] noneValid);
            Assert.IsFalse(noneValid[0]);
            Assert.IsTrue(double.IsNaN(none[0, 0]));
        }

        /// <summary>
        /// The GEV override handles empty inputs as the base member does, never modifies the calling instance
        /// (even one whose own parameters are invalid), and leaves the scalar InverseCDF on its known values.
        /// </summary>
        [TestMethod]
        public void GeneralizedExtremeValueOverrideEdgeCases()
        {
            var distribution = new GeneralizedExtremeValue(4.8294, 1.2409, -0.1136);
            var noRows = distribution.ParameterSetQuantiles(new double[0, 3], Probabilities, out bool[] noRowsValid);
            Assert.AreEqual(0, noRows.GetLength(0));
            Assert.AreEqual(Probabilities.Length, noRows.GetLength(1));
            Assert.IsEmpty(noRowsValid);
            var noProbabilities = distribution.ParameterSetQuantiles(new double[,] { { 4.0, 1.0, 0.1 }, { 4.0, 0.0, 0.1 } }, new double[0], out bool[] valid);
            Assert.AreEqual(0, noProbabilities.GetLength(1));
            CollectionAssert.AreEqual(new[] { true, false }, valid);
            CollectionAssert.AreEqual(new[] { 4.8294, 1.2409, -0.1136 }, distribution.GetParameters);

            var invalidSelf = new GeneralizedExtremeValue(4.0, -1.0, 0.0);
            var fromInvalid = invalidSelf.ParameterSetQuantiles(new double[,] { { 4.0, 1.0, 0.0 } }, new[] { 0.5 }, out bool[] fromInvalidValid);
            Assert.IsTrue(fromInvalidValid[0]);
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(new GeneralizedExtremeValue(4.0, 1.0, 0.0).InverseCDF(0.5)), BitConverter.DoubleToInt64Bits(fromInvalid[0, 0]));
            Assert.IsFalse(invalidSelf.ParametersValid);

            // Gumbel limit and a closed form: at κ = 0 the quantile is ξ − α ln(−ln p).
            var gumbel = new GeneralizedExtremeValue(10.0, 2.0, 0.0);
            Assert.AreEqual(10.0 - 2.0 * Math.Log(-Math.Log(0.9)), gumbel.InverseCDF(0.9), 1E-12);
            // κ ≠ 0: ξ + α (1 − (−ln p)^κ) / κ.
            var gev = new GeneralizedExtremeValue(10.0, 2.0, 0.2);
            Assert.AreEqual(10.0 + 2.0 * (1.0 - Math.Pow(-Math.Log(0.9), 0.2)) / 0.2, gev.InverseCDF(0.9), 1E-12);
        }

        /// <summary>Log-Pearson Type III rows with the clone's non-default base carried to every row.</summary>
        [TestMethod]
        public void LogPearsonTypeIIIMatchesPerRowLoopAndCarriesBase()
        {
            var parameterSets = new double[,]
            {
                { 3.0, 0.25, 0.4 },
                { 3.1, 0.3, -0.6 },
                { 2.9, 0.2, 0.0 },
                { 3.0, 0.0, 0.4 },
            };
            var distribution = new LogPearsonTypeIII(3.0, 0.25, 0.4) { Base = Math.E };
            AssertMatchesPerRowLoop(distribution, parameterSets, 1);

            // The base is not a parameter: a base-e row differs from the base-10 row.
            var baseTen = new LogPearsonTypeIII(3.0, 0.25, 0.4);
            var natural = distribution.ParameterSetQuantiles(parameterSets, new[] { 0.5 }, out _);
            var common = baseTen.ParameterSetQuantiles(parameterSets, new[] { 0.5 }, out _);
            Assert.AreNotEqual(common[0, 0], natural[0, 0]);
            Assert.AreEqual(distribution.InverseCDF(0.5), natural[0, 0]);
        }

        /// <summary>Normal rows, including a zero standard deviation.</summary>
        [TestMethod]
        public void NormalMatchesPerRowLoop()
        {
            var parameterSets = new double[,]
            {
                { 100.0, 15.0 },
                { -3.0, 0.5 },
                { 0.0, 1.0 },
                { 1.0, -2.0 },
            };
            AssertMatchesPerRowLoop(new Normal(100.0, 15.0), parameterSets, 1);
        }

        /// <summary>Gumbel rows, including an invalid scale.</summary>
        [TestMethod]
        public void GumbelMatchesPerRowLoop()
        {
            var parameterSets = new double[,]
            {
                { 10.0, 2.0 },
                { 12.5, 3.25 },
                { 10.0, 0.0 },
            };
            AssertMatchesPerRowLoop(new Gumbel(10.0, 2.0), parameterSets, 1);
        }

        /// <summary>Kappa-four rows across both shape signs, and an invalid scale.</summary>
        [TestMethod]
        public void KappaFourMatchesPerRowLoop()
        {
            var parameterSets = new double[,]
            {
                { 100.0, 10.0, 0.1, -0.2 },
                { 100.0, 10.0, -0.1, 0.3 },
                { 90.0, 12.0, 0.05, 0.0 },
                { 100.0, -10.0, 0.1, 0.2 },
            };
            AssertMatchesPerRowLoop(new KappaFour(100.0, 10.0, 0.1, -0.2), parameterSets, 1);
        }

        /// <summary>Zero rows return empty results; zero probabilities still report validity.</summary>
        [TestMethod]
        public void EmptyInputsReturnEmptyResults()
        {
            var distribution = new Normal(0.0, 1.0);

            var none = distribution.ParameterSetQuantiles(new double[0, 2], Probabilities, out bool[] noneValid);
            var noProbabilities = distribution.ParameterSetQuantiles(new double[,] { { 0.0, 1.0 }, { 0.0, -1.0 } }, new double[0], out bool[] valid);

            Assert.AreEqual(0, none.GetLength(0));
            Assert.IsEmpty(noneValid);
            Assert.AreEqual(0, noProbabilities.GetLength(1));
            CollectionAssert.AreEqual(new[] { true, false }, valid);
        }

        /// <summary>Null arguments and a wrong column count are rejected.</summary>
        [TestMethod]
        public void InvalidArgumentsThrow()
        {
            var distribution = new GeneralizedExtremeValue(4.0, 1.0, -0.1);

            Assert.ThrowsExactly<ArgumentNullException>(() => distribution.ParameterSetQuantiles(null!, Probabilities, out _));
            Assert.ThrowsExactly<ArgumentNullException>(() => distribution.ParameterSetQuantiles(new double[1, 3], null!, out _));
            Assert.ThrowsExactly<ArgumentException>(() => distribution.ParameterSetQuantiles(new double[1, 2], Probabilities, out _));
        }

        /// <summary>
        /// Table-backed distributions, which do not support setting parameters from a vector, report that when a
        /// row is valid; their instances are not modified.
        /// </summary>
        [TestMethod]
        public void TableBackedDistributionsDoNotSupportParameterRows()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.1d, 0.5d, 0.9d });
            var kernel = new KernelDensity(new[] { 1.2d, 3.4d, 2.2d, 5.6d, 4.4d });

            Assert.ThrowsExactly<NotImplementedException>(() => empirical.ParameterSetQuantiles(new double[,] { { 0.0, 1.0 } }, new[] { 0.5 }, out _));
            Assert.ThrowsExactly<NotImplementedException>(() => kernel.ParameterSetQuantiles(new double[,] { { 0.0, 1.0, 2.0 } }, new[] { 0.5 }, out _));
            Assert.AreEqual(2d, empirical.InverseCDF(0.5d), 1E-12);
        }
    }
}
