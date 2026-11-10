using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the Multinomial distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Reference values verified analytically and against R dmultinom.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Multinomial
    {

        /// <summary>
        /// Test basic construction.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            Assert.AreEqual(3, M.Dimension);
            Assert.AreEqual(10, M.NumberOfTrials);
            Assert.IsTrue(M.ParametersValid);
        }

        /// <summary>
        /// Test invalid parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(0, new double[] { 0.5, 0.5 }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { 0.5 })); // too few categories
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { -0.1, 0.6, 0.5 })); // negative prob
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(10, new double[] { 0.3, 0.3, 0.3 })); // don't sum to 1
        }

        /// <summary>
        /// Test type and display names.
        /// </summary>
        [TestMethod]
        public void Test_TypeAndName()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            Assert.AreEqual(MultivariateDistributionType.Multinomial, M.Type);
            Assert.AreEqual("Multinomial", M.DisplayName);
            Assert.AreEqual("Mult", M.ShortDisplayName);
        }

        /// <summary>
        /// Test mean vector. Mean[i] = N * p[i].
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var mean = M.Mean;
            Assert.AreEqual(2.0, mean[0], 1e-10);
            Assert.AreEqual(3.0, mean[1], 1e-10);
            Assert.AreEqual(5.0, mean[2], 1e-10);
        }

        /// <summary>
        /// Test variance vector. Var[i] = N * p[i] * (1 - p[i]).
        /// </summary>
        [TestMethod]
        public void Test_Variance()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var variance = M.Variance;
            Assert.AreEqual(10 * 0.2 * 0.8, variance[0], 1e-10);
            Assert.AreEqual(10 * 0.3 * 0.7, variance[1], 1e-10);
            Assert.AreEqual(10 * 0.5 * 0.5, variance[2], 1e-10);
        }

        /// <summary>
        /// Test covariance. Cov(Xi, Xj) = -N * pi * pj.
        /// </summary>
        [TestMethod]
        public void Test_Covariance()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            Assert.AreEqual(-10 * 0.2 * 0.3, M.Covariance(0, 1), 1e-10);
            Assert.AreEqual(M.Variance[0], M.Covariance(0, 0), 1e-10);
        }

        /// <summary>
        /// Test PMF for a simple fair coin case: Mult(10, (0.5, 0.5)) = Binomial(10, 0.5).
        /// </summary>
        [TestMethod]
        public void Test_PMF_Binomial()
        {
            // Mult(10, (0.5, 0.5)) at x = (5, 5) should equal C(10,5) * 0.5^10
            // = 252 * 0.0009765625 = 0.24609375
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            double pmf = M.PDF(new double[] { 5, 5 });
            Assert.AreEqual(0.24609375, pmf, 1e-8);
        }

        /// <summary>
        /// Test PMF for a 3-category case.
        /// </summary>
        [TestMethod]
        public void Test_PMF_ThreeCategory()
        {
            // Mult(4, (0.2, 0.3, 0.5)) at x = (1, 1, 2)
            // = 4! / (1! 1! 2!) * 0.2^1 * 0.3^1 * 0.5^2
            // = 24 / (1*1*2) * 0.2 * 0.3 * 0.25
            // = 12 * 0.015 = 0.18
            var M = new Multinomial(4, new double[] { 0.2, 0.3, 0.5 });
            double pmf = M.PDF(new double[] { 1, 1, 2 });
            Assert.AreEqual(0.18, pmf, 1e-10);
        }

        /// <summary>
        /// Test that PMF returns 0 for invalid count vectors.
        /// </summary>
        [TestMethod]
        public void Test_PMF_Invalid()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            // Counts don't sum to N
            Assert.AreEqual(0.0, M.PDF(new double[] { 3, 3 }), 1e-10);
            // Negative count
            Assert.AreEqual(0.0, M.PDF(new double[] { -1, 11 }), 1e-10);
            // Wrong dimension
            Assert.AreEqual(0.0, M.PDF(new double[] { 5, 3, 2 }), 1e-10);
        }

        /// <summary>
        /// Test LogPMF consistency with PMF.
        /// </summary>
        [TestMethod]
        public void Test_LogPMF()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var x = new double[] { 2, 3, 5 };
            Assert.AreEqual(Math.Log(M.PDF(x)), M.LogPMF(x), 1e-10);
        }

        /// <summary>
        /// Test that invalid count vectors return the canonical log(0) sentinel.
        /// </summary>
        [TestMethod]
        public void Test_LogPMF_InvalidInput_ReturnsNegativeInfinity()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });

            Assert.AreEqual(double.NegativeInfinity, M.LogPMF(new double[] { 3, 3 }));
            Assert.AreEqual(double.NegativeInfinity, M.LogPMF(new double[] { -1, 11 }));
            Assert.AreEqual(double.NegativeInfinity, M.LogPMF(new double[] { 5, 3, 2 }));
        }

        /// <summary>
        /// Test random sampling: each sample sums to N and moments converge.
        /// </summary>
        [TestMethod]
        public void Test_Sampling()
        {
            var M = new Multinomial(100, new double[] { 0.2, 0.3, 0.5 });
            var samples = M.GenerateRandomValues(5000, seed: 42);

            // All samples should sum to N
            for (int i = 0; i < 5000; i++)
            {
                double sum = 0;
                for (int j = 0; j < 3; j++)
                {
                    Assert.IsGreaterThanOrEqualTo(0, samples[i, j], $"Sample [{i},{j}] is negative");
                    sum += samples[i, j];
                }
                Assert.AreEqual(100.0, sum, 1e-10, $"Sample {i} does not sum to N");
            }

            // Sample means should converge to theoretical means
            var sampleMeans = new double[3];
            for (int i = 0; i < 5000; i++)
                for (int j = 0; j < 3; j++)
                    sampleMeans[j] += samples[i, j];

            var theoreticalMeans = M.Mean;
            for (int j = 0; j < 3; j++)
            {
                sampleMeans[j] /= 5000;
                Assert.AreEqual(theoreticalMeans[j], sampleMeans[j], 1.0); // within 1 count of N*p
            }
        }

        /// <summary>
        /// Test the static Sample method for weighted categorical sampling.
        /// </summary>
        [TestMethod]
        public void Test_WeightedSample()
        {
            var rng = new MersenneTwister(12345);
            var weights = new double[] { 1.0, 3.0, 6.0 }; // 10%, 30%, 60%
            var counts = new int[3];

            int n = 10000;
            for (int i = 0; i < n; i++)
            {
                int idx = Multinomial.Sample(weights, rng);
                Assert.IsGreaterThanOrEqualTo(0, idx);
                Assert.IsLessThan(3, idx);
                counts[idx]++;
            }

            // Check proportions
            Assert.AreEqual(0.1, counts[0] / (double)n, 0.02);
            Assert.AreEqual(0.3, counts[1] / (double)n, 0.02);
            Assert.AreEqual(0.6, counts[2] / (double)n, 0.02);
        }

        /// <summary>
        /// Test the static Sample method with edge cases.
        /// </summary>
        [TestMethod]
        public void Test_WeightedSample_EdgeCases()
        {
            var rng = new MersenneTwister(12345);
            // All weight on one category
            Assert.AreEqual(1, Multinomial.Sample(new double[] { 0.0, 1.0, 0.0 }, rng));

            // Single positive weight
            Assert.AreEqual(0, Multinomial.Sample(new double[] { 5.0, 0.0 }, rng));

            // Invalid: all zeros
            Assert.Throws<ArgumentException>(() => Multinomial.Sample(new double[] { 0.0, 0.0 }, rng));
        }

        /// <summary>
        /// Test Clone produces an independent copy.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var M = new Multinomial(10, new double[] { 0.2, 0.3, 0.5 });
            var M2 = M.Clone() as Multinomial;
            Assert.IsNotNull(M2);
            Assert.AreEqual(M.Dimension, M2.Dimension);
            Assert.AreEqual(M.NumberOfTrials, M2.NumberOfTrials);
            Assert.IsTrue(M2.ParametersValid);
        }

        /// <summary>
        /// Test CDF throws NotImplementedException.
        /// </summary>
        [TestMethod]
        public void Test_CDF_Throws()
        {
            var M = new Multinomial(10, new double[] { 0.5, 0.5 });
            Assert.Throws<NotImplementedException>(() => M.CDF(new double[] { 5, 5 }));
        }

        /// <summary>
        /// The categorical index map reproduces the reference vectors and numpy
        /// clip(searchsorted(c, u, side="right"), 0, n - 1): boundaries select the next category, a zero-weight
        /// category is skipped, and variates at or above the last entry clip to the last category.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_ReferenceVectors()
        {
            var uniforms = new[] { 0.0, 0.1999999999999999, 0.2, 0.2000000000000001, 0.5, 0.75, 0.9999999999999999, 1.0, 1.5 };
            CollectionAssert.AreEqual(new[] { 0, 0, 1, 1, 2, 2, 2, 2, 2 }, Multinomial.CategoricalIndices(new[] { 0.2, 0.5, 1.0 }, uniforms));
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 0, 2, 2, 2, 2, 2 }, Multinomial.CategoricalIndices(new[] { 0.5, 0.5, 1.0 }, uniforms));

            // numpy oracle values, including infinities, negative variates and repeated entries.
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 1, 2, 2, 2 },
                Multinomial.CategoricalIndices(new[] { 0.2, 0.5, 1.0 }, new[] { double.NegativeInfinity, -1.0, 0.0, 0.2, 0.5, 1.0, double.PositiveInfinity }));
            CollectionAssert.AreEqual(new[] { 0, 2, 2, 5, 5, 5, 6, 6 },
                Multinomial.CategoricalIndices(new[] { 0.0, 0.0, 0.3, 0.3, 0.3, 0.9, 0.9 }, new[] { -0.1, 0.0, 0.1, 0.3, 0.29999999999999999, 0.5, 0.9, 0.95 }));
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, Multinomial.CategoricalIndices(new[] { 2.5 }, new[] { 0.0, 2.5, 3.0 }));
            Assert.IsEmpty(Multinomial.CategoricalIndices(new[] { 0.2, 1.0 }, new double[0]));
        }

        /// <summary>
        /// On random unnormalized tables with ties, every index equals the definition (the first k with
        /// c[k] &gt; u, else n - 1) found by a linear scan, and the index frequencies match the weights.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_MatchesLinearScanAndFrequencies()
        {
            var rng = new MersenneTwister(2026);
            for (int trial = 0; trial < 20; trial++)
            {
                int n = 1 + trial * 7;
                var cumulative = new double[n];
                double running = 0;
                for (int k = 0; k < n; k++)
                {
                    running += rng.NextDouble() < 0.2 ? 0.0 : 3.0 * rng.NextDouble();
                    cumulative[k] = running;
                }
                var uniforms = new double[500];
                for (int i = 0; i < uniforms.Length; i++)
                    uniforms[i] = i % 10 == 0 ? cumulative[rng.Next(0, n)] : (1.2 * rng.NextDouble() - 0.1) * running;
                var indices = Multinomial.CategoricalIndices(cumulative, uniforms);
                for (int i = 0; i < uniforms.Length; i++)
                {
                    int expected = n - 1;
                    for (int k = 0; k < n; k++)
                        if (cumulative[k] > uniforms[i]) { expected = k; break; }
                    Assert.AreEqual(expected, indices[i]);
                }
            }

            var weights = new[] { 0.1, 0.0, 0.25, 0.4, 0.25 };
            var table = new double[weights.Length];
            double sum = 0;
            for (int k = 0; k < weights.Length; k++) { sum += weights[k]; table[k] = sum; }
            var draws = new double[200000];
            for (int i = 0; i < draws.Length; i++) draws[i] = rng.NextDouble();
            var counts = new int[weights.Length];
            foreach (int k in Multinomial.CategoricalIndices(table, draws)) counts[k]++;
            for (int k = 0; k < weights.Length; k++)
                Assert.AreEqual(weights[k], (double)counts[k] / draws.Length, 0.005);
            Assert.AreEqual(0, counts[1]);
        }

        /// <summary>
        /// Null arrays, an empty table, non-finite or decreasing table entries, and NaN variates throw the
        /// documented exceptions naming the argument.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_InvalidInputs()
        {
            var table = new[] { 0.2, 0.5, 1.0 };
            var u = new[] { 0.3 };
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentNullException>(() => Multinomial.CategoricalIndices(null!, u)).ParamName);
            Assert.AreEqual("uniforms", Assert.Throws<ArgumentNullException>(() => Multinomial.CategoricalIndices(table, null!)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentException>(() => Multinomial.CategoricalIndices(new double[0], u)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentException>(() => Multinomial.CategoricalIndices(new[] { 0.5, 0.4, 1.0 }, u)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentOutOfRangeException>(() => Multinomial.CategoricalIndices(new[] { 0.5, double.NaN }, u)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentOutOfRangeException>(() => Multinomial.CategoricalIndices(new[] { 0.5, double.PositiveInfinity }, u)).ParamName);
            Assert.AreEqual("uniforms", Assert.Throws<ArgumentOutOfRangeException>(() => Multinomial.CategoricalIndices(table, new[] { 0.1, double.NaN })).ParamName);
        }

        /// <summary>
        /// The inclusive form reproduces numpy clip(searchsorted(c, u, side="left"), 0, n - 1): a variate on a boundary
        /// c[k] selects category k (the lowest k of a repeated entry), u = 0 selects category 0 even when c[0] is zero,
        /// variates near and at one clip as the definition requires, and the strict default is unchanged.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_InclusiveReferenceVectors()
        {
            var uniforms = new[] { 0.0, 0.1999999999999999, 0.2, 0.2000000000000001, 0.5, 0.75, 0.9999999999999999, 1.0, 1.5 };
            var table = new[] { 0.2, 0.5, 1.0 };
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 1, 1, 2, 2, 2, 2 }, Multinomial.CategoricalIndices(table, uniforms, true));
            CollectionAssert.AreEqual(new[] { 0, 0, 1, 1, 2, 2, 2, 2, 2 }, Multinomial.CategoricalIndices(table, uniforms, false));
            CollectionAssert.AreEqual(Multinomial.CategoricalIndices(table, uniforms), Multinomial.CategoricalIndices(table, uniforms, false));

            // numpy side="left" oracle values on a table with zero weights, repeated entries and saturation at one.
            CollectionAssert.AreEqual(new[] { 0, 0, 0, 2, 2, 2, 5, 5, 6 },
                Multinomial.CategoricalIndices(new[] { 0.0, 0.0, 0.3, 0.3, 0.3, 0.9, 0.9 },
                    new[] { double.NegativeInfinity, -0.1, 0.0, 0.1, 0.3, 0.29999999999999999, 0.5, 0.9, 0.95 }, true));
            CollectionAssert.AreEqual(new[] { 0, 2, 2, 3, 4 },
                Multinomial.CategoricalIndices(new[] { 0.25, 0.5, 0.9999999999999999, 1.0, 1.0 },
                    new[] { 0.0, 0.75, 0.9999999999999999, 1.0, double.PositiveInfinity }, true));
            CollectionAssert.AreEqual(new[] { 0, 0, 0 }, Multinomial.CategoricalIndices(new[] { 2.5 }, new[] { 0.0, 2.5, 3.0 }, true));
            Assert.IsEmpty(Multinomial.CategoricalIndices(new[] { 0.2, 1.0 }, new double[0], true));
        }

        /// <summary>
        /// On a Poisson CDF table the inclusive form is the inverse-transform count min{k : F(k) ≥ u}: at every exact
        /// tie u = F(k) it returns k where the strict form returns k + 1, and on random variates both forms equal
        /// their linear-scan definitions.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_InclusivePoissonTiesAndLinearScan()
        {
            var poisson = new Poisson(3.5);
            var cdf = new double[40];
            for (int k = 0; k < cdf.Length; k++) cdf[k] = poisson.CDF(k);
            int last = cdf.Length - 1;

            var ties = new double[last];
            Array.Copy(cdf, ties, last);
            var inclusive = Multinomial.CategoricalIndices(cdf, ties, true);
            var strict = Multinomial.CategoricalIndices(cdf, ties, false);
            for (int k = 0; k < last; k++)
            {
                int first = Array.IndexOf(cdf, cdf[k]);
                Assert.AreEqual(first, inclusive[k], $"tie at k = {k}");
                int expectedStrict = last;
                for (int j = 0; j < last; j++) if (cdf[j] > cdf[k]) { expectedStrict = j; break; }
                Assert.AreEqual(expectedStrict, strict[k], $"tie at k = {k}");
            }
            Assert.AreEqual(0, Multinomial.CategoricalIndices(cdf, new[] { 0.0 }, true)[0]);

            var rng = new MersenneTwister(7);
            var uniforms = new double[2000];
            for (int i = 0; i < uniforms.Length; i++) uniforms[i] = i % 5 == 0 ? cdf[rng.Next(0, cdf.Length)] : rng.NextDouble();
            uniforms[0] = 1d - Math.Pow(2, -53);
            var mapped = Multinomial.CategoricalIndices(cdf, uniforms, true);
            for (int i = 0; i < uniforms.Length; i++)
            {
                int expected = last;
                for (int k = 0; k < last; k++)
                    if (cdf[k] >= uniforms[i]) { expected = k; break; }
                Assert.AreEqual(expected, mapped[i]);
            }
        }

        /// <summary>
        /// The inclusive overload validates its arguments exactly as the strict form does.
        /// </summary>
        [TestMethod]
        public void Test_CategoricalIndices_InclusiveInvalidInputs()
        {
            var table = new[] { 0.2, 0.5, 1.0 };
            var u = new[] { 0.3 };
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentNullException>(() => Multinomial.CategoricalIndices(null!, u, true)).ParamName);
            Assert.AreEqual("uniforms", Assert.Throws<ArgumentNullException>(() => Multinomial.CategoricalIndices(table, null!, true)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentException>(() => Multinomial.CategoricalIndices(new double[0], u, true)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentException>(() => Multinomial.CategoricalIndices(new[] { 0.5, 0.4, 1.0 }, u, true)).ParamName);
            Assert.AreEqual("cumulativeWeights", Assert.Throws<ArgumentOutOfRangeException>(() => Multinomial.CategoricalIndices(new[] { 0.5, double.NegativeInfinity }, u, true)).ParamName);
            Assert.AreEqual("uniforms", Assert.Throws<ArgumentOutOfRangeException>(() => Multinomial.CategoricalIndices(table, new[] { double.NaN }, true)).ParamName);
        }

    }
}
