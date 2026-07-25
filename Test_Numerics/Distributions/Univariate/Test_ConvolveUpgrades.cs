using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Unit tests for the N8 convolution upgrades: the atom-aware exact-lattice
    /// <see cref="EmpiricalDistribution.ConvolveDiscrete"/> (point masses — e.g. the
    /// zero-inflation atom of a defective risk curve — have no representation in the
    /// continuous-PDF pipeline) and the opt-in log-spaced output grid.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_ConvolveUpgrades
    {
        /// <summary>
        /// Test the discrete convolution against exact enumeration: total mass is the product
        /// of the input totals, the mean is the sum of the input means (the moment-preserving
        /// two-node split makes this exact), and the cumulative mass at probes between the
        /// enumerated sum atoms matches the enumerated CDF.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveDiscrete_ExactEnumeration()
        {
            double[] values1 = { 0d, 10d };
            double[] masses1 = { 0.6d, 0.4d };
            double[] values2 = { 0d, 5d, 20d };
            double[] masses2 = { 0.5d, 0.3d, 0.2d };

            EmpiricalDistribution.ConvolveDiscrete(values1, masses1, values2, masses2, 4096, out var values, out var masses);

            double total = 0d, mean = 0d;
            for (int i = 0; i < values.Length; i++)
            {
                total += masses[i];
                mean += masses[i] * values[i];
            }
            Assert.AreEqual(1d, total, 1E-12, "Total mass is the product of the input totals.");

            // Exact input means: 4.0 and 5.5 → the convolved mean is 9.5 exactly (two-node
            // splits preserve each input mean).
            Assert.AreEqual(9.5d, mean, 1E-9, "The convolved mean is the sum of the input means.");

            // Enumerated sums: {0: 0.30, 5: 0.18, 10: 0.20, 15: 0.12, 20: 0.12, 30: 0.08}.
            // Cumulative probes at least one lattice node away from every atom are exact up to
            // the split smear.
            Assert.AreEqual(0.48d, CumulativeAt(values, masses, 7.5d), 1E-9, "P(sum ≤ 7.5).");
            Assert.AreEqual(0.68d, CumulativeAt(values, masses, 12.5d), 1E-9, "P(sum ≤ 12.5).");
            Assert.AreEqual(0.92d, CumulativeAt(values, masses, 25d), 1E-9, "P(sum ≤ 25).");
        }

        /// <summary>
        /// Test the motivating case: zero-inflation atoms of defective risk curves convolve
        /// exactly — the joint zero atom is the product of the input zero masses, a point mass
        /// the continuous-PDF pipeline cannot represent.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveDiscrete_ZeroAtoms()
        {
            double[] values1 = { 0d, 100d };
            double[] masses1 = { 0.7d, 0.3d };
            double[] values2 = { 0d, 50d };
            double[] masses2 = { 0.8d, 0.2d };

            EmpiricalDistribution.ConvolveDiscrete(values1, masses1, values2, masses2, 4096, out var values, out var masses);

            // P(both zero) = 0.7·0.8 = 0.56 sits at the origin node (both inputs' atoms land
            // exactly on lattice nodes here only at the ends; probe just above the origin).
            Assert.AreEqual(0.56d, CumulativeAt(values, masses, 1d), 1E-9, "The joint zero atom is exact.");
            Assert.AreEqual(0.56d + 0.7d * 0.2d + 0.3d * 0.8d, CumulativeAt(values, masses, 120d), 1E-9, "P(sum ≤ 120).");
            Assert.AreEqual(1d, CumulativeAt(values, masses, 151d), 1E-12, "The ladder is exhaustive.");

            // Guards: mismatched atom lists are loud.
            Assert.Throws<ArgumentException>(() =>
                EmpiricalDistribution.ConvolveDiscrete(new[] { 0d }, new[] { 0.5d, 0.5d }, values2, masses2, 4096, out _, out _));
            Assert.Throws<ArgumentException>(() =>
                EmpiricalDistribution.ConvolveDiscrete(new[] { 0d }, new[] { -0.5d }, values2, masses2, 4096, out _, out _));
        }

        /// <summary>
        /// Test the log-spaced output option: false delegates to the linear-grid method
        /// unchanged (the same instance), true re-reads the same convolved CDF on a log-spaced
        /// ladder (agreeing at probes) and requires a strictly positive summed support.
        /// </summary>
        [TestMethod]
        public void Test_Convolve_LogSpacedOutput()
        {
            var dist1 = new EmpiricalDistribution(new[] { 10d, 40d, 100d }, new[] { 0.1d, 0.6d, 0.95d });
            var dist2 = new EmpiricalDistribution(new[] { 5d, 20d, 60d }, new[] { 0.15d, 0.5d, 0.9d });

            var linear = EmpiricalDistribution.Convolve(dist1, dist2, 1024, logSpacedOutput: false);
            var logSpaced = EmpiricalDistribution.Convolve(dist1, dist2, 1024, logSpacedOutput: true);

            // The log-spaced ladder reads the same convolved CDF.
            foreach (double probe in new[] { 40d, 80d, 140d })
            {
                Assert.AreEqual(linear.CDF(probe), logSpaced.CDF(probe), 5E-3, $"CDF agreement at {probe}.");
            }
            Assert.IsGreaterThan(0d, logSpaced.Minimum, "The log-spaced ladder lives on positive support.");

            // A non-positive summed support cannot be log-spaced.
            var negative = new EmpiricalDistribution(new[] { -5d, 0d, 10d }, new[] { 0.1d, 0.5d, 0.9d });
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.Convolve(negative, dist2, 1024, logSpacedOutput: true));
        }

        /// <summary>
        /// Accumulates the lattice mass at and below a probe value.
        /// </summary>
        /// <param name="values">The lattice values.</param>
        /// <param name="masses">The lattice masses.</param>
        /// <param name="probe">The probe value.</param>
        /// <returns>The cumulative mass.</returns>
        private static double CumulativeAt(double[] values, double[] masses, double probe)
        {
            double sum = 0d;
            for (int i = 0; i < values.Length && values[i] <= probe; i++) sum += masses[i];
            return sum;
        }
    }
}
