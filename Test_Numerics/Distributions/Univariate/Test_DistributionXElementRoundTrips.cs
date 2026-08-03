using System;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Unit tests for empirical and kernel-density XML round trips, including their data tables,
    /// transforms, kernels, bandwidths, and optional sample weights.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_DistributionXElementRoundTrips
    {
        /// <summary>
        /// Test the empirical distribution round-trip: the X/probability tables, the
        /// transforms, and evaluation restore exactly, through both the direct FromXElement and
        /// the distribution factory.
        /// </summary>
        [TestMethod]
        public void Test_EmpiricalDistribution_RoundTrip()
        {
            var original = new EmpiricalDistribution(
                new[] { 100d, 250d, 600d, 2000d },
                new[] { 0.05d, 0.4d, 0.8d, 0.99d })
            {
                XTransform = Transform.Logarithmic,
                ProbabilityTransform = Transform.NormalZ,
            };

            var element = original.ToXElement();
            var restored = (EmpiricalDistribution)UnivariateDistributionFactory.CreateDistribution(element);

            Assert.HasCount(original.XValues.Count, restored.XValues);
            for (int i = 0; i < original.XValues.Count; i++)
            {
                Assert.AreEqual(original.XValues[i], restored.XValues[i], 0d);
                Assert.AreEqual(original.ProbabilityValues[i], restored.ProbabilityValues[i], 0d);
            }
            Assert.AreEqual(original.XTransform, restored.XTransform);
            Assert.AreEqual(original.ProbabilityTransform, restored.ProbabilityTransform);
            Assert.AreEqual(original.CDF(500d), restored.CDF(500d), 1E-12);
            Assert.AreEqual(original.InverseCDF(0.5d), restored.InverseCDF(0.5d), 1E-12);

            // A serialized empirical distribution requires both data tables.
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.FromXElement(new System.Xml.Linq.XElement("Distribution")));
        }

        /// <summary>
        /// Test that a descending (exceedance-convention) probability ladder round-trips with
        /// its order preserved.
        /// </summary>
        [TestMethod]
        public void Test_EmpiricalDistribution_DescendingProbabilities_RoundTrip()
        {
            var original = new EmpiricalDistribution(
                new[] { 1d, 10d, 100d },
                new[] { 0.999d, 0.5d, 0.001d },
                SortOrder.Ascending,
                SortOrder.Descending);

            var restored = (EmpiricalDistribution)UnivariateDistributionFactory.CreateDistribution(original.ToXElement());
            for (int i = 0; i < original.ProbabilityValues.Count; i++)
            {
                Assert.AreEqual(original.ProbabilityValues[i], restored.ProbabilityValues[i], 0d, "The exceedance ladder must restore in its stored order.");
            }
            Assert.AreEqual(original.CDF(10d), restored.CDF(10d), 1E-12);
        }

        /// <summary>
        /// Test the kernel density round-trip: the sample, kernel type, bandwidth, transforms,
        /// and the optional per-sample weights restore exactly, through both the direct
        /// FromXElement and the distribution factory.
        /// </summary>
        [TestMethod]
        public void Test_KernelDensity_RoundTrip()
        {
            var sample = new[] { 1.2d, 3.4d, 2.2d, 5.6d, 4.4d, 3.1d, 2.7d, 4.9d };
            var original = new KernelDensity(sample, KernelDensity.KernelType.Epanechnikov, 0.75d);

            var element = original.ToXElement();
            var restored = (KernelDensity)UnivariateDistributionFactory.CreateDistribution(element);

            Assert.HasCount(original.SampleData.Count, restored.SampleData);
            for (int i = 0; i < original.SampleData.Count; i++)
                Assert.AreEqual(original.SampleData[i], restored.SampleData[i], 0d);
            Assert.AreEqual(original.KernelDistribution, restored.KernelDistribution);
            Assert.AreEqual(original.Bandwidth, restored.Bandwidth, 0d);
            Assert.AreEqual(original.PDF(3d), restored.PDF(3d), 1E-12);
            Assert.AreEqual(original.CDF(3d), restored.CDF(3d), 1E-12);

            // The weighted variant restores its weights (weighted evaluation must match).
            var weights = new[] { 1d, 2d, 1d, 3d, 1d, 2d, 1d, 1d };
            var weighted = new KernelDensity(sample, weights, KernelDensity.KernelType.Gaussian, 0.6d);
            var weightedRestored = (KernelDensity)UnivariateDistributionFactory.CreateDistribution(weighted.ToXElement());
            Assert.AreEqual(weighted.PDF(3d), weightedRestored.PDF(3d), 1E-12, "Weighted kernel evaluation must survive the round-trip.");
            Assert.AreEqual(weighted.Mean, weightedRestored.Mean, 1E-12);

            Assert.Throws<ArgumentException>(() => KernelDensity.FromXElement(new XElement("Distribution")));
        }

        /// <summary>
        /// Test that malformed serialized payloads are rejected: undefined enum ordinals on the
        /// empirical table, zero-sum kernel weights, non-finite bandwidths, and a missing
        /// distribution type discriminator.
        /// </summary>
        [TestMethod]
        public void Test_MalformedXElements_AreRejected()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d }, new[] { 0d, 1d });
            XElement invalidOrder = empirical.ToXElement();
            invalidOrder.SetAttributeValue("ProbabilityOrder", "999");
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(invalidOrder));

            XElement invalidTransform = empirical.ToXElement();
            invalidTransform.SetAttributeValue(nameof(EmpiricalDistribution.XTransform), "999");
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.FromXElement(invalidTransform));

            var kernel = new KernelDensity(new[] { 1d, 2d, 3d }, new[] { 1d, 1d, 1d }, KernelDensity.KernelType.Gaussian, 0.5d);
            XElement zeroWeights = kernel.ToXElement();
            zeroWeights.SetAttributeValue("Weights", "0|0|0");
            Assert.Throws<ArgumentException>(() => KernelDensity.FromXElement(zeroWeights));

            XElement invalidBandwidth = kernel.ToXElement();
            invalidBandwidth.SetAttributeValue(nameof(KernelDensity.Bandwidth), "NaN");
            Assert.Throws<ArgumentException>(() => KernelDensity.FromXElement(invalidBandwidth));

            XElement missingType = new Normal().ToXElement();
            missingType.Attribute(nameof(UnivariateDistributionBase.Type))!.Remove();
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(missingType));
        }
    }
}
