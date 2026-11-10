using System;
using System.Linq;
using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Unit tests for distribution XML round trips: empirical and kernel-density data tables, transforms,
    /// kernels, bandwidths, and optional sample weights; non-parameter settings of the log families; and
    /// competing-risks components that need their own serialized elements.
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
        /// Test that a probability ladder explicitly configured without a sort order serializes
        /// and restores that order rather than inferring one from its endpoints.
        /// </summary>
        [TestMethod]
        public void Test_EmpiricalDistribution_NoneProbabilityOrder_RoundTrip()
        {
            var original = new EmpiricalDistribution(
                new[] { 1d, 2d, 3d },
                new[] { 0.1d, 0.9d, 0.4d },
                SortOrder.Ascending,
                SortOrder.None);

            XElement element = original.ToXElement();
            Assert.AreEqual(nameof(SortOrder.None),
                element.Attribute(nameof(EmpiricalDistribution.ProbabilityOrder))?.Value);

            var restored = EmpiricalDistribution.FromXElement(element);
            Assert.AreEqual(SortOrder.None, restored.ProbabilityOrder);
            Assert.HasCount(original.ProbabilityValues.Count, restored.ProbabilityValues);
            for (int i = 0; i < original.ProbabilityValues.Count; i++)
                Assert.AreEqual(original.ProbabilityValues[i], restored.ProbabilityValues[i], 0d);
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
            invalidOrder.SetAttributeValue(nameof(EmpiricalDistribution.ProbabilityOrder), "999");
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

        /// <summary>
        /// A logarithm base other than 10 round-trips on the log-Normal and log-Pearson Type III families, and the
        /// restored distributions evaluate the natural-log model. Reference values from
        /// scipy.stats.norm.cdf((ln 3 - 1) / 0.2) and scipy.stats.pearson3.cdf(ln 3, 0.1, loc=1, scale=0.2).
        /// </summary>
        [TestMethod]
        public void Test_LogFamilies_NonDefaultBase_RoundTrip()
        {
            var logNormal = new LogNormal(1d, 0.2d) { Base = Math.E };
            var logPearson = new LogPearsonTypeIII(1d, 0.2d, 0.1d) { Base = Math.E };

            var restoredNormal = (LogNormal)UnivariateDistributionFactory.CreateDistribution(logNormal.ToXElement());
            var restoredPearson = (LogPearsonTypeIII)UnivariateDistributionFactory.CreateDistribution(logPearson.ToXElement());

            Assert.AreEqual(Math.E, restoredNormal.Base, 0d);
            Assert.AreEqual(Math.E, restoredPearson.Base, 0d);
            Assert.AreEqual(logNormal.CDF(3d), restoredNormal.CDF(3d), 0d);
            Assert.AreEqual(logPearson.CDF(3d), restoredPearson.CDF(3d), 0d);
            Assert.AreEqual(0.6890154133836028, restoredNormal.CDF(3d), 1E-12);
            Assert.AreEqual(0.6934665774884826, restoredPearson.CDF(3d), 1E-9);
        }

        /// <summary>
        /// Default-configured log families serialize without the optional attributes, in exactly the form
        /// written before those attributes existed, and an invalid serialized base is rejected.
        /// </summary>
        [TestMethod]
        public void Test_LogFamilies_DefaultBase_SerializedFormUnchanged()
        {
            Assert.AreEqual("<Distribution Type=\"LogPearsonTypeIII\" Mu=\"1\" Sigma=\"0.20000000000000001\" Gamma=\"0.10000000000000001\" />",
                new LogPearsonTypeIII(1d, 0.2d, 0.1d).ToXElement().ToString(SaveOptions.DisableFormatting));
            Assert.AreEqual("<Distribution Type=\"LogNormal\" Mu=\"1\" Sigma=\"0.20000000000000001\" />",
                new LogNormal(1d, 0.2d).ToXElement().ToString(SaveOptions.DisableFormatting));
            Assert.AreEqual("<Distribution Type=\"LnNormal\" Mean=\"10\" StandardDeviation=\"2\" />",
                new LnNormal(10d, 2d).ToXElement().ToString(SaveOptions.DisableFormatting));

            XElement legacy = XElement.Parse("<Distribution Type=\"LogNormal\" Mu=\"1\" Sigma=\"0.2\" />");
            Assert.AreEqual(10d, ((LogNormal)UnivariateDistributionFactory.CreateDistribution(legacy)).Base, 0d);

            XElement invalidBase = new LogNormal(1d, 0.2d) { Base = Math.E }.ToXElement();
            invalidBase.SetAttributeValue(nameof(LogNormal.Base), "1");
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(invalidBase));
            invalidBase.SetAttributeValue(nameof(LogNormal.Base), "NaN");
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(invalidBase));
        }

        /// <summary>
        /// An Ln-Normal whose log coordinates were set directly round-trips those coordinates bit for bit and stays
        /// off the stored physical-moment surface; its mean matches the closed form exp(mu + sigma^2/2).
        /// </summary>
        [TestMethod]
        public void Test_LnNormal_LogCoordinates_RoundTrip()
        {
            var original = new LnNormal { Mu = 2.2829747364174051d, Sigma = 0.198042200435365d };

            XElement element = original.ToXElement();
            var restored = (LnNormal)UnivariateDistributionFactory.CreateDistribution(element);

            Assert.IsNotNull(element.Attribute(nameof(LnNormal.Mu)));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(original.Mu), BitConverter.DoubleToInt64Bits(restored.Mu));
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(original.Sigma), BitConverter.DoubleToInt64Bits(restored.Sigma));
            Assert.AreEqual(original.Mean, restored.Mean, 0d);
            Assert.AreEqual(original.StandardDeviation, restored.StandardDeviation, 0d);
            Assert.AreEqual(Math.Exp(2.2829747364174051d + 0.198042200435365d * 0.198042200435365d / 2d), restored.Mean, 1E-13);

            XElement onlyMu = original.ToXElement();
            onlyMu.Attribute(nameof(LnNormal.Sigma))!.Remove();
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(onlyMu));
        }

        /// <summary>
        /// A competing-risks distribution with an empirical component serializes, restores the empirical table and
        /// the parametric component, and evaluates identically; its flattened parameters round-trip through
        /// SetParameters.
        /// </summary>
        [TestMethod]
        public void Test_CompetingRisks_EmpiricalComponent_RoundTrip()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.1d, 0.5d, 0.9d }) { XTransform = Transform.Logarithmic };
            var original = new CompetingRisks(new UnivariateDistributionBase[] { new Normal(2d, 1d), empirical });

            XElement element = original.ToXElement();
            var restored = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(element);

            var restoredEmpirical = (EmpiricalDistribution)restored.Distributions[1];
            CollectionAssert.AreEqual(empirical.XValues.ToArray(), restoredEmpirical.XValues.ToArray());
            CollectionAssert.AreEqual(empirical.ProbabilityValues.ToArray(), restoredEmpirical.ProbabilityValues.ToArray());
            Assert.AreEqual(Transform.Logarithmic, restoredEmpirical.XTransform);
            CollectionAssert.AreEqual(original.GetParameters, restored.GetParameters);
            foreach (double x in new[] { 1.2d, 1.8d, 2.5d })
                Assert.AreEqual(original.CDF(x), restored.CDF(x), 0d);

            original.SetParameters(new[] { 2.5d, 1.5d });
            Assert.AreEqual(2.5d, ((Normal)original.Distributions[0]).Mu, 0d);
            Assert.AreEqual(1.5d, ((Normal)original.Distributions[0]).Sigma, 0d);
            Assert.ThrowsExactly<ArgumentException>(() => original.SetParameters(new[] { 2.5d, 1.5d, 0d, 1d }));
        }

        /// <summary>
        /// Competing-risks components keep a non-default logarithm base through XML, while a configuration that
        /// needs no nested components serializes exactly as before and a component list that disagrees with the
        /// listed types is rejected.
        /// </summary>
        [TestMethod]
        public void Test_CompetingRisks_NestedComponents_RoundTrip()
        {
            var original = new CompetingRisks(new UnivariateDistributionBase[]
            {
                new LogPearsonTypeIII(1d, 0.2d, 0.1d) { Base = Math.E },
                new LogPearsonTypeIII(1.1d, 0.2d, 0.1d) { Base = Math.E }
            });
            var restored = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(original.ToXElement());
            Assert.AreEqual(Math.E, ((LogPearsonTypeIII)restored.Distributions[0]).Base, 0d);
            Assert.AreEqual(Math.E, ((LogPearsonTypeIII)restored.Distributions[1]).Base, 0d);
            Assert.AreEqual(original.CDF(3d), restored.CDF(3d), 0d);

            var plain = new CompetingRisks(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            Assert.AreEqual("<Distribution Type=\"CompetingRisks\" XTransform=\"None\" ProbabilityTransform=\"NormalZ\" MinimumOfRandomVariables=\"True\" "
                + "Dependency=\"Independent\" PRNGSeed=\"12345\" Distributions=\"Weibull|Weibull\" Parameters=\"5|1.3999999999999999|12|2.2999999999999998\">"
                + "<CorrelationMatrix /></Distribution>", plain.ToXElement().ToString(SaveOptions.DisableFormatting));

            XElement mismatched = original.ToXElement();
            mismatched.Element("Components")!.Elements().First().Remove();
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(mismatched));
            XElement wrongType = original.ToXElement();
            wrongType.Element("Components")!.Elements().First().ReplaceWith(new Normal().ToXElement());
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(wrongType));
        }

        /// <summary>
        /// Mixture and nested competing-risks components, which cannot be created from their type alone, are restored
        /// from their nested elements with identical evaluated values.
        /// </summary>
        [TestMethod]
        public void Test_CompetingRisks_CompositeComponents_RoundTrip()
        {
            var mixture = new Mixture(new[] { 0.3d, 0.7d }, new UnivariateDistributionBase[] { new Normal(10d, 2d), new Normal(15d, 3d) });
            var inner = new CompetingRisks(new UnivariateDistributionBase[] { new Weibull(5d, 1.4d), new Weibull(12d, 2.3d) });
            var original = new CompetingRisks(new UnivariateDistributionBase[] { mixture, inner, new Normal(12d, 4d) });

            var restored = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(original.ToXElement());

            Assert.IsInstanceOfType<Mixture>(restored.Distributions[0]);
            Assert.IsInstanceOfType<CompetingRisks>(restored.Distributions[1]);
            foreach (double x in new[] { 2d, 6d, 11d, 14d, 20d })
            {
                Assert.AreEqual(original.CDF(x), restored.CDF(x), 0d);
                Assert.AreEqual(original.PDF(x), restored.PDF(x), 0d);
            }
            Assert.AreEqual(original.ToXElement().ToString(SaveOptions.DisableFormatting), restored.ToXElement().ToString(SaveOptions.DisableFormatting));
        }

        /// <summary>
        /// A truncated distribution keeps its bounds through the XElement round trip, including a base distribution with
        /// parameters named Min and Max, and inconsistent bounds are rejected.
        /// </summary>
        [TestMethod]
        public void Test_TruncatedDistribution_RoundTrip()
        {
            foreach (var original in new[]
            {
                new TruncatedDistribution(new GammaDistribution(2d, 3d), 0.5d, 10d),
                new TruncatedDistribution(new Uniform(0d, 20d), 2d, 8d),
                new TruncatedDistribution(new Normal(10d, 2d), 6d, 12d),
            })
            {
                var restored = (TruncatedDistribution)UnivariateDistributionFactory.CreateDistribution(original.ToXElement());
                Assert.AreEqual(original.BaseDistribution.GetType(), restored.BaseDistribution.GetType());
                Assert.AreEqual(original.Min, restored.Min);
                Assert.AreEqual(original.Max, restored.Max);
                foreach (double p in new[] { 0.05d, 0.5d, 0.95d })
                {
                    double x = original.InverseCDF(p);
                    Assert.AreEqual(original.CDF(x), restored.CDF(x), 0d);
                    Assert.AreEqual(original.PDF(x), restored.PDF(x), 0d);
                }
                Assert.AreEqual(original.ToXElement().ToString(SaveOptions.DisableFormatting), restored.ToXElement().ToString(SaveOptions.DisableFormatting));
            }

            XElement oneBound = new TruncatedDistribution(new GammaDistribution(2d, 3d), 0.5d, 10d).ToXElement();
            oneBound.Attribute("TruncationMax")!.Remove();
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(oneBound));
            XElement reversed = new TruncatedDistribution(new GammaDistribution(2d, 3d), 0.5d, 10d).ToXElement();
            reversed.SetAttributeValue("TruncationMin", "20");
            Assert.ThrowsExactly<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(reversed));
        }
    }
}
