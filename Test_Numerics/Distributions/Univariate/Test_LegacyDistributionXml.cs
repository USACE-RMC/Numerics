using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System;
using System.Linq;

namespace Distributions.Univariate
{
    /// <summary>
    /// Golden-fixture tests reading distribution XML exactly as version 2.1.4 wrote it.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// The fixture strings below were captured verbatim from the 2.1.4 writers. They pin the
    /// factory's contract for files produced by earlier releases: well-formed payloads load with
    /// their exact values, while payloads carrying non-finite parameters, missing attributes,
    /// unknown types, or damaged correlation matrices are rejected loudly instead of loading
    /// degraded values. The 2.1.4 scalar writer could not serialize the table-bearing
    /// EmpiricalDistribution or KernelDensity at all (it faulted on their parameter arrays), so no
    /// legacy scalar-writer files exist for those types; the tableless-element tests cover
    /// externally authored payloads.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_LegacyDistributionXml
    {
        private const string NormalAllAttributes =
            "<Distribution Type=\"Normal\" Mu=\"3.5\" Sigma=\"0.75\" />";

        private const string GevAllAttributes =
            "<Distribution Type=\"GeneralizedExtremeValue\" Xi=\"100\" Alpha=\"25\" Kappa=\"-0.10000000000000001\" />";

        private const string NormalNaNParameters =
            "<Distribution Type=\"Normal\" Mu=\"NaN\" Sigma=\"NaN\" />";

        private const string CompetingRisksCorrelationMatrix =
            "<Distribution Type=\"CompetingRisks\" XTransform=\"None\" ProbabilityTransform=\"NormalZ\" " +
            "MinimumOfRandomVariables=\"True\" Dependency=\"CorrelationMatrix\" Distributions=\"Normal|Normal\" " +
            "Parameters=\"10|2|14|3\"><CorrelationMatrix><Correlation_Row>1|0.5</Correlation_Row>" +
            "<Correlation_Row>0.5|1</Correlation_Row></CorrelationMatrix></Distribution>";

        /// <summary>
        /// A well-formed 2.1.4 scalar payload loads with its exact parameter values.
        /// </summary>
        [TestMethod]
        public void LegacyScalarXml_AllAttributes_Loads()
        {
            var normal = (Normal)UnivariateDistributionFactory.CreateDistribution(XElement.Parse(NormalAllAttributes));
            Assert.AreEqual(3.5, normal.Mu, 0d);
            Assert.AreEqual(0.75, normal.Sigma, 0d);

            var gev = (GeneralizedExtremeValue)UnivariateDistributionFactory.CreateDistribution(XElement.Parse(GevAllAttributes));
            Assert.AreEqual(100d, gev.Xi, 0d);
            Assert.AreEqual(25d, gev.Alpha, 0d);
            Assert.AreEqual(-0.1, gev.Kappa, 0d);
        }

        /// <summary>
        /// A 2.1.4 payload carrying NaN parameters (a failed-fit artifact the old writer emitted
        /// silently) is rejected instead of loading a degraded distribution.
        /// </summary>
        [TestMethod]
        public void LegacyScalarXml_NaNParameters_Rejected()
        {
            Assert.Throws<ArgumentException>(
                () => UnivariateDistributionFactory.CreateDistribution(XElement.Parse(NormalNaNParameters)));
        }

        /// <summary>
        /// A legacy payload missing a parameter attribute is rejected instead of defaulting the
        /// parameter to zero.
        /// </summary>
        [TestMethod]
        public void LegacyScalarXml_MissingParameterAttribute_Rejected()
        {
            var element = XElement.Parse(NormalAllAttributes);
            element.Attribute("Sigma")!.Remove();
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(element));
        }

        /// <summary>
        /// A legacy payload whose type is not a defined distribution is rejected instead of
        /// defaulting to a deterministic distribution.
        /// </summary>
        [TestMethod]
        public void LegacyScalarXml_UnknownType_Rejected()
        {
            var element = XElement.Parse(NormalAllAttributes);
            element.SetAttributeValue("Type", "NotADistribution");
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(element));
        }

        /// <summary>
        /// A well-formed 2.1.4 competing-risks payload loads with its children, dependency, and
        /// correlation matrix intact.
        /// </summary>
        [TestMethod]
        public void LegacyCompetingRisksXml_WellFormed_Loads()
        {
            var element = XElement.Parse(CompetingRisksCorrelationMatrix);
            var risks = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(element);
            Assert.HasCount(2, risks.Distributions);
            Assert.AreEqual(Probability.DependencyType.CorrelationMatrix, risks.Dependency);
            Assert.AreEqual(0.5, risks.CorrelationMatrix[0, 1], 0d);
            Assert.AreEqual(10d, ((Normal)risks.Distributions[0]).Mu, 0d);
            Assert.AreEqual(3d, ((Normal)risks.Distributions[1]).Sigma, 0d);
        }

        /// <summary>
        /// A competing-risks payload with a truncated or unparseable correlation matrix is rejected
        /// instead of loading NaN matrix entries.
        /// </summary>
        [TestMethod]
        public void LegacyCompetingRisksXml_DamagedCorrelationMatrix_Rejected()
        {
            var truncated = XElement.Parse(CompetingRisksCorrelationMatrix);
            truncated.Element("CorrelationMatrix")!.Elements("Correlation_Row").Last().Remove();
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(truncated));

            var corrupted = XElement.Parse(CompetingRisksCorrelationMatrix);
            corrupted.Element("CorrelationMatrix")!.Elements("Correlation_Row").First().Value = "1|abc";
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(corrupted));
        }

        /// <summary>
        /// Tableless empirical and kernel-density elements are rejected: the table attributes are
        /// required, and no legacy writer ever produced a valid tableless payload.
        /// </summary>
        [TestMethod]
        public void LegacyTableBearingXml_TablelessElements_Rejected()
        {
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(
                XElement.Parse("<Distribution Type=\"Empirical\" />")));
            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(
                XElement.Parse("<Distribution Type=\"KernelDensity\" />")));
        }
    }
}
