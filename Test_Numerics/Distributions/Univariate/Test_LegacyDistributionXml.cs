using System.Xml.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
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

        private const string CompetingRisksCompatibilityTemplate =
            "<Distribution Type=\"CompetingRisks\" XTransform=\"Logarithmic\" ProbabilityTransform=\"NormalZ\" " +
            "MinimumOfRandomVariables=\"False\" Dependency=\"{0}\" PRNGSeed=\"24680\" " +
            "Distributions=\"Normal|Gumbel\" Parameters=\"10|2|14|3\">{1}</Distribution>";

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

            var roundTrip = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(risks.ToXElement());
            Assert.AreEqual(Probability.DependencyType.CorrelationMatrix, roundTrip.Dependency);
            Assert.IsTrue(roundTrip.MinimumOfRandomVariables);
            Assert.AreEqual(1d, roundTrip.CorrelationMatrix[0, 0], 0d);
            Assert.AreEqual(0.5d, roundTrip.CorrelationMatrix[0, 1], 0d);
            Assert.AreEqual(0.5d, roundTrip.CorrelationMatrix[1, 0], 0d);
            Assert.AreEqual(1d, roundTrip.CorrelationMatrix[1, 1], 0d);
            CollectionAssert.AreEqual(new[] { 10d, 2d, 14d, 3d }, roundTrip.GetParameters);
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
        /// A version-2 competing-risks payload may omit its optional correlation matrix or carry
        /// the empty/all-zero placeholders written before dependency-aware matrix persistence.
        /// Import preserves every other field and leaves the matrix unconfigured.
        /// </summary>
        /// <param name="dependency">The dependency value stored by the legacy file.</param>
        /// <param name="matrixXml">The legacy optional matrix representation.</param>
        [TestMethod]
        [DataRow("Independent", "")]
        [DataRow("Independent", "<CorrelationMatrix />")]
        [DataRow("Independent", "<CorrelationMatrix>  \r\n  </CorrelationMatrix>")]
        [DataRow("Independent", "<CorrelationMatrix><Correlation_Row>0|0</Correlation_Row><Correlation_Row>0|0</Correlation_Row></CorrelationMatrix>")]
        [DataRow("CorrelationMatrix", "")]
        [DataRow("CorrelationMatrix", "<CorrelationMatrix />")]
        [DataRow("CorrelationMatrix", "<CorrelationMatrix>  \r\n  </CorrelationMatrix>")]
        [DataRow("CorrelationMatrix", "<CorrelationMatrix><Correlation_Row>0|0</Correlation_Row><Correlation_Row>0|0</Correlation_Row></CorrelationMatrix>")]
        public void Version2CompetingRisksXml_MissingOrPlaceholderMatrix_LoadsWithoutChangingConfiguration(
            string dependency,
            string matrixXml)
        {
            string xml = string.Format(
                System.Globalization.CultureInfo.InvariantCulture,
                CompetingRisksCompatibilityTemplate,
                dependency,
                matrixXml);

            var risks = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(XElement.Parse(xml));

            Assert.AreEqual((Probability.DependencyType)Enum.Parse(typeof(Probability.DependencyType), dependency), risks.Dependency);
            Assert.IsNull(risks.CorrelationMatrix);
            Assert.AreEqual(Transform.Logarithmic, risks.XTransform);
            Assert.AreEqual(Transform.NormalZ, risks.ProbabilityTransform);
            Assert.IsFalse(risks.MinimumOfRandomVariables);
            Assert.AreEqual(24680, risks.PRNGSeed);
            Assert.AreEqual(UnivariateDistributionType.Normal, risks.Distributions[0].Type);
            Assert.AreEqual(UnivariateDistributionType.Gumbel, risks.Distributions[1].Type);
            CollectionAssert.AreEqual(new[] { 10d, 2d, 14d, 3d }, risks.GetParameters);
        }

        /// <summary>
        /// A correlation-dependent legacy import without configured correlation data remains
        /// unusable for Gaussian-copula simulation until the caller supplies a valid matrix.
        /// </summary>
        [TestMethod]
        public void Version2CompetingRisksXml_CorrelationDependencyWithoutMatrix_FailsWhenNumericallyUsed()
        {
            string xml = string.Format(
                System.Globalization.CultureInfo.InvariantCulture,
                CompetingRisksCompatibilityTemplate,
                "CorrelationMatrix",
                "<CorrelationMatrix />");
            var risks = (CompetingRisks)UnivariateDistributionFactory.CreateDistribution(XElement.Parse(xml));

            ArgumentException exception = Assert.Throws<ArgumentException>(() => risks.GenerateRandomValues(1, 12345));
            StringAssert.Contains(exception.Message, "requires a correlation matrix");
        }

        /// <summary>
        /// Populated correlation matrices that are truncated, out of range, asymmetric, or lack
        /// a unit diagonal remain invalid and cannot be mistaken for legacy placeholders.
        /// </summary>
        /// <param name="matrixXml">The malformed populated matrix.</param>
        [TestMethod]
        [DataRow("<CorrelationMatrix><Correlation_Row>1|0.2</Correlation_Row><Correlation_Row>0.2</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>1|1.2</Correlation_Row><Correlation_Row>1.2|1</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>1|0.2</Correlation_Row><Correlation_Row>0.3|1</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>0.9|0.2</Correlation_Row><Correlation_Row>0.2|1</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>0|0</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>0</Correlation_Row><Correlation_Row>0</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row>1|NaN</Correlation_Row><Correlation_Row>NaN|1</Correlation_Row></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix>not a matrix</CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Unexpected>0|0</Unexpected></CorrelationMatrix>")]
        [DataRow("<CorrelationMatrix><Correlation_Row><Unexpected>0|0</Unexpected></Correlation_Row><Correlation_Row>0|0</Correlation_Row></CorrelationMatrix>")]
        public void Version2CompetingRisksXml_PopulatedMalformedMatrix_RemainsRejected(string matrixXml)
        {
            string xml = string.Format(
                System.Globalization.CultureInfo.InvariantCulture,
                CompetingRisksCompatibilityTemplate,
                "Independent",
                matrixXml);

            Assert.Throws<ArgumentException>(() => UnivariateDistributionFactory.CreateDistribution(XElement.Parse(xml)));
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
