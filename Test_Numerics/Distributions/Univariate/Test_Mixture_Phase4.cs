using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Regression tests for the corrected physical-simplex and positive-hurdle mixture behavior.
    /// </summary>
    [TestClass]
    public class Test_Mixture_Phase4
    {
        /// <summary>
        /// Verifies that every public setter family copies caller-owned arrays and lists.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_Setters_DoNotMutateOrAliasCallerArrays()
        {
            var constructorWeights = new[] { 0.4, 0.6 };
            var constructorDistributions = new UnivariateDistributionBase[]
            {
                new Normal(0.0, 1.0),
                new Normal(3.0, 1.0)
            };
            var mixture = new Mixture(constructorWeights, constructorDistributions);
            constructorWeights[0] = 0.9;
            constructorDistributions[0] = new Normal(10.0, 1.0);
            Assert.AreEqual(0.4, mixture.Weights[0], 0.0);
            Assert.AreEqual(0.0, mixture.Distributions[0].Mean, 0.0);

            var baseWeights = new[] { 0.25, 0.75 };
            var baseDistributions = new UnivariateDistributionBase[]
            {
                new Normal(1.0, 1.0),
                new Normal(4.0, 1.0)
            };
            mixture.SetParameters(baseWeights, baseDistributions);
            baseWeights[0] = 0.8;
            baseDistributions[0] = new Normal(20.0, 1.0);
            Assert.AreEqual(0.25, mixture.Weights[0], 0.0);
            Assert.AreEqual(1.0, mixture.Distributions[0].Mean, 0.0);

            var interfaceWeights = new[] { 0.3, 0.7 };
            var interfaceDistributions = new IUnivariateDistribution[]
            {
                new Normal(2.0, 1.0),
                new Normal(5.0, 1.0)
            };
            mixture.SetParameters(interfaceWeights, interfaceDistributions);
            interfaceWeights[0] = 0.6;
            interfaceDistributions[0] = new Normal(30.0, 1.0);
            Assert.AreEqual(0.3, mixture.Weights[0], 0.0);
            Assert.AreEqual(2.0, mixture.Distributions[0].Mean, 0.0);

            var parameterWeights = new[] { 0.2, 0.8 };
            var distributionParameters = new[] { 6.0, 1.0, 9.0, 2.0 };
            mixture.SetParameters(parameterWeights, distributionParameters);
            parameterWeights[0] = 0.5;
            distributionParameters[0] = 60.0;
            Assert.AreEqual(0.2, mixture.Weights[0], 0.0);
            Assert.AreEqual(6.0, mixture.Distributions[0].GetParameters[0], 0.0);

            var listParameters = new List<double> { 0.35, 0.65, 7.0, 1.0, 10.0, 2.0 };
            double[] listSnapshot = listParameters.ToArray();
            mixture.SetParameters(listParameters);
            CollectionAssert.AreEqual(listSnapshot, listParameters);
            listParameters[0] = 0.9;
            Assert.AreEqual(0.35, mixture.Weights[0], 0.0);

            var referencedParameters = new[] { 2.0, 6.0, 8.0, 1.0, 11.0, 2.0 };
            double[] referencedSnapshot = referencedParameters.ToArray();
            mixture.SetParameters(ref referencedParameters);
            CollectionAssert.AreEqual(referencedSnapshot, referencedParameters);
            Assert.AreEqual(0.25, mixture.Weights[0], 1E-15);
            Assert.AreEqual(0.75, mixture.Weights[1], 1E-15);
        }

        /// <summary>
        /// Verifies ordinary and zero-inflated physical-simplex validation.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ValidatesConfiguredPhysicalSimplex()
        {
            var ordinary = new Mixture(
                new[] { 0.4, 0.4 },
                new UnivariateDistributionBase[] { new Normal(0.0, 1.0), new Normal(3.0, 1.0) });
            Assert.IsFalse(ordinary.ParametersValid);
            Assert.Throws<ArgumentOutOfRangeException>(() => ordinary.PDF(1.0));

            var zeroInflated = new Mixture(
                new[] { 0.5, 0.5 },
                new UnivariateDistributionBase[] { new Normal(0.0, 1.0), new Normal(3.0, 1.0) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.1
            };
            Assert.IsTrue(zeroInflated.ParametersValid);
            zeroInflated.SetParameters(
                new[] { 0.4, 0.4 },
                new UnivariateDistributionBase[] { new Normal(0.0, 1.0), new Normal(3.0, 1.0) });
            Assert.IsFalse(zeroInflated.ParametersValid);
            Assert.Throws<ArgumentOutOfRangeException>(() => zeroInflated.CDF(1.0));

            zeroInflated.ZeroWeight = 1.0;
            Assert.IsFalse(zeroInflated.ParametersValid);
        }

        /// <summary>
        /// Verifies the analytical positive-hurdle Normal density, distribution, log, and quantile identities.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ZeroInflatedNormal_UsesPositiveHurdleIdentities()
        {
            var normal = new Normal(0.0, 1.0);
            var mixture = new Mixture(new[] { 1.0 }, new UnivariateDistributionBase[] { normal })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.2
            };

            const double x = 1.3;
            double positiveMass = normal.CCDF(0.0);
            double expectedPdf = 0.8 * normal.PDF(x) / positiveMass;
            double expectedCdf = 0.2 + 0.8 * (normal.CDF(x) - normal.CDF(0.0)) / positiveMass;
            double expectedCcdf = 0.8 * normal.CCDF(x) / positiveMass;

            Assert.AreEqual(expectedPdf, mixture.PDF(x), 1E-14);
            Assert.AreEqual(Math.Log(expectedPdf), mixture.LogPDF(x), 1E-14);
            Assert.AreEqual(expectedCdf, mixture.CDF(x), 1E-14);
            Assert.AreEqual(Math.Log(expectedCdf), mixture.LogCDF(x), 1E-14);
            Assert.AreEqual(Math.Log(expectedCcdf), mixture.LogCCDF(x), 1E-14);

            const double probability = 0.6;
            double conditionalProbability = (probability - 0.2) / 0.8;
            double expectedQuantile = normal.InverseCDF(normal.CDF(0.0) + conditionalProbability * positiveMass);
            Assert.AreEqual(expectedQuantile, mixture.InverseCDF(probability), 1E-6);
        }

        /// <summary>
        /// Verifies the atom at zero and absence of negative support under the hurdle model.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ZeroInflatedModel_HasZeroJumpAndNoNegativeSupport()
        {
            var mixture = new Mixture(
                new[] { 1.0 },
                new UnivariateDistributionBase[] { new Normal(0.0, 1.0) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.2
            };

            Assert.AreEqual(0.0, mixture.PDF(-1.0), 0.0);
            Assert.AreEqual(double.NegativeInfinity, mixture.LogPDF(-1.0));
            Assert.AreEqual(0.0, mixture.CDF(-1.0), 0.0);
            Assert.AreEqual(0.2, mixture.PDF(0.0), 0.0);
            Assert.AreEqual(Math.Log(0.2), mixture.LogPDF(0.0), 0.0);
            Assert.AreEqual(0.2, mixture.CDF(0.0), 0.0);
            Assert.AreEqual(Math.Log(0.8), mixture.LogCCDF(0.0), 1E-15);
            Assert.AreEqual(0.0, mixture.InverseCDF(0.2), 0.0);
        }

        /// <summary>
        /// Verifies simulation uses the atom and strictly positive conditioned components.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ZeroInflatedSimulation_UsesAtomAndPositiveConditioning()
        {
            var mixture = new Mixture(
                new[] { 1.0 },
                new UnivariateDistributionBase[] { new Normal(0.0, 1.0) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.2
            };

            double[] sample = mixture.GenerateRandomValues(20000, 12345);
            double atomFrequency = sample.Count(value => value == 0.0) / (double)sample.Length;

            Assert.AreEqual(0.2, atomFrequency, 0.015);
            Assert.IsFalse(sample.Any(value => value < 0.0));
            Assert.IsTrue(sample.Any(value => value > 0.0));
        }

        /// <summary>
        /// Verifies a hurdle component must have finite, nonzero probability above zero.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ZeroInflatedModel_RejectsComponentWithoutPositiveMass()
        {
            var mixture = new Mixture(
                new[] { 1.0 },
                new UnivariateDistributionBase[] { new Deterministic(0.0) });
            mixture.ZeroWeight = 0.2;
            mixture.IsZeroInflated = true;

            Assert.IsFalse(mixture.ParametersValid);
            ArgumentOutOfRangeException exception =
                Assert.Throws<ArgumentOutOfRangeException>(() => mixture.PDF(1.0));
            StringAssert.Contains(exception.Message, "positive probability above zero");
        }

        /// <summary>
        /// Verifies zero-inflated EM rejects negative exact observations with row context.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_ZeroInflatedEM_RejectsNegativeExactObservationWithRowContext()
        {
            var mixture = new Mixture(
                new[] { 1.0 },
                new UnivariateDistributionBase[] { new Normal(1.0, 1.0) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.1
            };

            InvalidOperationException exception =
                Assert.Throws<InvalidOperationException>(() => mixture.MLE(new[] { 1.0, -2.5, 2.0 }));
            StringAssert.Contains(exception.Message, "row 1");
            StringAssert.Contains(exception.Message, "-2.5");
        }

        /// <summary>
        /// Verifies EM reports an impossible row instead of silently continuing.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_EM_ImpossibleRowThrowsWithRowContext()
        {
            var mixture = new Mixture(
                new[] { 1.0 },
                new UnivariateDistributionBase[] { new Normal(1.0, 1.0) })
            {
                IsZeroInflated = true,
                ZeroWeight = 0.0
            };

            InvalidOperationException exception =
                Assert.Throws<InvalidOperationException>(() => mixture.MLE(new[] { 0.0, 1.0, 2.0 }));
            StringAssert.Contains(exception.Message, "row 0");
            StringAssert.Contains(exception.Message, "value 0");
            StringAssert.Contains(exception.Message, "zero or nonfinite");
        }
    }
}
