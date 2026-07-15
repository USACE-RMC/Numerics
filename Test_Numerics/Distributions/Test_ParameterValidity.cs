using System;
using System.Collections.Generic;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions
{
    /// <summary>
    /// Regression tests for parameter-validity state across univariate, multivariate, and bivariate distributions.
    /// </summary>
    [TestClass]
    public class Test_ParameterValidity
    {
        /// <summary>
        /// Verifies every factory-supported flat univariate parameter vector transitions from valid to invalid and back.
        /// </summary>
        [TestMethod]
        public void UnivariateFlatSettersTrackFinalParameterTuple()
        {
            foreach (UnivariateDistributionType type in Enum.GetValues(typeof(UnivariateDistributionType)))
            {
                if (!SupportsFlatParameterReplacement(type)) continue;

                var distribution = UnivariateDistributionFactory.CreateDistribution(type);
                Assert.IsTrue(distribution.ParametersValid, $"{type} default parameters must be valid.");
                double[] validParameters = distribution.GetParameters.ToArray();
                double[] invalidFirst = validParameters.ToArray();
                invalidFirst[0] = double.NaN;
                distribution.SetParameters(invalidFirst);
                Assert.IsFalse(distribution.ParametersValid, $"{type} must reject NaN in its first parameter.");

                invalidFirst[0] = double.PositiveInfinity;
                distribution.SetParameters(invalidFirst);
                Assert.IsFalse(distribution.ParametersValid, $"{type} must reject infinity in its first parameter.");
                distribution.SetParameters(validParameters);
                Assert.IsTrue(distribution.ParametersValid, $"{type} must recover after restoring its first parameter.");

                double[] invalidFinal = validParameters.ToArray();
                invalidFinal[invalidFinal.Length - 1] = double.NaN;
                distribution.SetParameters(invalidFinal);
                Assert.IsFalse(distribution.ParametersValid, $"{type} must reject NaN in its final parameter.");

                invalidFinal[invalidFinal.Length - 1] = double.PositiveInfinity;
                distribution.SetParameters(invalidFinal);
                Assert.IsFalse(distribution.ParametersValid, $"{type} must reject infinity in its final parameter.");
                distribution.SetParameters(validParameters);
                Assert.IsTrue(distribution.ParametersValid, $"{type} must recover after restoring its final parameter.");
                CollectionAssert.AreEqual(validParameters, distribution.GetParameters, $"{type} must store the restored tuple.");
            }
        }

        /// <summary>
        /// Verifies data-backed and composite univariate distributions refresh validity after every replacement path.
        /// </summary>
        [TestMethod]
        public void SpecializedUnivariateSettersTrackFinalState()
        {
            var empirical = new EmpiricalDistribution(new[] { 0d, 1d }, new[] { 0d, 1d });
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 2d });
            Assert.IsFalse(empirical.ParametersValid);
            Assert.Throws<ArgumentOutOfRangeException>(() => empirical.CDF(0.5d));
            Assert.Throws<ArgumentException>(() => empirical.SetParameters(new[] { 0d }, new[] { 0d, 1d }));
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 1d });
            Assert.IsTrue(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, 0d, 1d }, new[] { 0d, 0.5d, 1d });
            Assert.IsTrue(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, -1d }, new[] { 0d, 1d });
            Assert.IsFalse(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 1d });
            Assert.IsTrue(empirical.ParametersValid);


            var kernel = new KernelDensity(new[] { -1d, 0d, 1d });
            kernel.Bandwidth = double.NaN;
            Assert.IsFalse(kernel.ParametersValid);
            kernel.Bandwidth = double.PositiveInfinity;
            Assert.IsFalse(kernel.ParametersValid);
            kernel.Bandwidth = 0.5d;
            Assert.IsTrue(kernel.ParametersValid);

            var mixture = new Mixture(
                new[] { 0.8d, 0.8d },
                new UnivariateDistributionBase[] { new Normal(), new Exponential() });
            Assert.IsFalse(mixture.ParametersValid);
            mixture.SetParameters(
                new[] { 0.5d, 0.5d },
                new UnivariateDistributionBase[] { new Normal(), new Exponential() });
            Assert.IsTrue(mixture.ParametersValid);
            mixture.SetParameters(new[] { 0.5d, 0.5d }, new[] { 0d, -1d, 0d, 1d });
            Assert.IsFalse(mixture.ParametersValid);

            var emptyMixture = new Mixture(Array.Empty<double>(), Array.Empty<UnivariateDistributionBase>());
            Assert.IsFalse(emptyMixture.ParametersValid);

            var competingRisks = new CompetingRisks(Array.Empty<UnivariateDistributionBase>());
            var zeroInflatedMixture = new Mixture(
                new[] { 0.4d, 0.5d },
                new UnivariateDistributionBase[] { new Normal(), new Exponential() });
            zeroInflatedMixture.IsZeroInflated = true;
            zeroInflatedMixture.ZeroWeight = 0.1d;
            Assert.IsTrue(zeroInflatedMixture.ParametersValid);
            zeroInflatedMixture.IsZeroInflated = false;
            Assert.IsFalse(zeroInflatedMixture.ParametersValid);
            Assert.IsFalse(competingRisks.ParametersValid);
            competingRisks.SetParameters(new UnivariateDistributionBase[] { new Normal(0d, -1d) });
            Assert.IsFalse(competingRisks.ParametersValid);
            competingRisks.SetParameters(new UnivariateDistributionBase[] { new Normal() });
            Assert.IsTrue(competingRisks.ParametersValid);
            Assert.Throws<ArgumentException>(() => competingRisks.SetParameters(new[] { 0d }));

            var truncated = new TruncatedDistribution(new Normal(), -1d, 1d);
            Assert.IsTrue(truncated.ParametersValid);
            truncated.SetParameters(new[] { 0d, 1d, 2d, 1d });
            Assert.IsFalse(truncated.ParametersValid);
            Assert.Throws<ArgumentOutOfRangeException>(() => truncated.CDF(1.5d));
            truncated.SetParameters(new[] { 2d, 3d, -1d, 1d });
            Assert.IsTrue(truncated.ParametersValid);
            CollectionAssert.AreEqual(new[] { 2d, 3d, -1d, 1d }, truncated.GetParameters);
            truncated.SetParameters(new[] { 2d, 3d, 100d, 101d });
            Assert.Throws<ArgumentException>(() => truncated.SetParameters(new[] { 0d, 1d, -1d }));
            Assert.IsFalse(truncated.ParametersValid);
        }

        /// <summary>
        /// Verifies multivariate setters reject non-finite values before mutation and bivariate empirical data can recover.
        /// </summary>
        [TestMethod]
        public void MultivariateSettersRejectNonFiniteValuesTransactionally()
        {
            var covariance = new[,] { { 1d, 0.25d }, { 0.25d, 2d } };
            var normal = new MultivariateNormal(new[] { 1d, 2d }, covariance);
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.SetParameters(new[] { double.NaN, 2d }, covariance));
            CollectionAssert.AreEqual(new[] { 1d, 2d }, normal.Mean);
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.SetParameters(new[] { 1d, 2d }, new[,] { { 1d, double.PositiveInfinity }, { 0.25d, 2d } }));
            CollectionAssert.AreEqual(new[] { 1d, 2d }, normal.Mean);
            Assert.AreEqual(0.25d, normal.Covariance[0, 1]);
            Assert.Throws<ArgumentOutOfRangeException>(() => normal.SetParameters(new[] { 1d }, covariance));
            CollectionAssert.AreEqual(new[] { 1d, 2d }, normal.Mean);
            normal.SetParameters(new[] { 3d, 4d }, covariance);
            CollectionAssert.AreEqual(new[] { 3d, 4d }, normal.Mean);

            var student = new MultivariateStudentT(5d, new[] { 1d, 2d }, covariance);
            Assert.Throws<ArgumentOutOfRangeException>(() => student.SetParameters(double.NaN, new[] { 0d, 0d }, covariance));
            Assert.AreEqual(5d, student.DegreesOfFreedom);
            Assert.Throws<ArgumentOutOfRangeException>(() => student.SetParameters(5d, new[] { 0d, double.PositiveInfinity }, covariance));
            Assert.AreEqual(5d, student.DegreesOfFreedom);
            CollectionAssert.AreEqual(new[] { 1d, 2d }, student.Location);
            Assert.Throws<ArgumentOutOfRangeException>(() => student.SetParameters(5d, new[] { 0d }, covariance));
            CollectionAssert.AreEqual(new[] { 1d, 2d }, student.Location);
            student.SetParameters(7d, new[] { 3d, 4d }, covariance);
            Assert.AreEqual(7d, student.DegreesOfFreedom);
            CollectionAssert.AreEqual(new[] { 3d, 4d }, student.Location);

            var probabilities = new[,] { { 0.1d, 0.2d }, { 0.3d, 0.8d } };
            var empirical = new BivariateEmpirical(new[] { 0d, 1d }, new[] { 0d, 1d }, probabilities);
            empirical.SetParameters(new[] { 0d, double.NaN }, new[] { 0d, 1d }, probabilities);
            Assert.IsFalse(empirical.ParametersValid);
            Assert.Throws<ArgumentOutOfRangeException>(() => empirical.CDF(0.5d, 0.5d));
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, double.NegativeInfinity }, probabilities);
            Assert.IsFalse(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 1d }, new[,] { { 0.1d, double.PositiveInfinity }, { 0.3d, 0.8d } });
            Assert.IsFalse(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 1d }, new[,] { { 0.1d }, { 0.3d } });
            Assert.IsFalse(empirical.ParametersValid);
            empirical.SetParameters(new[] { 0d, 1d }, new[] { 0d, 1d }, probabilities);
            Assert.IsTrue(empirical.ParametersValid);

            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(new[] { 1d, double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Dirichlet(new[] { 1d, double.PositiveInfinity }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(2, new[] { 0.5d, double.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => new Multinomial(2, new[] { 0.5d, double.PositiveInfinity }));
        }

        /// <summary>
        /// Verifies every copula rejects non-finite dependency parameters and recovers when restored.
        /// </summary>
        [TestMethod]
        public void CopulasRejectNonFiniteParametersAndRecover()
        {
            var copulas = new BivariateCopula[]
            {
                new AMHCopula(0.5d),
                new ClaytonCopula(2d),
                new FrankCopula(4d),
                new GumbelCopula(2d),
                new JoeCopula(2d),
                new NormalCopula(0.5d),
                new StudentTCopula(0.5d, 5d)
            };

            foreach (BivariateCopula copula in copulas)
            {
                double validTheta = copula.Theta;
                copula.Theta = double.NaN;
                Assert.IsFalse(copula.ParametersValid, $"{copula.Type} must reject NaN.");
                copula.Theta = double.PositiveInfinity;
                Assert.IsFalse(copula.ParametersValid, $"{copula.Type} must reject infinity.");
                copula.Theta = validTheta;
                Assert.IsTrue(copula.ParametersValid, $"{copula.Type} must recover after restoring theta.");
            }

            var student = new StudentTCopula(0.5d, 5d);
            student.DegreesOfFreedom = double.PositiveInfinity;
            Assert.IsFalse(student.ParametersValid);
            student.DegreesOfFreedom = 5d;
            Assert.IsTrue(student.ParametersValid);
        }

        /// <summary>
        /// Determines whether a factory-created distribution supports flattened numeric parameter replacement.
        /// </summary>
        /// <param name="type">The univariate distribution type.</param>
        /// <returns><see langword="true"/> when <see cref="UnivariateDistributionBase.SetParameters"/> is supported without external components.</returns>
        private static bool SupportsFlatParameterReplacement(UnivariateDistributionType type)
        {
            return type != UnivariateDistributionType.CompetingRisks &&
                   type != UnivariateDistributionType.Empirical &&
                   type != UnivariateDistributionType.KernelDensity &&
                   type != UnivariateDistributionType.Mixture &&
                   type != UnivariateDistributionType.UserDefined;
        }
    }
}
