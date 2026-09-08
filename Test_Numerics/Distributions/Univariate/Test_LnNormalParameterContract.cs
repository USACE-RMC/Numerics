using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Physical-parameter preservation across conversion, cloning and serialization.</summary>
    [TestClass]
    public class Test_LnNormalParameterContract
    {
        /// <summary>Persisted physical moments are not repeatedly rounded through logarithmic coordinates.</summary>
        [TestMethod]
        public void PhysicalParameterVectorSurvivesRoundTripExactly()
        {
            double[] parameters = { 9.9999999999999982, 9.9999999999999964 };
            var distribution = new LnNormal(parameters[0], parameters[1]);
            CollectionAssert.AreEqual(parameters, distribution.GetParameters);
            CollectionAssert.AreEqual(parameters, distribution.Clone().GetParameters);
            var clone = new LnNormal();
            clone.SetParameters(distribution.GetParameters);
            CollectionAssert.AreEqual(parameters, clone.GetParameters);
        }

        /// <summary>Direct mutations of log coordinates invalidate the preserved physical representation.</summary>
        [TestMethod]
        public void LogParameterMutationRefreshesPhysicalMoments()
        {
            var distribution = new LnNormal(10, 2);
            distribution.Mu = 0;
            distribution.Sigma = 1;
            Assert.AreEqual(Math.Exp(.5), distribution.Mean, 2E-15);
            Assert.AreEqual(Math.Sqrt(Math.E * (Math.E - 1)), distribution.StandardDeviation, 2E-15);
        }
    }
}
