using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Guards bulk parameter update allocations and state for GEV and generalized Pareto fitting.</summary>
    [TestClass]
    public class Test_GeneralizedBulkSetterRegressions
    {
#if NET8_0_OR_GREATER
        /// <summary>Valid repeated updates reuse their inputs without constructing temporary validation lists.</summary>
        /// <param name="pareto">Whether to exercise generalized Pareto instead of GEV.</param>
        /// <param name="list">Whether to use the list overload instead of the scalar overload.</param>
        [TestMethod]
        [DataRow(false, false)]
        [DataRow(false, true)]
        [DataRow(true, false)]
        [DataRow(true, true)]
        public void ValidBulkUpdates_DoNotAllocate(bool pareto, bool list)
        {
            UnivariateDistributionBase distribution = Create(pareto);
            double[] parameters = { -2d, 1E-20, 0.25d };
            Action update = () => Set(distribution, parameters, list);
            for (int i = 0; i < 1000; i++) update();

            long before = GC.GetAllocatedBytesForCurrentThread();
            for (int i = 0; i < 1000; i++) update();
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.IsTrue(distribution.ParametersValid);
            CollectionAssert.AreEqual(parameters, distribution.GetParameters);
            Assert.AreEqual(0L, allocated, $"1000 bulk updates allocated {allocated} bytes ({allocated / 1000d} bytes/call).");
        }
#endif

        /// <summary>Bulk updates retain exact values, validation precedence, and recovery from every invalid parameter.</summary>
        /// <param name="pareto">Whether to exercise generalized Pareto instead of GEV.</param>
        /// <param name="list">Whether to use the list overload instead of the scalar overload.</param>
        [TestMethod]
        [DataRow(false, false)]
        [DataRow(false, true)]
        [DataRow(true, false)]
        [DataRow(true, true)]
        public void BulkUpdates_PreserveValuesValidityAndValidationExceptions(bool pareto, bool list)
        {
            UnivariateDistributionBase distribution = Create(pareto);
            double[][] invalid = {
                new[] { double.NaN, 2d, 0.25d },
                new[] { 1d, 0d, 0.25d },
                new[] { 1d, -2d, 0.25d },
                new[] { 1d, double.PositiveInfinity, 0.25d },
                new[] { 1d, 2d, double.NegativeInfinity },
                new[] { double.PositiveInfinity, -2d, double.NaN }
            };
            foreach (double[] parameters in invalid)
            {
                var expected = distribution.ValidateParameters(parameters, false);
                Assert.IsNotNull(expected);
                Set(distribution, parameters, list);
                Assert.IsFalse(distribution.ParametersValid);
                CollectionAssert.AreEqual(parameters, distribution.GetParameters);
                var actual = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.LogPDF(3d));
                Assert.AreEqual(expected.ParamName, actual.ParamName);
                Assert.AreEqual(expected.Message, actual.Message);
                Assert.AreEqual(expected.ActualValue, actual.ActualValue);

                double[] valid = { -0d, 1E-20, -0.25d };
                Set(distribution, valid, list);
                Assert.IsTrue(distribution.ParametersValid);
                double[] retained = distribution.GetParameters;
                for (int i = 0; i < valid.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(valid[i]), BitConverter.DoubleToInt64Bits(retained[i]));
                Assert.IsNull(distribution.ValidateParameters(retained, false));
            }
        }

        /// <summary>Rejected list lengths retain the existing parameter exception and leave distribution state unchanged.</summary>
        /// <param name="pareto">Whether to exercise generalized Pareto instead of GEV.</param>
        [TestMethod]
        [DataRow(false)]
        [DataRow(true)]
        public void InvalidListLength_PreservesExceptionAndState(bool pareto)
        {
            UnivariateDistributionBase distribution = Create(pareto);
            double[] original = distribution.GetParameters;
            foreach (double[] parameters in new double[][] { null, Array.Empty<double>(), new[] { 1d, 2d }, new[] { 1d, 2d, 3d, 4d } })
            {
                var expected = distribution.ValidateParameters(parameters, false);
                var actual = Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => distribution.SetParameters(parameters));
                Assert.AreEqual(expected.ParamName, actual.ParamName);
                Assert.AreEqual(expected.Message, actual.Message);
                CollectionAssert.AreEqual(original, distribution.GetParameters);
                Assert.IsTrue(distribution.ParametersValid);
            }
        }

        private static UnivariateDistributionBase Create(bool pareto)
        {
            return pareto ? new GeneralizedPareto(1d, 2d, 0.25d) : new GeneralizedExtremeValue(1d, 2d, 0.25d);
        }

        private static void Set(UnivariateDistributionBase distribution, double[] parameters, bool list)
        {
            if (list) distribution.SetParameters(parameters);
            else if (distribution is GeneralizedPareto pareto) pareto.SetParameters(parameters[0], parameters[1], parameters[2]);
            else ((GeneralizedExtremeValue)distribution).SetParameters(parameters[0], parameters[1], parameters[2]);
        }
    }
}
