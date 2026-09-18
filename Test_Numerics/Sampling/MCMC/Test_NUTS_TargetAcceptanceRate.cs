using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for the NUTS dual-averaging target acceptance rate.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Raising the target acceptance rate is the standard response to divergent transitions, so the
    /// value has to be reachable. These tests pin the default at the Hoffman and Gelman value of 0.80,
    /// pin the accepted range to the open interval (0, 1), and check that a raised target does what it
    /// is raised for: a smaller adapted step size.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Hoffman, M.D. and Gelman, A. (2014). "The No-U-Turn Sampler: Adaptively Setting Path Lengths
    /// in Hamiltonian Monte Carlo." Journal of Machine Learning Research, 15, 1593-1623.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_NUTS_TargetAcceptanceRate
    {
        /// <summary>
        /// Builds a small, fully deterministic single-chain sampler on a standard normal target.
        /// </summary>
        /// <returns>An unsampled sampler.</returns>
        private static NUTS BuildSampler()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-50d, 50d), new Uniform(-50d, 50d) };

            static double logLH(double[] x) => -0.5d * (x[0] * x[0] + x[1] * x[1]);

            return new NUTS(priors, logLH, maxTreeDepth: 6)
            {
                NumberOfChains = 1,
                ParallelizeChains = false,
                InitialIterations = 1,
                ThinningInterval = 1,
                WarmupIterations = 60,
                Iterations = 150,
                OutputLength = 100,
                PRNGSeed = 12345
            };
        }

        /// <summary>
        /// The target acceptance rate defaults to the Hoffman and Gelman value of 0.80.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_TargetAcceptanceRate_DefaultsToPointEight()
        {
            var sampler = BuildSampler();
            Assert.AreEqual(0.80d, sampler.TargetAcceptanceRate, 0d,
                "The target acceptance rate default must remain 0.80 so that existing results are unchanged.");
        }

        /// <summary>
        /// The target acceptance rate is settable and reads back exactly what was written.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_TargetAcceptanceRate_IsSettable()
        {
            var sampler = BuildSampler();
            sampler.TargetAcceptanceRate = 0.95d;
            Assert.AreEqual(0.95d, sampler.TargetAcceptanceRate, 0d);
        }

        /// <summary>
        /// Values outside the open interval (0, 1), and non-finite values, are rejected when the
        /// sampler validates its settings.
        /// </summary>
        /// <remarks>
        /// The sampler validates configuration at the start of <see cref="MCMCSampler.Sample"/> rather
        /// than at assignment, which is the convention every other setting in this class and in
        /// <see cref="MCMCSampler"/> follows. Validation runs before any chain work, so these calls
        /// throw immediately.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_TargetAcceptanceRate_RejectsValuesOutsideTheUnitInterval()
        {
            double[] invalid =
            {
                0d, 1d, -0.5d, -1e-16d, 1.5d, double.NaN,
                double.PositiveInfinity, double.NegativeInfinity
            };

            foreach (double value in invalid)
            {
                var sampler = BuildSampler();
                sampler.TargetAcceptanceRate = value;
                Assert.ThrowsExactly<ArgumentException>(() => sampler.Sample(),
                    $"A target acceptance rate of {value:R} should have been rejected.");
            }
        }

        /// <summary>
        /// Values strictly inside (0, 1) are accepted.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_TargetAcceptanceRate_AcceptsValuesInsideTheUnitInterval()
        {
            double[] valid = { 1e-12d, 0.6d, 0.8d, 0.95d, 0.99d, 1d - 1e-12d };

            foreach (double value in valid)
            {
                var sampler = BuildSampler();
                sampler.TargetAcceptanceRate = value;
                sampler.Sample();
                Assert.AreEqual(value, sampler.TargetAcceptanceRate, 0d);
            }
        }

        /// <summary>
        /// Raising the target acceptance rate must shorten the adapted step size, which is the whole
        /// reason the setting exists.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_TargetAcceptanceRate_RaisingTheTargetShortensTheStepSize()
        {
            var baseline = BuildSampler();
            baseline.Sample();

            var raised = BuildSampler();
            raised.TargetAcceptanceRate = 0.99d;
            raised.Sample();

            Assert.IsLessThan(baseline.StepSizes[0], raised.StepSizes[0],
                $"Raising the target from 0.80 to 0.99 left the step size at {raised.StepSizes[0]:R}, " +
                $"which is not shorter than the default run's {baseline.StepSizes[0]:R}.");
        }
    }
}
