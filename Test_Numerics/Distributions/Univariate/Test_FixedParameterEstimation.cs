using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for <see cref="IFixedParameterEstimation"/> on <see cref="GeneralizedPareto"/>: the general form of the
    /// known-location fit.
    /// </summary>
    /// <remarks>
    /// The fixed-location maximum-likelihood estimate must set bitwise the parameters of the existing
    /// <see cref="GeneralizedPareto.MLE(IList{double}, double)"/>, and with no parameter fixed the estimate must be
    /// bitwise that of <see cref="GeneralizedPareto.Estimate(IList{double}, ParameterEstimationMethod)"/>. The SciPy
    /// reference is the one of <c>Test_GPA_MLE_FixedLocation_Reference</c> (scipy.stats.genpareto log-likelihood
    /// maximized at a fixed location of 0; tolerance from the Numerics Nelder-Mead convergence).
    /// </remarks>
    [TestClass]
    public class Test_FixedParameterEstimation
    {
        /// <summary>The excesses of the SciPy reference case.</summary>
        private static readonly double[] Excesses = { 1.2, 0.4, 3.9, 2.2, 0.15, 5.6, 1.05, 0.8, 7.3, 2.9, 0.33, 1.7, 4.4, 0.62, 2.5, 11.2, 0.09, 1.33, 3.1, 0.95 };

        /// <summary>Asserts two parameter vectors are bitwise equal.</summary>
        private static void AssertBitwise(double[] expected, double[] actual, string label)
        {
            Assert.HasCount(expected.Length, actual);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[i]), BitConverter.DoubleToInt64Bits(actual[i]), $"{label}[{i}]");
        }

        /// <summary>A fixed location with maximum likelihood sets bitwise the parameters of MLE(sample, location), and matches SciPy.</summary>
        [TestMethod]
        public void FixedLocationMaximumLikelihoodEqualsTheKnownLocationFit()
        {
            var random = new Random(31);
            foreach (double location in new[] { 0d, 10d, 50400d })
            {
                var sample = new double[60];
                var source = new GeneralizedPareto(location, 7d, -0.1);
                for (int i = 0; i < sample.Length; i++) sample[i] = source.InverseCDF(random.NextDouble());
                var gpa = new GeneralizedPareto();
                ((IFixedParameterEstimation)gpa).Estimate(sample, ParameterEstimationMethod.MaximumLikelihood, new[] { location, double.NaN, double.NaN });
                AssertBitwise(new GeneralizedPareto().MLE(sample, location), gpa.GetParameters, $"location {location}");
            }

            var fit = new GeneralizedPareto();
            fit.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, new[] { 0d, double.NaN, double.NaN });
            Assert.AreEqual(0d, fit.Xi);
            Assert.AreEqual(2.388710825197821, fit.Alpha, 1E-3 * 2.388710825197821);
            Assert.AreEqual(-0.07721168825622311, fit.Kappa, 1E-3);
            Assert.AreEqual(-38.959310020261206, fit.LogLikelihood(Excesses), 1E-6);
        }

        /// <summary>With no parameter fixed, every method gives bitwise the two-argument estimate.</summary>
        [TestMethod]
        public void NoFixedParameterIsTheOrdinaryEstimate()
        {
            var sample = new GeneralizedPareto(5d, 3d, 0.1).GenerateRandomValues(80, 12345);
            var heavy = new GeneralizedPareto(5d, 3d, -0.1).GenerateRandomValues(80, 777);
            foreach (var data in new[] { sample, heavy })
            foreach (var method in new[] { ParameterEstimationMethod.MethodOfMoments, ParameterEstimationMethod.MethodOfLinearMoments, ParameterEstimationMethod.MaximumLikelihood })
            {
                var ordinary = new GeneralizedPareto();
                var general = new GeneralizedPareto();
                Exception failure = null!;
                try { ordinary.Estimate(data, method); }
                catch (Exception e) { failure = e; }
                if (failure != null)
                {
                    // The ordinary estimate fails for this sample (for example the moment solve); the general form fails alike.
                    var again = Assert.Throws<Exception>(() => general.Estimate(data, method, new[] { double.NaN, double.NaN, double.NaN }));
                    Assert.AreEqual(failure.GetType(), again.GetType(), method.ToString());
                    continue;
                }
                general.Estimate(data, method, new[] { double.NaN, double.NaN, double.NaN });
                AssertBitwise(ordinary.GetParameters, general.GetParameters, method.ToString());
            }
        }

        /// <summary>Unsupported combinations and malformed fixed parameters are refused; the known-location checks still apply.</summary>
        [TestMethod]
        public void Refusals()
        {
            var gpa = new GeneralizedPareto();
            double nan = double.NaN;
            Assert.ThrowsExactly<ArgumentNullException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, null!));
            Assert.ThrowsExactly<ArgumentException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, new[] { 0d, nan }));
            Assert.ThrowsExactly<ArgumentException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, new[] { double.PositiveInfinity, nan, nan }));
            Assert.ThrowsExactly<NotSupportedException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MethodOfLinearMoments, new[] { 0d, nan, nan }));
            Assert.ThrowsExactly<NotSupportedException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, new[] { 0d, 2d, nan }));
            Assert.ThrowsExactly<NotSupportedException>(() => gpa.Estimate(Excesses, ParameterEstimationMethod.MaximumLikelihood, new[] { nan, nan, 0.1 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gpa.Estimate(new[] { 1d, -0.5 }, ParameterEstimationMethod.MaximumLikelihood, new[] { 0d, nan, nan }));
            Assert.IsInstanceOfType<IEstimation>(gpa);
        }
    }
}
