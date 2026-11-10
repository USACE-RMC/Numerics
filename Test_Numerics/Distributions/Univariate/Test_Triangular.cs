using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Triangular distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Triangular
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_TriangularDist()
        {
            double true_mean = 3.3333333333333335d;
            double true_median = 3.2613872124741694d;
            double true_mode = 3.0d;
            double true_stdDev = Math.Sqrt(1.0555555555555556d);
            double true_pdf = 0.20000000000000001d;
            double true_cdf = 0.10000000000000001d;
            double true_icdf = 2.0d;
            var T = new Triangular(1, 3, 6);

            Assert.AreEqual(true_mean, T.Mean, 0.0001d);
            Assert.AreEqual(true_median, T.Median, 0.0001d);
            Assert.AreEqual(true_mode, T.Mode, 0.0001d);
            Assert.AreEqual(true_stdDev, T.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_pdf, T.PDF(2.0d), 0.0001d);
            Assert.AreEqual(true_cdf, T.CDF(2.0d), 0.0001d);
            Assert.AreEqual(true_icdf, T.InverseCDF(true_cdf), 0.0001d);
        }

        /// <summary>
        /// Verified using R 'mc2d'
        /// </summary>
        [TestMethod]
        public void Test_Triangular_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.0400000, 0.2500000, 0.6666667, 0.9166667, 0.9866667 };
            var tri = new Triangular(0, 0.25, 1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = tri.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], tri.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// Test estimation with moments.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_MOM()
        {
            var dist = new Triangular(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 1);
        }

        /// <summary>
        /// Test estimation with likelihood.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_MLE()
        {
            var dist = new Triangular(-2, 10, 35);
            var sample = dist.GenerateRandomValues(1000, 12345);
            dist.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(-2, dist.Min, 1);
            Assert.AreEqual(10, dist.MostLikely, 1);
            Assert.AreEqual(35, dist.Max, 1);
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.Min);
            Assert.AreEqual(0.5, T.Mode);
            Assert.AreEqual(1, T.Max);

            var T2 = new Triangular(-1,1,2);
            Assert.AreEqual(-1, T2.Min);
            Assert.AreEqual(1, T2.Mode);
            Assert.AreEqual(2, T2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var T = new Triangular(double.NaN,double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(T.ParametersValid);

            var T2 = new Triangular(double.PositiveInfinity,double.NaN,double.NaN);
            Assert.IsFalse(T2.ParametersValid);

            var T3 = new Triangular(4, 1, 0);
            Assert.IsFalse(T3.ParametersValid);

            var T4 = new Triangular(1, 0, -1);
            Assert.IsFalse(T4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var T = new Triangular();
            Assert.AreEqual("Min (a)", T.ParametersToString[0, 0]);
            Assert.AreEqual("Most Likely (c)", T.ParametersToString[1, 0]);
            Assert.AreEqual("Max (b)", T.ParametersToString[2, 0]);
            Assert.AreEqual("0", T.ParametersToString[0, 1]);
            Assert.AreEqual("0.5", T.ParametersToString[1, 1]);
            Assert.AreEqual("1", T.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Triangular(1, 3, 6);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var T = new Triangular();
            Assert.AreEqual(0.5, T.Mean);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(3.3333, T2.Mean, 1e-04);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var T = new Triangular();
            Assert.AreEqual(0.5, T.Median);

            var T2 = new Triangular(1,3,6);
            Assert.AreEqual(3.26138, T2.Median, 1e-05);
        }

        /// <summary>
        /// Testing mode
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var T = new Triangular();
            Assert.AreEqual(0.5, T.Mode);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(3, T2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            var T = new Triangular();
            Assert.AreEqual(0.20412, T.StandardDeviation, 1e-04);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(1.02739, T2.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.Skewness);
        }

        /// <summary>
        /// Testing kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var T = new Triangular();
            Assert.AreEqual(12d / 5d, T.Kurtosis);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(12d / 5d, T2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.Minimum);
            Assert.AreEqual(1, T.Maximum);

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(1, T2.Minimum);
            Assert.AreEqual(6, T2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.PDF(-1));
            Assert.AreEqual(1.6, T.PDF(0.4));
            Assert.AreEqual(1.6, T.PDF(0.6));
            Assert.AreEqual(2, T.PDF(0.5));

            var T2 = new Triangular(1, 3, 6);
            Assert.AreEqual(0.2, T2.PDF(2), 1e-04);
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.CDF(-1));
            Assert.AreEqual(1, T.CDF(2));
            Assert.AreEqual(0.32, T.CDF(0.4), 1e-04);
            Assert.AreEqual(0.68, T.CDF(0.6), 1e-04);

            var T2 = new Triangular(1,3, 6);
            Assert.AreEqual(0.1, T2.CDF(2), 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var T = new Triangular();
            Assert.AreEqual(0, T.InverseCDF(0));
            Assert.AreEqual(1, T.InverseCDF(1));
            Assert.AreEqual(0.31622, T.InverseCDF(0.2), 1e-04);
            Assert.AreEqual(0.5, T.InverseCDF(0.5));
        }

        /// <summary>
        /// MLE succeeds when the sample minimum is exactly zero, where the established bounds collapse to zero.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_MLE_ZeroSampleMinimum()
        {
            var sample = new double[] { 0, 1, 2, 2.5, 3, 4, 5, 6, 8, 9 };
            var constraints = new Triangular().GetParameterConstraints(sample);
            // S = max(|0|, |9|, 9) = 9, bound = 10^(ceil(log10 9) + 1) = 100
            CollectionAssert.AreEqual(new double[] { -100, 0, 9 }, constraints.Item2);
            CollectionAssert.AreEqual(new double[] { 0, 9, 100 }, constraints.Item3);

            var fit = new Triangular();
            fit.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.IsTrue(fit.ParametersValid);
            // Multi-start optimum from scipy.optimize.minimize(Nelder-Mead, 200 starts) of the negative
            // log-likelihood sum(log(scipy.stats.triang((c-a)/(b-a), loc=a, scale=b-a).pdf(x))):
            // a = 0, c = 0, b = 11.2239153, log-likelihood -22.762864272955646.
            Assert.AreEqual(-22.762864272955646, fit.LogLikelihood(sample), 1E-6);
            Assert.AreEqual(0d, fit.Min, 1E-6);
            Assert.AreEqual(11.2239153, fit.Max, 1E-2);
        }

        /// <summary>
        /// Samples whose established bounds contain the initial values keep exactly those constraints.
        /// </summary>
        [TestMethod]
        public void Test_Triangular_ParameterConstraints_NonZeroSampleMatchesFormula()
        {
            var sample = new double[] { 1, 2, 3, 3.5, 4, 5, 6, 7, 9, 10 };
            double min = 1, max = 10, mean = Numerics.Tools.Mean(sample);
            var constraints = new Triangular().GetParameterConstraints(sample);
            CollectionAssert.AreEqual(new double[] { min - Numerics.Tools.DoubleMachineEpsilon, mean * 3 - max - min, max + Numerics.Tools.DoubleMachineEpsilon }, constraints.Item1);
            CollectionAssert.AreEqual(new double[] { -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(min)) + 1d)), min, max }, constraints.Item2);
            CollectionAssert.AreEqual(new double[] { min, max, Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(max)) + 1d)) }, constraints.Item3);
        }
    }
}
