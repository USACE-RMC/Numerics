using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Uniform distribution algorithm.
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
    public class Test_Uniform
    {

        /// <summary>
        /// Verified using Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_UniformDist()
        {
            double true_mean = 0.76d;
            double true_median = 0.76d;
            double true_stdDev = Math.Sqrt(0.03853333333333335d);
            double true_pdf = 1.4705882352941173d;
            double true_cdf = 0.70588235294117641d;
            double true_icdf = 0.9d;
            var U = new Uniform(0.42d, 1.1d);
            Assert.AreEqual(true_mean, U.Mean, 0.0001d);
            Assert.AreEqual(true_median, U.Median, 0.0001d);
            Assert.AreEqual(true_stdDev, U.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_pdf, U.PDF(0.9d), 0.0001d);
            Assert.AreEqual(true_cdf, U.CDF(0.9d), 0.0001d);
            Assert.AreEqual(true_icdf, U.InverseCDF(true_cdf), 0.0001d);
        }

        /// <summary>
        /// Verified using R 'stats'
        /// </summary>
        [TestMethod]
        public void Test_Uniform_R()
        {
            var x = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var p = new double[5];
            var true_p = new double[] { 0.1, 0.25, 0.5, 0.75, 0.9 };
            var u = new Uniform(0, 1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = u.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], u.InverseCDF(p[i]), 1E-7);
            }
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Min);
            Assert.AreEqual(1, U.Max);

            var U2 = new Uniform(2,10);
            Assert.AreEqual(2, U2.Min);
            Assert.AreEqual(10, U2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var U = new Uniform(1, 0);
            Assert.IsFalse(U.ParametersValid);

            var U2 = new Uniform(double.NaN, 0);
            Assert.IsFalse(U2.ParametersValid);

            var U3 = new Uniform(0,double.NaN);
            Assert.IsFalse(U3.ParametersValid);

            var U4 = new Uniform(0,double.PositiveInfinity);
            Assert.IsFalse(U4.ParametersValid);

            var U5 = new Uniform(double.PositiveInfinity, 0);
            Assert.IsFalse(U5.ParametersValid);
        }

        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var U = new Uniform();
            Assert.AreEqual("Min", U.ParametersToString[0, 0]);
            Assert.AreEqual("Max", U.ParametersToString[1, 0]);
            Assert.AreEqual("0", U.ParametersToString[0, 1]);
            Assert.AreEqual("1", U.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Uniform(2, 10);
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
            var U = new Uniform();
            Assert.AreEqual(0.5, U.Mean);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(6, U2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var U = new Uniform();
            Assert.AreEqual(0.5, U.Median);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(6, U2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var U = new Uniform();
            Assert.AreEqual(double.NaN, U.Mode);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(double.NaN, U2.Mode);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            var U = new Uniform();
            Assert.AreEqual(0.288675, U.StandardDeviation, 1e-05);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(2.3094, U2.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Skewness);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(0, U2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var U = new Uniform();
            Assert.AreEqual(9d / 5d, U.Kurtosis);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(9d / 5d, U2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.Minimum);
            Assert.AreEqual(1, U.Maximum);

            var U2 = new Uniform(2, 10);
            Assert.AreEqual(2, U2.Minimum);
            Assert.AreEqual(10, U2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.PDF(-1));
            Assert.AreEqual(0, U.PDF(2));
            Assert.AreEqual(1, U.PDF(1));
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.CDF(0));
            Assert.AreEqual(1, U.CDF(1));
            Assert.AreEqual(0.5, U.CDF(0.5));
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var U = new Uniform();
            Assert.AreEqual(0, U.InverseCDF(0));
            Assert.AreEqual(1, U.InverseCDF(1));
            Assert.AreEqual(0.3, U.InverseCDF(0.3));
        }

        /// <summary>
        /// For Min &lt; Max the PDF, CDF, quantile and mean are bit-identical to the established expressions,
        /// which are inlined here.
        /// </summary>
        [TestMethod]
        public void Test_OrdinaryParameters_BitIdentical()
        {
            (double min, double max)[] parameters = { (0, 1), (-3.7, 12.25), (1E-8, 2E-8), (-1E6, 1E6), (2, 2.0000001) };
            double[] fractions = { 0, 1E-12, 0.1, 0.3333333333333333, 0.5, 0.9, 0.999999, 1 };
            foreach (var (min, max) in parameters)
            {
                var U = new Uniform(min, max);
                Assert.AreEqual((min + max) / 2.0d, U.Mean, 0d);
                foreach (double f in fractions)
                {
                    double x = min + f * (max - min);
                    Assert.AreEqual(x < min || x > max ? 0d : 1.0d / (max - min), U.PDF(x), 0d);
                    double cdf = x <= min ? 0d : x >= max ? 1d : (x - min) / (max - min);
                    Assert.AreEqual(cdf, U.CDF(x), 0d);
                    double q = f == 0 ? min : f == 1 ? max : min + f * (max - min);
                    Assert.AreEqual(q, U.InverseCDF(f), 0d);
                }
                Assert.AreEqual(0d, U.PDF(min - 1));
                Assert.AreEqual(0d, U.CDF(min - 1));
                Assert.AreEqual(1d, U.CDF(max + 1));
            }
        }

        /// <summary>
        /// Min = Max is a valid point mass: the CDF steps from 0 to 1 at Min and every quantile is Min.
        /// </summary>
        [TestMethod]
        public void Test_PointMass()
        {
            var U = new Uniform(2, 2);
            Assert.IsTrue(U.ParametersValid);
            Assert.AreEqual(0d, U.CDF(1));
            Assert.AreEqual(1d, U.CDF(2));
            Assert.AreEqual(1d, U.CDF(3));
            Assert.AreEqual(0d, U.PDF(2));
            Assert.AreEqual(0d, U.PDF(1));
            Assert.AreEqual(2d, U.InverseCDF(0));
            Assert.AreEqual(2d, U.InverseCDF(0.4));
            Assert.AreEqual(2d, U.InverseCDF(1));
            Assert.AreEqual(2d, U.Mean);

            var Z = new Uniform(0, 0);
            Assert.AreEqual(0d, Z.InverseCDF(0.7));
            Assert.AreEqual(0d, Z.Mean);
        }

        /// <summary>
        /// A not-a-number argument returns not-a-number.
        /// </summary>
        [TestMethod]
        public void Test_NaNArgument()
        {
            Assert.IsTrue(double.IsNaN(new Uniform(0, 1).PDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(new Uniform(0, 1).CDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(new Uniform(2, 2).PDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(new Uniform(2, 2).CDF(double.NaN)));
        }

        /// <summary>
        /// Maximum likelihood gives (min, max); the method of moments gives mean ∓ √3·s.
        /// </summary>
        [TestMethod]
        public void Test_Estimation()
        {
            var data = new[] { 2.3, 5.1, 3.3, 4.8, 2.9, 4.1 };
            var U = new Uniform();
            U.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
            // scipy.stats.uniform.fit(data) = (loc=2.3, scale=2.8)
            Assert.AreEqual(2.3, U.Min, 0d);
            Assert.AreEqual(5.1, U.Max, 0d);
            U.Estimate(data, ParameterEstimationMethod.MethodOfMoments);
            // numpy: mean(data) -/+ sqrt(3)*std(data, ddof=1)
            Assert.AreEqual(1.8408116907963217, U.Min, 1E-14);
            Assert.AreEqual(5.659188309203678, U.Max, 1E-14);

            var moments = U.MomentsFromParameters(new[] { 0.0, 12.0 });
            Assert.AreEqual(6.0, moments[0], 1E-15);
            Assert.AreEqual(12.0 / Math.Sqrt(12.0), moments[1], 1E-15);
            Assert.AreEqual(0.0, moments[2]);
            Assert.AreEqual(1.8, moments[3], 1E-15);
            var parms = U.ParametersFromMoments(moments);
            Assert.AreEqual(0.0, parms[0], 1E-14);
            Assert.AreEqual(12.0, parms[1], 1E-14);

            var constraints = U.GetParameterConstraints(data);
            Assert.AreEqual(2.3, constraints.Item1[0], 0d);
            Assert.AreEqual(5.1, constraints.Item1[1], 0d);
            Assert.AreEqual(2.3, constraints.Item3[0], 0d);
            Assert.AreEqual(5.1, constraints.Item2[1], 0d);
            Assert.IsTrue(constraints.Item2[0] < 2.3 && constraints.Item3[1] > 5.1);

            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.Estimate(new[] { 1.0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.Estimate(new[] { 1.0, double.NaN }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.Estimate(new[] { 1.0, 1.0 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentNullException>(() => U.Estimate(null!, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<NotImplementedException>(() => U.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments));
        }

        /// <summary>
        /// The parametric bootstrap re-estimates a valid distribution from a seeded sample.
        /// </summary>
        [TestMethod]
        public void Test_Bootstrap()
        {
            var U = new Uniform(1, 4);
            var b = (Uniform)U.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 200, 12345);
            Assert.IsTrue(b.ParametersValid);
            Assert.IsTrue(b.Min >= 1 && b.Min < 1.1);
            Assert.IsTrue(b.Max <= 4 && b.Max > 3.9);
        }

        /// <summary>
        /// The quantile function validates the parameters before returning the endpoints.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_ValidatesParametersFirst()
        {
            var U = new Uniform(0, 1);
            U.Min = 2;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.InverseCDF(0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.InverseCDF(1));
        }
    }
}
