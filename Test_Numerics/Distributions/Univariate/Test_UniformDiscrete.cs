using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Discrete Uniform distribution algorithm.
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
    public class Test_UniformDiscrete
    {

        /// <summary>
        /// Verified using Palisade's @Risk and Accord.Net
        /// </summary>
        [TestMethod()]
        public void Test_UniformDiscreteDist()
        {
            double true_mean = 4.0d;
            double true_median = 4.0d;
            // scipy.stats.randint(2, 7).std() = 1.4142135624; N=5, sqrt((N^2-1)/12)
            double true_stdDev = Math.Sqrt((5.0d * 5.0d - 1.0d) / 12.0d);
            int true_skew = 0;
            double true_kurt = 1.7d;
            double true_pdf = 0.2d;
            double true_cdf = 0.2d;
            double true_icdf05 = 2.0d;
            double true_icdf95 = 6.0d;
            var U = new UniformDiscrete(2d, 6d);
            Assert.AreEqual(true_mean, U.Mean, 0.0001d);
            Assert.AreEqual(true_median, U.Median, 0.0001d);
            Assert.AreEqual(true_stdDev, U.StandardDeviation, 0.0001d);
            Assert.AreEqual(true_skew, U.Skewness, 0.0001d);
            Assert.AreEqual(true_kurt, U.Kurtosis, 0.0001d);
            Assert.AreEqual(true_pdf, U.PDF(4.0d), 0.0001d);
            Assert.AreEqual(true_cdf, U.CDF(2.0d), 0.0001d);
            Assert.AreEqual(true_icdf05, U.InverseCDF(0.17d), 0.0001d);
            Assert.AreEqual(true_icdf95, U.InverseCDF(0.87d), 0.0001d);
        }
        
        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0, U.Min);
            Assert.AreEqual(1, U.Max);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(2, U2.Min);
            Assert.AreEqual(10, U2.Max);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var U = new UniformDiscrete(1, 0);
            Assert.IsFalse(U.ParametersValid);

            var U2 = new UniformDiscrete(double.NaN, 0);
            Assert.IsFalse(U2.ParametersValid);

            var U3 = new UniformDiscrete(0, double.NaN);
            Assert.IsFalse(U3.ParametersValid);

            var U4 = new UniformDiscrete(0, double.PositiveInfinity);
            Assert.IsFalse(U4.ParametersValid);

            var U5 = new UniformDiscrete(double.PositiveInfinity, 0);
            Assert.IsFalse(U5.ParametersValid);
        }
        
        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual("Min", U.ParametersToString[0, 0]);
            Assert.AreEqual("Max", U.ParametersToString[1, 0]);
            Assert.AreEqual("0", U.ParametersToString[0, 1]);
            Assert.AreEqual("1", U.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0.5, U.Mean);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(6, U2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0.5, U.Median);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(6, U2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod]
        public void Test_Mode()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(double.NaN, U.Mode);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(double.NaN, U2.Mode);
        }

        /// <summary>
        /// Testing Standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            // scipy.stats.randint(0, 2).std() = 0.5; N=2, sqrt((4-1)/12) = 0.5
            var U = new UniformDiscrete();
            Assert.AreEqual(0.5, U.StandardDeviation, 1e-05);

            // scipy.stats.randint(2, 11).std() = 2.5819888975; N=9, sqrt((81-1)/12)
            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(2.5819888975, U2.StandardDeviation, 1e-04);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod]
        public void Test_Skewness()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0, U.Skewness);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(0, U2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod]
        public void Test_Kurtosis()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(1, U.Kurtosis);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(1.77, U2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0, U.Minimum);
            Assert.AreEqual(1, U.Maximum);

            var U2 = new UniformDiscrete(2, 10);
            Assert.AreEqual(2, U2.Minimum);
            Assert.AreEqual(10, U2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0, U.PDF(-1));
            Assert.AreEqual(0, U.PDF(2));
            Assert.AreEqual(0.5, U.PDF(1));
        }

        /// <summary>
        /// Testing CDF.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0.5, U.CDF(0));
            Assert.AreEqual(1, U.CDF(1));
            // The CDF is a step function: scipy.stats.randint(0, 2).cdf(0.5) = 0.5
            Assert.AreEqual(0.5, U.CDF(0.5));
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var U = new UniformDiscrete();
            Assert.AreEqual(0, U.InverseCDF(0));
            Assert.AreEqual(1, U.InverseCDF(1));
            Assert.AreEqual(0, U.InverseCDF(0.3));
        }

        /// <summary>
        /// Verify discrete uniform standard deviation uses the correct formula: sqrt((N^2-1)/12).
        /// Reference: scipy.stats.randint(low, high+1).std()
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation_Discrete()
        {
            // randint(0, 2).std() = 0.5000000000  {0,1}
            var u1 = new UniformDiscrete(0, 1);
            Assert.AreEqual(0.5, u1.StandardDeviation, 1E-6);

            // randint(0, 6).std() = 1.7078251277  {0,...,5}
            var u2 = new UniformDiscrete(0, 5);
            Assert.AreEqual(1.7078251277, u2.StandardDeviation, 1E-6);

            // randint(0, 11).std() = 3.1622776602  {0,...,10}
            var u3 = new UniformDiscrete(0, 10);
            Assert.AreEqual(3.1622776602, u3.StandardDeviation, 1E-6);

            // randint(1, 7).std() = 1.7078251277  {1,...,6}
            var u4 = new UniformDiscrete(1, 6);
            Assert.AreEqual(1.7078251277, u4.StandardDeviation, 1E-6);

            // randint(3, 8).std() = 1.4142135624  {3,...,7}
            var u5 = new UniformDiscrete(3, 7);
            Assert.AreEqual(1.4142135624, u5.StandardDeviation, 1E-6);
        }

        /// <summary>
        /// The PMF is zero away from the integers and the CDF floors its argument, as in scipy.stats.randint.
        /// </summary>
        [TestMethod()]
        public void Test_NonIntegerArguments()
        {
            var U = new UniformDiscrete(0, 1);
            // scipy.stats.randint(0, 2).pmf(0.5) = 0, .cdf(0.5) = 0.5
            Assert.AreEqual(0d, U.PDF(0.5));
            Assert.AreEqual(0.5, U.CDF(0.5));
            var D = new UniformDiscrete(1, 6);
            // scipy.stats.randint(1, 7).cdf(3.5) = 0.5
            Assert.AreEqual(0.5, D.CDF(3.5));
            Assert.AreEqual(1d / 6d, D.PDF(3));
            Assert.AreEqual(1d, D.CDF(double.PositiveInfinity));
            Assert.AreEqual(0d, D.CDF(double.NegativeInfinity));
            Assert.IsTrue(double.IsNaN(U.PDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(U.CDF(double.NaN)));
        }

        /// <summary>
        /// The quantile is the smallest k with CDF(k) ≥ p, including where p·N is an integer.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_SmallestK()
        {
            // Adjacent doubles for nonnegative arguments.
            static double NextUp(double x) => BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(x) + 1);
            static double NextDown(double x) => BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(x) - 1);

            // scipy.stats.randint(0, 2).ppf(0.5) = 0
            Assert.AreEqual(0d, new UniformDiscrete(0, 1).InverseCDF(0.5));
            var D = new UniformDiscrete(1, 6);
            // scipy.stats.randint(1, 7).ppf(0.5) = 3, .ppf(1/6) = 1
            Assert.AreEqual(3d, D.InverseCDF(0.5));
            Assert.AreEqual(1d, D.InverseCDF(1d / 6d));
            Assert.AreEqual(1d, D.InverseCDF(NextDown(1d / 6d)));
            // CDF(1) = 1/6 lies one unit in the last place below p, so the smallest k with CDF(k) >= p is 2
            // (scipy.stats.randint(1, 7).ppf returns 1 here because p*6 rounds down to exactly 1).
            Assert.AreEqual(2d, D.InverseCDF(NextUp(1d / 6d)));

            foreach (var (min, max) in new (double, double)[] { (0, 1), (1, 6), (-3, 4), (10, 109), (5, 5) })
            {
                var U = new UniformDiscrete(min, max);
                double n = max - min + 1;
                for (int j = 0; j <= n; j++)
                {
                    foreach (double p in new[] { j / n, NextDown(j / n), NextUp(j / n) })
                    {
                        if (!(p > 0 && p < 1)) continue;
                        double k = U.InverseCDF(p);
                        Assert.IsTrue(k >= min && k <= max);
                        Assert.IsGreaterThanOrEqualTo(p, U.CDF(k));
                        Assert.IsTrue(k == min || U.CDF(k - 1) < p);
                    }
                }
            }
        }

        /// <summary>
        /// Non-integer bounds are invalid; Min = Max is a valid point mass whose kurtosis is undefined.
        /// </summary>
        [TestMethod()]
        public void Test_IntegerParametersAndPointMass()
        {
            Assert.IsFalse(new UniformDiscrete(0.5, 2.2).ParametersValid);
            Assert.IsFalse(new UniformDiscrete(0, 2.2).ParametersValid);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new UniformDiscrete(0.5, 2).PDF(1));

            var P = new UniformDiscrete(3, 3);
            Assert.IsTrue(P.ParametersValid);
            Assert.AreEqual(1d, P.PDF(3));
            Assert.AreEqual(0d, P.PDF(2));
            Assert.AreEqual(0d, P.CDF(2));
            Assert.AreEqual(1d, P.CDF(3));
            Assert.AreEqual(3d, P.InverseCDF(0.4));
            Assert.AreEqual(3d, P.Mean);
            Assert.IsTrue(double.IsNaN(P.Kurtosis));
        }

        /// <summary>
        /// Maximum likelihood gives the sample minimum and maximum.
        /// </summary>
        [TestMethod()]
        public void Test_MLE()
        {
            var data = new double[] { 4, 2, 9, 7, 3, 3 };
            var U = new UniformDiscrete();
            U.Estimate(data, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(2d, U.Min);
            Assert.AreEqual(9d, U.Max);

            var constraints = U.GetParameterConstraints(data);
            Assert.AreEqual(2d, constraints.Item1[0]);
            Assert.AreEqual(9d, constraints.Item1[1]);
            Assert.AreEqual(2d, constraints.Item3[0]);
            Assert.AreEqual(9d, constraints.Item2[1]);
            Assert.IsTrue(constraints.Item2[0] < 2d && constraints.Item3[1] > 9d);

            // A constant sample estimates a point mass.
            U.Estimate(new double[] { 5, 5, 5 }, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(5d, U.Min);
            Assert.AreEqual(5d, U.Max);
            Assert.IsTrue(U.ParametersValid);

            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.Estimate(new[] { 1.0, 2.5 }, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.Estimate(new double[0], ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<ArgumentNullException>(() => U.Estimate(null!, ParameterEstimationMethod.MaximumLikelihood));
            Assert.ThrowsExactly<NotImplementedException>(() => U.Estimate(data, ParameterEstimationMethod.MethodOfMoments));
        }

        /// <summary>
        /// The parametric bootstrap re-estimates a valid distribution from a seeded sample.
        /// </summary>
        [TestMethod()]
        public void Test_Bootstrap()
        {
            var U = new UniformDiscrete(1, 6);
            var b = (UniformDiscrete)U.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 200, 12345);
            Assert.IsTrue(b.ParametersValid);
            Assert.AreEqual(1d, b.Min);
            Assert.AreEqual(6d, b.Max);
        }

        /// <summary>
        /// The quantile function validates the parameters before returning the endpoints.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF_ValidatesParametersFirst()
        {
            var U = new UniformDiscrete(0, 1);
            U.Min = 2;
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.InverseCDF(0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => U.InverseCDF(1));
        }
    }

}
