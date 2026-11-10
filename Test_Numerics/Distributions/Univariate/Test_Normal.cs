using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Normal distribution algorithm.
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
    public class Test_Normal
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
        private double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

        /// <summary>
        /// Verification of Normal Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.1 page 88.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_MOM_Fit()
        {
            var z = Normal.StandardZ(1 - 1E-16);

            var norm = new Normal();
            norm.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 12665d;
            double true_u2 = 4710d;
            Assert.IsLessThan( 0.01d, (u1 - true_u1) / true_u1 );
            Assert.IsLessThan( 0.01d, (u2 - true_u2) / true_u2 );

        }

        /// <summary>
        /// Verification of Normal Distribution fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_Normal_LMOM_Fit()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var norm = new Normal();
            norm.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 9.957516d;
            double true_u2 = 3.513431d;
            Assert.AreEqual(true_u1, u1, 0.0001d);
            Assert.AreEqual(true_u2, u2, 0.0001d);
            var lmom = norm.LinearMomentsFromParameters(norm.GetParameters);
            Assert.AreEqual(9.9575163d, lmom[0],  0.0001d);
            Assert.AreEqual(1.9822411d, lmom[1], 0.0001d);
            Assert.AreEqual(0.0000000d, lmom[2], 0.0001d);
            Assert.AreEqual(0.1226017d, lmom[3], 0.0001d);
        }

        /// <summary>
        /// Verification of Normal Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// Sigma should be the population standard deviation.
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_MLE_Fit()
        {
            var norm = new Normal();
            norm.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double u1 = norm.Mu;
            double u2 = norm.Sigma;
            double true_u1 = 12665d;
            double true_u2 = 4660d;
            Assert.IsLessThan(0.01d, (u1 - true_u1) / true_u1);
            Assert.IsLessThan(0.01d, (u2 - true_u2) / true_u2);

            var lh = norm.LogLikelihood(sample);
        }

        /// <summary>
        /// Test the quantile function for the Normal Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.2 page 91.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_Quantile()
        {
            var N = new Normal(12665d, 4710d);
            double q100 = N.InverseCDF(0.99d);
            double true_q100 = 23624d;
            Assert.IsLessThan(0.01d, (q100 - true_q100) / true_q100);
            double p = N.CDF(q100);
            double true_p = 0.99d;
            Assert.IsLessThan(0.01d, (p - true_p) / true_p);
        }

        /// <summary>
        /// Test the standard error for the Normal Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 5.1.3 page 94.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_Normal_StandardError()
        {
            var z = Normal.StandardZ(1-double.Epsilon);

            // Maximum Likelihood
            var N = new Normal(12665d, 4710d);
            double qVar99 = Math.Sqrt(N.QuantileVariance(0.99d, 48, ParameterEstimationMethod.MaximumLikelihood));
            double true_qVar99 = 1309d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99);
        }

        /// <summary>
        /// Verifying that input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var N = new Normal();
            Assert.AreEqual(0,N.Mu);
            Assert.AreEqual(1,N.Sigma);

            var N2 = new Normal(1, 1);
            Assert.AreEqual(1,N2.Mu);
            Assert.AreEqual(1,N2.Sigma);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var N = new Normal(double.NaN, 1);
            Assert.IsFalse(N.ParametersValid);

            var N2 = new Normal(1,double.NaN);
            Assert.IsFalse(N2.ParametersValid);

            var N3 = new Normal(double.PositiveInfinity, 1);
            Assert.IsFalse(N3.ParametersValid);

            var N4 = new Normal(1,double.PositiveInfinity);
            Assert.IsFalse(N4.ParametersValid);
        }

        /// <summary>
        /// Testing parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var N = new Normal();
            Assert.AreEqual("Mean (µ)",N.ParametersToString[0, 0]);
            Assert.AreEqual("Std Dev (σ)",N.ParametersToString[1, 0] );
            Assert.AreEqual("0",N.ParametersToString[0, 1]);
            Assert.AreEqual("1",N.ParametersToString[1,1] );
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new Normal(5, 9);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var N = new Normal();
            Assert.AreEqual(0,N.Mean);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(5, N2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var N = new Normal();
            Assert.AreEqual(0, N.Median);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(5, N2.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var N = new Normal();
            Assert.AreEqual(0, N.Mode);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(5, N2.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var N = new Normal();
            Assert.AreEqual(1, N.StandardDeviation);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(9, N2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var N = new Normal();
            Assert.AreEqual(0, N.Skewness);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(0, N2.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var N = new Normal();
            Assert.AreEqual(3, N.Kurtosis);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(3, N2.Kurtosis);
        }

        /// <summary>
        /// Testing minimum and maximum functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var N = new Normal();
            Assert.AreEqual(double.NegativeInfinity,N.Minimum);
            Assert.AreEqual(double.PositiveInfinity, N.Maximum);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(double.NegativeInfinity, N2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, N2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var N = new Normal();
            Assert.AreEqual(0.39894, N.PDF(0),  1e-04);
            Assert.AreEqual(0.24197, N.PDF(1), 1e-04);

            var N2 = new Normal(5, 9);
            Assert.AreEqual(0.03549, N2.PDF(-1), 1e-04);
        }

        /// <summary>
        /// Verifies standard normal quantiles retain their bitwise values across approximation branches and tails.
        /// </summary>
        /// <param name="probability">The nonexceedance probability.</param>
        /// <param name="expectedBits">The IEEE 754 bits of the existing quantile result.</param>
        /// <remarks>Compatibility values were captured from the pre-refactor c0d67b9 implementation.</remarks>
        [TestMethod]
        [DataRow(0d, -4592761413495173886L)]
        [DataRow(double.Epsilon, -4592761413495173886L)]
        [DataRow(1e-300, -4592961304261246354L)]
        [DataRow(1e-20, -4601968184877893045L)]
        [DataRow(1e-12, -4603765893743795414L)]
        [DataRow(0.075d, -4614210144286022365L)]
        [DataRow(0.5d, 0L)]
        [DataRow(0.925d, 4609161892568753444L)]
        [DataRow(0.999999999999d, 4619606146584096420L)]
        [DataRow(0.9999999999999999d, 4620811176048912977L)]
        [DataRow(1d, 9218868437227405312L)]
        public void Test_StandardZ_PreservesBranchAndTailValues(double probability, long expectedBits)
        {
            Assert.AreEqual(expectedBits, BitConverter.DoubleToInt64Bits(Normal.StandardZ(probability)));
        }

        /// <summary>
        /// Verifies distribution quantile endpoints and the existing NaN and invalid-probability contracts.
        /// </summary>
        [TestMethod]
        public void Test_Quantiles_ProbabilityEdges()
        {
            var normal = new Normal();
            Assert.AreEqual(double.NegativeInfinity, normal.InverseCDF(0d));
            Assert.AreEqual(double.PositiveInfinity, normal.InverseCDF(1d));
            Assert.IsTrue(double.IsNaN(normal.InverseCDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(Normal.StandardZ(double.NaN)));

            foreach (double probability in new[] { -1d, 2d, double.NegativeInfinity, double.PositiveInfinity })
            {
                var distributionError = Assert.Throws<ArgumentOutOfRangeException>(() => normal.InverseCDF(probability));
                var standardError = Assert.Throws<ArgumentOutOfRangeException>(() => Normal.StandardZ(probability));
                Assert.AreEqual("probability", distributionError.ParamName);
                Assert.AreEqual("probability", standardError.ParamName);
            }
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var N = new Normal(5,2);
            Assert.AreEqual(0.006209, N.CDF(0),  1e-04);
            Assert.AreEqual(0.30853, N.CDF(4),  1e-04);
            Assert.AreEqual(0.5, N.CDF(5));
        }

        /// <summary>
        /// Verifies the closed-form maximum likelihood estimates: the sample mean and the population standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_ClosedForm()
        {
            var norm = new Normal();
            // scipy.stats.norm.fit([1, 2, 3, 4, 5.5]) = (3.1, 1.5620499351813308)
            double[] parameters = norm.MLE(new[] { 1d, 2d, 3d, 4d, 5.5d });
            Assert.AreEqual(3.1d, parameters[0], 1E-15);
            Assert.AreEqual(1.5620499351813308d, parameters[1], 1E-15);

            norm.Estimate(new[] { 1d, 2d, 3d, 4d, 5.5d }, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(3.1d, norm.Mu, 1E-15);
            Assert.AreEqual(1.5620499351813308d, norm.Sigma, 1E-15);

            // scipy.stats.norm.fit(sample) = (12665.208333333334, 4660.423795814842)
            norm.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(12665.208333333334d, norm.Mu, 1E-12 * 12665.208333333334d);
            Assert.AreEqual(4660.423795814842d, norm.Sigma, 1E-12 * 4660.423795814842d);

            // A large common offset does not cancel the dispersion.
            // scipy.stats.norm.fit([1e9 + 1, 1e9 + 2, 1e9 + 3, 1e9 + 4, 1e9 + 5.5]) = (1000000003.1, 1.562049935181331)
            parameters = norm.MLE(new[] { 1E9 + 1d, 1E9 + 2d, 1E9 + 3d, 1E9 + 4d, 1E9 + 5.5d });
            Assert.AreEqual(1000000003.1d, parameters[0], 1E-15 * 1E9);
            Assert.AreEqual(1.562049935181331d, parameters[1], 1E-12);
        }

        /// <summary>
        /// Verifies that maximum likelihood estimation keeps the sample validation of the parameter constraints.
        /// </summary>
        [TestMethod]
        public void Test_MLE_InvalidSamples()
        {
            var norm = new Normal();
            Assert.Throws<ArgumentNullException>(() => norm.MLE(null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => norm.MLE(new[] { 1d, 2d, 3d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => norm.MLE(new[] { 2d, 2d, 2d, 2d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => norm.MLE(new[] { 1d, 2d, double.NaN, 4d }));
            Assert.Throws<ArgumentOutOfRangeException>(() => norm.MLE(new[] { 1d, 2d, double.PositiveInfinity, 4d }));
        }

        /// <summary>
        /// Verifies conditional moments with infinite limits against truncated Normal moments.
        /// </summary>
        /// <remarks>
        /// The moments are the conditional mean and the conditional moments about µ, which for a standard Normal equal
        /// scipy.stats.truncnorm.moment(k, a, b) for k = 1..4.
        /// </remarks>
        [TestMethod]
        public void Test_ConditionalMoments_InfiniteLimits()
        {
            var standard = new Normal(0d, 1d);
            // scipy.stats.truncnorm.moment(k, -inf, 0), k = 1..4
            AssertMoments(new[] { -0.7978845608028654d, 1d, -1.5957691216057308d, 3d }, standard.ConditionalMoments(double.NegativeInfinity, 0d), 1E-15);
            // scipy.stats.truncnorm.moment(k, 0, inf), k = 1..4
            AssertMoments(new[] { 0.7978845608028654d, 1d, 1.5957691216057308d, 3d }, standard.ConditionalMoments(0d, double.PositiveInfinity), 1E-15);
            // scipy.stats.truncnorm.moment(k, -inf, inf), k = 1..4
            AssertMoments(new[] { 0d, 1d, 0d, 3d }, standard.ConditionalMoments(double.NegativeInfinity, double.PositiveInfinity), 1E-15);
            // mpmath.quad(t^k exp(-t^2/2), [5, inf]) / mpmath.quad(exp(-t^2/2), [5, inf]) at 50 digits, k = 1..4
            AssertMoments(new[] { 5.1865039671258421156d, 26.932519835629210578d, 140.03560711239773712d, 729.11055539761789619d }, standard.ConditionalMoments(5d, double.PositiveInfinity), 1E-13);
            // Location and scale: mean 10 + 2 m1, moments about 10 scaled by 2^k, from scipy.stats.truncnorm.moment(k, -inf, 0.5)
            // = [-0.5091604338370335, 0.7454197830814833, -1.1456109761333253, 2.1726142950148204]
            var normal = new Normal(10d, 2d);
            AssertMoments(new[] { 10d + 2d * -0.5091604338370335d, 4d * 0.7454197830814833d, 8d * -1.1456109761333253d, 16d * 2.1726142950148204d },
                normal.ConditionalMoments(double.NegativeInfinity, 11d), 1E-14);
            AssertMoments(new[] { 10d + 2d * 0.7978845608028654d, 4d, 8d * 1.5957691216057308d, 48d }, normal.ConditionalMoments(10d, double.PositiveInfinity), 1E-15);
        }

        /// <summary>
        /// Verifies conditional moments on intervals in the upper tail, where the interval mass is differenced as survival probabilities.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_UpperTail()
        {
            var standard = new Normal(0d, 1d);
            // mpmath.quad(t^k exp(-t^2/2), [a, b]) / mpmath.quad(exp(-t^2/2), [a, b]) at 50 digits, k = 1..4
            AssertMoments(new[] { 8.1211889929797971226d, 65.967859202478959388d, 535.97037693355651102d, 4355.5936971380993874d }, standard.ConditionalMoments(8d, 9d), 1E-13);
            AssertMoments(new[] { 3.2604542855900211539d, 10.679839942199728351d, 35.154336739497886921d, 116.31543769842735713d }, standard.ConditionalMoments(3d, 4d), 1E-13);
            AssertMoments(new[] { 0.96355397941640390687d, 1.2408884948541009767d, 1.9873300825463330579d, 3.7377210154906842412d }, standard.ConditionalMoments(0.25d, 40d), 1E-13);
            AssertMoments(new[] { 30.033259667433607515d, 901.99779002300608895d, 27090.000220024984343d, 813604.0043907705208d }, standard.ConditionalMoments(30d, 31d), 1E-12);
            // Location and scale.
            var normal = new Normal(-50d, 4d);
            AssertMoments(new[] { -50d + 4d * 8.1211889929797971226d, 16d * 65.967859202478959388d, 64d * 535.97037693355651102d, 256d * 4355.5936971380993874d },
                normal.ConditionalMoments(-18d, -14d), 1E-13);
        }

        /// <summary>
        /// Verifies conditional moments on intervals in the lower tail, where the interval mass is differenced as lower-tail probabilities.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_LowerTail()
        {
            var standard = new Normal(0d, 1d);
            // Mirror images of the upper-tail mpmath values: odd moments change sign.
            AssertMoments(new[] { -8.1211889929797971226d, 65.967859202478959388d, -535.97037693355651102d, 4355.5936971380993874d }, standard.ConditionalMoments(-9d, -8d), 1E-13);
            AssertMoments(new[] { -3.2604542855900211539d, 10.679839942199728351d, -35.154336739497886921d, 116.31543769842735713d }, standard.ConditionalMoments(-4d, -3d), 1E-13);
            AssertMoments(new[] { -30.033259667433607515d, 901.99779002300608895d, -27090.000220024984343d, 813604.0043907705208d }, standard.ConditionalMoments(-31d, -30d), 1E-12);
            var normal = new Normal(-50d, 4d);
            AssertMoments(new[] { -50d - 4d * 8.1211889929797971226d, 16d * 65.967859202478959388d, -64d * 535.97037693355651102d, 256d * 4355.5936971380993874d },
                normal.ConditionalMoments(-86d, -82d), 1E-13);
        }

        /// <summary>
        /// Verifies the not-a-number contract for empty or reversed intervals.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_InvalidInterval()
        {
            var normal = new Normal(0d, 1d);
            foreach (var moments in new[] { normal.ConditionalMoments(1d, 1d), normal.ConditionalMoments(2d, -2d) })
                foreach (double m in moments) Assert.IsTrue(double.IsNaN(m));
        }

        /// <summary>
        /// Verifies that finite intervals containing the mean keep the closed-form erf-based evaluation bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_FiniteLowerOrCentralIntervals_SameBits()
        {
            foreach (double mu in new[] { 0d, 3.7d, -120.25d })
                foreach (double sigma in new[] { 1d, 0.35d, 42d })
                    foreach (double alpha in new[] { -6d, -2.5d, -1d, -0.3d, 0d })
                        foreach (double width in new[] { 0.01d, 0.5d, 2d, 7d })
                        {
                            double a = mu + alpha * sigma, b = a + width * sigma;
                            if ((a - mu) / sigma > 0d || (b - mu) / sigma < 0d) continue;
                            var expected = ErfConditionalMoments(mu, sigma, a, b);
                            var actual = new Normal(mu, sigma).ConditionalMoments(a, b);
                            for (int k = 0; k < 4; k++) Assert.AreEqual(expected[k], actual[k], 0d);
                        }
        }

        /// <summary>The closed-form truncated Normal moments with erf-based probabilities, for finite limits.</summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation.</param>
        /// <param name="a">The finite lower limit.</param>
        /// <param name="b">The finite upper limit.</param>
        /// <returns>The conditional mean and the conditional moments about the mean.</returns>
        private static double[] ErfConditionalMoments(double mu, double sigma, double a, double b)
        {
            double α = (a - mu) / sigma;
            double β = (b - mu) / sigma;
            double Φα = 0.5 * (1 + Erf.Function(α / Math.Sqrt(2)));
            double Φβ = 0.5 * (1 + Erf.Function(β / Math.Sqrt(2)));
            double Z = Φβ - Φα;
            double φα = Math.Exp(-0.5 * α * α) / Math.Sqrt(2 * Math.PI);
            double φβ = Math.Exp(-0.5 * β * β) / Math.Sqrt(2 * Math.PI);
            double λ = (φα - φβ) / Z;
            double δ = (α * φα - β * φβ) / Z;
            double τ = ((α * α + 2) * φα - (β * β + 2) * φβ) / Z;
            double κ = (α * α * α * φα - β * β * β * φβ) / Z;
            return new[] { mu + sigma * λ, sigma * sigma * (1 + δ), sigma * sigma * sigma * τ, sigma * sigma * sigma * sigma * (3 + κ + 3 * δ) };
        }

        /// <summary>
        /// StandardZ keeps a finite lower end at p = 0 and returns positive infinity at p = 1, as its remarks state.
        /// </summary>
        /// <remarks>Reference: mpmath 1.3.0 at 40 digits, the roots of ncdf(z) = 2^-1074 (-38.467405617144346) and of ncdf(-z) = 2^-53 (8.209536151601387).</remarks>
        [TestMethod]
        public void Test_StandardZ_Ends()
        {
            Assert.AreEqual(-38.467405617144346, Normal.StandardZ(0d), 1E-12);
            Assert.AreEqual(Normal.StandardZ(double.Epsilon), Normal.StandardZ(0d));
            Assert.IsTrue(double.IsPositiveInfinity(Normal.StandardZ(1d)));
            Assert.AreEqual(8.209536151601387, Normal.StandardZ(1d - Math.Pow(2d, -53)), 1E-12);
        }

        /// <summary>Asserts four moments to a relative tolerance, with the same tolerance as an absolute floor.</summary>
        /// <param name="expected">The reference moments.</param>
        /// <param name="actual">The computed moments.</param>
        /// <param name="tolerance">The relative tolerance.</param>
        private static void AssertMoments(double[] expected, double[] actual, double tolerance)
        {
            for (int k = 0; k < 4; k++)
                Assert.AreEqual(expected[k], actual[k], tolerance * Math.Max(1d, Math.Abs(expected[k])));
        }
    }
}
