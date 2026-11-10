using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Gamma distribution algorithm.
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
    public class Test_GammaDistribution
    {
        // Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
        private double[] sample = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

        /// <summary>
        /// Method of Moments.
        /// </summary>
        /// <remarks>
        /// Solution found in the book “The Gamma family and derived distributions applied in hydrology” by Bobée and Ashkar, 1991.
        /// </remarks>
        [TestMethod()]
        public void Test_GammaDist_MOM()
        {
            var G = new GammaDistribution();
            G.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double alpha = 1d / G.Theta;
            double lambda = G.Kappa;
            double trueA = 0.08317d;
            double trueL = 15.91188d;
            Assert.IsLessThan(0.01d, (alpha - trueA) / trueA);
            Assert.IsLessThan(0.01d, (lambda - trueL) / trueL);
        }

        /// <summary>
        /// Verifies gamma distribution fitting with linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_GammaDist_LMOM_Fit()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var G = new GammaDistribution();
            G.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double scale = G.Theta;
            double shape = G.Kappa;
            double true_scale = 1.280143d;
            double true_shape = 7.778442d;
            Assert.AreEqual(true_scale, scale, 0.0001d);
            Assert.AreEqual(true_shape, shape, 0.0001d);
            var lmom = G.LinearMomentsFromParameters(G.GetParameters);
            Assert.AreEqual(9.9575163d, lmom[0],  0.0001d);
            Assert.AreEqual(1.9822363d, lmom[1],  0.0001d);
            Assert.AreEqual(0.1175059d, lmom[2],  0.0001d);
            Assert.AreEqual(0.1268391d, lmom[3],  0.0001d);
        }

        /// <summary>
        /// Verification of Gamma fit with maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 7.4 page 93.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GammaDist_MLE()
        {
            var G = new GammaDistribution();
            G.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double alpha = 1d / G.Theta;
            double lambda = G.Kappa;
            double trueA = 0.08833d;
            double trueL = 16.89937d;
            Assert.IsLessThan(0.01d, (alpha - trueA) / trueA);
            Assert.IsLessThan(0.01d, (lambda - trueL) / trueL);
        }

        /// <summary>
        /// Test the quantile function for the Gamma Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 5.3 page 52.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GammaDist_Quantile()
        {
            var G = new GammaDistribution(1d / 0.08833d, 16.89937d);
            double q1000 = G.InverseCDF(0.99d);
            double true_1000 = 315.87d;
            Assert.IsLessThan(0.01d, (q1000 - true_1000) / true_1000);
            double p = G.CDF(q1000);
            double true_p = 0.99d;
            Assert.AreEqual(true_p, p);
        }

        /// <summary>
        /// Test the standard error for the Gamma Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 6.1.3 page 138.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GammaDist_StandardError()
        {

            // Method of Moments
            var G = new GammaDistribution(1d / 0.08317d, 15.9118d);
            double se1000 = Math.Sqrt(G.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MethodOfMoments));
            double true_se1000 = 16.024d;
            Assert.IsLessThan(0.01d, (se1000 - true_se1000) / true_se1000);

            // Maximum Likelihood
            G = new GammaDistribution(1d / 0.08833d, 16.89937d);
            se1000 = Math.Sqrt(G.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MaximumLikelihood));
            true_se1000 = 15.022d;
            Assert.IsLessThan(0.01d, (se1000 - true_se1000) / true_se1000);
        }

        /// <summary>
        /// Testing parameters are creating gamma distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var G = new GammaDistribution(2, 10);
            Assert.AreEqual(2,G.Theta);
            Assert.AreEqual(10, G.Kappa);

            var G2 = new GammaDistribution(-1, 4);
            Assert.AreEqual(-1,G2.Theta);
            Assert.AreEqual(4,G2.Kappa);

        }

        /// <summary>
        /// Validating inverse scale parameter.
        /// </summary>
        [TestMethod()]
        public void Test_Rate()
        {
            var G = new GammaDistribution(2, 2);
            Assert.AreEqual(0.5,G.Rate);

            var G2 = new GammaDistribution();
            Assert.AreEqual(0.1,G2.Rate);
        }

        /// <summary>
        /// Checking parameters to string function.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var G = new GammaDistribution();
            Assert.AreEqual("Scale (θ)", G.ParametersToString[0, 0]);
            Assert.AreEqual("Shape (κ)",G.ParametersToString[1, 0]);
            Assert.AreEqual("10",G.ParametersToString[0, 1]);
            Assert.AreEqual("2",G.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Checking bad parameters for Gamma distribution.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var G = new GammaDistribution(-1, 0);
            Assert.IsFalse(G.ParametersValid);

            var G2 = new GammaDistribution(0, 0);
            Assert.IsFalse(G2.ParametersValid);

            var G3 = new GammaDistribution(0,1);
            Assert.IsFalse(G3.ParametersValid);

            var G4 = new GammaDistribution(1, 0);
            Assert.IsFalse(G4.ParametersValid);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GammaDistribution();
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Checking mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(20,G.Mean);
        }

        /// <summary>
        /// Checking median function output. (There is no closed form solution)
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(G.InverseCDF(0.5), G.Median);
        }

        /// <summary>
        /// Checking Mode function with different scales. For shape at or below 1 the mode is 0, the lower bound of the support.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(10, G.Mode);

            var G2 = new GammaDistribution(10, 1);
            Assert.AreEqual(0d, G2.Mode);

            var G3 = new GammaDistribution(10, 0.5);
            Assert.AreEqual(0d, G3.Mode);
        }

        /// <summary>
        /// Checking standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(14.142135, G.StandardDeviation,1e-04);

            var G2 = new GammaDistribution(1, 2);
            Assert.AreEqual(1.4142135, G2.StandardDeviation,  1e-04);
        }

        /// <summary>
        /// Checking skew function with different shapes.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(1.4142135, G.Skewness, 1e-04);

            var G2 = new GammaDistribution(10, 100);
            Assert.AreEqual(0.2, G2.Skewness);
        }

        /// <summary>
        /// Checking Kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(6, G.Kurtosis);

            var G2 = new GammaDistribution(10, 6);
            Assert.AreEqual(4, G2.Kurtosis);

            var G3 = new GammaDistribution(10, 2.5);
            Assert.AreEqual(5.4,G3.Kurtosis);
        }

        /// <summary>
        /// Checking minimum is 0.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(0,G.Minimum);
        }

        /// <summary>
        /// Checking maximum is infinity
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var G = new GammaDistribution();
            Assert.AreEqual(double.PositiveInfinity,G.Maximum);
        }

        /// <summary>
        /// Verification of Gamma fit with maximum likelihood using Newton-Raphson.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 7.4 page 93.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void ValidateMLE_NR()
        {
            var G = new GammaDistribution();
            G.MLE_NR(sample);
            double alpha = 1d / G.Theta;
            double lambda = G.Kappa;
            double trueA = 0.08833d;
            double trueL = 16.89937d;
            Assert.IsLessThan( 0.2d, (alpha - trueA) / trueA );
            Assert.IsLessThan(0.01d, (lambda - trueL) / trueL);
        }

        /// <summary>
        /// The Newton-Raphson fit sets the scale to the sample mean over the fitted shape and agrees with the
        /// general maximum likelihood fit.
        /// </summary>
        /// <remarks>
        /// The gamma likelihood equations give scale = mean / shape at the maximum, so the Newton-Raphson solution
        /// and the Nelder-Mead maximum likelihood solution must agree.
        /// </remarks>
        [TestMethod()]
        public void Test_MLE_NR_ScaleIsMeanOverFittedShape()
        {
            var sample = new GammaDistribution(2d, 4d).GenerateRandomValues(300, 12345);
            double mean = 0d;
            for (int i = 0; i < sample.Length; i++) mean += sample[i];
            mean /= sample.Length;

            var nr = new GammaDistribution(2d, 4d);
            nr.MLE_NR(sample);
            Assert.AreEqual(mean / nr.Kappa, nr.Theta, 1E-12 * nr.Theta);

            // The fitted shape solves the profile likelihood equation ln(k) - digamma(k) = ln(mean) - mean(ln x).
            double lnsum = 0d;
            for (int i = 0; i < sample.Length; i++) lnsum += Math.Log(sample[i]);
            double s = Math.Log(mean) - lnsum / sample.Length;
            Assert.AreEqual(s, Math.Log(nr.Kappa) - Numerics.Mathematics.SpecialFunctions.Gamma.Digamma(nr.Kappa), 1E-10);

            // The general maximum likelihood fit (Nelder-Mead) reaches the same optimum, within its own
            // convergence tolerance, and never a higher likelihood.
            var mle = new GammaDistribution();
            mle.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(mle.Kappa, nr.Kappa, 1E-3 * mle.Kappa);
            Assert.AreEqual(mle.Theta, nr.Theta, 1E-3 * mle.Theta);
            Assert.IsGreaterThanOrEqualTo(mle.LogLikelihood(sample) - 1E-9, nr.LogLikelihood(sample));
            Assert.AreEqual(mle.LogLikelihood(sample), nr.LogLikelihood(sample), 1E-5);

            // The fit must not depend on the shape held before the call.
            var other = new GammaDistribution(10d, 0.5d);
            other.MLE_NR(sample);
            Assert.AreEqual(nr.Theta, other.Theta, 1E-12 * nr.Theta);
            Assert.AreEqual(nr.Kappa, other.Kappa, 1E-12 * nr.Kappa);
        }

        /// <summary>
        /// Verification of Gamma fit with maximum likelihood using Newton-Raphson.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 7.4 page 93.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void ValidateMLE_Bobee()
        {
            var G = new GammaDistribution();
            G.MLE_Bobee(sample);
            double alpha = 1d / G.Theta;
            double lambda = G.Kappa;
            double trueA = 0.08833d;
            double trueL = 16.89937d;
            Assert.IsLessThan(0.2d, (alpha - trueA) / trueA);
            Assert.IsLessThan(0.01d, (lambda - trueL) / trueL);
        }

        /// <summary>
        /// Checking PDF with different locations and parameters.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var G = new GammaDistribution(10,1);
            Assert.AreEqual(0.090483, G.PDF(1),  1e-04);
            Assert.AreEqual(0.036787, G.PDF(10), 1e-04);

            var G2 = new GammaDistribution(1,1);
            Assert.AreEqual(0.367879, G2.PDF(1), 1e-04);
            Assert.AreEqual(0.0000453999, G2.PDF(10),  1e-10);
        }

        /// <summary>
        /// A Gamma shape below one has a genuine density singularity at zero.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF_ShapeBelowOneIsPositiveInfinityAtZero()
        {
            Assert.AreEqual(double.PositiveInfinity, new GammaDistribution(2.0, 0.5).LogPDF(0.0));
        }

        /// <summary>
        /// Checking CDF function with different parameters at different locations.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var G = new GammaDistribution(10, 1);
            Assert.AreEqual(0.09516258, G.CDF(1),  1e-04);
            Assert.AreEqual(0.63212, G.CDF(10), 1e-04);

            var G2 = new GammaDistribution(1, 1);
            Assert.AreEqual(0.999954, G2.CDF(10), 1e-04);

            var G3 = new GammaDistribution(0.1, 10);
            Assert.AreEqual(0.54207028, G3.CDF(1), 1e-04);
            Assert.AreEqual(0.999999, G3.CDF(10), 1e-04);
        }

        /// <summary>
        /// Checking inverse CDF function with different parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var G = new GammaDistribution(10,1);
            Assert.AreEqual(0,G.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,G.InverseCDF(1) );

            var G2 = new GammaDistribution(9.30149316e-07,1082.2442991605726);
            Assert.AreEqual(1.0792e-03, G2.InverseCDF(0.99),1e-04);
            Assert.AreEqual(1.0817e-03, G2.InverseCDF(0.9919),1e-04);
            Assert.AreEqual(1.0834e-03, G2.InverseCDF(0.993),1e-04);
        }

        /// <summary>
        /// Validating Wilson Hilferty Inverse CDF function.
        /// <see href="extension://efaidnbmnnnibpcajpcglclefindmkaj/https://faculty.washington.edu/smukherj/mypubs/Gamma_normal_methods_2008.pdf"/>
        /// </summary>
        [TestMethod()]
        public void ValidateWilsonHilfertyInverseCDF()
        {
            var G = new GammaDistribution(1, 1);
            Assert.AreEqual(4.62111, G.WilsonHilfertyInverseCDF(.99),1e-04);
            Assert.AreEqual(6.92202, G.WilsonHilfertyInverseCDF(0.999),1e-04);

            var G2 = new GammaDistribution(1, 2);
            Assert.AreEqual(0.3566877, G2.WilsonHilfertyInverseCDF(0.05),1e-04);
            Assert.AreEqual(0.5326, G2.WilsonHilfertyInverseCDF(0.10), 1e-04);
        }

        /// <summary>
        /// Verifies the frequency-factor skew derivative at and around zero skew, where it is the smooth derivative of the
        /// Cornish-Fisher series with no switch.
        /// </summary>
        /// <remarks>
        /// Reference values are the derivative and value of the exact order-6 skew series of the standardized gamma quantile,
        /// derived with SymPy 1.14 and evaluated with mpmath 1.3.0 at 60 digits.
        /// </remarks>
        [TestMethod]
        public void Test_PartialKp_ZeroSkewLimitAndContinuity()
        {
            const double probability = 0.9d;
            double z = Normal.StandardZ(probability);
            double expected = (z * z - 1.0d) / 6.0d;

            Assert.AreEqual(expected, GammaDistribution.PartialKp(0.0d, probability), 1E-15);
            Assert.AreEqual(0.10705296165768363, GammaDistribution.PartialKp(9.9E-5, probability), 1E-15);
            Assert.AreEqual(0.10707184335968709, GammaDistribution.PartialKp(-9.9E-5, probability), 1E-15);
            Assert.AreEqual(0.10705277093275638, GammaDistribution.PartialKp(1.01E-4, probability), 1E-15);
            Assert.AreEqual(0.10707203408328514, GammaDistribution.PartialKp(-1.01E-4, probability), 1E-15);
            Assert.AreEqual(1.281562164255128, GammaDistribution.FrequencyFactorKp(9.9E-5, probability), 1E-15);
            Assert.AreEqual(1.2815623783608607, GammaDistribution.FrequencyFactorKp(1.01E-4, probability), 1E-15);
            Assert.AreEqual(1.2815409658994292, GammaDistribution.FrequencyFactorKp(-9.9E-5, probability), 1E-15);
            Assert.AreEqual(1.2815407517555517, GammaDistribution.FrequencyFactorKp(-1.01E-4, probability), 1E-15);

            const double h = 1E-3;
            double finiteDifference = (GammaDistribution.FrequencyFactorKp(h, probability) -
                                       GammaDistribution.FrequencyFactorKp(-h, probability)) / (2.0d * h);
            Assert.AreEqual(expected, finiteDifference, 1E-4);
        }

        /// <summary>
        /// Near zero skew the frequency factor and its derivative agree with the exact standardized Pearson III quantile and
        /// its skew derivative, with no switch at small skew.
        /// </summary>
        /// <remarks>
        /// Columns: skew, probability, the exact standardized quantile and its skew derivative, computed with mpmath 1.3.0 at
        /// 60 digits by Newton iteration on the quadrature of the gamma density and by implicit differentiation. At these
        /// skews the fifth- and higher-order terms of the polynomial contribute below 1E-15. Returning the Normal variate here
        /// instead would be in error by up to 3.6E-4 at the 1E-6 quantile.
        /// </remarks>
        [TestMethod]
        [DataRow(9.9e-05, 1e-06, -4.753067995664727, 3.5990718404725808)]
        [DataRow(0.000101, 1e-06, -4.753060797523106, 3.599069780999969)]
        [DataRow(-9.9e-05, 0.99, 2.326275077531289, 0.7353208179293219)]
        [DataRow(0.0003, 1e-06, -4.752344603028048, 3.598864832474004)]
        [DataRow(-0.0003, 0.99, 2.3261272770106918, 0.7353311276274129)]
        [DataRow(0.0003, 0.5, -4.9999999933333326e-05, -0.166666666)]
        [DataRow(-0.0003, 0.01, -2.326568466452883, 0.735300340244945)]
        public void Test_FrequencyFactorKp_ExactQuantileNearZeroSkew(double skew, double probability, double quantile, double derivative)
        {
            Assert.AreEqual(quantile, GammaDistribution.FrequencyFactorKp(skew, probability), 2E-15 * (1d + Math.Abs(quantile)));
            Assert.AreEqual(derivative, GammaDistribution.PartialKp(skew, probability), 2E-15 * (1d + Math.Abs(derivative)));
        }

        /// <summary>
        /// At the probability endpoints the frequency factor is the Normal variate at and near zero skew: positive infinity
        /// at probability one and <see cref="Normal.StandardZ(double)"/> of zero at probability zero.
        /// </summary>
        [TestMethod]
        public void Test_FrequencyFactorKp_EndpointsNearZeroSkew()
        {
            foreach (double c in new[] { 0d, 1E-12, -1E-12, 5E-5, -5E-5, 9.99E-5, -9.99E-5 })
            {
                Assert.AreEqual(double.PositiveInfinity, GammaDistribution.FrequencyFactorKp(c, 1d), $"skew {c}");
                Assert.AreEqual(Normal.StandardZ(0d), GammaDistribution.FrequencyFactorKp(c, 0d), 0d, $"skew {c}");
            }
        }

        /// <summary>
        /// The frequency factor is continuous through zero skew: it returns the Normal variate at zero, and the first-order
        /// term (U² - 1)C/6 at small skew, with no switch.
        /// </summary>
        [TestMethod]
        public void Test_FrequencyFactorKp_SmoothThroughZeroSkew()
        {
            foreach (double p in new[] { 1E-6, 0.01, 0.5, 0.99 })
            {
                double z = Normal.StandardZ(p);
                Assert.AreEqual(z, GammaDistribution.FrequencyFactorKp(0d, p), 0d);
                foreach (double c in new[] { 1E-12, -1E-12, 1E-8, -1E-8 })
                    Assert.AreEqual(z + c * (z * z - 1d) / 6d, GammaDistribution.FrequencyFactorKp(c, p), 1E-15 * (1d + Math.Abs(z)));
            }
        }

        /// <summary>
        /// The L-scale stays finite and accurate for large shapes, where Γ(κ) and Γ(κ + 1/2) overflow separately
        /// beyond κ ≈ 171 (skews at or below 0.1527).
        /// </summary>
        /// <remarks>
        /// Reference values λ2 = θ Γ(κ + 1/2) / (√π Γ(κ)) with θ = 2, computed with Python mpmath 1.3.0 at 50
        /// significant digits from loggamma. The tolerance is 2E-15 relative.
        /// </remarks>
        [TestMethod]
        [DataRow(0.5, 0.6366197723675814)]
        [DataRow(10d, 3.5239410400390625)]
        [DataRow(171d, 14.744690070861031)]
        [DataRow(171.5, 14.766262369935228)]
        [DataRow(172d, 14.787803199723196)]
        [DataRow(1000d, 35.67802229170864)]
        [DataRow(1E6, 1128.3790260481255)]
        public void Test_LinearMomentsFromParameters_LargeShape(double shape, double expectedL2)
        {
            var gamma = new GammaDistribution(2d, shape);
            double[] lmom = gamma.LinearMomentsFromParameters(new[] { 2d, shape });
            Assert.AreEqual(2d * shape, lmom[0], 1E-15 * 2d * shape);
            Assert.AreEqual(expectedL2, lmom[1], 2E-15 * expectedL2);
            Assert.IsTrue(lmom[2] > 0d && lmom[2] < 1d);
            Assert.IsTrue(lmom[3] > 0d && lmom[3] < 1d);
        }
    }
}
