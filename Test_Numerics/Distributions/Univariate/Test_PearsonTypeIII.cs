using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Pearson Type III distribution algorithm.
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
    public class Test_PearsonTypeIII
    {

        // Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
        private double[] sample = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

        /// <summary>
        /// Verification of PIII fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.3 page 70.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_MOM()
        {
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double xi = P3.Xi;
            double beta = P3.Beta;
            double alpha = P3.Alpha;
            double mu = P3.Mu;
            double sigma = P3.Sigma;
            double gamma = P3.Gamma;
            double mean = P3.Mean;
            double stDev = P3.StandardDeviation;
            double skew = P3.Skewness;
            double true_xi = 79.84941d;
            double true_beta = 1d / 0.04846d;
            double true_alpha = 5.40148d;
            double true_mu = 191.31739d;
            double true_sigma = 47.96161d;
            double true_gamma = 0.86055d;
            double true_mean = 191.31739d;
            double true_stDev = 47.96161d;
            double true_skew = 0.86055d;
            Assert.IsLessThan(0.01d, (xi - true_xi) / true_xi );
            Assert.IsLessThan(0.01d, (beta - true_beta) / true_beta );
            Assert.IsLessThan(0.01d, (alpha - true_alpha) / true_alpha);
            Assert.IsLessThan(0.01d, (mu - true_mu) / true_mu);
            Assert.IsLessThan(0.01d, (sigma - true_sigma) / true_sigma);
            Assert.IsLessThan(0.01d, (gamma - true_gamma) / true_gamma);
            Assert.IsLessThan(0.01d, (mean - true_mean) / true_mean);
            Assert.IsLessThan(0.01d, (stDev - true_stDev) / true_stDev);
            Assert.IsLessThan(0.01d, (skew - true_skew) / true_skew);
        }

        /// <summary>
        /// Verification of PIII fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_P3_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = P3.Xi;
            double a = P3.Alpha;
            double b = P3.Beta;
            double true_x = 863.4104d;
            double true_a = 10.02196d;
            double true_b = 78.36751d;
            Assert.AreEqual(true_x, x, 0.001d);
            Assert.AreEqual(true_a, a, 0.001d);
            Assert.AreEqual(true_b, b, 0.001d);
            var lmom = P3.LinearMomentsFromParameters(P3.GetParameters);
            Assert.AreEqual(1648.806d, lmom[0],  0.001d);
            Assert.AreEqual(138.2366d, lmom[1],  0.001d);
            Assert.AreEqual(0.1033889d, lmom[2],  0.001d);
            Assert.AreEqual(0.1258521d, lmom[3],  0.001d);
        }

        /// <summary>
        /// Verification of PIII fit with maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 page 64.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_MLE()
        {
            var P3 = new PearsonTypeIII();
            P3.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double xi = P3.Xi;
            double beta = P3.Beta;
            double alpha = P3.Alpha;
            double mu = P3.Mu;
            double sigma = P3.Sigma;
            double gamma = P3.Gamma;
            double mean = P3.Mean;
            double stDev = P3.StandardDeviation;
            double skew = P3.Skewness;
            double true_xi = 39.38903d;
            double true_beta = 1d / 0.06872d;
            double true_alpha = 10.44062d;
            double true_mu = 191.31739d;
            double true_sigma = 47.01925d;
            double true_gamma = 0.61897d;
            double true_mean = 191.31739d;
            double true_stDev = 47.01925d;
            double true_skew = 0.61897d;
            Assert.IsLessThan(0.01d, (xi - true_xi) / true_xi);
            Assert.IsLessThan(0.01d, (beta - true_beta) / true_beta);
            Assert.IsLessThan(0.01d, (alpha - true_alpha) / true_alpha);
            Assert.IsLessThan(0.01d, (mu - true_mu) / true_mu);
            Assert.IsLessThan(0.01d, (sigma - true_sigma) / true_sigma);
            Assert.IsLessThan(0.01d, (gamma - true_gamma) / true_gamma);
            Assert.IsLessThan(0.01d, (mean - true_mean) / true_mean);
            Assert.IsLessThan(0.01d, (stDev - true_stDev) / true_stDev);
            Assert.IsLessThan(0.01d, (skew - true_skew) / true_skew);
        }

        /// <summary>
        /// Test the quantile function for the Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 page 64.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_Quantile()
        {
            var P3 = new PearsonTypeIII(191.31739d, 47.01925d, -0.61897d);
            double q999 = P3.InverseCDF(0.99d);
            double true_q999 = 321.48d;
            Assert.IsLessThan(0.01d, (q999 - true_q999) / true_q999);
        }

        /// <summary>
        /// Test the standard error for the Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 6.1 &amp; 6.3 page 64-70.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_P3_StandardError()
        {

            // Method of Moments
            var P3 = new PearsonTypeIII(191.31739d, 47.96161d, 0.86055d);
            double qVar99 = Math.Sqrt(P3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MethodOfMoments));
            double true_qVar99 = 27.175d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99);

            // Maximum Likelihood
            P3 = new PearsonTypeIII(191.31739d, 47.01925d, 0.61897d);
            qVar99 = Math.Sqrt(P3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MaximumLikelihood));
            true_qVar99 = 20.045d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99 );
        }

        /// <summary>
        /// Verifying input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(100,P3.Mu);
            Assert.AreEqual(10,P3.Sigma);
            Assert.AreEqual(0, P3.Gamma);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(1,P3ii.Mu);
            Assert.AreEqual(1, P3ii.Sigma);
            Assert.AreEqual(1, P3ii.Gamma);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var P3 = new PearsonTypeIII(double.PositiveInfinity, double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(P3.ParametersValid);

            var P3ii = new PearsonTypeIII(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(P3ii.ParametersValid);

            var P3iii = new PearsonTypeIII(1, 0, 1);
            Assert.IsFalse(P3iii.ParametersValid);
        }

        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual("Mean (µ)",P3.ParametersToString[0, 0]);
            Assert.AreEqual("Std Dev (σ)", P3.ParametersToString[1, 0]);
            Assert.AreEqual("Skew (γ)", P3.ParametersToString[2, 0]);
            Assert.AreEqual("100", P3.ParametersToString[0, 1]);
            Assert.AreEqual("10", P3.ParametersToString[1, 1]);
            Assert.AreEqual("0", P3.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new PearsonTypeIII();
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
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(100, P3.Mean);

            var P3ii = new PearsonTypeIII(100, 1, 1);
            Assert.AreEqual(100, P3ii.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(100, P3.Median);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(100, P3.Mode);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(0.5, P3ii.Mode);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(10, P3.StandardDeviation);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(1, P3ii.StandardDeviation);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(0, P3.Skewness);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(0, P3.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(3, P3.Kurtosis);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(4.5, P3ii.Kurtosis);
        }

        /// <summary>
        /// Testing minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(double.NegativeInfinity,P3.Minimum);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(-1, P3ii.Minimum);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var P3 = new PearsonTypeIII();
            Assert.AreEqual(double.PositiveInfinity, P3.Maximum);

            var P3ii = new PearsonTypeIII(1, 1, 1);
            Assert.AreEqual(double.PositiveInfinity, P3ii.Maximum);
        }

        /// <summary>
        /// Verifies signed, finite Pearson III L-moments and their Normal limit.
        /// </summary>
        [TestMethod]
        public void Test_LinearMoments_SignedSmallAndZeroSkew()
        {
            var distribution = new PearsonTypeIII();
            double[] normalMoments = distribution.LinearMomentsFromParameters([10.0d, 2.0d, 0.0d]);
            Assert.AreEqual(10.0d, normalMoments[0]);
            Assert.AreEqual(2.0d / Math.Sqrt(Math.PI), normalMoments[1], 1E-14);
            Assert.AreEqual(0.0d, normalMoments[2]);
            Assert.AreEqual(0.12260172d, normalMoments[3], 1E-14);

            double[] normalParameters = distribution.ParametersFromLinearMoments(normalMoments);
            Assert.AreEqual(10.0d, normalParameters[0]);
            Assert.AreEqual(2.0d, normalParameters[1], 1E-14);
            Assert.AreEqual(0.0d, normalParameters[2]);

            double[] positive = distribution.LinearMomentsFromParameters([10.0d, 2.0d, 0.1d]);
            double[] negative = distribution.LinearMomentsFromParameters([10.0d, 2.0d, -0.1d]);
            Assert.AreEqual(10.0d, positive[0]);
            Assert.AreEqual(positive[0], negative[0]);
            Assert.AreEqual(positive[1], negative[1], 1E-14);
            Assert.AreEqual(positive[2], -negative[2], 1E-14);
            Assert.IsGreaterThan(0.0d, positive[2]);

            foreach (double moment in positive)
            {
                Assert.IsFalse(double.IsNaN(moment));
                Assert.IsFalse(double.IsInfinity(moment));
            }

            double[] recovered = distribution.ParametersFromLinearMoments(negative);
            Assert.AreEqual(10.0d, recovered[0]);
            Assert.AreEqual(2.0d, recovered[1], 1E-4);
            Assert.AreEqual(-0.1d, recovered[2], 1E-4);

            double expectedDerivative = 2.0d * GammaDistribution.PartialKp(0.0d, 0.9d);
            var zeroSkew = new PearsonTypeIII(10.0d, 2.0d, 0.0d);
            Assert.AreEqual(expectedDerivative, zeroSkew.QuantileGradientForMoments(0.9d)[2], 1E-14);
        }

        /// <summary>
        /// Pearson III shapes below one have a genuine density singularity at the support
        /// boundary for both skew directions.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF_ShapeBelowOneIsPositiveInfinityAtLocation()
        {
            var positiveSkew = new PearsonTypeIII(1.0, 1.5, 3.0);
            var negativeSkew = new PearsonTypeIII(-1.0, 1.5, -3.0);

            Assert.AreEqual(0.0, positiveSkew.Xi, 0.0);
            Assert.AreEqual(0.0, negativeSkew.Xi, 0.0);
            Assert.AreEqual(double.PositiveInfinity, positiveSkew.LogPDF(positiveSkew.Xi));
            Assert.AreEqual(double.PositiveInfinity, negativeSkew.LogPDF(negativeSkew.Xi));
        }


        /// <summary>
        /// The L-scale and the L-moment scale estimate use the exact gamma ratio for every shape, so they stay accurate
        /// across alpha = 4/skew^2 = 100 and finite beyond alpha = 171 (skews at or below 0.1527).
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.3.0 at 50 significant digits:
        /// L2 = sigma Gamma(alpha + 1/2) / (sqrt(pi alpha) Gamma(alpha)) with sigma = 0.3, and the scale
        /// sigma = L2 sqrt(pi alpha) Gamma(alpha) / Gamma(alpha + 1/2) at the fitted skew with L2 = 0.17. The tolerance
        /// is 3E-15 relative.
        /// </remarks>
        [TestMethod]
        [DataRow(2.0, 0.15)]
        [DataRow(0.5, 0.16793992090970278)]
        [DataRow(0.2, 0.16904543702776925)]
        [DataRow(0.19, 0.1690660404626627)]
        [DataRow(0.1527, 0.16913358854788108)]
        [DataRow(0.15, 0.16913790831003775)]
        [DataRow(-0.1, 0.16920399056827415)]
        [DataRow(0.01, 0.16925634613741877)]
        [DataRow(0.002, 0.1692568539072188)]
        public void Test_LinearMomentsFromParameters_ExactGammaRatio(double skew, double expectedL2)
        {
            double[] lmom = new PearsonTypeIII().LinearMomentsFromParameters(new[] { 3d, 0.3, skew });
            Assert.AreEqual(expectedL2, lmom[1], 3E-15 * expectedL2);
        }

        /// <summary>
        /// The L-moment scale estimate uses the exact gamma ratio at the fitted skew.
        /// </summary>
        /// <remarks>Reference values as for <see cref="Test_LinearMomentsFromParameters_ExactGammaRatio"/>.</remarks>
        [TestMethod]
        [DataRow(0.05, 0.30663370214036545, 0.3022037834482887)]
        [DataRow(0.03, 0.1841192827432352, 0.30163652985276324)]
        [DataRow(0.0249, 0.1528394551987608, 0.3015371952184198)]
        [DataRow(0.01, 0.0613966415922153, 0.3013526514107132)]
        [DataRow(-0.005, -0.030699430945504754, 0.3013260290932623)]
        [DataRow(0.001, 0.006139957284868968, 0.30131750963471404)]
        public void Test_ParametersFromLinearMoments_ExactGammaRatio(double t3, double expectedSkew, double expectedSigma)
        {
            double[] p = new PearsonTypeIII().ParametersFromLinearMoments(new[] { 3d, 0.17, t3 });
            Assert.AreEqual(3d, p[0], 0d);
            Assert.AreEqual(expectedSkew, p[2], 1E-15 * Math.Abs(expectedSkew));
            Assert.AreEqual(expectedSigma, p[1], 3E-15 * expectedSigma);
        }

        /// <summary>
        /// The log tails and log density carry the rounding residual of the centered unit gamma coordinate below absolute
        /// skew 0.1, so they keep near full accuracy where the uncentered coordinate is in error by up to 1.5E-13.
        /// </summary>
        /// <remarks>
        /// Columns: skew, x (the 1E-6, 0.5 and 0.99 quantiles of P3(0, 1, skew)), ln F(x), ln(1 - F(x)) and ln f(x),
        /// computed with mpmath 1.3.0 at 60 digits by quadrature of the gamma density (tails) and closed form (density).
        /// </remarks>
        [TestMethod]
        [DataRow(0.05, -4.574785358853457, -13.815510557964272, -1.0000005000003347e-06, -12.135095893353434)]
        [DataRow(0.05, -0.008333024620929786, -0.6931471805599453, -0.6931471805599453, -0.9188169936914243)]
        [DataRow(0.05, 2.363047431519595, -0.010050335853501454, -4.60517018598809, -3.663076425606502)]
        [DataRow(-0.05, -4.934638161479374, -13.815510557964272, -1.0000005000003353e-06, -12.293655407637814)]
        [DataRow(-0.05, 0.008333024620929786, -0.6931471805599453, -0.6931471805599453, -0.9188169936914243)]
        [DataRow(-0.05, 2.289520092650801, -0.010050335853501456, -4.60517018598809, -3.5855121536914183)]
        [DataRow(0.03, -4.6459194427364725, -13.815510557964274, -1.0000005000003344e-06, -12.168163508231066)]
        [DataRow(0.03, -0.004999933327857031, -0.6931471805599453, -0.6931471805599453, -0.9188947816838093)]
        [DataRow(0.03, 2.3483838020219627, -0.010050335853501452, -4.60517018598809, -3.64793842524562)]
        [DataRow(-0.03, -4.8618559211873125, -13.815510557964274, -1.000000500000333e-06, -12.263256188295664)]
        [DataRow(-0.03, 0.004999933327857031, -0.6931471805599453, -0.6931471805599453, -0.9188947816838093)]
        [DataRow(-0.03, 2.304265772351977, -0.010050335853501445, -4.605170185988091, -3.6014072902489267)]
        public void Test_LogTailsAndDensity_CenteredBelowSkew0_1(double skew, double x, double logCdf, double logCcdf, double logPdf)
        {
            var d = new PearsonTypeIII(0d, 1d, skew);
            Assert.AreEqual(logCdf, d.LogCDF(x), 1E-14);
            Assert.AreEqual(logCcdf, d.LogCCDF(x), 1E-14 * Math.Abs(logCcdf) + 1E-20);
            Assert.AreEqual(logPdf, d.LogPDF(x), 1E-14);
        }

        /// <summary>
        /// The maximum-likelihood skew information uses its Bernoulli series through B36 below skew²/4 = 0.15 and the
        /// trigamma form above it.
        /// </summary>
        /// <remarks>
        /// The reference is 1/T with T = a³[ψ'(a) - 1/(a-1) + 1/(2(a-1)²)], a = 4/skew², from mpmath 1.3.0 at 50 digits;
        /// the covariance of the skew at sample size one is 1/T. The tolerance is 1E-13 relative; the trigamma form alone is
        /// in error by up to 3E-13 for skew²/4 between 0.05 and 0.15.
        /// </remarks>
        [TestMethod]
        [DataRow(0.3, 5.604637825817938)]
        [DataRow(0.5, 4.948232207631228)]
        [DataRow(0.7, 4.0697304381664265)]
        [DataRow(0.77, 3.7300037226093092)]
        [DataRow(0.78, 3.680673250778186)]
        [DataRow(1.0, 2.584704719215949)]
        public void Test_ParameterCovariance_SkewInformation(double skew, double inverseInformation)
        {
            var c = new PearsonTypeIII(0d, 1d, skew).ParameterCovariance(1, ParameterEstimationMethod.MaximumLikelihood);
            Assert.AreEqual(inverseInformation, c[2, 2], 1E-13 * inverseInformation);
        }

        /// <summary>
        /// Within absolute skew 0.03 and |skew|(1 + |z|) at most 0.2 the quantile and its skew derivative come from the
        /// order-8 skew series, which removes the cancellation of the gamma quantile against its shape.
        /// </summary>
        /// <remarks>
        /// Columns: skew, probability, the standardized quantile and its skew derivative, from mpmath 1.3.0 at 60 digits by
        /// Newton iteration on the quadrature of the gamma density and by implicit differentiation of the tail. On this range
        /// the gamma path alone is in error by up to 2.7E-13 in the quantile and 3.4E-8 relative in the derivative.
        /// </remarks>
        [TestMethod]
        [DataRow(1e-12, 1e-06, -4.7534243088193, 3.599173776617046)]
        [DataRow(1e-06, 0.01, -2.326347138725077, 0.7353157898213117)]
        [DataRow(-0.0001, 0.99, 2.3262743422104455, 0.735320869231516)]
        [DataRow(0.000999, 1e-06, -4.74982924823925, 3.5981444512326153)]
        [DataRow(0.001001, 1e-06, -4.7498220519524095, 3.5981423889685757)]
        [DataRow(-0.002, 0.5, 0.0003333333135802397, -0.166666637037019)]
        [DataRow(0.005, 0.999999, 4.77143301522681, 3.604302348820344)]
        [DataRow(-0.01, 1e-06, -4.789467268116546, 3.6093924042486107)]
        [DataRow(0.02, 0.01, -2.311631433034089, 0.7363215141260572)]
        [DataRow(-0.0299, 0.99, 2.304339452999121, 0.7368040603980558)]
        [DataRow(-0.016, 1e-30, -11.814324379694858, 22.048097400863778)]
        public void Test_Quantile_SkewSeries(double skew, double probability, double quantile, double derivative)
        {
            var d = new PearsonTypeIII(0d, 1d, skew);
            Assert.AreEqual(quantile, d.InverseCDF(probability), 4E-16 * (1d + Math.Abs(quantile)));
            Assert.AreEqual(derivative, d.QuantileGradient(probability)[2], 1E-13 * Math.Abs(derivative));
            var shifted = new PearsonTypeIII(3d, 0.3, skew);
            Assert.AreEqual(3d + 0.3 * quantile, shifted.InverseCDF(probability), 1E-15 * (1d + Math.Abs(3d + 0.3 * quantile)));
            Assert.AreEqual(0.3 * derivative, shifted.QuantileGradient(probability)[2], 1E-13 * Math.Abs(0.3 * derivative));
        }

        /// <summary>
        /// Just outside the skew-series region the gamma path takes over, continuous with the series to about 1E-14 in the
        /// quantile and 3E-12 relative in the skew derivative.
        /// </summary>
        /// <remarks>Reference values as for <see cref="Test_Quantile_SkewSeries"/>.</remarks>
        [TestMethod]
        [DataRow(0.0301, 0.5, -0.005016599325540457, -0.16665995455652624)]
        [DataRow(0.03, 1e-10, -6.1653599779701125, 6.486885714447091)]
        [DataRow(0.05, 1e-06, -4.574785358853457, 3.5457194876203704)]
        public void Test_Quantile_GammaPathBeyondSkewSeries(double skew, double probability, double quantile, double derivative)
        {
            var d = new PearsonTypeIII(0d, 1d, skew);
            Assert.AreEqual(quantile, d.InverseCDF(probability), 1E-14 * (1d + Math.Abs(quantile)));
            Assert.AreEqual(derivative, d.QuantileGradient(probability)[2], 5E-12 * Math.Abs(derivative));
        }
    }


}
