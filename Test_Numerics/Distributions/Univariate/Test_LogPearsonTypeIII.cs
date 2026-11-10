using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Log-Pearson Type III distribution algorithm.
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
    public class Test_LogPearsonTypeIII
    {

        /// <summary>
        /// Verifies cloning preserves the configured logarithm base and resulting distribution.
        /// </summary>
        [TestMethod]
        public void Test_Clone_PreservesBase()
        {
            var source = new LogPearsonTypeIII(4.2d, 0.4d, 0.25d) { Base = Math.E };

            var clone = (LogPearsonTypeIII)source.Clone();

            Assert.AreNotSame(source, clone);
            Assert.AreEqual(source.Mu, clone.Mu, 0d);
            Assert.AreEqual(source.Sigma, clone.Sigma, 0d);
            Assert.AreEqual(source.Gamma, clone.Gamma, 0d);
            Assert.AreEqual(source.Base, clone.Base, 0d);
            Assert.AreEqual(source.CDF(75d), clone.CDF(75d), 0d);
        }

        /// <summary>
        /// Verifies bootstrap distributions retain the configured logarithm base.
        /// </summary>
        [TestMethod]
        public void Test_Bootstrap_PreservesBase()
        {
            var source = new LogPearsonTypeIII(4.2d, 0.4d, 0.25d) { Base = Math.E };

            var bootstrap = (LogPearsonTypeIII)source.Bootstrap(
                ParameterEstimationMethod.MethodOfMoments, 40, 12345);

            Assert.AreEqual(source.Base, bootstrap.Base, 0d);
        }

        // Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee & F. Ashkar, Water Resources Publications, 1991.
        // Table 1.2 Maximum annual peak discharge values in cms, observed at the Harricana River at Amos (Quebec, Canada)
        private double[] sample = new double[] { 122d, 244d, 214d, 173d, 229d, 156d, 212d, 263d, 146d, 183d, 161d, 205d, 135d, 331d, 225d, 174d, 98.8d, 149d, 238d, 262d, 132d, 235d, 216d, 240d, 230d, 192d, 195d, 172d, 173d, 172d, 153d, 142d, 317d, 161d, 201d, 204d, 194d, 164d, 183d, 161d, 167d, 179d, 185d, 117d, 192d, 337d, 125d, 166d, 99.1d, 202d, 230d, 158d, 262d, 154d, 164d, 182d, 164d, 183d, 171d, 250d, 184d, 205d, 237d, 177d, 239d, 187d, 180d, 173d, 174d };

        /// <summary>
        /// Verification of LPIII fit with indirect method of moments.
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
        public void Test_LP3_IndirectMOM()
        {
            var LP3 = new LogPearsonTypeIII();
            LP3.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double xi = LP3.Xi;
            double beta = LP3.Beta;
            double alpha = LP3.Alpha;
            double meanOfLog = LP3.Mu;
            double stDevOfLog = LP3.Sigma;
            double skewOfLog = LP3.Gamma;
            double mean = LP3.Mean;
            double stDev = LP3.StandardDeviation;
            double skew = LP3.Skewness;
            double true_xi = 7.53821d;
            double true_beta = 1d / -460.31089d;
            double true_alpha = 2425.57481d;
            double true_meanOfLog = 2.26878d;
            double true_stDevOfLog = 0.10699d;
            double true_skewOfLog = -0.04061d;
            double true_mean = 191.38768d;
            double true_stDev = 47.62977d;
            double true_skew = 0.71589d;
            Assert.IsLessThan(0.01d, (xi - true_xi) / true_xi );
            Assert.IsLessThan(0.01d, (beta - true_beta) / true_beta );
            Assert.IsLessThan(0.01d, (alpha - true_alpha) / true_alpha);
            Assert.IsLessThan(0.01d, (meanOfLog - true_meanOfLog) / true_meanOfLog);
            Assert.IsLessThan(0.01d, (stDevOfLog - true_stDevOfLog) / true_stDevOfLog);
            Assert.IsLessThan(0.01d, (skewOfLog - true_skewOfLog) / true_skewOfLog);
            Assert.IsLessThan(0.01d, (mean - true_mean) / true_mean);
            Assert.IsLessThan(0.01d, (stDev - true_stDev) / true_stDev);
            Assert.IsLessThan(0.01d, (skew - true_skew) / true_skew);
        }

        /// <summary>
        /// Verification of LPIII fit with maximum likelihood.
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
        public void Test_LP3_MLE()
        {
            var LP3 = new LogPearsonTypeIII();
            LP3.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double xi = LP3.Xi;
            double beta = LP3.Beta;
            double alpha = LP3.Alpha;
            double meanOfLog = LP3.Mu;
            double stDevOfLog = LP3.Sigma;
            double skewOfLog = LP3.Gamma;
            double mean = LP3.Mean;
            double stDev = LP3.StandardDeviation;
            double skew = LP3.Skewness;
            double true_xi = 9.53033d;
            double true_beta = 1d / -643.69408d;
            double true_alpha = 4674.2179d;
            double true_meanOfLog = 2.26878d;
            double true_stDevOfLog = 0.10621d;
            double true_skewOfLog = -0.02925d;
            double true_mean = 191.30891d;
            double true_stDev = 47.32124d;
            double true_skew = 0.72396d;
            Assert.IsLessThan(0.01d, (xi - true_xi) / true_xi);
            Assert.IsLessThan(0.01d, (beta - true_beta) / true_beta);
            Assert.IsLessThan(0.01d, (alpha - true_alpha) / true_alpha);
            Assert.IsLessThan(0.01d, (meanOfLog - true_meanOfLog) / true_meanOfLog);
            Assert.IsLessThan(0.01d, (stDevOfLog - true_stDevOfLog) / true_stDevOfLog);
            Assert.IsLessThan(0.01d, (skewOfLog - true_skewOfLog) / true_skewOfLog);
            Assert.IsLessThan(0.01d, (mean - true_mean) / true_mean);
            Assert.IsLessThan(0.01d, (stDev - true_stDev) / true_stDev);
            Assert.IsLessThan(0.01d, (skew - true_skew) / true_skew);
        }

        /// <summary>
        /// Test the quantile function for the Log-Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 7.1 page 87.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LP3_Quantile()
        {
            var LP3 = new LogPearsonTypeIII(2.26878d, 0.10621d, -0.02925d);
            double q1000 = LP3.InverseCDF(0.99d);
            double true_q1000 = 326.25d;
            Assert.IsLessThan(0.01d, (q1000 - true_q1000) / true_q1000);
            double p = LP3.CDF(q1000);
            double true_p = 0.99d;
            Assert.IsLessThan(0.01d, (p - true_p) / true_p);
        }

        /// <summary>
        /// Test the standard error for the Log-Pearson Type III Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "The Gamma Family and Derived Distributions Applied in Hydrology", B. Bobee &amp; F. Ashkar, Water Resources Publications, 1991.
        /// </para>
        /// <para>
        /// Example 7.1 &amp; 6.4 page 87-93.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_LP3_StandardError()
        {

            // Method of Moments
            var LP3 = new LogPearsonTypeIII(2.26878d, 0.10699d, -0.04061d);
            double qVar99 = Math.Sqrt(LP3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MethodOfMoments));
            double true_qVar99 = 25.053d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99);

            // Maximum Likelihood
            LP3 = new LogPearsonTypeIII(2.26878d, 0.10621d, -0.02925d);
            qVar99 = Math.Sqrt(LP3.QuantileVariance(0.99d, 69, ParameterEstimationMethod.MaximumLikelihood));
            true_qVar99 = 25d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99 );

        }

        /// <summary>
        /// Verifying that input parameters can create distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(3,LP3.Mu);
            Assert.AreEqual(0.5,LP3.Sigma);
            Assert.AreEqual(0, LP3.Gamma);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var LP3 = new LogPearsonTypeIII(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(LP3.ParametersValid);

            var LP3ii = new LogPearsonTypeIII(double.PositiveInfinity, double.NegativeInfinity, double.PositiveInfinity);
            Assert.IsFalse(LP3ii.ParametersValid);
        }

        /// <summary>
        /// Testing parameter to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual("Mean (of log) (µ)", LP3.ParametersToString[0, 0]);
            Assert.AreEqual("Std Dev (of log) (σ)", LP3.ParametersToString[1, 0]);
            Assert.AreEqual("Skew (of log) (γ)", LP3.ParametersToString[2, 0]);
            Assert.AreEqual("3", LP3.ParametersToString[0, 1]);
            Assert.AreEqual("0.5", LP3.ParametersToString[1, 1]);
            Assert.AreEqual("0", LP3.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var trueDist = new LnNormal(10, 5);
            var dist = new LogPearsonTypeIII(trueDist.Mu, trueDist.Sigma, 0.0) { Base = Math.E };
            Assert.AreEqual(trueDist.Mean, dist.Mean, 1E-2);
            Assert.AreEqual(trueDist.StandardDeviation, dist.StandardDeviation, 1E-2);
            Assert.AreEqual(trueDist.Skewness, dist.Skewness, 1E-2);
            Assert.AreEqual(trueDist.Kurtosis, dist.Kurtosis, 1E-2);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(1000, LP3.Median,  1e-04);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var LP3 = new LogPearsonTypeIII();
            // Frozen defining-formula evidence: docs/distributions/oracles/normal-pearson.R.
            Assert.AreEqual(265.67685816501944, LP3.Mode, 1e-10);

            var LP3ii = new LogPearsonTypeIII(1, 1, 1);
            Assert.AreEqual(0.4980296964878923, LP3ii.Mode, 1e-13);
        }

        /// <summary>
        /// Testing minimum.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(0, LP3.Minimum );

            var LP3ii = new LogPearsonTypeIII(1,1,1);
            Assert.AreEqual(0.1, LP3ii.Minimum,  1e-05);

            // A bounded upper tail requires negative skew and a positive standard deviation.
            var LP3iii = new LogPearsonTypeIII(1, 1, -1);
            Assert.AreEqual(0,LP3iii.Minimum);
        }

        /// <summary>
        /// Testing maximum.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(double.PositiveInfinity,LP3.Maximum );

            var LP3ii = new LogPearsonTypeIII(1,1,1);
            Assert.AreEqual(double.PositiveInfinity, LP3ii.Maximum);

            var LP3iii = new LogPearsonTypeIII(1, 1, -1);
            Assert.AreEqual(1000, LP3iii.Maximum, 1e-04);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(0, LP3.PDF(-1));
            Assert.AreEqual(5.2774e-09, LP3.PDF(1),1e-13);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(0,LP3.CDF(-1));
            Assert.AreEqual(9.8658e-10, LP3.CDF(1), 1e-13);
        }

        /// <summary>
        /// Testing inverse CDF.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(0,LP3.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,LP3.InverseCDF(1));
            Assert.AreEqual(546.7637, LP3.InverseCDF(0.3), 1e-04);
        }

        /// <summary>
        /// Testing Wilson-Hilferty Inverse CDF.
        /// </summary>
        [TestMethod()]
        public void ValidateWilsonHilfertyInverseCDF()
        {
            var LP3 = new LogPearsonTypeIII();
            Assert.AreEqual(0, LP3.WilsonHilfertyInverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity, LP3.WilsonHilfertyInverseCDF(1));
            Assert.AreEqual(747.01005, LP3.WilsonHilfertyInverseCDF(0.4), 1e-05);
        }

        /// <summary>
        /// Verifies signed, finite Log-Pearson III L-moments and their Normal limit.
        /// </summary>
        [TestMethod]
        public void Test_LinearMoments_SignedSmallAndZeroSkew()
        {
            var distribution = new LogPearsonTypeIII();
            double[] normalMoments = distribution.LinearMomentsFromParameters([2.0d, 0.3d, 0.0d]);
            Assert.AreEqual(2.0d, normalMoments[0]);
            Assert.AreEqual(0.3d / Math.Sqrt(Math.PI), normalMoments[1], 1E-14);
            Assert.AreEqual(0.0d, normalMoments[2]);
            Assert.AreEqual(0.12260172d, normalMoments[3], 1E-14);

            double[] normalParameters = distribution.ParametersFromLinearMoments(normalMoments);
            Assert.AreEqual(2.0d, normalParameters[0]);
            Assert.AreEqual(0.3d, normalParameters[1], 1E-14);
            Assert.AreEqual(0.0d, normalParameters[2]);

            double[] positive = distribution.LinearMomentsFromParameters([2.0d, 0.3d, 0.1d]);
            double[] negative = distribution.LinearMomentsFromParameters([2.0d, 0.3d, -0.1d]);
            Assert.AreEqual(2.0d, positive[0]);
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
            Assert.AreEqual(2.0d, recovered[0]);
            Assert.AreEqual(0.3d, recovered[1], 1E-5);
            Assert.AreEqual(-0.1d, recovered[2], 1E-4);
        }
        /// <summary>
        /// Verify the parameter constraints admit a negative log10-space mean.
        /// </summary>
        /// <remarks>
        /// See the matching LogNormal test: the location parameter is the mean of the
        /// log10-transformed data, negative whenever the data are mostly below 1; a machine-epsilon
        /// lower bound would reject any such sample before a fit could start.
        /// </remarks>
        [TestMethod]
        public void Test_LP3_ParameterConstraints_AllowNegativeLogMean()
        {
            var sample = new double[] { 0.12, 0.31, 0.45, 0.08, 0.90, 1.4, 0.25, 0.6, 0.5, 0.75, 0.2, 0.33 };
            var constraints = new LogPearsonTypeIII().GetParameterConstraints(sample);
            var initials = constraints.Item1;
            var lowers = constraints.Item2;
            var uppers = constraints.Item3;
            Assert.IsLessThan(0d, initials[0], "Fixture precondition: the log10 mean is negative.");
            Assert.IsLessThan(uppers[0], lowers[0], "The mean bounds must not be inverted.");
            Assert.IsTrue(initials[0] >= lowers[0] && initials[0] <= uppers[0],
                "The initial mean must sit inside its own bounds.");
            Assert.AreEqual(double.NegativeInfinity, new LogPearsonTypeIII().MinimumOfParameters[0]);
        }

        /// <summary>
        /// Log-Pearson III shapes below one have a genuine density singularity at the
        /// transformed support boundary for both skew directions.
        /// </summary>
        [TestMethod]
        public void Test_LogPDF_ShapeBelowOneIsPositiveInfinityAtLocation()
        {
            var positiveSkew = new LogPearsonTypeIII(1.0, 1.5, 3.0);
            var negativeSkew = new LogPearsonTypeIII(-1.0, 1.5, -3.0);

            Assert.AreEqual(0.0, positiveSkew.Xi, 0.0);
            Assert.AreEqual(0.0, negativeSkew.Xi, 0.0);
            Assert.AreEqual(double.PositiveInfinity, positiveSkew.LogPDF(1.0));
            Assert.AreEqual(double.PositiveInfinity, negativeSkew.LogPDF(1.0));
        }

        /// <summary>
        /// The support endpoints and the log density, log CDF and log CCDF equal the expressions evaluated with
        /// Math.Log(Base), Math.Log(Math.Log(Base)) and Math.Exp(Xi * Math.Log(Base)) on every call, bit for bit,
        /// including the endpoints, after every parameter setter, a base change, SetParameters and Clone.
        /// </summary>
        [TestMethod()]
        public void Test_CachedBaseLogarithms()
        {
            var dist = new LogPearsonTypeIII(1.2, 0.3, 0.5);
            var steps = new System.Action[]
            {
                () => { },
                () => dist.Gamma = -0.8,
                () => dist.Mu = 2.1,
                () => dist.Sigma = 0.05,
                () => dist.Base = Math.E,
                () => dist.Gamma = 0d,
                () => dist.SetParameters(0.4, 0.2, 3.5),
                () => dist.Base = 2d,
                () => dist.SetParameters(new[] { 3d, 0.4, -2.5 }),
                () => dist.Base = 10d,
            };
            foreach (var step in steps)
            {
                step();
                foreach (var d in new[] { dist, (LogPearsonTypeIII)dist.Clone() })
                {
                    double logBase = Math.Log(d.Base);
                    double boundary = d.Gamma == 0d ? double.NaN : Math.Exp(d.Xi * logBase);
                    double minimum = d.Gamma > 0d ? Math.Exp(d.Xi * Math.Log(d.Base)) : 0d;
                    double maximum = d.Gamma < 0d ? Math.Exp(d.Xi * Math.Log(d.Base)) : double.PositiveInfinity;
                    Assert.AreEqual(minimum, d.Minimum, 0.0);
                    Assert.AreEqual(maximum, d.Maximum, 0.0);
                    var grid = new System.Collections.Generic.List<double> { 0d, 1E-300, 0.5, 1d, 3d, 15d, 100d, 1E4, 1E300, double.PositiveInfinity };
                    if (!double.IsNaN(boundary)) { grid.Add(boundary); grid.Add(boundary * (1 + 1E-12)); grid.Add(boundary * (1 - 1E-12)); }
                    foreach (double x in grid)
                    {
                        double expected;
                        if (x < minimum || x > maximum || double.IsPositiveInfinity(x)) expected = double.NegativeInfinity;
                        else if (x == 0d)
                        {
                            double rate = -1d / (d.Beta * logBase);
                            expected = d.Gamma >= 0d || rate > 1d ? double.NegativeInfinity : rate < 1d || d.Alpha > 1d ? double.PositiveInfinity
                                : d.Alpha == 1d ? -d.Xi * logBase : double.NegativeInfinity;
                        }
                        else
                        {
                            double logX = Math.Log(x);
                            double transformed = (d.Gamma > 0d && x == boundary) || (d.Gamma < 0d && x == boundary) ? d.Xi : logX / logBase;
                            expected = PearsonTypeIII.LogPDF(d.Mu, d.Sigma, d.Gamma, transformed) - Math.Log(logBase) - logX;
                        }
                        Assert.AreEqual(expected, d.LogPDF(x), 0.0, $"LogPDF x {x:R}");
                        double expectedCdf = x <= minimum ? double.NegativeInfinity : x >= maximum ? 0d : PearsonTypeIII.LogTail(d.Mu, d.Sigma, d.Gamma, Math.Log(x, d.Base), false);
                        Assert.AreEqual(expectedCdf, d.LogCDF(x), 0.0, $"LogCDF x {x:R}");
                        double expectedCcdf = x <= minimum ? 0d : x >= maximum ? double.NegativeInfinity : PearsonTypeIII.LogTail(d.Mu, d.Sigma, d.Gamma, Math.Log(x, d.Base), true);
                        Assert.AreEqual(expectedCcdf, d.LogCCDF(x), 0.0, $"LogCCDF x {x:R}");
                    }
                }
            }
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
            double[] lmom = new LogPearsonTypeIII().LinearMomentsFromParameters(new[] { 3d, 0.3, skew });
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
            double[] p = new LogPearsonTypeIII().ParametersFromLinearMoments(new[] { 3d, 0.17, t3 });
            Assert.AreEqual(3d, p[0], 0d);
            Assert.AreEqual(expectedSkew, p[2], 1E-15 * Math.Abs(expectedSkew));
            Assert.AreEqual(expectedSigma, p[1], 3E-15 * expectedSigma);
        }

        /// <summary>
        /// The base-log L-moment members are those of Pearson type III.
        /// </summary>
        [TestMethod]
        public void Test_LinearMoments_DelegateToPearsonTypeIII()
        {
            foreach (double skew in new[] { -2.5, -0.4, 0d, 1E-9, 0.3, 1.7 })
            {
                var parameters = new[] { 3d, 0.25, skew };
                CollectionAssert.AreEqual(new PearsonTypeIII().LinearMomentsFromParameters(parameters), new LogPearsonTypeIII().LinearMomentsFromParameters(parameters));
            }
            foreach (double t3 in new[] { -0.6, -0.1, 0d, 0.05, 0.5 })
            {
                var moments = new[] { 3d, 0.15, t3, 0.15 };
                CollectionAssert.AreEqual(new PearsonTypeIII().ParametersFromLinearMoments(moments), new LogPearsonTypeIII().ParametersFromLinearMoments(moments));
            }
        }

        /// <summary>
        /// The real-space skewness and kurtosis use the cancellation-free form below
        /// sigma·ln(10)·(1 + |skew|) = 0.15, where the closed form is in error by up to 1.7E-9 relative.
        /// </summary>
        /// <remarks>
        /// Columns: sigma, skew, skewness and kurtosis of LogPearsonTypeIII(3, sigma, skew) in base 10, from the raw moments
        /// E[10^(rX)] = exp(r ln10 µ + α[-ln(1 - u) - u]), u = r ln10 sigma skew / 2, with mpmath 1.3.0 at 60 digits.
        /// </remarks>
        [TestMethod]
        [DataRow(0.00178, 2.0, 2.0248465835835114, 9.200633577802794)]
        [DataRow(0.00178, -2.0, -1.9756574439817294, 8.807028992529881)]
        [DataRow(0.0237, 1.5, 1.785720434623108, 8.318876023644016)]
        [DataRow(0.01, 0.5, 0.5744908860077079, 3.5347293989594415)]
        [DataRow(0.03, -0.3, -0.0928908070637434, 2.9640568957578286)]
        [DataRow(0.05, 0.1, 0.45409947825904245, 3.3829676148063212)]
        [DataRow(0.004, 0.0, 0.027632388480010317, 3.0013574514431385)]
        public void Test_SkewnessKurtosis_SmallLogScale(double sigma, double skew, double skewness, double kurtosis)
        {
            var d = new LogPearsonTypeIII(3d, sigma, skew);
            Assert.AreEqual(skewness, d.Skewness, 1E-14 * Math.Abs(skewness));
            Assert.AreEqual(kurtosis, d.Kurtosis, 1E-14 * kurtosis);
        }
    }
}
