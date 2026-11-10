using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using System.Diagnostics;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Empirical distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_EmpiricalDistribution
    {
        /// <summary>
        /// This test creates an empirical distribution that is designed to closely match a Normal distribution. 
        /// The test compares the moments of the two distributions.
        /// </summary>
       [TestMethod]
        public void Test_Empirical_Normal_Moments()
        {
            var norm = new Normal(100, 15);
            var eCDF = norm.CreateCDFGraph();
            var emp  = new EmpiricalDistribution(eCDF.GetColumn(0), eCDF.GetColumn(1));

            // Test moments
            Assert.AreEqual(norm.Mean, emp.Mean, 1E-1);
            Assert.AreEqual(norm.StandardDeviation, emp.StandardDeviation, 1E-1);
            Assert.AreEqual(norm.Skewness, emp.Skewness, 1E-1);
            Assert.AreEqual(norm.Kurtosis, emp.Kurtosis, 1E-1);

        }

        /// <summary>
        /// This test creates an empirical distribution that is designed to closely match a Normal distribution. 
        /// The test compares the distribution functions of the two distributions.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Normal_Dist()
        {
            var norm = new Normal(100, 15);
            var eCDF = norm.CreateCDFGraph();
            var emp = new EmpiricalDistribution(eCDF.GetColumn(0), eCDF.GetColumn(1));

            // Test distribution functions
            Assert.AreEqual(norm.PDF(80), emp.PDF(80), 1E-4);
            Assert.AreEqual(norm.CDF(80), emp.CDF(80), 1E-4);
            Assert.AreEqual(norm.CCDF(80), emp.CCDF(80), 1E-4);
            Assert.AreEqual(80, emp.InverseCDF(emp.CDF(80)), 1E-4);

        }

        /// <summary>
        /// This test compares against Palisades '@Risk' empirical CDF function. 
        /// This test is described in the RMC-BestFit verification report. 
        /// </summary>
        [TestMethod]
        public void Test_Empirical_PalisadesAtRisk()
        {
            // nonparametric distribution for USGS 01562000 from Bulletin 17C test sites
            var xValues = new double[] { 3180, 4340, 4670, 4720, 5020, 6180, 6270, 7410, 7800, 8130, 8320, 8400, 8450, 8640, 8690, 8900, 8990, 9040, 9220, 9640, 9830, 10200, 10300, 10600, 10800, 10800, 11100, 11100, 11300, 11600, 11700, 11700, 11800, 11800, 12000, 12200, 12200, 12300, 12500, 12600, 12700, 12700, 12900, 13200, 13200, 13400, 13400, 13600, 13800, 14000, 14100, 14500, 14500, 14600, 15100, 15100, 15200, 15600, 16200, 17200, 17400, 17700, 17700, 17800, 18000, 18300, 18400, 18400, 18400, 18500, 18500, 18600, 18900, 19100, 19200, 19400, 19900, 20400, 20900, 21000, 21200, 21500, 21800, 22100, 22300, 22400, 22500, 22700, 22800, 23600, 26800, 29000, 31300, 39200, 40200, 42900, 45800, 71300, 80500 };
            var pValues = new double[] { 0.010036801605888, 0.020073603211777, 0.030110404817665, 0.040147206423553, 0.050184008029441, 0.0602208096353291, 0.070257611241218, 0.080294412847106, 0.090331214452994, 0.100368016058882, 0.110404817664771, 0.120441619270659, 0.130478420876547, 0.140515222482436, 0.150552024088324, 0.160588825694212, 0.1706256273001, 0.180662428905989, 0.190699230511877, 0.200736032117765, 0.210772833723653, 0.220809635329542, 0.23084643693543, 0.240883238541318, 0.250920040147206, 0.260956841753095, 0.270993643358983, 0.281030444964871, 0.291067246570759, 0.301104048176648, 0.311140849782536, 0.321177651388424, 0.331214452994312, 0.341251254600201, 0.351288056206089, 0.361324857811977, 0.371361659417865, 0.381398461023754, 0.391435262629642, 0.40147206423553, 0.411508865841419, 0.421545667447307, 0.431582469053195, 0.441619270659083, 0.451656072264972, 0.46169287387086, 0.471729675476748, 0.481766477082636, 0.491803278688525, 0.501840080294413, 0.511876881900301, 0.521913683506189, 0.531950485112078, 0.541987286717966, 0.552024088323854, 0.562060889929742, 0.572097691535631, 0.582134493141519, 0.592171294747407, 0.602208096353295, 0.612244897959184, 0.622281699565072, 0.63231850117096, 0.642355302776848, 0.652392104382737, 0.662428905988625, 0.672465707594513, 0.682502509200401, 0.69253931080629, 0.702576112412178, 0.712612914018066, 0.722649715623955, 0.732686517229843, 0.742723318835731, 0.752760120441619, 0.762796922047508, 0.772833723653396, 0.782870525259284, 0.792907326865172, 0.802944128471061, 0.812980930076949, 0.823017731682837, 0.833054533288725, 0.843091334894614, 0.853128136500502, 0.86316493810639, 0.873201739712278, 0.883238541318167, 0.893275342924055, 0.903312144529943, 0.913348946135831, 0.92338574774172, 0.933422549347608, 0.943459350953496, 0.953496152559384, 0.963532954165273, 0.973569755771161, 0.989071038251366, 0.994535519125683 };
            
            // @Risk does not do any interpolation transforms
            var emp = new EmpiricalDistribution(xValues, pValues) { ProbabilityTransform = Transform.None };

            // Test moments
            // The method included in Numerics is more accurate than the method in @Risk
            // Compare at 10% relative difference
            Assert.AreEqual(16763.82, emp.Mean, 1E-1 * 16763.82);
            Assert.AreEqual(11405.12, emp.StandardDeviation, 1E-1 * 11405.12);
            Assert.AreEqual(3.0303, emp.Skewness, 1E-1 * 3.0303);
            Assert.AreEqual(15.4169, emp.Kurtosis, 1E-1 * 15.4169);

            // Test percentiles
            Assert.AreEqual(5014.50, emp.InverseCDF(0.05), 1E-2);
            Assert.AreEqual(10781.67, emp.InverseCDF(0.25), 1E-2);
            Assert.AreEqual(13963.33, emp.InverseCDF(0.50), 1E-2);
            Assert.AreEqual(19172.50, emp.InverseCDF(0.75), 1E-2);
            Assert.AreEqual(39851.67, emp.InverseCDF(0.95), 1E-2);

        }

        /// <summary>
        /// Test: Convolve two uniform distributions.
        /// The convolution of two uniform distributions should approximate a triangular distribution.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveTwoUniformDistributions()
        {
            // Create two uniform distributions U(0,1)
            // Using sample data to create empirical distributions
            var sample1 = new double[100];
            var sample2 = new double[100];

            for (int i = 0; i < 100; i++)
            {
                sample1[i] = i / 99.0; // Values from 0 to 1
                sample2[i] = i / 99.0;
            }

            var dist1 = new EmpiricalDistribution(sample1, PlottingPositions.PlottingPostionType.Weibull);
            var dist2 = new EmpiricalDistribution(sample2, PlottingPositions.PlottingPostionType.Weibull);

            // Convolve with 1000 points
            var convolved = EmpiricalDistribution.Convolve(dist1, dist2, 1000);

            // Assert number of points
            Assert.HasCount(1000, convolved.XValues, "Should have exactly 1000 points");

            // Expected: Min ≈ 0, Max ≈ 2, Mean ≈ 1
            Assert.AreEqual(0.0, convolved.Minimum, 0.05, "Minimum should be approximately 0");
            Assert.AreEqual(2.0, convolved.Maximum, 0.05, "Maximum should be approximately 2");
            Assert.AreEqual(1.0, convolved.Mean, 0.1, "Mean should be approximately 1");

            // Median should be near 1 for symmetric distribution
            Assert.AreEqual(1.0, convolved.Median, 0.1, "Median should be approximately 1");

            // Verify CDF properties
            Assert.IsTrue(convolved.CDF(convolved.Minimum) >= 0.0 && convolved.CDF(convolved.Minimum) <= 0.01,
                "CDF at minimum should be close to 0");
            Assert.IsTrue(convolved.CDF(convolved.Maximum) >= 0.99 && convolved.CDF(convolved.Maximum) <= 1.0,
                "CDF at maximum should be close to 1");
            Assert.AreEqual(0.5, convolved.CDF(convolved.Mean), 0.1, "CDF at mean should be approximately 0.5");
        }

        /// <summary>
        /// Test: Convolve two normal-like distributions.
        /// Creates empirical distributions from normal samples.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveTwoNormalDistributions()
        {
            // Create sample data that approximates normal distributions
            // N(0, 1) and N(0, 1)
            var sample1 = new List<double>();
            var sample2 = new List<double>();

            // Generate more points for better accuracy (500 points)
            for (int i = 1; i < 500; i++)
            {
                double p = i / 500.0;
                double z = Normal.StandardZ(p);
                sample1.Add(z);
                sample2.Add(z);
            }

            var dist1 = new EmpiricalDistribution(sample1, PlottingPositions.PlottingPostionType.Weibull);
            var dist2 = new EmpiricalDistribution(sample2, PlottingPositions.PlottingPostionType.Weibull);

            // Convolve with 2048 points
            var convolved = EmpiricalDistribution.Convolve(dist1, dist2, 2048);

            // Assert number of points
            Assert.HasCount(2048, convolved.XValues, "Should have exactly 2048 points");

            // Expected: For N(0,1) + N(0,1) = N(0, sqrt(2)) 
            // Mean ≈ 0, StdDev ≈ 1.414
            // Use more generous tolerances due to empirical approximation
            Assert.AreEqual(0.0, convolved.Mean, 0.3, "Mean should be approximately 0");
            Assert.AreEqual(1.414, convolved.StandardDeviation, 0.15, "StdDev should be approximately 1.414");
        }

        /// <summary>
        /// Test: Convolve two distributions with different ranges.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveDifferentRanges()
        {
            // Distribution 1: Uniform on [0, 10]
            var sample1 = new double[50];
            for (int i = 0; i < 50; i++)
            {
                sample1[i] = 10.0 * i / 49.0;
            }

            // Distribution 2: Uniform on [5, 15]
            var sample2 = new double[50];
            for (int i = 0; i < 50; i++)
            {
                sample2[i] = 5.0 + 10.0 * i / 49.0;
            }

            var dist1 = new EmpiricalDistribution(sample1);
            var dist2 = new EmpiricalDistribution(sample2);

            // Convolve with 500 points
            var convolved = EmpiricalDistribution.Convolve(dist1, dist2, 500);

            // Assert number of points
            Assert.HasCount(500, convolved.XValues, "Should have exactly 500 points");

            // Expected: Range ≈ [5, 25], Mean ≈ 15
            Assert.AreEqual(5.0, convolved.Minimum, 0.5, "Minimum should be approximately 5");
            Assert.AreEqual(25.0, convolved.Maximum, 0.5, "Maximum should be approximately 25");
            Assert.AreEqual(15.0, convolved.Mean, 0.5, "Mean should be approximately 15");
        }

        /// <summary>
        /// Test: Convolve five identical distributions.
        /// Tests the list convolution method with multiple distributions.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveFiveDistributions()
        {
            // Create 5 identical uniform distributions U(0, 2)
            var distributions = new List<EmpiricalDistribution>();

            for (int j = 0; j < 5; j++)
            {
                var sample = new double[100];
                for (int i = 0; i < 100; i++)
                {
                    sample[i] = 2.0 * i / 99.0; // Values from 0 to 2
                }
                distributions.Add(new EmpiricalDistribution(sample));
            }

            // Convolve all five with 1024 points
            var convolved = EmpiricalDistribution.Convolve(distributions, 1024);

            // Assert number of points
            Assert.HasCount(1024, convolved.XValues);

            // Expected: Range ≈ [0, 10], Mean ≈ 5
            Assert.AreEqual(0.0, convolved.Minimum, 0.2, "Minimum should be approximately 0");
            Assert.AreEqual(10.0, convolved.Maximum, 0.2, "Maximum should be approximately 10");
            Assert.AreEqual(5.0, convolved.Mean, 0.3, "Mean should be approximately 5");

            // Should be roughly symmetric, so median near mean
            Assert.AreEqual(convolved.Mean, convolved.Median, 0.5, "Median should be close to mean for symmetric distribution");

            // Verify CDF properties
            Assert.IsLessThanOrEqualTo(0.01,convolved.CDF(convolved.Minimum), "CDF at minimum should be close to 0");
            Assert.IsGreaterThanOrEqualTo(0.99,convolved.CDF(convolved.Maximum), "CDF at maximum should be close to 1");
        }

        /// <summary>
        /// Test: Convolve five different distributions.
        /// </summary>
        [TestMethod]
        public void Test_ConvolveFiveDifferentDistributions()
        {
            var distributions = new List<EmpiricalDistribution>();

            // Distribution 1: U(0, 1)
            var sample1 = Enumerable.Range(0, 50).Select(i => i / 49.0).ToArray();
            distributions.Add(new EmpiricalDistribution(sample1));

            // Distribution 2: U(1, 3)
            var sample2 = Enumerable.Range(0, 50).Select(i => 1.0 + 2.0 * i / 49.0).ToArray();
            distributions.Add(new EmpiricalDistribution(sample2));

            // Distribution 3: U(0, 2)
            var sample3 = Enumerable.Range(0, 50).Select(i => 2.0 * i / 49.0).ToArray();
            distributions.Add(new EmpiricalDistribution(sample3));

            // Distribution 4: U(0.5, 1.5)
            var sample4 = Enumerable.Range(0, 50).Select(i => 0.5 + 1.0 * i / 49.0).ToArray();
            distributions.Add(new EmpiricalDistribution(sample4));

            // Distribution 5: U(2, 4)
            var sample5 = Enumerable.Range(0, 50).Select(i => 2.0 + 2.0 * i / 49.0).ToArray();
            distributions.Add(new EmpiricalDistribution(sample5));

            // Calculate expected sum properties
            double expectedMean = distributions.Sum(d => d.Mean);
            double expectedVariance = distributions.Sum(d => Math.Pow(d.StandardDeviation, 2));
            double expectedStdDev = Math.Sqrt(expectedVariance);

            // Convolve all five
            var convolved = EmpiricalDistribution.Convolve(distributions, 1000);

            // Assert number of points
            Assert.HasCount(1000, convolved.XValues, "Should have exactly 1000 points");

            // Compare with expected (allow for numerical error)
            double meanError = Math.Abs(convolved.Mean - expectedMean) / expectedMean;
            double stdDevError = Math.Abs(convolved.StandardDeviation - expectedStdDev) / expectedStdDev;

            Assert.IsLessThan(0.05,meanError, $"Mean error {meanError:P2} should be less than 5%");
            Assert.IsLessThan(0.15,stdDevError, $"StdDev error {stdDevError:P2} should be less than 15%");

            // Verify range is reasonable
            double expectedMin = distributions.Sum(d => d.Minimum);
            double expectedMax = distributions.Sum(d => d.Maximum);
            Assert.AreEqual(expectedMin, convolved.Minimum, 0.5, "Minimum should be close to sum of minimums");
            Assert.AreEqual(expectedMax, convolved.Maximum, 0.5, "Maximum should be close to sum of maximums");
        }

        /// <summary>
        /// Test: Test with non-power-of-2 point counts.
        /// </summary>
        [TestMethod]
        public void Test_NonPowerOfTwoPoints()
        {
            var sample1 = Enumerable.Range(0, 50).Select(i => i / 49.0).ToArray();
            var sample2 = Enumerable.Range(0, 50).Select(i => i / 49.0).ToArray();

            var dist1 = new EmpiricalDistribution(sample1);
            var dist2 = new EmpiricalDistribution(sample2);

            // Test various non-power-of-2 sizes
            int[] testSizes = { 100, 500, 1000, 1500, 3000 };

            foreach (var size in testSizes)
            {
                var convolved = EmpiricalDistribution.Convolve(dist1, dist2, size);

                // Assert correct number of points
                Assert.HasCount(size, convolved.XValues, $"Should have exactly {size} points");

                // Assert reasonable properties
                Assert.AreEqual(1.0, convolved.Mean, 0.1, $"Mean should be approximately 1 for size {size}");
                Assert.IsTrue(convolved.Minimum >= 0 && convolved.Minimum <= 0.1, $"Minimum should be near 0 for size {size}");
                Assert.IsTrue(convolved.Maximum >= 1.9 && convolved.Maximum <= 2.1, $"Maximum should be near 2 for size {size}");
            }
        }

        /// <summary>
        /// Duplicate sample values are valid empirical ordinates for bootstrap samples.
        /// </summary>
        [TestMethod]
        public void Test_DuplicateSampleValues_AreValidForCdfAndInverseCdf()
        {
            var sample = new[] { 100d, 100d, 125d, 150d, 150d, 200d };
            var distribution = new EmpiricalDistribution(sample);

            double cdf = distribution.CDF(150d);
            double quantile = distribution.InverseCDF(0.5d);

            Assert.IsFalse(double.IsNaN(cdf) || double.IsInfinity(cdf));
            Assert.IsTrue(cdf >= 0d && cdf <= 1d);
            Assert.IsFalse(double.IsNaN(quantile) || double.IsInfinity(quantile));
            Assert.IsTrue(quantile >= distribution.Minimum && quantile <= distribution.Maximum);
        }

        /// <summary>
        /// Duplicate x-values are valid when the empirical ordinate set is explicitly non-strict on x.
        /// </summary>
        [TestMethod]
        public void Test_NonStrictDuplicateXValues_AreValidForCdfAndInverseCdf()
        {
            var xValues = new[] { 100d, 100d, 125d, 150d, 200d };
            var pValues = new[] { 0.1d, 0.2d, 0.45d, 0.7d, 0.95d };
            var distribution = new EmpiricalDistribution(xValues, pValues);

            double cdf = distribution.CDF(125d);
            double quantile = distribution.InverseCDF(0.45d);

            Assert.IsFalse(double.IsNaN(cdf) || double.IsInfinity(cdf));
            Assert.IsTrue(cdf >= 0d && cdf <= 1d);
            Assert.AreEqual(125d, quantile, 1E-12);
        }

        /// <summary>
        /// Duplicate x-values are rejected when the empirical ordinate set requires strict x ordering.
        /// </summary>
        [TestMethod]
        public void Test_StrictDuplicateXValues_ThrowWhenUsed()
        {
            var orderedPairedData = new OrderedPairedData(
                new[] { 100d, 100d, 125d },
                new[] { 0.1d, 0.2d, 0.9d },
                true,
                SortOrder.Ascending,
                true,
                SortOrder.Ascending);
            var distribution = new EmpiricalDistribution(orderedPairedData);

            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.CDF(110d));
        }

        /// <summary>
        /// Decreasing x-values remain invalid for empirical distributions configured as ascending.
        /// </summary>
        [TestMethod]
        public void Test_DecreasingXValues_ThrowWhenUsed()
        {
            var distribution = new EmpiricalDistribution(
                new[] { 100d, 150d, 125d },
                new[] { 0.1d, 0.5d, 0.9d });

            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.CDF(125d));
        }

        /// <summary>
        /// Empirical X values require an explicit ascending order even when the values happen to be ordered.
        /// </summary>
        [TestMethod]
        public void Test_NoneXOrder_ThrowsForCdfAndInverseCdf()
        {
            var distribution = new EmpiricalDistribution(
                new[] { 100d, 125d, 150d },
                new[] { 0.1d, 0.5d, 0.9d },
                SortOrder.None,
                SortOrder.Ascending);

            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.CDF(125d));
            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.InverseCDF(0.5d));
        }

        /// <summary>
        /// Descending empirical X values are rejected because empirical interpolation requires ascending X values.
        /// </summary>
        [TestMethod]
        public void Test_DescendingXOrder_ThrowsForCdfAndInverseCdf()
        {
            var distribution = new EmpiricalDistribution(
                new[] { 200d, 150d, 100d },
                new[] { 0.9d, 0.5d, 0.1d },
                SortOrder.Descending,
                SortOrder.Descending);
            var orderedPairedData = new OrderedPairedData(
                new[] { 200d, 150d, 100d },
                new[] { 0.9d, 0.5d, 0.1d },
                true,
                SortOrder.Descending,
                true,
                SortOrder.Descending);

            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.CDF(150d));
            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.InverseCDF(0.5d));
            Assert.Throws<ArgumentOutOfRangeException>(() => new EmpiricalDistribution(orderedPairedData));
        }

        /// <summary>
        /// Descending probabilities remain valid when ascending X values define a survival curve.
        /// </summary>
        [TestMethod]
        public void Test_DescendingProbabilityOrder_SupportsCdfAndInverseCdf()
        {
            var distribution = new EmpiricalDistribution(
                new[] { 100d, 150d, 200d },
                new[] { 0.9d, 0.5d, 0.1d },
                SortOrder.Ascending,
                SortOrder.Descending);

            Assert.AreEqual(0.5d, distribution.CDF(150d), 1E-12);
            Assert.AreEqual(150d, distribution.InverseCDF(0.5d), 1E-12);
        }

        /// <summary>
        /// Non-finite probabilities remain invalid.
        /// </summary>
        [TestMethod]
        public void Test_NonFiniteProbability_ThrowsWhenUsed()
        {
            var distribution = new EmpiricalDistribution(
                new[] { 100d, 125d, 150d },
                new[] { 0.1d, double.NaN, 0.9d });

            Assert.Throws<ArgumentOutOfRangeException>(() => distribution.CDF(125d));
        }

        /// <summary>
        /// The extrapolation default (None) must reproduce the endpoint hold exactly,
        /// including the far-tail inverse guards, so a model without extrapolation is bit-identical. This
        /// is the bit-identity gate for the distribution's extrapolation property.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_DefaultNone_MatchesEndpointHold()
        {
            var dist = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d });
            Assert.AreEqual(ExtrapolationSides.None, dist.Extrapolation);

            Assert.AreEqual(0.1d, dist.CDF(50d));
            Assert.AreEqual(0.9d, dist.CDF(400d));
            Assert.AreEqual(100d, dist.InverseCDF(0.01d));
            Assert.AreEqual(300d, dist.InverseCDF(0.99d));
            Assert.AreEqual(100d, dist.InverseCDF(0d));
            Assert.AreEqual(300d, dist.InverseCDF(1d));
            Assert.AreEqual(100d, dist.InverseCDF(1E-17));
            Assert.AreEqual(300d, dist.InverseCDF(1d - 1E-17));
            Assert.AreEqual(0d, dist.PDF(50d));
            Assert.AreEqual(0d, dist.PDF(400d));
        }

        /// <summary>
        /// Sided extension on the default normal-Z probability axis, hand-computed from the
        /// boundary segments in z-space, both directions, including the one-sided holds. The
        /// extrapolation sides are defined in X space.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_SidedExtension_NormalZ()
        {
            var dist = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d })
            {
                Extrapolation = ExtrapolationSides.Both
            };
            double z01 = Normal.StandardZ(0.1);
            double z05 = Normal.StandardZ(0.5);
            double z09 = Normal.StandardZ(0.9);

            // CDF extends the boundary segments linearly in z and back-transforms into (0, 1).
            double validLow = Normal.StandardCDF(z01 + (z05 - z01) / 100d * (50d - 100d));
            double validHigh = Normal.StandardCDF(z09 + (z09 - z05) / 100d * (400d - 300d));
            Assert.AreEqual(validLow, dist.CDF(50d), 1E-12);
            Assert.AreEqual(validHigh, dist.CDF(400d), 1E-12);
            Assert.IsTrue(dist.CDF(50d) > 0d && dist.CDF(400d) < 1d);

            // InverseCDF extends the same segments beyond the table span.
            double validXLow = 100d + (200d - 100d) / (z05 - z01) * (Normal.StandardZ(0.01) - z01);
            double validXHigh = 200d + (300d - 200d) / (z09 - z05) * (Normal.StandardZ(0.99) - z05);
            Assert.AreEqual(validXLow, dist.InverseCDF(0.01d), 1E-9);
            Assert.AreEqual(validXHigh, dist.InverseCDF(0.99d), 1E-9);

            // The extended pair stays consistent.
            Assert.AreEqual(0.01d, dist.CDF(dist.InverseCDF(0.01d)), 1E-10);
            Assert.AreEqual(0.99d, dist.CDF(dist.InverseCDF(0.99d)), 1E-10);

            // One-sided policies hold the other side, in X space.
            var below = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d })
            {
                Extrapolation = ExtrapolationSides.Below
            };
            Assert.AreEqual(validLow, below.CDF(50d), 1E-12);
            Assert.AreEqual(0.9d, below.CDF(400d));
            Assert.AreEqual(validXLow, below.InverseCDF(0.01d), 1E-9);
            Assert.AreEqual(300d, below.InverseCDF(0.99d));

            var above = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d })
            {
                Extrapolation = ExtrapolationSides.Above
            };
            Assert.AreEqual(0.1d, above.CDF(50d));
            Assert.AreEqual(validHigh, above.CDF(400d), 1E-12);
            Assert.AreEqual(100d, above.InverseCDF(0.01d));
            Assert.AreEqual(validXHigh, above.InverseCDF(0.99d), 1E-9);
        }

        /// <summary>
        /// A probability axis under a None transform extends linearly and can leave the unit
        /// interval; the CDF must clamp its output to [0, 1] under extrapolation.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_CdfClamp_LinearProbabilityAxis()
        {
            var dist = new EmpiricalDistribution(new[] { 0d, 1d }, new[] { 0.2d, 0.8d })
            {
                ProbabilityTransform = Transform.None,
                Extrapolation = ExtrapolationSides.Both
            };
            Assert.AreEqual(1d, dist.CDF(2d));
            Assert.AreEqual(0d, dist.CDF(-1d));
            Assert.AreEqual(0.5d, dist.CDF(0.5d), 1E-12);
        }

        /// <summary>
        /// A survival-form table (descending stored probabilities) must extend the same X sides
        /// as its ascending twin: the complemented CDF branch and the X-side to lookup-side
        /// mapping in the inverse are both exercised.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_DescendingOrientation_MatchesAscending()
        {
            var ascending = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d })
            {
                Extrapolation = ExtrapolationSides.Both
            };
            var survival = new EmpiricalDistribution(
                new[] { 100d, 200d, 300d },
                new[] { 0.9d, 0.5d, 0.1d },
                SortOrder.Ascending,
                SortOrder.Descending)
            {
                Extrapolation = ExtrapolationSides.Both
            };

            // Same distribution, either storage orientation: extended values agree.
            Assert.AreEqual(ascending.CDF(50d), survival.CDF(50d), 1E-12);
            Assert.AreEqual(ascending.CDF(400d), survival.CDF(400d), 1E-12);
            Assert.AreEqual(ascending.InverseCDF(0.01d), survival.InverseCDF(0.01d), 1E-8);
            Assert.AreEqual(ascending.InverseCDF(0.99d), survival.InverseCDF(0.99d), 1E-8);

            // One-sided X-space semantics on the survival form: Below extends only the low-X end.
            var survivalBelow = new EmpiricalDistribution(
                new[] { 100d, 200d, 300d },
                new[] { 0.9d, 0.5d, 0.1d },
                SortOrder.Ascending,
                SortOrder.Descending)
            {
                Extrapolation = ExtrapolationSides.Below
            };
            Assert.IsLessThan(100d, survivalBelow.InverseCDF(0.01d));
            Assert.AreEqual(300d, survivalBelow.InverseCDF(0.99d));
            Assert.IsLessThan(0.1d, survivalBelow.CDF(50d));
            Assert.AreEqual(0.9d, survivalBelow.CDF(400d));
        }

        /// <summary>
        /// The far-tail inverse guards stay total under extrapolation (probabilities at or
        /// beyond the 1e-16 floors evaluate the extended lookup at the floor instead of holding
        /// the endpoint), and the inverse stays monotone across both boundaries.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_InverseMonotoneAndGuards()
        {
            var dist = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d })
            {
                Extrapolation = ExtrapolationSides.Both
            };

            double x0 = dist.InverseCDF(0d);
            double x1 = dist.InverseCDF(1d);
            Assert.IsTrue(Tools.IsFinite(x0) && x0 < 100d);
            Assert.IsTrue(Tools.IsFinite(x1) && x1 > 300d);
            Assert.AreEqual(dist.InverseCDF(1E-17), x0, 0d);
            Assert.AreEqual(dist.InverseCDF(1d - 1E-17), x1, 0d);

            var grid = new[] { 0d, 1E-16, 1E-10, 1E-4, 0.01d, 0.1d, 0.3d, 0.5d, 0.7d, 0.9d, 0.99d, 1d - 1E-4, 1d - 1E-10, 1d - 1E-16, 1d };
            for (int i = 1; i < grid.Length; i++)
            {
                Assert.IsGreaterThanOrEqualTo(dist.InverseCDF(grid[i - 1]), dist.InverseCDF(grid[i]),
                    $"InverseCDF must be monotone across the extended tail at u = {grid[i]}.");
            }
        }

        /// <summary>
        /// The policy survives Clone, serializes only when non-default (conditional presence),
        /// restores by name, defaults when absent, and rejects an unparseable value.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Extrapolation_CloneAndSerialization()
        {
            var dist = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d });

            // Default: no attribute, and set-then-clear restores the exact serialized form.
            var defaultXml = dist.ToXElement();
            Assert.IsNull(defaultXml.Attribute(nameof(EmpiricalDistribution.Extrapolation)));
            string baseline = defaultXml.ToString();
            dist.Extrapolation = ExtrapolationSides.Both;
            dist.Extrapolation = ExtrapolationSides.None;
            Assert.AreEqual(baseline, dist.ToXElement().ToString());

            // Clone carries the policy.
            dist.Extrapolation = ExtrapolationSides.Both;
            var clone = (EmpiricalDistribution)dist.Clone();
            Assert.AreEqual(ExtrapolationSides.Both, clone.Extrapolation);

            // Non-default round-trips by enum name.
            var xml = dist.ToXElement();
            Assert.AreEqual(nameof(ExtrapolationSides.Both), xml.Attribute(nameof(EmpiricalDistribution.Extrapolation))?.Value);
            var restored = EmpiricalDistribution.FromXElement(xml);
            Assert.AreEqual(ExtrapolationSides.Both, restored.Extrapolation);
            Assert.AreEqual(dist.InverseCDF(0.01d), restored.InverseCDF(0.01d), 1E-12);

            // Absent attribute reads as the default; an invalid value throws (the strict pattern).
            Assert.AreEqual(ExtrapolationSides.None, EmpiricalDistribution.FromXElement(defaultXml).Extrapolation);
            var invalid = dist.ToXElement();
            invalid.SetAttributeValue(nameof(EmpiricalDistribution.Extrapolation), "Sideways");
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.FromXElement(invalid));
        }

        /// <summary>
        /// A survival table (strictly descending probabilities) keeps its orientation through
        /// Clone and XML: CDF = 1 − S and InverseCDF reads the exceedance axis on every copy.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_SurvivalTable_CloneAndSerialization()
        {
            var dist = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.5d, 0.1d, 0.01d }, SortOrder.Ascending, SortOrder.Descending);
            Assert.IsTrue(dist.ParametersValid);

            var clone = (EmpiricalDistribution)dist.Clone();
            var restored = EmpiricalDistribution.FromXElement(dist.ToXElement());
            foreach (var copy in new[] { clone, restored })
            {
                Assert.AreEqual(SortOrder.Descending, copy.ProbabilityOrder);
                Assert.IsTrue(copy.StrictProbabilities);
                Assert.IsFalse(copy.StrictXValues);
                Assert.IsTrue(copy.ParametersValid);
                Assert.AreEqual(1d - 0.1d, copy.CDF(200d), 0d);
                Assert.AreEqual(dist.CDF(250d), copy.CDF(250d), 0d);
                Assert.AreEqual(200d, copy.InverseCDF(1d - 0.1d), 1E-10);
                Assert.AreEqual(dist.InverseCDF(0.7d), copy.InverseCDF(0.7d), 0d);
            }
        }

        /// <summary>
        /// A tied-probability (non-strict) table and a strict-X table built from ordered paired
        /// data keep their strictness through Clone and XML. The XML writes strictness only when
        /// it departs from the default, so a default distribution serializes exactly as before.
        /// </summary>
        [TestMethod]
        public void Test_Empirical_Strictness_CloneAndSerialization()
        {
            // Default: no strictness attributes; the serialized form is unchanged.
            var standard = new EmpiricalDistribution(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d });
            var standardXml = standard.ToXElement();
            Assert.IsNull(standardXml.Attribute(nameof(EmpiricalDistribution.StrictXValues)));
            Assert.IsNull(standardXml.Attribute(nameof(EmpiricalDistribution.StrictProbabilities)));
            Assert.AreEqual(
                "<Distribution Type=\"Empirical\" XTransform=\"None\" ProbabilityTransform=\"NormalZ\" XValues=\"100|200|300\" ProbabilityValues=\"0.10000000000000001|0.5|0.90000000000000002\" ProbabilityOrder=\"Ascending\" />",
                standardXml.ToString());

            // Tied probabilities: valid only with non-strict probabilities.
            var tied = new EmpiricalDistribution(new OrderedPairedData(new[] { 100d, 200d, 300d, 400d }, new[] { 0.1d, 0.5d, 0.5d, 0.9d }, false, SortOrder.Ascending, false, SortOrder.Ascending));
            Assert.IsTrue(tied.ParametersValid);
            Assert.IsFalse(tied.StrictProbabilities);
            var tiedClone = (EmpiricalDistribution)tied.Clone();
            Assert.IsFalse(tiedClone.StrictProbabilities);
            Assert.IsTrue(tiedClone.ParametersValid);
            Assert.AreEqual(tied.CDF(350d), tiedClone.CDF(350d), 0d);
            var tiedXml = tied.ToXElement();
            Assert.AreEqual(bool.FalseString, tiedXml.Attribute(nameof(EmpiricalDistribution.StrictProbabilities))?.Value);
            Assert.IsNull(tiedXml.Attribute(nameof(EmpiricalDistribution.StrictXValues)));
            var tiedRestored = EmpiricalDistribution.FromXElement(tiedXml);
            Assert.IsFalse(tiedRestored.StrictProbabilities);
            Assert.IsTrue(tiedRestored.ParametersValid);
            Assert.AreEqual(tied.CDF(350d), tiedRestored.CDF(350d), 0d);
            Assert.AreEqual(tiedXml.ToString(), tiedRestored.ToXElement().ToString());

            // Strict X values round-trip as well.
            var strictX = new EmpiricalDistribution(new OrderedPairedData(new[] { 100d, 200d, 300d }, new[] { 0.1d, 0.5d, 0.9d }, true, SortOrder.Ascending, true, SortOrder.Ascending));
            Assert.IsTrue(strictX.StrictXValues);
            Assert.IsTrue(((EmpiricalDistribution)strictX.Clone()).StrictXValues);
            var strictXml = strictX.ToXElement();
            Assert.AreEqual(bool.TrueString, strictXml.Attribute(nameof(EmpiricalDistribution.StrictXValues))?.Value);
            Assert.IsTrue(EmpiricalDistribution.FromXElement(strictXml).StrictXValues);

            // An unparseable strictness value is rejected.
            tiedXml.SetAttributeValue(nameof(EmpiricalDistribution.StrictProbabilities), "sometimes");
            Assert.Throws<ArgumentException>(() => EmpiricalDistribution.FromXElement(tiedXml));
        }

        /// <summary>
        /// A not-a-number probability is rejected like any probability outside [0, 1], while the endpoints keep
        /// returning the table bounds.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_NaNProbabilityThrows()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.1d, 0.5d, 0.9d });
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => empirical.InverseCDF(double.NaN));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => empirical.InverseCDF(-0.1d));
            Assert.AreEqual(1d, empirical.InverseCDF(0d), 0d);
            Assert.AreEqual(3d, empirical.InverseCDF(1d), 0d);
            Assert.AreEqual(2d, empirical.InverseCDF(0.5d), 1E-12);
        }

        /// <summary>
        /// A table whose probabilities have no sort order cannot be inverted: the inverse CDF, median, and moments
        /// report an invalid operation carrying the interpolation failure, while the CDF still interpolates.
        /// </summary>
        [TestMethod]
        public void Test_UnorderedProbabilities_InverseThrowsInvalidOperation()
        {
            var empirical = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.1d, 0.9d, 0.4d }, SortOrder.Ascending, SortOrder.None);

            var inverse = Assert.ThrowsExactly<InvalidOperationException>(() => empirical.InverseCDF(0.5d));
            Assert.IsNotNull(inverse.InnerException);
            Assert.AreEqual(typeof(Exception), inverse.InnerException!.GetType());
            Assert.ThrowsExactly<InvalidOperationException>(() => empirical.Median);
            Assert.ThrowsExactly<InvalidOperationException>(() => empirical.Mean);
            Assert.AreEqual(0.9d, empirical.CDF(2d), 1E-12);
        }

        /// <summary>
        /// Minimum and Maximum equal the first and last entries of the public X-value collection bit for bit,
        /// including signed zeros and single-ordinate tables, and an empty table keeps the sequence exception.
        /// </summary>
        [TestMethod]
        public void Test_MinimumMaximum_MatchXValueEndpoints()
        {
            var tables = new[]
            {
                new EmpiricalDistribution(new[] { -0d, 2d, 4d, 8d }, new[] { 0.05d, 0.4d, 0.8d, 0.99d }),
                new EmpiricalDistribution(new[] { -3d, -1d, 0d }, new[] { 0.2d, 0.6d, 0.7d }),
                new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.9d, 0.5d, 0.1d }, SortOrder.Ascending, SortOrder.Descending),
                new EmpiricalDistribution(new[] { 5d }, new[] { 0.5d })
            };
            foreach (var table in tables)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(table.XValues.First()), BitConverter.DoubleToInt64Bits(table.Minimum));
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(table.XValues.Last()), BitConverter.DoubleToInt64Bits(table.Maximum));
            }
            var empty = new EmpiricalDistribution(new double[0], new double[0]);
            Assert.ThrowsExactly<InvalidOperationException>(() => empty.Minimum);
            Assert.ThrowsExactly<InvalidOperationException>(() => empty.Maximum);
        }

        /// <summary>
        /// The descending survival table used by the survival-quantile tests: x = 1..8 with exceedance
        /// probabilities from 0.5 down to 1E-18, interpolated linearly in normal-Z space.
        /// </summary>
        private static readonly double[] SurvivalX = { 1d, 2d, 3d, 4d, 5d, 6d, 7d, 8d };

        /// <summary>
        /// The exceedance probabilities of <see cref="SurvivalX"/>.
        /// </summary>
        private static readonly double[] SurvivalP = { 0.5, 0.1, 1E-2, 1E-4, 1E-6, 1E-9, 1E-13, 1E-18 };

        /// <summary>
        /// Survival probability at x by independent linear interpolation of the standard normal Z of the table.
        /// </summary>
        /// <param name="x">A value inside the table range.</param>
        /// <returns>The interpolated exceedance probability.</returns>
        private static double TableSurvival(double x)
        {
            int k = Math.Min((int)Math.Floor(x) - 1, SurvivalX.Length - 2);
            double z0 = Normal.StandardZ(SurvivalP[k]), z1 = Normal.StandardZ(SurvivalP[k + 1]);
            return Normal.StandardCDF(z0 + (x - SurvivalX[k]) * (z1 - z0));
        }

        /// <summary>
        /// A descending (survival) table is read at q directly: the survival quantile at q = 1E-8, 1E-12, 1E-15 and
        /// 1E-17 round-trips through the table to q within 1E-12 relative, while InverseCDF(1 - q) loses up to
        /// 8E-4 relative at 1E-15 and returns the table maximum at 1E-17, below its far-tail floor.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_SurvivalTableRoundTrip()
        {
            var emp = new EmpiricalDistribution(SurvivalX, SurvivalP, SortOrder.Ascending, SortOrder.Descending);
            foreach (double q in new[] { 1E-8, 1E-12, 1E-15, 1E-17 })
            {
                double x = emp.InverseCCDF(q);
                Assert.IsLessThan(8d, x);
                Assert.AreEqual(q, TableSurvival(x), 1E-12 * q, $"q {q}");
            }
            Assert.AreEqual(1E-15, TableSurvival(emp.InverseCDF(1d - 1E-15)), 1E-3 * 1E-15);
            Assert.AreNotEqual(1E-15, TableSurvival(emp.InverseCDF(1d - 1E-15)), 1E-6 * 1E-15);
            Assert.AreEqual(8d, emp.InverseCDF(1d - 1E-17), 0d);
            // Table knots are returned exactly, and the survival quantile decreases with q.
            Assert.AreEqual(6d, emp.InverseCCDF(1E-9), 1E-13);
            Assert.AreEqual(7d, emp.InverseCCDF(1E-13), 1E-13);
            Assert.IsGreaterThan(emp.InverseCCDF(1E-12), emp.InverseCCDF(1E-13));
        }

        /// <summary>
        /// Below the smallest stored exceedance probability the survival quantile holds the table maximum, or, with
        /// upper extrapolation, extends the last normal-Z segment; the result equals the independent extension.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_SurvivalTableTailAndExtrapolation()
        {
            var held = new EmpiricalDistribution(SurvivalX, SurvivalP, SortOrder.Ascending, SortOrder.Descending);
            Assert.AreEqual(8d, held.InverseCCDF(1E-20), 0d);
            var extended = new EmpiricalDistribution(SurvivalX, SurvivalP, SortOrder.Ascending, SortOrder.Descending)
            {
                Extrapolation = ExtrapolationSides.Above
            };
            double z7 = Normal.StandardZ(1E-13), z8 = Normal.StandardZ(1E-18);
            double expected = 7d + (Normal.StandardZ(1E-20) - z7) / (z8 - z7);
            Assert.AreEqual(expected, extended.InverseCCDF(1E-20), 1E-12);
            Assert.IsGreaterThan(8d, extended.InverseCCDF(1E-20));
        }

        /// <summary>
        /// Where the base form loses nothing the result is InverseCDF(1 - q) bit for bit: q = 0, q of at least 0.5,
        /// and ascending tables. Invalid probabilities and an unordered table throw.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_BaseCasesAndInvalidInputs()
        {
            var survival = new EmpiricalDistribution(SurvivalX, SurvivalP, SortOrder.Ascending, SortOrder.Descending);
            foreach (double q in new[] { 0d, 0.5, 0.75, 1d - 1E-17, 1d })
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(survival.InverseCDF(1d - q)), BitConverter.DoubleToInt64Bits(survival.InverseCCDF(q)));
            var ascending = new EmpiricalDistribution(SurvivalX, new[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999 }, SortOrder.Ascending, SortOrder.Ascending);
            foreach (double q in new[] { 0d, 1E-17, 1E-6, 0.2, 0.5, 0.9, 1d })
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(ascending.InverseCDF(1d - q)), BitConverter.DoubleToInt64Bits(ascending.InverseCCDF(q)));
            foreach (double bad in new[] { -0.1, 1.1, double.NaN })
                Assert.AreEqual("probability", Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => survival.InverseCCDF(bad)).ParamName);
            var unordered = new EmpiricalDistribution(new[] { 1d, 2d, 3d }, new[] { 0.1d, 0.9d, 0.4d }, SortOrder.Ascending, SortOrder.None);
            Assert.ThrowsExactly<InvalidOperationException>(() => unordered.InverseCCDF(0.3));
        }

        /// <summary>
        /// Mode is the midpoint of the highest bin of the density grid: for x = {0, 1, 2, 4} and probabilities
        /// {0, 0.2, 0.7, 1} the bin densities are 0.2, 0.5 and 0.15, so the mode is 1.5. Ties go to the first bin.
        /// </summary>
        [TestMethod]
        public void Test_ModeIsHighestBin()
        {
            Assert.AreEqual(1.5, new EmpiricalDistribution(new[] { 0d, 1d, 2d, 4d }, new[] { 0d, 0.2, 0.7, 1d }).Mode, 1E-15);
            Assert.AreEqual(0.5, new EmpiricalDistribution(new[] { 0d, 1d, 2d }, new[] { 0d, 0.5, 1d }).Mode, 1E-15);
            // A tie between the first two bins (densities 0.4, 0.4 and 0.1) goes to the first.
            Assert.AreEqual(0.5, new EmpiricalDistribution(new[] { 0d, 1d, 2d, 3d }, new[] { 0d, 0.4, 0.8, 0.9 }).Mode, 1E-15);
        }

        /// <summary>
        /// With the NormalZ probability transform a table that reaches p = 1 has a segment whose end is infinite in z
        /// (StandardZ(1) is +∞): the CDF at the node next to it is the node's probability rather than NaN, and inside the
        /// segment it is the limit 1. StandardZ(0) is finite, so a table that starts at p = 0 was already finite; its
        /// nodes are checked too. Both x transforms behave the same.
        /// </summary>
        [TestMethod]
        public void Test_NormalZNodeNextToCertainProbability()
        {
            foreach (var xTransform in new[] { Transform.None, Transform.Logarithmic })
            {
                var upper = new EmpiricalDistribution(new[] { 0d, 1d, 2d, 4d }, new[] { 0d, 0.2, 0.7, 1d }) { XTransform = xTransform, ProbabilityTransform = Transform.NormalZ };
                Assert.AreEqual(0.2, upper.CDF(1d), 1E-15);
                Assert.AreEqual(0.7, upper.CDF(2d), 1E-15);
                Assert.AreEqual(1d, upper.CDF(3d), 0d);
                Assert.AreEqual(1d, upper.CDF(4d), 0d);
                Assert.AreEqual(2d, upper.InverseCDF(0.7), 1E-12);

                var lower = new EmpiricalDistribution(new[] { 0d, 1d, 2d, 4d }, new[] { 0d, 0.3, 0.8, 0.9 }) { XTransform = xTransform, ProbabilityTransform = Transform.NormalZ };
                Assert.AreEqual(0d, lower.CDF(0d), 0d);
                Assert.AreEqual(0.3, lower.CDF(1d), 1E-15);
                Assert.AreEqual(0.8, lower.CDF(2d), 1E-15);
                Assert.AreEqual(1d, lower.InverseCDF(0.3), 1E-12);
            }
        }

        /// <summary>
        /// A point mass is the mode: x = {0, 1, 1, 2} with probabilities {0, 0.2, 0.7, 1} has an atom of 0.5 at 1, and a
        /// sample with ties has its atoms at the tied values, the largest at the most repeated one.
        /// </summary>
        [TestMethod]
        public void Test_ModeIsLargestAtom()
        {
            Assert.AreEqual(1d, new EmpiricalDistribution(new[] { 0d, 1d, 1d, 2d }, new[] { 0d, 0.2, 0.7, 1d }).Mode, 0d);
            Assert.AreEqual(2d, new EmpiricalDistribution(new[] { 1d, 2d, 2d, 2d, 3d, 4d, 4d }).Mode, 0d);
        }

        /// <summary>
        /// Under the NormalZ probability transform and the Logarithmic x transform the density varies inside a bin, and
        /// Mode still returns the midpoint of the bin with the highest average density, |Δp|/|Δx|: for x = {1, 2, 3, 5}
        /// and probabilities {0.1, 0.3, 0.8, 0.9} that is the bin [2, 3], so 2.5.
        /// </summary>
        [TestMethod]
        public void Test_ModeIsBinMidpointUnderTransforms()
        {
            foreach (var x in new[] { Transform.None, Transform.Logarithmic })
            {
                foreach (var p in new[] { Transform.None, Transform.NormalZ, Transform.Logarithmic })
                {
                    var distribution = new EmpiricalDistribution(new[] { 1d, 2d, 3d, 5d }, new[] { 0.1, 0.3, 0.8, 0.9 }) { XTransform = x, ProbabilityTransform = p };
                    Assert.AreEqual(2.5, distribution.Mode, 0d, $"{x}, {p}");
                }
            }
        }

        /// <summary>
        /// Long tables take a vectorized scan of the bins; it returns the same bin as a direct scan. The expected values
        /// are exact by construction, with integer ordinates and probabilities that are multiples of 1/512: on an
        /// equally spaced table every bin has the same density and the first bin wins; a single heavier bin anywhere
        /// wins, including the last bins, which fall past the last full vector; equal heavier bins go to the first; and
        /// a repeated x value anywhere makes the largest atom the mode.
        /// </summary>
        [TestMethod]
        public void Test_ModeOnLongTables()
        {
            const int n = 259;
            double[] Grid() { var x = new double[n]; for (int i = 0; i < n; i++) x[i] = i; return x; }
            double[] Linear() { var p = new double[n]; for (int i = 0; i < n; i++) p[i] = i / 512d; return p; }

            // The ordinates and probabilities i / 512 are exact, so every bin has density 1/512 exactly: the first bin.
            Assert.AreEqual(0.5, new EmpiricalDistribution(Grid(), Linear()).Mode, 0d);

            // One bin twice as wide holds the same mass: every other bin is denser, still the first bin wins.
            var wide = Grid();
            for (int i = 101; i < n; i++) wide[i] += 1d;
            Assert.AreEqual(0.5, new EmpiricalDistribution(wide, Linear()).Mode, 0d);

            foreach (int heavy in new[] { 0, 1, 2, 3, 4, 5, 7, 8, 9, 57, 120, n - 4, n - 3, n - 2 })
            {
                // Narrowing bin [heavy, heavy + 1] to half its width doubles its density.
                var x = Grid();
                for (int i = heavy + 1; i < n; i++) x[i] -= 0.5;
                Assert.AreEqual(heavy + 0.25, new EmpiricalDistribution(x, Linear()).Mode, 0d, $"heavy bin {heavy}");
            }

            // Two equally narrowed bins: the first one wins.
            var twin = Grid();
            for (int i = 31; i < n; i++) twin[i] -= 0.5;
            for (int i = 150; i < n; i++) twin[i] -= 0.5;
            Assert.AreEqual(30.25, new EmpiricalDistribution(twin, Linear()).Mode, 0d);

            foreach (int tie in new[] { 1, 64, n - 1 })
            {
                // A repeated x value is a point mass, the mode wherever it sits.
                var x = Grid();
                x[tie] = x[tie - 1];
                Assert.AreEqual(x[tie], new EmpiricalDistribution(x, Linear()).Mode, 0d, $"tie at {tie}");
            }
        }
    }
}
