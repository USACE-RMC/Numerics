/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
            Assert.AreEqual(1000, convolved.XValues.Count, "Should have exactly 1000 points");

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
            Assert.AreEqual(2048, convolved.XValues.Count, "Should have exactly 2048 points");

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
            Assert.AreEqual(500, convolved.XValues.Count, "Should have exactly 500 points");

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
            Assert.AreEqual(1024, convolved.XValues.Count, "Should have exactly 1024 points");

            // Expected: Range ≈ [0, 10], Mean ≈ 5
            Assert.AreEqual(0.0, convolved.Minimum, 0.2, "Minimum should be approximately 0");
            Assert.AreEqual(10.0, convolved.Maximum, 0.2, "Maximum should be approximately 10");
            Assert.AreEqual(5.0, convolved.Mean, 0.3, "Mean should be approximately 5");

            // Should be roughly symmetric, so median near mean
            Assert.AreEqual(convolved.Mean, convolved.Median, 0.5, "Median should be close to mean for symmetric distribution");

            // Verify CDF properties
            Assert.IsTrue(convolved.CDF(convolved.Minimum) <= 0.01, "CDF at minimum should be close to 0");
            Assert.IsTrue(convolved.CDF(convolved.Maximum) >= 0.99, "CDF at maximum should be close to 1");
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
            Assert.AreEqual(1000, convolved.XValues.Count, "Should have exactly 1000 points");

            // Compare with expected (allow for numerical error)
            double meanError = Math.Abs(convolved.Mean - expectedMean) / expectedMean;
            double stdDevError = Math.Abs(convolved.StandardDeviation - expectedStdDev) / expectedStdDev;

            Assert.IsTrue(meanError < 0.05, $"Mean error {meanError:P2} should be less than 5%");
            Assert.IsTrue(stdDevError < 0.15, $"StdDev error {stdDevError:P2} should be less than 15%");

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
                Assert.AreEqual(size, convolved.XValues.Count, $"Should have exactly {size} points");

                // Assert reasonable properties
                Assert.AreEqual(1.0, convolved.Mean, 0.1, $"Mean should be approximately 1 for size {size}");
                Assert.IsTrue(convolved.Minimum >= 0 && convolved.Minimum <= 0.1, $"Minimum should be near 0 for size {size}");
                Assert.IsTrue(convolved.Maximum >= 1.9 && convolved.Maximum <= 2.1, $"Maximum should be near 2 for size {size}");
            }
        }

    }
}
