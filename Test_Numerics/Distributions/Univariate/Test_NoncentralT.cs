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
using Numerics.Distributions;
using System;
using System.Diagnostics;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Noncentral T distribution algorithm.
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
    public class Test_NoncentralT
    {
        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_PDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double pdf = t.PDF(1.4d);
            double result = 0.23552141805184526d;
            Assert.AreEqual(pdf, result, 1E-6);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_CDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double result = 0.15955740661144721d;
            Assert.AreEqual(cdf, result, 1E-6);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_InverseCDF()
        {
            var t = new NoncentralT(4d, 2.42d);
            double cdf = t.CDF(1.4d);
            double invcdf = t.InverseCDF(cdf);
            double result = 1.4d;
            Assert.AreEqual(invcdf, result, 1E-6);
            var table = new[,] { { 3.0d, 0.0d, 1d, 0.89758361765043326d }, { 3.0d, 0.0d, 2d, 0.9522670169d }, { 3.0d, 0.0d, 3d, 0.97116555718878128d }, { 3.0d, 0.5d, 1d, 0.8231218864d }, { 3.0d, 0.5d, 2d, 0.904902151d }, { 3.0d, 0.5d, 3d, 0.9363471834d }, { 3.0d, 1.0d, 1d, 0.7301025986d }, { 3.0d, 1.0d, 2d, 0.8335594263d }, { 3.0d, 1.0d, 3d, 0.8774010255d }, { 3.0d, 2.0d, 1d, 0.5248571617d }, { 3.0d, 2.0d, 2d, 0.6293856597d }, { 3.0d, 2.0d, 3d, 0.6800271741d }, { 3.0d, 4.0d, 1d, 0.20590131975d }, { 3.0d, 4.0d, 2d, 0.2112148916d }, { 3.0d, 4.0d, 3d, 0.2074730718d }, { 15.0d, 7.0d, 15d, 0.9981130072d }, { 15.0d, 7.0d, 20d, 0.999487385d }, { 15.0d, 7.0d, 25d, 0.9998391562d }, { 0.05d, 1.0d, 1d, 0.168610566972d }, { 0.05d, 1.0d, 2d, 0.16967950985d }, { 0.05d, 1.0d, 3d, 0.1701041003d }, { 4.0d, 2.0d, 10d, 0.9247683363d }, { 4.0d, 3.0d, 10d, 0.7483139269d }, { 4.0d, 4.0d, 10d, 0.4659802096d }, { 5.0d, 2.0d, 10d, 0.9761872541d }, { 5.0d, 3.0d, 10d, 0.8979689357d }, { 5.0d, 4.0d, 10d, 0.7181904627d }, { 6.0d, 2.0d, 10d, 0.9923658945d }, { 6.0d, 3.0d, 10d, 0.9610341649d }, { 6.0d, 4.0d, 10d, 0.868800735d } };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double x = table[i, 0];
                double delta = table[i, 1];
                double degF = table[i, 2];
                var target = new NoncentralT(degF, delta);
                double expected = table[i, 3];
                double actual = target.CDF(x);
                Assert.AreEqual(expected, actual, 0.000001d);
                double expectedX = target.InverseCDF(actual);
                Assert.AreEqual(expectedX, x, 0.000001d);
            }

        }

        /// <summary>
        /// Verifying input parameters can create Noncentral T distribution.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var t = new NoncentralT();
            Assert.AreEqual(10,t.DegreesOfFreedom);
            Assert.AreEqual(0, t.Noncentrality);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(1, t2.DegreesOfFreedom);
            Assert.AreEqual(1, t2.Noncentrality);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var t = new NoncentralT(0, 1);
            Assert.IsFalse(t.ParametersValid);

            var t2 = new NoncentralT(1,double.PositiveInfinity);
            Assert.IsFalse(t2.ParametersValid);

            var t3 = new NoncentralT(1,double.NaN);
            Assert.IsFalse(t3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var t = new NoncentralT();
            Assert.AreEqual("Degrees of Freedom (ν)", t.ParametersToString[0, 0]);
            Assert.AreEqual("Noncentrality (μ)", t.ParametersToString[1, 0]);
            Assert.AreEqual("10", t.ParametersToString[0, 1]);
            Assert.AreEqual("0", t.ParametersToString[1, 1]);
        }

        /// <summary>
        /// Testing mean.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Mean);

            var t2 = new NoncentralT(0, 1);
            Assert.AreEqual(double.NaN,t2.Mean);
        }

        /// <summary>
        /// Testing median.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Median, 1e-04);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(1.3202, t2.Median, 1e-04);
        }

        /// <summary>
        /// Testing mode.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0, t.Mode,  1E-4);

            var t3 = new NoncentralT(10, 1);
            Assert.AreEqual(0.9329, t3.Mode,  1e-04);
        }

        /// <summary>
        /// Testing standard deviation.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var t = new NoncentralT();
            Assert.AreEqual(1.1180, t.StandardDeviation,1e-04);

            var t2 = new NoncentralT(1, 0);
            Assert.AreEqual(double.NaN,t2.StandardDeviation);
        }

        /// <summary>
        /// Testing skew.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.0, t.Skewness,  1E-4);
        }

        /// <summary>
        /// Testing Kurtosis
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var t = new NoncentralT();
            Assert.AreEqual(4.0, t.Kurtosis, 1E-4);
        }

        /// <summary>
        /// Testing min and max functions.
        /// </summary>
        [TestMethod()]
        public void Test_MinMax()
        {
            var t = new NoncentralT();
            Assert.AreEqual(double.NegativeInfinity,t.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t.Maximum);

            var t2 = new NoncentralT(1, 1);
            Assert.AreEqual(double.NegativeInfinity, t2.Minimum);
            Assert.AreEqual(double.PositiveInfinity, t2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.38910, t.PDF(0), 1e-04);
            Assert.AreEqual(0.23036, t.PDF(1),1e-04);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(0.82955, t.CDF(1), 1e-04);
        }

        /// <summary>
        /// Testing inverse CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var t = new NoncentralT();
            Assert.AreEqual(double.NegativeInfinity,t.InverseCDF(0) );
            Assert.AreEqual(double.PositiveInfinity, t.InverseCDF(1));
            Assert.AreEqual(-0.26018, t.InverseCDF(0.4), 1e-04);
        }

        /// <summary>
        /// Validates NoncentralT InverseCDF against Stedinger (1983) Table 1.
        /// Table 1: Percentage points of ζ(0.90)-distribution for the 10-year event.
        /// Relationship: ζ_α(p) = NCT.InverseCDF(α; df=n-1, δ=z_p·√n) / √n
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_StedsingerTable1()
        {
            double zp = Normal.StandardZ(0.90);

            // Exact values from Stedinger (1983) Table 1 (p=0.90)
            // { sampleSize, alpha, expected zeta }
            var table = new[,]
            {
                { 10d, 0.005, 0.436 }, { 10d, 0.050, 0.712 }, { 10d, 0.250, 1.043 },
                { 10d, 0.750, 1.671 }, { 10d, 0.950, 2.355 }, { 10d, 0.995, 3.368 },
                { 20d, 0.005, 0.651 }, { 20d, 0.050, 0.858 }, { 20d, 0.250, 1.104 },
                { 20d, 0.750, 1.528 }, { 20d, 0.950, 1.926 }, { 20d, 0.995, 2.423 },
                { 50d, 0.005, 0.857 }, { 50d, 0.050, 1.000 }, { 50d, 0.250, 1.164 },
                { 50d, 0.750, 1.426 }, { 50d, 0.950, 1.646 }, { 50d, 0.995, 1.890 },
                { 100d, 0.005, 0.970 }, { 100d, 0.050, 1.077 }, { 100d, 0.250, 1.196 },
                { 100d, 0.750, 1.380 }, { 100d, 0.950, 1.527 }, { 100d, 0.995, 1.683 },
            };

            for (int i = 0; i < table.GetLength(0); i++)
            {
                int n = (int)table[i, 0];
                double alpha = table[i, 1];
                double expectedZeta = table[i, 2];

                int df = n - 1;
                double delta = zp * Math.Sqrt(n);
                var nct = new NoncentralT(df, delta);
                double computedZeta = nct.InverseCDF(alpha) / Math.Sqrt(n);

                Assert.AreEqual(expectedZeta, computedZeta, 0.002,
                    $"Table 1 mismatch: n={n}, α={alpha}, expected={expectedZeta}, computed={computedZeta:F4}");
            }
        }

        /// <summary>
        /// Validates NoncentralT InverseCDF against Stedinger (1983) Table 3.
        /// Table 3: Percentage points of ζ(0.99)-distribution for the 100-year event.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_StedsingerTable3()
        {
            double zp = Normal.StandardZ(0.99);

            // Exact values from Stedinger (1983) Table 3 (p=0.99)
            var table = new[,]
            {
                { 10d, 0.005, 1.232 }, { 10d, 0.050, 1.562 }, { 10d, 0.250, 2.008 },
                { 10d, 0.750, 2.927 }, { 10d, 0.950, 3.981 }, { 10d, 0.995, 5.582 },
                { 20d, 0.005, 1.483 }, { 20d, 0.050, 1.749 }, { 20d, 0.250, 2.085 },
                { 20d, 0.750, 2.697 }, { 20d, 0.950, 3.295 }, { 20d, 0.995, 4.059 },
                { 50d, 0.005, 1.744 }, { 50d, 0.050, 1.936 }, { 50d, 0.250, 2.163 },
                { 50d, 0.750, 2.538 }, { 50d, 0.950, 2.862 }, { 50d, 0.995, 3.230 },
                { 100d, 0.005, 1.894 }, { 100d, 0.050, 2.040 }, { 100d, 0.250, 2.207 },
                { 100d, 0.750, 2.470 }, { 100d, 0.950, 2.684 }, { 100d, 0.995, 2.915 },
            };

            for (int i = 0; i < table.GetLength(0); i++)
            {
                int n = (int)table[i, 0];
                double alpha = table[i, 1];
                double expectedZeta = table[i, 2];

                int df = n - 1;
                double delta = zp * Math.Sqrt(n);
                var nct = new NoncentralT(df, delta);
                double computedZeta = nct.InverseCDF(alpha) / Math.Sqrt(n);

                Assert.AreEqual(expectedZeta, computedZeta, 0.002,
                    $"Table 3 mismatch: n={n}, α={alpha}, expected={expectedZeta}, computed={computedZeta:F4}");
            }
        }

        /// <summary>
        /// Tests NoncentralT CDF and InverseCDF with large noncentrality parameters
        /// encountered in real confidence interval calculations (δ > 20).
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_LargeNoncentrality()
        {
            double[] probs = { 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99 };

            // δ = z_0.99 * √100 ≈ 23.3
            var nct1 = new NoncentralT(99, Normal.StandardZ(0.99) * Math.Sqrt(100));
            foreach (double p in probs)
            {
                double x = nct1.InverseCDF(p);
                Assert.IsFalse(double.IsNaN(x), $"InverseCDF({p}) returned NaN for δ≈23.3");
                Assert.IsFalse(double.IsInfinity(x), $"InverseCDF({p}) returned Infinity for δ≈23.3");
                double roundtrip = nct1.CDF(x);
                Assert.AreEqual(p, roundtrip, 1e-4, $"CDF(InverseCDF({p})) roundtrip failed for δ≈23.3");
            }

            // δ = z_0.99 * √200 ≈ 32.9
            var nct2 = new NoncentralT(199, Normal.StandardZ(0.99) * Math.Sqrt(200));
            foreach (double p in probs)
            {
                double x = nct2.InverseCDF(p);
                Assert.IsFalse(double.IsNaN(x), $"InverseCDF({p}) returned NaN for δ≈32.9");
                Assert.IsFalse(double.IsInfinity(x), $"InverseCDF({p}) returned Infinity for δ≈32.9");
                double roundtrip = nct2.CDF(x);
                Assert.AreEqual(p, roundtrip, 1e-4, $"CDF(InverseCDF({p})) roundtrip failed for δ≈32.9");
            }
        }

        /// <summary>
        /// Verifies that NoncentralTConfidenceIntervals produces results equivalent to
        /// MonteCarloConfidenceIntervals. The noncentral-t method is exact (Stedinger 1983);
        /// the Monte Carlo method is approximate.
        /// </summary>
        [TestMethod()]
        public void Test_NoncentralT_CI_vs_MonteCarlo_CI()
        {
            var dist = new Normal(100, 15);
            int sampleSize = 30;
            var quantiles = new double[] { 0.01, 0.10, 0.50, 0.90, 0.99 };
            var percentiles = new double[] { 0.05, 0.25, 0.50, 0.75, 0.95 };

            var nctCI = dist.NoncentralTConfidenceIntervals(sampleSize, quantiles, percentiles);
            var mcCI = dist.MonteCarloConfidenceIntervals(sampleSize, 100000, quantiles, percentiles);

            for (int i = 0; i < quantiles.Length; i++)
            {
                for (int j = 0; j < percentiles.Length; j++)
                {
                    double nctVal = nctCI[i, j];
                    double mcVal = mcCI[i, j];
                    double tol = Math.Abs(nctVal) > 1.0 ? 0.02 * Math.Abs(nctVal) : 0.5;
                    Assert.AreEqual(nctVal, mcVal, tol,
                        $"NCT vs MC: quantile={quantiles[i]}, percentile={percentiles[j]}: NCT={nctVal:F4}, MC={mcVal:F4}");
                }
            }
        }
    }
}
