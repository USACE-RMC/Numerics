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
using Numerics.Distributions;

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

    }
}
