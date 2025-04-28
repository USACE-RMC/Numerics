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
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Joe Copula. All tests are compared against the R 'copula' package. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_JoeCopula
    {

        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new JoeCopula(10);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 4.342727e-05, 1E-6);

            copula = new JoeCopula(20);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 8.645443e-11, 1E-6);

            copula = new JoeCopula(1.1);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.9512265, 1E-6);

            copula = new JoeCopula(5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.02110735, 1E-6);

            // Test symmetry
            copula = new JoeCopula(10);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 4.342727e-05, 1E-6);

            copula = new JoeCopula(20);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 8.645443e-11, 1E-6);

            copula = new JoeCopula(1.1);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.9512265, 1E-6);

            copula = new JoeCopula(5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.02110735, 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new JoeCopula(10);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1999999, 1E-6);

            copula = new JoeCopula(20);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.2, 1E-6);

            copula = new JoeCopula(1.1);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1656223, 1E-6);

            copula = new JoeCopula(5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.199895, 1E-6);

            // Test symmetry
            copula = new JoeCopula(10);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1999999, 1E-6);

            copula = new JoeCopula(20);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.2, 1E-6);

            copula = new JoeCopula(1.1);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1656223, 1E-6);

            copula = new JoeCopula(5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.199895, 1E-6);

        }

        private double[] data1 = new double[] { 105.080527801019, 115.026257454374, 117.358129524541, 84.3469768155338, 84.5362454182564, 115.832013303087, 103.698610148398, 88.8411298970013, 94.7422294382998, 109.79929825752, 94.0976045955494, 112.417576846617, 110.180042813334, 105.093201216588, 136.275040650675, 91.0470309608862, 84.2474349276179, 100.535501865457, 80.9058461473664, 119.213220466183, 100.572669469592, 103.131258346566, 100.810317630538, 114.061624071483, 83.5975402855976, 89.0115333472513, 66.3268165070727, 84.7329927826247, 69.7680817872407, 115.650963119626, 92.5472957212446, 106.114531288404, 108.032885546296, 110.38969090409, 99.6051228521885, 71.2293045424982, 93.7583750323159, 88.56905726009, 101.688577044379, 88.1152932378495, 102.179951148463, 147.661448558484, 130.057283273608, 109.854007734265, 81.2705409811028, 117.66598298105, 88.8769701106725, 102.802427495432, 93.4008884560194, 90.2277365987588, 103.036719364918, 77.5184782680141, 91.5324658931925, 101.835989095901, 101.800043749889, 109.764650578377, 79.7834884356497, 80.6891470677998, 116.183753310048, 97.7251444184232, 101.93962603208, 91.479747920405, 106.616497717471, 88.9732889869284, 100.347155748138, 64.2776497298067, 97.8630230693096, 115.039994756502, 100.160509753111, 99.0744120656348, 111.633096916686, 99.3446937251762, 113.010250612283, 91.3689738704099, 93.0857066194397, 112.761306730353, 111.792394529326, 51.8176166632693, 107.422232145344, 102.659904000871, 82.1117724410236, 99.3163646335942, 113.555145795822, 116.501871333875, 95.4376498231323, 118.137910178593, 101.824522891831, 98.675229136211, 87.3054670597499, 107.50132629382, 88.6217672049816, 108.587609616577, 117.070978807494, 100.829377608782, 85.9716443690607, 87.1615186321506, 101.929033409522, 129.887046688029, 106.165881910737, 110.123691041869 };
        private double[] data2 = new double[] { 65.3627311420053, 109.338093441946, 110.187130564663, 72.6572491437794, 66.1830652216622, 106.917958610739, 89.4544458566931, 34.2121510489786, 87.5612581119189, 60.5637237767687, 33.8034654209032, 99.2844647305164, 104.97303205735, 6.04617859719163, 138.811261697691, 71.2893972974008, 59.417366566882, 85.355111220543, 50.1413604581927, 102.920060652383, 93.0914947236245, 70.1509627285045, 11.6186828960635, 117.810957879192, 79.1966502763724, 70.2131986632539, 72.5040290039477, 63.6963355882566, 48.4038689280691, 105.297596037192, 87.9566889097355, 87.2650426632532, 92.1659510753509, 51.9553483976797, 69.1607656479071, 69.8590431888949, 63.1696994924036, 83.3814974744098, 77.8342636732521, 54.3557740053106, 75.9564039272592, 153.295740162045, 125.192966195081, 84.9856299114578, 86.854354257246, 98.0159599592548, 3.77678567039993, 74.3992868304184, 61.7221737834943, 75.0682654590025, 107.067397020988, 66.3698849697294, 78.3453997857911, 83.9149525731405, 69.5344335662347, 105.933049927905, 26.5327793191245, 54.4009262051242, 89.3584519067692, 77.5365813897511, 66.5148556530152, 75.574562907648, 69.8755523955032, 75.15359174989, 55.0871908145974, 76.301388721602, 84.6034786632286, 86.412293254073, 89.8760641259544, 73.4531945546952, 103.905973398154, 76.6907578549793, 81.9676998718682, 78.9370335717157, 20.7797553590814, 96.0968795553399, 86.6087907592008, 70.4507347850656, 107.294695066862, 79.0020795665946, 61.8228011685283, 70.9975241569947, 76.1327812702326, 74.3139042973535, 56.2911048478489, 121.512535634237, 92.9142895046758, 88.7950084218728, 63.705929356528, 84.3367159226979, 66.7354953690805, 102.711710544246, 97.0527842334408, 69.4933842704733, 84.4978321017289, 74.067511177418, 75.1669385237322, 134.226620060869, 107.962687902957, 97.4269493483909 };

        /// <summary>
        /// Test fitting with the method of maximum pseudo likelihood.
        /// </summary>
        [TestMethod]
        public void Test_MPL_Fit()
        {
            // get ranks of data
            var rank1 = Statistics.RanksInPlace(data1);
            var rank2 = Statistics.RanksInPlace(data2);
            // get plotting positions
            for (int i = 0; i < data1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }
            // Fit copula
            BivariateCopula copula = new JoeCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(2.664326, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new JoeCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(2.965269, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new JoeCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(3.032809, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(99.85959, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(15.47272, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(79.50469, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(25.11105, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }

    }
}
