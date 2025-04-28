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
    /// Unit tests for the Clayton Copula. All tests are compared against the R 'copula' package. 
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
    public class Test_ClaytonCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new ClaytonCopula(10);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.3113E-05, 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 2.387424E-11, 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.25, 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.00731365, 1E-6);

            // Test symmetry
            copula = new ClaytonCopula(10);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.3113E-05, 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 2.387424E-11, 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.25, 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.00731365, 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new ClaytonCopula(10);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.2, 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.2, 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1167184, 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1999737, 1E-6);

            // Test symmetry
            copula = new ClaytonCopula(10);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.2, 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.2, 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1167184, 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1999737, 1E-6);

        }

        private double[] data1 = new double[] { 135.852695757514, 104.082298038859, 108.737560538974, 99.2685553595344, 134.734517358047, 90.96671487189, 77.2933084429689, 115.40695335195, 108.98103448672, 79.030140378764, 119.344309239779, 114.057958622873, 121.376980352823, 102.545934287125, 116.171688724822, 90.2737860301091, 118.384838925781, 78.5650043184329, 122.502611362073, 119.21456701156, 122.901928011128, 79.7224942508733, 103.26308186608, 134.856147348969, 87.2443704188882, 107.484484922021, 102.369642440458, 106.295354390286, 103.333259725228, 114.944097101105, 101.947395165995, 99.3368989207062, 107.600964955177, 63.8796678112047, 122.606004825966, 96.2971006035691, 101.928179675011, 98.6811421254304, 136.414461816426, 94.4488702685262, 118.178109208696, 99.3447054496099, 73.6895736287131, 111.543756772209, 94.1765631568512, 108.378057968653, 82.0682187268786, 93.2316741233556, 92.1848572845123, 101.62266736115, 90.9160104057016, 97.6045250944967, 123.668316008262, 91.5790391025517, 118.507497005719, 116.028361864995, 106.431274472801, 88.1845431623773, 100.593879884171, 106.477111540899, 125.501606137489, 120.558072482172, 93.9426575106378, 61.491148690679, 81.7991169923448, 70.0617929783816, 124.844586043561, 97.6329727079425, 98.4538343555991, 78.7796746466831, 93.4697095163469, 94.4204666650014, 121.959329690142, 86.5743545885648, 119.317279401324, 78.9585491320436, 115.814312639539, 102.440872193275, 124.220180616412, 88.2386416291415, 111.637544083286, 74.24445197777, 139.719855097958, 111.842726007709, 106.593685525292, 95.7079419062174, 100.768445085655, 89.9934014730302, 89.7346792611329, 90.5836289523697, 95.6118805500478, 98.0110615835651, 104.685614771651, 105.840599667288, 110.492236632241, 115.285805520348, 112.406028227278, 97.3674746819745, 77.7088711413582, 102.843691800664 };
        private double[] data2 = new double[] { 106.247204356978, 119.634538634978, 72.5075040420413, 96.4624345242022, 93.7887427218662, 88.0710430691384, 54.3965822738992, 133.939640413503, 147.313192863588, 39.596862389244, 113.685219500925, 65.661791172481, 125.494512719209, 122.353356065013, 86.4850623431384, 56.2427898884044, 112.239490666668, 48.2002798009377, 112.836791609701, 70.306705761088, 59.4614291231979, 44.3647320315126, 91.7637292129578, 95.6179940736635, 60.4815479415048, 58.3786398648642, 94.1724393628596, 87.9631903747507, 91.5086027340312, 114.720801315278, 89.8113482015865, 82.2111665540629, 74.6453503599532, 43.9664396683432, 74.8202107513964, 98.3796818314241, 105.279101457561, 97.8452040931319, 99.7804654911693, 75.7974672384221, 139.913874523229, 81.687017386353, 34.7670214588405, 68.5923527448465, 93.0469463410792, 90.1397446225404, 51.1775088910909, 55.8146471548922, 77.1273200963725, 66.6745681621479, 95.5181914247239, 61.5430213312507, 76.4727428009896, 64.5965344713655, 75.1297368849306, 117.072661743052, 79.6250066166689, 37.2231394673498, 71.1212261596989, 76.6759488420322, 107.428580054692, 110.536313019896, 87.680895839141, 30.0011658085059, 48.6775096662048, 21.1838108219624, 95.5909448097503, 83.8207210333862, 116.670732017431, 40.4019864082416, 91.469151821309, 70.7404242103512, 45.6819094938675, 86.6712212708917, 121.431171230518, 38.4195748844222, 143.699998594758, 64.2498614632223, 75.9822870091017, 59.6123708381618, 113.933833528108, 42.0682009456935, 73.2020173163658, 51.1044236832251, 98.8788896427372, 46.7142903461546, 80.8155805864023, 87.4214895184672, 62.741943393305, 80.361655652615, 70.7251173704377, 89.8903637671805, 74.0531075745106, 72.4690834186899, 107.394378394182, 89.3201273021353, 54.0735001346141, 66.3099470947684, 55.5401474570749, 63.1916628514916 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new ClaytonCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(1.326613, copula.Theta, 1E-4);
        }

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
            BivariateCopula copula = new ClaytonCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(1.534016, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new ClaytonCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(1.485167, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new ClaytonCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(1.512075, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(102.4226, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(16.70995, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(80.54073, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(26.11687, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }
    }
}
