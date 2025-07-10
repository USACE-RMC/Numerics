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
    /// Unit tests for the Frank Copula. All tests are compared against the R 'copula' package. 
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
    public class Test_FrankCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new FrankCopula(10);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.02469702, 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.0001228828, 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.089996, 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.2408784, 1E-6);

            // Test symmetry
            copula = new FrankCopula(10);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.02469702, 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.0001228828, 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.089996, 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.2408784, 1E-6);
        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new FrankCopula(10);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1998148, 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1999997, 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1534309, 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1960338, 1E-6);

            // Test symmetry
            copula = new FrankCopula(10);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1998148, 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1999997, 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1534309, 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1960338, 1E-6);
        }


        private double[] data1 = new double[] { 89.4465116086853, 91.7257825797634, 102.113218235982, 80.9258507082031, 83.4871920273991, 115.695292225839, 82.9627093424983, 88.6038265767081, 128.772275696669, 102.984538236351, 121.578135793465, 108.431406395934, 80.0293061353146, 111.431178502881, 79.2589288833566, 118.158376357189, 80.7620921129063, 108.376156927961, 95.8476352445177, 114.311862750292, 94.7747731262201, 118.875199050171, 93.2349488601981, 98.3959955105056, 108.649741186885, 122.182781815476, 77.0562107320218, 107.540344227161, 117.370406779471, 112.182445319761, 93.0289332852279, 90.529599567024, 120.201762251141, 114.278074578765, 104.568760050718, 99.7636591326382, 110.046851177132, 94.1493960120251, 104.412586917358, 98.7607116866476, 98.7463881469725, 105.1085034505, 121.269695545648, 113.731813064378, 107.800720302073, 106.36957073507, 91.6744161473404, 116.374681447127, 81.77133664558, 85.6156796128899, 107.491138655012, 102.732572344842, 121.586784481517, 116.993480194217, 124.032257210449, 99.1570143649983, 71.7128277313518, 54.6418351602958, 120.014404396034, 96.6470795565503, 107.428703111402, 99.3256773796286, 96.3098320017686, 115.238337503403, 72.0162782248935, 80.5473603316434, 85.8099029975343, 134.325558667819, 102.547186178411, 109.9543113498, 112.218655560109, 103.637934442487, 105.599530792916, 104.374711756055, 117.597527250648, 108.696922523792, 95.7635882353413, 113.556079933849, 91.7625855813032, 101.738578831165, 63.8120024793359, 63.7998925020047, 100.584112662855, 103.281311868372, 87.666605980906, 116.010424048388, 120.789094918399, 115.413149279569, 93.3814849565503, 91.1936357063422, 95.8835908172229, 104.01051720411, 118.117029822129, 101.092380522772, 105.948133484193, 101.511409332803, 91.6808235778427, 112.274907117709, 113.335593237504, 109.199879710188 };
        private double[] data2 = new double[] { 90.2725830142295, 50.9675107814588, 72.8091017436628, 37.834394607052, 66.0298908678443, 131.543666094728, 23.4646646758453, 38.387445258362, 90.0045371339763, 64.0836922901644, 113.713441516048, 86.2362874241136, 68.9809017958243, 83.2228976474751, 40.2100658509047, 126.376008778581, 79.489380523177, 105.171916388245, 92.2811081494569, 111.661707493053, 76.4955819232933, 64.4896701009032, 59.6898595730439, 95.893807000091, 84.416144203673, 119.239043275835, 41.0210587766723, 80.3403252114675, 98.987503709067, 117.127667708957, 75.0679115951468, 47.6762994744396, 98.631369967578, 145.07181315581, 73.404950592094, 77.2596264433046, 85.5789153296777, 58.5136745231103, 90.0353584106097, 88.7890539339681, 82.6531147843648, 85.5315028451663, 107.486065585877, 115.188926317826, 88.5805804644887, 94.0736037456068, 70.9824270336269, 104.904710824485, 56.6492033519975, 49.5998637369, 87.9079770735456, 85.6153352921363, 109.84039662402, 109.70332488243, 109.013245187645, 95.3425794566713, 37.0734004156272, 68.2576366613219, 83.0475074812988, 76.0326431319754, 99.6616980295915, 73.6251299875544, 58.5130416264008, 125.024141331457, 20.7835530562455, 55.6776014119347, 40.421171413744, 109.884427624024, 93.0410150569105, 98.1134492367469, 103.110412789979, 70.3696410672021, 87.9378660295918, 88.828217701758, 145.785916448632, 124.799405281758, 61.6454234189007, 74.1550632448535, 53.7379563429746, 103.133087738044, 64.6916624574593, 21.1435390667273, 77.6876521796421, 75.6302369353579, 71.5323052064531, 160.885271053561, 95.1756368878511, 119.852069907319, 68.0321254425055, 74.3625711976756, 79.8206997832619, 84.3152823622181, 113.17252314479, 54.5297187685205, 82.1517568341444, 87.0940228919095, 69.4555843222865, 81.7696161911127, 73.1591225199505, 93.8908358293613 };

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
            BivariateCopula copula = new FrankCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(7.718761, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new FrankCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(8.130259, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new FrankCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(8.231434, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(101.4276, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(14.88317, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(81.08905, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(26.93537, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }
    }
}
