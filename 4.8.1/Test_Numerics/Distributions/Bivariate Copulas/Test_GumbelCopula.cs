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
    /// Unit tests for the Gumbel Copula. All tests are compared against the R 'copula' package. 
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
    public class Test_GumbelCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new GumbelCopula(10);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.559981e-07, 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 7.954156e-16, 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.9004251, 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.001609717, 1E-6);

            // Test symmetry
            copula = new GumbelCopula(10);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.559981e-07, 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 7.954156e-16, 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.9004251, 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.001609717, 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new GumbelCopula(10);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.2, 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.2, 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1694668, 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1999967, 1E-6);

            // Test symmetry
            copula = new GumbelCopula(10);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.2, 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.2, 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1694668, 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1999967, 1E-6);

        }

        private double[] data1 = new double[] { 112.472889925526, 108.197032461095, 113.64666688938, 113.965212991957, 117.333562148345, 100.559345907851, 102.180053944561, 110.97067894608, 104.077049999534, 102.473973961948, 107.113945888296, 106.082540016705, 87.0978719961949, 75.1150155332284, 112.802882878613, 138.120743213339, 79.1789109173273, 106.623115297302, 110.92063481489, 109.543820203885, 73.138959971343, 98.3431545632918, 102.952511132798, 93.0465814076382, 111.688341287489, 79.2581834425453, 103.249260938183, 85.8316965436919, 84.6132137548954, 110.962644400059, 84.0087454428379, 107.223748328677, 121.506255972552, 123.73645847787, 79.0481183395772, 146.264725342603, 95.3276983471542, 107.603210240998, 99.9974503587899, 90.8379276433018, 104.67679932997, 73.4569986162982, 113.260163491206, 97.0605626966906, 89.0374440090811, 133.729158410007, 78.9190177023254, 108.345190856042, 108.953076977455, 97.0704392084539, 86.2765535028879, 102.219521216094, 63.2511232136125, 76.7529568677982, 74.1684311723453, 84.3685879177063, 94.7059136752836, 92.5916995697557, 107.923148783474, 105.261642868879, 84.7155767681601, 105.132628194701, 102.500956437438, 88.5763924520429, 95.9492723947193, 131.375734529198, 123.260221112319, 73.5681907570129, 93.1887735334748, 106.007277431115, 100.356675515534, 91.4930966213857, 122.640033750305, 107.143857399285, 91.0914154936143, 69.1482456631603, 122.434236175527, 111.906044824347, 111.632367993799, 84.6813208809395, 99.4489764877588, 106.884486176956, 95.3274691248045, 107.839259276592, 90.3967465603013, 76.5485833533036, 98.4847265047609, 96.2342050639354, 91.963420241494, 101.351308748786, 95.4725609428008, 145.010119874994, 88.3665385050844, 109.704283131981, 81.2884314113531, 109.836099370176, 102.386450396465, 93.0092098477109, 83.7158456419459, 89.7486244037052 };
        private double[] data2 = new double[] { 101.222248773938, 81.7313758478523, 123.411631899952, 93.2413907814405, 102.50494017697, 96.5306123268502, 18.1227322903384, 74.307009995493, 82.2705528425937, 89.6693353083137, 87.2256413770307, 81.1744664866471, 81.2165207417765, 81.4927985171313, 100.651709224046, 130.599135003894, 45.4442163124226, 66.6919540534553, 83.5881809451671, 88.4770416926437, 75.1039197272015, 52.1553833221993, 49.528146406319, 66.3900605051234, 97.2816608378485, 26.5831021290638, 88.5107383079526, 78.4085272967254, 70.7569274155408, 109.703373315445, 87.4165170640218, 85.6182893046897, 110.738685438064, 119.768351328997, 21.6161419679033, 154.357497576718, 53.8795395966417, 75.943071862762, 105.35244272571, 86.5868773757059, 76.4752235253919, 48.5832730671902, 84.6605006576596, 77.1067893953142, 46.6141185543174, 126.161273608663, 67.6305173666572, 94.9907772124782, 112.133285294602, 72.9506266499651, 37.128394091496, 63.4492833651505, 59.0759448691672, 56.0401720089048, 74.4223335518519, 43.9584062661215, 64.1266855983013, 103.095921282224, 107.649317301385, 72.9193080139833, 91.8847445373411, 75.1454442126996, 76.9540119546376, 64.9699411254905, 55.4012215444807, 118.83302426895, 106.302264714345, 48.1415318738189, 62.2914940564682, 84.5812336961348, 80.3688979491351, 47.5308700972038, 107.759633266225, 95.6723525891828, 27.2716316191787, 99.1054148896251, 82.2568018516306, 96.0432032275742, 96.3923733319187, 62.287551336599, 47.541023973346, 95.9802695805108, 80.2482799857865, 88.743034379118, 51.9583360652995, 73.0816610190799, 33.729130931116, 55.4771104609277, 93.3789037599054, 90.0234978870232, 62.3635712815712, 146.062306846747, 54.4992333698561, 89.1013264737916, 80.1329566564616, 123.235472480686, 73.4015267834993, 49.6295804684941, 74.6326460427444, 64.7088795224836 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new GumbelCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(1.978417, copula.Theta, 1E-4);
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
            BivariateCopula copula = new GumbelCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(2.097953, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new GumbelCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(2.031034, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new GumbelCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(1.999167, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(99.90953, ((Normal)copula.MarginalDistributionX).Mu, 1E-1);
            Assert.AreEqual(15.80668, ((Normal)copula.MarginalDistributionX).Sigma, 1E-1);
            // Marginal-Y
            Assert.AreEqual(78.36136, ((Normal)copula.MarginalDistributionY).Mu, 1E-1);
            Assert.AreEqual(24.83673, ((Normal)copula.MarginalDistributionY).Sigma, 1E-1);
        }

    }
}
