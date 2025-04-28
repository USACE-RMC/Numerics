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
    /// Unit tests for the Normal Copula. All tests are compared against the R 'copula' package. 
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
    public class Test_NormalCopula
    {

        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new NormalCopula(0.8);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.09803021, 1E-6);

            copula = new NormalCopula(0.1);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.928971, 1E-6);

            copula = new NormalCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.462211, 1E-6);

            copula = new NormalCopula(-0.9);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 3.208773, 1E-6);

            // Test symmetry
            copula = new NormalCopula(0.8);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.09803021, 1E-6);

            copula = new NormalCopula(0.1);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.928971, 1E-6);

            copula = new NormalCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.462211, 1E-6);

            copula = new NormalCopula(-0.9);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 3.208773, 1E-6);
        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {

            var copula = new NormalCopula(0.8);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1996831, 1E-2);

            copula = new NormalCopula(0.1);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1675602, 1E-2);

            copula = new NormalCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1128494, 1E-2);

            copula = new NormalCopula(-0.9);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.05006756, 1E-2);

            // Test symmetry
            copula = new NormalCopula(0.8);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1996831, 1E-2);

            copula = new NormalCopula(0.1);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1675602, 1E-2);

            copula = new NormalCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1128494, 1E-2);

            copula = new NormalCopula(-0.9);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.05006756, 1E-2);
        }

        private double[] data1 = new double[] { 122.094066003419, 92.8321267206161, 86.4920318705377, 87.6183663113541, 102.558777787492, 103.627475117762, 127.084948716539, 105.908684131013, 110.065795957654, 105.924647125867, 110.009738155469, 126.490833800772, 64.1264871206211, 81.3150800229481, 92.0780134395721, 106.040322550555, 113.158086143066, 117.051057784044, 127.110531266645, 108.907371862136, 105.476247114194, 108.629403495407, 98.7803988364997, 93.217925588845, 97.7219451830075, 109.178093756809, 137.69504856252, 106.884615327674, 112.139177456202, 85.7416217661797, 71.0610938629716, 112.644166631765, 119.545871678548, 70.5169833274982, 99.6896817997206, 100.987892854545, 103.659280253554, 75.6075621013066, 118.810868919796, 109.113664695226, 113.636425353944, 100.008375355612, 113.178917359795, 80.4269472604342, 88.3638384448237, 90.2905074656314, 98.7995143316863, 98.4698060067802, 108.279297570816, 86.1578437055905, 101.183725242941, 85.5531148952956, 111.024195253862, 121.934506174556, 104.169993666179, 84.4652994609478, 99.6099259747033, 95.3130792386208, 115.45680252817, 120.213139478586, 95.5691788140058, 92.7950300448044, 102.58430893827, 86.7105161576407, 82.8059368562185, 107.335705516294, 112.603259240932, 102.780778760832, 128.958090528336, 105.139162595628, 118.272661482198, 99.8275937885748, 94.2856024560543, 108.48679008009, 100.147734981682, 88.7006383425785, 89.6441478272035, 112.24266306884, 99.8184811468069, 120.592090049738, 124.023170133661, 101.250961381805, 90.0000027551006, 108.781064635426, 94.9203320035987, 99.9491821782837, 88.7473944659517, 94.3643253649856, 105.814317118952, 92.6866900633813, 111.020330544613, 111.676189456988, 115.70235103978, 124.659106152655, 81.3866270495082, 120.178528245778, 93.6511977805724, 114.099368762143, 119.062045395294, 74.1998497412903 };
        private double[] data2 = new double[] { 127.869024514059, 53.5970265830273, 35.6871183968043, 77.5937820397885, 84.619117510857, 110.477376636164, 114.679535976765, 109.338354392258, 88.5987759167264, 72.6695216679034, 111.932652280673, 86.3677960278751, 23.9336347978345, 51.2377830227977, 82.4565771813309, 92.9162733515069, 117.465381827514, 104.862362549521, 131.059266136887, 67.2743851584176, 100.263235166171, 113.734275000025, 73.1582387829997, 78.4353197703676, 60.0180359279642, 106.709991071405, 123.175455301514, 98.7006449949188, 99.860486991242, 55.7603096813567, 53.7716423706874, 104.659445447656, 119.899401349887, 59.8670226375024, 94.0117104763717, 101.424610891155, 114.256354904191, 53.5051841563538, 118.35993465227, 73.1605008375787, 87.4677698350712, 75.4031529479113, 105.404958657365, 53.336411944238, 61.2731445424292, 72.377272009744, 88.959659863884, 80.1301183393358, 98.624093971352, 81.9603074727622, 52.0788199186743, 75.49358652998, 90.2428259997917, 101.326931349259, 48.1343463500222, 56.9295918059918, 89.0348875829931, 69.0012535890253, 100.355241744174, 74.00820280539, 63.9482913881998, 64.4973782209222, 95.8934144135508, 85.4028102356618, 37.8958459664423, 99.2194777630975, 126.581868541047, 91.8287794302242, 143.543939198862, 108.751405708845, 100.951567564812, 73.5051068155712, 83.419507788205, 84.9090133796832, 59.3886126411711, 84.0348703304947, 78.1503115396303, 104.953483626903, 77.6450718557069, 117.615613165515, 118.131904013699, 76.3190144944821, 62.0183469143453, 97.4729901076061, 49.3396925267253, 58.6790714873228, 45.0596168059506, 85.3857426310419, 65.0772008397323, 58.8836242438228, 79.2838406333912, 102.608529398935, 83.7509120512927, 103.106132785215, 52.8403092456187, 88.4802383528401, 64.2906982187616, 93.0489784548541, 116.065815369284, 26.2779209375887 };

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
            BivariateCopula copula = new NormalCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(0.800082, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new NormalCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(0.7871479, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new NormalCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(0.7871358, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(102.5959, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(14.25443, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(84.17168, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(24.58027, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }
    }
}
