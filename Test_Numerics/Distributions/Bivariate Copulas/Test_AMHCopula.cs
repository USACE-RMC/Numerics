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
    /// Unit tests for the AMH Copula. All tests are compared against the R 'copula' package. 
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
    public class Test_AMHCopula
    {

        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new AMHCopula(0.8);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.6491167, 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 0.9630927, 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.158995, 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(copula.PDF(0.2, 0.8), 1.259423, 1E-6);

            // Test symmetry
            copula = new AMHCopula(0.8);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.6491167, 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 0.9630927, 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.158995, 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(copula.PDF(0.8, 0.2), 1.259423, 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new AMHCopula(0.8);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1834862, 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1626016, 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1481481, 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(copula.CDF(0.2, 0.8), 0.1398601, 1E-6);

            // Test symmetry
            copula = new AMHCopula(0.8);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1834862, 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1626016, 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1481481, 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(copula.CDF(0.8, 0.2), 0.1398601, 1E-6);

        }

        private double[] data1 = new double[] { 82.6207861517446, 122.99162719591, 85.5464857510852, 88.271076232573, 75.5957634925482, 115.666296571075, 88.1535968625871, 89.0486583921505, 95.6520223040253, 93.8704340350363, 94.8983409930905, 98.190034452137, 119.062097667798, 101.396326234446, 102.044102658515, 122.09278272123, 87.5930294017727, 89.467574539934, 68.759249687252, 97.3170891237622, 115.23726463545, 88.7249972292411, 104.808913056224, 81.4376985296622, 88.5496448465636, 85.9566777828142, 124.159708173032, 82.3989677011281, 81.0754908215474, 96.9227210098483, 81.9598285215311, 106.023194937944, 98.9429238885091, 99.6105440441767, 76.8967482790596, 93.8733402407522, 97.1593258552729, 83.8476766091024, 83.3372066079669, 84.2378120351813, 130.294157330837, 92.0183206123361, 91.0344208437507, 62.6399930777225, 100.724150023017, 105.381851336153, 106.052861517738, 105.695937045651, 98.7180560328115, 103.659278801543, 135.355695309546, 112.126374548073, 75.8894228126445, 113.303033282488, 108.572361018471, 84.4697433912218, 102.158309989318, 104.169634603778, 126.258771553477, 92.0530458784438, 121.851579456857, 112.009972090868, 82.566074288035, 100.197743002857, 109.653015549203, 103.193458190759, 120.818259910061, 109.584519216444, 102.469182012814, 107.385973965207, 97.2671206491186, 112.417149113877, 97.8701780074462, 86.9632866191297, 129.722705545036, 118.480752744013, 109.568435955988, 122.488504266404, 112.534825234329, 110.178551830369, 85.8393034892458, 128.148787670707, 88.8620889850668, 100.212554271105, 111.449083807286, 97.9922242053636, 96.0208140774823, 118.595431696388, 96.5101956950509, 89.0070522979225, 111.690497497928, 98.2256095188502, 109.751617939794, 113.262405920361, 89.1969936327915, 126.836051984505, 94.7671769565446, 101.621377036869, 85.9837125336272, 80.6830712615443 };
        private double[] data2 = new double[] { 23.3689941895676, 139.075605092248, 125.391972129225, 68.0466559628945, 21.3849501824245, 72.623739980515, 118.431110203811, 55.6433779771001, 61.7394973031962, 117.820425024964, 51.4915091268943, 90.6318678923286, 90.2225515717848, 82.9583115041602, 58.8208994489536, 89.3934020879786, 84.2687389199291, 52.8123343271364, 133.922245003475, 65.7460273014489, 94.0306797439919, 92.7622201239545, 67.4956202291755, 52.1748570092584, 76.289463467535, 30.6442121096422, 80.367744128015, 84.589322268581, 89.4931756387178, 99.8444599387594, 147.360635147879, 26.7990606303491, 75.9359921708942, 132.224006705134, 43.5713891672347, 69.4331898107351, 113.826594948864, 44.856181127981, 59.6818830690231, 34.8777185962398, 124.522970067359, 55.0531003117469, 112.641824278124, 53.3587463398865, 82.3818569319419, 131.306495042214, 92.5500215682676, 96.1266469742884, 126.22260773411, 105.118018987749, 98.5541508572697, 81.6528854407349, 77.1433525995117, 76.8338526759968, 89.723133164602, 55.8177007161922, 94.0035174340863, 164.923630526842, 133.423371650793, 89.0043059347403, 89.7587258229626, 88.4640664624267, 104.247452450998, 67.9322685363551, 83.4200659472935, 73.8431913350094, 17.7019731704449, 114.710483451006, 101.065426716036, 55.4556570293316, 67.8514034006699, 127.73989553804, 125.34525003517, 57.0971025319617, 85.1768621826704, 101.705832867004, 85.3021580666208, 89.4450158150175, 75.5941517929671, 124.424804607655, 83.3449371155293, 96.7642165625312, 57.9643876742987, 126.024606562933, 66.5425894700683, 125.964000028259, 45.1213432721552, 86.2256615507381, 101.231675185045, 27.2201147115773, 95.015298927854, 63.4674026223287, 86.224525433146, 116.818569536018, 63.082998559889, 73.7218123340297, 92.935473959462, 84.2370209154411, 49.9905902403315, 61.0021219814154 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new AMHCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(0.7965838, copula.Theta, 1E-4);
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
            BivariateCopula copula = new AMHCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(0.8321504, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new AMHCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(0.8392506, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new AMHCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(0.8367644, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(100.2479, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(14.74416, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(83.76367, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(29.76202, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }
    }

}
