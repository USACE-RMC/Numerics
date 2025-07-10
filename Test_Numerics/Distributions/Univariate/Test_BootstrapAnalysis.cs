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

namespace Distributions.Univariate
{

    /// <summary>
    /// Unit test for bootstrap analysis of univariate distributions. 
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
    public class Test_BootstrapAnalysis
    {

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Normal' method for the Normal distribution with those from the 'boot' package. 
        /// </summary>
        [TestMethod]
        public void Test_NormalCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.NormalQuantileCI(probabilities);

            /* Below are the results from 'boot' using comparable settings:
            *  Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var true05 = new double[] { 5.451061, 5.3807, 5.284459, 5.208962, 5.130887, 5.023248, 4.938038, 4.849109, 4.724954, 4.625182, 4.519401, 4.368283, 4.243235, 4.106139, 3.899265, 3.713675, 3.485485, 3.317521, 3.031074, 2.732176, 2.546064, 2.283231, 2.063337, 1.813848, 1.646703 };
            var true95 = new double[] { 6.098854, 6.010689, 5.890184, 5.795726, 5.698123, 5.563704, 5.457429, 5.34666, 5.192298, 5.068532, 4.937633, 4.751333, 4.597948, 4.430811, 4.181338, 3.961466, 3.698674, 3.512604, 3.213796, 2.927438, 2.759516, 2.531367, 2.345799, 2.138941, 2.001853 };

            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(true05[i], CIs[i, 0], 0.01 * true05[i]);
                Assert.AreEqual(true95[i], CIs[i, 1], 0.01 * true95[i]);
            }
        }

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Percentile' method for the Normal distribution with those from the 'boot' package. 
        /// </summary>
        [TestMethod]
        public void Test_PercentileCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.PercentileQuantileCI(probabilities);

            /* Below are the results from 'boot' using comparable settings:
            *  Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var true05 = new double[] { 5.448065, 5.378149, 5.280887, 5.205002, 5.127415, 5.020843, 4.935679, 4.847712, 4.723488, 4.623268, 4.517452, 4.367232, 4.24338, 4.106147, 3.900035, 3.714043, 3.486811, 3.318758, 3.032698, 2.734326, 2.548027, 2.284871, 2.06533, 1.816327, 1.648123 };
            var true95 = new double[] { 6.095642, 6.00775, 5.887441, 5.792823, 5.696038, 5.562841, 5.457508, 5.346516, 5.191728, 5.067822, 4.9371, 4.750642, 4.596812, 4.428688, 4.178886, 3.95975, 3.69779, 3.51265, 3.215176, 2.930193, 2.761301, 2.533379, 2.3472, 2.139633, 2.003285 };

            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(true05[i], CIs[i, 0], 0.01 * true05[i]);
                Assert.AreEqual(true95[i], CIs[i, 1], 0.01 * true95[i]);
            }

        }

        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'Bootstrap-t' or 'Student' method for the Normal distribution with the "true" Noncentral-t intervals. 
        /// </summary>
        [TestMethod]
        public void Test_BootstrapTCI()
        {
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.BootstrapTQuantileCI(probabilities);

            /* Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var trueCIs = dist.MonteCarloConfidenceIntervals(100, 10000, probabilities, new double[] { 0.05, 0.95 });
            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], CIs[i, 0], 0.01 * trueCIs[i, 0]);
                Assert.AreEqual(trueCIs[i, 1], CIs[i, 1], 0.01 * trueCIs[i, 1]);
            }

        }


        /// <summary>
        /// This test compares the quantile confidence intervals obtained by the 'BCa' method for the Normal distribution with the "true" Noncentral-t intervals. 
        /// </summary>
        [TestMethod]
        public void Test_BCaCI()
        {
            var sampleData = new double[] { 3.292764, 3.354733, 2.945348, 2.773251, 3.302944, 2.091022, 3.315049, 2.861908, 2.85792, 2.540339, 2.941876, 3.908656, 3.185314, 3.260108, 2.624734, 3.40845, 2.556821, 2.834211, 3.560356, 3.149362, 3.389811, 3.727893, 2.677836, 2.223431, 2.201145, 3.902549, 2.759176, 3.31019, 3.306062, 2.918845, 3.405937, 4.098417, 4.024595, 3.816223, 3.127136, 3.245594, 2.837957, 2.168975, 3.883867, 3.012901, 3.564255, 1.809821, 2.469867, 3.46857, 3.427226, 3.730365, 2.293451, 3.283702, 3.291594, 2.346601, 2.729807, 3.973846, 3.026795, 3.175831, 2.664512, 3.138977, 3.345586, 3.411898, 4.072533, 1.826528, 3.074796, 2.328734, 3.276652, 3.794981, 2.70656, 2.083811, 3.44407, 3.796744, 3.258427, 2.352164, 3.027308, 2.607675, 2.475324, 4.165256, 3.701353, 3.4713, 3.413129, 2.59423, 3.238124, 3.510629, 3.322692, 3.521572, 2.847815, 4.238555, 3.48561, 3.93355, 3.336021, 2.846023, 3.268262, 3.412435, 2.518049, 2.572459, 3.943473, 2.80409, 2.509684, 3.343666, 2.747478, 4.07886, 2.700101, 2.652727 };
            var probabilities = new double[] { 0.999999, 0.999998, 0.999995, 0.99999, 0.99998, 0.99995, 0.9999, 0.9998, 0.9995, 0.999, 0.998, 0.995, 0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01 };
            var dist = new Normal(3.122599, 0.5573654);
            var boot = new BootstrapAnalysis(dist, ParameterEstimationMethod.MethodOfMoments, 100);
            var CIs = boot.BCaQuantileCI(sampleData, probabilities);

            /* Since bootstrap methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 1% of 'boot' results. 
            */

            var trueCIs = dist.MonteCarloConfidenceIntervals(100, 10000, probabilities, new double[] { 0.05, 0.95 });
            for (int i = 0; i < probabilities.Length; i++)
            {
                Assert.AreEqual(trueCIs[i, 0], CIs[i, 0], 0.01 * trueCIs[i, 0]);
                Assert.AreEqual(trueCIs[i, 1], CIs[i, 1], 0.01 * trueCIs[i, 1]);
            }

        }

    }
}
