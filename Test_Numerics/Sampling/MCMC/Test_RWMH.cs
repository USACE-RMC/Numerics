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

using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit test for the Random Walk Metropolis-Hastings (RWMH) sampler. 
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
    public class Test_RWMH
    {

        /// <summary>
        /// This test compares the results obtained using RWMH with those from the 'rstan' package. 
        /// </summary>
        [TestMethod]
        public void Test_RWMH_NormalDist_RStan()
        {

            // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
            // Table 5.1.1 Tippecanoe River Near Delphi, Indiana (Station 43) Data
            double[] sample = new double[] { 6290d, 2700d, 13100d, 16900d, 14600d, 9600d, 7740d, 8490d, 8130d, 12000d, 17200d, 15000d, 12400d, 6960d, 6500d, 5840d, 10400d, 18800d, 21400d, 22600d, 14200d, 11000d, 12800d, 15700d, 4740d, 6950d, 11800d, 12100d, 20600d, 14600d, 14600d, 8900d, 10600d, 14200d, 14100d, 14100d, 12500d, 7530d, 13400d, 17600d, 13400d, 19200d, 16900d, 15500d, 14500d, 21900d, 10400d, 7460d };

            // Create uniform priors
            var normDist = new Normal();
            var constraints = normDist.GetParameterConstraints(sample);
            var muPrior = new Uniform(constraints.Item2[0], constraints.Item3[0]);
            var sigmaPrior = new Uniform(constraints.Item2[1], constraints.Item3[1]);
            var priors = new List<IUnivariateDistribution> { muPrior, sigmaPrior };

            // Create log-likelihood function
            double logLH(double[] x)
            {
                var dist = new Normal(x[0], x[1]);
                return dist.LogLikelihood(sample);
            }

            // Create and run sampler
            var sampler = new RWMH(priors, logLH, new Matrix(2));
            sampler.Initialize = MCMCSampler.InitializationType.MAP;
            sampler.Sample();
            var results = new MCMCResults(sampler);

            /* Below are the results from 'rstan' using comparable MCMC settings:
            *            mean se_mean     sd       5%      50%      95% n_eff Rhat
            *  mu    12663.69    7.10 706.60 11488.50 12671.08 13801.45  9897    1
            *  sigma  4844.09    5.22 519.08  4077.80  4796.63  5771.81  9880    1
            *  lp__   -466.13    0.01   1.03  -468.17  -465.81  -465.15  9958    1
            * 
            *  Since MCMC methods rely on random number generation, results will not be 
            *  exactly the same as those produced by other samplers. Therefore, these 
            *  comparisons aim to verify whether the results are within 5% of 'rstan' results. 
            */

            // Mu 
            Assert.AreEqual(12663.69, results.ParameterResults[0].SummaryStatistics.Mean, 0.05 * 12663.69);
            Assert.AreEqual(706.60, results.ParameterResults[0].SummaryStatistics.StandardDeviation, 0.05 * 706.60);
            Assert.AreEqual(11488.50, results.ParameterResults[0].SummaryStatistics.LowerCI, 0.05 * 11488.50);
            Assert.AreEqual(12671.08, results.ParameterResults[0].SummaryStatistics.Median, 0.05 * 12671.08);
            Assert.AreEqual(13801.45, results.ParameterResults[0].SummaryStatistics.UpperCI, 0.05 * 13801.45);
            // Sigma 
            Assert.AreEqual(4844.09, results.ParameterResults[1].SummaryStatistics.Mean, 0.05 * 4844.09);
            Assert.AreEqual(519.08, results.ParameterResults[1].SummaryStatistics.StandardDeviation, 0.05 * 519.08);
            Assert.AreEqual(4077.80, results.ParameterResults[1].SummaryStatistics.LowerCI, 0.05 * 4077.80);
            Assert.AreEqual(4796.63, results.ParameterResults[1].SummaryStatistics.Median, 0.05 * 4796.63);
            Assert.AreEqual(5771.81, results.ParameterResults[1].SummaryStatistics.UpperCI, 0.05 * 5771.81);
        }

    }
}
