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

using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Competing Risks distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://reliability.readthedocs.io/en/latest/Competing%20risk%20models.html" />
    /// </para>
    /// </remarks>

    [TestClass]
    public class Test_CompetingRisks
    {

        /// <summary>
        /// Test the moments of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_Moments()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            // Numerics computes the moments using numerical integration
            Assert.AreEqual(27.0445, cr.Mean, 1E-2);
            Assert.AreEqual(25.0845, cr.Median, 1E-2);
            Assert.AreEqual(16.6581, cr.Mode, 1E-2);
            Assert.AreEqual(15.60225, cr.StandardDeviation, 1E-2);
            Assert.AreEqual(0.3371, cr.Skewness, 1E-2);
            Assert.AreEqual(-0.8719, cr.Kurtosis - 3, 1E-2); // Package reports excess kurtosis.

        }

        /// <summary>
        /// Test the CDF and Inverse CDF of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_CDF()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            Assert.AreEqual(4.6431, cr.InverseCDF(0.05), 1E-3);
            Assert.AreEqual(25.0845, cr.InverseCDF(0.50), 1E-3);
            Assert.AreEqual(54.2056, cr.InverseCDF(0.95), 1E-3);

        }

    }
}
