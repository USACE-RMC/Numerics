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
    /// Testing the Pert Percentile distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_PertPercentileDists
    {

        /// <summary>
        /// Verified using Palisade's @Risk
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentile_AtRisk()
        {
            double true_min = -2.6453d;
            double true_max = 2.6453d;
            int true_mean = 0;
            int true_mode = 0;
            int true_median = 0;
            double true_stdDev = 0.9998d;
            int true_skew = 0;
            double true_kurt = 2.3333d;
            double true_icdf05 = -1.644;
            double true_icdf50 = 0.0;
            double true_icdf95 = 1.644;
            var pert = new PertPercentile(true_icdf05, true_icdf50, true_icdf95);
            pert.SolveParameters();
            Assert.AreEqual(pert.Minimum, true_min, 1E-3);
            Assert.AreEqual(pert.Maximum, true_max, 1E-3);
            Assert.AreEqual(pert.Mean, true_mean, 1E-3);
            Assert.AreEqual(pert.Median, true_median, 1E-3);
            Assert.AreEqual(pert.Mode, true_mode, 1E-3);
            Assert.AreEqual(pert.StandardDeviation, true_stdDev, 1E-3);
            Assert.AreEqual(pert.Skewness, true_skew, 1E-3);
            Assert.AreEqual(pert.Kurtosis, true_kurt, 1E-3);
            Assert.AreEqual(pert.CDF(true_icdf05), 0.05, 1E-3);
            Assert.AreEqual(pert.CDF(true_icdf50), 0.5, 1E-3);
            Assert.AreEqual(pert.CDF(true_icdf95), 0.95, 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.05d), true_icdf05, 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.5d), true_icdf50, 1E-3);
            Assert.AreEqual(pert.InverseCDF(0.95d), true_icdf95, 1E-3);
        }

        /// <summary>
        /// Verified using Palisade's @Risk
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentileZ_AtRisk()
        {
            double true_min = Normal.StandardCDF(-2.6453d);
            double true_max = Normal.StandardCDF(2.6453d);
            double true_mean = Normal.StandardCDF(0);
            double true_mode = Normal.StandardCDF(0);
            double true_median = Normal.StandardCDF(0);
            double true_icdf05 = Normal.StandardCDF(-1.644);
            double true_icdf50 = Normal.StandardCDF(0.0);
            double true_icdf95 = Normal.StandardCDF(1.644);
            var pert = new PertPercentileZ(true_icdf05, true_icdf50, true_icdf95);
            pert.SolveParameters();
            Assert.AreEqual(pert.Minimum, true_min, 1E-2);
            Assert.AreEqual(pert.Maximum, true_max, 1E-2);
            Assert.AreEqual(pert.Mean, true_mean, 1E-2);
            Assert.AreEqual(pert.Median, true_median, 1E-2);
            Assert.AreEqual(pert.Mode, true_mode, 1E-2);
            Assert.AreEqual(pert.CDF(true_icdf05), 0.05, 1E-2);
            Assert.AreEqual(pert.CDF(true_icdf50), 0.5, 1E-2);
            Assert.AreEqual(pert.CDF(true_icdf95), 0.95, 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.05d), true_icdf05, 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.5d), true_icdf50, 1E-2);
            Assert.AreEqual(pert.InverseCDF(0.95d), true_icdf95, 1E-2);
        }

        /// <summary>
        /// Test several Pert Percentile configurations.
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentile()
        {
            double fifth = 12;
            double fiftieth = 21;
            double ninetyfifth = 32;
            var pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            pert.SolveParameters();
            var p1 = pert.CDF(fifth);
            var p2 = pert.CDF(fiftieth);
            var p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0;
            fiftieth = 3;
            ninetyfifth = 6;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0;
            fiftieth = 2;
            ninetyfifth = 3;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);


            fifth = 14;
            fiftieth = 29;
            ninetyfifth = 36;
            pert = new PertPercentile(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

        }


        /// <summary>
        /// Test several Pert PercentileZ configurations.
        /// </summary>
        [TestMethod()]
        public void Test_PertPercentileZ()
        {
            double fifth = 0.12;
            double fiftieth = 0.21;
            double ninetyfifth = 0.32;
            var pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            pert.SolveParameters();
            var p1 = pert.CDF(fifth);
            var p2 = pert.CDF(fiftieth);
            var p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0.01;
            fiftieth = 0.3;
            ninetyfifth = 0.6;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

            fifth = 0.01;
            fiftieth = 0.02;
            ninetyfifth = 0.03;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);


            fifth = 0.014;
            fiftieth = 0.029;
            ninetyfifth = 0.036;
            pert = new PertPercentileZ(fifth, fiftieth, ninetyfifth);
            p1 = pert.CDF(fifth);
            p2 = pert.CDF(fiftieth);
            p3 = pert.CDF(ninetyfifth);

            Assert.AreEqual(pert.CDF(fifth), 0.05, 1E-1);
            Assert.AreEqual(pert.CDF(fiftieth), 0.5, 1E-1);
            Assert.AreEqual(pert.CDF(ninetyfifth), 0.95, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.05d), fifth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.5d), fiftieth, 1E-1);
            Assert.AreEqual(pert.InverseCDF(0.95d), ninetyfifth, 1E-1);

        }

    }
}
