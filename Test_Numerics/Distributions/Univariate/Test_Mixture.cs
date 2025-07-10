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
using Numerics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    ///<summary>
    /// Unit testing for the Mixture distribution. 
    ///</summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Mixture
    {
        /// <summary>
        /// Test the PDF and Log PDF with the R 'mistr' package.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_PDF()
        {
            var n1 = new Normal(10, 2);
            var n2 = new Normal(20, 1);
            var n3 = new Normal(30, 5);
            var weights = new double[] { 0.3, 0.2, 0.5 };
            var dists = new IUnivariateDistribution[] { n1, n2, n3 };
            var mix = new Mixture(weights, dists);

            var xVals = Tools.Sequence(0d, 60d, 1d);
            var truePDF = new double[] { 0.0000002236155, 0.000002399538, 0.00002008072, 0.000130921, 0.0006648309, 0.002629394, 0.008099041, 0.01942865, 0.0362981, 0.05281569, 0.05985473, 0.05283899, 0.0363568, 0.01955086, 0.008337055, 0.003072727, 0.001483088, 0.002375569, 0.01305772, 0.051944, 0.08518778, 0.05628918, 0.02189028, 0.01585912, 0.01944537, 0.02419737, 0.02896916, 0.03332246, 0.03682701, 0.03910427, 0.03989423, 0.03910427, 0.03682701, 0.03332246, 0.02896916, 0.02419707, 0.01941861, 0.01497275, 0.01109208, 0.007895016, 0.005399097, 0.003547459, 0.002239453, 0.001358297, 0.0007915452, 0.0004431848, 0.0002384088, 0.0001232219, 0.00006119019, 0.00002919469, 0.00001338302, 0.000005894307, 0.000002494247, 0.000001014085, 0.0000003961299, 0.000000148672, 0.00000005361035, 0.00000001857362, 0.000000006182621, 0.00000000197732, 0.0000000006075883 };
            var trueLogPDF = new double[] { -15.313338, -12.940234, -10.815751, -8.940917, -7.315978, -5.941002, -4.81601, -3.941006, -3.31599, -2.940947, -2.815835, -2.940506, -3.314374, -3.934736, -4.787045, -5.78519, -6.513629, -6.042518, -4.338376, -2.957589, -2.462897, -2.877253, -3.821713, -4.144011, -3.940146, -3.721511, -3.541524, -3.401524, -3.301524, -3.241524, -3.221524, -3.241524, -3.301524, -3.401524, -3.541524, -3.721524, -3.941524, -4.201524, -4.501524, -4.841524, -5.221524, -5.641524, -6.101524, -6.601524, -7.141524, -7.721524, -8.341524, -9.001524, -9.701524, -10.441524, -11.221524, -12.041524, -12.901524, -13.801524, -14.741524, -15.721524, -16.741524, -17.801524, -18.901524, -20.041524, -21.221524 };
        
            for (int i = 0; i < xVals.Length; i++)
            {
                Assert.AreEqual(truePDF[i], mix.PDF(xVals[i]), 1E-5);
                Assert.AreEqual(trueLogPDF[i], mix.LogPDF(xVals[i]), 1E-5);
            }
        
        }

        /// <summary>
        /// Test the CDF and Log CDF with the R 'mistr' package.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_CDF()
        {
            var n1 = new Normal(10, 2);
            var n2 = new Normal(20, 1);
            var n3 = new Normal(30, 5);
            var weights = new double[] { 0.3, 0.2, 0.5 };
            var dists = new IUnivariateDistribution[] { n1, n2, n3 };
            var mix = new Mixture(weights, dists);

            var xVals = Tools.Sequence(0d, 60d, 1d);
            var trueCDF = new double[] { 0.00000008648877, 0.00000102096, 0.000009506731, 0.00006980538, 0.0004050192, 0.001863043, 0.006825436, 0.02004322, 0.04759928, 0.09256793, 0.1500158, 0.2074749, 0.252483, 0.2801263, 0.2935185, 0.2988121, 0.3008789, 0.3025308, 0.3086393, 0.3386818, 0.411375, 0.4862341, 0.5228496, 0.5401083, 0.5575285, 0.5793276, 0.6059277, 0.6371266, 0.6722891, 0.7103701, 0.75, 0.7896299, 0.8277109, 0.8628734, 0.8940723, 0.9206724, 0.9424652, 0.9596217, 0.9726004, 0.9820348, 0.9886249, 0.9930483, 0.9959012, 0.9976694, 0.9987224, 0.9993251, 0.9996564, 0.9998315, 0.9999204, 0.9999638, 0.9999842, 0.9999933, 0.9999973, 0.9999989, 0.9999996, 0.9999999, 1, 1, 1, 1, 1 };
            var trueLogCDF = new double[] { -16.26325, -13.79477, -11.56351, -9.569799, -7.811576, -6.285544, -4.987099, -3.909865, -3.044938, -2.379812, -1.897014, -1.572745, -1.376411, -1.272515, -1.225815, -1.20794, -1.201047, -1.195572, -1.175582, -1.082694, -0.8882501, -0.7210651, -0.6484614, -0.6159855, -0.5842417, -0.5458872, -0.5009946, -0.450787, -0.3970668, -0.3419691, -0.2876821, -0.236191, -0.1890914, -0.1474872, -0.1119686, -0.08265104, -0.05925632, -0.04121617, -0.02778202, -0.01812849, -0.01144026, -0.006976, -0.004107191, -0.002333314, -0.001278382, -0.0006751769, -0.000343628, -0.0001684788, -0.00007955746, -0.00003617468, -0.00001583575, -0.000006672897, -0.000002706276, -0.000001056228, -0.0000003966642, -0.0000001433258, -0.00000004982213, -0.00000001666022, -0.000000005358795, -0.000000001657873, -0.0000000004932939 };
            
            for (int i = 0; i < xVals.Length; i++)
            {
                Assert.AreEqual(trueCDF[i], mix.CDF(xVals[i]), 1E-5);
                Assert.AreEqual(trueLogCDF[i], mix.LogCDF(xVals[i]), 1E-5);
            }

        }

        /// <summary>
        /// Test the inverse CDF of the mixture.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_InverseCDF()
        {
            var n1 = new Normal(10, 2);
            var n2 = new Normal(20, 1);
            var n3 = new Normal(30, 5);
            var weights = new double[] { 0.3, 0.2, 0.5 };
            var dists = new IUnivariateDistribution[] { n1, n2, n3 };
            var mix = new Mixture(weights, dists);

            var xVals = Tools.Sequence(0d, 60d, 1d);
            for (int i = 0; i < xVals.Length; i++)
            {
                Assert.AreEqual(xVals[i], mix.InverseCDF(mix.CDF(xVals[i])), 1E-6);
            }

        }

        /// <summary>
        /// Test the first two moments of the mixture.
        /// </summary>
        [TestMethod]
        public void Test_Mixture_Moments()
        {
            var n1 = new Normal(10, 2);
            var n2 = new Normal(20, 1);
            var n3 = new Normal(30, 5);
            var weights = new double[] { 0.3, 0.2, 0.5 };
            var dists = new IUnivariateDistribution[] { n1, n2, n3 };
            var mix = new Mixture(weights, dists);

            // There are known solutions to the mean and standard deviation of a mixture
            // μ_XYZ= ω_X∙μ_X + ω_Y∙μ_Y + ω_Z∙μ_Z
            // σ_XYZ= √(ω_X∙(μ_X^2 + σ_X^2) + ω_Y∙(μ_Y^2 + σ_Y^2) + ω_Z∙(μ_Z^2 + σ_Z^2) - μ_XYZ^2)
            var w = weights;
            var u = new double[] { n1.Mean, n2.Mean, n3.Mean };
            var s = new double[] { n1.Variance, n2.Variance, n3.Variance };
            var trueMean = w[0] * u[0] + w[1] * u[1] + w[2] * u[2];
            var trueSD = Math.Sqrt(w[0] * (u[0] * u[0] + s[0]) + w[1] * (u[1] * u[1] + s[1]) + w[2] * (u[2] * u[2] + s[2]) - trueMean * trueMean);

            // Numerics computes the moments using numerical integration
            Assert.AreEqual(trueMean, mix.Mean, 1E-4);
            Assert.AreEqual(trueSD, mix.StandardDeviation, 1E-4);

        }

        /// <summary>
        /// Test the MLE fit for a mixture of 2 Normal distributions.
        /// </summary>
        [TestMethod]
        public void Test_2D_Mixture_MLE()
        {
            var mix = new Mixture(new[] { 0.3, 0.7 }, new[] { new Normal(0, 1), new Normal(3, 0.1) });
            var sample = mix.GenerateRandomValues(1000, 12345);
            mix.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);

            var w = mix.Weights;
            var d = mix.Distributions;
            Array.Sort(w, d);

            // Weights
            Assert.AreEqual(0.3, w[0], 1E-1);
            Assert.AreEqual(0.7, w[1], 1E-1);
            // N1 parameters
            Assert.AreEqual(0.0, d[0].GetParameters[0], 1E-1);
            Assert.AreEqual(1.0, d[0].GetParameters[1], 1E-1);
            // N2 parameters
            Assert.AreEqual(3.0, d[1].GetParameters[0], 1E-1);
            Assert.AreEqual(0.1, d[1].GetParameters[1], 1E-1);

        }

        /// <summary>
        /// Test the MLE fit for a zero-inflated mixture of 2 Normal distributions.
        /// </summary>
        [TestMethod]
        public void Test_2D_Mixture_MLE_ZeroInflated()
        {
            var mix = new Mixture(new[] { 0.3, 0.6 }, new[] { new Normal(3, 0.1), new Normal(5, 2) });
            mix.IsZeroInflated = true;
            mix.ZeroWeight = 0.1;
            var sample = mix.GenerateRandomValues(1000, 12345);
            mix.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);

            var w = mix.Weights;
            var d = mix.Distributions;
            Array.Sort(w, d);

            // Weights
            Assert.AreEqual(0.3, w[0], 1E-1);
            Assert.AreEqual(0.6, w[1], 1E-1);
            // N1 parameters
            Assert.AreEqual(3.0, d[0].GetParameters[0], 1E-1);
            Assert.AreEqual(0.1, d[0].GetParameters[1], 1E-1);
            // N2 parameters
            Assert.AreEqual(5.0, d[1].GetParameters[0], 1E-1);
            Assert.AreEqual(2.0, d[1].GetParameters[1], 1E-1);

        }

        /// <summary>
        /// Test the MLE fit for a mixture of 3 Normal distributions.
        /// </summary>
        [TestMethod]
        public void Test_3D_Mixture_MLE()
        {
            var mix = new Mixture(new[] { 0.2, 0.3, 0.5 }, new[] { new Normal(0, 1), new Normal(3, 0.1), new Normal(5, 2) });
            var sample = mix.GenerateRandomValues(1000, 12345);
            mix.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);

            var w = mix.Weights;
            var d = mix.Distributions;
            Array.Sort(w, d);

            // Weights
            Assert.AreEqual(0.2, w[0], 1E-1);
            Assert.AreEqual(0.3, w[1], 1E-1);
            Assert.AreEqual(0.5, w[2], 1E-1);
            // N1 parameters
            Assert.AreEqual(0.0, d[0].GetParameters[0], 1E-1);
            Assert.AreEqual(1.0, d[0].GetParameters[1], 1E-1);
            // N2 parameters
            Assert.AreEqual(3.0, d[1].GetParameters[0], 1E-1);
            Assert.AreEqual(0.1, d[1].GetParameters[1], 1E-1);
            // N3 parameters
            Assert.AreEqual(5.0, d[2].GetParameters[0], 1E-1);
            Assert.AreEqual(2.0, d[2].GetParameters[1], 1E-1);

        }

    }
}
