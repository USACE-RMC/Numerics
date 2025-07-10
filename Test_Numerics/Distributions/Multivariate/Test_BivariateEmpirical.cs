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

namespace Distributions.Multivariate
{
    /// <summary>
    /// Unit tests for the bivariate empirical distribution. 
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
    public class Test_BivariateEmpirical
    {
        /// <summary>
        /// Test the Bivariate empirical distribution CDF.
        /// </summary>
        [TestMethod()]
        public void Test_BivariateEmp()
        {
            var X1vals = new double[] { 1d, 2d, 3d, 4d, 5d };
            var X2vals = new double[] { 6d, 7d, 8d, 9d, 10d };
            var Pvals = new double[,] { { 0.00405294623516298d, 0.0132662170105167d, 0.0207236588305977d, 0.022603272182165d, 0.0227468879771625d }, { 0.0132662170105167d, 0.0625140947096637d, 0.127398206576625d, 0.154872951858602d, 0.158508394165442d }, { 0.0207236588305977d, 0.127398206576625d, 0.333333333333333d, 0.468742952645168d, 0.497973526882419d }, { 0.022603272182165d, 0.154872951858602d, 0.468742952645168d, 0.745203586846751d, 0.831860831130881d }, { 0.0227468879771625d, 0.158508394165442d, 0.497973526882419d, 0.831860831130881d, 0.958552682338805d } };
            var bv = new BivariateEmpirical(X1vals, X2vals, Pvals);

            // Real space
            double true_cdf1 = 0.008659582d;
            double true_cdf2 = 0.388089576d;
            double true_cdf3 = 0.497681221d;
            double cdf1 = bv.CDF(1.5d, 6d);
            double cdf2 = bv.CDF(3.22d, 8.15d);
            double cdf3 = bv.CDF(3d, 9.99d);
            Assert.AreEqual(cdf1, true_cdf1, 1E-6);
            Assert.AreEqual(cdf2, true_cdf2, 1E-6);
            Assert.AreEqual(cdf3, true_cdf3, 1E-6);

            // Log 10 space
            double true_cdf4 = 0.379680109d;
            bv.ProbabilityTransform = Numerics.Data.Transform.Logarithmic;
            double cdf4 = bv.CDF(3.22d, 8.15d);
            Assert.AreEqual(cdf4, true_cdf4, 1E-6);

            // Z space
            double true_cdf5 = 0.386806419d;
            bv.ProbabilityTransform = Numerics.Data.Transform.NormalZ;
            double cdf5 = bv.CDF(3.22d, 8.15d);
            Assert.AreEqual(cdf5, true_cdf5, 1E-6);
        }

    }
}
