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
using Numerics.Data;
using System;

namespace Data.Interpolation
{
    /// <summary>
    /// Unit tests for the CubicSpling class. The output values of the class were tested against
    /// R's  "cubicspline( )" function of the "pracma" package.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> References: </b>
    /// Hans W. Borchers (2011). pracma: Practical Numerical Math Functions. R package version 2.4.4,
    /// https://cran.r-project.org/web/packages/pracma
    /// </remarks>
    [TestClass]
    public class Test_CubicSpline
    {
        /// <summary>
        /// Test the base class Interpolater's sequential search function
        /// </summary>
        [TestMethod]
        public void Test_Sequential()
        {
            var values = new double[1000];
            for (int i = 1; i <= 1000; i++)
                values[i - 1] = i;
            var spline = new CubicSpline(values, values);
            var lo = spline.SequentialSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var spline2 = new CubicSpline(values, values, SortOrder.Descending);
            lo = spline2.SequentialSearch(872.5);
            Assert.AreEqual(lo, 127);
        }

        /// <summary>
        /// Test the base class Interpolater's bisection search function
        /// </summary>
        [TestMethod]
        public void Test_Bisection()
        {
            var values = new double[1000];
            for (int i = 1; i <= 1000; i++)
                values[i - 1] = i;
            var spline = new CubicSpline(values, values);
            var lo = spline.BisectionSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var spline2 = new CubicSpline(values, values, SortOrder.Descending);
            lo = spline2.BisectionSearch(872.5);
            Assert.AreEqual(lo, 127);
        }

        /// <summary>
        /// Test the base class Interpolater's Hunt search function
        /// </summary>
        [TestMethod]
        public void Test_Hunt()
        {
            var values = new double[1000];
            for (int i = 1; i <= 1000; i++)
                values[i - 1] = i;
            var spline = new CubicSpline(values, values);
            var lo = spline.HuntSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var spline2 = new CubicSpline(values, values, SortOrder.Descending);
            lo = spline2.HuntSearch(872.5);
            Assert.AreEqual(lo, 127);
        }

        /// <summary>
        /// Test the interpolate function with one value
        /// </summary>
        [TestMethod]
        public void Test_Interpolate()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var spline = new CubicSpline(XArray, YArray);
            double X = 8d;
            double Y = spline.Interpolate(X);
            Assert.AreEqual(Y, 11.4049889205445d, 1E-6);
        }

        /// <summary>
        /// Test the interpolate function with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Interpolate_R()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var spline = new CubicSpline(XArray, YArray);

            var Xout = new double[] { 8, 20, 30, 56 };
            var trueYout = new double[] { 11.40499, 19.68530, 25.35189, 34.35289 };
            for(int i = 0; i < Xout.Length; i++)
            {
                double Y = spline.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout[i], Y, 1E-4);
            }
        }

        /// <summary>
        /// Test the interpolate function with another list of values
        /// </summary>
        [TestMethod]
        public void Test_Interpolate_List()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var spline = new CubicSpline(XArray, YArray);
            var X = new double[34];
            X[0] = 6;
            for (int i = 1; i < X.Length; i++)
                X[i] = X[i - 1] + 2;

            var Y = spline.Interpolate(X);
            var true_Y = new double[] { 9.96, 11.4049889205445, 12.8430203387148, 14.2671367521368, 15.6703806584362, 17.045794555239, 18.3864209401709, 19.6853023108579, 20.9354811649256, 22.13, 23.2634521159782, 24.3366340218424, 25.3518930288462, 26.3115764482431, 27.2180315912868, 28.0736057692308, 28.8806462933286, 29.6415004748338, 30.358515625, 31.0340390550807, 31.6704180763295, 32.27, 32.8352193880579, 33.3688598053181, 33.8737920673077, 34.3528869895537, 34.8090153875831, 35.2450480769231, 35.6638558731007, 36.0683095916429, 36.4612800480769, 36.8456380579297, 37.2242544367284, 37.6};
            for (int i = 0; i < X.Length; i++)
                Assert.AreEqual(Y[i], true_Y[i], 1E-6);
        }

    }
}
