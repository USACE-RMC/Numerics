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
    /// Unit tests for the Polynomial class. All interpolation values were validated against R's
    /// "polyinterp( )" function of the "cmna" package
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
    /// Howard JP (2017). Computational Methods for Numerical Analysis with R, series Numerical Analysis and Scientific 
    /// Computing Series. Chapman and Hall/CRC, New York. https://jameshoward.us/cmna/.
    /// </remarks>
    [TestClass]
    public class Test_Polynomial
    {
        /// <summary>
        /// Test the base class Interpolater's sequential search function
        /// </summary>
        [TestMethod()]
        public void Test_Sequential()
        {
            var values = new double[1000];
            for (int i = 1; i <= 1000; i++)
                values[i - 1] = i;
            var poly = new Polynomial(3, values, values);
            var lo = poly.SequentialSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var poly2 = new Polynomial(3, values, values, SortOrder.Descending);
            lo = poly2.SequentialSearch(872.5);
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
            var poly = new Polynomial(3, values, values);
            var lo = poly.BisectionSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var poly2 = new Polynomial(3, values, values, SortOrder.Descending);
            lo = poly2.BisectionSearch(872.5);
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
            var poly = new Polynomial(3, values, values);
            var lo = poly.HuntSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            var poly2 = new Polynomial(3, values, values, SortOrder.Descending);
            lo = poly2.HuntSearch(872.5);
            Assert.AreEqual(lo, 127);
        }

        /// <summary>
        /// Test the polynomial class and interpolation function with an order 3 polynomial
        /// </summary>
        [TestMethod]
        public void Test_Interpolate_Order3()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var poly = new Polynomial(3, XArray, YArray);
            double X = 8d;
            double Y = poly.Interpolate(X);
            Assert.AreEqual(Y, 11.5415808882467, 1E-6);
        }

        /// <summary>
        /// Test with an order 3 polynomial interpolation with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Interpolate_Order3_R()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var poly = new Polynomial(3, XArray, YArray);
            
            var Xout = new double[] { 8, 20, 30, 56 };
            var trueYout = new double[] { 11.54158, 19.80796, 25.24398, 34.46549 };
            for(int i = 0; i < Xout.Length; i++)
            {
                double Y = poly.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout[i], Y, 1E-5);
            }
        }

        /// <summary>
        /// Test with an order 3 polynomial interpolation with a large list of values
        /// </summary>
        [TestMethod]
        public void Test_Interpolate_List()
        {
            var XArray = new double[] { 6d, 24d, 48d, 72d };
            var YArray = new double[] { 9.96d, 22.13d, 32.27d, 37.60d };
            var poly = new Polynomial(3, XArray, YArray);
            var X = new double[34];
            X[0] = 6;
            for (int i = 1; i < X.Length; i++)
                X[i] = X[i - 1] + 2;

            var Y = poly.Interpolate(X);
            var true_Y = new double[] { 9.95999999999906, 11.5415808882467, 13.0626606341182, 14.5245941558435, 15.9287363716526, 17.2764421997752, 18.5690665584413, 19.807964365881, 20.994490540324, 22.1300000000003, 23.2158476631398, 24.2533884479724, 25.2439772727281, 26.1889690556368, 27.0897187149283, 27.9475811688326, 28.7639113355797, 29.5400641333994, 30.2773944805217, 30.9772572951765, 31.6410074955936, 32.2700000000031, 32.8655897266348, 33.4291315937186, 33.9619805194845, 34.4654914221624, 34.9410192199823, 35.3899188311739, 35.8135451739673, 36.2132531665923, 36.5903977272789, 36.9463337742571, 37.2824162257566, 37.6000000000075, 37.9004400152396, 38.1850911896829, 38.4553084415673, 38.7124466891227, 38.9578608505791, 39.1929058441663 };
            for (int i = 0; i < X.Length; i++)
                Assert.AreEqual(Y[i], true_Y[i], 1E-6);
        }
    }
}
