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
using Numerics.Data;

namespace Data.Interpolation
{
    /// <summary>
    /// Unit tests for the Linear class. All interpolation values were validated with R's "approx( )" function
    /// and extrapolation values were validated with R's "approxExtrap( )" function on the "Hmisc" package.
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
    /// Harrell FE (2014): Hmisc: A package of miscellaneous R functions. Programs available from https://hbiostat.org/R/Hmisc/.
    /// </remarks>
    [TestClass]
    public class Test_Linear
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
            var LI = new Linear(values, values);
            var lo = LI.SequentialSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            LI = new Linear(values, values, SortOrder.Descending);
            lo = LI.SequentialSearch(872.5);
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
            var LI = new Linear(values, values);
            var lo = LI.BisectionSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            LI = new Linear(values, values, SortOrder.Descending);
            lo = LI.BisectionSearch(872.5);
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
            var LI = new Linear(values, values);
            var lo = LI.HuntSearch(872.5d);
            Assert.AreEqual(lo, 871);

            Array.Reverse(values);
            LI = new Linear(values, values, SortOrder.Descending);
            lo = LI.HuntSearch(872.5);
            Assert.AreEqual(lo, 127);
        }

        /// <summary>
        /// Test the most basic implementation of the linear class and its interpolation function with one value
        /// </summary>
        [TestMethod]
        public void Test_Lin()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            var LI = new Linear(XArray, YArray);
            double X = 75d;
            double Y = LI.Interpolate(X);
            Assert.AreEqual(Y, 150.0d, 1E-6);
        }

        /// <summary>
        /// Test the most basic implementation of the linear class and its interpolation function with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Lin_List1()
        {
            var XArray = new double[] { 50, 100, 150, 200, 250, 300 };
            var YArray = new double[] { 0.001, 0.01, 0.1, 0.7, 0.95, 0.99 };
            var LI = new Linear(XArray, YArray);
            var Xout = new double[] { 76, 80, 96, 162, 170, 216 };
            var trueYout = new double[] { 0.00568, 0.0064, 0.00928, 0.244, 0.34, 0.78 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = LI.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout[i], Y, 1E-6);
            }
        }

        /// <summary>
        /// Test the  all of the possible log transforms within the linear class and its interpolation function
        /// </summary>
        [TestMethod]
        public void Test_Log()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            double X = 75d;

            var LinLog = new Linear(XArray, YArray) { YTransform = Transform.Logarithmic };
            double Y1 = LinLog.Interpolate(X);
            Assert.AreEqual(Y1, 141.42135623731d, 1E-6);

            var LogLin = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic };
            double Y2 = LogLin.Interpolate(X);
            Assert.AreEqual(Y2, 158.496250072116d, 1E-6);
;
            var LogLog = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            double Y3 = LogLog.Interpolate(X);
            Assert.AreEqual(Y3, 150.0d, 1E-6);
        }

        /// <summary>
        /// Test the  all of the possible log transforms within the linear class and its interpolation function
        /// with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Log_List()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            var Xout = new double[] { 76, 80, 96, 162, 170, 216 };

            var LinLog = new Linear(XArray, YArray) { YTransform = Transform.Logarithmic };
            var trueYout1 = new double[] { 143.3955, 151.5717, 189.2115, 321.4449, 336.5865, 429.6068 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y1 = LinLog.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout1[i], Y1, 1E-4);
            }

            var LogLin = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic };
            var trueYout2 = new double[] { 160.4071d, 167.8072d, 194.1106d, 326.7521d, 343.5075d, 434.4895d };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = LogLin.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout2[i], Y, 1E-4);
            }

            var LogLog = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            var trueYout3 = new double[] { 152, 160, 192, 324, 340, 432 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = LogLog.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout3[i], Y, 1E-6);
            }
        }

        /// <summary>
        /// Test the  all of the possible normal z transforms within the linear class and its interpolation function
        /// </summary>
        [TestMethod]
        public void Test_Z()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            double X = 0.18d;

            var LinZ = new Linear(XArray, YArray) { YTransform = Transform.NormalZ };
            double Y1 = LinZ.Interpolate(X);
            Assert.AreEqual(Y1, 0.358762529d, 1E-6);

            var ZLin = new Linear(XArray, YArray) { XTransform = Transform.NormalZ };
            double Y2 = ZLin.Interpolate(X);
            Assert.AreEqual(Y2, 0.362146174d, 1E-6);

            var ZZ = new Linear(XArray, YArray) { XTransform = Transform.NormalZ, YTransform = Transform.NormalZ };
            double Y3 = ZZ.Interpolate(X);
            Assert.AreEqual(Y3, 0.36093855992815d, 1E-6);
        }

        /// <summary>
        /// Test the  all of the possible normal z transforms within the linear class and its interpolation function
        /// with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Z_R()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.001, 0.01, 0.7, 0.95, 0.99 };
            var Xout = new double[] { 0.07, 0.09, 0.12, 0.18, 0.21, 0.23 };

            var LinZ = new Linear(XArray, YArray) { YTransform = Transform.NormalZ };
            var trueYout1 = new double[] { 0.002679041, 0.006585261, 0.117801570, 0.884282862, 0.962556228, 0.980000061 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = LinZ.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout1[i], Y, 1E-6);
            }

            var ZLin = new Linear(XArray, YArray) { XTransform = Transform.NormalZ };
            var trueYout2 = new double[] { 0.005188150, 0.008533366, 0.309976505, 0.855365435, 0.958424502, 0.974597253 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = ZLin.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout2[i], Y, 1E-6);
            }

            var ZZ = new Linear(XArray, YArray) { XTransform = Transform.NormalZ, YTransform = Transform.NormalZ };
            var trueYout3 = new double[] { 0.003121301, 0.007126363, 0.138520808, 0.888903833, 0.963143032, 0.980487623 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = ZZ.Interpolate(Xout[i]);
                Assert.AreEqual(trueYout3[i], Y, 1E-6);
            }
        }

        /// <summary>
        /// Test the linear class and its interpolation function with the arrays reversed
        /// </summary>
        [TestMethod]
        public void Test_RevLin()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };         
            Array.Reverse(XArray);
            Array.Reverse(YArray);
            var LI = new Linear(XArray, YArray, SortOrder.Descending);
            double X = 75d;
            double Y = LI.Interpolate(X);
            Assert.AreEqual(Y, 150.0d, 1E-6);
        }

        /// <summary>
        /// Test the linear class and its interpolation function with the arrays reversed and with all possible
        /// log transforms
        /// </summary>
        [TestMethod]
        public void Test_Rev_Log()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);
            double X = 75d;

            var LinLog = new Linear(XArray, YArray, SortOrder.Descending) { YTransform = Transform.Logarithmic };
            double Y1 = LinLog.Interpolate(X);
            Assert.AreEqual(Y1, 141.42135623731d, 1E-6);

            var LogLin = new Linear(XArray, YArray, SortOrder.Descending) { XTransform = Transform.Logarithmic };
            double Y2 = LogLin.Interpolate(X);
            Assert.AreEqual(Y2, 158.496250072116d, 1E-6);

            var LogLog = new Linear(XArray, YArray, SortOrder.Descending) { XTransform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            double Y3 = LogLog.Interpolate(X);
            Assert.AreEqual(Y3, 150.0d, 1E-6);
        }

        /// <summary>
        /// Test the linear class and its interpolation function with the arrays reversed and with all possible
        /// normal z transforms
        /// </summary>
        [TestMethod]
        public void Test_Rev_Z()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);
            double X = 0.18d;

            var LinZ = new Linear(XArray, YArray, SortOrder.Descending) { YTransform = Transform.NormalZ };
            double Y1 = LinZ.Interpolate(X);
            Assert.AreEqual(Y1, 0.358762529d, 1E-6);

            var ZLin = new Linear(XArray, YArray, SortOrder.Descending) { XTransform = Transform.NormalZ };
            double Y2 = ZLin.Interpolate(X);
            Assert.AreEqual(Y2, 0.362146174d, 1E-6);

            var ZZ = new Linear(XArray, YArray, SortOrder.Descending) { XTransform = Transform.NormalZ, YTransform = Transform.NormalZ };
            double Y3 = ZZ.Interpolate(X);
            Assert.AreEqual(Y3, 0.36093855992815d, 1E-6);
        }

        // ???
        [TestMethod]
        public void Test_Lin_List()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            var LI = new Linear(XArray, YArray) { UseSmartSearch = true };

            int N = 10;
            double delta = (260 - 40) / ((double)N - 1);
            var xVals = new double[N];
            xVals[0] = 40;
            for (int i = 1; i < N; i++)
                xVals[i] = xVals[i - 1] + delta;

            var yVals = LI.Interpolate(xVals);
            var trueVals = new double[] { 100, 128.888888888889, 177.777777777778, 226.666666666667, 275.555555555556, 324.444444444444, 373.333333333333, 422.222222222222, 471.111111111111, 500 };
            for (int i = 1; i < N; i++)
                Assert.AreEqual(yVals[i], trueVals[i], 1E-6);
        }

        /// <summary>
        /// Test the extrapolate function with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Extrapolate()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            var LI = new Linear(XArray, YArray);

            var Xout = new double[] { 300, 310, 350, 400, 440 };
            var trueYout = new double[] { 600, 620, 700, 800, 880 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y = LI.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout[i], Y, 1E-6);
            }
        }

        /// <summary>
        /// Test all of the possible log transforms within the linear class and its extrapolation function
        /// with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Extrapolate_Log()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            var Xout = new double[] { 300, 310, 350, 400, 440 };

            var LogLin = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic };
            var trueYout1 = new double[] { 581.7059, 596.4004, 650.7873, 710.6284, 753.3409 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y1 = LogLin.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout1[i], Y1, 1E-4);
            }

            var LinLog = new Linear(XArray, YArray) { YTransform = Transform.Logarithmic };
            var trueYout2 = new double[] { 625.0000,  653.5247,  781.2500,  976.5625, 1167.4225 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y2 = LinLog.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout2[i], Y2, 1E-4);
            }

            var LogLog = new Linear(XArray, YArray) { XTransform = Transform.Logarithmic, YTransform = Transform.Logarithmic };
            var trueYout3 = new double[] { 600, 620, 700, 800, 880 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y3 = LogLog.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout3[i], Y3, 1E-6);
            }
        }

        /// <summary>
        /// Test all of the possible normal z transforms within the linear class and its extrapolation function
        /// with a list of values
        /// </summary>
        [TestMethod]
        public void Test_Extrapolate_Z()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            var Xout = new double[] { 0.3, 0.45, 0.5, 0.6, 0.77 };

            var ZLin = new Linear(XArray, YArray) { XTransform = Transform.NormalZ };
            var trueYout1 = new double[] { 0.5898031, 0.8283812, 0.9035683, 1.0551538, 1.3456435 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y1 = ZLin.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout1[i], Y1, 1E-6);
            }

            var LinZ = new Linear(XArray, YArray) { YTransform = Transform.NormalZ };
            var trueYout2 = new double[] { 0.6000000, 0.8445627, 0.8973751, 0.9619212, 0.9957908 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y2 = LinZ.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout2[i], Y2, 1E-6);
            }

            var ZZ = new Linear(XArray, YArray) { XTransform = Transform.NormalZ, YTransform = Transform.NormalZ };
            var trueYout3 = new double[] { 0.5899878, 0.7972798, 0.8467110, 0.9202071, 0.9839199 };
            for (int i = 0; i < Xout.Length; i++)
            {
                double Y3 = ZZ.Extrapolate(Xout[i]);
                Assert.AreEqual(trueYout3[i], Y3, 1E-6);
            }
        }

    }
}
