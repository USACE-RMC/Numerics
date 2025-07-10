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

namespace Data.PairedData
{
    /// <summary>
    /// Unit tests for the OrderedPairedData interpolation methods
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_PairedDataInterpolation
    {
        /// <summary>
        /// Test the sequential search method within the OrderedPairedData class
        /// </summary>
        [TestMethod()]
        public void Test_Sequential()
        {
            // ASC
            var opd = new OrderedPairedData(true, SortOrder.Ascending, false, SortOrder.Ascending);
            for (int i = 1; i <= 1000; i++)
                opd.Add(new Ordinate(i, i));
            // X
            var lo = opd.SequentialSearchX(872.5d);
            Assert.AreEqual(lo, 871);
            // Y
            lo = opd.SequentialSearchY(872.5d);
            Assert.AreEqual(lo, 871);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.SequentialSearchX(872.5d);
            Assert.AreEqual(lo, 127);
            // Y
            lo = opd.SequentialSearchY(872.5d);
            Assert.AreEqual(lo, 127);

        }

        /// <summary>
        /// Test the bisection search method within the OrderedPairedData class
        /// </summary>
        [TestMethod()]
        public void Test_Bisection()
        {
            // ASC
            var opd = new OrderedPairedData(true, SortOrder.Ascending, false, SortOrder.Ascending);
            for (int i = 1; i <= 1000; i++)
                opd.Add(new Ordinate(i, i));
            // X
            var lo = opd.BisectionSearchX(872.5d);
            Assert.AreEqual(lo, 871);
            // Y
            lo = opd.BisectionSearchY(872.5d);
            Assert.AreEqual(lo, 871);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.BisectionSearchX(872.5d);
            Assert.AreEqual(lo, 127);
            // Y
            lo = opd.BisectionSearchY(872.5d);
            Assert.AreEqual(lo, 127);

        }

        /// <summary>
        /// Test the Hunt search method within the OrderedPairedData class
        /// </summary>
        [TestMethod()]
        public void Test_Hunt()
        {
            // ASC
            var opd = new OrderedPairedData(true, SortOrder.Ascending, false, SortOrder.Ascending);
            for (int i = 1; i <= 1000; i++)
                opd.Add(new Ordinate(i, i));
            // X
            var lo = opd.HuntSearchX(872.5d);
            Assert.AreEqual(lo, 871);
            // Y
            lo = opd.HuntSearchY(872.5d);
            Assert.AreEqual(lo, 871);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.HuntSearchX(872.5d);
            Assert.AreEqual(lo, 127);
            // Y
            lo = opd.HuntSearchY(872.5d);
            Assert.AreEqual(lo, 127);

        }

        /// <summary>
        /// Test linear interpolation
        /// </summary>
        [TestMethod()]
        public void Test_Lin()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 75d;
            double Y = opd.GetYFromX(X);
            Assert.AreEqual(Y, 150.0d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y);
            Assert.AreEqual(X, xFromY, 1E-6);

        }

        /// <summary>
        /// Test linear interpolation with a logarithmic transform on y
        /// </summary>
        [TestMethod()]
        public void Test_LinLog()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.None, Transform.Logarithmic);
            Assert.AreEqual(Y, 141.42135623731d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.None, Transform.Logarithmic);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with a logarithmic transform on x
        /// </summary>
        [TestMethod()]
        public void Test_LogLin()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.Logarithmic, Transform.None);
            Assert.AreEqual(Y, 158.496250072116d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.Logarithmic, Transform.None);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with a logarithmic transform on x and y
        /// </summary>
        [TestMethod()]
        public void Test_LogLog()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.Logarithmic, Transform.Logarithmic);
            Assert.AreEqual(Y, 150.0d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.Logarithmic, Transform.Logarithmic);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with a normal z transform on y
        /// </summary>
        [TestMethod()]
        public void Test_LinZ()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.None, Transform.NormalZ);
            Assert.AreEqual(Y, 0.358762529d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.None, Transform.NormalZ);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with a normal z transform on x
        /// </summary>
        [TestMethod()]
        public void Test_ZLin()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.NormalZ, Transform.None);
            Assert.AreEqual(Y, 0.362146174d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.NormalZ, Transform.None);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with a normal z transform on both x and y
        /// </summary>
        [TestMethod()]
        public void Test_ZZ()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.NormalZ, Transform.NormalZ);
            Assert.AreEqual(Y, 0.36093855992815d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.NormalZ, Transform.NormalZ);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed
        /// </summary>
        [TestMethod()]
        public void Test_RevLinear()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 75d;
            double Y = opd.GetYFromX(X);
            Assert.AreEqual(Y, 150.0d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a logarithmic transform on y
        /// </summary>
        [TestMethod()]
        public void Test_RevLinLog()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.None, Transform.Logarithmic);
            Assert.AreEqual(Y, 141.42135623731d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.None, Transform.Logarithmic);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a logarithmic transform on x
        /// </summary>
        [TestMethod()]
        public void Test_RevLogLin()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.Logarithmic, Transform.None);
            Assert.AreEqual(Y, 158.496250072116d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.Logarithmic, Transform.None);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a logarithmic transform on x and y
        /// </summary>
        [TestMethod()]
        public void Test_RevLogLog()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 75d;
            double Y = opd.GetYFromX(X, Transform.Logarithmic, Transform.Logarithmic);
            Assert.AreEqual(Y, 150.0d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.Logarithmic, Transform.Logarithmic);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a normal z transform on y
        /// </summary>
        [TestMethod()]
        public void Test_RevLinZ()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.None, Transform.NormalZ);
            Assert.AreEqual(Y, 0.358762529d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.None, Transform.NormalZ);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a normal z transform on x
        /// </summary>
        [TestMethod()]
        public void Test_RevZLin()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.NormalZ, Transform.None);
            Assert.AreEqual(Y, 0.362146174d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.NormalZ, Transform.None);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation with the values reversed and a normal z transform on x and y
        /// </summary>
        [TestMethod()]
        public void Test_RevZZ()
        {
            var XArray = new double[] { 0.05d, 0.1d, 0.15d, 0.2d, 0.25d };
            var YArray = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d };
            Array.Reverse(XArray);
            Array.Reverse(YArray);

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Descending, true, SortOrder.Descending);
            double X = 0.18d;
            double Y = opd.GetYFromX(X, Transform.NormalZ, Transform.NormalZ);
            Assert.AreEqual(Y, 0.36093855992815d, 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(Y, Transform.NormalZ, Transform.NormalZ);
            Assert.AreEqual(X, xFromY, 1E-6);
        }

        /// <summary>
        /// Test linear interpolation on a list of input values
        /// </summary>
        [TestMethod()]
        public void Test_Lin_List()
        {
            var XArray = new double[] { 50d, 100d, 150d, 200d, 250d };
            var YArray = new double[] { 100d, 200d, 300d, 400d, 500d };

            // Given X
            var opd = new OrderedPairedData(XArray, YArray, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double X = 75d;
            double Y = opd.GetYFromX(X);
            Assert.AreEqual(Y, 150.0d, 1E-6);

            var yVals = opd.GetYFromX(XArray);
            for (int i = 1; i < YArray.Length; i++)
                Assert.AreEqual(YArray[i], yVals[i], 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(yVals);
            for (int i = 1; i < YArray.Length; i++)
                Assert.AreEqual(XArray[i], xFromY[i], 1E-6);
        }

    }
}
