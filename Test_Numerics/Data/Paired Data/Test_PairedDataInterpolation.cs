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
            Assert.AreEqual(871, lo);
            // Y
            lo = opd.SequentialSearchY(872.5d);
            Assert.AreEqual(871, lo);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.SequentialSearchX(872.5d);
            Assert.AreEqual(127, lo);
            // Y
            lo = opd.SequentialSearchY(872.5d);
            Assert.AreEqual(127, lo);

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
            Assert.AreEqual(871, lo);
            // Y
            lo = opd.BisectionSearchY(872.5d);
            Assert.AreEqual(871, lo);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.BisectionSearchX(872.5d);
            Assert.AreEqual(127, lo);
            // Y
            lo = opd.BisectionSearchY(872.5d);
            Assert.AreEqual(127, lo);

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
            Assert.AreEqual(871, lo);
            // Y
            lo = opd.HuntSearchY(872.5d);
            Assert.AreEqual(871, lo);

            // DSC
            opd = new OrderedPairedData(true, SortOrder.Descending, false, SortOrder.Descending);
            for (int i = 1000; i >= 1; i--)
                opd.Add(new Ordinate(i, i));
            // X
            lo = opd.HuntSearchX(872.5d);
            Assert.AreEqual(127, lo);
            // Y
            lo = opd.HuntSearchY(872.5d);
            Assert.AreEqual(127, lo);

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
            Assert.AreEqual(150.0d, Y, 1E-6);

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
            Assert.AreEqual(141.42135623731d, Y, 1E-6);

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
            Assert.AreEqual(158.496250072116d, Y, 1E-6);

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
            Assert.AreEqual(150.0d, Y, 1E-6);

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
            Assert.AreEqual(0.358762529d, Y, 1E-6);

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
            Assert.AreEqual(0.362146174d, Y, 1E-6);

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
            Assert.AreEqual(0.36093855992815d, Y, 1E-6);

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
            Assert.AreEqual(150.0d, Y, 1E-6);

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
            Assert.AreEqual(141.42135623731d, Y, 1E-6);

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
            Assert.AreEqual(158.496250072116d, Y, 1E-6);

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
            Assert.AreEqual(150.0d, Y, 1E-6);

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
            Assert.AreEqual(0.358762529d, Y, 1E-6);

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
            Assert.AreEqual(0.362146174d, Y, 1E-6);

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
            Assert.AreEqual(0.36093855992815d, Y, 1E-6);

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
            Assert.AreEqual(150.0d, Y, 1E-6);

            var yVals = opd.GetYFromX(XArray);
            for (int i = 1; i < YArray.Length; i++)
                Assert.AreEqual(YArray[i], yVals[i], 1E-6);

            // Given Y
            var xFromY = opd.GetXFromY(yVals);
            for (int i = 1; i < YArray.Length; i++)
                Assert.AreEqual(XArray[i], xFromY[i], 1E-6);
        }

        /// <summary>
        /// Regression pin for the default out-of-range behavior: with no extrapolation requested,
        /// lookups hold the boundary ordinate on both sides, at and beyond the endpoints, for
        /// ascending and descending tables and for transformed axes. This is the bit-identity gate
        /// for ExtrapolationSides.None.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_DefaultNone_MatchesEndpointHold()
        {
            // Ascending, no transforms.
            var asc = new OrderedPairedData(new double[] { 1d, 2d, 4d, 8d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(30d, asc.GetYFromX(3d), 1E-12);
            Assert.AreEqual(10d, asc.GetYFromX(1d), 0d);
            Assert.AreEqual(80d, asc.GetYFromX(8d), 0d);
            Assert.AreEqual(10d, asc.GetYFromX(0.5d), 0d);
            Assert.AreEqual(80d, asc.GetYFromX(16d), 0d);
            Assert.AreEqual(1d, asc.GetXFromY(5d), 0d);
            Assert.AreEqual(8d, asc.GetXFromY(100d), 0d);

            // Descending X, no transforms.
            var desc = new OrderedPairedData(new double[] { 8d, 4d, 2d, 1d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Descending, true, SortOrder.Ascending);
            Assert.AreEqual(10d, desc.GetYFromX(16d), 0d);
            Assert.AreEqual(80d, desc.GetYFromX(0.5d), 0d);

            // Logarithmic axes hold the raw boundary ordinates as well.
            var logs = new OrderedPairedData(new double[] { 1d, 10d, 100d }, new double[] { 2d, 20d, 200d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(2d, logs.GetYFromX(0.1d, Transform.Logarithmic, Transform.Logarithmic), 0d);
            Assert.AreEqual(200d, logs.GetYFromX(1000d, Transform.Logarithmic, Transform.Logarithmic), 0d);

            // Normal-z probability axis holds too.
            var normz = new OrderedPairedData(new double[] { 1d, 2d, 3d }, new double[] { 0.2d, 0.5d, 0.8d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(0.2d, normz.GetYFromX(0d, Transform.None, Transform.NormalZ), 0d);
            Assert.AreEqual(0.8d, normz.GetYFromX(4d, Transform.None, Transform.NormalZ), 0d);
        }

        /// <summary>
        /// Linear extrapolation with no transforms: each end extends the boundary segment's secant,
        /// the side flags act independently, and ascending and descending tables agree on the same
        /// underlying line.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_LinearBothEnds_NoTransform()
        {
            var asc = new OrderedPairedData(new double[] { 1d, 2d, 4d, 8d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            // Below: through (1,10)-(2,20) at x = 0.5 -> 5. Above: through (4,40)-(8,80) at x = 10 -> 100.
            Assert.AreEqual(5d, asc.GetYFromX(0.5d, extrapolation: ExtrapolationSides.Both), 1E-12);
            Assert.AreEqual(100d, asc.GetYFromX(10d, extrapolation: ExtrapolationSides.Both), 1E-12);
            // Side selectivity: the unrequested side still holds.
            Assert.AreEqual(5d, asc.GetYFromX(0.5d, extrapolation: ExtrapolationSides.Below), 1E-12);
            Assert.AreEqual(80d, asc.GetYFromX(10d, extrapolation: ExtrapolationSides.Below), 0d);
            Assert.AreEqual(10d, asc.GetYFromX(0.5d, extrapolation: ExtrapolationSides.Above), 0d);
            Assert.AreEqual(100d, asc.GetYFromX(10d, extrapolation: ExtrapolationSides.Above), 1E-12);
            // Exactly at the endpoints the boundary ordinate is returned unchanged.
            Assert.AreEqual(10d, asc.GetYFromX(1d, extrapolation: ExtrapolationSides.Both), 0d);
            Assert.AreEqual(80d, asc.GetYFromX(8d, extrapolation: ExtrapolationSides.Both), 0d);

            // The same line authored descending extrapolates to the same values: below the minimum
            // x of 1 through (2,40)-(1,80) -> 100 at x = 0.5; above the maximum x of 8 through
            // (8,10)-(4,20) -> 5 at x = 10.
            var desc = new OrderedPairedData(new double[] { 8d, 4d, 2d, 1d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Descending, true, SortOrder.Ascending);
            Assert.AreEqual(100d, desc.GetYFromX(0.5d, extrapolation: ExtrapolationSides.Both), 1E-12);
            Assert.AreEqual(5d, desc.GetYFromX(10d, extrapolation: ExtrapolationSides.Both), 1E-12);
            Assert.AreEqual(100d, desc.GetYFromX(0.5d, extrapolation: ExtrapolationSides.Below), 1E-12);
            Assert.AreEqual(10d, desc.GetYFromX(10d, extrapolation: ExtrapolationSides.Below), 0d);
        }

        /// <summary>
        /// Logarithmic-axis extrapolation extends log-linearly: the table y = 2x on log-log axes
        /// continues to y = 2x beyond both ends.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_LogTransform()
        {
            var logs = new OrderedPairedData(new double[] { 1d, 10d, 100d }, new double[] { 2d, 20d, 200d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(2000d, logs.GetYFromX(1000d, Transform.Logarithmic, Transform.Logarithmic, ExtrapolationSides.Both), 1E-9 * 2000d);
            Assert.AreEqual(0.2d, logs.GetYFromX(0.1d, Transform.Logarithmic, Transform.Logarithmic, ExtrapolationSides.Both), 1E-9 * 0.2d);
        }

        /// <summary>
        /// Normal-z extrapolation extends linearly in standard-normal z and therefore always
        /// back-transforms into (0, 1): the probability ladder 0.2/0.5/0.8 is exactly linear in z,
        /// and the far extension stays a valid probability.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_NormalZTransform()
        {
            var normz = new OrderedPairedData(new double[] { 1d, 2d, 3d }, new double[] { 0.2d, 0.5d, 0.8d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            double z05 = Numerics.Distributions.Normal.StandardZ(0.5d);
            double z08 = Numerics.Distributions.Normal.StandardZ(0.8d);
            double expected = Numerics.Distributions.Normal.StandardCDF(z08 + (z08 - z05));
            Assert.AreEqual(expected, normz.GetYFromX(4d, Transform.None, Transform.NormalZ, ExtrapolationSides.Both), 1E-12);
            // Extreme extension stays within [0, 1], saturating at the floating-point boundaries.
            double far = normz.GetYFromX(100d, Transform.None, Transform.NormalZ, ExtrapolationSides.Both);
            Assert.IsGreaterThan(0.999d, far);
            Assert.IsLessThanOrEqualTo(1d, far);
            double low = normz.GetYFromX(-100d, Transform.None, Transform.NormalZ, ExtrapolationSides.Both);
            Assert.IsGreaterThanOrEqualTo(0d, low);
            Assert.IsLessThan(0.001d, low);
        }

        /// <summary>
        /// GetXFromY extrapolates on the y-range sides with the same semantics as GetYFromX.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_GetXFromY()
        {
            var asc = new OrderedPairedData(new double[] { 1d, 2d, 4d, 8d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(0.5d, asc.GetXFromY(5d, extrapolation: ExtrapolationSides.Both), 1E-12);
            Assert.AreEqual(10d, asc.GetXFromY(100d, extrapolation: ExtrapolationSides.Both), 1E-12);
            Assert.AreEqual(1d, asc.GetXFromY(5d, extrapolation: ExtrapolationSides.Above), 0d);
            Assert.AreEqual(8d, asc.GetXFromY(100d, extrapolation: ExtrapolationSides.Below), 0d);
        }

        /// <summary>
        /// Plateau (non-strict) boundary segments extend at slope zero, which coincides with the
        /// endpoint hold, and a single-point table always holds regardless of the flags.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_PlateauAndSinglePoint()
        {
            var plateau = new OrderedPairedData(new double[] { 1d, 2d, 3d, 4d }, new double[] { 5d, 5d, 10d, 10d }, true, SortOrder.Ascending, false, SortOrder.Ascending);
            Assert.AreEqual(5d, plateau.GetYFromX(0d, extrapolation: ExtrapolationSides.Both), 0d);
            Assert.AreEqual(10d, plateau.GetYFromX(6d, extrapolation: ExtrapolationSides.Both), 0d);

            var single = new OrderedPairedData(new double[] { 2d }, new double[] { 7d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            Assert.AreEqual(7d, single.GetYFromX(100d, extrapolation: ExtrapolationSides.Both), 0d);
            Assert.AreEqual(7d, single.GetYFromX(-100d, extrapolation: ExtrapolationSides.Both), 0d);
        }

        /// <summary>
        /// The list overloads match the scalar overload point-for-point under extrapolation.
        /// </summary>
        [TestMethod]
        public void Test_Extrapolation_ListOverloads_MatchScalar()
        {
            var asc = new OrderedPairedData(new double[] { 1d, 2d, 4d, 8d }, new double[] { 10d, 20d, 40d, 80d }, true, SortOrder.Ascending, true, SortOrder.Ascending);
            var xs = new double[] { 0.5d, 1d, 3d, 8d, 10d };
            var ys = asc.GetYFromX(xs, Transform.None, Transform.None, ExtrapolationSides.Both);
            for (int i = 0; i < xs.Length; i++)
            {
                Assert.AreEqual(asc.GetYFromX(xs[i], Transform.None, Transform.None, ExtrapolationSides.Both), ys[i], 0d);
            }

            var targets = new double[] { 5d, 10d, 30d, 80d, 100d };
            var backs = asc.GetXFromY(targets, Transform.None, Transform.None, ExtrapolationSides.Both);
            for (int i = 0; i < targets.Length; i++)
            {
                Assert.AreEqual(asc.GetXFromY(targets[i], Transform.None, Transform.None, ExtrapolationSides.Both), backs[i], 0d);
            }
        }

    }
}
