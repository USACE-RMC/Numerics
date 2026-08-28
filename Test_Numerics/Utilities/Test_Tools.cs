using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;

namespace Utilities
{
    /// <summary>
    /// Testing the Tools class
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b> Authors: </b>
    /// <list type="bullet"> 
    /// <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Tools
    {
        /// <summary>
        /// Test Sign function with varying inputs.
        /// </summary>
        [TestMethod]
        public void Test_Sign()
        {
            double a = 1;
            double b = -2;
            var result = Tools.Sign(a, b);

            Assert.AreEqual(-1, result);

            double c = -1;
            double d = -2;
            var result2 = Tools.Sign(c, d);
            Assert.AreEqual(-1, result2);

            double e = -1;
            double f = 2;
            var result3 = Tools.Sign(e, f);
            Assert.AreEqual(1, result3);
        }

        /// <summary>
        /// Testing the squared value function.
        /// </summary>
        [TestMethod]
        public void Test_Sqr()
        {
            double a = 9;
            var result = Tools.Sqr(a);
            Assert.AreEqual(81, result);
        }

        /// <summary>
        /// Testing Pow function.
        /// </summary>
        [TestMethod]
        public void Test_Pow()
        {
            double a = 2;
            int b = 3;
            var result = Tools.Pow(a, b);
            Assert.AreEqual(8, result);

            double c = 123;
            int d = 0;
            var result2 = Tools.Pow(c, d);
            Assert.AreEqual(1, result2);
        }

        /// <summary>
        /// Testing Log10 function with different inputs.
        /// </summary>
        [TestMethod]
        public void Test_Log10()
        {
            double x = 100;
            var result = Tools.Log10(x);
            Assert.AreEqual(2, result);

            double x2 = 1;
            var result2 = Tools.Log10(x2);
            Assert.AreEqual(0, result2);
        }

        /// <summary>
        /// Testing natural log with different inputs.
        /// </summary>
        [TestMethod]
        public void Test_Log()
        {
            double x = 1;
            var result = Tools.Log(x);
            Assert.AreEqual(0, result);

            double x2 = 2.9;
            var result2 = Tools.Log(x2);
            Assert.AreEqual(1.0647, result2, 1E-04);
        }

        /// <summary>
        /// Testing Euclidean Distance between two points.
        /// </summary>
        [TestMethod]
        public void Test_Distance()
        {
            double x1 = 3;
            double y1 = 0;
            double x2 = 0;
            double y2 = 4;

            var result = Tools.Distance(x1, y1, x2, y2);
            Assert.AreEqual(5, result);

            List<double> point1 = new List<double> { 3d, 0 };
            List<double> point2 = new List<double> { 0, 4d };
            var result2 = Tools.Distance(point1, point2);
            Assert.AreEqual(5, result2);
        }

        /// <summary>
        /// Completing Min-Max normalization and denormalization on data
        /// </summary>
        [TestMethod]
        public void Test_Normalize_Denormalize()
        {
            var trueData = new double[] { 0, 250, 500, 750, 1000 };
            var trueNorm = new double[] { 0, 0.25, 0.5, 0.75, 1.0 };
            // Test normalization
            var dataNorm = Tools.Normalize(trueData);
            for (int i = 0; i < trueData.Length; i++)
            {
                Assert.AreEqual(trueNorm[i], dataNorm[i]);
            }
            // Test de-normalization
            var dataDenorm = Tools.Denormalize(dataNorm, Tools.Min(trueData), Tools.Max(trueData));
            for (int i = 0; i < trueData.Length; i++)
            {
                Assert.AreEqual(trueData[i], dataDenorm[i]);
            }
        }

        /// <summary>
        /// Standardizing (Z-score normalization) each value in a list.
        /// </summary>
        [TestMethod]
        public void Test_Standardize()
        {
            List<double> values = new List<double> { 3, 3, 4, 4, 6 };
            var result = Tools.Standardize(values);
            var true_result = new double[] { -0.8164, -0.8164, 0, 0, 1.63299 };
            for(int i = 0;i < values.Count; i++)
            {
                Assert.AreEqual(true_result[i], result[i], 1E-04);
            }
            
        }

        /// <summary>
        /// Destandardizing each value in a list.
        /// </summary>
        [TestMethod]
        public void Test_Destandardize()
        {
            List<double> values = new List<double> { -0.8164, -0.8164, 0, 0, 1.63299 };
            var result = Tools.Destandardize(values, 4, 1.224745);
            var true_result = new double[] { 3, 3, 4, 4, 6 };
            for (int i = 0; i < values.Count; i++)
            {
                Assert.AreEqual(true_result[i], result[i], 1E-03);
            }
        }

        /// <summary>
        /// Summing the int values in a list.
        /// </summary>
        [TestMethod]
        public void Test_SumInt()
        {
            List<int> values = new List<int> { 1, 2, 3 };
            var result = Tools.Sum(values);
            Assert.AreEqual(6, result);
        }

        /// <summary>
        /// Summing the double values in a list.
        /// </summary>
        [TestMethod]
        public void Test_SumDouble()
        {
            List<double> values = new List<double> { 1.4d, 2.3, 3.3 };
            var result = Tools.Sum(values);
            Assert.AreEqual(7, result);
        }

        /// <summary>
        /// Testing the sum of a list of values with indicator variables.
        /// </summary>
        [TestMethod]
        public void Test_SumIndicator()
        {
            List<double> values = new List<double> { 1.4, 2.3, 3.2, 3.3d };
            List<int> predictors = new List<int> { 1, 1, 1, 0 };
            var result = Tools.Sum(values,predictors);
            Assert.AreEqual(6.9, result);
        }

        /// <summary>
        /// Testing the sum of the product of two lists' values.
        /// </summary>
        [TestMethod]
        public void Test_SumProduct()
        {
            List<double> list1 = new List<double> { 1, 2, 3 };
            List<double> list2 = new List<double> { 4, 5, 6 };

            var result = Tools.SumProduct(list1, list2);
            Assert.AreEqual(32, result);
        }

        /// <summary>
        /// Testing the average of values in a list.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Mean(values);
            Assert.AreEqual(2, result);
        }

        /// <summary>
        /// Testing average of values in a list with indicator variables.
        /// </summary>
        [TestMethod]
        public void Test_MeanIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3,4 };
            List<int> indicators = new List<int> { 1, 1, 1, 0 };
            var result = Tools.Mean(values,indicators);
            Assert.AreEqual(2, result);
        }

        /// <summary>
        /// Testing product of a list of values.
        /// </summary>
        [TestMethod]
        public void Test_Product()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Product(values);
            Assert.AreEqual(6, result);
        }

        /// <summary>
        /// Testing product of a list of values with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_ProductIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3,4 };
            List<int> indicators = new List<int> { 1,1,1,0 };
            var result = Tools.Product(values,indicators);
            Assert.AreEqual(6, result);
        }

        /// <summary>
        /// Testing minimum value in a list
        /// </summary>
        [TestMethod]
        public void Test_Min()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Min(values);
            Assert.AreEqual(1, result);
        }

        /// <summary>
        /// Testing arg min value in a list
        /// </summary>
        [TestMethod]
        public void Test_ArgMin()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.ArgMin(values);
            Assert.AreEqual(0, result);
        }

        /// <summary>
        /// Testing max value in a list
        /// </summary>
        [TestMethod]
        public void Test_Max()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.Max(values);
            Assert.AreEqual(3, result);
        }


        /// <summary>
        /// Testing arg max value in a list
        /// </summary>
        [TestMethod]
        public void Test_ArgMax()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            var result = Tools.ArgMax(values);
            Assert.AreEqual(2, result);
        }

        /// <summary>
        /// Testing minimum value in a list with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_MinIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3, 4};
            List<int> indicators = new List<int> { 0, 1, 1, 1 };
            var result = Tools.Min(values,indicators);
            Assert.AreEqual(2, result);
        }

        /// <summary>
        /// Testing maximum value in a list with indicator variables
        /// </summary>
        [TestMethod]
        public void Test_MaxIndicator()
        {
            List<double> values = new List<double> { 1, 2, 3 };
            List<int> indicators = new List<int> { 1, 0, 0 };
            var result = Tools.Max(values, indicators);
            Assert.AreEqual(1, result);
        }

        /// <summary>
        /// Verifies concurrent additions are committed exactly once and return their serialization positions.
        /// </summary>
        [TestMethod]
        public void Test_ParallelAdd_CommitsEveryConcurrentUpdate()
        {
            const int updateCount = 25000;
            double total = 0d;
            var committedValues = new double[updateCount];

            Parallel.For(
                0,
                updateCount,
                index => committedValues[index] = Tools.ParallelAdd(ref total, 1d));

            Assert.AreEqual((double)updateCount, total);
            Array.Sort(committedValues);
            for (int i = 0; i < committedValues.Length; i++)
            {
                Assert.AreEqual(i + 1d, committedValues[i]);
            }
        }

        /// <summary>
        /// Verifies the compare-and-swap loop terminates when the accumulator contains NaN.
        /// </summary>
        [TestMethod]
        public void Test_ParallelAdd_NaNAccumulatorTerminates()
        {
            double total = double.NaN;

            double result = Tools.ParallelAdd(ref total, 1d);

            Assert.IsTrue(double.IsNaN(total));
            Assert.IsTrue(double.IsNaN(result));
        }

        /// <summary>
        /// Verifies signed zero, NaN, and infinity follow double-precision addition semantics.
        /// </summary>
        [TestMethod]
        public void Test_ParallelAdd_HandlesSpecialValues()
        {
            double negativeZero = BitConverter.Int64BitsToDouble(long.MinValue);

            double total = negativeZero;
            double result = Tools.ParallelAdd(ref total, negativeZero);

            Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(total));
            Assert.AreEqual(long.MinValue, BitConverter.DoubleToInt64Bits(result));

            total = negativeZero;
            result = Tools.ParallelAdd(ref total, 0d);

            Assert.AreEqual(0L, BitConverter.DoubleToInt64Bits(total));
            Assert.AreEqual(0L, BitConverter.DoubleToInt64Bits(result));

            total = 1d;
            result = Tools.ParallelAdd(ref total, double.NaN);

            Assert.IsTrue(double.IsNaN(total));
            Assert.IsTrue(double.IsNaN(result));

            total = double.PositiveInfinity;
            result = Tools.ParallelAdd(ref total, 1d);

            Assert.AreEqual(double.PositiveInfinity, total);
            Assert.AreEqual(double.PositiveInfinity, result);

            total = double.NegativeInfinity;
            result = Tools.ParallelAdd(ref total, -1d);

            Assert.AreEqual(double.NegativeInfinity, total);
            Assert.AreEqual(double.NegativeInfinity, result);

            total = double.PositiveInfinity;
            result = Tools.ParallelAdd(ref total, double.NegativeInfinity);

            Assert.IsTrue(double.IsNaN(total));
            Assert.IsTrue(double.IsNaN(result));
        }

        /// <summary>
        /// Verifies concurrent mixed-sign additions are all committed when the accumulator revisits values.
        /// </summary>
        [TestMethod]
        public void Test_ParallelAdd_CommitsMixedSignUpdates()
        {
            const int updateCount = 25000;
            double total = 0d;

            Parallel.For(
                0,
                updateCount,
                index => Tools.ParallelAdd(ref total, index % 2 == 0 ? 2d : -1d));

            Assert.AreEqual(updateCount / 2d, total);
        }

        /// <summary>
        /// Verifies a negative-to-positive zero race cannot cause finite concurrent updates to be lost.
        /// </summary>
        [TestMethod]
        public void Test_ParallelAdd_SignedZeroContentionDoesNotLoseUpdates()
        {
            const int rounds = 256;
            double total = 0d;
            bool everyUpdateCommitted = true;
            double negativeZero = BitConverter.Int64BitsToDouble(long.MinValue);

            using (var startBarrier = new Barrier(5))
            using (var finishBarrier = new Barrier(5))
            {
                Func<double, Task> createWorker = valueToAdd => Task.Factory.StartNew(
                    () =>
                    {
                        for (int round = 0; round < rounds; round++)
                        {
                            startBarrier.SignalAndWait();
                            Tools.ParallelAdd(ref total, valueToAdd);
                            finishBarrier.SignalAndWait();
                        }
                    },
                    CancellationToken.None,
                    TaskCreationOptions.LongRunning,
                    TaskScheduler.Default);

                Task[] workers =
                {
                    createWorker(0d),
                    createWorker(0d),
                    createWorker(1d),
                    createWorker(1d)
                };

                for (int round = 0; round < rounds; round++)
                {
                    total = negativeZero;
                    startBarrier.SignalAndWait();
                    finishBarrier.SignalAndWait();
                    if (total != 2d) everyUpdateCommitted = false;
                }

                Task.WaitAll(workers);
            }

            Assert.IsTrue(everyUpdateCommitted, "A finite update was lost during signed-zero contention.");
        }

        /// <summary>
        /// Testing the Log-sum-exponential function.
        /// <para> 
        /// <see href="https://www.rdocumentation.org/packages/matrixStats/versions/1.3.0/topics/logSumExp"/>
        /// </para>
        /// </summary>
        [TestMethod]
        public void Test_LogSumExp()
        {
            double u = 1000.01;
            double v = 1000.02;
            List<double> values = new List<double> { 1000.01, 1000.02 };

            var result = Tools.LogSumExp(u, v);
            var result2 = Tools.LogSumExp(values);
            Assert.AreEqual(1000.70815,result, 1E-04);
            Assert.AreEqual(1000.70815,result2, 1E-04);
        }

        /// <summary>
        /// Testing min helpers when every candidate is positive infinity.
        /// </summary>
        /// <remarks>
        /// Reference behavior verified against Python numpy 2.4.2: np.min([inf, inf]) returns inf
        /// and np.argmin([inf, inf]) returns 0.
        /// </remarks>
        [TestMethod]
        public void Test_MinHelpers_AllPositiveInfinity()
        {
            List<double> values = new List<double> { double.PositiveInfinity, double.PositiveInfinity };
            List<int> indicators = new List<int> { 1, 1 };

            Tools.MinMax(values, out double min, out double max);

            Assert.AreEqual(double.PositiveInfinity, min);
            Assert.AreEqual(double.PositiveInfinity, max);
            Assert.AreEqual(0, Tools.ArgMin(values));
            Assert.AreEqual(double.PositiveInfinity, Tools.Min(values));
            Assert.AreEqual(double.PositiveInfinity, Tools.Min(values, indicators));
        }

        /// <summary>
        /// Testing max helpers when every candidate is negative infinity.
        /// </summary>
        [TestMethod]
        public void Test_MaxHelpers_AllNegativeInfinity()
        {
            List<double> values = new List<double> { double.NegativeInfinity, double.NegativeInfinity };
            List<int> indicators = new List<int> { 1, 1 };

            Tools.MinMax(values, out double min, out double max);

            Assert.AreEqual(double.NegativeInfinity, min);
            Assert.AreEqual(double.NegativeInfinity, max);
            Assert.AreEqual(0, Tools.ArgMax(values));
            Assert.AreEqual(double.NegativeInfinity, Tools.Max(values));
            Assert.AreEqual(double.NegativeInfinity, Tools.Max(values, indicators));
        }

        /// <summary>
        /// Testing log-sum-exponential when every log input is negative infinity.
        /// </summary>
        [TestMethod]
        public void Test_LogSumExp_AllNegativeInfinity()
        {
            List<double> values = new List<double> { double.NegativeInfinity, double.NegativeInfinity };

            Assert.AreEqual(double.NegativeInfinity, Tools.LogSumExp(double.NegativeInfinity, double.NegativeInfinity));
            Assert.AreEqual(double.NegativeInfinity, Tools.LogSumExp(values));
        }

        /// <summary>
        /// Testing integer sequence.
        /// </summary>
        [TestMethod]
        public void Test_IntegerSequence()
        {
            int start = 0;
            int end = 3;
            var result = Tools.Sequence(start, end);
            var true_result = new int[] { 0, 1, 2, 3 };
            for (int i = 0;  i < true_result.Length; i++)
            {
                Assert.AreEqual(true_result[i], result[i]);
            }
        }

        /// <summary>
        /// Testing compress into byte array.
        /// </summary>
        [TestMethod]
        public void Test_Compress()
        {
            byte[] data = new byte[3];
            data[0] = 0;
            data[1] = 128;
            data[2] = 255;
            var result = Tools.Compress(data);
            Assert.IsLessThanOrEqualTo(result.Length, data.Length);
        }

        /// <summary>
        /// Testing decompress a byte array.
        /// </summary>
        [TestMethod]
        public void Test_Decompress()
        {
            byte[] data = new byte[3];
            data[0] = 0;
            data[1] = 128;
            data[2] = 255;
            var result = Tools.Decompress(data);
            Assert.IsGreaterThanOrEqualTo(result.Length, data.Length);
        }

        /// <summary>
        /// Expm1 computes exp(x) - 1 without cancellation: tiny arguments return themselves
        /// exactly, small arguments match the series exp(x) - 1 = x + x^2/2 + x^3/6 to full
        /// precision, deep negatives saturate at exactly -1, arguments large enough that the
        /// exponential itself overflows return positive infinity while the band just below, where
        /// only the compensation's intermediate product overflows, stays finite, and the round trip
        /// with Log1p closes.
        /// </summary>
        [TestMethod]
        public void Test_Expm1()
        {
            Assert.AreEqual(0d, Tools.Expm1(0d), 0d);
            Assert.AreEqual(1E-18, Tools.Expm1(1E-18), 0d);
            Assert.AreEqual(Math.E - 1d, Tools.Expm1(1d), 3E-16);
            // Series references: 1e-8 + 0.5e-16 + 1.667e-25 and its negative-argument mirror.
            Assert.AreEqual(1.0000000050000000167E-8, Tools.Expm1(1E-8), 1E-24);
            Assert.AreEqual(-9.9999999500000002E-9, Tools.Expm1(-1E-8), 1E-24);
            // Deep negatives saturate at -1: within 1E-15 at -746, and full underflow returns exactly -1.
            Assert.AreEqual(-1d, Tools.Expm1(-746d), 1E-15);
            Assert.AreEqual(-1d, Tools.Expm1(-800d), 0d);
            Assert.AreEqual(-1d, Tools.Expm1(double.NegativeInfinity), 0d);
            Assert.AreEqual(double.PositiveInfinity, Tools.Expm1(800d));
            Assert.AreEqual(double.PositiveInfinity, Tools.Expm1(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(Tools.Expm1(double.NaN)));
            // Inside the band where exp(x) is finite but the compensation's intermediate
            // (u - 1) * x overflows, the exact difference is returned: exp(x) - 1 and exp(x) agree
            // to the last unit there. Unguarded, these answered positive infinity.
            Assert.AreEqual(Math.Exp(705d), Tools.Expm1(705d), 0d);
            Assert.AreEqual(Math.Exp(709d), Tools.Expm1(709d), 0d);
            // The exponential itself first overflows just above 709.78.
            Assert.AreEqual(double.PositiveInfinity, Tools.Expm1(709.8d));

            foreach (double x in new[] { 1E-12, 0.5d, 3d })
            {
                Assert.AreEqual(x, Tools.Expm1(Tools.Log1p(x)), 1E-14 * x);
            }
        }

        /// <summary>
        /// Log1p pins for the companion helper: tiny arguments return themselves exactly, small
        /// arguments match the series log(1 + x) = x - x^2/2 + x^3/3 to full precision, and the
        /// domain edges produce negative infinity at -1 and NaN below it.
        /// </summary>
        [TestMethod]
        public void Test_Log1p()
        {
            Assert.AreEqual(0d, Tools.Log1p(0d), 0d);
            Assert.AreEqual(1E-18, Tools.Log1p(1E-18), 0d);
            Assert.AreEqual(Math.Log(2d), Tools.Log1p(1d), 3E-16);
            Assert.AreEqual(9.9999999500000003E-9, Tools.Log1p(1E-8), 1E-24);
            Assert.AreEqual(-1.00000000500000003E-8, Tools.Log1p(-1E-8), 1E-24);
            Assert.AreEqual(double.NegativeInfinity, Tools.Log1p(-1d));
            Assert.IsTrue(double.IsNaN(Tools.Log1p(-1.5d)));
        }
    }

}
