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
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the RunningStatistics class. Most of these methods were tested against the same validation methods used in testing the Statistics
    /// class.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_RunningStatistics
    {
        /// <summary>
        /// Test the construction of the running covariance matrix.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix()
        {
            var matrix = new RunningCovarianceMatrix(5);
            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 0 }, { 0 }, { 0 }, { 0 }, { 0 } });
            var trueCovariance = new Matrix(new double[,] { { 1, 0, 0, 0, 0 },
                                                            { 0, 1, 0, 0, 0 },
                                                            { 0, 0, 1, 0, 0 },
                                                            { 0, 0, 0, 1, 0 },
                                                            { 0, 0, 0, 0, 1 } });

            for(int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// Test the push method of the running covariance matrix.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix_Push()
        {
            var matrix = new RunningCovarianceMatrix(5);
            var values = new double[] { 1, 2, 3, 4, 5 };

            matrix.Push(values);
            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 1 }, { 2 }, { 3 }, { 4 }, { 5 } });
            var trueCovariance = new Matrix(new double[,] { { 1, 0, 0, 0, 0 },
                                                            { 0, 1, 0, 0, 0 },
                                                            { 0, 0, 1, 0, 0 },
                                                            { 0, 0, 0, 1, 0 },
                                                            { 0, 0, 0, 0, 1 } });

            for (int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// Test the push method of the running covariance matrix with mulitple pushes.
        /// </summary>
        [TestMethod]
        public void Test_RunningMatrix_MultiplePush()
        {
            var matrix = new RunningCovarianceMatrix(5);
            matrix.Push(new double[] { 1, 1, 1, 1, 1 });
            matrix.Push(new double[] { 2, 2, 2, 2, 2 });
            matrix.Push(new double[] { 3, 3, 3, 3, 3 });
            matrix.Push(new double[] { 4, 4, 4, 4, 4 });
            matrix.Push(new double[] { 5, 5, 5, 5, 5 });

            var mean = matrix.Mean;
            var covariance = matrix.Covariance;

            var trueMean = new Matrix(new double[,] { { 3 }, { 3 }, { 3 }, { 3 }, { 3 } });
            var trueCovariance = new Matrix(new double[,] { { 11, 10, 10, 10, 10 },
                                                            { 10, 11, 10, 10, 10 },
                                                            { 10, 10, 11, 10, 10 },
                                                            { 10, 10, 10, 11, 10 },
                                                            { 10, 10, 10, 10, 11 } });

            for (int i = 0; i < mean.NumberOfRows; i++)
            {
                CollectionAssert.AreEqual(trueMean.Row(i), mean.Row(i));
                CollectionAssert.AreEqual(trueCovariance.Row(i), covariance.Row(i));
            }
        }

        /// <summary>
        /// Test the construction of the running statistic class.
        /// </summary>
        [TestMethod]
        public void Test_RunningStat_Construction()
        {
            var runningstat1 = new RunningStatistics();

            var values = new double[] { 1, 2, 3, 4, 5 };
            var runningstat2 = new RunningStatistics(values);

            var test1 = runningstat1.Count;
            Assert.AreEqual(0, test1);

            var test2 = runningstat2.Count;
            Assert.AreEqual(5, test2);
        }

        /// <summary>
        /// Test the sample methods of the running statistics class. These methods were validated against the values found in the testing of the
        /// Statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Stats()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat = new RunningStatistics(sample1);

            var min = runningstat.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = runningstat.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = runningstat.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = runningstat.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = runningstat.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = runningstat.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = runningstat.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = runningstat.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }

        /// <summary>
        /// Test the population method of the running statistics class. These methods were validated against the values found in the testing of the
        /// Statistics class, or through the direct equation for the value calculated using R.
        /// </summary>
        [TestMethod]
        public void Test_PopStats()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat = new RunningStatistics(sample1);

            var popVariance = runningstat.PopulationVariance;
            Assert.AreEqual(2266.97824826717, popVariance, 1E-10);

            var popSD = runningstat.PopulationStandardDeviation;
            Assert.AreEqual(47.6127950058298, popSD, 1E-10);

            // sqrt(n)* sum((sample1-mean(sample1))**3) / sum((sample1-mean(sample1))**2)**(3/2)  -- R eqm
            var popSkewness = runningstat.PopulationSkewness;
            Assert.AreEqual(0.84172348076203, popSkewness, 1E-10);

            // (sum((sample1-mean(sample1))**4)/n) / (sum((sample1-mean(sample1))**2)/n)**2 - 3  -- R eqn
            var popKurtosis = runningstat.PopulationKurtosis;
            Assert.AreEqual(1.16237367377557, popKurtosis, 1E-10);
        }

        /// <summary>
        /// Test the push methods of the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Push()
        {
            var values = new double[] { 1, 2, 3, 4, 5 };
            var runningstat = new RunningStatistics(values);

            // push new min
            runningstat.Push(-1);
            var newMin = runningstat.Minimum;
            Assert.AreEqual(-1, newMin);

            // push new max
            runningstat.Push(6);
            var newMax = runningstat.Maximum;
            Assert.AreEqual(6, newMax);

            // push new max and min in a list
            runningstat.Push(new double[] { -5, 2.5, 3.33, 12, 20 });
            var newMinList = runningstat.Minimum;
            Assert.AreEqual(-5, newMinList);
            var newMaxList = runningstat.Maximum;
            Assert.AreEqual(20, newMaxList);
        }

        /// <summary>
        /// Test the combine method of the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Combine()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166 };
            var runningstat1 = new RunningStatistics(sample1);

            var sample2 = new double[] { 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat2 = new RunningStatistics(sample2);

            var complete = RunningStatistics.Combine(runningstat1, runningstat2);

            // values from previous tests
            var min = complete.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = complete.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = complete.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = complete.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = complete.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = complete.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = complete.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = complete.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }

        /// <summary>
        /// Test the overload of the + operator for the running statistics class.
        /// </summary>
        [TestMethod]
        public void Test_Add()
        {
            var sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166 };
            var runningstat1 = new RunningStatistics(sample1);

            var sample2 = new double[] { 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
            var runningstat2 = new RunningStatistics(sample2);

            var complete = runningstat1 + runningstat2;

            // values from previous tests
            var min = complete.Minimum;
            Assert.AreEqual(98.8, min, 1E-10);

            var max = complete.Maximum;
            Assert.AreEqual(337.0, max, 1E-10);

            var mean = complete.Mean;
            Assert.AreEqual(191.317391304348, mean, 1E-10);

            var variance = complete.Variance;
            Assert.AreEqual(2300.31616368286, variance, 1E-10);

            var sd = complete.StandardDeviation;
            Assert.AreEqual(47.9616113541118, sd, 1E-10);

            var coeffOfVariation = complete.CoefficientOfVariation;
            Assert.AreEqual(0.250691330396694, coeffOfVariation, 1E-10);

            var skewness = complete.Skewness;
            Assert.AreEqual(0.8605451107461, skewness, 1E-10);

            var kurtosis = complete.Kurtosis;
            Assert.AreEqual(1.3434868130194, kurtosis, 1E-10);
        }
    }
}
