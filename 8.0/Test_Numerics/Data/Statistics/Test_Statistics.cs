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
using System.Diagnostics;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the Statistics class. These methods were test against various R methods mostly from the "base" package. The specific functions and any other packages
    /// used are documented in the individual tests.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// R Core Team (2024). _R: A Language and Environment for Statistical Computing_.R Foundation for Statistical Computing, Vienna,
    /// Austria. <see href="https://www.R-project.org/"/>
    /// </remarks>
    [TestClass]
    public class Test_Statistics
    {
        /// <summary>
        /// Datasets to test the methods on.
        /// </summary>
        private double[] _sample1 = new double[] { 122, 244, 214, 173, 229, 156, 212, 263, 146, 183, 161, 205, 135, 331, 225, 174, 98.8, 149, 238, 262, 132, 235, 216, 240, 230, 192, 195, 172, 173, 172, 153, 142, 317, 161, 201, 204, 194, 164, 183, 161, 167, 179, 185, 117, 192, 337, 125, 166, 99.1, 202, 230, 158, 262, 154, 164, 182, 164, 183, 171, 250, 184, 205, 237, 177, 239, 187, 180, 173, 174 };
        private double[] _sample2 = new double[] { 279d, 105d, 171d, 171d, 129d, 127d, 194d, 234d, 251d, 152d, 207d, 205d, 183d, 137d, 148d, 189d, 182d, 236d, 148d, 150d, 207d, 252d, 237d, 209d, 225d, 137d, 207d, 129d, 148d, 192d, 95d, 231d, 255d, 220d, 205d, 163d, 265d, 190d, 226d, 123d, 108d, 145d, 197d, 233d, 133d, 177d, 211d, 180d, 200d, 197d, 142d, 166d, 251d, 254d, 226d, 197d, 250d, 194d, 190d, 181d, 290d, 185d, 123d, 208d, 238d, 179d, 189d, 225d, 236d };
        private double[] _jackKnifeSample = new double[] { 3.35602585719312d, 3.07554696139253d, 3.04921802267018d, 3.07554696139253d, 3.13033376849501d, 3.13033376849501d, 3.31597034545692d, 3.15533603746506d, 3.17897694729317d, 3.29446622616159d, 3.37657695705651d, 3.72015930340596d, 3.25042000230889d, 3.25042000230889d, 3.10380372095596d, 3.45178643552429d, 3.15533603746506d, 3.37657695705651d, 3.45178643552429d, 3.29446622616159d };

        /// <summary>
        /// Test the Minimum method against R's "min()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_Minimum()
        {
            double val = Numerics.Data.Statistics.Statistics.Minimum(_sample1);
            double trueVal = 98.8d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the Maximum method against R's "max()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_Maximum()
        {
            double val = Numerics.Data.Statistics.Statistics.Maximum(_sample1);
            double trueVal = 337.0d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the Mean method against R's "mean()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            double val = Numerics.Data.Statistics.Statistics.Mean(_sample1);
            double trueVal = 191.317391304348d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the ParallelMean method with the direct equation. Should also be the same as the arithmetic mean in this case.
        /// </summary>
        [TestMethod]
        public void Test_ParallelMean()
        {
            // basic equation for parallel mean
            var parallel = _sample1.AsParallel();
            var valid = parallel.Sum() / parallel.Count();

            double test = Numerics.Data.Statistics.Statistics.ParallelMean(_sample1);
            double regMean = Numerics.Data.Statistics.Statistics.Mean(_sample1);
            Assert.AreEqual(valid, test, 1E-6);
            Assert.AreEqual(test, regMean, 1E-6);
        }

        /// <summary>
        /// Test the GeometricMean method against R's "geometric.mean()" method from the "psych" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// William Revelle (2024). psych: Procedures for Psychological, Psychometric, and Personality Research. Northwestern University,
        /// Evanston, Illinois. R package version 2.4.6, <see href="https://CRAN.R-project.org/package=psych"/>
        /// </remarks>
        [TestMethod]
        public void Test_GeometricMean()
        {
            double val = Numerics.Data.Statistics.Statistics.GeometricMean(_sample1);
            double trueVal = 185.685629284673d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the HarmonicMean method against R's "harmonic.mean()" method from the "psych" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// William Revelle (2024). psych: Procedures for Psychological, Psychometric, and Personality Research. Northwestern University,
        /// Evanston, Illinois. R package version 2.4.6, <see href="https://CRAN.R-project.org/package=psych"/>
        /// </remarks>
        [TestMethod]
        public void Test_HarmonicMean()
        {
            double val = Numerics.Data.Statistics.Statistics.HarmonicMean(_sample1);
            double trueVal = 180.183870381546d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the Variance method against R's "var()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_Variance()
        {
            double val = Numerics.Data.Statistics.Statistics.Variance(_sample1);
            double trueVal = 2300.31616368286d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the PopulationVariance method against R's "var()" method buy multiplying the sample variance calculated from "var()" by
        /// (n-1) / n as per the definition of population variance.
        /// </summary>
        [TestMethod]
        public void Test_PopulationVariance()
        {
            double val = Numerics.Data.Statistics.Statistics.PopulationVariance(_sample1);
            double trueVal = 2266.97824826717d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the StandardDeviation method against R's "sd()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_StandardDeviation()
        {
            double val = Numerics.Data.Statistics.Statistics.StandardDeviation(_sample1);
            double trueVal = 47.9616113541118d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the PopulationStandardDeviation method against R by taking the square root of the population variance, calculated by multiplying the 
        /// "var()" method by (n-1) / n 
        /// </summary>
        [TestMethod]
        public void Test_PopulationStandardDeviation()
        {
            double val = Numerics.Data.Statistics.Statistics.PopulationStandardDeviation(_sample1);
            double trueVal = 47.6127950058298d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the MeanVariance method against values from R's "mean()" and "var()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_MeanVariance()
        {
            var test = Numerics.Data.Statistics.Statistics.MeanVariance(_sample1);
            double trueMean = 191.317391304348d;
            double trueVar = 2300.31616368286d;

            Assert.AreEqual(trueMean, test.Item1, 1E-10);
            Assert.AreEqual(trueVar, test.Item2, 1E-10);
        }

        /// <summary>
        /// Test the MeanStandardDeviation method against R's "mean()" and "sd()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_MeanStandardDeviation()
        {
            var test = Numerics.Data.Statistics.Statistics.MeanStandardDeviation(_sample1);
            double trueMean = 191.317391304348d;
            double trueSD = 47.9616113541118d;

            Assert.AreEqual(trueMean, test.Item1, 1E-10);
            Assert.AreEqual(trueSD, test.Item2, 1E-10);
        }

        /// <summary>
        /// Test the CoefficientOfVariation method against R's "cv()" method from the "goeveg" package.
        /// </summary>
        /// <remarks>
        /// <b> Refernces: </b>
        /// von Lampe F, Schellenberg J (2024). goeveg: Functions for Community Data and Ordinations. 
        /// R package version 0.7.5, <see href="https://CRAN.R-project.org/package=goeveg"/>
        /// </remarks>
        [TestMethod]
        public void Test_CoefficientOfCVariation()
        {
            double test = Numerics.Data.Statistics.Statistics.CoefficientOfVariation(_sample1);
            double valid = 0.250691330396694;
            Assert.AreEqual(valid, test, 1E-10);
        }

        /// <summary>
        /// Test the Skewness method against R's "skewness()" method from the "EnvStats" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Millard SP (2013). EnvStats: An R Package for Environmental Statistics. Springer, New York.
        /// ISBN 978-1-4614-8455-4, <see href="https://www.springer.com"/>
        /// </remarks>
        [TestMethod]
        public void Test_Skewness()
        {
            double val = Numerics.Data.Statistics.Statistics.Skewness(_sample1);
            double trueVal = 0.8605451107461d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the JackKnifeStandardError method against R's "jackknife()" method from the "bootstrap" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// original S, StatLib f, Leisch. bRTRpbF (2019). bootstrap: Functions for the Book "An Introduction to the Bootstrap". 
        /// R package version 2019.6, <see href="https://CRAN.R-project.org/package=bootstrap"/>
        /// </remarks>
        [TestMethod]
        public void Test_JackKnifeStandardError()
        {
            double val = Numerics.Data.Statistics.Statistics.JackKnifeStandardError(_jackKnifeSample, Numerics.Data.Statistics.Statistics.Mean);
            double trueVal = 0.0372182162668589d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the JackKnifeSample method against R's "jackknife()" method from the "bootstrap" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// original S, StatLib f, Leisch. bRTRpbF (2019). bootstrap: Functions for the Book "An Introduction to the Bootstrap". 
        /// R package version 2019.6, <see href="https://CRAN.R-project.org/package=bootstrap"/>
        /// </remarks>
        [TestMethod]
        public void Test_JackKnifeSample()
        {
            var test = Numerics.Data.Statistics.Statistics.JackKnifeSample(_jackKnifeSample, Numerics.Data.Statistics.Statistics.Mean);
            var valid = new double[] { 3.254582, 3.269344, 3.270730, 3.269344, 3.266461, 3.266461, 3.256690, 3.265145, 3.263901, 3.257822, 3.253501, 3.235417, 3.260140, 3.260140, 3.267857, 3.249542, 3.265145, 3.253501, 3.249542, 3.257822 };
            for (int i = 0; i < valid.Length; i++)
            {
                Assert.AreEqual(valid[i], test[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the Kurtosis method with R's "kurtosis()" method from the "PerformanceAnalytics" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// Peterson BG, Carl P (2020). PerformanceAnalytics: Econometric Tools for Performance and Risk Analysis. 
        /// R package version 2.0.4, <see href="https://CRAN.R-project.org/package=PerformanceAnalytics"/>
        /// </remarks>
        [TestMethod]
        public void Test_Kurtosis()
        {
            double val = Numerics.Data.Statistics.Statistics.Kurtosis(_sample1);
            double trueVal = 1.3434868130194d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the Covariance method against R's "cov()" method from the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_Covariance()
        {
            double val = Numerics.Data.Statistics.Statistics.Covariance(_sample1, _sample2);
            double trueVal = -253.54405370844d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test thr PopulationCovariance method against R's "cov" method by multiplying the sample covariance by (n-1) / n.
        /// </summary>
        [TestMethod]
        public void Test_PopulationCovariance()
        {
            double val = Numerics.Data.Statistics.Statistics.PopulationCovariance(_sample1, _sample2);
            double trueVal = -249.869502205419d;
            Assert.AreEqual(val, trueVal, 1E-10);
        }

        /// <summary>
        /// Test the ProductMoments methods that returns {Mean, Standard Deviation, Skew, and Kurtosis} against the validation values 
        /// attained in the previous tests of each individual statistic.
        /// </summary>
        [TestMethod]
        public void Test_ComputeProductMoments()
        {
            var vals = Numerics.Data.Statistics.Statistics.ProductMoments(_sample1);
            double trueVal1 = 191.317391304348d;
            double trueVal2 = 47.9616113541118d;
            double trueVal3 = 0.8605451107461d;
            double trueVal4 = 1.3434868130194d;

            Assert.AreEqual(vals[0], trueVal1, 1E-10);
            Assert.AreEqual(vals[1], trueVal2, 1E-10);
            Assert.AreEqual(vals[2], trueVal3, 1E-10);
            Assert.AreEqual(vals[3], trueVal4, 1E-10);
        }

        /// <summary>
        /// Test the LinearMoments method against the "samlmu()" method of the "lmom" package.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// J. R. M. Hosking (2023). L-Moments. R package, version 3.0. URL: <see href="https://CRAN.R-project.org/package=lmom."/>
        /// </remarks>
        [TestMethod]
        public void Test_ComputeLinearMoments()
        {
            var data = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var lmoms = Numerics.Data.Statistics.Statistics.LinearMoments(data);
            double trueVal1 = 1648.8064516d;
            double trueVal2 = 138.2365591d;
            double trueVal3 = 0.1033903d;
            double trueVal4 = 0.1940943d;

            Assert.AreEqual(lmoms[0], trueVal1, 1E-7);
            Assert.AreEqual(lmoms[1], trueVal2, 1E-7);
            Assert.AreEqual(lmoms[2], trueVal3, 1E-7);
            Assert.AreEqual(lmoms[3], trueVal4, 1E-7);
        }

        /// <summary>
        /// Test the Percentile method against R's "quantile()" method from the "stats" package.
        /// </summary>
        [TestMethod]
        public void Test_Percentiles()
        {
            double val0 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0);
            double trueVal0 = 98.8d;
            Assert.AreEqual(val0, trueVal0, 1E-2);

            double val1 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.01d);
            double trueVal1 = 99.004d;
            Assert.AreEqual(val1, trueVal1, 1E-2);

            double val2 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.05d);
            double trueVal2 = 123.2d;
            Assert.AreEqual(val2, trueVal2, 1E-2);

            double val3 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.25d);
            double trueVal3 = 164.0d;
            Assert.AreEqual(val3, trueVal3, 1E-2);

            double val4 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.5d);
            double trueVal4 = 183.0d;
            Assert.AreEqual(val4, trueVal4, 1E-2);

            double val5 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.75d);
            double trueVal5 = 216.0d;
            Assert.AreEqual(val5, trueVal5, 1E-2);

            double val6 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.95d);
            double trueVal6 = 262.6d;
            Assert.AreEqual(val6, trueVal6, 1E-2);

            double val7 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 0.99d);
            double trueVal7 = 332.92d;
            Assert.AreEqual(val7, trueVal7, 1E-2);

            double val8 = Numerics.Data.Statistics.Statistics.Percentile(_sample1, 1);
            double trueVal8 = 337d;
            Assert.AreEqual(val8, trueVal8, 1E-2);
        }

        /// <summary>
        /// Test the FiveNumberSummary methods that returns {min, 25th-percentile, 50th-percentile, 75th-percentile, max} against the validation values 
        /// attained in the previous tests of each individual statistic.
        /// </summary>
        [TestMethod]
        public void Test_FiveNumberSummary()
        {
            var summary = Numerics.Data.Statistics.Statistics.FiveNumberSummary(_sample1);
            double trueMin = 98.8d;
            double p25 = 164.0d;
            double p50 = 183.0d;
            double p75 = 216.0d;
            double trueMax = 337.0d;

            Assert.AreEqual(trueMin, summary[0], 1E-2);
            Assert.AreEqual(p25, summary[1], 1E-2);
            Assert.AreEqual(p50, summary[2], 1E-2);
            Assert.AreEqual(p75, summary[3], 1E-2);
            Assert.AreEqual(trueMax, summary[4], 1E-2);
        }

        /// <summary>
        /// Test the SevenNumberSummary methods that returns {min, 5th percentile, 25th-percentile, 50th-percentile, 75th-percentile, 95th-percentile, max} against the 
        /// validation values attained in the previous tests of each individual statistic.
        /// </summary>
        [TestMethod]
        public void Test_SevenNumberSummary()
        {
            var summary = Numerics.Data.Statistics.Statistics.SevenNumberSummary(_sample1);
            double trueMin = 98.8d;
            double p5 = 123.2d;
            double p25 = 164.0d;
            double p50 = 183.0d;
            double p75 = 216.0d;
            double p95 = 262.6d;
            double trueMax = 337.0d;

            Assert.AreEqual(trueMin, summary[0], 1E-2);
            Assert.AreEqual(p5, summary[1], 1E-2);
            Assert.AreEqual(p25, summary[2], 1E-2);
            Assert.AreEqual(p50, summary[3], 1E-2);
            Assert.AreEqual(p75, summary[4], 1E-2);
            Assert.AreEqual(p95, summary[5], 1E-2);
            Assert.AreEqual(trueMax, summary[6], 1E-2);
        }

        /// <summary>
        /// Test the RanksInplace() method that returns the rank of each entry against R's "rank()" method of the "base" package.
        /// </summary>
        [TestMethod]
        public void Test_RanksInPlace()
        {
            var test = Numerics.Data.Statistics.Statistics.RanksInPlace(_sample1);
            var valid = new double[] { 4.0, 62.0, 51.0, 27.0, 54.0, 13.0, 50.0, 66.0, 9.0, 36.0, 16.0, 48.5, 7.0, 68.0, 53.0, 29.5, 1.0, 10.0, 59.0, 64.5, 6.0, 57.0, 52.0, 61.0, 55.5, 41.5, 44.0, 24.5, 27.0, 24.5, 11.0, 8.0, 67.0, 16.0, 45.0, 47.0, 43.0, 19.0, 36.0, 16.0, 22.0, 32.0, 39.0, 3.0, 41.5, 69.0, 5.0, 21.0, 2.0, 46.0, 55.5, 14.0, 64.5, 12.0, 19.0, 34.0, 19.0, 36.0, 23.0, 63.0, 38.0, 48.5, 58.0, 31.0, 60.0, 40.0, 33.0, 27.0, 29.5 };

            for (int i = 0; i < valid.Length; i++)
            {
                Assert.AreEqual(valid[i], test[i]);
            }
        }

        /// <summary>
        /// Test the RanksInPlace() list of ties. 
        /// </summary>
        [TestMethod]
        public void Test_RanksInPlace_Ties()
        {
            var data = _sample1.ToArray();
            Array.Sort(data);
            Numerics.Data.Statistics.Statistics.RanksInPlace(data, out var ties);
            var valid = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };
            for (int i = 0; i < valid.Length; i++)
            {
                Assert.AreEqual(valid[i], ties[i]);
            }
        }

        /// <summary>
        /// Test the Entropy function against the value derived from the direct function: sum of p*ln(p)
        /// </summary>
        [TestMethod]
        public void Test_Entropy()
        {
            var norm = new Normal(100, 15);
            var data = new double[30];
            for (int i = 1; i <= 30; i++)
            {
                data[i - 1] = norm.InverseCDF((double)i / 31);
            }

            Func<double, double> func = (double x) =>
            {
                return norm.PDF(x);
            };

            var test = Numerics.Data.Statistics.Statistics.Entropy(data, func);
            double valid = 2.252937012209;
            Assert.AreEqual(valid, test, 1E-10);
        }
    }
}

