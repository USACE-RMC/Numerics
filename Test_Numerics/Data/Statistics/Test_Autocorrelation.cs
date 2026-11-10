using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Distributions;

namespace Data.Statistics
{
    /// <summary>
    /// Unit testing for the autocorrelation class. These methods were tested against values attained from R's "acf()" method from the "stats" package.
    /// function from the "stats" package.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    /// <b> References: </b>
    ///  R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, 
    ///  Vienna, Austria. ISBN 3-900051-07-0, URL <see href="http://www.R-project.org/."/>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Autocorrelation
    {
        /// <summary>
        /// Test the covariance method with an array of data
        /// </summary>
        [TestMethod()]
        public void Test_Covariance()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var acvf = Autocorrelation.Function(sample, type: Autocorrelation.Type.Covariance);
            var true_acvf = new double[] { 22.105148d, 21.129089d, 20.268354d, 19.478147d, 18.595148d, 17.568448d, 16.893142d, 16.297683d, 15.644932d, 15.187918d, 14.655771d, 13.893153d, 13.017684d, 12.294179d, 11.298568d, 10.43392d, 9.836657d, 9.290689d, 8.714833d, 8.114298d, 7.537443d, 7.01491d, 6.621986d, 5.987118d };
            for (int i = 0; i < true_acvf.Length; i++)
                Assert.AreEqual(acvf[i, 1], true_acvf[i], 0.0001d);
        }

        /// <summary>
        /// The uncentered functions of a zero-mean process, against exact rational values (fractions.Fraction): the
        /// autocovariance Σ x_t x_(t-k) / n, the autocorrelation Σ x_t x_(t-k) / Σ x_t² (numpy agrees), and the partial
        /// autocorrelation from the exact Yule-Walker solutions. The centered default is unchanged.
        /// </summary>
        [TestMethod()]
        public void Test_Uncentered()
        {
            var sample = new[] { 0.3, -1.2, 0.8, 0.5, -0.7, 1.1, -0.4, -0.9, 0.6, 0.2 };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1, 0, 0, 0), sample);
            var cov = new[] { 0.549, -0.254, -0.15, 0.288 };
            var acf = new[] { 1.0, -0.46265938069216755, -0.273224043715847, 0.5245901639344263 };
            var pacf = new[] { -0.46265938069216755, -0.6199886020642927, 0.06813084015755634 };
            foreach (var (c, r, p) in new[]
            {
                (Autocorrelation.Function(sample, 3, Autocorrelation.Type.Covariance, false), Autocorrelation.Function(sample, 3, Autocorrelation.Type.Correlation, false), Autocorrelation.Function(sample, 3, Autocorrelation.Type.Partial, false)),
                (Autocorrelation.Function(ts, 3, Autocorrelation.Type.Covariance, false), Autocorrelation.Function(ts, 3, Autocorrelation.Type.Correlation, false), Autocorrelation.Function(ts, 3, Autocorrelation.Type.Partial, false)),
            })
            {
                for (int k = 0; k <= 3; k++)
                {
                    Assert.AreEqual(k, c![k, 0]);
                    Assert.AreEqual(cov[k], c[k, 1], 1e-15);
                    Assert.AreEqual(acf[k], r![k, 1], 1e-15);
                }
                for (int k = 0; k < 3; k++)
                {
                    Assert.AreEqual(k + 1, p![k, 0]);
                    Assert.AreEqual(pacf[k], p[k, 1], 1e-14);
                }
            }

            var centered = Autocorrelation.Function(sample, 3, Autocorrelation.Type.Correlation, true)!;
            var byDefault = Autocorrelation.Function(sample, 3)!;
            for (int k = 0; k <= 3; k++)
                Assert.AreEqual(byDefault[k, 1], centered[k, 1], 0d);
        }

        /// <summary>
        /// The autocorrelation about zero to lag 20 of a 200-value series, against the correctly rounded exact
        /// ratios S_k / S_0 of the double inputs (Python fractions.Fraction, Python 3 with numpy 2.5.4;
        /// statsmodels 0.15.0 acovf(demean=False, adjusted=False, fft=False) / acovf[0] is within 3 ulp of the
        /// same values). The autocovariance about zero keeps the bits of the centered loop evaluated about zero,
        /// the correlation is each lag sum divided once by S_0, and the list, array and time-series inputs agree
        /// bitwise.
        /// </summary>
        [TestMethod()]
        public void Test_Uncentered_ExactRatios()
        {
            const int n = 200, lagMax = 20;
            var x = new double[n];
            for (int t = 0; t < n; t++)
                x[t] = ((t * 37) % 101 - 50) / 10.0 + 0.3 * ((t * 53) % 17 - 8) / 7.0;
            var exact = new[] { 1.0, -0.3787176557080698, -0.1830780409782735, 0.46746662717808185, -0.4772807834319895, 0.1338066896813549, 0.055205326588923714, -0.4526233763951284, 0.5690599179503525, -0.23206700297672522, -0.32468748660343033, 0.7943325964057234, -0.400576463671353, -0.10203088730921354, 0.32548062290471713, -0.4542297178755157, 0.23414858498380062, -0.034700514748903295, -0.4042519183294049, 0.6732011249496646, -0.2785445052303177 };

            var sums = new double[lagMax + 1];
            for (int lag = 0; lag <= lagMax; lag++)
                for (int t = lag; t < n; t++)
                    sums[lag] += (x[t] - 0d) * (x[t - lag] - 0d);

            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1, 0, 0, 0), x);
            var inputs = new[]
            {
                (Autocorrelation.Function(x, lagMax, Autocorrelation.Type.Correlation, false)!, Autocorrelation.Function(x, lagMax, Autocorrelation.Type.Covariance, false)!),
                (Autocorrelation.Function(new List<double>(x), lagMax, Autocorrelation.Type.Correlation, false)!, Autocorrelation.Function(new List<double>(x), lagMax, Autocorrelation.Type.Covariance, false)!),
                (Autocorrelation.Function(ts, lagMax, Autocorrelation.Type.Correlation, false)!, Autocorrelation.Function(ts, lagMax, Autocorrelation.Type.Covariance, false)!),
            };
            foreach (var (acf, acvf) in inputs)
            {
                Assert.AreEqual(lagMax + 1, acf.GetLength(0));
                for (int k = 0; k <= lagMax; k++)
                {
                    Assert.AreEqual(k, acf[k, 0]);
                    Assert.AreEqual(k, acvf[k, 0]);
                    Assert.AreEqual(exact[k], acf[k, 1], 2e-15);
                    Assert.AreEqual(sums[k] / sums[0], acf[k, 1], 0d);
                    Assert.AreEqual(sums[k] / n, acvf[k, 1], 0d);
                }
            }

            Assert.IsNull(Autocorrelation.Function(new double[10], 3, Autocorrelation.Type.Correlation, false));
            Assert.IsNull(Autocorrelation.Function(new[] { 1.0 }, 1, Autocorrelation.Type.Correlation, false));
        }

        /// <summary>
        /// The centered default keeps the bits of the mean-removed loop S_k / n divided by S_0 / n.
        /// </summary>
        [TestMethod()]
        public void Test_Centered_Bitwise()
        {
            const int n = 200, lagMax = 20;
            var x = new double[n];
            for (int t = 0; t < n; t++)
                x[t] = ((t * 37) % 101 - 50) / 10.0 + 0.3 * ((t * 53) % 17 - 8) / 7.0 + 140.0;
            double mean = Numerics.Data.Statistics.Statistics.Mean(x);
            var cov = new double[lagMax + 1];
            for (int lag = 0; lag <= lagMax; lag++)
            {
                double c = 0d;
                for (int t = lag; t < n; t++)
                    c += (x[t] - mean) * (x[t - lag] - mean);
                cov[lag] = c / n;
            }
            var acf = Autocorrelation.Function(x, lagMax)!;
            var acvf = Autocorrelation.Function(x, lagMax, Autocorrelation.Type.Covariance)!;
            for (int k = 0; k <= lagMax; k++)
            {
                Assert.AreEqual(cov[k], acvf[k, 1], 0d);
                Assert.AreEqual(cov[k] / cov[0], acf[k, 1], 0d);
            }
        }

        /// <summary>
        /// Test the covariance method with TimeSeries data
        /// </summary>
        [TestMethod()]
        public void Test_Covariance_TimeSeries()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1, 0, 0, 0), sample);
            var acvf = Autocorrelation.Function(ts, type: Autocorrelation.Type.Covariance);
            var true_acvf = new double[] { 22.105148d, 21.129089d, 20.268354d, 19.478147d, 18.595148d, 17.568448d, 16.893142d, 16.297683d, 15.644932d, 15.187918d, 14.655771d, 13.893153d, 13.017684d, 12.294179d, 11.298568d, 10.43392d, 9.836657d, 9.290689d, 8.714833d, 8.114298d, 7.537443d, 7.01491d, 6.621986d, 5.987118d };
            for (int i = 0; i < true_acvf.Length; i++)
                Assert.AreEqual(acvf[i, 1], true_acvf[i], 0.0001d);
        }

        /// <summary>
        /// Test the autocorrelation method with an array of data
        /// </summary>
        [TestMethod()]
        public void Test_Correlation()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var acf = Autocorrelation.Function(sample, type: Autocorrelation.Type.Correlation);
            var true_acf = new double[] { 1d, 0.9558447d, 0.9169065d, 0.8811588d, 0.8412135d, 0.7947673d, 0.7642176d, 0.73728d, 0.7077506d, 0.6870761d, 0.6630026d, 0.6285031d, 0.5888983d, 0.5561681d, 0.5111283d, 0.4720131d, 0.4449939d, 0.4202953d, 0.3942445d, 0.3670773d, 0.3409814d, 0.3173428d, 0.2995676d, 0.2708472d };
            for (int i = 0; i < true_acf.Length; i++)
                Assert.AreEqual(acf[i, 1], true_acf[i], 0.0001d);
            // Test FFT Autocorrelation
            var acf2 = Fourier.Autocorrelation(sample);
            for (int i = 0; i < true_acf.Length; i++)
                Assert.AreEqual(acf2[i, 1], true_acf[i], 0.0001d);
        }

        /// <summary>
        /// Test the autocorrelation method with TimeSeries data
        /// </summary>
        [TestMethod()]
        public void Test_Correlation_TimeSeries()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1, 0, 0, 0), sample);
            var acf = Autocorrelation.Function(ts, type: Autocorrelation.Type.Correlation);
            var true_acf = new double[] { 1d, 0.9558447d, 0.9169065d, 0.8811588d, 0.8412135d, 0.7947673d, 0.7642176d, 0.73728d, 0.7077506d, 0.6870761d, 0.6630026d, 0.6285031d, 0.5888983d, 0.5561681d, 0.5111283d, 0.4720131d, 0.4449939d, 0.4202953d, 0.3942445d, 0.3670773d, 0.3409814d, 0.3173428d, 0.2995676d, 0.2708472d };
            for (int i = 0; i < true_acf.Length; i++)
                Assert.AreEqual(acf[i, 1], true_acf[i], 0.0001d);
        }

        /// <summary>
        /// Test the partial autocorrelation method with an array of data
        /// </summary>
        [TestMethod()]
        public void Test_Partial()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var pacf = Autocorrelation.Function(sample, type: Autocorrelation.Type.Partial);
            var true_pacf = new double[] { 0.955844726d, 0.037833688d, 0.020103563d, -0.063283074d, -0.101338145d, 0.148214184d, 0.042511772d, -0.023952331d, 0.07833128d, -0.075887833d, -0.117610117d, -0.085556279d, 0.03336607d, -0.118322713d, 0.051528628d, 0.089399743d, 0.0019337d, -0.005678822d, -0.104080153d, -0.032036409d, 0.083778031d, 0.078051358d, -0.120023159d };
            for (int i = 0; i < true_pacf.Length; i++)
                Assert.AreEqual(pacf[i, 1], true_pacf[i], 0.0001d);
        }

        /// <summary>
        /// Test the partial autocorrelation method with TimeSeries data
        /// </summary>
        [TestMethod()]
        public void Test_Partial_TimeSeries()
        {
            var sample = new[] { 142.25d, 141.23d, 141.33d, 140.82d, 141.31d, 140.58d, 141.58d, 142.15d, 143.07d, 142.85d, 143.17d, 142.54d, 143.07d, 142.26d, 142.97d, 143.86d, 142.57d, 142.19d, 142.35d, 142.63d, 144.15d, 144.73d, 144.7d, 144.97d, 145.12d, 144.78d, 145.06d, 143.94d, 143.77d, 144.8d, 145.67d, 145.44d, 145.56d, 145.61d, 146.05d, 145.74d, 145.83d, 143.88d, 140.39d, 139.34d, 140.05d, 137.93d, 138.78d, 139.59d, 140.54d, 141.31d, 140.42d, 140.18d, 138.43d, 138.97d, 139.31d, 139.26d, 140.08d, 141.1d, 143.48d, 143.28d, 143.5d, 143.12d, 142.14d, 142.54d, 142.24d, 142.16d, 142.97d, 143.69d, 143.67d, 144.65d, 144.33d, 144.82d, 143.74d, 144.9d, 145.83d, 146.97d, 146.6d, 146.55d, 148.22d, 148.37d, 148.23d, 148.73d, 149.49d, 149.09d, 149.64d, 148.42d, 148.9d, 149.97d, 150.75d, 150.88d, 150.58d, 150.64d, 150.73d, 149.75d, 150.86d, 150.7d, 150.8d, 151.38d, 152.01d, 152.58d, 152.7d, 152.95d, 152.53d, 151.5d, 151.94d, 151.46d, 153.67d, 153.88d, 153.54d, 153.74d, 152.86d, 151.56d, 149.58d, 150.93d, 150.67d, 150.5d, 152.06d, 153.14d, 153.38d, 152.55d, 153.58d, 151.08d, 151.52d, 150.24d, 150.21d, 148.13d, 150.38d, 150.9d, 150.87d, 152.18d, 152.4d, 152.38d, 153.16d, 152.29d, 150.75d, 152.37d, 154.57d, 154.99d, 154.93d, 154.23d, 155.2d, 154.89d, 154.18d, 153.12d, 152.02d, 150.19d, 148.21d, 145.93d, 148.33d, 145.18d, 146.76d, 147.28d, 144.21d, 145.94d, 148.41d, 147.43d, 144.39d, 146.5d, 145.7d, 142.72d, 139.79d, 145.5d, 145.17d, 144.6d, 146.01d, 147.34d, 146.48d, 147.85d, 146.16d, 144.37d, 145.45d, 147.65d, 147.45d, 148.2d, 147.95d, 146.48d, 146.52d, 146.24d, 147.29d, 148.55d, 147.96d, 148.31d, 148.83d, 153.41d, 153.34d, 152.71d, 152.42d, 150.81d, 152.25d, 152.91d, 152.85d, 152.6d, 154.61d, 153.81d, 154.11d, 155.03d, 155.39d, 155.6d, 156.04d, 156.93d, 155.46d, 156.27d, 154.41d, 154.98d };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1, 0, 0, 0), sample);
            var pacf = Autocorrelation.Function(ts, type: Autocorrelation.Type.Partial);
            var true_pacf = new double[] { 0.955844726d, 0.037833688d, 0.020103563d, -0.063283074d, -0.101338145d, 0.148214184d, 0.042511772d, -0.023952331d, 0.07833128d, -0.075887833d, -0.117610117d, -0.085556279d, 0.03336607d, -0.118322713d, 0.051528628d, 0.089399743d, 0.0019337d, -0.005678822d, -0.104080153d, -0.032036409d, 0.083778031d, 0.078051358d, -0.120023159d };
            for (int i = 0; i < true_pacf.Length; i++)
                Assert.AreEqual(pacf[i, 1], true_pacf[i], 0.0001d);
        }

        /// <summary>
        /// Test the method for attaining a confidence interval for the correlation function
        /// </summary>
        [TestMethod]
        public void Test_CorrelationConfidenceInterval()
        {
            var sizes = new int[] { 20, 30, 40, 50, 50, 60 };
            var intervals = new double[] { 0.90, 0.95, 0.99 };

            for(int i = 0; i < sizes.Length;i++)
            {
                for(int j = 0; j < intervals.Length; j++)
                {
                    var cci = Autocorrelation.CorrelationConfidenceInterval(sizes[i], intervals[j]);
                    var high = Normal.StandardZ(0.5*(1 - intervals[j])) / Math.Sqrt(sizes[i]);
                    var low = -high;

                    Assert.AreEqual(high, cci[0], 1E-6);
                    Assert.AreEqual(low, cci[1], 1E-6);
                }
            }

        }

        /// <summary>
        /// A maximum lag at or beyond the number of observations gives one row per lag, with zero for every lag of n or
        /// more, centered and uncentered, for a list and a time series; the lower lags match the call with lagMax = n - 1.
        /// </summary>
        [TestMethod]
        public void Test_LagMaxBeyondSampleSize()
        {
            var d = new double[] { 1, -2, 3, -1, 2 };
            var ts = new TimeSeries(TimeInterval.OneDay, new DateTime(2020, 1, 1), d);
            foreach (bool centered in new[] { true, false })
            {
                foreach (var type in new[] { Autocorrelation.Type.Covariance, Autocorrelation.Type.Correlation })
                {
                    foreach (var acf in new[] { Autocorrelation.Function(d, 10, type, centered), Autocorrelation.Function(ts, 10, type, centered) })
                    {
                        var reference = Autocorrelation.Function(d, 4, type, centered);
                        Assert.AreEqual(11, acf!.GetLength(0));
                        for (int lag = 0; lag <= 10; lag++)
                        {
                            Assert.AreEqual(lag, acf[lag, 0]);
                            Assert.AreEqual(lag <= 4 ? reference![lag, 1] : 0d, acf[lag, 1], 1E-15, $"{type} centered={centered} lag {lag}");
                        }
                    }
                }
                var pacf = Autocorrelation.Function(d, 10, Autocorrelation.Type.Partial, centered);
                var pacfReference = Autocorrelation.Function(d, 4, Autocorrelation.Type.Partial, centered);
                Assert.AreEqual(10, pacf!.GetLength(0));
                for (int lag = 0; lag < 4; lag++)
                    Assert.AreEqual(pacfReference![lag, 1], pacf[lag, 1], 1E-15, $"Partial centered={centered} row {lag}");
            }
        }
    }
}
