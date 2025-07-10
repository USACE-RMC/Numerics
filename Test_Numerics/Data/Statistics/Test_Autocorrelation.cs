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
    }
}
