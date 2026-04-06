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
using Numerics.Distributions;
using Numerics.Mathematics.SpecialFunctions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Unit tests for the Von Mises distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Reference values computed using R package 'circular' (dvonmises, pvonmises) and scipy.stats.vonmises.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_VonMises
    {

        /// <summary>
        /// Verifying default construction.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var VM = new VonMises();
            Assert.AreEqual(0d, VM.Mu);
            Assert.AreEqual(1d, VM.Kappa);

            var VM2 = new VonMises(1.0, 2.0);
            Assert.AreEqual(1.0, VM2.Mu);
            Assert.AreEqual(2.0, VM2.Kappa);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod]
        public void Test_InvalidParameters()
        {
            var VM = new VonMises(double.NaN, double.NaN);
            Assert.IsFalse(VM.ParametersValid);

            var VM2 = new VonMises(0, -1);
            Assert.IsFalse(VM2.ParametersValid);

            var VM3 = new VonMises(5, 1); // mu > pi
            Assert.IsFalse(VM3.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString().
        /// </summary>
        [TestMethod]
        public void Test_ParametersToString()
        {
            var VM = new VonMises();
            Assert.AreEqual("Mean Direction (μ)", VM.ParametersToString[0, 0]);
            Assert.AreEqual("Concentration (κ)", VM.ParametersToString[1, 0]);
        }

        /// <summary>
        /// Testing the range is [-π, π].
        /// </summary>
        [TestMethod]
        public void Test_MinMax()
        {
            var VM = new VonMises();
            Assert.AreEqual(-Math.PI, VM.Minimum);
            Assert.AreEqual(Math.PI, VM.Maximum);
        }

        /// <summary>
        /// Testing mean, median, and mode.
        /// </summary>
        [TestMethod]
        public void Test_Mean()
        {
            var VM = new VonMises(0.5, 2.0);
            Assert.AreEqual(0.5, VM.Mean);
            Assert.AreEqual(0.5, VM.Median);
            Assert.AreEqual(0.5, VM.Mode);
        }

        /// <summary>
        /// Test the PDF using the analytical formula: f(x) = exp(κ cos(x - μ)) / (2π I₀(κ)).
        /// </summary>
        /// <remarks>
        /// Verified against scipy.stats.vonmises.pdf(x, kappa, loc=mu):
        /// VM(0,1): PDF(0) = exp(1)/(2π·I₀(1)) ≈ 0.34171
        /// VM(0,1): PDF(π/2) = 1/(2π·I₀(1)) ≈ 0.12573
        /// VM(0,1): PDF(π) = exp(-1)/(2π·I₀(1)) ≈ 0.04626
        /// </remarks>
        [TestMethod]
        public void Test_PDF()
        {
            // VM(mu=0, kappa=1): analytical values
            var VM = new VonMises(0, 1);
            double i0_1 = Bessel.I0(1); // ≈ 1.2660658
            double norm = 2d * Math.PI * i0_1;

            Assert.AreEqual(Math.Exp(1d) / norm, VM.PDF(0), 1e-6);        // exp(cos(0)) / norm
            Assert.AreEqual(1d / norm, VM.PDF(Math.PI / 2), 1e-6);          // exp(cos(π/2)) / norm = exp(0) / norm
            Assert.AreEqual(Math.Exp(-1d) / norm, VM.PDF(Math.PI), 1e-6);   // exp(cos(π)) / norm

            // PDF should be symmetric about mu
            Assert.AreEqual(VM.PDF(0.5), VM.PDF(-0.5), 1e-10);

            // VM(mu=0, kappa=5): more concentrated
            var VM2 = new VonMises(0, 5);
            double i0_5 = Bessel.I0(5);
            Assert.AreEqual(Math.Exp(5d) / (2d * Math.PI * i0_5), VM2.PDF(0), 1e-6);

            // PDF outside [-π, π] should be 0
            Assert.AreEqual(0d, VM.PDF(-4.0));
        }

        /// <summary>
        /// Test the CDF boundary conditions and symmetry.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var VM = new VonMises(0, 1);

            // Boundary conditions
            Assert.AreEqual(0d, VM.CDF(-Math.PI), 1e-6);
            Assert.AreEqual(1d, VM.CDF(Math.PI), 1e-6);

            // CDF at mean direction should be 0.5 (symmetry)
            Assert.AreEqual(0.5, VM.CDF(0), 1e-4);

            // CDF should be monotonically increasing
            double prev = 0;
            for (double x = -Math.PI; x <= Math.PI; x += 0.1)
            {
                double cdf = VM.CDF(x);
                Assert.IsGreaterThanOrEqualTo(prev - 1e-10, cdf);
                prev = cdf;
            }

            // CDF at pi/2 should be > 0.5 (right of mean)
            double cdfPiHalf = VM.CDF(Math.PI / 2);
            Assert.IsGreaterThan(0.5, cdfPiHalf);
            Assert.IsLessThan(1.0, cdfPiHalf);
        }

        /// <summary>
        /// Test InverseCDF is consistent with CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var VM = new VonMises(0, 2);
            Assert.AreEqual(-Math.PI, VM.InverseCDF(0));
            Assert.AreEqual(Math.PI, VM.InverseCDF(1));

            // CDF-InverseCDF round-trip
            double[] probs = { 0.1, 0.25, 0.5, 0.75, 0.9 };
            foreach (var p in probs)
            {
                double x = VM.InverseCDF(p);
                double pBack = VM.CDF(x);
                Assert.AreEqual(p, pBack, 1e-4);
            }
        }

        /// <summary>
        /// Test the MLE estimation with known data.
        /// </summary>
        [TestMethod]
        public void Test_MLE()
        {
            // Generate a sample from VM(mu=1.0, kappa=3.0) and verify MLE recovers parameters
            var VM = new VonMises(1.0, 3.0);
            var sample = VM.GenerateRandomValues(5000, seed: 12345);

            var fitted = new VonMises();
            fitted.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);

            Assert.AreEqual(1.0, fitted.Mu, 0.1);
            Assert.AreEqual(3.0, fitted.Kappa, 0.3);
        }

        /// <summary>
        /// Test the circular variance.
        /// </summary>
        /// <remarks>
        /// Circular variance = 1 - I1(kappa)/I0(kappa).
        /// For kappa=0: variance = 1 (uniform on circle).
        /// For kappa→∞: variance → 0 (all mass at mu).
        /// For kappa=1: A(1) = I1(1)/I0(1) ≈ 0.44606, so V ≈ 0.55394.
        /// </remarks>
        [TestMethod]
        public void Test_Variance()
        {
            // kappa = 0 → uniform on circle → variance = 1
            var VM0 = new VonMises(0, 0);
            Assert.AreEqual(1.0, VM0.Variance, 1e-6);

            // kappa = 1 → variance ≈ 0.55394
            var VM1 = new VonMises(0, 1);
            Assert.AreEqual(0.55394, VM1.Variance, 1e-3);
        }

        /// <summary>
        /// Test the Bessel function I0 against known values.
        /// </summary>
        [TestMethod]
        public void Test_BesselI0()
        {
            // I0(0) = 1
            Assert.AreEqual(1.0, Bessel.I0(0), 1e-6);
            // I0(1) ≈ 1.2660658
            Assert.AreEqual(1.2660658, Bessel.I0(1), 1e-5);
            // I0(2) ≈ 2.2795853
            Assert.AreEqual(2.2795853, Bessel.I0(2), 1e-5);
            // I0(5) ≈ 27.239872
            Assert.AreEqual(27.239872, Bessel.I0(5), 1e-3);
        }

        /// <summary>
        /// Test the Bessel function I1 against known values.
        /// </summary>
        [TestMethod]
        public void Test_BesselI1()
        {
            // I1(0) = 0
            Assert.AreEqual(0.0, Bessel.I1(0), 1e-6);
            // I1(1) ≈ 0.5651591
            Assert.AreEqual(0.5651591, Bessel.I1(1), 1e-5);
            // I1(2) ≈ 1.5906369
            Assert.AreEqual(1.5906369, Bessel.I1(2), 1e-5);
        }

        /// <summary>
        /// Test random sampling produces values in [-π, π] and recovers mean direction.
        /// </summary>
        [TestMethod]
        public void Test_RandomSampling()
        {
            var VM = new VonMises(0.5, 5.0);
            var samples = VM.GenerateRandomValues(10000, seed: 42);

            // All values should be in [-π, π]
            foreach (var s in samples)
            {
                Assert.IsGreaterThanOrEqualTo(-Math.PI, s);
                Assert.IsLessThanOrEqualTo(Math.PI, s);
            }

            // Mean direction should be close to mu
            double sumSin = 0, sumCos = 0;
            foreach (var s in samples)
            {
                sumSin += Math.Sin(s);
                sumCos += Math.Cos(s);
            }
            double meanDir = Math.Atan2(sumSin, sumCos);
            Assert.AreEqual(0.5, meanDir, 0.05);
        }

        /// <summary>
        /// Test bootstrap method.
        /// </summary>
        [TestMethod]
        public void Test_Bootstrap()
        {
            var VM = new VonMises(0, 2);
            var bootstrapped = VM.Bootstrap(ParameterEstimationMethod.MaximumLikelihood, 500, seed: 123);
            Assert.IsTrue(bootstrapped.ParametersValid);
        }
    }
}
