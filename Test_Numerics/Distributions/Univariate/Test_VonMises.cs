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
        /// The MLE rejects null, too-small, nonfinite and degenerate samples.
        /// </summary>
        [TestMethod]
        public void Test_MLE_InvalidSamples()
        {
            var VM = new VonMises();
            Assert.ThrowsExactly<ArgumentNullException>(() => VM.MLE(null!));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => VM.MLE(Array.Empty<double>()));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => VM.MLE(new[] { 0.7 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => VM.MLE(new[] { 0.1, double.NaN, 0.3 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => VM.MLE(new[] { 0.1, double.PositiveInfinity, 0.3 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => VM.MLE(new[] { double.NegativeInfinity, 0.2, 0.3 }));

            // Constant samples have no angular dispersion: R̄ = 1 and the likelihood grows without bound in κ.
            Assert.ThrowsExactly<ArgumentException>(() => VM.MLE(new[] { 0.5, 0.5, 0.5, 0.5 }));
            Assert.ThrowsExactly<ArgumentException>(() => VM.MLE(new[] { -2.9, -2.9 }));
            Assert.ThrowsExactly<ArgumentException>(() => VM.Estimate(new[] { 1.1, 1.1, 1.1 }, ParameterEstimationMethod.MaximumLikelihood));
        }

        /// <summary>
        /// A sample with zero mean resultant length gives κ = 0, up to the rounding of the trigonometric sums.
        /// </summary>
        [TestMethod]
        public void Test_MLE_ZeroResultantLength()
        {
            var VM = new VonMises();
            var parameters = VM.MLE(new[] { 0d, Math.PI / 2d, Math.PI, -Math.PI / 2d });
            Assert.AreEqual(0d, parameters[1], 1E-15);
        }

        /// <summary>
        /// The MLE matches the exact solution of A(κ) = R̄ on fixed samples.
        /// </summary>
        /// <remarks>
        /// Oracle: mpmath 1.4.1 (40 digits) atan2(Σ sin, Σ cos) and findroot(besseli(1, k)/besseli(0, k) - R̄);
        /// scipy 1.17.1 stats.vonmises.fit(sample, fscale=1) agrees (κ = 3.6654827931049607 and 1.4431177910033508).
        /// The mean direction does not depend on the Bessel ratio and is asserted bit for bit.
        /// </remarks>
        [TestMethod]
        public void Test_MLE_KnownSamples()
        {
            var VM = new VonMises();
            var sample = new double[] { 0.82, 1.35, 0.47, 1.91, 1.12, 0.66, 1.58, 0.95, 2.24, 0.31, 1.07, 1.44, 0.73, 1.76, 0.18, 1.29, 0.88, 2.05, 1.18, 0.54 };
            var parameters = VM.MLE(sample);
            Assert.AreEqual(1.117841234661942, parameters[0], 0d);
            Assert.AreEqual(1.1178412346619419, parameters[0], 1E-15);
            Assert.AreEqual(3.6654827931049603, parameters[1], 1E-14 * 3.6654827931049603);

            var sample2 = new double[] { -1.3860957223901167, -2.443019534443491, -1.4660779414560707, -1.7161501499710763, -3.1115876144062966, 2.7418430215623495, -0.8905788747480656, -2.1008128542337676, -3.0116976139367626, -2.929734295378788, 0.9356487978942987, -2.1965807819198453, 2.8849365389829087, -1.421790802788399, -2.9537664029979656 };
            var parameters2 = VM.MLE(sample2);
            Assert.AreEqual(-2.358164404828794, parameters2[0], 0d);
            Assert.AreEqual(1.4431177910033538, parameters2[1], 1E-14 * 1.4431177910033538);
        }

        /// <summary>
        /// A tightly concentrated sample gives a finite κ above the overflow range of I₀.
        /// </summary>
        [TestMethod]
        public void Test_MLE_LargeKappa()
        {
            var sample = new VonMises(0.3, 2000d).GenerateRandomValues(50, 7);
            var fitted = new VonMises();
            fitted.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            Assert.IsTrue(fitted.ParametersValid);
            Assert.AreEqual(0.3, fitted.Mu, 0.01);
            Assert.AreEqual(2000d, fitted.Kappa, 800d);

            // The fitted κ solves A(κ) = R̄ to double precision.
            double sumSin = 0d, sumCos = 0d;
            foreach (var x in sample) { sumSin += Math.Sin(x); sumCos += Math.Cos(x); }
            double rBar = Math.Sqrt(sumSin * sumSin + sumCos * sumCos) / sample.Length;
            Assert.AreEqual(rBar, Bessel.I1e(fitted.Kappa) / Bessel.I0e(fitted.Kappa), 1E-15);
        }

        /// <summary>
        /// PDF, CDF, InverseCDF and StandardDeviation over κ from 0.01 to 5000.
        /// </summary>
        /// <remarks>
        /// <para>
        /// PDF oracle: scipy 1.17.1 stats.vonmises.pdf(x, k). The relative tolerance 1E-14 max(1, κ) reflects the
        /// conditioning of exp(κ (cos x - 1)), whose relative error is κ times that of cos.
        /// </para>
        /// <para>
        /// CDF oracle: mpmath 1.4.1 (40 digits) quad(exp(k (cos t - 1)), [-π, ..., x]) / (2π besseli(0, k) exp(-k)),
        /// with breakpoints around the mode. scipy's vonmises.cdf uses a normal approximation for κ &gt; 50 that
        /// differs from mpmath by up to 4E-7 (for example 0.15906007481932388 at κ = 100, x = -0.1), so it is not
        /// used here.
        /// </para>
        /// <para>
        /// InverseCDF oracle: mpmath findroot of the CDF above. The 5E-9 absolute tolerance is that of the Brent
        /// root finder.
        /// </para>
        /// <para>
        /// StandardDeviation oracle: sqrt(1 - scipy.special.i1e(k) / scipy.special.i0e(k)) at 1E-15 relative, and
        /// mpmath sqrt(1 - besseli(1, k) / besseli(0, k)) at 1E-11 relative (the subtraction 1 - A(κ) loses about
        /// log10(κ) digits for large κ). scipy's vonmises.std is the linear standard deviation and does not apply.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_Oracle_KappaGrid()
        {
            double[] kappas = { 0.01, 1.0, 10.0, 100.0, 500.0, 800.0, 5000.0 };
            double[] xs = { -1.0, -0.1, 0.1, 0.5, 2.0 };
            double[] ps = { 0.001, 0.1, 0.9, 0.999 };
            // scipy.stats.vonmises.pdf(x, k)
            double[,] pdf = { { 0.16001318783208981, 0.16074242744907072, 0.16074242744907072, 0.16055379188456195, 0.15849003875385817 },
                              { 0.21578146511029628, 0.3400076166746798, 0.3400076166746798, 0.3023382476531484, 0.08291508547317152 },
                              { 0.012552617230965511, 1.1843481078434028, 1.1843481078434028, 0.36603602164656857, 8.809037084481075E-07 },
                              { 4.324613145875511E-20, 2.4176760562675517, 2.4176760562675517, 1.9223951790410252E-05, 1.252820314712858E-61 },
                              { 1.3433848090714657E-99, 0.733592123039229, 0.733592123039229, 2.3317256797193976E-26, 2.740993760844653E-307 },
                              { 2.1729495047091729E-159, 0.20732725444587974, 0.20732725444587974, 3.312934660253534E-42, 0.0 },
                              { 0.0, 4.0000640763224805E-10, 4.0000640763224805E-10, 4.210124228805052E-265, 0.0 } };
            // mpmath CDF; the value 4.38E-1001 at κ = 5000, x = -1 underflows to zero
            double[,] cdf = { { 0.33950402171138178, 0.48392522197192228, 0.51607477802807772, 0.58034216768359747, 0.81975555369517935 },
                              { 0.20564469256531652, 0.46588578919208413, 0.53411421080791587, 0.66407647460499303, 0.93424095588998317 },
                              { 0.0014108017064043275, 0.37754136491149819, 0.62245863508850181, 0.93864493191090481, 0.99999989481852285 },
                              { 5.1016782326113297E-22, 0.15906045738451988, 0.84093954261548012, 0.99999961271360566, 1.0 },
                              { 3.1932020762915065E-102, 0.012722585013892834, 0.98727741498610717, 1.0, 1.0 },
                              { 3.240777964596878E-162, 0.0023507385260872993, 0.9976492614739127, 1.0, 1.0 },
                              { 0.0, 7.8622442488496296E-13, 0.99999999999921378, 1.0, 1.0 } };
            // mpmath findroot(CDF(q) - p)
            double[,] ppf = { { -3.1352461629871933, -2.5073605251989272, 2.5073605251989274, 3.1352461629871933 },
                              { -3.1199706062242234, -1.6515617183319102, 1.6515617183319104, 3.1199706062242234 },
                              { -1.0383338432211791, -0.41371903020106368, 0.41371903020106373, 1.038333843221179 },
                              { -0.31066424322376332, -0.1284047905917765, 0.12840479059177651, 0.3106642432237633 },
                              { -0.13834436025024223, -0.057334931460261367, 0.057334931460261374, 0.13834436025024222 },
                              { -0.10932775946460887, -0.045320655086987514, 0.04532065508698752, 0.10932775946460886 },
                              { -0.043707056157497312, -0.018124577294593212, 0.018124577294593214, 0.043707056157497308 } };
            // sqrt(1 - scipy.special.i1e(k) / scipy.special.i0e(k))
            double[] sd = { 0.9974968984908967, 0.7440497524382798, 0.22671606481490028, 0.07079990815552865, 0.03163069716490751, 0.02500391083648131, 0.010000250046898407 };
            // mpmath sqrt(1 - besseli(1, k) / besseli(0, k))
            double[] sdExact = { 0.99749689849089674, 0.74404975243827982, 0.22671606481490023, 0.070799908155528242, 0.031630697164907186, 0.025003910836481279, 0.010000250046889459 };

            for (int i = 0; i < kappas.Length; i++)
            {
                var VM = new VonMises(0d, kappas[i]);
                for (int j = 0; j < xs.Length; j++)
                {
                    Assert.AreEqual(pdf[i, j], VM.PDF(xs[j]), 1E-14 * Math.Max(1d, kappas[i]) * pdf[i, j], $"PDF κ={kappas[i]}, x={xs[j]}");
                    Assert.AreEqual(cdf[i, j], VM.CDF(xs[j]), 1E-12, $"CDF κ={kappas[i]}, x={xs[j]}");
                }
                for (int j = 0; j < ps.Length; j++)
                    Assert.AreEqual(ppf[i, j], VM.InverseCDF(ps[j]), 5E-9, $"InverseCDF κ={kappas[i]}, p={ps[j]}");
                Assert.AreEqual(sd[i], VM.StandardDeviation, 1E-15 * sd[i], $"SD κ={kappas[i]}");
                Assert.AreEqual(sdExact[i], VM.StandardDeviation, 1E-11 * sdExact[i], $"SD exact κ={kappas[i]}");
            }
        }

        /// <summary>
        /// The CDF is accurate with a nonzero mean direction, and its tail values keep relative accuracy.
        /// </summary>
        /// <remarks>
        /// Oracle: mpmath 1.4.1 quad of exp(k (cos(t - 1) - 1)) from -π to x, split at μ = 1, divided by
        /// 2π besseli(0, k) exp(-k); the tail values use μ = 0.
        /// </remarks>
        [TestMethod]
        public void Test_CDF_NonzeroMeanAndTails()
        {
            var VM = new VonMises(1.0, 2.0);
            Assert.AreEqual(0.0099712386632325323, VM.CDF(-2.5), 1E-14);
            Assert.AreEqual(0.4620903193977252, VM.CDF(0.9), 1E-14);
            Assert.AreEqual(0.5392802929733679, VM.CDF(1.05), 1E-14);
            Assert.AreEqual(0.99619770517346191, VM.CDF(3.0), 1E-14);

            var VM2 = new VonMises(1.0, 800.0);
            Assert.AreEqual(0.0023507385260873051, VM2.CDF(0.9), 1E-14);
            Assert.AreEqual(0.92129630807929995, VM2.CDF(1.05), 1E-13);
            Assert.AreEqual(1.0, VM2.CDF(3.0), 1E-13);

            Assert.AreEqual(5.1016782326113297E-22, new VonMises(0d, 100d).CDF(-1.0), 1E-12 * 5.1016782326113297E-22);
            Assert.AreEqual(3.7585362033492754E-10, new VonMises(0d, 10d).CDF(-3.0), 1E-12 * 3.7585362033492754E-10);
        }

        /// <summary>
        /// For κ beyond the overflow range of exp(κ) and I₀(κ), every property is finite and accurate.
        /// </summary>
        /// <remarks>
        /// Unscaled exp(κ cos(x - μ)) / I₀(κ) overflows at κ = 800; the scaled evaluation keeps every property
        /// finite. Oracle: scipy.stats.vonmises.pdf(0.1, 800) =
        /// 0.20732725444587974; mpmath 1.4.1 for μ = -0.5, κ = 1E5: PDF(μ) = 1 / (2π besseli(0, k) exp(-k)) =
        /// 126.15646840453545, CDF(-0.499) = 0.6240850279713659 and CDF(-0.501) = 0.3759149720286341.
        /// </remarks>
        [TestMethod]
        public void Test_LargeKappa_Finite()
        {
            var VM = new VonMises(0d, 800d);
            Assert.AreEqual(0.20732725444587974, VM.PDF(0.1), 1E-12 * 0.20732725444587974);
            Assert.AreEqual(0.025003910836481279, VM.StandardDeviation, 1E-13);
            Assert.AreEqual(0.5, VM.CDF(0d), 1E-13);
            Assert.AreEqual(0d, VM.InverseCDF(0.5), 1E-9);

            var VM2 = new VonMises(-0.5, 1E5);
            Assert.AreEqual(126.15646840453545, VM2.PDF(-0.5), 1E-13 * 126.15646840453545);
            Assert.AreEqual(0.5, VM2.CDF(-0.5), 1E-11);
            Assert.AreEqual(0.6240850279713659, VM2.CDF(-0.499), 1E-11);
            Assert.AreEqual(0.3759149720286341, VM2.CDF(-0.501), 1E-11);
            Assert.AreEqual(1d, VM2.CDF(0.5), 1E-11);
            Assert.AreEqual(0d, VM2.CDF(-1.5), 1E-11);
            Assert.AreEqual(-0.5, VM2.InverseCDF(0.5), 5E-9);
        }

        /// <summary>
        /// The density keeps full relative accuracy at large κ.
        /// </summary>
        [TestMethod()]
        public void Test_PDF_LargeKappaAccuracy()
        {
            // mpmath (50 digits): exp(k*(cos(x) - 1)) / (2*pi*besseli(0, k)*exp(-k))
            Assert.AreEqual(21.969060954882484, new VonMises(0d, 5000d).PDF(0.01), 21.969060954882484 * 1E-14);
            Assert.AreEqual(2.9736846445490933, new VonMises(0d, 5000d).PDF(0.03), 2.9736846445490933 * 1E-14);
            Assert.AreEqual(103.28818726838194, new VonMises(0d, 1E5).PDF(0.002), 103.28818726838194 * 1E-14);
        }

        /// <summary>
        /// For κ = 0 the distribution is uniform on [-π, π].
        /// </summary>
        [TestMethod]
        public void Test_ZeroKappa_Uniform()
        {
            var VM = new VonMises(0.4, 0d);
            Assert.AreEqual(1d / (2d * Math.PI), VM.PDF(-2.0), 1E-15);
            Assert.AreEqual(0.75, VM.CDF(Math.PI / 2d), 1E-14);
            Assert.AreEqual(1d, VM.StandardDeviation, 0d);
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
