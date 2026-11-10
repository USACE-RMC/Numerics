using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Kappa-4 distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_KappaFour
    {

        /// <summary>
        /// Verifies Kappa Four fitting with linear moments.
        /// </summary>
        [TestMethod]
        public void Test_K4_LMOM()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double xi = kappa4.Xi;
            double a = kappa4.Alpha;
            double k = kappa4.Kappa;
            double h= kappa4.Hondo;
            double true_xi = 8.68360234;
            double true_a = 3.10384972;
            double true_k = 0.14470737;
            double true_h = -0.07348014;

            Assert.AreEqual(xi, true_xi, 0.0001d);
            Assert.AreEqual(a, true_a, 0.0001d);
            Assert.AreEqual(k, true_k, 0.0001d);
            Assert.AreEqual(h, true_h, 0.0001d);

            var lmom = kappa4.LinearMomentsFromParameters(kappa4.GetParameters);
            Assert.AreEqual(9.95751634d, lmom[0],  0.0001d);
            Assert.AreEqual(1.98224114d, lmom[1],  0.0001d);
            Assert.AreEqual(0.06380885d, lmom[2],  0.0001d);
            Assert.AreEqual(0.12442297d, lmom[3],  0.0001d);
        }

        /// <summary>
        /// Verification using R lmom package.
        /// </summary>
        [TestMethod]
        public void Test_K4_CDF()
        {
            var x = new double[] { 5, 10, 12, 15, 18 };
            var p = new double[5];
            var true_p = new double[] { 0.07168831, 0.53317660, 0.73279234, 0.91293987, 0.97980084 };
            var k4 = new KappaFour(8.7, 3.1, 0.14, -0.1);
            for (int i = 0; i < 5; i++)
            {
                p[i] = k4.CDF(x[i]);
                Assert.AreEqual(true_p[i], p[i], 1E-7);
            }
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(x[i], k4.InverseCDF(p[i]), 1E-7);
            }

        }


        /// <summary>
        /// Verifies Kappa Four distribution values.
        /// </summary>
        [TestMethod]
        public void Test_K4_Dist()
        {
            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);

            double pdf = kappa4.PDF(14.9);
            double cdf = kappa4.CDF(14.9);
            double invcdf = kappa4.InverseCDF(cdf);
            double true_pdf = 0.0385446;
            double true_cdf = 0.9106253;
            double true_invcdf = 14.9;

            Assert.AreEqual(pdf, true_pdf, 0.0001d);
            Assert.AreEqual(cdf, true_cdf, 0.0001d);
            Assert.AreEqual(invcdf, true_invcdf, 0.0001d);

        }

        /// <summary>
        /// Pins the implemented zero-kappa closed form: the density must equal
        /// exp(-(x - xi) / alpha) / alpha * F(x)^(1 - h), evaluated with the distribution's own
        /// CDF. The independent evidence for the zero-kappa branch lives in the companion
        /// normalization and continuity tests.
        /// </summary>
        [TestMethod]
        public void Test_K4_ZeroKappa_PDFMatchesClosedForm()
        {
            const double xi = 2.5d;
            const double alpha = 1.75d;
            const double x = 4.0d;
            foreach (double hondo in new[] { -0.2d, 0d, 0.2d })
            {
                var distribution = new KappaFour(xi, alpha, 0d, hondo);
                double standardizedValue = (x - xi) / alpha;
                double expected = Math.Exp(-standardizedValue) / alpha
                    * Math.Pow(distribution.CDF(x), 1d - hondo);

                Assert.AreEqual(expected, distribution.PDF(x), 1E-12d);
            }
        }

        /// <summary>
        /// Verifies the zero-kappa inverse CDF algebra and CDF round trip.
        /// </summary>
        [TestMethod]
        public void Test_K4_ZeroKappa_InverseCDFMatchesAnalyticalSolution()
        {
            const double xi = 2.5d;
            const double alpha = 1.75d;
            foreach (double hondo in new[] { -0.2d, 0.2d })
            {
                var distribution = new KappaFour(xi, alpha, 0d, hondo);
                foreach (double probability in new[] { 0.1d, 0.5d, 0.9d })
                {
                    double expected = xi - alpha
                        * Math.Log((1d - Math.Pow(probability, hondo)) / hondo);
                    double quantile = distribution.InverseCDF(probability);

                    Assert.AreEqual(expected, quantile, 1E-12d);
                    Assert.AreEqual(probability, distribution.CDF(quantile), 1E-12d);
                }
            }
        }

        /// <summary>
        /// Verifies normalization and support for the positive-hondo zero-kappa case.
        /// </summary>
        [TestMethod]
        public void Test_K4_ZeroKappa_PositiveHondoHasNormalizedSupportedDensity()
        {
            var distribution = new KappaFour(2.5d, 1.75d, 0d, 0.2d);
            const double lowerProbability = 1E-9d;
            const double upperProbability = 1d - lowerProbability;
            double lower = distribution.InverseCDF(lowerProbability);
            double upper = distribution.InverseCDF(upperProbability);
            var integrator = new AdaptiveGaussKronrod(distribution.PDF, lower, upper);

            integrator.Integrate();

            Assert.AreEqual(2.5d + 1.75d * Math.Log(0.2d), distribution.Minimum, 1E-12d);
            Assert.AreEqual(0d, distribution.PDF(distribution.Minimum - 1E-6d));
            Assert.AreEqual(upperProbability - lowerProbability, integrator.Result, 1E-8d);
        }

        /// <summary>
        /// Verifies continuity of the general Kappa Four formulas as kappa approaches zero.
        /// </summary>
        [TestMethod]
        public void Test_K4_GeneralFormulaIsContinuousAtZeroKappa()
        {
            const double xi = 2.5d;
            const double alpha = 1.75d;
            const double hondo = 0.2d;
            const double x = 4.0d;
            const double probability = 0.7d;
            var limit = new KappaFour(xi, alpha, 0d, hondo);

            foreach (double kappa in new[] { -1E-7d, 1E-7d })
            {
                var nearby = new KappaFour(xi, alpha, kappa, hondo);
                Assert.AreEqual(limit.CDF(x), nearby.CDF(x), 1E-7d);
                Assert.AreEqual(limit.PDF(x), nearby.PDF(x), 1E-7d);
                Assert.AreEqual(limit.InverseCDF(probability), nearby.InverseCDF(probability), 1E-6d);
            }
        }

        /// <summary>
        /// Verifies that representative finite shape pairs are valid and define internally consistent supports.
        /// </summary>
        [TestMethod]
        public void Test_K4_FiniteShapePairsHaveConsistentSupport()
        {
            (double Kappa, double Hondo)[] shapePairs =
            [
                (-0.25d, -0.5d),
                (-0.25d, 0.5d),
                (0.25d, -0.5d),
                (0.25d, 0.5d)
            ];
            double[] probabilities = [1E-6d, 0.1d, 0.5d, 0.9d, 1d - 1E-6d];

            foreach ((double kappa, double hondo) in shapePairs)
            {
                var distribution = new KappaFour(2.5d, 1.75d, kappa, hondo);
                Assert.IsTrue(distribution.ParametersValid, $"Finite shape pair (kappa={kappa}, h={hondo}) must be admissible.");

                double previousQuantile = double.NegativeInfinity;
                foreach (double probability in probabilities)
                {
                    double quantile = distribution.InverseCDF(probability);
                    Assert.IsTrue(Tools.IsFinite(quantile));
                    Assert.IsGreaterThan(previousQuantile, quantile);
                    Assert.AreEqual(probability, distribution.CDF(quantile), 1E-10d);
                    Assert.IsTrue(quantile >= distribution.Minimum && quantile <= distribution.Maximum);
                    previousQuantile = quantile;
                }

                double medianDensity = distribution.PDF(distribution.InverseCDF(0.5d));
                Assert.IsTrue(Tools.IsFinite(medianDensity) && medianDensity > 0d);
                if (Tools.IsFinite(distribution.Minimum))
                {
                    Assert.AreEqual(0d, distribution.CDF(distribution.Minimum));
                    Assert.AreEqual(0d, distribution.PDF(distribution.Minimum - 1d));
                }
                if (Tools.IsFinite(distribution.Maximum))
                {
                    Assert.AreEqual(1d, distribution.CDF(distribution.Maximum));
                    Assert.AreEqual(0d, distribution.PDF(distribution.Maximum + 1d));
                }
            }
        }

        /// <summary>
        /// Verifies Kappa Four partial derivative calculations.
        /// </summary>
        [TestMethod]
        public void Test_K4_PartialDerivatives()
        {

            // Air quality - wind data from R
            var data = new double[] { 7.4, 8, 12.6, 11.5, 14.3, 14.9, 8.6, 13.8, 20.1, 8.6, 6.9, 9.7, 9.2, 10.9, 13.2, 11.5, 12, 18.4, 11.5, 9.7, 9.7, 16.6, 9.7, 12, 16.6, 14.9, 8, 12, 14.9, 5.7, 7.4, 8.6, 9.7, 16.1, 9.2, 8.6, 14.3, 9.7, 6.9, 13.8, 11.5, 10.9, 9.2, 8, 13.8, 11.5, 14.9, 20.7, 9.2, 11.5, 10.3, 6.3, 1.7, 4.6, 6.3, 8, 8, 10.3, 11.5, 14.9, 8, 4.1, 9.2, 9.2, 10.9, 4.6, 10.9, 5.1, 6.3, 5.7, 7.4, 8.6, 14.3, 14.9, 14.9, 14.3, 6.9, 10.3, 6.3, 5.1, 11.5, 6.9, 9.7, 11.5, 8.6, 8, 8.6, 12, 7.4, 7.4, 7.4, 9.2, 6.9, 13.8, 7.4, 6.9, 7.4, 4.6, 4, 10.3, 8, 8.6, 11.5, 11.5, 11.5, 9.7, 11.5, 10.3, 6.3, 7.4, 10.9, 10.3, 15.5, 14.3, 12.6, 9.7, 3.4, 8, 5.7, 9.7, 2.3, 6.3, 6.3, 6.9, 5.1, 2.8, 4.6, 7.4, 15.5, 10.9, 10.3, 10.9, 9.7, 14.9, 15.5, 6.3, 10.9, 11.5, 6.9, 13.8, 10.3, 10.3, 8, 12.6, 9.2, 10.3, 10.3, 16.6, 6.9, 13.2, 14.3, 8, 11.5 };
            var kappa4 = new KappaFour();
            kappa4.Estimate(data, ParameterEstimationMethod.MethodOfLinearMoments);
            double p = 0.999;

            var pd1 = kappa4.QuantileGradient(p);
            var pd2 = NumericalDerivative.Gradient(x =>
            {
                var K4 = new KappaFour();
                K4.SetParameters(x);
                return K4.InverseCDF(p);
            },kappa4.GetParameters);

            for (int i = 0; i < pd1.Length; i++)
            {
                Assert.AreEqual(pd1[i], pd2[i], 1E-2);
            }
        }

        /// <summary>
        /// Checking that Kappa-4 is being created with inputs.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var k4 = new KappaFour();
            Assert.AreEqual(100,k4.Xi);
            Assert.AreEqual(10, k4.Alpha);
            Assert.AreEqual(0, k4.Kappa);
            Assert.AreEqual(0, k4.Hondo);

            var k4ii = new KappaFour(100, 10, 1, 1);
            Assert.AreEqual(100, k4ii.Xi);
            Assert.AreEqual(10, k4ii.Alpha);
            Assert.AreEqual(1, k4ii.Kappa);
            Assert.AreEqual(1, k4ii.Hondo);
        }

        /// <summary>
        /// Testting Kappa-4 with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var k4 = new KappaFour(double.NaN,double.NaN,double.NaN, double.NaN);
            Assert.IsFalse(k4.ParametersValid);

            var k4ii = new KappaFour(double.PositiveInfinity,double.PositiveInfinity,double.PositiveInfinity,double.PositiveInfinity);
            Assert.IsFalse(k4ii.ParametersValid);

            var k4iii = new KappaFour(100, 0, 1, 1);
            Assert.IsFalse(k4iii.ParametersValid);
        }

        /// <summary>
        /// Testing ParametersToString
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var k4 = new KappaFour();
            Assert.AreEqual("Location (ξ)",k4.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (α)", k4.ParametersToString[1, 0]);
            Assert.AreEqual("Shape (κ)", k4.ParametersToString[2, 0]);
            Assert.AreEqual("Shape (h)", k4.ParametersToString[3, 0]);
            Assert.AreEqual("100", k4.ParametersToString[0, 1]);
            Assert.AreEqual("10", k4.ParametersToString[1, 1]);
            Assert.AreEqual("0", k4.ParametersToString[2, 1]);
            Assert.AreEqual("0", k4.ParametersToString[3, 1]);
        }

        /// <summary>
        /// A NaN probability is rejected by the inverse CDF.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_NaNProbability()
        {
            var k4 = new KappaFour(100, 10, 0.1, 0.2);
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => k4.InverseCDF(double.NaN));
        }

        /// <summary>
        /// The density matches the support checks written against <see cref="KappaFour.Minimum"/> and
        /// <see cref="KappaFour.Maximum"/> on a grid that includes points outside the support, the upper endpoint
        /// and interior points.
        /// </summary>
        [TestMethod]
        public void Test_PDF_MatchesSupportExpression()
        {
            var cases = new[] { new KappaFour(100, 10, 0.1, 0.2), new KappaFour(100, 10, -0.2, -0.5), new KappaFour(100, 10, 1.5, 0.3), new KappaFour(0, 1, 2, 1E-300) };
            foreach (var k4 in cases)
            {
                double lower = k4.Minimum, upper = k4.Maximum;
                var grid = new System.Collections.Generic.List<double> { double.NegativeInfinity, double.PositiveInfinity, upper };
                if (!double.IsInfinity(lower)) { grid.Add(lower - 1); grid.Add(lower + 1E-6); }
                if (!double.IsInfinity(upper)) { grid.Add(upper + 1); grid.Add(upper - 1E-6); }
                double start = double.IsInfinity(lower) ? k4.InverseCDF(0.01) : lower, end = double.IsInfinity(upper) ? k4.InverseCDF(0.99) : upper;
                for (int i = 1; i < 10; i++) grid.Add(start + (end - start) * i / 10d);
                foreach (double x in grid)
                {
                    double expected;
                    if (double.IsInfinity(x) || x < k4.Minimum || x > k4.Maximum) expected = 0d;
                    else if (x == k4.Maximum) expected = k4.Kappa < 1d ? 0d : k4.Kappa == 1d ? 1d / k4.Alpha : double.PositiveInfinity;
                    else if (x == k4.Minimum) continue;
                    else expected = System.Math.Exp(k4.LogPDF(x));
                    Assert.AreEqual(expected, k4.PDF(x), 0.0, $"x {x:R}");
                    double expectedLog = double.IsInfinity(x) || x < k4.Minimum || x > k4.Maximum ? double.NegativeInfinity
                        : x == k4.Maximum ? (expected > 0d ? System.Math.Log(expected) : double.NegativeInfinity) : k4.LogPDF(x);
                    Assert.AreEqual(expectedLog, k4.LogPDF(x), 0.0, $"x {x:R}");
                }
            }
        }

        /// <summary>
        /// When the compensated lower endpoint cannot be evaluated (|κ log h| above 700) the analytical endpoint is
        /// returned on every access, and changing a parameter recomputes the endpoint.
        /// </summary>
        [TestMethod]
        public void Test_Minimum_CompensatedEndpointFailure()
        {
            Assert.IsFalse(KappaFourBoundary.TryLowerEndpoint(0, 1, 2, 1E-300, out _));
            var k4 = new KappaFour(0, 1, 2, 1E-300);
            double first = k4.Minimum;
            Assert.AreEqual(first, k4.Minimum, 0.0);
            Assert.AreEqual(new KappaFour(0, 1, 2, 1E-300).Minimum, first, 0.0);
            k4.Kappa = 0.1;
            Assert.IsTrue(KappaFourBoundary.TryLowerEndpoint(0, 1, 0.1, 1E-300, out double compensated));
            Assert.AreEqual(compensated, k4.Minimum, 0.0);
            k4.Kappa = 2;
            Assert.AreEqual(first, k4.Minimum, 0.0);
            k4.Hondo = 0.5;
            Assert.AreEqual(new KappaFour(0, 1, 2, 0.5).Minimum, k4.Minimum, 0.0);
        }

        /// <summary>
        /// The closed-form L-moments agree with the integrated branch on both sides of the |κ| = 1E-3 and |h| = 1E-3
        /// switch, so the switch has no visible jump.
        /// </summary>
        /// <remarks>
        /// Reference values for ξ = 10 and α = 2 computed with Python mpmath 1.3.0 at 60 significant digits from the
        /// probability-weighted moments g_r = r Γ(1 + κ)Γ(r/h)/(h^(1 + κ) Γ(1 + κ + r/h)) and the corresponding h &lt; 0
        /// form. The tolerances are 1E-14 relative for λ1 and λ2 and 5E-14 absolute for τ3 and τ4.
        /// </remarks>
        [TestMethod]
        [DataRow(1.0001E-3, 0.5, 11.611806369381824, 1.1652130285604996, 0.2566167063892183, 0.15274541345776865)]
        [DataRow(0.9999E-3, 0.5, 11.611806748811835, 1.1652133189146059, 0.25661681157679006, 0.1527454765504469)]
        [DataRow(-1.0001E-3, -0.3, 10.841733433984263, 1.5484856166284047, 0.11614608171186738, 0.1539005969063767)]
        [DataRow(-0.9999E-3, -0.3, 10.841732999219193, 1.5484854037514884, 0.11614593481821361, 0.1539005580784055)]
        [DataRow(0.1, 1.0001E-3, 10.974030955066072, 1.2736241449291492, 0.10745707129714385, 0.1289016707921116)]
        [DataRow(0.1, 0.9999E-3, 10.974030745864114, 1.273624256512736, 0.10745702942238328, 0.12890167578315836)]
        [DataRow(-0.2, -1.0001E-3, 11.64136560448079, 1.7315869564248136, 0.3049652393972515, 0.21800793124654994)]
        [DataRow(-0.2, -0.9999E-3, 11.641365790779878, 1.731586877120042, 0.30496526492759024, 0.21800793510072788)]
        [DataRow(1.0001E-3, 1.0001E-3, 11.153455189286937, 1.3845148644476606, 0.16946425165295032, 0.1501168358775242)]
        [DataRow(-1.0001E-3, -1.0001E-3, 11.155411613302801, 1.3880765113072797, 0.17038661393273602, 0.15063364195259785)]
        [DataRow(0.3, 0.01, 10.695153466140328, 1.1164247518990507, -0.006343883233203591, 0.10549200453289603)]
        [DataRow(-0.5, 20.0, 25.09910570460915, 1.1562480487707685, 0.9161754131378103, 0.8118370212481122)]
        public void Test_LinearMomentsFromParameters_ContinuousAcrossIntegrationSwitch(double kappa, double hondo, double l1, double l2, double t3, double t4)
        {
            var k4 = new KappaFour(10d, 2d, kappa, hondo);
            double[] lmom = k4.LinearMomentsFromParameters(new[] { 10d, 2d, kappa, hondo });
            Assert.AreEqual(l1, lmom[0], 1E-14 * l1);
            Assert.AreEqual(l2, lmom[1], 1E-14 * l2);
            Assert.AreEqual(t3, lmom[2], 5E-14);
            Assert.AreEqual(t4, lmom[3], 5E-14);
        }
     }
}
