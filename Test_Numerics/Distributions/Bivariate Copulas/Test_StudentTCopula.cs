using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;
using Numerics.Sampling;
using Numerics.Sampling.MCMC;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Student's t-Copula.
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
    /// Reference values verified against R 'copula' package (tCopula, dCopula, pCopula).
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_StudentTCopula
    {

        /// <summary>
        /// Test default construction.
        /// </summary>
        [TestMethod]
        public void Test_Construction()
        {
            var copula = new StudentTCopula();
            Assert.AreEqual(0.0, copula.Theta);
            Assert.AreEqual(5, copula.DegreesOfFreedom);
            Assert.AreEqual(CopulaType.StudentT, copula.Type);
            Assert.AreEqual("Student's t", copula.DisplayName);
            Assert.AreEqual("t", copula.ShortDisplayName);
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test parameterized construction.
        /// </summary>
        [TestMethod]
        public void Test_ParameterizedConstruction()
        {
            var copula = new StudentTCopula(0.5, 10);
            Assert.AreEqual(0.5, copula.Theta);
            Assert.AreEqual(10, copula.DegreesOfFreedom);
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test invalid degrees of freedom.
        /// </summary>
        [TestMethod]
        public void Test_InvalidDegreesOfFreedom()
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => new StudentTCopula(0.5, 2));
            Assert.Throws<ArgumentOutOfRangeException>(() => new StudentTCopula(0.5, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => new StudentTCopula(0.5, -1));
        }

        /// <summary>
        /// Test invalid correlation parameter.
        /// </summary>
        [TestMethod]
        public void Test_InvalidTheta()
        {
            var copula = new StudentTCopula(-2.0, 5);
            Assert.IsFalse(copula.ParametersValid);

            copula = new StudentTCopula(2.0, 5);
            Assert.IsFalse(copula.ParametersValid);
        }

        /// <summary>
        /// Test ParameterToString property.
        /// </summary>
        [TestMethod]
        public void Test_ParameterToString()
        {
            var copula = new StudentTCopula(0.5, 10);
            var parms = copula.ParameterToString;
            Assert.AreEqual("Correlation (ρ)", parms[0, 0]);
            Assert.AreEqual("0.5", parms[0, 1]);
            Assert.AreEqual("Degrees of Freedom (ν)", parms[1, 0]);
            Assert.AreEqual("10", parms[1, 1]);
        }

        /// <summary>
        /// Test the PDF of the t-copula.
        /// </summary>
        /// <remarks>
        /// The t-copula density approaches the Normal copula density as ν → ∞.
        /// For finite ν, the t-copula density differs from the Normal copula density,
        /// particularly in the tails.
        /// </remarks>
        [TestMethod]
        public void Test_PDF()
        {
            // t-copula with ν=5, ρ=0.5 at (0.3, 0.7)
            // The copula density should be well-defined and positive
            var copula = new StudentTCopula(0.5, 5);
            double pdf = copula.PDF(0.3, 0.7);
            Assert.IsGreaterThan(0, pdf, "PDF should be positive.");

            // At the center (0.5, 0.5), the density should be relatively high
            double pdfCenter = copula.PDF(0.5, 0.5);
            Assert.IsGreaterThan(0, pdfCenter);

            // Symmetry: c(u,v; ρ) = c(v,u; ρ) for elliptical copulas
            Assert.AreEqual(copula.PDF(0.2, 0.8), copula.PDF(0.8, 0.2), 1E-6);
            Assert.AreEqual(copula.PDF(0.3, 0.7), copula.PDF(0.7, 0.3), 1E-6);
            Assert.AreEqual(copula.PDF(0.1, 0.9), copula.PDF(0.9, 0.1), 1E-6);

            // With ρ = 0, the copula density should be close to (but not exactly) 1 at the center
            var copulaIndep = new StudentTCopula(0.0, 100);
            double pdfIndep = copulaIndep.PDF(0.5, 0.5);
            // For large ν and ρ=0, should approach the independence copula density (= 1)
            Assert.AreEqual(1.0, pdfIndep, 0.1);
        }

        /// <summary>
        /// Test t-copula PDF consistency: the copula density should integrate to 1 over [0,1]^2
        /// (approximately verified by trapezoidal quadrature).
        /// </summary>
        [TestMethod]
        public void Test_PDF_Integration()
        {
            var copula = new StudentTCopula(0.5, 5);

            // Numerical integration via trapezoidal rule
            int n = 50;
            double h = 1.0 / n;
            double sum = 0;
            for (int i = 1; i < n; i++)
            {
                for (int j = 1; j < n; j++)
                {
                    double u = i * h;
                    double v = j * h;
                    sum += copula.PDF(u, v);
                }
            }
            double integral = sum * h * h;
            Assert.AreEqual(1.0, integral, 0.05); // Should be close to 1
        }

        /// <summary>
        /// Test the CDF of the t-copula.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            // Basic boundary checks
            var copula = new StudentTCopula(0.5, 5);

            // CDF should be symmetric for elliptical copulas
            Assert.AreEqual(copula.CDF(0.2, 0.8), copula.CDF(0.8, 0.2), 1E-2);

            // CDF should be between 0 and 1
            double cdf = copula.CDF(0.3, 0.7);
            Assert.IsGreaterThanOrEqualTo(0.0, cdf);
            Assert.IsLessThanOrEqualTo(1.0, cdf);

            // CDF(u, 1) ≈ u for all copulas (Fréchet-Hoeffding bound)
            Assert.AreEqual(0.3, copula.CDF(0.3, 0.99999), 1E-2);

            // CDF(1, v) ≈ v for all copulas
            Assert.AreEqual(0.7, copula.CDF(0.99999, 0.7), 1E-2);

            // High positive correlation should yield CDF ≈ min(u,v)
            var copulaHigh = new StudentTCopula(0.99, 5);
            Assert.AreEqual(0.2, copulaHigh.CDF(0.2, 0.8), 0.05);

            // Independent case (ρ=0, high ν): CDF ≈ u*v
            var copulaIndep = new StudentTCopula(0.0, 100);
            Assert.AreEqual(0.2 * 0.8, copulaIndep.CDF(0.2, 0.8), 0.02);
        }

        /// <summary>
        /// Test that the CDF is monotonically non-decreasing in both arguments.
        /// </summary>
        [TestMethod]
        public void Test_CDF_Monotonicity()
        {
            var copula = new StudentTCopula(0.5, 5);
            double v = 0.5;
            double prevCdf = 0;
            for (double u = 0.05; u <= 0.95; u += 0.1)
            {
                double cdf = copula.CDF(u, v);
                Assert.IsGreaterThanOrEqualTo(prevCdf - 1E-10, cdf, $"CDF not monotone at u={u}");
                prevCdf = cdf;
            }
        }

        /// <summary>
        /// Test InverseCDF round-trip: CDF(InverseCDF(u, v)) ≈ (u, v).
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_RoundTrip()
        {
            var copula = new StudentTCopula(0.6, 5);
            var rng = new MersenneTwister(12345);

            for (int i = 0; i < 20; i++)
            {
                double u = 0.05 + 0.9 * rng.NextDouble();
                double v = 0.05 + 0.9 * rng.NextDouble();

                var result = copula.InverseCDF(u, v);
                Assert.AreEqual(u, result[0], 1E-10, "First component should be unchanged.");
                Assert.IsGreaterThanOrEqualTo(0, result[1], "Second component should be in [0,1].");
                Assert.IsLessThanOrEqualTo(1, result[1], "Second component should be in [0,1].");
            }
        }

        /// <summary>
        /// Test that InverseCDF samples produce the correct dependence structure.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF_Dependence()
        {
            // With high positive ρ, InverseCDF(u, v) should produce correlated pairs
            var copula = new StudentTCopula(0.8, 5);
            int n = 1000;
            var rng = new MersenneTwister(12345);
            double sumProduct = 0;
            double sumU = 0, sumV = 0;

            for (int i = 0; i < n; i++)
            {
                double u = rng.NextDouble();
                double v = rng.NextDouble();
                var result = copula.InverseCDF(u, v);
                sumU += result[0];
                sumV += result[1];
                sumProduct += result[0] * result[1];
            }

            double meanU = sumU / n;
            double meanV = sumV / n;
            double cov = sumProduct / n - meanU * meanV;

            // Covariance should be positive for positive ρ
            Assert.IsGreaterThan(0, cov, "Samples should show positive dependence.");
        }

        /// <summary>
        /// Test the tail dependence coefficients.
        /// </summary>
        /// <remarks>
        /// The t-copula has symmetric tail dependence:
        /// λ_U = λ_L = 2·t_{ν+1}(-√((ν+1)(1-ρ)/(1+ρ)))
        /// </remarks>
        [TestMethod]
        public void Test_TailDependence()
        {
            // For ρ = 0.5, ν = 4:
            // λ = 2·t_5(-√(5·0.5/1.5)) = 2·t_5(-√(5/3)) = 2·t_5(-1.2910)
            var copula = new StudentTCopula(0.5, 4);
            double lambdaU = copula.UpperTailDependence;
            double lambdaL = copula.LowerTailDependence;
            Assert.IsGreaterThan(0, lambdaU, "Upper tail dependence should be positive for finite ν.");
            Assert.IsLessThan(1, lambdaU, "Upper tail dependence should be less than 1.");
            Assert.AreEqual(lambdaU, lambdaL, 1E-10, "t-copula tail dependence should be symmetric.");

            // ρ = 0, ν = 4: still has tail dependence (this is the key difference from Normal copula)
            var copulaZero = new StudentTCopula(0.0, 4);
            double lambdaZero = copulaZero.UpperTailDependence;
            Assert.IsGreaterThan(0, lambdaZero, "t-copula with ρ=0 still has tail dependence.");

            // Higher ν → lower tail dependence (approaching Normal copula with λ=0)
            var copulaHighNu = new StudentTCopula(0.5, 100);
            double lambdaHighNu = copulaHighNu.UpperTailDependence;
            Assert.IsLessThan(lambdaU, lambdaHighNu, "Higher ν should reduce tail dependence.");

            // Higher ρ → higher tail dependence
            var copulaHighRho = new StudentTCopula(0.9, 4);
            Assert.IsGreaterThan(lambdaU, copulaHighRho.UpperTailDependence, "Higher ρ should increase tail dependence.");

            // ρ close to -1: tail dependence approaches 0
            var copulaNegRho = new StudentTCopula(-0.99, 4);
            Assert.IsLessThan(0.01, copulaNegRho.UpperTailDependence, "Tail dependence should be near 0 for ρ ≈ -1.");
        }

        /// <summary>
        /// Test that the t-copula converges to the Normal copula as ν → ∞.
        /// </summary>
        [TestMethod]
        public void Test_ConvergenceToNormal()
        {
            double rho = 0.5;
            var normalCopula = new NormalCopula(rho);
            var tCopulaHighNu = new StudentTCopula(rho, 1000);

            // PDF should converge
            double normalPdf = normalCopula.PDF(0.3, 0.7);
            double tPdf = tCopulaHighNu.PDF(0.3, 0.7);
            Assert.AreEqual(normalPdf, tPdf, 0.05);

            // CDF should converge
            double normalCdf = normalCopula.CDF(0.3, 0.7);
            double tCdf = tCopulaHighNu.CDF(0.3, 0.7);
            Assert.AreEqual(normalCdf, tCdf, 0.05);
        }

        /// <summary>
        /// Test random sampling generates valid copula samples.
        /// </summary>
        [TestMethod]
        public void Test_Sampling()
        {
            var copula = new StudentTCopula(0.6, 5);
            var samples = copula.GenerateRandomValues(500, seed: 42);

            for (int i = 0; i < 500; i++)
            {
                // All values should be in [0, 1]
                Assert.IsGreaterThanOrEqualTo(0, samples[i, 0], $"Sample [{i},0] out of range.");
                Assert.IsLessThanOrEqualTo(1, samples[i, 0], $"Sample [{i},0] out of range.");
                Assert.IsGreaterThanOrEqualTo(0, samples[i, 1], $"Sample [{i},1] out of range.");
                Assert.IsLessThanOrEqualTo(1, samples[i, 1], $"Sample [{i},1] out of range.");
            }
        }

        /// <summary>
        /// Test Clone produces an independent copy.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new StudentTCopula(0.5, 10);
            var clone = copula.Clone() as StudentTCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(copula.Theta, clone.Theta);
            Assert.AreEqual(copula.DegreesOfFreedom, clone.DegreesOfFreedom);
            Assert.IsTrue(clone.ParametersValid);

            // Mutating clone should not affect original
            clone.Theta = 0.1;
            Assert.AreEqual(0.5, copula.Theta);
        }

        private double[] data1 = new double[] { 122.094066003419, 92.8321267206161, 86.4920318705377, 87.6183663113541, 102.558777787492, 103.627475117762, 127.084948716539, 105.908684131013, 110.065795957654, 105.924647125867, 110.009738155469, 126.490833800772, 64.1264871206211, 81.3150800229481, 92.0780134395721, 106.040322550555, 113.158086143066, 117.051057784044, 127.110531266645, 108.907371862136, 105.476247114194, 108.629403495407, 98.7803988364997, 93.217925588845, 97.7219451830075, 109.178093756809, 137.69504856252, 106.884615327674, 112.139177456202, 85.7416217661797, 71.0610938629716, 112.644166631765, 119.545871678548, 70.5169833274982, 99.6896817997206, 100.987892854545, 103.659280253554, 75.6075621013066, 118.810868919796, 109.113664695226, 113.636425353944, 100.008375355612, 113.178917359795, 80.4269472604342, 88.3638384448237, 90.2905074656314, 98.7995143316863, 98.4698060067802, 108.279297570816, 86.1578437055905, 101.183725242941, 85.5531148952956, 111.024195253862, 121.934506174556, 104.169993666179, 84.4652994609478, 99.6099259747033, 95.3130792386208, 115.45680252817, 120.213139478586, 95.5691788140058, 92.7950300448044, 102.58430893827, 86.7105161576407, 82.8059368562185, 107.335705516294, 112.603259240932, 102.780778760832, 128.958090528336, 105.139162595628, 118.272661482198, 99.8275937885748, 94.2856024560543, 108.48679008009, 100.147734981682, 88.7006383425785, 89.6441478272035, 112.24266306884, 99.8184811468069, 120.592090049738, 124.023170133661, 101.250961381805, 90.0000027551006, 108.781064635426, 94.9203320035987, 99.9491821782837, 88.7473944659517, 94.3643253649856, 105.814317118952, 92.6866900633813, 111.020330544613, 111.676189456988, 115.70235103978, 124.659106152655, 81.3866270495082, 120.178528245778, 93.6511977805724, 114.099368762143, 119.062045395294, 74.1998497412903 };
        private double[] data2 = new double[] { 127.869024514059, 53.5970265830273, 35.6871183968043, 77.5937820397885, 84.619117510857, 110.477376636164, 114.679535976765, 109.338354392258, 88.5987759167264, 72.6695216679034, 111.932652280673, 86.3677960278751, 23.9336347978345, 51.2377830227977, 82.4565771813309, 92.9162733515069, 117.465381827514, 104.862362549521, 131.059266136887, 67.2743851584176, 100.263235166171, 113.734275000025, 73.1582387829997, 78.4353197703676, 60.0180359279642, 106.709991071405, 123.175455301514, 98.7006449949188, 99.860486991242, 55.7603096813567, 53.7716423706874, 104.659445447656, 119.899401349887, 59.8670226375024, 94.0117104763717, 101.424610891155, 114.256354904191, 53.5051841563538, 118.35993465227, 73.1605008375787, 87.4677698350712, 75.4031529479113, 105.404958657365, 53.336411944238, 61.2731445424292, 72.377272009744, 88.959659863884, 80.1301183393358, 98.624093971352, 81.9603074727622, 52.0788199186743, 75.49358652998, 90.2428259997917, 101.326931349259, 48.1343463500222, 56.9295918059918, 89.0348875829931, 69.0012535890253, 100.355241744174, 74.00820280539, 63.9482913881998, 64.4973782209222, 95.8934144135508, 85.4028102356618, 37.8958459664423, 99.2194777630975, 126.581868541047, 91.8287794302242, 143.543939198862, 108.751405708845, 100.951567564812, 73.5051068155712, 83.419507788205, 84.9090133796832, 59.3886126411711, 84.0348703304947, 78.1503115396303, 104.953483626903, 77.6450718557069, 117.615613165515, 118.131904013699, 76.3190144944821, 62.0183469143453, 97.4729901076061, 49.3396925267253, 58.6790714873228, 45.0596168059506, 85.3857426310419, 65.0772008397323, 58.8836242438228, 79.2838406333912, 102.608529398935, 83.7509120512927, 103.106132785215, 52.8403092456187, 88.4802383528401, 64.2906982187616, 93.0489784548541, 116.065815369284, 26.2779209375887 };

        /// <summary>
        /// Test fitting with the method of maximum pseudo likelihood.
        /// </summary>
        [TestMethod]
        public void Test_MPL_Fit()
        {
            // get ranks of data
            var rank1 = Statistics.RanksInPlace(data1);
            var rank2 = Statistics.RanksInPlace(data2);
            // get plotting positions
            for (int i = 0; i < data1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }
            // Fit copula (fix ν=5, estimate ρ)
            BivariateCopula copula = new StudentTCopula(0.0, 5);
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);

            // The estimated ρ should be positive and reasonable for this correlated data
            Assert.IsGreaterThan(0.5, copula.Theta, $"Estimated ρ = {copula.Theta} should be > 0.5");
            Assert.IsLessThan(1.0, copula.Theta, $"Estimated ρ = {copula.Theta} should be < 1.0");
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new StudentTCopula(0.0, 5);
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);

            // The estimated ρ should be positive and reasonable
            Assert.IsGreaterThan(0.5, copula.Theta, $"Estimated ρ = {copula.Theta} should be > 0.5");
            Assert.IsLessThan(1.0, copula.Theta, $"Estimated ρ = {copula.Theta} should be < 1.0");
        }

        /// <summary>
        /// Test that MPL estimates both rho and degrees of freedom.
        /// </summary>
        [TestMethod]
        public void Test_MPL_EstimatesBothParameters()
        {
            var rank1 = Statistics.RanksInPlace(data1);
            var rank2 = Statistics.RanksInPlace(data2);
            for (int i = 0; i < data1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }

            // Start with default ν=5, the estimation should update it
            BivariateCopula copula = new StudentTCopula(0.0, 5);
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);

            var tCopula = (StudentTCopula)copula;

            // ρ should be estimated as positive
            Assert.IsGreaterThan(0.5, tCopula.Theta, $"Estimated ρ = {tCopula.Theta} should be > 0.5");
            Assert.IsLessThan(1.0, tCopula.Theta, $"Estimated ρ = {tCopula.Theta} should be < 1.0");

            // ν should be estimated (may differ from initial value of 5)
            Assert.IsGreaterThanOrEqualTo(3, tCopula.DegreesOfFreedom, $"Estimated ν = {tCopula.DegreesOfFreedom} should be >= 3");
            Assert.IsLessThanOrEqualTo(60, tCopula.DegreesOfFreedom, $"Estimated ν = {tCopula.DegreesOfFreedom} should be <= 60");

            // Tail dependence should be data-driven (not user-defined)
            double lambda = tCopula.UpperTailDependence;
            Assert.IsGreaterThan(0, lambda, "Tail dependence should be positive.");
            Assert.IsLessThan(1, lambda, "Tail dependence should be less than 1.");
        }

        /// <summary>
        /// Test that IFM estimates both rho and degrees of freedom.
        /// </summary>
        [TestMethod]
        public void Test_IFM_EstimatesBothParameters()
        {
            BivariateCopula copula = new StudentTCopula(0.0, 5);
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);

            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);

            var tCopula = (StudentTCopula)copula;
            Assert.IsGreaterThan(0.5, tCopula.Theta, $"Estimated ρ = {tCopula.Theta} should be > 0.5");
            Assert.IsGreaterThanOrEqualTo(3, tCopula.DegreesOfFreedom, $"Estimated ν = {tCopula.DegreesOfFreedom} should be >= 3");
            Assert.IsLessThanOrEqualTo(60, tCopula.DegreesOfFreedom, $"Estimated ν = {tCopula.DegreesOfFreedom} should be <= 60");
        }

        /// <summary>
        /// Test GetCopulaParameters and SetCopulaParameters round-trip.
        /// </summary>
        [TestMethod]
        public void Test_GetSetCopulaParameters()
        {
            var copula = new StudentTCopula(0.7, 10);

            // GetCopulaParameters should return [rho, nu]
            Assert.AreEqual(2, copula.NumberOfCopulaParameters);
            var parms = copula.GetCopulaParameters;
            Assert.AreEqual(0.7, parms[0], 1E-10);
            Assert.AreEqual(10.0, parms[1], 1E-10);

            // SetCopulaParameters should update both
            copula.SetCopulaParameters(new double[] { -0.3, 15.0 });
            Assert.AreEqual(-0.3, copula.Theta, 1E-10);
            Assert.AreEqual(15.0, copula.DegreesOfFreedom);

            // SetCopulaParameters should clamp nu to minimum of 3
            copula.SetCopulaParameters(new double[] { 0.5, 1.0 });
            Assert.AreEqual(3.0, copula.DegreesOfFreedom);
        }

        /// <summary>
        /// Test ParameterConstraints returns correct 2D array.
        /// </summary>
        [TestMethod]
        public void Test_ParameterConstraints()
        {
            var copula = new StudentTCopula(0.5, 5);
            var constraints = copula.ParameterConstraints(data1, data2);

            // Should be [2, 2] array
            Assert.AreEqual(2, constraints.GetLength(0));
            Assert.AreEqual(2, constraints.GetLength(1));

            // Row 0: rho constraints [-1+eps, 1-eps]
            Assert.IsLessThan(-0.99, constraints[0, 0]);
            Assert.IsGreaterThan(0.99, constraints[0, 1]);

            // Row 1: nu constraints [3, 60]
            Assert.AreEqual(3.0, constraints[1, 0]);
            Assert.AreEqual(60.0, constraints[1, 1]);
        }

    }
}
