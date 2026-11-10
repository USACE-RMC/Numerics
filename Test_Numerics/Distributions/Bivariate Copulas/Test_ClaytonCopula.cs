using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Clayton Copula. Reference values come from the R 'copula' package for the original methods and from pyvinecopulib 0.7.6 and mpmath 1.4.1 for the conditional-distribution methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_ClaytonCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new ClaytonCopula(10);
            Assert.AreEqual(1.3113E-05, copula.PDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(2.387424E-11, copula.PDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(1.25, copula.PDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(0.00731365, copula.PDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new ClaytonCopula(10);
            Assert.AreEqual(1.3113E-05, copula.PDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(2.387424E-11, copula.PDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(1.25, copula.PDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(0.00731365, copula.PDF(0.8, 0.2), 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new ClaytonCopula(10);
            Assert.AreEqual(0.2, copula.CDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(0.2, copula.CDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(0.1167184, copula.CDF(0.2, 0.8), 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(0.1999737, copula.CDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new ClaytonCopula(10);
            Assert.AreEqual(0.2, copula.CDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(20);
            Assert.AreEqual(0.2, copula.CDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(-0.5);
            Assert.AreEqual(0.1167184, copula.CDF(0.8, 0.2), 1E-6);

            copula = new ClaytonCopula(5);
            Assert.AreEqual(0.1999737, copula.CDF(0.8, 0.2), 1E-6);

        }

        private double[] data1 = new double[] { 135.852695757514, 104.082298038859, 108.737560538974, 99.2685553595344, 134.734517358047, 90.96671487189, 77.2933084429689, 115.40695335195, 108.98103448672, 79.030140378764, 119.344309239779, 114.057958622873, 121.376980352823, 102.545934287125, 116.171688724822, 90.2737860301091, 118.384838925781, 78.5650043184329, 122.502611362073, 119.21456701156, 122.901928011128, 79.7224942508733, 103.26308186608, 134.856147348969, 87.2443704188882, 107.484484922021, 102.369642440458, 106.295354390286, 103.333259725228, 114.944097101105, 101.947395165995, 99.3368989207062, 107.600964955177, 63.8796678112047, 122.606004825966, 96.2971006035691, 101.928179675011, 98.6811421254304, 136.414461816426, 94.4488702685262, 118.178109208696, 99.3447054496099, 73.6895736287131, 111.543756772209, 94.1765631568512, 108.378057968653, 82.0682187268786, 93.2316741233556, 92.1848572845123, 101.62266736115, 90.9160104057016, 97.6045250944967, 123.668316008262, 91.5790391025517, 118.507497005719, 116.028361864995, 106.431274472801, 88.1845431623773, 100.593879884171, 106.477111540899, 125.501606137489, 120.558072482172, 93.9426575106378, 61.491148690679, 81.7991169923448, 70.0617929783816, 124.844586043561, 97.6329727079425, 98.4538343555991, 78.7796746466831, 93.4697095163469, 94.4204666650014, 121.959329690142, 86.5743545885648, 119.317279401324, 78.9585491320436, 115.814312639539, 102.440872193275, 124.220180616412, 88.2386416291415, 111.637544083286, 74.24445197777, 139.719855097958, 111.842726007709, 106.593685525292, 95.7079419062174, 100.768445085655, 89.9934014730302, 89.7346792611329, 90.5836289523697, 95.6118805500478, 98.0110615835651, 104.685614771651, 105.840599667288, 110.492236632241, 115.285805520348, 112.406028227278, 97.3674746819745, 77.7088711413582, 102.843691800664 };
        private double[] data2 = new double[] { 106.247204356978, 119.634538634978, 72.5075040420413, 96.4624345242022, 93.7887427218662, 88.0710430691384, 54.3965822738992, 133.939640413503, 147.313192863588, 39.596862389244, 113.685219500925, 65.661791172481, 125.494512719209, 122.353356065013, 86.4850623431384, 56.2427898884044, 112.239490666668, 48.2002798009377, 112.836791609701, 70.306705761088, 59.4614291231979, 44.3647320315126, 91.7637292129578, 95.6179940736635, 60.4815479415048, 58.3786398648642, 94.1724393628596, 87.9631903747507, 91.5086027340312, 114.720801315278, 89.8113482015865, 82.2111665540629, 74.6453503599532, 43.9664396683432, 74.8202107513964, 98.3796818314241, 105.279101457561, 97.8452040931319, 99.7804654911693, 75.7974672384221, 139.913874523229, 81.687017386353, 34.7670214588405, 68.5923527448465, 93.0469463410792, 90.1397446225404, 51.1775088910909, 55.8146471548922, 77.1273200963725, 66.6745681621479, 95.5181914247239, 61.5430213312507, 76.4727428009896, 64.5965344713655, 75.1297368849306, 117.072661743052, 79.6250066166689, 37.2231394673498, 71.1212261596989, 76.6759488420322, 107.428580054692, 110.536313019896, 87.680895839141, 30.0011658085059, 48.6775096662048, 21.1838108219624, 95.5909448097503, 83.8207210333862, 116.670732017431, 40.4019864082416, 91.469151821309, 70.7404242103512, 45.6819094938675, 86.6712212708917, 121.431171230518, 38.4195748844222, 143.699998594758, 64.2498614632223, 75.9822870091017, 59.6123708381618, 113.933833528108, 42.0682009456935, 73.2020173163658, 51.1044236832251, 98.8788896427372, 46.7142903461546, 80.8155805864023, 87.4214895184672, 62.741943393305, 80.361655652615, 70.7251173704377, 89.8903637671805, 74.0531075745106, 72.4690834186899, 107.394378394182, 89.3201273021353, 54.0735001346141, 66.3099470947684, 55.5401474570749, 63.1916628514916 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new ClaytonCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(1.326613, copula.Theta, 1E-4);
        }

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
            // Fit copula
            BivariateCopula copula = new ClaytonCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(1.534016, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new ClaytonCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(1.485167, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new ClaytonCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(1.512075, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(102.4226, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(16.70995, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(80.54073, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(26.11687, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }

        /// <summary>
        /// Test the tail dependence coefficients.
        /// Clayton copula: λ_L = 2^(-1/θ), λ_U = 0.
        /// </summary>
        [TestMethod]
        public void Test_TailDependence()
        {
            // θ = 2: λ_L = 2^(-1/2) ≈ 0.70711
            var copula = new ClaytonCopula(2.0);
            Assert.AreEqual(0.0, copula.UpperTailDependence, 1E-10, "Clayton copula should have no upper tail dependence.");
            Assert.AreEqual(Math.Pow(2.0, -0.5), copula.LowerTailDependence, 1E-6, "λ_L should equal 2^(-1/θ) for θ=2.");

            // θ = 1: λ_L = 2^(-1) = 0.5
            var copula1 = new ClaytonCopula(1.0);
            Assert.AreEqual(0.5, copula1.LowerTailDependence, 1E-10, "λ_L should equal 0.5 for θ=1.");

            // θ → ∞: λ_L → 1
            var copulaLarge = new ClaytonCopula(100.0);
            Assert.IsGreaterThan(0.99, copulaLarge.LowerTailDependence, "λ_L should approach 1 for large θ.");

            // θ ≤ 0: λ_L = 0 (no dependence or negative dependence)
            var copulaZero = new ClaytonCopula(0.0);
            Assert.AreEqual(0.0, copulaZero.LowerTailDependence, 1E-10, "λ_L should be 0 for θ ≤ 0.");
        }

        /// <summary>
        /// Test that ParametersValid tracks the dependency parameter's valid range.
        /// ValidateParameter returns null for a valid dependence parameter and a non-null
        /// exception otherwise, so ParametersValid follows the parameter.
        /// </summary>
        [TestMethod]
        public void Test_ParametersValid()
        {
            var copula = new ClaytonCopula(2.0);
            Assert.IsTrue(copula.ParametersValid);

            copula.Theta = -2.0; // below the minimum of -1
            Assert.IsFalse(copula.ParametersValid);

            copula.Theta = 3.0;
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test Clone produces an independent copy with deep-copied marginals.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new ClaytonCopula(2.0, new Normal(100, 10), new Gumbel(50, 5));
            var clone = copula.Clone() as ClaytonCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(copula.Theta, clone.Theta);
            Assert.IsTrue(clone.ParametersValid);

            // Marginals are deep-copied (distributions memoize lazily, so clones must
            // not share marginal instances), with identical quantiles
            Assert.AreNotSame(copula.MarginalDistributionX, clone.MarginalDistributionX);
            Assert.AreNotSame(copula.MarginalDistributionY, clone.MarginalDistributionY);
            Assert.AreEqual(copula.MarginalDistributionY.InverseCDF(0.9), clone.MarginalDistributionY.InverseCDF(0.9));

            // Mutating the clone should not affect the original
            clone.Theta = 4.0;
            Assert.AreEqual(2.0, copula.Theta);
        }

        /// <summary>
        /// Test the forward conditional CDF (h-function), h(v|u) = ∂C(u,v)/∂u.
        /// Reference values were computed with the Python package pyvinecopulib 0.7.6
        /// (Bicop.hfunc1, family clayton), cross-checked at generation against a 50-digit
        /// numerical ∂C/∂u of the Clayton CDF computed with mpmath 1.4.1 (agreement ≤ 1E-15).
        /// The analytic conditional is also cross-checked against a central finite difference
        /// of the closed-form CDF (ε = 1E-6; the difference-quotient noise floor is ~1E-9, so
        /// 1E-7 is asserted) and pinned at the exact upper edge h(1|u) = 1.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDF()
        {
            var copula = new ClaytonCopula(2d);
            Assert.AreEqual(0.8743161176077272, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.06882371771256163, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.9991210147197616, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(0.0001713170464197122, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            copula = new ClaytonCopula(5d);
            Assert.AreEqual(0.985754642959101, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.006108127860986364, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.9999997399342684, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(2.9401186465615086E-08, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            foreach (double theta in new[] { 2d, 5d })
            {
                var c = new ClaytonCopula(theta);
                foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    foreach (double v in new[] { 0.1, 0.5, 0.9 })
                    {
                        double fd = (c.CDF(u + 1E-6, v) - c.CDF(u - 1E-6, v)) / 2E-6;
                        Assert.AreEqual(fd, c.ConditionalCDF(u, v), 1E-7);
                    }
                    Assert.AreEqual(1d, c.ConditionalCDF(u, 1d), 1E-12);
                }
            }
        }

        /// <summary>
        /// Test the scalar inverse conditional CDF. Reference values were computed with the
        /// Python package pyvinecopulib 0.7.6 (Bicop.hinv1, family clayton). Also verifies the
        /// closed-form round trip h(InverseConditionalCDF(u, t) | u) = t (both directions are
        /// closed forms, so 1E-10 is asserted) and that InverseCDF(u, t)[1] recomposes the
        /// scalar bit-for-bit.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDF()
        {
            var copula = new ClaytonCopula(2d);
            Assert.AreEqual(0.5010908594248752, copula.InverseConditionalCDF(0.3, 0.7), 1E-9);
            Assert.AreEqual(0.3359223320941678, copula.InverseConditionalCDF(0.9, 0.05), 1E-9);

            copula = new ClaytonCopula(5d);
            Assert.AreEqual(0.37039687871242977, copula.InverseConditionalCDF(0.3, 0.7), 1E-9);
            Assert.AreEqual(0.5500280933001612, copula.InverseConditionalCDF(0.9, 0.05), 1E-9);

            foreach (double theta in new[] { 2d, 5d })
            {
                var c = new ClaytonCopula(theta);
                foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    foreach (double t in new[] { 0.1, 0.5, 0.9 })
                    {
                        double v = c.InverseConditionalCDF(u, t);
                        Assert.AreEqual(t, c.ConditionalCDF(u, v), 1E-10);
                        Assert.AreEqual(v, c.InverseCDF(u, t)[1], 0d);
                    }
                }
            }
        }

        /// <summary>
        /// LogPDF and PDF against the closed-form log density evaluated in mpmath 1.3.0 at 120 digits, for θ from -0.9 to 100
        /// and u, v from 1E-12 to 1 - 1E-12, where the direct density underflows or overflows. The log density is within 1E-12
        /// (relative above 1 in magnitude) and PDF is its exponential.
        /// </summary>
        [TestMethod]
        public void Test_LogPDFMatchesMpmath()
        {
            double[,] cases =
            {
            { -0.9, 1e-12, 1e-12, double.NegativeInfinity },
            { -0.9, 1e-12, 0.5, double.NegativeInfinity },
            { -0.9, 1e-12, 0.999999999999, 22.61729917857869 },
            { -0.9, 0.01, 0.99, 2.5894456155757988 },
            { -0.9, 0.2, 0.8, 0.4921553035203186 },
            { -0.9, 0.5, 0.5, 0.1775914309458633 },
            { -0.9, 0.999999999999, 0.5, -1.6787526304885023 },
            { -0.9, 0.999999999999, 0.999999999999, -2.302585092992246 },
            { -0.5, 1e-12, 1e-12, double.NegativeInfinity },
            { -0.5, 1e-12, 0.5, double.NegativeInfinity },
            { -0.5, 1e-12, 0.999999999999, 13.12236337740483 },
            { -0.5, 0.01, 0.99, 1.614463080360851 },
            { -0.5, 0.2, 0.8, 0.2231435513142097 },
            { -0.5, 0.5, 0.5, 0.0 },
            { -0.5, 0.999999999999, 0.5, -0.34657359027947265 },
            { -0.5, 0.999999999999, 0.999999999999, -0.6931471805589453 },
            { 2.0, 1e-12, 1e-12, 26.996765453196794 },
            { 2.0, 1e-12, 0.5, -52.08398840150915 },
            { 2.0, 1e-12, 0.999999999999, -54.16342994318599 },
            { 2.0, 0.01, 0.99, -8.081582151755068 },
            { 2.0, 0.2, 0.8, -1.5064594045945112 },
            { 2.0, 0.5, 0.5, 0.3927199993894983 },
            { 2.0, 0.999999999999, 0.5, -0.28768207245003097 },
            { 2.0, 0.999999999999, 0.999999999999, 1.0986122886641099 },
            { 100.0, 1e-12, 1e-12, 30.852915799844318 },
            { 100.0, 1e-12, 0.5, -2688.479125839459 },
            { 100.0, 1e-12, 0.999999999999, -2758.4869910759126 },
            { 100.0, 0.01, 0.99, -454.88681416076423 },
            { 100.0, 0.2, 0.8, -133.7911720438336 },
            { 100.0, 0.5, 0.5, 3.9150418644757146 },
            { 100.0, 0.999999999999, 0.5, -64.69959753905228 },
            { 100.0, 0.999999999999, 0.999999999999, 4.615120516641264 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new ClaytonCopula(cases[i, 0]);
                double u = cases[i, 1], v = cases[i, 2], expected = cases[i, 3];
                double actual = copula.LogPDF(u, v);
                if (double.IsNegativeInfinity(expected))
                {
                    Assert.IsTrue(double.IsNegativeInfinity(actual), $"case {i}");
                    Assert.AreEqual(0d, copula.PDF(u, v), 0d, $"case {i}");
                    continue;
                }
                Assert.AreEqual(expected, actual, 1E-12 * Math.Max(1d, Math.Abs(expected)), $"case {i}: theta={cases[i, 0]}, u={u:R}, v={v:R}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Math.Exp(actual)), BitConverter.DoubleToInt64Bits(copula.PDF(u, v)), $"case {i}");
            }
        }

        /// <summary>
        /// t = 1 returns 1 at every u, including u = 0 where the closed form multiplies an infinite power by zero; the other
        /// edges keep the closed-form limits.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDFEdges()
        {
            var copula = new ClaytonCopula(2d);
            Assert.AreEqual(1d, copula.InverseConditionalCDF(0d, 1d), 0d);
            Assert.AreEqual(1d, copula.InverseConditionalCDF(1E-300, 1d), 0d);
            Assert.AreEqual(0d, copula.InverseConditionalCDF(0d, 0.5), 0d);
            Assert.AreEqual(Math.Pow(0.5, 1d / 3d), copula.InverseConditionalCDF(1d, 0.5), 1E-15);
            Assert.AreEqual(1d, new ClaytonCopula(-0.5).InverseConditionalCDF(0d, 0.5), 0d);
        }

        /// <summary>
        /// θ = 0 is the independence limit: PDF 1, LogPDF 0, CDF uv, ConditionalCDF v and InverseConditionalCDF t, on the
        /// closed unit square.
        /// </summary>
        [TestMethod]
        public void Test_ThetaZeroIsIndependence()
        {
            var copula = new ClaytonCopula(0d);
            double[] grid = { 0d, 1E-12, 0.01, 0.3, 0.5, 0.7, 0.99, 1d - 1E-12, 1d };
            foreach (double u in grid)
            {
                foreach (double v in grid)
                {
                    Assert.AreEqual(1d, copula.PDF(u, v), 0d);
                    Assert.AreEqual(0d, copula.LogPDF(u, v), 0d);
                    Assert.AreEqual(u * v, copula.CDF(u, v), 0d);
                    Assert.AreEqual(v, copula.ConditionalCDF(u, v), 0d);
                    Assert.AreEqual(v, copula.InverseConditionalCDF(u, v), 0d);
                }
            }
        }

        /// <summary>
        /// The members are continuous through θ = 0: at θ = ±1E-12 and ±1E-9 they agree with mpmath 1.3.0 at 80 digits to
        /// 5E-14 (relative for CDF, ConditionalCDF and InverseConditionalCDF, absolute for LogPDF), and so with the θ = 0
        /// limits to within θ times a bounded factor.
        /// </summary>
        [TestMethod]
        public void Test_ContinuousThroughThetaZero()
        {
            // theta, u, v, then mpmath CDF(u, v), LogPDF(u, v), ConditionalCDF(u, v) and InverseConditionalCDF(u, v).
            double[,] cases =
            {
            { -1e-09, 1e-12, 0.5, 4.999999904238178e-13, 8.171803680557638e-09, 0.49999999077039137, 0.5000000092296085 },
            { -1e-09, 0.3, 0.7, 0.20999999990982032, 1.3122081580171854e-10, 0.6999999999490736, 0.7000000000509263 },
            { -1e-09, 0.99, 0.01, 0.009899999999541793, 3.5689370142956436e-09, 0.010000000045588868, 0.009999999954411134 },
            { -1e-09, 0.5, 1e-06, 4.999999952119089e-07, 3.932475496458715e-09, 1.0000000042393283e-06, 9.999999957606716e-07 },
            { -1e-12, 1e-12, 0.5, 4.999999999904238e-13, 8.171803913756873e-12, 0.4999999999907704, 0.5000000000092296 },
            { -1e-12, 0.3, 0.7, 0.2099999999999098, 1.3122081577798051e-13, 0.699999999999949, 0.7000000000000509 },
            { -1e-12, 0.99, 0.01, 0.009899999999999541, 3.5689370148093665e-12, 0.010000000000045589, 0.009999999999954412 },
            { -1e-12, 0.5, 1e-06, 4.999999999952119e-07, 3.9324755472243095e-12, 1.0000000000042394e-06, 9.999999999957606e-07 },
            { 1e-12, 1e-12, 0.5, 5.000000000095762e-13, -8.17180391422374e-12, 0.5000000000092296, 0.4999999999907704 },
            { 1e-12, 0.3, 0.7, 0.21000000000009017, -1.3122081577793298e-13, 0.7000000000000509, 0.699999999999949 },
            { 1e-12, 0.99, 0.01, 0.009900000000000459, -3.568937014810395e-12, 0.009999999999954412, 0.010000000000045589 },
            { 1e-12, 0.5, 1e-06, 5.00000000004788e-07, -3.932475547325942e-12, 9.999999999957606e-07, 1.0000000000042394e-06 },
            { 1e-09, 1e-12, 0.5, 5.000000095761821e-13, -8.171804147422972e-09, 0.5000000092296085, 0.49999999077039137 },
            { 1e-09, 0.3, 0.7, 0.21000000009017963, -1.3122081575419498e-10, 0.7000000000509263, 0.6999999999490736 },
            { 1e-09, 0.99, 0.01, 0.009900000000458207, -3.5689370153241182e-09, 0.009999999954411134, 0.010000000045588866 },
            { 1e-09, 0.5, 1e-06, 5.000000047880911e-07, -3.932475598091537e-09, 9.999999957606716e-07, 1.0000000042393285e-06 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new ClaytonCopula(cases[i, 0]);
                double u = cases[i, 1], v = cases[i, 2];
                Assert.AreEqual(cases[i, 3], copula.CDF(u, v), 5E-14 * cases[i, 3], $"case {i}: CDF");
                Assert.AreEqual(cases[i, 4], copula.LogPDF(u, v), 5E-14, $"case {i}: LogPDF");
                Assert.AreEqual(cases[i, 5], copula.ConditionalCDF(u, v), 5E-14 * cases[i, 5], $"case {i}: ConditionalCDF");
                Assert.AreEqual(cases[i, 6], copula.InverseConditionalCDF(u, v), 5E-14 * cases[i, 6], $"case {i}: InverseConditionalCDF");
                // Within 1E-12 of independence the values sit within about 30·|θ| of the limits.
                Assert.AreEqual(u * v, copula.CDF(u, v), 1E-7 * u * v);
                Assert.AreEqual(0d, copula.LogPDF(u, v), 1E-7);
            }
        }

        /// <summary>
        /// Every generator member against mpmath 1.3.0 at 150 digits, for φ(t) = (t^(−θ) − 1)/θ (Nelsen 2006, Table 4.1), with the independence limit −ln t at θ = 0: Generator, GeneratorInverse (at
        /// the rounded generator value), GeneratorPrime, GeneratorPrime2 and GeneratorPrimeInverse (at the rounded slope),
        /// on t from 1E-12 to 1 - 1E-12 and θ including the independence region, to 1E-13 relative.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorMembersMatchMpmath()
        {
            // theta, t, phi(t), phi^-1(phi), phi'(t), phi''(t), (phi')^-1(phi').
            double[,] cases =
            {
            { -0.9, 1e-12, 1.1111111110935012, 9.999944734812045e-13, -15.848931924611126, 1584893192461.1123, 9.999999999999996e-13 },
            { -0.9, 0.01, 1.093501186750432, 0.00999999999999999, -1.5848931924611134, 15.84893192461113, 0.009999999999999997 },
            { -0.9, 0.5, 0.5156814097020593, 0.5, -1.0717734625362931, 0.21435469250725858, 0.5000000000000001 },
            { -0.9, 0.99, 0.0100050184301842, 0.99, -1.0010055388008428, 0.10111167058594368, 0.989999999999999 },
            { -0.9, 0.999999999999, 9.999778782799286e-13, 0.999999999999, -1.0000000000001, 0.10000000000010997, 0.9999999999990008 },
            { -1e-12, 1e-12, 27.63102111554681, 9.999999999999998e-13, -999999999972.369, 9.99999999971369e+23, 1e-12 },
            { -1e-12, 0.01, 4.605170185977488, 0.009999999999999997, -99.99999999953948, 9999.999999943948, 0.01 },
            { -1e-12, 0.5, 0.693147180559705, 0.5, -1.9999999999986138, 3.999999999993227, 0.5 },
            { -1e-12, 0.99, 0.0100503358535014, 0.99, -1.010101010101, 1.0203040506060503, 0.99 },
            { -1e-12, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.0000000000009999, 1.0000000000009999, 0.9999999999990001 },
            { 0.0, 1e-12, 27.631021115928547, 1.000000000000001e-12, -1000000000000.0, 1e+24, 1e-12 },
            { 0.0, 0.01, 4.605170185988091, 0.010000000000000004, -100.0, 10000.0, 0.01 },
            { 0.0, 0.5, 0.6931471805599453, 0.5, -2.0, 4.0, 0.5 },
            { 0.0, 0.99, 0.01005033585350145, 0.99, -1.0101010101010102, 1.0203040506070808, 0.99 },
            { 0.0, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.000000000001, 1.000000000002, 0.9999999999989999 },
            { 1e-12, 1e-12, 27.631021116310286, 9.999999999999986e-13, -1000000000027.631, 1.000000000028631e+24, 1e-12 },
            { 1e-12, 0.01, 4.605170185998695, 0.010000000000000004, -100.00000000046052, 10000.000000056052, 0.01 },
            { 1e-12, 0.5, 0.6931471805601855, 0.5, -2.0000000000013864, 4.000000000006772, 0.49999999999999994 },
            { 1e-12, 0.99, 0.010050335853501501, 0.99, -1.0101010101010202, 1.0203040506081116, 0.9900000000000001 },
            { 1e-12, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.000000000001, 1.000000000003, 0.9999999999989999 },
            { 2.0, 1e-12, 5e+23, 1e-12, -1e+36, 3e+48, 1e-12 },
            { 2.0, 0.01, 4999.5, 0.01, -999999.9999999999, 300000000.0, 0.01 },
            { 2.0, 0.5, 1.5, 0.5, -8.0, 48.0, 0.5 },
            { 2.0, 0.99, 0.010152025303540464, 0.99, -1.0306101521283646, 3.1230610670556502, 0.99 },
            { 2.0, 0.999999999999, 9.999778782813784e-13, 0.999999999999, -1.000000000003, 3.0000000000119997, 0.999999999999 },
            { 50.0, 1e-12, double.PositiveInfinity, 0.0, double.NegativeInfinity, double.PositiveInfinity, 0.0 },
            { 50.0, 0.01, 1.999999999999998e+98, 0.01, -9.99999999999999e+101, 5.099999999999995e+105, 0.01 },
            { 50.0, 0.5, 22517998136852.46, 0.5, -2251799813685248.0, 2.296835809958953e+17, 0.5 },
            { 50.0, 0.99, 0.013057519728068095, 0.99, -1.6695717034377826, 86.00823926800699, 0.99 },
            { 50.0, 0.999999999999, 9.999778783053773e-13, 0.999999999999, -1.000000000050999, 51.000000002651944, 0.999999999999 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new ClaytonCopula(cases[i, 0]);
                double t = cases[i, 1];
                AssertClose(cases[i, 2], copula.Generator(t), $"case {i}: Generator");
                if (cases[i, 2] > 0d && !double.IsInfinity(cases[i, 2]))
                    AssertClose(cases[i, 3], copula.GeneratorInverse(cases[i, 2]), $"case {i}: GeneratorInverse");
                AssertClose(cases[i, 4], copula.GeneratorPrime(t), $"case {i}: GeneratorPrime");
                AssertClose(cases[i, 5], copula.GeneratorPrime2(t), $"case {i}: GeneratorPrime2");
                if (cases[i, 4] < 0d && !double.IsInfinity(cases[i, 4]))
                    AssertClose(cases[i, 6], copula.GeneratorPrimeInverse(cases[i, 4]), $"case {i}: GeneratorPrimeInverse");
            }
        }

        private static void AssertClose(double expected, double actual, string message)
        {
            if (double.IsInfinity(expected) || expected == 0d)
                Assert.AreEqual(expected, actual, 1E-300, message);
            else
                Assert.AreEqual(expected, actual, 1E-13 * Math.Abs(expected), message);
        }
    }
}
