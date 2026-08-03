using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Independence (product) copula, Π(u,v) = u·v. The copula is exact,
    /// so every statistical function is pinned against the closed form directly.
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
    public class Test_IndependenceCopula
    {
        /// <summary>
        /// Test the zero-parameter surface: no parameters exist, the parameter accessors are
        /// empty, and SetCopulaParameters ignores its input entirely.
        /// </summary>
        [TestMethod]
        public void Test_ZeroParameterSurface()
        {
            var copula = new IndependenceCopula();
            Assert.AreEqual(CopulaType.Independence, copula.Type);
            Assert.AreEqual("Independence", copula.DisplayName);
            Assert.AreEqual("Π", copula.ShortDisplayName);
            Assert.AreEqual(0, copula.NumberOfCopulaParameters);
            Assert.IsEmpty(copula.GetCopulaParameters);
            Assert.AreEqual(0, copula.ParameterToString.GetLength(0));
            Assert.AreEqual("", copula.ParameterNameShortForm);
            Assert.AreEqual(0d, copula.ThetaMinimum);
            Assert.AreEqual(0d, copula.ThetaMaximum);
            Assert.AreEqual(0, copula.ParameterConstraints(new[] { 1d, 2d }, new[] { 3d, 4d }).GetLength(0));

            // SetCopulaParameters is a no-op regardless of input
            copula.SetCopulaParameters(new double[] { 5d });
            Assert.AreEqual(0d, copula.Theta);
            Assert.IsTrue(copula.ParametersValid);
            copula.SetCopulaParameters(System.Array.Empty<double>());
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test that the copula is permanently valid: ValidateParameter always returns null
        /// and never throws, and no Theta assignment — not even NaN, which every
        /// parameterized family rejects — can invalidate it, because the dependency
        /// parameter carries no meaning for a zero-parameter copula.
        /// </summary>
        [TestMethod]
        public void Test_PermanentValidity()
        {
            var copula = new IndependenceCopula();
            Assert.IsTrue(copula.ParametersValid);
            Assert.IsNull(copula.ValidateParameter(double.NaN, false));
            Assert.IsNull(copula.ValidateParameter(double.PositiveInfinity, true));

            copula.Theta = double.NaN;
            Assert.IsTrue(copula.ParametersValid);
            copula.Theta = double.PositiveInfinity;
            Assert.IsTrue(copula.ParametersValid);

            // Statistical functions remain usable in every state
            Assert.AreEqual(0.35, copula.CDF(0.5, 0.7), 0d);
        }

        /// <summary>
        /// Test PDF: the product copula density is identically 1, so LogPDF is identically 0.
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new IndependenceCopula();
            Assert.AreEqual(1d, copula.PDF(0.2, 0.8), 0d);
            Assert.AreEqual(1d, copula.PDF(0.5, 0.5), 0d);
            Assert.AreEqual(1d, copula.PDF(0.99, 0.01), 0d);
            Assert.AreEqual(0d, copula.LogPDF(0.2, 0.8), 0d);
        }

        /// <summary>
        /// Test CDF: C(u,v) = u·v exactly, with the copula boundary conditions
        /// C(0,v) = 0, C(u,1) = u, and C(1,v) = v.
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new IndependenceCopula();
            Assert.AreEqual(0.2 * 0.8, copula.CDF(0.2, 0.8), 0d);
            Assert.AreEqual(0.25, copula.CDF(0.5, 0.5), 0d);
            Assert.AreEqual(0.045000000000000005, copula.CDF(0.9, 0.05), 0d);
            Assert.AreEqual(0d, copula.CDF(0d, 0.7), 0d);
            Assert.AreEqual(0.7, copula.CDF(0.7, 1d), 0d);
            Assert.AreEqual(0.7, copula.CDF(1d, 0.7), 0d);

            // Joint exceedance identities: OR = 1 - uv, AND = (1-u)(1-v)
            Assert.AreEqual(1d - 0.16, copula.ORJointExceedanceProbability(0.2, 0.8), 1E-15);
            Assert.AreEqual(0.8 * 0.2, copula.ANDJointExceedanceProbability(0.2, 0.8), 1E-15);
        }

        /// <summary>
        /// Test the forward conditional CDF: under independence h(v|u) = v exactly, for any
        /// conditioning value.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDF()
        {
            var copula = new IndependenceCopula();
            foreach (double u in new[] { 0.05, 0.3, 0.5, 0.7, 0.95 })
            {
                foreach (double v in new[] { 0d, 0.1, 0.5, 0.9, 1d })
                {
                    Assert.AreEqual(v, copula.ConditionalCDF(u, v), 0d);
                }
            }
        }

        /// <summary>
        /// Test conditional inversion: the scalar inverse is the identity in t, the array
        /// form returns the pair unchanged, and the array form recomposes the scalar
        /// bit-for-bit.
        /// </summary>
        [TestMethod]
        public void Test_InverseCDF()
        {
            var copula = new IndependenceCopula();
            foreach (double u in new[] { 0.05, 0.3, 0.5, 0.7, 0.95 })
            {
                foreach (double t in new[] { 0.1, 0.5, 0.9 })
                {
                    Assert.AreEqual(t, copula.InverseConditionalCDF(u, t), 0d);
                    var pair = copula.InverseCDF(u, t);
                    Assert.AreEqual(u, pair[0], 0d);
                    Assert.AreEqual(t, pair[1], 0d);
                }
            }
        }

        /// <summary>
        /// Test the tail dependence coefficients: independent variables have none.
        /// </summary>
        [TestMethod]
        public void Test_TailDependence()
        {
            var copula = new IndependenceCopula();
            Assert.AreEqual(0d, copula.UpperTailDependence, 0d);
            Assert.AreEqual(0d, copula.LowerTailDependence, 0d);
        }

        /// <summary>
        /// Test random generation: a seeded sample has the right shape, stays inside the unit
        /// square, reproduces bit-for-bit under the same seed, and carries no rank
        /// dependence. With n = 2,000 the standard error of Kendall's tau under independence
        /// is √(2(2n+5)/(9n(n−1))) ≈ 0.0149, so |τ| &lt; 0.05 is a ≈3.4σ bound on the fixed
        /// seed.
        /// </summary>
        [TestMethod]
        public void Test_GenerateRandomValues()
        {
            var copula = new IndependenceCopula();
            int n = 2000;
            var sample = copula.GenerateRandomValues(n, seed: 12345);
            Assert.AreEqual(n, sample.GetLength(0));
            Assert.AreEqual(2, sample.GetLength(1));

            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = sample[i, 0];
                y[i] = sample[i, 1];
                Assert.IsTrue(sample[i, 0] >= 0d && sample[i, 0] <= 1d);
                Assert.IsTrue(sample[i, 1] >= 0d && sample[i, 1] <= 1d);
            }

            var repeat = copula.GenerateRandomValues(n, seed: 12345);
            for (int i = 0; i < n; i++)
            {
                Assert.AreEqual(sample[i, 0], repeat[i, 0], 0d);
                Assert.AreEqual(sample[i, 1], repeat[i, 1], 0d);
            }

            double tau = Correlation.KendallsTau(x, y);
            Assert.AreEqual(0d, tau, 0.05);
        }

        /// <summary>
        /// Test Clone produces an independent copy with deep-copied marginals.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new IndependenceCopula(new Normal(100, 10), new Gumbel(50, 5));
            var clone = copula.Clone() as IndependenceCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(CopulaType.Independence, clone.Type);
            Assert.IsTrue(clone.ParametersValid);

            // Marginals are deep-copied (distributions memoize lazily, so clones must
            // not share marginal instances), with identical quantiles
            Assert.AreNotSame(copula.MarginalDistributionX, clone.MarginalDistributionX);
            Assert.AreNotSame(copula.MarginalDistributionY, clone.MarginalDistributionY);
            Assert.AreEqual(copula.MarginalDistributionY.InverseCDF(0.9), clone.MarginalDistributionY.InverseCDF(0.9));
        }

        /// <summary>
        /// Test that parameter estimation is a benign no-op for the zero-parameter copula
        /// under all three estimation methods: no exception is thrown (the full-likelihood
        /// path would otherwise demand estimable marginals) and the copula state is
        /// untouched.
        /// </summary>
        [TestMethod]
        public void Test_Estimation_NoOp()
        {
            var dataX = new double[] { 1d, 2d, 3d, 4d, 5d };
            var dataY = new double[] { 2d, 1d, 4d, 3d, 5d };
            foreach (CopulaEstimationMethod method in new[] { CopulaEstimationMethod.PseudoLikelihood, CopulaEstimationMethod.InferenceFromMargins, CopulaEstimationMethod.FullLikelihood })
            {
                BivariateCopula copula = new IndependenceCopula();
                BivariateCopulaEstimation.Estimate(ref copula, dataX, dataY, method);
                Assert.IsInstanceOfType(copula, typeof(IndependenceCopula));
                Assert.AreEqual(0d, copula.Theta);
                Assert.AreEqual(0, copula.NumberOfCopulaParameters);
                Assert.IsTrue(copula.ParametersValid);
                Assert.IsNull(copula.MarginalDistributionX);
                Assert.IsNull(copula.MarginalDistributionY);
            }
        }
    }
}
