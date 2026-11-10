using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.MachineLearning;
using System.Collections.Generic;
using System.Reflection;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;

namespace MachineLearning
{
    /// <summary>
    /// Unit tests for Gaussian Mixture Model (GMM) classification.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_GMM
    {
        /// <summary>
        /// Gaussian Mixture Model (GMM) is tested against the Iris dataset in R using the 'mclust' package.
        /// </summary>
        [TestMethod]
        public void Test_GMM_Iris()
        {
            var sepalLength = new double[] { 5.1, 4.9, 4.7, 4.6, 5, 5.4, 4.6, 5, 4.4, 4.9, 5.4, 4.8, 4.8, 4.3, 5.8, 5.7, 5.4, 5.1, 5.7, 5.1, 5.4, 5.1, 4.6, 5.1, 4.8, 5, 5, 5.2, 5.2, 4.7, 4.8, 5.4, 5.2, 5.5, 4.9, 5, 5.5, 4.9, 4.4, 5.1, 5, 4.5, 4.4, 5, 5.1, 4.8, 5.1, 4.6, 5.3, 5, 7, 6.4, 6.9, 5.5, 6.5, 5.7, 6.3, 4.9, 6.6, 5.2, 5, 5.9, 6, 6.1, 5.6, 6.7, 5.6, 5.8, 6.2, 5.6, 5.9, 6.1, 6.3, 6.1, 6.4, 6.6, 6.8, 6.7, 6, 5.7, 5.5, 5.5, 5.8, 6, 5.4, 6, 6.7, 6.3, 5.6, 5.5, 5.5, 6.1, 5.8, 5, 5.6, 5.7, 5.7, 6.2, 5.1, 5.7, 6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2, 6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6, 6.9, 5.6, 7.7, 6.3, 6.7, 7.2, 6.2, 6.1, 6.4, 7.2, 7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6, 6.9, 6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9 };
            var sepalWidth = new double[] { 3.5, 3, 3.2, 3.1, 3.6, 3.9, 3.4, 3.4, 2.9, 3.1, 3.7, 3.4, 3, 3, 4, 4.4, 3.9, 3.5, 3.8, 3.8, 3.4, 3.7, 3.6, 3.3, 3.4, 3, 3.4, 3.5, 3.4, 3.2, 3.1, 3.4, 4.1, 4.2, 3.1, 3.2, 3.5, 3.6, 3, 3.4, 3.5, 2.3, 3.2, 3.5, 3.8, 3, 3.8, 3.2, 3.7, 3.3, 3.2, 3.2, 3.1, 2.3, 2.8, 2.8, 3.3, 2.4, 2.9, 2.7, 2, 3, 2.2, 2.9, 2.9, 3.1, 3, 2.7, 2.2, 2.5, 3.2, 2.8, 2.5, 2.8, 2.9, 3, 2.8, 3, 2.9, 2.6, 2.4, 2.4, 2.7, 2.7, 3, 3.4, 3.1, 2.3, 3, 2.5, 2.6, 3, 2.6, 2.3, 2.7, 3, 2.9, 2.9, 2.5, 2.8, 3.3, 2.7, 3, 2.9, 3, 3, 2.5, 2.9, 2.5, 3.6, 3.2, 2.7, 3, 2.5, 2.8, 3.2, 3, 3.8, 2.6, 2.2, 3.2, 2.8, 2.8, 2.7, 3.3, 3.2, 2.8, 3, 2.8, 3, 2.8, 3.8, 2.8, 2.8, 2.6, 3, 3.4, 3.1, 3, 3.1, 3.1, 3.1, 2.7, 3.2, 3.3, 3, 2.5, 3, 3.4, 3 };
            var petalLength = new double[] { 1.4, 1.4, 1.3, 1.5, 1.4, 1.7, 1.4, 1.5, 1.4, 1.5, 1.5, 1.6, 1.4, 1.1, 1.2, 1.5, 1.3, 1.4, 1.7, 1.5, 1.7, 1.5, 1, 1.7, 1.9, 1.6, 1.6, 1.5, 1.4, 1.6, 1.6, 1.5, 1.5, 1.4, 1.5, 1.2, 1.3, 1.4, 1.3, 1.5, 1.3, 1.3, 1.3, 1.6, 1.9, 1.4, 1.6, 1.4, 1.5, 1.4, 4.7, 4.5, 4.9, 4, 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4, 4.7, 3.6, 4.4, 4.5, 4.1, 4.5, 3.9, 4.8, 4, 4.9, 4.7, 4.3, 4.4, 4.8, 5, 4.5, 3.5, 3.8, 3.7, 3.9, 5.1, 4.5, 4.5, 4.7, 4.4, 4.1, 4, 4.4, 4.6, 4, 3.3, 4.2, 4.2, 4.2, 4.3, 3, 4.1, 6, 5.1, 5.9, 5.6, 5.8, 6.6, 4.5, 6.3, 5.8, 6.1, 5.1, 5.3, 5.5, 5, 5.1, 5.3, 5.5, 6.7, 6.9, 5, 5.7, 4.9, 6.7, 4.9, 5.7, 6, 4.8, 4.9, 5.6, 5.8, 6.1, 6.4, 5.6, 5.1, 5.6, 6.1, 5.6, 5.5, 4.8, 5.4, 5.6, 5.1, 5.1, 5.9, 5.7, 5.2, 5, 5.2, 5.4, 5.1 };
            var petalWidth = new double[] { 0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.3, 0.2, 0.2, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.4, 0.2, 0.5, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.4, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.3, 0.3, 0.2, 0.6, 0.4, 0.3, 0.2, 0.2, 0.2, 0.2, 1.4, 1.5, 1.5, 1.3, 1.5, 1.3, 1.6, 1, 1.3, 1.4, 1, 1.5, 1, 1.4, 1.3, 1.4, 1.5, 1, 1.5, 1.1, 1.8, 1.3, 1.5, 1.2, 1.3, 1.4, 1.4, 1.7, 1.5, 1, 1.1, 1, 1.2, 1.6, 1.5, 1.6, 1.5, 1.3, 1.3, 1.3, 1.2, 1.4, 1.2, 1, 1.3, 1.2, 1.3, 1.3, 1.1, 1.3, 2.5, 1.9, 2.1, 1.8, 2.2, 2.1, 1.7, 1.8, 1.8, 2.5, 2, 1.9, 2.1, 2, 2.4, 2.3, 1.8, 2.2, 2.3, 1.5, 2.3, 2, 2, 1.8, 2.1, 1.8, 1.8, 1.8, 2.1, 1.6, 1.9, 2, 2.2, 1.5, 1.4, 2.3, 2.4, 1.8, 1.8, 2.1, 2.4, 2.3, 1.9, 2.3, 2.5, 2.3, 1.9, 2, 2.3, 1.8 };
            var featureList = new List<double[]> { sepalLength, sepalWidth, petalLength, petalWidth };
            var features = new Matrix(featureList) { Header = new string[] { "Sepal Length", "Sepal Width", "Petal Length", "Petal Width" } };

            var gmm = new GaussianMixtureModel(features, 3);
            gmm.Train(12345);

            // Test cluster counts
            Assert.AreEqual(0.3005423, gmm.Weights[0], 1E-2);
            Assert.AreEqual(0.3661243, gmm.Weights[1], 1E-2);
            Assert.AreEqual(0.3333333, gmm.Weights[2], 1E-2);

            // Test cluster means
            var trueMean1 = new double[] { 5.915044, 2.777451, 4.204002, 1.298935 };
            var trueMean2 = new double[] { 6.546807, 2.949613, 5.482252, 1.985523 };
            var trueMean3 = new double[] { 5.006000, 3.428000, 1.462000, 0.246000 };
            for (int i = 0; i < 4; i++)
            {
                Assert.AreEqual(trueMean1[i], gmm.Means[0, i], 1E-2);
                Assert.AreEqual(trueMean2[i], gmm.Means[1, i], 1E-2);
                Assert.AreEqual(trueMean3[i], gmm.Means[2, i], 1E-2);
            }

            // The log-likelihood carries the full Gaussian normalizing constant. R mclust
            // (Mclust, G = 3, VVV) reports -180.185 and Python sklearn 1.8.0
            // (GaussianMixture(3, covariance_type='full'), score(X) * n) reports -180.18547759250404
            // for this fixture; the window covers EM stopping differences between implementations.
            Assert.AreEqual(-180.185, gmm.LogLikelihood, 0.5);
        }

        /// <summary>
        /// A training run that exhausts its iteration budget reports the log-likelihood of its
        /// final iteration rather than a placeholder.
        /// </summary>
        [TestMethod]
        public void Test_GMM_LogLikelihood_AtIterationCap()
        {
            var x = new double[100];
            var y = new double[100];
            var rnd = new MersenneTwister(4242);
            for (int i = 0; i < 100; i++)
            {
                double t = i < 50 ? 0.0 : 8.0;
                x[i] = t + rnd.NextDouble();
                y[i] = t + rnd.NextDouble();
            }
            var gmm = new GaussianMixtureModel(new Matrix(new List<double[]> { x, y }), 2) { MaxIterations = 2 };
            gmm.Train(12345);

            Assert.IsLessThan(0, gmm.LogLikelihood, "the iteration-capped run must report its final evaluated log-likelihood");
            Assert.IsFalse(double.IsNaN(gmm.LogLikelihood) || double.IsInfinity(gmm.LogLikelihood));
        }
        /// <summary>
        /// A mixture initialized from a degenerate clustering keeps finite parameters: a component
        /// that receives no responsibility retains its previous state instead of spreading NaN.
        /// </summary>
        [TestMethod]
        public void Test_GMM_DegenerateFixture_StaysFinite()
        {
            var x = new double[12];
            var y = new double[12];
            for (int i = 6; i < 12; i++) { x[i] = 10; y[i] = 10; }
            var gmm = new GaussianMixtureModel(new Matrix(new List<double[]> { x, y }), 3);
            gmm.Train(12345);

            Assert.IsFalse(double.IsNaN(gmm.LogLikelihood));
            for (int k = 0; k < 3; k++)
                for (int d = 0; d < 2; d++)
                    Assert.IsFalse(double.IsNaN(gmm.Means[k, d]), $"mean [{k},{d}] is NaN");
        }

        /// <summary>
        /// A singular initial component covariance remains a loud training failure, but the
        /// public workflow identifies the component and preserves the factorization context.
        /// </summary>
        [TestMethod]
        public void Test_GMM_SingularInitialComponentReportsContextualFailure()
        {
            var data = new double[,]
            {
                { 0d, 0d }, { 1d, 1d }, { 2d, 2d },
                { 10d, 10d }, { 11d, 11d }, { 12d, 12d }
            };
            var gmm = new GaussianMixtureModel(data, 2);

            var exception = Assert.Throws<InvalidOperationException>(() => gmm.Train(12345));

            StringAssert.Contains(exception.Message, "component");
            Assert.IsNotNull(exception.InnerException);
        }

        /// <summary>
        /// Verifies that the M-step does not ridge covariance matrices that are already positive definite.
        /// </summary>
        /// <remarks>
        /// The helper returns the symmetrized input unchanged when Cholesky accepts it. This test
        /// recomputes each M-step covariance externally from the public responsibilities after one EM
        /// iteration and requires exact agreement after applying only the established diagonal floor.
        /// </remarks>
        [TestMethod]
        public void Test_GMM_MStep_WellConditionedCovarianceIsNotRidged()
        {
            var data = new double[,]
            {
                { 1.0, 2.1 }, { 1.2, 1.9 }, { 0.8, 2.3 }, { 1.1, 2.0 }, { 0.9, 1.8 }, { 1.3, 2.2 },
                { 8.0, 9.1 }, { 8.2, 8.9 }, { 7.8, 9.3 }, { 8.1, 9.0 }, { 7.9, 8.8 }, { 8.3, 9.2 }
            };
            var gmm = new GaussianMixtureModel(data, 2) { MaxIterations = 1 };
            gmm.Train(seed: 42);

            int n = data.GetLength(0);
            int dims = data.GetLength(1);

            // Per-dimension population variance of the whole sample, matching the M-step's floor.
            var colVar = new double[dims];
            for (int d = 0; d < dims; d++)
            {
                double colMean = 0;
                for (int i = 0; i < n; i++)
                    colMean += data[i, d];
                colMean /= n;
                double v = 0;
                for (int i = 0; i < n; i++)
                    v += (data[i, d] - colMean) * (data[i, d] - colMean);
                colVar[d] = v / n;
            }

            for (int k = 0; k < 2; k++)
            {
                double wgt = 0d;
                for (int i = 0; i < n; i++)
                    wgt += gmm.LikelihoodMatrix[i, k];
                Assert.IsGreaterThan(0d, wgt, "Fixture precondition: both components carry responsibility.");

                // Recompute the raw M-step covariance from the responsibilities and stored means.
                var expected = new double[dims, dims];
                for (int d = 0; d < dims; d++)
                {
                    for (int j = 0; j < dims; j++)
                    {
                        double sum = 0;
                        for (int i = 0; i < n; i++)
                            sum += gmm.LikelihoodMatrix[i, k] * (data[i, d] - gmm.Means[k, d]) * (data[i, j] - gmm.Means[k, j]);
                        expected[d, j] = sum / wgt;
                    }
                }
                // Apply the established diagonal floor. No ridge is needed for these covariances.
                for (int d = 0; d < dims; d++)
                    expected[d, d] = System.Math.Max(expected[d, d], 1E-6 * colVar[d]);

                for (int d = 0; d < dims; d++)
                    for (int j = 0; j < dims; j++)
                        Assert.AreEqual(expected[d, j], gmm.Sigmas[k][d, j], 0d,
                            $"Sigma[{k}][{d},{j}] must preserve the un-ridged covariance.");
            }
        }

        /// <summary>
        /// Verifies that a ridge required by a rank-deficient M-step covariance is stored for the next E-step.
        /// </summary>
        /// <remarks>
        /// The full fixture has a usable initial covariance. The test then assigns responsibility only to
        /// the first three collinear observations and invokes one M-step, producing the exactly rank-one
        /// covariance <c>[[2/3, 2/3], [2/3, 2/3]]</c>. The returned trace-scaled ridge must be assigned to
        /// <see cref="GaussianMixtureModel.Sigmas"/>; discarding the pure helper's return value leaves the
        /// stored covariance singular.
        /// </remarks>
        [TestMethod]
        public void Test_GMM_MStep_RequiredPositiveDefiniteRepairIsStored()
        {
            var data = new double[,]
            {
                { 0d, 0d }, { 1d, 1d }, { 2d, 2d },
                { 0d, 2d }, { 1d, 0d }, { 2d, 0d }
            };
            var gmm = new GaussianMixtureModel(data, 1) { MaxIterations = 1 };
            gmm.Train(seed: 42);

            for (int i = 0; i < data.GetLength(0); i++)
                gmm.LikelihoodMatrix[i, 0] = i < 3 ? 1d : 0d;

            MethodInfo mStep = typeof(GaussianMixtureModel).GetMethod(
                "MStep",
                BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.IsNotNull(mStep);
            mStep.Invoke(gmm, null);

            double rawVariance = 2d / 3d;
            double baseRidge = 1E-10d * rawVariance;
            Assert.AreEqual(rawVariance + baseRidge, gmm.Sigmas[0][0, 0], 0d);
            Assert.AreEqual(rawVariance + baseRidge, gmm.Sigmas[0][1, 1], 0d);
            Assert.AreEqual(rawVariance, gmm.Sigmas[0][0, 1], 0d);
            Assert.AreEqual(rawVariance, gmm.Sigmas[0][1, 0], 0d);
            Assert.IsTrue(new CholeskyDecomposition(gmm.Sigmas[0]).IsPositiveDefinite);
        }

        /// <summary>
        /// Verifies one M-step against scikit-learn's GaussianMixture with full covariances.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Oracle: scikit-learn 1.9.1, <c>GaussianMixture(2, covariance_type='full', reg_covar=0.0, max_iter=1,
        /// weights_init=[0.5, 0.5], means_init=[[0.2, 0.3], [2.0, 2.0]], precisions_init=inv([[[1, 0.3], [0.3, 1]],
        /// [[1, -0.2], [-0.2, 1.5]]])).fit(X)</c>, which runs one E-step at the initial parameters and one M-step.
        /// The responsibilities below are that E-step's, from <c>_estimate_log_prob_resp</c> at the initial
        /// parameters, and the expected values are the fitted <c>weights_</c>, <c>means_</c> and
        /// <c>covariances_</c>.
        /// </para>
        /// <para>
        /// The test sets these responsibilities and runs one M-step. The covariance floor does not bind on this
        /// fixture, so the result is the plain EM update. scikit-learn adds 10 machine epsilon to each component's
        /// responsibility sum, which moves the values by about 1E-15.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_GMM_MStep_MatchesScikitLearn()
        {
            var data = new double[,]
            {
                { 0.0, 0.2 }, { 0.5, -0.3 }, { 1.1, 0.9 }, { -0.4, 0.1 }, { 0.8, 1.4 },
                { 1.6, 1.2 }, { 2.3, 2.0 }, { 1.9, 2.6 }, { 2.8, 1.7 }, { 1.2, 2.1 }
            };
            var responsibilities = new double[]
            {
                0.9796283429489151, 0.020371657051085058, 0.9667488260898032, 0.03325117391019694,
                0.6764783226423449, 0.3235216773576551, 0.9924304009977198, 0.007569599002279958,
                0.6348226281157242, 0.3651773718842756, 0.37168045591962245, 0.6283195440803774,
                0.07232849603643313, 0.9276715039635671, 0.054737405112947946, 0.9452625948870519,
                0.04617142321376257, 0.9538285767862374, 0.23529370279634693, 0.7647062972036534
            };
            var weights = new[] { 0.5030320003873621, 0.4969679996126379 };
            var means = new[] { 0.519860421693883, 0.5582796355392851, 1.8481945977789889, 1.8294286128438793 };
            var covariances = new[]
            {
                0.4966605576296129, 0.37085605948641437, 0.3708560594864144, 0.5419387205323211,
                0.435968002661546, 0.10691915725472677, 0.10691915725472677, 0.2985049441544398
            };

            var gmm = new GaussianMixtureModel(data, 2) { MaxIterations = 1 };
            gmm.Train(12345);
            for (int i = 0; i < data.GetLength(0); i++)
                for (int k = 0; k < 2; k++)
                    gmm.LikelihoodMatrix[i, k] = responsibilities[2 * i + k];

            MethodInfo mStep = typeof(GaussianMixtureModel).GetMethod(
                "MStep",
                BindingFlags.Instance | BindingFlags.NonPublic);
            Assert.IsNotNull(mStep);
            mStep.Invoke(gmm, null);

            const double tolerance = 1E-12;
            for (int k = 0; k < 2; k++)
            {
                Assert.AreEqual(weights[k], gmm.Weights[k], tolerance);
                for (int d = 0; d < 2; d++)
                {
                    Assert.AreEqual(means[2 * k + d], gmm.Means[k, d], tolerance);
                    for (int j = 0; j < 2; j++)
                        Assert.AreEqual(covariances[4 * k + 2 * d + j], gmm.Sigmas[k][d, j], tolerance);
                }
            }
        }

        /// <summary>
        /// Test that expectation-maximization training reproduces fixed reference means, weights
        /// and covariance exactly for a fixed seed.
        /// </summary>
        /// <remarks>
        /// Floating-point math functions on .NET Framework can differ from .NET in the last bit, and
        /// the iterative fit carries such differences into the last bits of the estimates, so that
        /// runtime has its own reference values.
        /// </remarks>
        [TestMethod]
        public void Test_GMM_ReferenceValuesExact()
        {
            var gmm = new GaussianMixtureModel(MachineLearningFixtures.ReferenceFixture(), 3);
            gmm.Train(12345);
#if NETFRAMEWORK
            MachineLearningFixtures.AssertExact(new[] { 1.6836975512647632, 0.15252813217254019, 6.5311314435756849, 4.9214154281214215, -1.5193501223541142, 0.50401330706019221 }, gmm.Means);
            CollectionAssert.AreEqual(new[] { 0.19183326669715589, 0.60927329841690037, 0.19889343488594363 }, gmm.Weights);
            MachineLearningFixtures.AssertExact(new[] { 0.896469916487836, -1.0502016999406476, -1.0502016999406476, 1.5091520124954216 }, gmm.Sigmas[0].ToArray());
#else
            MachineLearningFixtures.AssertExact(new[] { 1.6836975512647632, 0.1525281321725402, 6.531131443575685, 4.9214154281214215, -1.5193501223541142, 0.5040133070601922 }, gmm.Means);
            CollectionAssert.AreEqual(new[] { 0.1918332666971559, 0.6092732984169004, 0.19889343488594363 }, gmm.Weights);
            MachineLearningFixtures.AssertExact(new[] { 0.896469916487836, -1.0502016999406476, -1.0502016999406476, 1.5091520124954216 }, gmm.Sigmas[0].ToArray());
#endif
        }

        /// <summary>
        /// Test that a fit large enough for the parallel expectation and maximization passes reproduces
        /// fixed reference means, log-likelihood and iteration count exactly, independent of the
        /// processor count.
        /// </summary>
        [TestMethod]
        public void Test_GMM_LargeReferenceValuesExact()
        {
            var gmm = new GaussianMixtureModel(MachineLearningFixtures.LargeFixture(6000), 5) { MaxIterations = 50 };
            gmm.Train(11);
            MachineLearningFixtures.AssertExact(new[] { 9.5066066066066774, 14.000600600600507, 18.5021021021021, 2.9885885885885908, 13.494834834834855, 18.011351351351379, 2.5003403403403386, 7.0118518518518052, 3.4998987182798067, 7.9987795552719803, 12.517916736297826, 17.016797573289967, 17.499779779779786, 1.991271271271271, 6.4977777777777463, 10.991771771771708, 3.6312741312741319, 7.9257829257829249, 10.933290433290434, 15.656799656799659 }, gmm.Means);
            Assert.AreEqual(19299.806045657064, gmm.LogLikelihood, 0d);
            Assert.AreEqual(7, gmm.Iterations);
        }

    }
}
