using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.MachineLearning;
using System.Collections.Generic;
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
        /// Verify that the M-step's symmetric positive-definite repair is actually applied to the
        /// stored covariance matrices.
        /// </summary>
        /// <remarks>
        /// MatrixRegularization.MakeSymmetricPositiveDefinite is pure — it returns a symmetrized copy
        /// with a trace-scaled ridge — so the M-step must store that return value for the repair to
        /// reach the covariance the next E-step's Cholesky factorization consumes. This test recomputes
        /// the M-step covariance externally from the public responsibilities after a single EM iteration
        /// and asserts the stored matrices carry the repair's base ridge (1E-10 of the mean diagonal)
        /// exactly; a covariance stored WITHOUT the ridge would miss these assertions by exactly that
        /// amount.
        /// </remarks>
        [TestMethod]
        public void Test_GMM_MStep_PositiveDefiniteRepairIsApplied()
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
                // Apply the diagonal floor, then the repair's base ridge (the first Cholesky attempt
                // succeeds for these well-separated clusters, so exactly one base ridge is added).
                double trace = 0;
                for (int d = 0; d < dims; d++)
                {
                    expected[d, d] = System.Math.Max(expected[d, d], 1E-6 * colVar[d]);
                    trace += expected[d, d];
                }
                double baseRidge = 1e-10 * trace / dims;
                for (int d = 0; d < dims; d++)
                    expected[d, d] += baseRidge;

                for (int d = 0; d < dims; d++)
                    for (int j = 0; j < dims; j++)
                        Assert.AreEqual(expected[d, j], gmm.Sigmas[k][d, j], 0d,
                            $"Sigma[{k}][{d},{j}] must carry the positive-definite repair.");
            }
        }

    }
}
