using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.MachineLearning;
using System.Collections.Generic;
using Numerics.Mathematics.LinearAlgebra;
using System.Linq;

namespace MachineLearning
{
    /// <summary>
    /// Unit tests for K-Means classification.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_KMeans
    {
        /// <summary>
        /// K-Means is tested against the Iris dataset in R.
        /// </summary>
        [TestMethod]
        public void Test_KMeans_Iris()
        {
            var sepalLength = new double[] { 5.1, 4.9, 4.7, 4.6, 5, 5.4, 4.6, 5, 4.4, 4.9, 5.4, 4.8, 4.8, 4.3, 5.8, 5.7, 5.4, 5.1, 5.7, 5.1, 5.4, 5.1, 4.6, 5.1, 4.8, 5, 5, 5.2, 5.2, 4.7, 4.8, 5.4, 5.2, 5.5, 4.9, 5, 5.5, 4.9, 4.4, 5.1, 5, 4.5, 4.4, 5, 5.1, 4.8, 5.1, 4.6, 5.3, 5, 7, 6.4, 6.9, 5.5, 6.5, 5.7, 6.3, 4.9, 6.6, 5.2, 5, 5.9, 6, 6.1, 5.6, 6.7, 5.6, 5.8, 6.2, 5.6, 5.9, 6.1, 6.3, 6.1, 6.4, 6.6, 6.8, 6.7, 6, 5.7, 5.5, 5.5, 5.8, 6, 5.4, 6, 6.7, 6.3, 5.6, 5.5, 5.5, 6.1, 5.8, 5, 5.6, 5.7, 5.7, 6.2, 5.1, 5.7, 6.3, 5.8, 7.1, 6.3, 6.5, 7.6, 4.9, 7.3, 6.7, 7.2, 6.5, 6.4, 6.8, 5.7, 5.8, 6.4, 6.5, 7.7, 7.7, 6, 6.9, 5.6, 7.7, 6.3, 6.7, 7.2, 6.2, 6.1, 6.4, 7.2, 7.4, 7.9, 6.4, 6.3, 6.1, 7.7, 6.3, 6.4, 6, 6.9, 6.7, 6.9, 5.8, 6.8, 6.7, 6.7, 6.3, 6.5, 6.2, 5.9 };
            var sepalWidth = new double[] { 3.5, 3, 3.2, 3.1, 3.6, 3.9, 3.4, 3.4, 2.9, 3.1, 3.7, 3.4, 3, 3, 4, 4.4, 3.9, 3.5, 3.8, 3.8, 3.4, 3.7, 3.6, 3.3, 3.4, 3, 3.4, 3.5, 3.4, 3.2, 3.1, 3.4, 4.1, 4.2, 3.1, 3.2, 3.5, 3.6, 3, 3.4, 3.5, 2.3, 3.2, 3.5, 3.8, 3, 3.8, 3.2, 3.7, 3.3, 3.2, 3.2, 3.1, 2.3, 2.8, 2.8, 3.3, 2.4, 2.9, 2.7, 2, 3, 2.2, 2.9, 2.9, 3.1, 3, 2.7, 2.2, 2.5, 3.2, 2.8, 2.5, 2.8, 2.9, 3, 2.8, 3, 2.9, 2.6, 2.4, 2.4, 2.7, 2.7, 3, 3.4, 3.1, 2.3, 3, 2.5, 2.6, 3, 2.6, 2.3, 2.7, 3, 2.9, 2.9, 2.5, 2.8, 3.3, 2.7, 3, 2.9, 3, 3, 2.5, 2.9, 2.5, 3.6, 3.2, 2.7, 3, 2.5, 2.8, 3.2, 3, 3.8, 2.6, 2.2, 3.2, 2.8, 2.8, 2.7, 3.3, 3.2, 2.8, 3, 2.8, 3, 2.8, 3.8, 2.8, 2.8, 2.6, 3, 3.4, 3.1, 3, 3.1, 3.1, 3.1, 2.7, 3.2, 3.3, 3, 2.5, 3, 3.4, 3 };
            var petalLength = new double[] { 1.4, 1.4, 1.3, 1.5, 1.4, 1.7, 1.4, 1.5, 1.4, 1.5, 1.5, 1.6, 1.4, 1.1, 1.2, 1.5, 1.3, 1.4, 1.7, 1.5, 1.7, 1.5, 1, 1.7, 1.9, 1.6, 1.6, 1.5, 1.4, 1.6, 1.6, 1.5, 1.5, 1.4, 1.5, 1.2, 1.3, 1.4, 1.3, 1.5, 1.3, 1.3, 1.3, 1.6, 1.9, 1.4, 1.6, 1.4, 1.5, 1.4, 4.7, 4.5, 4.9, 4, 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4, 4.7, 3.6, 4.4, 4.5, 4.1, 4.5, 3.9, 4.8, 4, 4.9, 4.7, 4.3, 4.4, 4.8, 5, 4.5, 3.5, 3.8, 3.7, 3.9, 5.1, 4.5, 4.5, 4.7, 4.4, 4.1, 4, 4.4, 4.6, 4, 3.3, 4.2, 4.2, 4.2, 4.3, 3, 4.1, 6, 5.1, 5.9, 5.6, 5.8, 6.6, 4.5, 6.3, 5.8, 6.1, 5.1, 5.3, 5.5, 5, 5.1, 5.3, 5.5, 6.7, 6.9, 5, 5.7, 4.9, 6.7, 4.9, 5.7, 6, 4.8, 4.9, 5.6, 5.8, 6.1, 6.4, 5.6, 5.1, 5.6, 6.1, 5.6, 5.5, 4.8, 5.4, 5.6, 5.1, 5.1, 5.9, 5.7, 5.2, 5, 5.2, 5.4, 5.1 };
            var petalWidth = new double[] { 0.2, 0.2, 0.2, 0.2, 0.2, 0.4, 0.3, 0.2, 0.2, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.4, 0.4, 0.3, 0.3, 0.3, 0.2, 0.4, 0.2, 0.5, 0.2, 0.2, 0.4, 0.2, 0.2, 0.2, 0.2, 0.4, 0.1, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.3, 0.3, 0.2, 0.6, 0.4, 0.3, 0.2, 0.2, 0.2, 0.2, 1.4, 1.5, 1.5, 1.3, 1.5, 1.3, 1.6, 1, 1.3, 1.4, 1, 1.5, 1, 1.4, 1.3, 1.4, 1.5, 1, 1.5, 1.1, 1.8, 1.3, 1.5, 1.2, 1.3, 1.4, 1.4, 1.7, 1.5, 1, 1.1, 1, 1.2, 1.6, 1.5, 1.6, 1.5, 1.3, 1.3, 1.3, 1.2, 1.4, 1.2, 1, 1.3, 1.2, 1.3, 1.3, 1.1, 1.3, 2.5, 1.9, 2.1, 1.8, 2.2, 2.1, 1.7, 1.8, 1.8, 2.5, 2, 1.9, 2.1, 2, 2.4, 2.3, 1.8, 2.2, 2.3, 1.5, 2.3, 2, 2, 1.8, 2.1, 1.8, 1.8, 1.8, 2.1, 1.6, 1.9, 2, 2.2, 1.5, 1.4, 2.3, 2.4, 1.8, 1.8, 2.1, 2.4, 2.3, 1.9, 2.3, 2.5, 2.3, 1.9, 2, 2.3, 1.8 };
            var featureList = new List<double[]> { sepalLength, sepalWidth, petalLength, petalWidth };
            var features = new Matrix(featureList) { Header = new string[] { "Sepal Length", "Sepal Width", "Petal Length", "Petal Width" } };

            var kMeans = new KMeans(features, 3);
            kMeans.Train(12345);

            // Test cluster counts
            Assert.AreEqual(62, kMeans.Labels.Where(x => x == 0).Count());
            Assert.AreEqual(38, kMeans.Labels.Where(x => x == 1).Count());
            Assert.AreEqual(50, kMeans.Labels.Where(x => x == 2).Count());

            // Test cluster means
            var trueMean1 = new double[] { 5.901613, 2.748387, 4.393548, 1.433871 };
            var trueMean2 = new double[] { 6.850000, 3.073684, 5.742105, 2.071053 };
            var trueMean3 = new double[] { 5.006000, 3.428000, 1.462000, 0.246000 };
            for (int i = 0; i < 4; i++)
            {
                Assert.AreEqual(trueMean1[i], kMeans.Means[0, i], 1E-6);
                Assert.AreEqual(trueMean2[i], kMeans.Means[1, i], 1E-6);
                Assert.AreEqual(trueMean3[i], kMeans.Means[2, i], 1E-6);                
            }
        }

        /// <summary>
        /// A cluster emptied during training relocates to the point farthest from its assigned
        /// centroid, so a degenerate fixture with fewer distinct points than clusters converges to
        /// finite centroids and zero inertia.
        /// </summary>
        /// <remarks>
        /// Reference behavior verified against Python scikit-learn 1.8.0 KMeans(3) on the same
        /// fixture: finite centers and an inertia of exactly zero, since every point coincides
        /// with a centroid at convergence.
        /// </remarks>
        [TestMethod]
        public void Test_KMeans_EmptyClusterRelocation()
        {
            var rows = new double[12, 2];
            for (int i = 6; i < 12; i++) { rows[i, 0] = 10; rows[i, 1] = 10; }
            var km = new KMeans(rows, 3);
            km.Train(12345);

            double inertia = 0;
            for (int i = 0; i < 12; i++)
            {
                Assert.IsFalse(double.IsNaN(km.Means[km.Labels[i], 0]) || double.IsNaN(km.Means[km.Labels[i], 1]));
                inertia += Math.Pow(rows[i, 0] - km.Means[km.Labels[i], 0], 2) + Math.Pow(rows[i, 1] - km.Means[km.Labels[i], 1], 2);
            }
            Assert.AreEqual(0.0, inertia);
        }

        /// <summary>
        /// The cluster count is validated against the data size.
        /// </summary>
        [TestMethod]
        public void Test_KMeans_ClusterCountValidation()
        {
            var rows = new double[12, 2];
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new KMeans(rows, 0));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new KMeans(rows, 13));
        }

        /// <summary>
        /// A single-cluster fit converges on its first E-step, so the M-step must run before that
        /// convergence exit: the reported mean is the mean of the assigned points, not whichever
        /// data point the k-means++ initializer happened to seed.
        /// </summary>
        [TestMethod]
        public void Test_KMeans_SingleCluster_ReportsTheSampleMean()
        {
            var km = new KMeans(new double[] { 1, 2, 3, 10, 11, 12 }, 1);
            km.Train(7);
            Assert.AreEqual(1, km.Iterations);
            Assert.AreEqual(6.5, km.Means[0, 0], 1E-12);
        }

    }
}
