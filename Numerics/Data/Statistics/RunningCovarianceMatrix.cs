using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// A class for keeping track of a running, regularized covariance matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The co-moment matrix starts at the identity, not at zero, and each pushed vector adds its Welford update,
    /// so after N vectors <see cref="Covariance"/> is M₂ + I, where M₂ is the sum of outer products of the deviations
    /// from the mean. <see cref="SampleCovariance"/> is therefore (M₂ + I) / (N − 1) and
    /// <see cref="PopulationCovariance"/> is (M₂ + I) / N: the sample statistics shrunk toward the identity by a
    /// prior weight that fades as 1 / N. The identity keeps the matrix positive definite from the first vectors,
    /// which the adaptive random-walk Metropolis-Hastings sampler relies on for its proposal covariance. For the
    /// plain sample covariance, subtract the identity from <see cref="Covariance"/> before scaling.
    /// </para>
    /// </remarks>
    [Serializable]
    public class RunningCovarianceMatrix
    {
        /// <summary>
        /// Construct a covariance matrix with a zero mean and an identity co-moment matrix.
        /// </summary>
        /// <param name="size">The number of rows and columns.</param>
        public RunningCovarianceMatrix(int size)
        {
            // The mean vector set at zero.
            Mean = new Matrix(size, 1);
            // The initial covariance  is the identity matrix
            Covariance = Matrix.Identity(size);
        }

        /// <summary>
        /// The sample size N.
        /// </summary>
        public int N { get; private set; }

        /// <summary>
        /// The mean vector
        /// </summary>
        public Matrix Mean { get; private set; }

        /// <summary>
        /// The co-moment matrix M₂ + I: the sum of outer products of the deviations from the running mean, plus the
        /// identity it starts from. This is unadjusted by the sample size.
        /// </summary>
        public Matrix Covariance { get; private set; }

        /// <summary>
        /// The regularized sample covariance matrix (M₂ + I) / (N − 1), which includes the identity prior.
        /// </summary>
        public Matrix SampleCovariance => (1d / (N - 1)) * Covariance;

        /// <summary>
        /// The regularized population covariance matrix (M₂ + I) / N, which includes the identity prior.
        /// </summary>
        public Matrix PopulationCovariance => (1d / N) * Covariance;

        /// <summary>
        /// The sample correlation matrix.
        /// </summary>
        public Matrix SampleCorrelation
        { 
            get
            {
                var D = Matrix.Diagonal(SampleCovariance);
                D.Sqrt();
                var invSqrtD = !D;
                return (invSqrtD * SampleCovariance) * invSqrtD;
            }       
        }

        /// <summary>
        /// The population correlation matrix.
        /// </summary>
        public Matrix PopulationCorrelation
        {
            get
            {
                var D = Matrix.Diagonal(PopulationCovariance);
                D.Sqrt();
                var invSqrtD = !D;
                return (invSqrtD * PopulationCovariance) * invSqrtD;
            }
        }

        /// <summary>
        /// Add a new vector to the running statistics. 
        /// </summary>
        /// <param name="values">Vector of data values. The length of the vector but be the same as 
        /// the number of rows</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="values"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the length differs from the matrix size.</exception>
        /// <remarks>
        /// The Welford update is written element by element: the new mean is m + (x − m)/N (or m + x for the first
        /// vector) and each co-moment entry becomes C[i, j] + (0 + (xᵢ − mᵢ)(xⱼ − m′ⱼ)), with m the previous and m′
        /// the new mean. These are the operations, in the same order, of the matrix form
        /// M + (x − M)/N and C + (x − M)(x − M′)ᵀ. New <see cref="Mean"/> and <see cref="Covariance"/> matrices are created on each
        /// call, so references obtained earlier keep their values.
        /// </remarks>
        public void Push(IList<double> values)
        {
            N += 1;
            if (values == null) throw new ArgumentNullException(nameof(values));
            int size = Mean.NumberOfRows;
            if (values.Count != size) throw new ArgumentException("The matrix must have the same number of rows.", nameof(values));

            var oldMean = Mean.Array;
            var newMatrix = new Matrix(size, 1);
            var mean = newMatrix.Array;
            double scale = 1d / N;
            var x = new double[size];
            for (int i = 0; i < size; i++)
            {
                x[i] = values[i];
                double m = oldMean[i, 0];
                mean[i, 0] = m + (N == 1 ? x[i] : (x[i] - m) * scale);
            }

            var oldCovariance = Covariance.Array;
            var newCovariance = new Matrix(size, size);
            var c = newCovariance.Array;
            var after = new double[size];
            for (int j = 0; j < size; j++) after[j] = x[j] - mean[j, 0];
            for (int i = 0; i < size; i++)
            {
                double before = x[i] - oldMean[i, 0];
                for (int j = 0; j < size; j++)
                    c[i, j] = oldCovariance[i, j] + (0.0d + before * after[j]);
            }
            Mean = newMatrix;
            Covariance = newCovariance;
        }

    }
}
