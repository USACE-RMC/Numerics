using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// A class for keeping track of a running covariance matrix.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class RunningCovarianceMatrix
    {
        /// <summary>
        /// Construct a covariance matrix.
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
        /// The covariance matrix. This is unadjusted by the sample size.
        /// </summary>
        public Matrix Covariance { get; private set; }

        /// <summary>
        /// The sample covariance matrix corrected by the sample size with the degrees of freedom adjustment.
        /// </summary>
        public Matrix SampleCovariance => (1d / (N - 1)) * Covariance;

        /// <summary>
        /// The population covariance matrix corrected by the total sample size.
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
        public void Push(IList<double> values)
        {
            N += 1;
            var x = new Matrix(values.ToArray());
            // Update mean
            var oldMean = Mean;
            Mean += N == 1 ? x : (x - Mean) * (1d / N);
            // Update covariance                 
            Covariance += (x - oldMean) * Matrix.Transpose(x - Mean);
        }

    }
}
