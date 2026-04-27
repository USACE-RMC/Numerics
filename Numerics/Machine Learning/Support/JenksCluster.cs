using System;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Supporting class for a Jenks natural breaks cluster.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Jenks optimization method, also called the Jenks natural breaks classification method, 
    /// is a data clustering method designed to determine the best arrangement of values into different classes.
    /// </para>
    /// </remarks>
    public class JenksCluster
    {
        /// <summary>
        /// Creates a new Jenks cluster. 
        /// </summary>
        /// <param name="data">The sorted input data array.</param>
        /// <param name="startIndex">The starting index of the cluster.</param>
        /// <param name="endIndex">The ending index of the cluster.</param>
        public JenksCluster(double[] data, int startIndex,  int endIndex)
        {
            StartIndex = startIndex;
            EndIndex = endIndex;
            MinValue = data[startIndex];
            MaxValue = data[endIndex];

            // Compute summary statistics using Welford's algorithm for numerical stability
            double sum = 0;
            double mean = 0;
            double m2 = 0;
            for (int i = startIndex; i <= endIndex; i++)
            {
                sum += data[i];
                int k = i - startIndex + 1;
                double delta = data[i] - mean;
                mean += delta / k;
                double delta2 = data[i] - mean;
                m2 += delta * delta2;
            }

            Sum = sum;
            Average = mean;
            SumOfSquaredDeviations = m2;
            Variance = Count <= 1 ? 0 : m2 / (Count - 1);
        }

        /// <summary>
        /// The starting index of the cluster.
        /// </summary>
        public int StartIndex { get; private set; }

        /// <summary>
        /// The ending index of the cluster.
        /// </summary>
        public int EndIndex { get; private set; }

        /// <summary>
        /// The number of data points in the cluster.
        /// </summary>
        public int Count { get { return EndIndex - StartIndex + 1; } }

        /// <summary>
        /// The minimum value of the cluster.
        /// </summary>
        public double MinValue { get; private set; }

        /// <summary>
        /// The maximum value of the cluster.
        /// </summary>
        public double MaxValue { get; private set; }

        /// <summary>
        /// The sum of the values in the cluster.
        /// </summary>
        public double Sum { get; private set; }

        /// <summary>
        /// The average value of the cluster. 
        /// </summary>
        public double Average { get; private set; }

        /// <summary>
        /// The variance of the cluster.
        /// </summary>
        public double Variance { get; private set; }

        /// <summary>
        /// The sum of squared deviations.
        /// </summary>
        public double SumOfSquaredDeviations { get; private set; }

    }
}
