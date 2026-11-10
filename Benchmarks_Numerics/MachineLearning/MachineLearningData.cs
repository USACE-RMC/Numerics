using System.Runtime.InteropServices;

namespace Benchmarks_Numerics.MachineLearning
{
    /// <summary>
    /// Loads the benchmark data written by docs/machine-learning/benchmarks/pca_elastic_net_generate_data.py.
    /// </summary>
    public static class MachineLearningData
    {
        /// <summary>
        /// The data folder, from the NUMERICS_ML_BENCHMARKS environment variable.
        /// </summary>
        public static string Folder =>
            Environment.GetEnvironmentVariable("NUMERICS_ML_BENCHMARKS")
            ?? throw new InvalidOperationException("Set NUMERICS_ML_BENCHMARKS to the folder written by pca_elastic_net_generate_data.py.");

        /// <summary>
        /// Reads an n × p row-major little-endian float64 matrix.
        /// </summary>
        /// <param name="name">The file name.</param>
        /// <param name="n">The number of rows.</param>
        /// <param name="p">The number of columns.</param>
        public static double[,] ReadMatrix(string name, int n, int p)
        {
            var bytes = File.ReadAllBytes(Path.Combine(Folder, name));
            var result = new double[n, p];
            MemoryMarshal.Cast<byte, double>(bytes).CopyTo(MemoryMarshal.CreateSpan(ref result[0, 0], n * p));
            return result;
        }

        /// <summary>
        /// Reads a little-endian float64 vector.
        /// </summary>
        /// <param name="name">The file name.</param>
        public static double[] ReadVector(string name)
        {
            return MemoryMarshal.Cast<byte, double>(File.ReadAllBytes(Path.Combine(Folder, name))).ToArray();
        }
    }
}
