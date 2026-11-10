using BenchmarkDotNet.Attributes;
using Numerics.MachineLearning;

namespace Benchmarks_Numerics.MachineLearning
{
    /// <summary>
    /// PrincipalComponentAnalysis.Train on the sizes timed by pca_elastic_net_time_python.py (scikit-learn PCA().fit) and
    /// pca_elastic_net_time_r.R (prcomp(X, retx = FALSE)).
    /// </summary>
    [MemoryDiagnoser]
    public class PcaBenchmarks
    {
        private double[,] _x = null!;

        /// <summary>
        /// The data size as "n x p".
        /// </summary>
        [Params("1000x20", "10000x50", "100000x100")]
        public string Size { get; set; } = "";

        /// <summary>
        /// Loads the data for the current size.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var parts = Size.Split('x');
            int n = int.Parse(parts[0]), p = int.Parse(parts[1]);
            _x = MachineLearningData.ReadMatrix($"pca_{n}x{p}.bin", n, p);
        }

        /// <summary>
        /// Fits all components on the covariance basis.
        /// </summary>
        [Benchmark]
        public double Fit()
        {
            var pca = new PrincipalComponentAnalysis(_x);
            pca.Train();
            return pca.ExplainedVariance[0];
        }
    }
}
