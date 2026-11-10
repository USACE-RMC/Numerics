using BenchmarkDotNet.Attributes;
using Numerics.MachineLearning;

namespace Benchmarks_Numerics.MachineLearning
{
    /// <summary>
    /// ElasticNetRegression on the sizes timed by pca_elastic_net_time_python.py (scikit-learn enet_path and ElasticNetCV) and
    /// pca_elastic_net_time_r.R (glmnet and cv.glmnet): α = 0.5, 100 penalties, default tolerance 1E-7, 10 folds.
    /// </summary>
    [MemoryDiagnoser]
    public class ElasticNetBenchmarks
    {
        private double[,] _x = null!;
        private double[] _y = null!;

        /// <summary>
        /// The data size as "n x p".
        /// </summary>
        [Params("100x200", "1000x50", "10000x100", "100000x500")]
        public string Size { get; set; } = "";

        /// <summary>
        /// Loads the data for the current size.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var parts = Size.Split('x');
            int n = int.Parse(parts[0]), p = int.Parse(parts[1]);
            _x = MachineLearningData.ReadMatrix($"enet_{n}x{p}_X.bin", n, p);
            _y = MachineLearningData.ReadVector($"enet_{n}x{p}_y.bin");
        }

        /// <summary>
        /// Fits the coefficient path.
        /// </summary>
        [Benchmark]
        public int Path()
        {
            var model = new ElasticNetRegression(_x, _y) { Alpha = 0.5 };
            model.Fit();
            return model.Lambdas.Length;
        }

        /// <summary>
        /// Fits the path and runs 10-fold cross-validation.
        /// </summary>
        [Benchmark]
        public int CrossValidated()
        {
            var model = new ElasticNetRegression(_x, _y) { Alpha = 0.5 };
            model.CrossValidate(10, 1);
            return model.LambdaMinIndex;
        }
    }
}
