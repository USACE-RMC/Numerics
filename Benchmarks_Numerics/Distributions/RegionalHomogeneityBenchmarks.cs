using BenchmarkDotNet.Attributes;
using Numerics.Distributions.Regional;

namespace Benchmarks_Numerics.Distributions
{
    /// <summary>
    /// The full regional test (discordancy, kappa fit, 500 simulated regions, H1 to H3 and Z) of the 19-site North
    /// Cascades region of Hosking's LMOMENTS package, serial and parallel.
    /// </summary>
    /// <remarks>
    /// The Python comparison is the Stochastics region/homogeneity.py regional_test on the same record lengths with
    /// correlated = False and 500 simulations.
    /// </remarks>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class RegionalHomogeneityBenchmarks
    {
        private RegionalData _data = null!;

        /// <summary>
        /// Builds the Cascades regional data.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            int[] n = [98, 59, 90, 61, 65, 86, 78, 72, 67, 99, 49, 61, 69, 73, 70, 66, 59, 74, 82];
            double[] mean = [19.685, 62.580, 40.852, 46.045, 45.021, 31.042, 80.143, 41.305, 30.585, 32.932,
                17.560, 69.518, 47.653, 102.501, 52.413, 79.696, 44.643, 58.655, 39.024];
            double[] lcv = [0.1209, 0.0915, 0.1124, 0.1032, 0.0967, 0.1328, 0.1008, 0.1143, 0.1107, 0.1179,
                0.1308, 0.1119, 0.1018, 0.1025, 0.1054, 0.1174, 0.1115, 0.1003, 0.1046];
            double[] t3 = [0.0488, 0.0105, 0.0614, 0.0417, -0.0134, -0.0176, 0.0943, 0.0555, 0.0478, 0.0492,
                0.0940, -0.0429, 0.0435, 0.0182, -0.0224, 0.0124, -0.0346, 0.0446, 0.0128];
            double[] t4 = [0.1433, 0.1569, 0.1541, 0.1429, 0.1568, 0.1206, 0.1967, 0.1210, 0.1371, 0.0900,
                0.1273, 0.0927, 0.1446, 0.1047, 0.1664, 0.1317, 0.1032, 0.1450, 0.1583];
            var names = new string[n.Length];
            for (int i = 0; i < names.Length; i++) names[i] = (i + 1).ToString();
            _data = new RegionalData(names, n, mean, lcv, t3, t4);
        }

        /// <summary>The test on one thread.</summary>
        /// <returns>The heterogeneity measures.</returns>
        [Benchmark(Baseline = true)]
        public IReadOnlyList<double> RegionalTestSerial() => new RegionalHomogeneityTest(_data, 500, 12345, 1).Heterogeneity;

        /// <summary>The test on the default scheduler.</summary>
        /// <returns>The heterogeneity measures.</returns>
        [Benchmark]
        public IReadOnlyList<double> RegionalTestParallel() => new RegionalHomogeneityTest(_data, 500, 12345).Heterogeneity;
    }
}
