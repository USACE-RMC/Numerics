using BenchmarkDotNet.Attributes;
using Numerics.Mathematics.Optimization;

namespace Benchmarks_Numerics.Optimization
{
    /// <summary>
    /// Multi-source Dijkstra on a random sparse directed graph in CSR form, the workload of
    /// scipy.sparse.csgraph.dijkstra(min_only=True), with and without the predecessor and nearest-source outputs.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class DijkstraBenchmarks
    {
        private int[] _rowStart = null!, _toNode = null!, _sources = null!, _predecessors = null!, _nearest = null!;
        private double[] _weight = null!;

        /// <summary>The number of nodes; each has five outgoing edges on average.</summary>
        [Params(10000, 100000, 1000000)]
        public int N { get; set; }

        /// <summary>
        /// Builds a graph with 5 * N edges between uniformly random nodes, weights uniform on [0.001, 1.001), and ten
        /// random sources.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(20261108);
            int m = 5 * N;
            var from = new int[m];
            var to = new int[m];
            var w = new double[m];
            for (int k = 0; k < m; k++)
            {
                from[k] = rng.Next(N);
                to[k] = rng.Next(N);
                w[k] = rng.NextDouble() + 1E-3;
            }
            _rowStart = new int[N + 1];
            for (int k = 0; k < m; k++) _rowStart[from[k] + 1]++;
            for (int i = 0; i < N; i++) _rowStart[i + 1] += _rowStart[i];
            var cursor = (int[])_rowStart.Clone();
            _toNode = new int[m];
            _weight = new double[m];
            for (int k = 0; k < m; k++)
            {
                int slot = cursor[from[k]]++;
                _toNode[slot] = to[k];
                _weight[slot] = w[k];
            }
            _sources = new int[10];
            for (int i = 0; i < _sources.Length; i++) _sources[i] = rng.Next(N);
            _predecessors = new int[N];
            _nearest = new int[N];
        }

        /// <summary>The distances only.</summary>
        /// <returns>The distances.</returns>
        [Benchmark(Baseline = true)]
        public double[] Distances()
        {
            return Dijkstra.ShortestDistances(N, _rowStart, _toNode, _weight, _sources);
        }

        /// <summary>The distances with predecessors and nearest sources, as scipy's return_predecessors=True.</summary>
        /// <returns>The distances.</returns>
        [Benchmark]
        public double[] DistancesPredecessorsAndSources()
        {
            return Dijkstra.ShortestDistances(N, _rowStart, _toNode, _weight, _sources, double.PositiveInfinity, _predecessors, _nearest);
        }
    }
}
