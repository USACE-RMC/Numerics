using BenchmarkDotNet.Attributes;
using Numerics.Data;
using Numerics.Data.Statistics;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The autocorrelation about zero to lag 20: the lag-sum form over a contiguous array against the list form it
    /// replaced, which divided the autocovariances S_k / n by S_0 / n.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class UncenteredAutocorrelationBenchmarks
    {
        private const int LagMax = 20;
        private double[] _x = null!;
        private List<double> _list = null!;
        private TimeSeries _series = null!;

        /// <summary>The series length.</summary>
        [Params(30, 100, 1000, 10000, 100000)]
        public int N { get; set; }

        /// <summary>
        /// Prepares a zero-mean AR(1) series.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var rng = new Random(1);
            _x = new double[N];
            double previous = 0d;
            for (int i = 0; i < N; i++)
            {
                previous = 0.6 * previous + rng.NextDouble() - 0.5;
                _x[i] = previous;
            }
            _list = new List<double>(_x);
            _series = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1), _x);
        }

        /// <summary>The replaced list form, written out.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark(Baseline = true)]
        public double[,]? ListForm()
        {
            return Replaced(_x);
        }

        /// <summary>The replaced list form on a <see cref="List{T}"/>.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark]
        public double[,]? ListFormOnList()
        {
            return Replaced(_list);
        }

        /// <summary>The library member on an array.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark]
        public double[,]? Array()
        {
            return Autocorrelation.Function(_x, LagMax, Autocorrelation.Type.Correlation, false);
        }

        /// <summary>The library member on a <see cref="List{T}"/>.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark]
        public double[,]? List()
        {
            return Autocorrelation.Function(_list, LagMax, Autocorrelation.Type.Correlation, false);
        }

        /// <summary>The library member on a time-series.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark]
        public double[,]? Series()
        {
            return Autocorrelation.Function(_series, LagMax, Autocorrelation.Type.Correlation, false);
        }

        /// <summary>The replaced time-series form, written out.</summary>
        /// <returns>The autocorrelation table.</returns>
        [Benchmark]
        public double[,]? SeriesReplaced()
        {
            int n = _series.Count;
            var acf = new double[LagMax + 1, 2];
            for (int lag = 0; lag <= LagMax; lag++)
            {
                double cov = 0d;
                for (int t = lag; t < n; t++)
                    cov += (_series[t].Value - 0d) * (_series[t - lag].Value - 0d);
                acf[lag, 0] = lag;
                acf[lag, 1] = cov / n;
            }
            double den = acf[0, 1];
            for (int i = 0; i <= LagMax; i++) acf[i, 1] /= den;
            return acf;
        }

        private static double[,] Replaced(IList<double> data)
        {
            int n = data.Count;
            var acf = new double[LagMax + 1, 2];
            for (int lag = 0; lag <= LagMax; lag++)
            {
                double cov = 0d;
                for (int t = lag; t < n; t++)
                    cov += (data[t] - 0d) * (data[t - lag] - 0d);
                acf[lag, 0] = lag;
                acf[lag, 1] = cov / n;
            }
            double den = acf[0, 1];
            for (int i = 0; i <= LagMax; i++) acf[i, 1] /= den;
            return acf;
        }
    }
}
