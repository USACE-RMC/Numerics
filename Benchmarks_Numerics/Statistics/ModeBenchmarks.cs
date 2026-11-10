using BenchmarkDotNet.Attributes;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using System.Reflection;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The global Mode of KernelDensity, TruncatedDistribution and EmpiricalDistribution against the single bounded Brent
    /// search of the density over [Q(0.001), Q(0.999)] that preceded it.
    /// </summary>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class ModeBenchmarks
    {
        private KernelDensity _kde = null!;
        private TruncatedDistribution _truncatedNormal = null!, _truncatedMixture = null!;
        private EmpiricalDistribution _empirical = null!;

        /// <summary>The kernel density sample size.</summary>
        [Params(200, 2000, 100000)]
        public int N { get; set; }

        /// <summary>
        /// Builds a bimodal sample of N/2 quantiles of N(-2, 0.6) and N/2 of N(2, 1.2), and the truncated and empirical
        /// distributions.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            var left = new Normal(-2d, 0.6);
            var right = new Normal(2d, 1.2);
            var sample = new List<double>();
            int half = N / 2;
            for (int i = 0; i < half; i++) sample.Add(left.InverseCDF((i + 0.5) / half));
            for (int i = 0; i < half; i++) sample.Add(right.InverseCDF((i + 0.5) / half));
            _kde = new KernelDensity(sample);
            _truncatedNormal = new TruncatedDistribution(new Normal(10d, 3d), 8d, 25d);
            var mixture = new Mixture(new[] { 0.5, 0.5 }, new IUnivariateDistribution[] { left, right });
            _truncatedMixture = new TruncatedDistribution(mixture, -3.8, 5.6);
            _empirical = new EmpiricalDistribution(sample);
        }

        private static double LocalMode(UnivariateDistributionBase distribution)
        {
            var brent = new BrentSearch(distribution.PDF, distribution.InverseCDF(0.001), distribution.InverseCDF(0.999));
            brent.Maximize();
            return brent.BestParameterSet.Values[0];
        }

        /// <summary>The previous local search on the kernel density.</summary>
        /// <returns>The mode.</returns>
        [Benchmark(Baseline = true)]
        public double KdeLocal() => LocalMode(_kde);

        /// <summary>The global kernel density mode.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double KdeGlobal() => _kde.Mode;

        /// <summary>The previous local search on a truncated normal.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double TruncatedNormalLocal() => LocalMode(_truncatedNormal);

        /// <summary>The global truncated normal mode.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double TruncatedNormalGlobal() => _truncatedNormal.Mode;

        /// <summary>The previous local search on a truncated bimodal mixture.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double TruncatedMixtureLocal() => LocalMode(_truncatedMixture);

        /// <summary>The global truncated mixture mode.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double TruncatedMixtureGlobal() => _truncatedMixture.Mode;

        /// <summary>The previous local search on the empirical distribution.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double EmpiricalLocal() => LocalMode(_empirical);

        /// <summary>The highest-bin empirical mode, served from its cache after the first call.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double EmpiricalHighestBin() => _empirical.Mode;

        private static readonly FieldInfo ModeCacheField =
            typeof(EmpiricalDistribution).GetField("_modeCache", BindingFlags.Instance | BindingFlags.NonPublic)!;

        /// <summary>The highest-bin empirical mode on a first call: the cached mode is cleared before each scan.</summary>
        /// <returns>The mode.</returns>
        [Benchmark]
        public double EmpiricalHighestBinFirstCall()
        {
            ModeCacheField.SetValue(_empirical, null);
            return _empirical.Mode;
        }
    }
}
