using BenchmarkDotNet.Attributes;
using Numerics.Distributions;

namespace Benchmarks_Numerics.Distributions
{
    /// <summary>
    /// The Pearson type III and log-Pearson type III members whose evaluation changes with the skew: the quantile and its
    /// gradient, the conditional moments, the log tails and density, the frequency factor and its derivative, the
    /// real-space log-Pearson skewness and kurtosis, and the maximum-likelihood covariance.
    /// </summary>
    /// <remarks>
    /// Each benchmark sums the member over fixed inputs at one skew: seven probabilities from 1E-6 to 1 - 1E-6, four
    /// conditioning intervals (two censored, two bounded) or nine standardized points from -4 to 4. The skews cover zero,
    /// the small-skew expansions and the gamma path.
    /// </remarks>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class PearsonTypeIIIBenchmarks
    {
        /// <summary>The skew of the standardized distribution.</summary>
        [Params(0.0, 5E-5, 0.002, 0.02, 0.3, 1.2)]
        public double Skew;

        private static readonly double[] Probabilities = { 1E-6, 0.01, 0.1, 0.5, 0.9, 0.99, 1 - 1E-6 };
        private static readonly (double Lower, double Upper)[] Intervals =
            { (double.NegativeInfinity, -1.5), (-1.0, 0.5), (1.5, double.PositiveInfinity), (2.0, 3.0) };

        private PearsonTypeIII _pearson = null!;
        private LogPearsonTypeIII _logPearson = null!;
        private LogPearsonTypeIII _smallScale = null!;

        /// <summary>
        /// Creates the distributions at the benchmark skew.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            _pearson = new PearsonTypeIII(0d, 1d, Skew);
            _logPearson = new LogPearsonTypeIII(3d, 0.3, Skew);
            _smallScale = new LogPearsonTypeIII(3d, 0.02, Skew);
        }

        /// <summary>The quantile at seven probabilities.</summary>
        /// <returns>The sum of the quantiles.</returns>
        [Benchmark]
        public double InverseCDF()
        {
            double sum = 0d;
            foreach (double p in Probabilities) sum += _pearson.InverseCDF(p);
            return sum;
        }

        /// <summary>The skew component of the quantile gradient at seven probabilities.</summary>
        /// <returns>The sum of the skew derivatives.</returns>
        [Benchmark]
        public double QuantileGradient()
        {
            double sum = 0d;
            foreach (double p in Probabilities) sum += _pearson.QuantileGradient(p)[2];
            return sum;
        }

        /// <summary>The conditional moments on four intervals.</summary>
        /// <returns>The sum of the conditional second moments.</returns>
        [Benchmark]
        public double ConditionalMoments()
        {
            double sum = 0d;
            foreach (var (lower, upper) in Intervals) sum += _pearson.ConditionalMoments(lower, upper)[1];
            return sum;
        }

        /// <summary>The log distribution function and log density at nine points.</summary>
        /// <returns>The sum of the logarithms.</returns>
        [Benchmark]
        public double LogCDFAndLogPDF()
        {
            double sum = 0d;
            for (double z = -4d; z <= 4d; z += 1d) sum += _pearson.LogCDF(z) + _pearson.LogPDF(z);
            return sum;
        }

        /// <summary>The frequency factor at seven probabilities.</summary>
        /// <returns>The sum of the frequency factors.</returns>
        [Benchmark]
        public double FrequencyFactorKp()
        {
            double sum = 0d;
            foreach (double p in Probabilities) sum += GammaDistribution.FrequencyFactorKp(Skew, p);
            return sum;
        }

        /// <summary>The frequency-factor skew derivative at seven probabilities.</summary>
        /// <returns>The sum of the derivatives.</returns>
        [Benchmark]
        public double PartialKp()
        {
            double sum = 0d;
            foreach (double p in Probabilities) sum += GammaDistribution.PartialKp(Skew, p);
            return sum;
        }

        /// <summary>The log-Pearson skewness and kurtosis at log scales 0.3 and 0.02.</summary>
        /// <returns>The sum of the moments.</returns>
        [Benchmark]
        public double LogPearsonSkewnessKurtosis()
            => _logPearson.Skewness + _logPearson.Kurtosis + _smallScale.Skewness + _smallScale.Kurtosis;

        /// <summary>The skew variance of the maximum-likelihood covariance at sample size 100.</summary>
        /// <returns>The skew variance.</returns>
        [Benchmark]
        public double MaximumLikelihoodCovariance()
            => _pearson.ParameterCovariance(100, ParameterEstimationMethod.MaximumLikelihood)[2, 2];
    }
}
