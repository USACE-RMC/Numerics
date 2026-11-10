using BenchmarkDotNet.Attributes;
using Numerics.Distributions;
using Numerics.Mathematics.RootFinding;

namespace Benchmarks_Numerics.Distributions
{
    /// <summary>
    /// The generalized Pareto shape solve of the method of moments on 1,000 skews, at the solver tolerances of
    /// <see cref="GeneralizedPareto.SolveForKappa(double)"/> (absolute 1E-15, the minimum relative tolerance) against
    /// the same Brent solve with an absolute tolerance of 1E-8.
    /// </summary>
    /// <remarks>
    /// The skews are those of shapes spaced evenly on [-0.25, 2], which covers both brackets of the solve and the
    /// neighborhood of the exponential shape zero.
    /// </remarks>
    [MemoryDiagnoser]
    [ShortRunJob]
    public class GeneralizedParetoMomentBenchmarks
    {
        private const int N = 1000;
        private double[] _skews = null!;

        /// <summary>
        /// Prepares the skews.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            _skews = new double[N];
            for (int i = 0; i < N; i++)
            {
                double k = -0.25 + 2.25 * i / (N - 1);
                _skews[i] = Skewness(k);
            }
        }

        /// <summary>The skewness of the generalized Pareto distribution at shape k.</summary>
        /// <param name="k">The shape.</param>
        /// <returns>The skewness.</returns>
        private static double Skewness(double k) => 2d * (1d - k) * Math.Sqrt(1d + 2d * k) / (1d + 3d * k);

        /// <summary>The same Brent solve with an absolute tolerance of 1E-8.</summary>
        /// <returns>The sum of the shapes.</returns>
        [Benchmark(Baseline = true)]
        public double AbsoluteTolerance1E8()
        {
            double sum = 0d, third = 1d / 3d, boundary = Skewness(third);
            for (int i = 0; i < N; i++)
            {
                double skew = _skews[i];
                sum += skew >= boundary
                    ? Brent.Solve(x => Skewness(x) - skew, -third, third)
                    : Brent.Solve(x => Skewness(x) - skew, third, 250d);
            }
            return sum;
        }

        /// <summary>The shape solve of <see cref="GeneralizedPareto.SolveForKappa(double)"/>.</summary>
        /// <returns>The sum of the shapes.</returns>
        [Benchmark]
        public double SolveForKappa()
        {
            var gpa = new GeneralizedPareto();
            double sum = 0d;
            for (int i = 0; i < N; i++) sum += gpa.SolveForKappa(_skews[i]);
            return sum;
        }
    }
}
