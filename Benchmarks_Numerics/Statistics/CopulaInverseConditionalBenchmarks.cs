using BenchmarkDotNet.Attributes;
using Numerics.Distributions.Copulas;
using Numerics.Mathematics.RootFinding;

namespace Benchmarks_Numerics.Statistics
{
    /// <summary>
    /// The Gumbel and Joe inverse conditional CDFs (safeguarded Newton on the closed-form conditional, full precision)
    /// against a reference implementation of the earlier Brent inversion at its default tolerance of 1E-8, over 64
    /// (u, t) pairs at θ = 2.
    /// </summary>
    [ShortRunJob]
    public class CopulaInverseConditionalBenchmarks
    {
        private readonly double[] _u = new double[64], _t = new double[64];
        private readonly GumbelCopula _gumbel = new GumbelCopula(2d);
        private readonly JoeCopula _joe = new JoeCopula(2d);

        /// <summary>Fills the (u, t) pairs on a scrambled grid of (0, 1).</summary>
        [GlobalSetup]
        public void Setup()
        {
            for (int i = 0; i < _u.Length; i++)
            {
                _u[i] = (i + 0.5) / _u.Length;
                _t[i] = ((i * 37) % _t.Length + 0.5) / _t.Length;
            }
        }

        private static double GumbelBrent(double theta, double u, double t)
        {
            Func<double, double> f = x =>
            {
                double vu = Math.Pow(-Math.Log(u), theta - 1d) * Math.Exp(-Math.Pow(Math.Pow(-Math.Log(u), theta) + Math.Pow(-Math.Log(x), theta), 1d / theta)) * Math.Pow(Math.Pow(-Math.Log(u), theta) + Math.Pow(-Math.Log(x), theta), 1d / theta - 1d) / u;
                return vu - t;
            };
            if (f(1d) <= 0d) return 1d;
            return Brent.Solve(f, 0d, 1d);
        }

        private static double JoeBrent(double theta, double u, double t)
        {
            Func<double, double> f = x =>
            {
                double vu = -(Math.Pow(1d - x, theta) - 1d) * Math.Pow(Math.Pow(1d - u, theta) - Math.Pow(1d - u, theta) * Math.Pow(1d - x, theta) + Math.Pow(1d - x, theta), (-theta + 1d) / theta) * Math.Pow(1d - u, theta - 1d);
                return vu - t;
            };
            if (f(1d) <= 0d) return 1d;
            return Brent.Solve(f, 0d, 1d);
        }

        /// <summary>The earlier Gumbel Brent inversion.</summary>
        /// <returns>The sum of the roots.</returns>
        [Benchmark(Baseline = true)]
        public double GumbelBrentPrevious()
        {
            double s = 0d;
            for (int i = 0; i < _u.Length; i++) s += GumbelBrent(2d, _u[i], _t[i]);
            return s;
        }

        /// <summary>The Gumbel Newton inversion.</summary>
        /// <returns>The sum of the roots.</returns>
        [Benchmark]
        public double GumbelNewton()
        {
            double s = 0d;
            for (int i = 0; i < _u.Length; i++) s += _gumbel.InverseConditionalCDF(_u[i], _t[i]);
            return s;
        }

        /// <summary>The earlier Joe Brent inversion.</summary>
        /// <returns>The sum of the roots.</returns>
        [Benchmark]
        public double JoeBrentPrevious()
        {
            double s = 0d;
            for (int i = 0; i < _u.Length; i++) s += JoeBrent(2d, _u[i], _t[i]);
            return s;
        }

        /// <summary>The Joe Newton inversion.</summary>
        /// <returns>The sum of the roots.</returns>
        [Benchmark]
        public double JoeNewton()
        {
            double s = 0d;
            for (int i = 0; i < _u.Length; i++) s += _joe.InverseConditionalCDF(_u[i], _t[i]);
            return s;
        }
    }
}
