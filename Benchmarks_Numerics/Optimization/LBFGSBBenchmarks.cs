using BenchmarkDotNet.Attributes;
using Numerics.Mathematics.Optimization;

namespace Benchmarks_Numerics.Optimization
{
    /// <summary>
    /// LBFGSB on the problems timed by docs/mathematics/oracles/lbfgsb/lbfgsb_scipy_timing.py, with analytic gradients and
    /// the SciPy defaults m = 10, gtol = 1e-5 and ftol = 2.220446049250313e-9.
    /// </summary>
    [MemoryDiagnoser]
    public class LBFGSBBenchmarks
    {
        private Func<double[], double> _f = null!;
        private Func<double[], double[]> _g = null!;
        private double[] _x0 = null!, _lower = null!, _upper = null!;

        /// <summary>
        /// The problem: Rosenbrock in 10 or 100 dimensions without bounds, Rosenbrock in 100 dimensions on [-2, 0.9],
        /// or a separable weighted quadratic in 1000 dimensions on [-1, 1] with most bounds active.
        /// </summary>
        [Params("Rosenbrock10", "Rosenbrock100", "BoxRosenbrock100", "Quadratic1000")]
        public string Problem { get; set; } = "";

        /// <summary>
        /// Builds the current problem.
        /// </summary>
        [GlobalSetup]
        public void Setup()
        {
            int n = Problem switch { "Rosenbrock10" => 10, "Quadratic1000" => 1000, _ => 100 };
            _x0 = new double[n];
            _lower = new double[n];
            _upper = new double[n];
            if (Problem == "Quadratic1000")
            {
                var w = new double[n];
                var c = new double[n];
                for (int i = 0; i < n; i++)
                {
                    w[i] = 1 + i % 10;
                    c[i] = 1.5 * Math.Sin(i);
                    _lower[i] = -1;
                    _upper[i] = 1;
                }
                _f = x => { double s = 0; for (int i = 0; i < x.Length; i++) s += w[i] * (x[i] - c[i]) * (x[i] - c[i]); return s; };
                _g = x => { var g = new double[x.Length]; for (int i = 0; i < x.Length; i++) g[i] = 2 * w[i] * (x[i] - c[i]); return g; };
                return;
            }
            bool boxed = Problem == "BoxRosenbrock100";
            for (int i = 0; i < n; i++)
            {
                _x0[i] = i % 2 == 0 ? -1.2 : 1;
                _lower[i] = boxed ? -2 : double.NegativeInfinity;
                _upper[i] = boxed ? 0.9 : double.PositiveInfinity;
            }
            _f = x =>
            {
                double s = 0;
                for (int i = 0; i < x.Length - 1; i++)
                {
                    double a = x[i + 1] - x[i] * x[i], b = 1 - x[i];
                    s += 100 * a * a + b * b;
                }
                return s;
            };
            _g = x =>
            {
                var g = new double[x.Length];
                for (int i = 0; i < x.Length - 1; i++)
                {
                    double a = x[i + 1] - x[i] * x[i];
                    g[i] += -400 * x[i] * a - 2 * (1 - x[i]);
                    g[i + 1] += 200 * a;
                }
                return g;
            };
        }

        /// <summary>
        /// Minimizes the current problem.
        /// </summary>
        [Benchmark]
        public double Minimize()
        {
            var solver = new LBFGSB(_f, _x0.Length, _x0, _lower, _upper, _g)
            {
                AbsoluteTolerance = 1E-5,
                RelativeTolerance = 2.220446049250313E-9,
                MaxIterations = 15000,
                RecordTraces = false,
                ComputeHessian = false,
                ReportFailure = false
            };
            solver.Minimize();
            return solver.BestParameterSet.Fitness;
        }
    }
}
