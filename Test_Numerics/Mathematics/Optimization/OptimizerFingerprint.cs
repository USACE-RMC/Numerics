using System;
using System.Linq;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// A bitwise fingerprint of an optimizer run, used to pin a seeded search path exactly.
    /// </summary>
    /// <remarks>
    /// <para>
    /// The fingerprint holds the bit patterns of the best parameter set and fitness, the iteration and
    /// evaluation counters, the status, and an FNV-1a hash over the bit patterns of every trace entry.
    /// Two runs with equal fingerprints followed the same search path, evaluation for evaluation.
    /// </para>
    /// <para>
    /// The objectives below use only addition and multiplication, so their values are identical on
    /// every target framework, and a test can pin a run with an expected string.
    /// </para>
    /// </remarks>
    internal static class OptimizerFingerprint
    {
        /// <summary>
        /// The Rosenbrock function written with products only.
        /// </summary>
        /// <param name="x">The point.</param>
        /// <returns>The function value.</returns>
        public static double Rosenbrock(double[] x)
        {
            double s = 0;
            for (int i = 0; i < x.Length - 1; i++)
            {
                double a = x[i + 1] - x[i] * x[i], b = 1 - x[i];
                s += 100 * a * a + b * b;
            }
            return s;
        }

        /// <summary>
        /// The six-hump camel function written with products only.
        /// </summary>
        /// <param name="v">The point.</param>
        /// <returns>The function value.</returns>
        public static double SixHumpCamel(double[] v)
        {
            double x = v[0], y = v[1], x2 = x * x, y2 = y * y;
            return (4 - 2.1 * x2 + x2 * x2 / 3) * x2 + x * y + (-4 + 4 * y2) * y2;
        }

        /// <summary>
        /// Returns the fingerprint of a finished run.
        /// </summary>
        /// <param name="optimizer">The optimizer after Minimize or Maximize.</param>
        /// <returns>The fingerprint string.</returns>
        public static string Of(Optimizer optimizer)
        {
            ulong hash = 14695981039346656037UL;
            void Mix(double d)
            {
                hash ^= (ulong)BitConverter.DoubleToInt64Bits(d);
                hash *= 1099511628211UL;
            }
            foreach (var p in optimizer.ParameterSetTrace)
            {
                Mix(p.Fitness);
                foreach (var v in p.Values) Mix(v);
            }
            string values = string.Join(",", optimizer.BestParameterSet.Values.Select(v => BitConverter.DoubleToInt64Bits(v).ToString("X16")));
            return $"x={values}|f={BitConverter.DoubleToInt64Bits(optimizer.BestParameterSet.Fitness):X16}|it={optimizer.Iterations}|fe={optimizer.FunctionEvaluations}|st={optimizer.Status}|trace={optimizer.ParameterSetTrace.Count}:{hash:X16}";
        }
    }
}
