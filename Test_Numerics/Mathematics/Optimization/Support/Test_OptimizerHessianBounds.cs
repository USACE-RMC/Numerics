using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the end-of-run finite-difference Hessian computed by <see cref="Optimizer"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// A successful <see cref="Optimizer.Minimize"/> or <see cref="Optimizer.Maximize"/> differences the
    /// objective function about the reported solution. When the solution lies on a bound, an unclamped
    /// perturbation evaluates the objective outside the declared feasible region, where it may be
    /// undefined. These tests isolate that pass: the solver below performs no search of its own, so every
    /// objective evaluation after the first is a Hessian probe.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_OptimizerHessianBounds
    {
        /// <summary>
        /// An optimizer that reports a caller-supplied solution and performs no search.
        /// </summary>
        /// <remarks>
        /// This exists so the end-of-run Hessian pass can be observed on its own. A real solver interleaves
        /// its own finite-difference gradients and its own iterate repair with the Hessian pass, which makes
        /// it impossible to attribute an out-of-bounds evaluation to one path or the other.
        /// </remarks>
        private sealed class FixedSolutionOptimizer : Optimizer
        {
            private readonly double[] _solution;

            /// <summary>
            /// Construct a solver that reports <paramref name="solution"/> as a successful result.
            /// </summary>
            /// <param name="objectiveFunction">The objective function to evaluate.</param>
            /// <param name="solution">The parameter set to report.</param>
            /// <param name="lowerBounds">The inclusive lower bounds, or null for an unbounded problem.</param>
            /// <param name="upperBounds">The inclusive upper bounds, or null for an unbounded problem.</param>
            public FixedSolutionOptimizer(Func<double[], double> objectiveFunction, double[] solution, double[] lowerBounds, double[] upperBounds)
                : base(objectiveFunction, solution.Length)
            {
                _solution = solution;
                LowerBounds = lowerBounds;
                UpperBounds = upperBounds;
            }

            /// <summary>
            /// An array of lower bounds (inclusive) of the interval containing the optimal point.
            /// </summary>
            public double[] LowerBounds { get; }

            /// <summary>
            /// An array of upper bounds (inclusive) of the interval containing the optimal point.
            /// </summary>
            public double[] UpperBounds { get; }

            /// <inheritdoc />
            protected override double[] ParameterLowerBounds => LowerBounds;

            /// <inheritdoc />
            protected override double[] ParameterUpperBounds => UpperBounds;

            /// <inheritdoc />
            protected override void Optimize()
            {
                bool cancel = false;
                Evaluate((double[])_solution.Clone(), ref cancel);
                UpdateStatus(OptimizationStatus.Success);
            }
        }

        /// <summary>
        /// A separable quadratic whose exact Hessian is the constant diagonal matrix diag(2a, 2b).
        /// </summary>
        private static Func<double[], double> Quadratic(double a, double b, List<string> violations, double[] lower, double[] upper)
        {
            return x =>
            {
                if (violations != null)
                {
                    for (int i = 0; i < x.Length; i++)
                    {
                        if (x[i] < lower[i] || x[i] > upper[i])
                            violations.Add("p" + i + "=" + x[i].ToString("R"));
                    }
                }
                return a * x[0] * x[0] + b * x[1] * x[1];
            };
        }

        /// <summary>
        /// Verifies that the end-of-run Hessian never evaluates the objective outside the declared bounds
        /// when the reported solution sits on a corner of the feasible region.
        /// </summary>
        /// <remarks>
        /// Both parameters are on a bound, so an unclamped central difference steps outside in every
        /// coordinate and in all four off-diagonal stencil directions.
        /// </remarks>
        [TestMethod]
        public void Test_HessianDoesNotProbeOutsideBounds()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var violations = new List<string>();
            var solver = new FixedSolutionOptimizer(Quadratic(3d, 5d, violations, lower, upper), new double[] { 1d, 0d }, lower, upper);

            solver.Minimize();

            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsNotNull(solver.Hessian, "The Hessian was not computed.");
            Assert.IsEmpty(violations, "The objective was evaluated outside the declared bounds: " + string.Join(", ", violations));
        }

        /// <summary>
        /// Verifies that clamping the Hessian probes into the feasible region still recovers the exact
        /// curvature of a quadratic, so the guard does not cost accuracy.
        /// </summary>
        /// <remarks>
        /// A one-sided second difference is exact for a quadratic up to floating point rounding, so the
        /// tolerance here only absorbs cancellation in the difference quotient.
        /// </remarks>
        [TestMethod]
        public void Test_BoundedHessianRecoversQuadraticCurvature()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var solver = new FixedSolutionOptimizer(Quadratic(3d, 5d, null, lower, upper), new double[] { 1d, 0d }, lower, upper);

            solver.Minimize();

            Assert.IsNotNull(solver.Hessian);
            Assert.AreEqual(6d, solver.Hessian[0, 0], 1E-4);
            Assert.AreEqual(10d, solver.Hessian[1, 1], 1E-4);
            Assert.AreEqual(0d, solver.Hessian[0, 1], 1E-4);
            Assert.AreEqual(0d, solver.Hessian[1, 0], 1E-4);
        }

        /// <summary>
        /// Verifies that an optimizer which declares no bounds still computes its Hessian.
        /// </summary>
        /// <remarks>
        /// <see cref="Optimizer"/> returns null bounds by default, which is the correct description of an
        /// unconstrained problem. The finite-difference routine treats null as unlimited room, so the
        /// unbounded result is a plain central difference.
        /// </remarks>
        [TestMethod]
        public void Test_UnboundedOptimizerComputesHessian()
        {
            var solver = new FixedSolutionOptimizer(x => 3d * x[0] * x[0] + 5d * x[1] * x[1], new double[] { 0.5d, 0.25d }, null, null);

            solver.Minimize();

            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsNotNull(solver.Hessian);
            Assert.AreEqual(6d, solver.Hessian[0, 0], 1E-4);
            Assert.AreEqual(10d, solver.Hessian[1, 1], 1E-4);
        }

        /// <summary>
        /// Verifies that the maximization path guards its Hessian probes as well as the minimization path.
        /// </summary>
        /// <remarks>
        /// <see cref="Optimizer.Maximize"/> carries its own copy of the end-of-run Hessian block, so it
        /// needs its own coverage; the two are easy to change out of step.
        /// </remarks>
        [TestMethod]
        public void Test_MaximizeHessianDoesNotProbeOutsideBounds()
        {
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 1d, 1d };
            var violations = new List<string>();
            var solver = new FixedSolutionOptimizer(Quadratic(-3d, -5d, violations, lower, upper), new double[] { 1d, 1d }, lower, upper);

            solver.Maximize();

            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsNotNull(solver.Hessian, "The Hessian was not computed.");
            Assert.IsEmpty(violations, "The objective was evaluated outside the declared bounds: " + string.Join(", ", violations));
        }
    }
}
