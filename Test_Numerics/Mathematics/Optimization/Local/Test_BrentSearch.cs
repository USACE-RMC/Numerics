using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// Unit tests for the Brent optimization algorithm
    /// </summary>
    [TestClass]
    public class Test_BrentSearch
    {
        /// <summary>
        /// Test to find the minimum of a one dimensional function using Brent's method
        /// </summary>
        [TestMethod]
        public void Test_Minimize()
        {
            double lower = -3d;
            double upper = 3d;
            var solver = new BrentSearch(TestFunctions.FX, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = 1.0d;
            Assert.AreEqual(X, trueX, 1E-4);
        }

        /// <summary>
        /// Test to find the maximum of a one dimensional function using Brent's method
        /// </summary>
        [TestMethod]
        public void Test_Maximize()
        {
            double lower = -3;
            double upper = 3d;
            var solver = new BrentSearch(TestFunctions.FX, lower, upper);
            solver.Maximize();
            double F = -1*solver.BestParameterSet.Fitness;
            double trueF = 9.4815;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = -1.6667d;
            Assert.AreEqual(X, trueX, 1E-4);
        }

        /// <summary>
        /// Test the Brent algorithm with De Jong's function in 1-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            double lower = -5.12d;
            double upper = 5.12d;
            var solver = new BrentSearch((x) => { return TestFunctions.DeJong(new double[] { x }); }, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = 0.0;
            Assert.AreEqual(X, trueX, 1E-4);
        }

        /// <summary>
        /// Verifies geometric expansion to the right and the exact trial sequence.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_GeometricExpansionRight()
        {
            var evaluations = new List<double>();
            var solver = new BrentSearch(x =>
            {
                evaluations.Add(x);
                return (x - 10d) * (x - 10d);
            }, 0d, 1d);

            solver.Bracket(1d, 2d);

            CollectionAssert.AreEqual(new[] { 0d, 1d, 2d, 4d, 8d, 16d }, evaluations.ToArray());
            Assert.AreEqual(4d, solver.LowerBound, 0d);
            Assert.AreEqual(16d, solver.UpperBound, 0d);
        }

        /// <summary>
        /// Verifies that a positive initial step reverses direction and expands geometrically when the minimum is to the left.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_GeometricExpansionLeft()
        {
            var evaluations = new List<double>();
            var solver = new BrentSearch(x =>
            {
                evaluations.Add(x);
                return (x + 10d) * (x + 10d);
            }, 0d, 1d);

            solver.Bracket(1d, 2d);

            CollectionAssert.AreEqual(new[] { 0d, 1d, -1d, -3d, -7d, -15d }, evaluations.ToArray());
            Assert.AreEqual(-15d, solver.LowerBound, 0d);
            Assert.AreEqual(-3d, solver.UpperBound, 0d);
        }

        /// <summary>
        /// Verifies that a caller-provided expansion factor controls the geometric trial sequence.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_UsesCustomExpansionFactor()
        {
            var evaluations = new List<double>();
            var solver = new BrentSearch(x =>
            {
                evaluations.Add(x);
                return (x - 10d) * (x - 10d);
            }, 0d, 1d);

            solver.Bracket(1d, 3d);

            CollectionAssert.AreEqual(new[] { 0d, 1d, 2d, 5d, 14d, 41d }, evaluations.ToArray());
            Assert.AreEqual(5d, solver.LowerBound, 0d);
            Assert.AreEqual(41d, solver.UpperBound, 0d);
        }

        /// <summary>
        /// Verifies that geometric bracketing finds a distant minimum with logarithmic objective work and remains compatible with minimization.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_DistantMinimumReducesEvaluations()
        {
            int evaluations = 0;
            var solver = new BrentSearch(x =>
            {
                evaluations++;
                return (x - 1000d) * (x - 1000d);
            }, 0d, 1d);

            solver.Bracket(0.1d, 2d);
            int bracketEvaluations = evaluations;

            Assert.IsLessThanOrEqualTo(20, bracketEvaluations, $"Expected at most 20 bracket evaluations, but observed {bracketEvaluations}.");
            Assert.IsTrue(solver.LowerBound <= 1000d && solver.UpperBound >= 1000d);
            solver.Minimize();
            Assert.AreEqual(1000d, solver.BestParameterSet.Values[0], 1E-4);
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-8);
        }

        /// <summary>
        /// Verifies that a non-strict bracket terminates immediately for a flat objective.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_PlateauTerminates()
        {
            int evaluations = 0;
            var solver = new BrentSearch(x =>
            {
                evaluations++;
                return 1d;
            }, 0d, 1d);

            solver.Bracket(1d, 2d);

            Assert.AreEqual(3, evaluations);
            Assert.AreEqual(0d, solver.LowerBound, 0d);
            Assert.AreEqual(2d, solver.UpperBound, 0d);
        }

        /// <summary>
        /// Verifies bounded failure for a monotone objective and transactional preservation of the original bounds.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_MonotoneObjectiveReachesIterationLimit()
        {
            var solver = new BrentSearch(x => -x, 0d, 1d)
            {
                MaxIterations = 10,
                ReportFailure = false
            };

            solver.Bracket(1d, 2d);

            Assert.AreEqual(OptimizationStatus.MaximumIterationsReached, solver.Status);
            Assert.AreEqual(0d, solver.LowerBound, 0d);
            Assert.AreEqual(1d, solver.UpperBound, 0d);

            var throwingSolver = new BrentSearch(x => -x, 0d, 1d) { MaxIterations = 10 };
            var exception = AssertThrows<ArgumentException>(() => throwingSolver.Bracket(1d, 2d));
            Assert.AreEqual(nameof(Optimizer.MaxIterations), exception.ParamName);
        }

        /// <summary>
        /// Verifies validation of the public bracketing step and expansion-factor contract.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_RejectsInvalidInputs()
        {
            foreach (double invalidStep in new[] { 0d, double.NaN, double.NegativeInfinity, double.PositiveInfinity })
            {
                var solver = new BrentSearch(x => x * x, 0d, 1d);
                var exception = AssertThrows<ArgumentOutOfRangeException>(() => solver.Bracket(invalidStep, 2d));
                Assert.AreEqual("s", exception.ParamName);
            }

            foreach (double invalidExpansion in new[] { -1d, 0d, 1d, double.NaN, double.NegativeInfinity, double.PositiveInfinity })
            {
                var solver = new BrentSearch(x => x * x, 0d, 1d);
                var exception = AssertThrows<ArgumentOutOfRangeException>(() => solver.Bracket(1d, invalidExpansion));
                Assert.AreEqual("k", exception.ParamName);
            }
        }

        /// <summary>
        /// Verifies deterministic failure for NaN objectives and coordinate overflow without changing the original bounds.
        /// </summary>
        [TestMethod]
        public void Test_Bracket_RejectsNonFiniteSearchState()
        {
            var nanSolver = new BrentSearch(x => double.NaN, 0d, 1d) { ReportFailure = false };
            nanSolver.Bracket();
            Assert.AreEqual(OptimizationStatus.Failure, nanSolver.Status);
            Assert.AreEqual(0d, nanSolver.LowerBound, 0d);
            Assert.AreEqual(1d, nanSolver.UpperBound, 0d);

            var overflowSolver = new BrentSearch(x => -x, double.MaxValue, double.MaxValue) { ReportFailure = false };
            overflowSolver.Bracket(double.MaxValue, 2d);
            Assert.AreEqual(OptimizationStatus.Failure, overflowSolver.Status);
            Assert.AreEqual(double.MaxValue, overflowSolver.LowerBound);
            Assert.AreEqual(double.MaxValue, overflowSolver.UpperBound);

            var throwingSolver = new BrentSearch(x => -x, double.MaxValue, double.MaxValue);
            AssertThrows<ArithmeticException>(() => throwingSolver.Bracket(double.MaxValue, 2d));
        }

        /// <summary>
        /// Executes an action and returns the expected exception.
        /// </summary>
        /// <typeparam name="TException">The expected exception type.</typeparam>
        /// <param name="action">The action expected to throw.</param>
        /// <returns>The exception thrown by <paramref name="action"/>.</returns>
        private static TException AssertThrows<TException>(Action action) where TException : Exception
        {
            try
            {
                action();
            }
            catch (TException exception)
            {
                return exception;
            }
            catch (Exception exception)
            {
                Assert.Fail($"Expected {typeof(TException).Name}, but observed {exception.GetType().Name}.");
            }

            Assert.Fail($"Expected {typeof(TException).Name}, but no exception was thrown.");
            throw new InvalidOperationException("The exception assertion did not terminate as expected.");
        }


        /// <summary>
        /// A run is reproducible bit for bit.
        /// </summary>
        /// <remarks>
        /// The fingerprint pins the exact search path; see <see cref="OptimizerFingerprint"/>. A different
        /// fingerprint means different results.
        /// </remarks>
        [TestMethod]
        public void Test_SeededRunFingerprint()
        {
            var solver = new BrentSearch(x => (x - 1) * (x - 1) * (x + 2), -1.5, 4) { ComputeHessian = false };
            solver.Minimize();
            Assert.AreEqual("x=3FF00000003F52A8|f=3C477EAE10345EBA|it=12|fe=13|st=Success|trace=13:D01A04FA0F958673", OptimizerFingerprint.Of(solver));
        }

        /// <summary>
        /// Iterations counts the search iterations: one evaluation each after the initial one.
        /// </summary>
        [TestMethod]
        public void Test_IterationsAreCounted()
        {
            var solver = new BrentSearch(x => (x - 1) * (x - 1), -3, 5);
            solver.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, solver.Status);
            Assert.IsGreaterThan(0, solver.Iterations);
            Assert.AreEqual(solver.FunctionEvaluations - 1, solver.Iterations);
        }

        /// <summary>
        /// At SciPy's settings (RelativeTolerance = sqrt(2.2E-16), AbsoluteTolerance = xatol = 1E-5, MaxFunctionEvaluations =
        /// maxiter = 500) the minimizer, its value and the evaluation count are bitwise those of
        /// scipy.optimize.minimize_scalar(f, bounds, method='bounded') (SciPy 1.18.1), for smooth, non-smooth, flat and
        /// maximization objectives.
        /// </summary>
        [TestMethod]
        public void Test_BitwiseScipyMinimizeScalarBounded()
        {
            var cases = new (Func<double, double> F, double Lower, double Upper, bool Maximize, double X, double Fx, int Evaluations)[]
            {
                (x => (x + 3.0) * ((x - 1.0) * (x - 1.0)), -3.0, 3.0, false, 1.0000002210646803, 1.95478382377075e-13, 11),
                (x => x * x, -5.12, 5.12, false, 0.0, 0.0, 6),
                (x => x * x * x * x - 3.0 * x * x * x + 2.0, 0.0, 5.0, false, 2.2499994695548797, -6.542968749997151, 12),
                (x => Math.Abs(x - 0.3), -1.0, 1.0, false, 0.29999957916974357, 4.2083025642014604e-07, 19),
                (x => (x - 1000.0) * (x - 1000.0), 0.0, 2000.0, false, 1000.0, 0.0, 6),
                (x => Math.Cos(x), 0.0, 6.0, false, 3.1415926891518526, -0.9999999999999993, 9),
                // Maximize is the minimization of the negated objective: SciPy minimizes -(x + 3)(x - 1)^2.
                (x => (x + 3.0) * ((x - 1.0) * (x - 1.0)), -3.0, 3.0, true, -1.6666659774370731, -9.48148148147958, 11),
            };
            for (int i = 0; i < cases.Length; i++)
            {
                var c = cases[i];
                var solver = new BrentSearch(c.F, c.Lower, c.Upper)
                {
                    RelativeTolerance = Math.Sqrt(2.2E-16),
                    AbsoluteTolerance = 1E-5,
                    MaxFunctionEvaluations = 500,
                    ComputeHessian = false
                };
                if (c.Maximize) solver.Maximize(); else solver.Minimize();
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(c.X), BitConverter.DoubleToInt64Bits(solver.BestParameterSet.Values[0]), $"case {i}: x");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(c.Fx), BitConverter.DoubleToInt64Bits(solver.BestParameterSet.Fitness), $"case {i}: f");
                Assert.AreEqual(c.Evaluations, solver.FunctionEvaluations, $"case {i}: evaluations");
                Assert.AreEqual(OptimizationStatus.Success, solver.Status, $"case {i}: status");
            }
        }

        /// <summary>
        /// SciPy's maxiter is the evaluation budget: with maxiter = 12 on |x - 0.3| SciPy stops at x = 0.3002178696074389 with
        /// status 1 after 12 evaluations, and so does the search with MaxFunctionEvaluations = 12.
        /// </summary>
        [TestMethod]
        public void Test_ScipyEvaluationBudget()
        {
            var solver = new BrentSearch(x => Math.Abs(x - 0.3), -1.0, 1.0)
            {
                RelativeTolerance = Math.Sqrt(2.2E-16),
                AbsoluteTolerance = 1E-5,
                MaxFunctionEvaluations = 12,
                ReportFailure = false,
                ComputeHessian = false
            };
            solver.Minimize();
            Assert.AreEqual(BitConverter.DoubleToInt64Bits(0.3002178696074389), BitConverter.DoubleToInt64Bits(solver.BestParameterSet.Values[0]));
            Assert.AreEqual(12, solver.FunctionEvaluations);
            Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, solver.Status);
        }

        /// <summary>
        /// The constructor sets the absolute tolerance so the step tolerance is 1E-8 * |x| + 1E-3 * DoubleMachineEpsilon, the
        /// guard against seeking relative accuracy at an exact zero; a minimum at zero is found exactly and one at 1E-12 to
        /// the absolute guard.
        /// </summary>
        [TestMethod]
        public void Test_DefaultTolerances()
        {
            var solver = new BrentSearch(x => x * x * x * x, -1d, 2d);
            Assert.AreEqual(3d * (Numerics.Tools.DoubleMachineEpsilon * 1E-3), solver.AbsoluteTolerance, 0d);
            Assert.AreEqual(1E-8, solver.RelativeTolerance, 0d);
            Assert.AreEqual(Numerics.Tools.DoubleMachineEpsilon * 1E-3, solver.AbsoluteTolerance / 3d, 0d);
            solver.Minimize();
            Assert.AreEqual(0d, solver.BestParameterSet.Values[0], 1E-20);
            var shifted = new BrentSearch(x => (x - 1E-12) * (x - 1E-12), -1d, 1d);
            shifted.Minimize();
            Assert.AreEqual(1E-12, shifted.BestParameterSet.Values[0], 1E-19);
        }
    }
}
