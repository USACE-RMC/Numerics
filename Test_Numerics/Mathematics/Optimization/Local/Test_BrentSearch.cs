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

    }
}
