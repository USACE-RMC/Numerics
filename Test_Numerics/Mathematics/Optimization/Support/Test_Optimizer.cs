using System.Collections.Generic;
using System.Threading;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;
using System;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the optimizer base class trace contract.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_Optimizer
    {
        /// <summary>
        /// The parameter set trace records one best-so-far entry per function evaluation with
        /// nonincreasing fitness, ends at the best parameter set, and shares one values array
        /// across the entries recorded between improvements.
        /// </summary>
        /// <remarks>
        /// The trace is read-only: entries recorded between improvements alias the same values
        /// array, so the trace stores the search history without one array allocation per
        /// function evaluation.
        /// </remarks>
        [TestMethod]
        public void Test_ParameterSetTrace_Contract()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new NelderMead(TestFunctions.FXYZ, 3, initial, lower, upper);
            solver.Minimize();

            var trace = solver.ParameterSetTrace;
            Assert.HasCount(solver.FunctionEvaluations, trace);

            for (int i = 1; i < trace.Count; i++)
            {
                Assert.IsLessThanOrEqualTo(trace[i - 1].Fitness, trace[i].Fitness, $"fitness increased at entry {i}");
                if (trace[i].Fitness == trace[i - 1].Fitness)
                {
                    Assert.IsTrue(ReferenceEquals(trace[i].Values, trace[i - 1].Values), $"entries {i - 1} and {i} hold separate copies of one best-so-far point");
                }
            }

            var last = trace[trace.Count - 1];
            Assert.AreEqual(solver.BestParameterSet.Fitness, last.Fitness);
            CollectionAssert.AreEqual(solver.BestParameterSet.Values, last.Values);
        }

        /// <summary>
        /// Disabling trace recording leaves the parameter set trace empty while the solution is
        /// unchanged.
        /// </summary>
        [TestMethod]
        public void Test_ParameterSetTrace_RecordTracesOff()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };

            var traced = new NelderMead(TestFunctions.FXYZ, 3, initial, lower, upper);
            traced.Minimize();
            var untraced = new NelderMead(TestFunctions.FXYZ, 3, initial, lower, upper) { RecordTraces = false };
            untraced.Minimize();

            Assert.IsEmpty(untraced.ParameterSetTrace);
            Assert.AreEqual(traced.BestParameterSet.Fitness, untraced.BestParameterSet.Fitness);
            CollectionAssert.AreEqual(traced.BestParameterSet.Values, untraced.BestParameterSet.Values);
        }

        /// <summary>
        /// A NaN objective value ranks as the worst fitness, so a search that starts on a NaN
        /// point still moves to the finite minimum instead of reporting the NaN point.
        /// </summary>
        [TestMethod]
        public void Test_NaNObjective_RanksWorst()
        {
            Func<double[], double> f = x => x[0] < 0.5 ? double.NaN : (x[0] - 2) * (x[0] - 2);
            var nm = new NelderMead(f, 1, new[] { 0.45 }, new[] { 0d }, new[] { 5d });
            nm.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, nm.Status);
            Assert.AreEqual(2d, nm.BestParameterSet.Values[0], 1E-3);
            Assert.IsFalse(double.IsNaN(nm.BestParameterSet.Fitness) || double.IsInfinity(nm.BestParameterSet.Fitness));

            // NaN members are replaced like any worse point, so the population spread becomes finite and DE converges.
            var de = new DifferentialEvolution(x => x[0] < -4.5 ? double.NaN : (x[0] - 2) * (x[0] - 2), 1, new[] { -5d }, new[] { 5d }) { PopulationSize = 10 };
            de.Minimize();
            Assert.AreEqual(OptimizationStatus.Success, de.Status);
            Assert.AreEqual(2d, de.BestParameterSet.Values[0], 1E-3);
            Assert.IsLessThan(10000, de.Iterations);
        }

        /// <summary>
        /// An exception from the objective is rethrown with its original stack trace.
        /// </summary>
        [TestMethod]
        public void Test_ObjectiveException_KeepsStackTrace()
        {
            var solver = new NelderMead(ThrowingObjective, 1, new[] { 1d }, new[] { 0d }, new[] { 2d });
            var ex = Assert.ThrowsExactly<InvalidOperationException>(() => solver.Minimize());
            Assert.Contains(nameof(ThrowingObjective), ex.StackTrace!);
            Assert.AreEqual(OptimizationStatus.Failure, solver.Status);
        }

        /// <summary>
        /// An objective that always fails.
        /// </summary>
        /// <param name="x">The point.</param>
        /// <returns>Never returns.</returns>
        private static double ThrowingObjective(double[] x)
        {
            throw new InvalidOperationException("The objective failed.");
        }

        /// <summary>
        /// Cancelling the token stops the run at the next evaluation without throwing, keeps the best point
        /// found so far, and reports Cancelled; this holds for a population method, for a quasi-Newton method
        /// whose finite-difference gradient retries failed evaluations, and for local searches nested in a
        /// global one.
        /// </summary>
        [TestMethod]
        public void Test_Cancellation()
        {
            foreach (string kind in new[] { "DE", "BFGS", "MultiStart", "DEPolish" })
            {
                using var source = new CancellationTokenSource();
                int calls = 0;
                Func<double[], double> f = x =>
                {
                    if (++calls == 300) source.Cancel();
                    return OptimizerFingerprint.Rosenbrock(x);
                };
                var lower = new[] { -5d, -5d, -5d };
                var upper = new[] { 5d, 5d, 5d };
                var start = new[] { -1.2, 1d, -1.2 };
                Optimizer solver = kind switch
                {
                    "DE" => new DifferentialEvolution(f, 3, lower, upper),
                    "BFGS" => new BFGS(f, 3, start, lower, upper) { MaxIterations = 100000 },
                    "MultiStart" => new MultiStart(f, 3, start, lower, upper),
                    _ => new DifferentialEvolution(f, 3, lower, upper) { Polish = true, MaxIterations = 10 },
                };
                solver.CancellationToken = source.Token;
                solver.Minimize();
                Assert.AreEqual(OptimizationStatus.Cancelled, solver.Status, kind);
                Assert.AreEqual(300, solver.FunctionEvaluations, kind);
                Assert.IsNull(solver.Hessian, kind);
                Assert.IsFalse(double.IsNaN(solver.BestParameterSet.Fitness), kind);
            }
        }

        /// <summary>
        /// Progress is reported every ProgressInterval evaluations and once at the end, with nonincreasing
        /// best fitness, and does not change the result.
        /// </summary>
        [TestMethod]
        public void Test_Progress()
        {
            var lower = new[] { -5d, -5d, -5d };
            var upper = new[] { 5d, 5d, 5d };
            var plain = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false };
            plain.Minimize();

            var reports = new List<OptimizationProgress>();
            var watched = new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ComputeHessian = false, ProgressInterval = 250, Progress = new Collector(reports) };
            watched.Minimize();
            Assert.AreEqual(OptimizerFingerprint.Of(plain), OptimizerFingerprint.Of(watched));
            Assert.HasCount(watched.FunctionEvaluations / 250 + 1, reports);
            for (int k = 1; k < reports.Count; k++)
                Assert.IsLessThanOrEqualTo(reports[k - 1].BestFitness, reports[k].BestFitness);
            Assert.AreEqual(watched.FunctionEvaluations, reports[reports.Count - 1].FunctionEvaluations);
            Assert.AreEqual(watched.BestParameterSet.Fitness, reports[reports.Count - 1].BestFitness);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new DifferentialEvolution(OptimizerFingerprint.Rosenbrock, 3, lower, upper) { ProgressInterval = 0 }.Minimize());
        }

        /// <summary>
        /// A synchronous progress sink.
        /// </summary>
        private sealed class Collector : IProgress<OptimizationProgress>
        {
            private readonly List<OptimizationProgress> _reports;

            /// <summary>
            /// Constructs the sink.
            /// </summary>
            /// <param name="reports">The list that receives the reports.</param>
            public Collector(List<OptimizationProgress> reports) { _reports = reports; }

            /// <inheritdoc/>
            public void Report(OptimizationProgress value) => _reports.Add(value);
        }
    }
}
