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
    }
}
