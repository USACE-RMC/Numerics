using System;
using System.Linq;
using System.Threading;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Differentiation
{
    /// <summary>
    /// Unit tests for how the bounded numerical derivatives treat exceptions from the function.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// Cancellation and an optimizer's evaluation or iteration budget propagate; any other failed evaluation is retried
    /// with a smaller step. Gradient-based optimizers therefore stop at exactly their evaluation budget.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_DerivativeStopRequests
    {
        private static readonly double[] Point = { 0.3, -0.7 };
        private static readonly double[] Lower = { -5, -5 };
        private static readonly double[] Upper = { 5, 5 };

        private static double Quadratic(double[] x) => x[0] * x[0] + 3 * x[1] * x[1] + x[0] * x[1];

        /// <summary>
        /// A function that succeeds for the first <paramref name="good"/> calls and then throws.
        /// </summary>
        private static Func<double[], double> ThrowAfter(int good, Func<Exception> exception)
        {
            int calls = 0;
            return x =>
            {
                if (++calls > good) throw exception();
                return Quadratic(x);
            };
        }

        /// <summary>
        /// Cancellation from the function propagates out of the gradient, Jacobian and Hessian.
        /// </summary>
        [TestMethod]
        public void Test_Cancellation_Propagates()
        {
            Func<Exception> cancel = () => new OperationCanceledException();
            Assert.ThrowsExactly<OperationCanceledException>(() => NumericalDerivative.Gradient(ThrowAfter(2, cancel), Point, Lower, Upper));
            Assert.ThrowsExactly<OperationCanceledException>(() => NumericalDerivative.Hessian(ThrowAfter(3, cancel), Point, Lower, Upper));
            var f = ThrowAfter(2, cancel);
            Assert.ThrowsExactly<OperationCanceledException>(() => NumericalDerivative.Jacobian(x => new[] { f(x), x[0] }, Point, Lower, Upper));
        }

        /// <summary>
        /// An optimizer's budget stop, an argument exception for the evaluation or iteration limit, propagates.
        /// </summary>
        [TestMethod]
        public void Test_BudgetStop_Propagates()
        {
            foreach (var name in new[] { nameof(Optimizer.MaxFunctionEvaluations), nameof(Optimizer.MaxIterations) })
            {
                Func<Exception> stop = () => new ArgumentException("The budget has been reached.", name);
                var ex = Assert.ThrowsExactly<ArgumentException>(() => NumericalDerivative.Gradient(ThrowAfter(1, stop), Point, Lower, Upper));
                Assert.AreEqual(name, ex.ParamName);
                Assert.ThrowsExactly<ArgumentException>(() => NumericalDerivative.Hessian(ThrowAfter(4, stop), Point, Lower, Upper));
            }
        }

        /// <summary>
        /// Any other exception is still treated as a failed evaluation: the step is reduced and the evaluation retried.
        /// </summary>
        [TestMethod]
        public void Test_OtherExceptions_Retried()
        {
            int calls = 0;
            // The second evaluation fails once; the derivative retries with a smaller step and still succeeds.
            double F(double[] x)
            {
                if (++calls == 2) throw new InvalidOperationException("transient");
                return Quadratic(x);
            }
            var g = NumericalDerivative.Gradient(F, Point, Lower, Upper);
            Assert.AreEqual(2 * 0.3 - 0.7, g[0], 1E-5);
            Assert.AreEqual(6 * -0.7 + 0.3, g[1], 1E-5);

            Func<Exception> other = () => new ArgumentException("bad value", "x");
            var zero = NumericalDerivative.Gradient(ThrowAfter(1, other), Point, Lower, Upper);
            Assert.IsTrue(zero.All(v => v == 0.0));
        }

        /// <summary>
        /// Gradient descent and ADAM stop at exactly the evaluation budget, and cancellation stops them.
        /// </summary>
        [TestMethod]
        public void Test_GradientOptimizers_StopAtBudget()
        {
            double Rosenbrock(double[] x)
            {
                double s = 0;
                for (int i = 0; i < x.Length - 1; i++) s += 100 * Math.Pow(x[i + 1] - x[i] * x[i], 2) + Math.Pow(1 - x[i], 2);
                return s;
            }
            var x0 = Enumerable.Repeat(-1.2, 5).ToArray();
            var lo = Enumerable.Repeat(-5.0, 5).ToArray();
            var hi = Enumerable.Repeat(5.0, 5).ToArray();
            foreach (int budget in new[] { 23, 101 })
            {
                var gd = new GradientDescent(Rosenbrock, 5, x0, lo, hi) { MaxFunctionEvaluations = budget, ComputeHessian = false };
                gd.Minimize();
                Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, gd.Status);
                Assert.AreEqual(budget, gd.FunctionEvaluations);
                var adam = new ADAM(Rosenbrock, 5, x0, lo, hi) { MaxFunctionEvaluations = budget, ComputeHessian = false };
                adam.Minimize();
                Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, adam.Status);
                Assert.AreEqual(budget, adam.FunctionEvaluations);
            }

            var source = new CancellationTokenSource();
            int calls = 0;
            var cancelled = new GradientDescent(x => { if (++calls == 30) source.Cancel(); return Rosenbrock(x); }, 5, x0, lo, hi)
            {
                CancellationToken = source.Token,
                ComputeHessian = false
            };
            cancelled.Minimize();
            Assert.AreEqual(OptimizationStatus.Cancelled, cancelled.Status);
            Assert.AreEqual(30, calls);
        }
    }
}
