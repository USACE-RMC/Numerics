using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization;

/// <summary>Regression tests for L-BFGS-B termination, bounds and evaluation contracts.</summary>
/// <remarks>The cases mirror <see cref="Test_BFGSRegression"/>.</remarks>
[TestClass]
public class Test_LBFGSBRegression
{
    /// <summary>A small objective change must not hide a nonzero gradient, including a warm start.</summary>
    /// <param name="x">Initial first coordinate.</param>
    /// <param name="y">Initial second coordinate.</param>
    [TestMethod]
    [DataRow(1e-9, 1e-4)]
    [DataRow(1e-6, 1e-3)]
    public void ScaledQuadratic_RequiresStationarity(double x, double y)
    {
        var solver = new LBFGSB(p => 0.001 + 0.5 * (1e6 * p[0] * p[0] + p[1] * p[1]), 2,
            new[] { x, y }, new[] { -10d, -10d }, new[] { 10d, 10d },
            p => new[] { 1e6 * p[0], p[1] }) { ComputeHessian = false, RelativeTolerance = 1e-15 };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.IsLessThanOrEqualTo(solver.AbsoluteTolerance, Math.Abs(1e6 * solver.BestParameterSet.Values[0]));
        Assert.IsLessThanOrEqualTo(solver.AbsoluteTolerance, Math.Abs(solver.BestParameterSet.Values[1]));
        Assert.IsGreaterThan(0, solver.Iterations);
    }

    /// <summary>A stationary initial point succeeds without an accepted step or redundant gradient.</summary>
    [TestMethod]
    public void StationaryStart_DoesNotSearch()
    {
        int gradients = 0;
        var solver = new LBFGSB(p => p[0] * p[0], 1, new[] { 0d }, new[] { -10d }, new[] { 10d },
            p => { gradients++; return new[] { 2 * p[0] }; }) { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(0, solver.Iterations);
        Assert.AreEqual(1, solver.FunctionEvaluations);
        Assert.AreEqual(1, gradients);
    }

    /// <summary>Maximization applies the objective sign to a supplied gradient without mutating it.</summary>
    [TestMethod]
    public void Maximize_UsesSuppliedGradientSign()
    {
        var buffer = new double[1];
        var solver = new LBFGSB(p => -(p[0] - 2) * (p[0] - 2), 1,
            new[] { 0d }, new[] { -10d }, new[] { 10d },
            p => { buffer[0] = -2 * (p[0] - 2); return buffer; }) { ComputeHessian = false };
        solver.Maximize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(2d, solver.BestParameterSet.Values[0], 1e-8);
    }

    /// <summary>The standard Rosenbrock minimum must satisfy the requested gradient tolerance.</summary>
    /// <remarks>The relative-reduction test is tightened so that only the gradient test can stop the run.</remarks>
    [TestMethod]
    public void Rosenbrock_ReportsStationarySolution()
    {
        var solver = RosenbrockSolver();
        solver.RelativeTolerance = 1e-15;
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        var p = solver.BestParameterSet.Values;
        Assert.AreEqual(1d, p[0], 1e-8);
        Assert.AreEqual(1d, p[1], 1e-8);
        foreach (double g in solver.Gradient!(p)) Assert.IsLessThanOrEqualTo(solver.AbsoluteTolerance, Math.Abs(g));
    }

    /// <summary>A blocked coordinate must leave the other coordinate free to reach its constrained optimum.</summary>
    [TestMethod]
    public void BoundaryOptimum_UsesProjectedGradient()
    {
        var solver = new LBFGSB(p => Math.Pow(p[0] - 2, 2) + Math.Pow(p[1] - 0.3, 2), 2,
            new[] { 0d, 0d }, new[] { 0d, 0d }, new[] { 1d, 1d },
            p => new[] { 2 * (p[0] - 2), 2 * (p[1] - 0.3) });
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(1d, solver.BestParameterSet.Values[0]);
        Assert.AreEqual(0.3, solver.BestParameterSet.Values[1], 1e-8);
        Assert.IsNotNull(solver.Hessian);
    }

    /// <summary>A steep objective far from the start converges to its minimizer.</summary>
    [TestMethod]
    public void LargeDirection_IsScaledConsistently()
    {
        var solver = new LBFGSB(p => 1e8 * Math.Pow(p[0] - 3, 2), 1, new[] { 0d },
            new[] { -1000d }, new[] { 1000d }, p => new[] { 2e8 * (p[0] - 3) }) { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(3d, solver.BestParameterSet.Values[0], 1e-12);
    }

    /// <summary>A wrong derivative cannot satisfy the line search for a constant objective.</summary>
    [TestMethod]
    public void InconsistentGradient_ReportsLineSearchFailure()
    {
        var solver = new LBFGSB(_ => 1d, 1, new[] { 0d }, new[] { -10d }, new[] { 10d },
            _ => new[] { 1d }) { ReportFailure = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.LineSearchFailed, solver.Status);
        Assert.AreEqual(0d, solver.BestParameterSet.Values[0]);
        Assert.IsNotNull(solver.Hessian);
    }

    /// <summary>A non-finite gradient at the start cannot be reported as converged.</summary>
    /// <param name="value">The invalid gradient component.</param>
    [TestMethod]
    [DataRow(double.NaN)]
    [DataRow(double.PositiveInfinity)]
    public void InvalidGradient_ReportsFailure(double value)
    {
        var solver = new LBFGSB(p => p[0] * p[0], 1, new[] { 1d }, new[] { -10d }, new[] { 10d },
            _ => new[] { value }) { ReportFailure = false, ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Failure, solver.Status);
    }

    /// <summary>An invalid initial objective must not produce a successful result.</summary>
    [TestMethod]
    public void InvalidInitialObjective_ReportsFailure()
    {
        var solver = new LBFGSB(_ => double.NaN, 1, new[] { 1d }, new[] { -10d }, new[] { 10d },
            _ => new[] { 0d }) { ReportFailure = false, ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Failure, solver.Status);
    }

    /// <summary>An invalid trial shortens the search rather than poisoning the iterate.</summary>
    [TestMethod]
    public void InvalidTrialObjective_Backtracks()
    {
        var solver = new LBFGSB(p => p[0] > 0 ? double.NaN : Math.Pow(p[0] + 0.1, 2), 1,
            new[] { -1d }, new[] { -10d }, new[] { 10d }, p => new[] { 2 * (p[0] + 0.1) })
        { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(-0.1, solver.BestParameterSet.Values[0], 1e-8);
    }

    /// <summary>Evaluation exhaustion remains distinct from line-search failure with either derivative source.</summary>
    /// <param name="analytic">Whether to use an analytical gradient.</param>
    [TestMethod]
    [DataRow(true)]
    [DataRow(false)]
    public void EvaluationBudget_IsPreserved(bool analytic)
    {
        var solver = new LBFGSB(p => -p[0], 1, new[] { 0d }, new[] { double.NegativeInfinity },
            new[] { double.PositiveInfinity }, analytic ? _ => new[] { -1d } : null)
        { MaxFunctionEvaluations = 10, ReportFailure = false, ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, solver.Status);
        Assert.AreEqual(10, solver.FunctionEvaluations);
    }

    /// <summary>Accepted iterations, including the last allowed step, are counted exactly.</summary>
    [TestMethod]
    public void IterationBudget_IsPreserved()
    {
        var solver = RosenbrockSolver();
        solver.MaxIterations = 10;
        solver.ReportFailure = false;
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.MaximumIterationsReached, solver.Status);
        Assert.AreEqual(10, solver.Iterations);
    }

    /// <summary>Rounded-equal function values must not prevent acceptance of a stationary point.</summary>
    /// <remarks>SciPy 1.17.1 L-BFGS-B with gtol = ftol = 1e-8 returns x = 1 after 1 iteration and 3 evaluations.</remarks>
    [TestMethod]
    public void RoundedObjective_StillChecksWolfeGradient()
    {
        var solver = new LBFGSB(p => 1e12 + Math.Pow(p[0] - 1, 2), 1,
            new[] { 1.0001 }, new[] { -10d }, new[] { 10d }, p => new[] { 2 * (p[0] - 1) })
        { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(1d, solver.BestParameterSet.Values[0], 1e-12);
        Assert.AreEqual(1, solver.Iterations);
        Assert.AreEqual(3, solver.FunctionEvaluations);
    }

    /// <summary>Safeguarded interpolation must resolve a very short step on a steep quadratic.</summary>
    [TestMethod]
    public void SteepQuadratic_ResolvesSmallWolfeStep()
    {
        var solver = new LBFGSB(p => 0.5e12 * p[0] * p[0], 1, new[] { 1e-8 },
            new[] { -10d }, new[] { 10d }, p => new[] { 1e12 * p[0] }) { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.IsLessThanOrEqualTo(solver.AbsoluteTolerance, Math.Abs(1e12 * solver.BestParameterSet.Values[0]));
    }

    /// <summary>Roundoff in a sum of squares must not reject its stationary minimum.</summary>
    /// <param name="offset">An additive objective constant, including a negative minimum.</param>
    /// <param name="maximize">Whether to maximize the negative quadratic.</param>
    [TestMethod]
    [DataRow(0d, false)]
    [DataRow(-1.25d, false)]
    [DataRow(0d, true)]
    public void RoundedUpMinimum_WithSuppliedGradient_Converges(double offset, bool maximize)
    {
        double sign = maximize ? -1d : 1d;
        Func<double[], double> objective = p => sign * (offset + 0.5 *
            ((p[0] - 1) * (p[0] - 1) + (p[0] + 1) * (p[0] + 1)));
        var solver = new LBFGSB(objective, 1, new[] { 6.8000000000000005e-9 },
            new[] { -10d }, new[] { 10d }, p => new[] { sign * 2 * p[0] })
        { ComputeHessian = false, ReportFailure = false };
        if (maximize) solver.Maximize(); else solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(0d, solver.BestParameterSet.Values[0], 1e-8);
        Assert.AreEqual(sign * objective(solver.BestParameterSet.Values), solver.BestParameterSet.Fitness);
    }

    /// <summary>A rounding-sized increase returns the start, which then passes the relative-reduction test.</summary>
    /// <remarks>
    /// The line search cannot find a decrease; when its interval collapses it returns the best step, zero, so the
    /// run stops at the start with no reduction. SciPy 1.17.1 L-BFGS-B gives the same result: x = 0 after 1
    /// iteration and 10 evaluations.
    /// </remarks>
    [TestMethod]
    public void RoundoffRise_WithoutStationarity_ReturnsStart()
    {
        var solver = new LBFGSB(p => p[0] == 0 ? 1d : 1d + 4 * Numerics.Tools.DoubleMachineEpsilon,
            1, new[] { 0d }, new[] { -10d }, new[] { 10d }, _ => new[] { 1e-7 })
        { ComputeHessian = false, ReportFailure = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(0d, solver.BestParameterSet.Values[0]);
        Assert.AreEqual(1d, solver.BestParameterSet.Fitness);
        Assert.AreEqual(1, solver.Iterations);
        Assert.AreEqual(10, solver.FunctionEvaluations);
    }

    /// <summary>A stationary trial with a resolvable objective increase must still be rejected.</summary>
    [TestMethod]
    public void StationaryUphillTrial_OutsideRoundoff_Backtracks()
    {
        var solver = new LBFGSB(p => 1 + 2 * p[0] * p[0] * p[0] - 2.5 * p[0] * p[0],
            1, new[] { 1d }, new[] { -2d }, new[] { 2d }, p => new[] { 6 * p[0] * p[0] - 5 * p[0] })
        { ComputeHessian = false, ReportFailure = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(5d / 6d, solver.BestParameterSet.Values[0], 1e-6);
        Assert.IsLessThan(0.5, solver.BestParameterSet.Fitness);
    }

    /// <summary>A non-finite derivative at a trial point steps back instead of ending the run.</summary>
    [TestMethod]
    public void InvalidTrialGradient_StepsBack()
    {
        var solver = new LBFGSB(p => 0.5 * ((p[0] - 1) * (p[0] - 1) + (p[0] + 1) * (p[0] + 1)),
            1, new[] { 0.5 }, new[] { -10d }, new[] { 10d },
            p => new[] { p[0] < 0 ? double.NaN : 2 * p[0] })
        { ComputeHessian = false, ReportFailure = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(0d, solver.BestParameterSet.Values[0], 1e-6);
    }

    /// <summary>Line search trials must not overwrite evaluation-budget termination.</summary>
    [TestMethod]
    public void RoundoffTrial_PreservesEvaluationBudget()
    {
        var solver = new LBFGSB(p => p[0] == 0 ? 1d : 1d + 4 * Numerics.Tools.DoubleMachineEpsilon,
            1, new[] { 0d }, new[] { -10d }, new[] { 10d }, _ => new[] { 1e-7 })
        { ComputeHessian = false, ReportFailure = false, MaxFunctionEvaluations = 10 };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.MaximumFunctionEvaluationsReached, solver.Status);
        Assert.AreEqual(10, solver.FunctionEvaluations);
        Assert.AreEqual(0, solver.Iterations);
    }

    /// <summary>Creates a Rosenbrock problem with its analytical derivative.</summary>
    /// <returns>The configured optimizer.</returns>
    private static LBFGSB RosenbrockSolver() => new LBFGSB(
        p => Math.Pow(1 - p[0], 2) + 100 * Math.Pow(p[1] - p[0] * p[0], 2), 2,
        new[] { -1.2, 1d }, new[] { -10d, -10d }, new[] { 10d, 10d },
        p => new[] { 2 * (p[0] - 1) - 400 * p[0] * (p[1] - p[0] * p[0]), 200 * (p[1] - p[0] * p[0]) })
        { ComputeHessian = false };
}
