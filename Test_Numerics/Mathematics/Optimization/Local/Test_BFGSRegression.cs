using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization;

/// <summary>Analytical regression tests for BFGS termination, bounds and evaluation contracts.</summary>
[TestClass]
public class Test_BFGSRegression
{
    /// <summary>A small objective change must not hide a nonzero gradient, including a warm start.</summary>
    /// <param name="x">Initial first coordinate.</param>
    /// <param name="y">Initial second coordinate.</param>
    /// <remarks>SciPy 1.16.2 BFGS with gtol=1e-8 independently returns the analytical minimum (0,0).</remarks>
    [TestMethod]
    [DataRow(1e-9, 1e-4)]
    [DataRow(1e-6, 1e-3)]
    public void ScaledQuadratic_RequiresStationarity(double x, double y)
    {
        var solver = new BFGS(p => 0.001 + 0.5 * (1e6 * p[0] * p[0] + p[1] * p[1]), 2,
            new[] { x, y }, new[] { -10d, -10d }, new[] { 10d, 10d },
            p => new[] { 1e6 * p[0], p[1] }) { ComputeHessian = false };
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
        var solver = new BFGS(p => p[0] * p[0], 1, new[] { 0d }, new[] { -10d }, new[] { 10d },
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
        var solver = new BFGS(p => -(p[0] - 2) * (p[0] - 2), 1,
            new[] { 0d }, new[] { -10d }, new[] { 10d },
            p => { buffer[0] = -2 * (p[0] - 2); return buffer; }) { ComputeHessian = false };
        solver.Maximize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(2d, solver.BestParameterSet.Values[0], 1e-8);
    }

    /// <summary>The standard Rosenbrock minimum must satisfy the requested gradient tolerance.</summary>
    /// <remarks>SciPy 1.16.2 BFGS at gtol=1e-8 reaches (1,1) from (-1.2,1).</remarks>
    [TestMethod]
    public void Rosenbrock_ReportsStationarySolution()
    {
        var solver = RosenbrockSolver();
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
        var solver = new BFGS(p => Math.Pow(p[0] - 2, 2) + Math.Pow(p[1] - 0.3, 2), 2,
            new[] { 0d, 0d }, new[] { 0d, 0d }, new[] { 1d, 1d },
            p => new[] { 2 * (p[0] - 2), 2 * (p[1] - 0.3) });
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(1d, solver.BestParameterSet.Values[0]);
        Assert.AreEqual(0.3, solver.BestParameterSet.Values[1], 1e-8);
        Assert.IsNotNull(solver.Hessian);
    }

    /// <summary>Strong-Wolfe slope calculations must use the direction after step scaling.</summary>
    [TestMethod]
    public void LargeDirection_IsScaledConsistently()
    {
        var solver = new BFGS(p => 1e8 * Math.Pow(p[0] - 3, 2), 1, new[] { 0d },
            new[] { -1000d }, new[] { 1000d }, p => new[] { 2e8 * (p[0] - 3) }) { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(3d, solver.BestParameterSet.Values[0], 1e-12);
    }

    /// <summary>A wrong derivative cannot satisfy strong Wolfe for a constant objective.</summary>
    [TestMethod]
    public void InconsistentGradient_ReportsLineSearchFailure()
    {
        var solver = new BFGS(_ => 1d, 1, new[] { 0d }, new[] { -10d }, new[] { 10d },
            _ => new[] { 1d }) { ReportFailure = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.LineSearchFailed, solver.Status);
        Assert.AreEqual(0, solver.Iterations);
        Assert.IsNotNull(solver.Hessian);
    }

    /// <summary>Non-finite gradients cannot be reported as converged.</summary>
    /// <param name="value">The invalid gradient component.</param>
    [TestMethod]
    [DataRow(double.NaN)]
    [DataRow(double.PositiveInfinity)]
    public void InvalidGradient_ReportsFailure(double value)
    {
        var solver = new BFGS(p => p[0] * p[0], 1, new[] { 1d }, new[] { -10d }, new[] { 10d },
            _ => new[] { value }) { ReportFailure = false, ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Failure, solver.Status);
    }

    /// <summary>An invalid initial objective must not produce a successful result.</summary>
    [TestMethod]
    public void InvalidInitialObjective_ReportsFailure()
    {
        var solver = new BFGS(_ => double.NaN, 1, new[] { 1d }, new[] { -10d }, new[] { 10d },
            _ => new[] { 0d }) { ReportFailure = false, ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Failure, solver.Status);
    }

    /// <summary>An invalid trial brackets the valid part of the search rather than poisoning the iterate.</summary>
    [TestMethod]
    public void InvalidTrialObjective_Backtracks()
    {
        var solver = new BFGS(p => p[0] > 0 ? double.NaN : Math.Pow(p[0] + 0.1, 2), 1,
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
        var solver = new BFGS(p => -p[0], 1, new[] { 0d }, new[] { double.NegativeInfinity },
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

    /// <summary>Rounded-equal function values must not prevent acceptance of a point satisfying both Wolfe conditions.</summary>
    [TestMethod]
    public void RoundedObjective_StillChecksWolfeGradient()
    {
        var solver = new BFGS(p => 1e12 + Math.Pow(p[0] - 1, 2), 1,
            new[] { 1.0001 }, new[] { -10d }, new[] { 10d }, p => new[] { 2 * (p[0] - 1) })
        { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.AreEqual(1d, solver.BestParameterSet.Values[0], 1e-12);
        Assert.AreEqual(1, solver.Iterations);
    }

    /// <summary>Safeguarded interpolation must resolve a very short step on a steep quadratic.</summary>
    [TestMethod]
    public void SteepQuadratic_ResolvesSmallWolfeStep()
    {
        var solver = new BFGS(p => 0.5e12 * p[0] * p[0], 1, new[] { 1e-8 },
            new[] { -10d }, new[] { 10d }, p => new[] { 1e12 * p[0] }) { ComputeHessian = false };
        solver.Minimize();
        Assert.AreEqual(OptimizationStatus.Success, solver.Status);
        Assert.IsLessThanOrEqualTo(solver.AbsoluteTolerance, Math.Abs(1e12 * solver.BestParameterSet.Values[0]));
    }

    /// <summary>Creates a Rosenbrock problem with its analytical derivative.</summary>
    /// <returns>The configured optimizer.</returns>
    private static BFGS RosenbrockSolver() => new BFGS(
        p => Math.Pow(1 - p[0], 2) + 100 * Math.Pow(p[1] - p[0] * p[0], 2), 2,
        new[] { -1.2, 1d }, new[] { -10d, -10d }, new[] { 10d, 10d },
        p => new[] { 2 * (p[0] - 1) - 400 * p[0] * (p[1] - p[0] * p[0]), 200 * (p[1] - p[0] * p[0]) })
        { ComputeHessian = false };
}
