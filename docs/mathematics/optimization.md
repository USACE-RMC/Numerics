# Optimization

[← Previous: Numerical Differentiation](differentiation.md) | [Back to Index](../index.md) | [Next: Root Finding →](root-finding.md)

Optimization is the process of finding the parameter set that minimizes (or maximizes) an objective function. The ***Numerics*** library provides a comprehensive suite of local and global optimization algorithms for both unconstrained and constrained problems. These methods are essential for parameter estimation, model calibration, machine learning, and engineering design optimization.

## Overview

| Method | Type | Best For | Requires Derivatives |
|--------|------|----------|---------------------|
| **Local Methods** | | | |
| BFGS | Local | Smooth, differentiable functions | No (numerical) |
| L-BFGS-B | Local | Smooth functions with simple bounds, many parameters | No (numerical) |
| Nelder-Mead | Local | General purpose, robust | No |
| Powell | Local | Smooth functions without derivatives | No |
| Golden Section | Local | 1D problems | No |
| Brent Search | Local | 1D problems, smooth functions | No |
| Gradient Descent | Local | Large-scale, smooth problems | No (numerical) |
| ADAM | Local | Machine learning, stochastic | No (numerical) |
| **Global Methods** | | | |
| Differential Evolution | Global | Multimodal, robust | No |
| Particle Swarm | Global | Multimodal, fast convergence | No |
| Shuffled Complex Evolution | Global | Hydrological calibration | No |
| Simulated Annealing | Global | Discrete/continuous, multimodal | No |
| Multi-Start | Global | Combines local search with random starts | No |
| MLSL | Global | Smooth multimodal functions | Requires local solver |
| Dynamically Dimensioned Search | Global | Expensive simulation calibration, fixed budget | No |
| CMA-ES | Global | Correlated, badly scaled continuous problems | No |
| **Constrained** | | | |
| Augmented Lagrangian | Constrained | Equality and inequality constraints | No |
| **Discrete and multi-objective** | | | |
| Dynamic Programming | Exact | Staged decisions, reservoir operation rules | No |
| Knapsack | Exact | Budget-constrained project selection | No |
| Pareto front utilities | Analysis | Non-dominated sorting, hypervolume | No |

## Problem Formulation

An optimization problem seeks to find:

```math
\min_{\mathbf{x}} f(\mathbf{x})
```

subject to:

```math
\mathbf{x}_L \leq \mathbf{x} \leq \mathbf{x}_U
```

where $f(\mathbf{x})$ is the objective function, $\mathbf{x}$ is the parameter vector, and $\mathbf{x}_L$, $\mathbf{x}_U$ are lower and upper bounds.

For maximization problems, minimize $-f(\mathbf{x})$ or use the `Maximize()` method.

## Common Interface

All optimizers in ***Numerics*** inherit from the `Optimizer` base class and share a common interface.

### Input Properties
- `MaxIterations`: Maximum iterations allowed (default = 10,000)
- `MaxFunctionEvaluations`: Maximum function evaluations (default = int.MaxValue)
- `AbsoluteTolerance`: Absolute convergence tolerance (default = 1E-8)
- `RelativeTolerance`: Relative convergence tolerance (default = 1E-8)
- `ReportFailure`: Rethrow an exception raised by the objective or the algorithm (default = true). Reaching `MaxIterations` or `MaxFunctionEvaluations` never throws from `Minimize()` or `Maximize()`; it is reported through `Status` only
- `RecordTraces`: Save optimization trace (default = true)
- `ComputeHessian`: Compute Hessian at solution (default = true)
- `CancellationToken`: Stop the run early; the run ends with `OptimizationStatus.Cancelled` and keeps its best point (default = none)
- `Progress`, `ProgressInterval`: Receive an `OptimizationProgress` snapshot every `ProgressInterval` evaluations and at the end (default = none, 100)

### Output Properties
- `BestParameterSet`: Optimal solution (Values, Fitness)
- `Iterations`: Number of iterations performed
- `FunctionEvaluations`: Number of function evaluations
- `Status`: Optimization status (Success, Failure, etc.)
- `ParameterSetTrace`: Read-only trace of the best-so-far parameter set and fitness at every function evaluation (entries recorded between improvements share one values array)
- `Hessian`: Numerically differentiated Hessian matrix (if computed)

### Methods
- `Minimize()`: Minimize the objective function
- `Maximize()`: Maximize the objective function (minimizes $-f(\mathbf{x})$)

After `Maximize()`, `BestParameterSet.Fitness` and every trace entry hold the value the solver minimized, which is $-f(\mathbf{x})$. The maximum of the objective is therefore `-BestParameterSet.Fitness`. The `Hessian` is that of $f$ itself and is not negated.

```cs
var brent = new BrentSearch(x => -(x - 1) * (x - 1) + 4, -5, 5);
brent.Maximize();
Console.WriteLine($"x = {brent.BestParameterSet.Values[0]:F4}, f(x) = {-brent.BestParameterSet.Fitness:F4}"); // x = 1.0000, f(x) = 4.0000
```

### Progress and Cancellation

Long calibrations can report progress and be stopped. The token is checked at every function evaluation, including the evaluations of a local search run inside a global one. A cancelled run does not throw: it keeps the best point found so far and ends with `OptimizationStatus.Cancelled`. A NaN objective value is ranked as the worst possible fitness (positive infinity), so a point where the objective is undefined is always replaced.

```cs
using System.Threading;
using Numerics.Mathematics.Optimization;

double Sphere(double[] x) => x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

using var stop = new CancellationTokenSource();
var de = new DifferentialEvolution(Sphere, 3, new[] { -5.0, -5, -5 }, new[] { 5.0, 5, 5 })
{
    CancellationToken = stop.Token,
    ProgressInterval = 500,
    Progress = new Progress<OptimizationProgress>(p =>
        Console.WriteLine($"{p.FunctionEvaluations} evaluations, best {p.BestFitness:E3}"))
};
de.Minimize();
Console.WriteLine($"Status: {de.Status}");
```

`Progress<T>` posts its callback to the thread pool or the calling synchronization context; implement `IProgress<OptimizationProgress>` directly for a synchronous callback. For `AugmentedLagrange`, also give the inner optimizer the token, because it evaluates the penalized objective itself.

## Local Optimization

Local optimization methods find the nearest local minimum from a starting point. They are fast and efficient for smooth, unimodal functions but may get trapped in local minima for multimodal problems.

### BFGS (Broyden-Fletcher-Goldfarb-Shanno)

BFGS is a quasi-Newton method that builds an approximation to the Hessian matrix using gradient information [[1]](#1). It's one of the most effective methods for smooth, unconstrained optimization.

```cs
using Numerics.Mathematics.Optimization;

// Rosenbrock function: f(x,y) = (1-x)² + 100(y-x²)²
double Rosenbrock(double[] x)
{
    return Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - x[0] * x[0], 2);
}

// Its analytic gradient
double[] RosenbrockGradient(double[] x)
{
    return new[] { -2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0] * x[0]), 200 * (x[1] - x[0] * x[0]) };
}

// Initial guess
var initial = new double[] { -1.2, 1.0 };
var lower = new double[] { -2, -2 };
var upper = new double[] { 2, 2 };

// Create and configure optimizer
var bfgs = new BFGS(Rosenbrock, 2, initial, lower, upper, RosenbrockGradient);
bfgs.RelativeTolerance = 1e-8;
bfgs.MaxIterations = 1000;

// Minimize
bfgs.Minimize();

// Results
Console.WriteLine($"Optimal solution: x = {bfgs.BestParameterSet.Values[0]:F6}, y = {bfgs.BestParameterSet.Values[1]:F6}");
Console.WriteLine($"Minimum value: {bfgs.BestParameterSet.Fitness:F10}");
Console.WriteLine($"Iterations: {bfgs.Iterations}");
Console.WriteLine($"Function evaluations: {bfgs.FunctionEvaluations}");
Console.WriteLine($"Status: {bfgs.Status}");
```

Supply the gradient when it is available. A finite-difference gradient has an error near $10^{-8}$, the default `AbsoluteTolerance`, so without the gradient this example ends with `LineSearchFailed` at the same point instead of `Success`; with it, the run takes 22 iterations and 25 evaluations.

The search starts from a unit inverse Hessian, so the first step is a steepest-descent step whose length depends on the scale of the objective. When the curvature is far from one, for example a log-likelihood summed over thousands of observations, set `ScaleInitialInverseHessian = true`. After the first accepted step, the identity is then multiplied by $y^\top s / y^\top y$, where $s$ is the step and $y$ the change in gradient, before the first update (Nocedal and Wright [[1]](#1), eq. 6.20). The same scaling is applied after each metric restart. The option is off by default, and default results are unchanged.

```cs
var scaled = new BFGS(Rosenbrock, 2, initial, lower, upper);
scaled.ScaleInitialInverseHessian = true;
scaled.Minimize();
Console.WriteLine($"Scaled: {scaled.FunctionEvaluations} evaluations, status {scaled.Status}");
```

**Advantages**: Fast convergence, memory efficient, works well for most smooth problems.

**Disadvantages**: Can fail on non-smooth functions, requires bounded search space.

### L-BFGS-B

`LBFGSB` is the limited-memory BFGS method for bound-constrained problems of Byrd, Lu, Nocedal and Zhu [[27]](#27), in the L-BFGS-B 3.0 form [[28]](#28) [[29]](#29). It keeps the last `MemorySize` (default 10) step and gradient-change pairs instead of a dense inverse Hessian. Each iteration finds the generalized Cauchy point along the projected steepest-descent path, which fixes the active bounds. It then minimizes the quadratic model over the free variables and takes a More-Thuente line search [[30]](#30). `MaxLineSearchSteps` (default 20) caps the function evaluations in one line search. The memory and the work per iteration grow linearly with the number of parameters, so it suits problems with hundreds or thousands of parameters.

Each bound may be infinite, so a parameter can be unbounded, bounded on one side, bounded on both sides, or fixed (`lower == upper`). A start outside the box is projected onto it. The run stops when the largest projected gradient component is at most `AbsoluteTolerance` or the relative reduction $(f_k - f_{k+1})/\max(|f_k|, |f_{k+1}|, 1)$ is at most `RelativeTolerance`. These are SciPy's `gtol` and `ftol`; set them to `1e-5` and `2.220446049250313e-9` to reproduce SciPy's defaults. A trial point with a non-finite objective or gradient makes the line search restart with half the step. When a line search fails with an empty memory, the status is `LineSearchFailed`; otherwise the memory is discarded and the iteration is retried.

```cs
using Numerics.Mathematics.Optimization;

// Rosenbrock in 10 dimensions with an analytic gradient; x[3] is held at or below 0.5
int n = 10;
double F(double[] x)
{
    double s = 0;
    for (int i = 0; i < n - 1; i++) s += 100 * Math.Pow(x[i + 1] - x[i] * x[i], 2) + Math.Pow(1 - x[i], 2);
    return s;
}
double[] G(double[] x)
{
    var g = new double[n];
    for (int i = 0; i < n - 1; i++)
    {
        double a = x[i + 1] - x[i] * x[i];
        g[i] += -400 * x[i] * a - 2 * (1 - x[i]);
        g[i + 1] += 200 * a;
    }
    return g;
}

var start = new double[n];
var lowerBounds = Enumerable.Repeat(double.NegativeInfinity, n).ToArray();
var upperBounds = Enumerable.Repeat(double.PositiveInfinity, n).ToArray();
upperBounds[3] = 0.5;

var lbfgsb = new LBFGSB(F, n, start, lowerBounds, upperBounds, G)
{
    AbsoluteTolerance = 1e-5,
    RelativeTolerance = 2.220446049250313e-9
};
lbfgsb.Minimize();

Console.WriteLine($"Status: {lbfgsb.Status}, iterations {lbfgsb.Iterations}, evaluations {lbfgsb.FunctionEvaluations}");
Console.WriteLine($"f = {lbfgsb.BestParameterSet.Fitness:F6}, x[3] = {lbfgsb.BestParameterSet.Values[3]:F6}");
```

**Agreement with SciPy.** On 18 problems the iteration and evaluation counts equal those of `scipy.optimize.minimize(method="L-BFGS-B")` 1.17.1 with the same settings, and the optimum agrees to $2 \times 10^{-9}$ or better. The problems are the BFGS test functions, Rosenbrock in 2, 5 and 10 dimensions, a quadratic with four active bounds, a problem whose solution has every variable at a bound, and one-sided and fixed bounds. The oracle scripts are in [oracles/lbfgsb](oracles/lbfgsb/).

**Speed.** Wall time against SciPy (BenchmarkDotNet mean; SciPy median of 21 runs, with NumPy objective and gradient functions that add Python overhead to each evaluation). Intel Core Ultra 9 185H, Windows 11, .NET 10 (BenchmarkDotNet 0.15.8, `Benchmarks_Numerics/Optimization/LBFGSBBenchmarks.cs`); Python 3.14.3, NumPy 2.4.2, SciPy 1.17.1 (`oracles/lbfgsb/lbfgsb_scipy_timing.py`). The iteration counts match except on the 100-dimensional problems, where dot-product rounding changes the path after a few hundred iterations; both reach the same optimum.

| Problem | Iterations, Numerics / SciPy | Numerics | SciPy |
|---|---|---|---|
| Rosenbrock, 10-D, unbounded | 71 / 71 | 0.36 ms | 6.7 ms |
| Rosenbrock, 100-D, unbounded | 513 / 517 | 6.9 ms | 48 ms |
| Rosenbrock, 100-D, on $[-2, 0.9]$ | 590 / 600 | 10.1 ms | 55 ms |
| Weighted quadratic, 1000-D, on $[-1, 1]$ | 14 / 14 | 1.25 ms | 7.5 ms |

**Advantages**: Handles bounds exactly, low memory, and stopping rules that match SciPy.

**Disadvantages**: Needs a smooth objective; a finite-difference gradient costs about $2D$ evaluations per iteration.

### Nelder-Mead Simplex

The Nelder-Mead method is a direct search algorithm that uses a simplex (geometric figure with n+1 vertices in n dimensions) to search the parameter space [[2]](#2). It's robust and doesn't require derivatives.

```cs
var nm = new NelderMead(Rosenbrock, 2, initial, lower, upper);
nm.RelativeTolerance = 1e-6;
nm.Minimize();

Console.WriteLine($"Optimal solution: [{nm.BestParameterSet.Values[0]:F6}, {nm.BestParameterSet.Values[1]:F6}]");
Console.WriteLine($"Minimum value: {nm.BestParameterSet.Fitness:F10}");
```

**Advantages**: Very robust, doesn't require derivatives, handles non-smooth functions.

**Disadvantages**: Slower convergence than gradient-based methods, can stagnate.

### Powell's Method

Powell's method is a conjugate direction algorithm that doesn't require derivatives [[1]](#1). It performs successive line searches along conjugate directions.

```cs
var powell = new Powell(Rosenbrock, 2, initial, lower, upper);
powell.RelativeTolerance = 1e-8;
powell.Minimize();

Console.WriteLine($"Solution: [{powell.BestParameterSet.Values[0]:F6}, {powell.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: No derivatives required, good for smooth functions.

**Disadvantages**: Can be slow in high dimensions, sensitive to scaling.

### Gradient Descent

Steepest descent with a fixed step size `Alpha` (there is no line search). The gradient is numerical unless one is supplied to the constructor:

```cs
var gd = new GradientDescent(Rosenbrock, 2, initial, lower, upper);
gd.Alpha = 0.001; // Step size (learning rate)
gd.Minimize();
```

### ADAM Optimizer

Adaptive Moment Estimation, popular in machine learning applications [[3]](#3):

```cs
var adam = new ADAM(Rosenbrock, 2, initial, lower, upper);
adam.Alpha = 0.001;
adam.Beta1 = 0.9;  // First moment decay
adam.Beta2 = 0.999; // Second moment decay
adam.Minimize();
```

**Advantages**: Adaptive learning rates, works well for noisy objectives.

**Best for**: Machine learning, stochastic optimization problems.

### One-Dimensional Methods

For single-parameter optimization, specialized 1D methods are more efficient:

#### Golden Section Search

Uses the golden ratio to bracket the minimum:

```cs
Func<double, double> f1d = x => Math.Pow(x - 2, 2) + 3;

var golden = new GoldenSection(f1d, 0.0, 5.0);
golden.Minimize();

Console.WriteLine($"Minimum at x = {golden.BestParameterSet.Values[0]:F6}");
```

#### Brent Search

Combines golden section search with parabolic interpolation for faster convergence. The search is a port of SciPy's bounded minimizer, `scipy.optimize.minimize_scalar(method='bounded')`. Its step tolerance is `RelativeTolerance * |x| + AbsoluteTolerance / 3`; the constructor sets `AbsoluteTolerance` to about $3.3 \times 10^{-19}$, so the default is $10^{-8}|x| + 1.1 \times 10^{-19}$. With SciPy's settings the result is bitwise SciPy's:

```cs
var brent = new BrentSearch(f1d, 0.0, 5.0);
brent.Minimize();

// SciPy's defaults: xatol = 1e-5, maxiter = 500 evaluations, and its sqrt(2.2e-16) relative term.
var scipyLike = new BrentSearch(f1d, 0.0, 5.0)
{
    RelativeTolerance = Math.Sqrt(2.2E-16),
    AbsoluteTolerance = 1E-5,
    MaxFunctionEvaluations = 500
};
scipyLike.Minimize();
```

## Global Optimization

Global optimization methods are designed to find the global minimum across the entire search space, avoiding local minima. They typically require more function evaluations but are more robust for multimodal problems.

### Polishing

A global search finds the right basin but converges slowly inside it. Differential Evolution, Particle Swarm, SCE-UA, Simulated Annealing, DDS, and CMA-ES can refine their best point with a bounded local search when the global search ends:

- `Polish`: Run the local search (default = false; true for Simulated Annealing)
- `PolishMethod`: `LocalMethod.NelderMead` (default), `LocalMethod.BFGS`, or `LocalMethod.Powell`

The polish starts at the best point, uses the optimizer's tolerances and remaining evaluation budget, and its evaluations count toward `FunctionEvaluations`. It runs after convergence or after `MaxIterations`, but not after the evaluation budget is spent or a cancellation, or when fewer than 10 evaluations remain in the budget. The reported status is the global search's. Multi-Start and MLSL always polish with their own `LocalMethod` (`Polish` = true by default).

```cs
double Rosen(double[] x) => Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - x[0] * x[0], 2);
var polished = new DifferentialEvolution(Rosen, 2, new[] { -5.0, -5 }, new[] { 5.0, 5 }) { Polish = true };
polished.Minimize();
Console.WriteLine($"Polished minimum: {polished.BestParameterSet.Fitness:E2}");
```

### Differential Evolution

Differential Evolution (DE) is a population-based evolutionary algorithm that's very robust for continuous optimization [[4]](#4). For each member $\mathbf{x}_i$ of the population, DE creates a **trial vector** $\mathbf{u}$ using mutation and crossover:

**Mutation** (DE/rand/1): Three distinct population members $\mathbf{x}_{r_0}$, $\mathbf{x}_{r_1}$, $\mathbf{x}_{r_2}$ are randomly selected, and a **mutant vector** is formed:

```math
\mathbf{v} = \mathbf{x}_{r_0} + G \cdot (\mathbf{x}_{r_1} - \mathbf{x}_{r_2})
```

where $G$ is a scale factor. In the ***Numerics*** implementation, $G$ is **dithered** with 90% probability: $G = 0.5 + \text{rand}() \times 0.5$, which helps avoid stagnation.

**Crossover** (binomial): The trial vector $\mathbf{u}$ is assembled component-by-component:

```math
u_j = \begin{cases} v_j & \text{if } \text{rand}() \leq CR \text{ or } j = j_{\text{rand}} \\ x_{i,j} & \text{otherwise} \end{cases}
```

where $CR$ is the crossover probability and $j_{\text{rand}}$ is a randomly chosen index that ensures at least one component comes from the mutant vector.

**Selection** (greedy): The trial vector replaces the target only if it has equal or better fitness:

```math
\mathbf{x}_i^{(g+1)} = \begin{cases} \mathbf{u} & \text{if } f(\mathbf{u}) \leq f(\mathbf{x}_i^{(g)}) \\ \mathbf{x}_i^{(g)} & \text{otherwise} \end{cases}
```

The algorithm converges when the standard deviation of fitness values across the population falls below the tolerance.

```cs
using Numerics.Mathematics.Optimization;

// Rastrigin function (highly multimodal)
double Rastrigin(double[] x)
{
    double sum = 10 * x.Length;
    for (int i = 0; i < x.Length; i++)
    {
        sum += x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]);
    }
    return sum;
}

var lower = new double[] { -5.12, -5.12 };
var upper = new double[] { 5.12, 5.12 };

var de = new DifferentialEvolution(Rastrigin, 2, lower, upper);
de.PopulationSize = 50;
de.CrossoverProbability = 0.9;
de.Mutation = 0.75;
de.MaxIterations = 1000;
de.Minimize();

Console.WriteLine($"Global minimum: [{de.BestParameterSet.Values[0]:F6}, {de.BestParameterSet.Values[1]:F6}]");
Console.WriteLine($"Function value: {de.BestParameterSet.Fitness:F10}");
```

**Advantages**: Very robust, handles discontinuities, good for difficult problems.

**Parameters**:
- `PopulationSize`: Number of candidate solutions (default = 10 × dimensions)
- `CrossoverProbability`: Probability of crossover (default = 0.9)
- `Mutation`: Scaling factor for mutation (default = 0.75)
- `InitialPopulation`: Starting members that replace the first Latin hypercube members (default = none)
- `IntegerParameters`: Flags for integer-valued parameters (default = none)
- `UseParallel`, `MaxDegreeOfParallelism`: Evaluate each generation in parallel (default = off, no limit)

**Warm start and integer parameters**: `InitialPopulation` seeds the search with known points, such as a previous calibration or a baseline schedule; the other members and the random stream are unchanged. `IntegerParameters` marks parameters that must take integer values, such as start years in a project schedule: they are rounded to the nearest integer (ties to even) within $[\lceil x_L \rceil, \lfloor x_U \rfloor]$ before every evaluation, as SciPy's `integrality` option does. Polishing is skipped when there are integer parameters, and the Hessian is not meaningful for them.

```cs
// Two integer start years and one continuous fraction.
double Schedule(double[] x) => Math.Pow(x[0] - 3.4, 2) + Math.Pow(x[1] - 7.6, 2) + Math.Pow(x[2] - 0.25, 2);
var plan = new DifferentialEvolution(Schedule, 3, new[] { 0.0, 0, 0 }, new[] { 10.0, 10, 1 })
{
    IntegerParameters = new[] { true, true, false },
    InitialPopulation = new List<double[]> { new[] { 3.0, 8.0, 0.5 } },
    ComputeHessian = false
};
plan.Minimize();
Console.WriteLine($"Start years {plan.BestParameterSet.Values[0]}, {plan.BestParameterSet.Values[1]}; fraction {plan.BestParameterSet.Values[2]:F4}"); // 3, 8; 0.2500
```

**Parallel evaluation**: With `UseParallel = true`, each generation's trial vectors are built from the population at the start of the generation, evaluated concurrently, and selected in index order (deferred updating, as SciPy's `updating='deferred'`). The random draws stay on one thread, so the result is reproducible and does not depend on the number of threads. It is a different search from the default immediate updating, where an accepted trial can serve as a donor later in the same generation. The objective must be thread safe. For a simulation objective costing about 1 ms per run, the speedup measured on a 22-core machine was about 9 times. Particle Swarm has the same option, with synchronous moves toward the previous iteration's global best.

### Particle Swarm Optimization

PSO simulates social behavior of bird flocking or fish schooling [[5]](#5). Each particle $i$ has a position $\mathbf{x}_i$ and velocity $\mathbf{v}_i$ that are updated at each iteration using:

```math
v_{i,j}^{(t+1)} = w \cdot v_{i,j}^{(t)} + c_1 \cdot r_1 \cdot (p_{i,j} - x_{i,j}^{(t)}) + c_2 \cdot r_2 \cdot (g_j - x_{i,j}^{(t)})
```
```math
x_{i,j}^{(t+1)} = x_{i,j}^{(t)} + v_{i,j}^{(t+1)}
```

where:
- $w$ is the **inertia weight**, which decreases linearly from $w_{\max} = 0.9$ to $w_{\min} = 0.4$ over the optimization, balancing exploration (high $w$) and exploitation (low $w$)
- $c_1 = c_2 = 2.05$ are the cognitive and social acceleration coefficients
- $r_1, r_2 \sim U(0,1)$ are independent random numbers (per component, per particle)
- $\mathbf{p}_i$ is particle $i$'s personal best position (best position it has visited)
- $\mathbf{g}$ is the global best position (best position any particle has visited)

The three terms represent **momentum** (continue in the current direction), **cognitive pull** (return toward personal best), and **social pull** (move toward the swarm's best).

```cs
var pso = new ParticleSwarm(Rastrigin, 2, lower, upper);
pso.PopulationSize = 40;
pso.Minimize();

Console.WriteLine($"Solution: [{pso.BestParameterSet.Values[0]:F6}, {pso.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Fast convergence, simple concept, works well for continuous problems.

**Parameters**:
- `PopulationSize`: Number of particles (default = 30)
- `UseParallel`, `MaxDegreeOfParallelism`: Synchronous swarm with parallel evaluation (default = off)
- `Polish`, `PolishMethod`: Refine the best point at the end (default = off)

### Shuffled Complex Evolution (SCE-UA)

SCE-UA was specifically developed for calibrating hydrological models [[6]](#6). The algorithm partitions the population into $p$ **complexes**, evolves each complex independently using a **Competitive Complex Evolution** (CCE) strategy, then shuffles members between complexes to share information.

The CCE step within each complex proceeds as follows:

1. **Sub-complex selection**: Select $q = D+1$ distinct points from the complex using a trapezoidal probability distribution that favors better-ranked points: $P(i) = \frac{2(N+1-i)}{N(N+1)}$

2. **Reflection**: Compute the centroid $\mathbf{g}$ of all sub-complex points except the worst point $\mathbf{x}_w$, then reflect:

```math
\mathbf{r} = 2\mathbf{g} - \mathbf{x}_w
```

3. **Contraction**: If the reflected point is infeasible or worse than $\mathbf{x}_w$, try contraction:

```math
\mathbf{c} = \frac{\mathbf{g} + \mathbf{x}_w}{2}
```

4. **Mutation**: If contraction also fails, generate a random point within the smallest bounding box of the complex. It replaces the worst point whatever its value.

The evolved sub-complex is written back to the positions it was drawn from, as in Duan et al. (1994) [[22]](#22).

The shuffling step redistributes points across complexes, preventing any single complex from stagnating. This combination of local competitive evolution within complexes and global information sharing between complexes makes SCE-UA particularly effective for the rugged, multimodal objective functions typical of hydrological model calibration.

```cs
var sce = new ShuffledComplexEvolution(Rastrigin, 2, lower, upper);
sce.Complexes = 5;
sce.MaxIterations = 1000;
sce.Minimize();

Console.WriteLine($"Calibrated parameters: [{sce.BestParameterSet.Values[0]:F6}, {sce.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Excellent for hydrological model calibration, balances exploration and exploitation.

**Best for**: Calibrating watershed models, water resources applications.

### Dynamically Dimensioned Search (DDS)

DDS is a greedy, single-solution global search developed for calibrating watershed models when each model run is expensive [[14]](#14). Instead of stopping on a convergence test, DDS is given a budget of $m$ function evaluations and scales its search to that budget: it is global at the start and becomes local as the budget is spent, by perturbing fewer and fewer parameters of the best solution found so far. The ***Numerics*** implementation follows Figure 1 of Tolson and Shoemaker (2007) exactly.

**Initialization**: The search starts from a supplied initial point, or from the best of $M = \max(5, \lfloor 0.005\,m \rfloor)$ uniform random points. These evaluations count toward $m$.

**Neighborhood selection**: With $i$ the number of evaluations made so far, each of the $D$ parameters is selected for perturbation with probability

```math
P(i) = 1 - \frac{\ln i}{\ln m}
```

If no parameter is selected, one is chosen at random, so at least one parameter always changes. $P$ falls from 1 at the start to nearly 0 at the end of the budget.

**Perturbation**: Each selected parameter of the best solution is perturbed with a normal step whose standard deviation is a fixed fraction $r$ of the parameter range:

```math
x_j^{\text{new}} = x_j^{\text{best}} + \sigma_j\, Z, \qquad \sigma_j = r\,(x_j^{\max} - x_j^{\min}), \qquad Z \sim N(0, 1)
```

The bounds are reflecting. A value below $x_j^{\min}$ becomes $x_j^{\min} + (x_j^{\min} - x_j^{\text{new}})$, and if that passes $x_j^{\max}$ it is set to $x_j^{\min}$; the upper bound is handled the same way.

**Selection** (greedy): The candidate replaces the best solution when $f(\mathbf{x}^{\text{new}}) \leq f(\mathbf{x}^{\text{best}})$. The search stops when the $m$ evaluations are spent.

```cs
using Numerics.Mathematics.Optimization;

// Rastrigin function (highly multimodal)
double Rastrigin(double[] x)
{
    double sum = 10 * x.Length;
    for (int i = 0; i < x.Length; i++)
    {
        sum += x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]);
    }
    return sum;
}

var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };

var dds = new DynamicallyDimensionedSearch(Rastrigin, 5, lower, upper);
dds.EvaluationBudget = 5000;   // m: the search is scaled to this budget
dds.PerturbationSize = 0.2;    // r: perturbation standard deviation as a fraction of the range
dds.Minimize();

Console.WriteLine($"Status: {dds.Status}, evaluations: {dds.FunctionEvaluations}");
Console.WriteLine($"Best value: {dds.BestParameterSet.Fitness:F4}");

// Warm start from a previous calibration: one initial evaluation instead of a random sample.
var start = new double[] { 1.0, -1.0, 0.5, 0.0, 2.0 };
var warm = new DynamicallyDimensionedSearch(Rastrigin, 5, start, lower, upper) { EvaluationBudget = 1000 };
warm.Minimize();
Console.WriteLine($"Warm start best value: {warm.BestParameterSet.Fitness:F4}");
```

Spending the budget is the normal end of a DDS run, so it reports `OptimizationStatus.Success` and, when `ComputeHessian` is true, computes the Hessian at the best point. `MaxFunctionEvaluations` and `MaxIterations` still act as hard caps: if either is reached before the budget, the run stops early with that status. `Iterations` counts the DDS iterations, which is the budget less the initial evaluations.

**Advantages**: Needs no tuning, makes good use of a small evaluation budget, and handles many parameters well because most iterations change only a few of them.

**Parameters**:
- `EvaluationBudget`: The planned number of evaluations $m$, including the initial ones (default = 1,000)
- `PerturbationSize`: The neighborhood size $r$, in $(0, 1]$ (default = 0.2, as recommended by the authors)
- `PRNGSeed`: The random seed (default = 12345)
- `Polish`, `PolishMethod`: Refine the best point after the budget (default = off)

**Best for**: Calibrating watershed and other simulation models with roughly 5 to 30 parameters, where each run is expensive and you can only afford a few hundred to a few thousand runs. Set `EvaluationBudget` to the number of model runs you can afford, not to a large safety value, because the search schedule depends on it. DDS is designed to find good solutions within the budget, not to converge precisely. When you need precise parameter values, set `Polish = true` (see [Polishing](#polishing)) or refine the DDS result with a local method such as `NelderMead` or `BFGS`. With a known good starting point, the authors suggest reducing $r$ to about 0.1. For narrow curved valleys such as Rosenbrock's, or for cheap objectives where tens of thousands of evaluations are affordable, `ShuffledComplexEvolution` usually converges more precisely.


### Covariance Matrix Adaptation Evolution Strategy (CMA-ES)

CMA-ES samples a population from a multivariate normal distribution and adapts the mean, the step size and the full covariance matrix from the ranked samples [[23]](#23), [[24]](#24). Each generation draws $\lambda$ points

```math
\mathbf{x}_k = \mathbf{m} + \sigma \mathbf{B}\mathbf{D}\mathbf{z}_k, \qquad \mathbf{z}_k \sim \mathcal{N}(\mathbf{0}, \mathbf{I}), \qquad \mathbf{C} = \mathbf{B}\mathbf{D}^2\mathbf{B}^\top
```

and ranks them by fitness. With $\mathbf{y}_{i:\lambda} = (\mathbf{x}_{i:\lambda} - \mathbf{m})/\sigma$ the $i$-th best step and $\langle \mathbf{y} \rangle_w = \sum_{i=1}^{\mu} w_i \mathbf{y}_{i:\lambda}$, the mean moves to the weighted average of the $\mu = \lfloor \lambda/2 \rfloor$ best points, $\mathbf{m} \leftarrow \mathbf{m} + \sigma \langle \mathbf{y} \rangle_w$.

**Step size** (cumulative step-size adaptation): the conjugate evolution path $\mathbf{p}_\sigma$ accumulates the normalized steps, and $\sigma$ grows when consecutive steps point the same way and shrinks when they cancel:

```math
\mathbf{p}_\sigma \leftarrow (1 - c_\sigma)\mathbf{p}_\sigma + \sqrt{c_\sigma(2 - c_\sigma)\mu_{\text{eff}}}\, \mathbf{C}^{-1/2} \langle \mathbf{y} \rangle_w, \qquad
\sigma \leftarrow \sigma \exp\!\left(\frac{c_\sigma}{d_\sigma}\left(\frac{\lVert \mathbf{p}_\sigma \rVert}{\mathrm{E}\lVert \mathcal{N}(\mathbf{0}, \mathbf{I}) \rVert} - 1\right)\right)
```

**Covariance matrix**: the rank-one update uses the evolution path $\mathbf{p}_c$, and the rank-$\mu$ update uses all $\lambda$ ranked steps, with negative weights for the worst half (active CMA):

```math
\mathbf{p}_c \leftarrow (1 - c_c)\mathbf{p}_c + h_\sigma \sqrt{c_c(2 - c_c)\mu_{\text{eff}}}\, \langle \mathbf{y} \rangle_w
```
```math
\mathbf{C} \leftarrow \left(1 + c_1\delta(h_\sigma) - c_1 - c_\mu \textstyle\sum_j w_j\right)\mathbf{C} + c_1 \mathbf{p}_c \mathbf{p}_c^\top + c_\mu \sum_{i=1}^{\lambda} w_i^\circ\, \mathbf{y}_{i:\lambda}\mathbf{y}_{i:\lambda}^\top
```

where $h_\sigma$ stalls $\mathbf{p}_c$ while $\lVert \mathbf{p}_\sigma \rVert$ is large, and $w_i^\circ = w_i\, n / \lVert \mathbf{C}^{-1/2}\mathbf{y}_{i:\lambda} \rVert^2$ for the negative weights. All learning rates and weights are the defaults of Table 1 in Hansen's tutorial [[23]](#23), with $\lambda = 4 + \lfloor 3 \ln n \rfloor$. The eigendecomposition of $\mathbf{C}$ is updated lazily, after more than $\lambda / (c_1 + c_\mu) / n / 10$ evaluations. The updates depend on the ranks of the fitness values only, and the covariance learns the local shape of the function, so away from the bounds the method is invariant to rotations of the search space and to monotone transformations of the objective.

**Bounds**: the search runs in the unit cube, where every parameter range has length one. The objective is evaluated at $T(\mathbf{y})$, where $T$ is pycma's periodic linear-quadratic `BoundTransform` [[26]](#26): it is the identity on $[0.05, 0.95]$, quadratic within 0.05 of each bound, and maps every real number into $[0, 1]$. Every evaluated point is inside the bounds, and the strategy itself is unchanged. As in pycma, each coordinate standard deviation $\sigma\sqrt{C_{ii}}$ is capped at one third of the bound range.

**Restarts (IPOP)**: a run stops on the tutorial's termination criteria (TolFun, TolX, EqualFunValues, NoEffectAxis, NoEffectCoord, ConditionCov, Stagnation, TolXUp). The solver then restarts from a uniformly random mean with the population size doubled [[25]](#25), up to `Restarts` times. Larger populations smooth out the local structure of multimodal functions such as Rastrigin. TolFun uses the common tolerances: a run converges when the best values of the last $10 + \lceil 30n/\lambda \rceil$ generations and all values of the current generation lie within `AbsoluteTolerance + RelativeTolerance` $\times |f_{\text{best}}|$.

```cs
using Numerics.Mathematics.Optimization;

// Rastrigin function shifted so that the minimum f = 0 is at x = (1, ..., 1)
double ShiftedRastrigin(double[] x)
{
    double sum = 10 * x.Length;
    for (int i = 0; i < x.Length; i++)
    {
        double z = x[i] - 1;
        sum += z * z - 10 * Math.Cos(2 * Math.PI * z);
    }
    return sum;
}

int n = 5;
var lower = Enumerable.Repeat(-5.12, n).ToArray();
var upper = Enumerable.Repeat(5.12, n).ToArray();

var cma = new CMAES(ShiftedRastrigin, n, lower, upper);
cma.Restarts = 9;              // IPOP restarts with doubled population (default)
cma.MaxFunctionEvaluations = 200000;
cma.Minimize();

Console.WriteLine($"Status: {cma.Status}, restarts: {cma.RestartsPerformed}, evaluations: {cma.FunctionEvaluations}");
Console.WriteLine($"Solution: [{string.Join(", ", cma.BestParameterSet.Values.Select(v => v.ToString("F6")))}]");
Console.WriteLine($"Function value: {cma.BestParameterSet.Fitness:E3}");

// A single run from a given starting point, for a unimodal problem
double Rosenbrock(double[] x) => 100 * Math.Pow(x[1] - x[0] * x[0], 2) + Math.Pow(1 - x[0], 2);
var local = new CMAES(Rosenbrock, 2, new[] { -1.2, 1.0 }, new[] { -5.0, -5.0 }, new[] { 5.0, 5.0 });
local.Restarts = 0;
local.InitialStepSize = 0.1;   // in units of the bound range
local.Minimize();
Console.WriteLine($"Rosenbrock: [{local.BestParameterSet.Values[0]:F6}, {local.BestParameterSet.Values[1]:F6}] after {local.FunctionEvaluations} evaluations");
```

The IPOP run finds the global minimum of the 5-D Rastrigin function. Once a run has converged, reaching `MaxFunctionEvaluations` or `MaxIterations` during a later restart ends the optimization normally with `Success` and the best point found, so the evaluation limit acts as the restart budget. Without a limit, all `Restarts` runs are performed; with the default of 9 restarts the last population is 512 times the first, and a 10-D problem then uses roughly 0.5 to 1 million evaluations in total.

**Parameters**:
- `PopulationSize`: $\lambda$ of the first run (default = $4 + \lfloor 3 \ln n \rfloor$)
- `InitialStepSize`: $\sigma_0$ in units of the bound range (default = 0.3)
- `Restarts`: maximum number of IPOP restarts (default = 9; 0 for a single run)
- `MaxFunctionEvaluations`: evaluation budget shared by all runs (default = 10,000 × $n$ for CMA-ES). Once a run has converged, reaching the budget during a later restart ends the search with `Success`
- `InitialValues`: optional initial mean of the first run (default = center of the bounds), set through the constructor overload
- `Polish`, `PolishMethod`: Refine the best point after the restarts (default = off)
- `PRNGSeed`: random seed (default = 12345)

**When to prefer CMA-ES over Differential Evolution**: CMA-ES learns correlations and scaling between parameters, so it is much faster on smooth, ill-conditioned or non-separable problems. In 10 dimensions it reached $f < 10^{-8}$ about 30 times faster than DE on the sphere and 38 times faster on a rotated ellipsoid with condition number $10^6$; on the 10-D Rosenbrock function it needed about 5,500 evaluations against 116,000, and on the 20-D Rosenbrock function DE did not reach the target within 400,000 evaluations. With IPOP restarts it also solved the 10-D Rastrigin function in every one of 30 runs, with a median of about 91,000 evaluations against 568,000 for DE (26 of 30 runs). DE remains a good choice when the objective is discontinuous or noisy at the scale of the search, when parameters are nearly independent, or when evaluations can be spread over a large, fixed population. Use `Restarts = 0` when the problem is known to be unimodal; each restart costs at least as much as the run before it.

### Simulated Annealing

SA mimics the physical process of annealing in metallurgy [[7]](#7). It accepts uphill moves with decreasing probability, allowing escape from local minima.

```cs
var sa = new SimulatedAnnealing(Rastrigin, 2, lower, upper);
sa.InitialTemperature = 10.0;
sa.CoolingRate = 0.95;
sa.MaxIterations = 10000;
sa.Minimize();

Console.WriteLine($"Solution: [{sa.BestParameterSet.Values[0]:F6}, {sa.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Can escape local minima, works for discrete and continuous problems.

**Stopping and polishing**: each temperature reduction makes `TemperatureCycles` × `UpdateCycles` sweeps over the $D$ coordinates (up to $40D$ evaluations with the defaults). The search stops with `OptimizationStatus.Success` when the best value has not changed, within the tolerances, over the last `ToleranceSteps` (default = 20) reductions, and then polishes the best point with Nelder-Mead (`Polish` = true by default). The stall test usually ends a search after a few thousand evaluations, far fewer than the full `MaxIterations` schedule.

**Parameters**:
- `InitialTemperature`: Starting temperature (default = 10)
- `CoolingRate`: Temperature reduction factor (default = 0.95)

### Multi-Start Optimization

Combines local search with multiple random starting points. The local search is `LocalMethod.BFGS` (default), `LocalMethod.NelderMead`, `LocalMethod.Powell` or `LocalMethod.LBFGSB`, which holds the search exactly on the box and also polishes the final result; MLSL accepts the same four.

```cs
var initial = new double[] { 0.0, 0.0 };
var ms = new MultiStart(Rastrigin, 2, initial, lower, upper, LocalMethod.BFGS);
ms.MaxIterations = 20; // Number of random starts
ms.Minimize();

Console.WriteLine($"Best solution: [{ms.BestParameterSet.Values[0]:F6}, {ms.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Simple to implement, leverages fast local search.

**Disadvantages**: May waste evaluations in same basin of attraction.

Completing all the starts is the normal end, so the run reports the status of the final polish (`Success` in a normal run) and computes the Hessian.

### MLSL (Multi-Level Single Linkage)

Clustering-based global optimization that avoids redundant local searches:

```cs
var initial = new double[] { 0.0, 0.0 };
var mlsl = new MLSL(Rastrigin, 2, initial, lower, upper, LocalMethod.BFGS);
mlsl.Minimize();
```

**Advantages**: More efficient than multi-start, avoids redundant searches.

A sample point starts a local search only if no better sample point or local minimum lies within the critical distance of Rinnooy Kan and Timmer [[21]](#21), measured in the box scaled to the unit cube:

```math
r_k = \pi^{-1/2} \left( \Gamma\left(1 + \frac{D}{2}\right) \sigma \frac{\ln kN}{kN} \right)^{1/D}
```

where $k$ is the iteration, $N$ = `SampleSize`, and $\sigma$ = `Sigma` (default = 2).

## Constrained Optimization

### Augmented Lagrangian

The Augmented Lagrangian method handles equality and inequality constraints by adding penalty terms to the objective function [[8]](#8).

```cs
using Numerics.Mathematics.Optimization;

// Objective: minimize x² + y²
double Objective(double[] x)
{
    return x[0] * x[0] + x[1] * x[1];
}

// Constraint: x + y >= 1
var constraint = new Constraint(
    x => x[0] + x[1], 2, 1.0,  // g(x) >= value form
    ConstraintType.GreaterThanOrEqualTo
);

var lower = new double[] { -5, -5 };
var upper = new double[] { 5, 5 };
var initial = new double[] { 0, 0 };

// Create the inner optimizer
var bfgs = new BFGS(Objective, 2, initial, lower, upper);

// Create the Augmented Lagrange optimizer with constraints
var al = new AugmentedLagrange(Objective, bfgs, new IConstraint[] { constraint });
al.MaxIterations = 100;
al.Minimize();

Console.WriteLine($"Optimal solution: [{al.BestParameterSet.Values[0]:F6}, {al.BestParameterSet.Values[1]:F6}]");
Console.WriteLine($"Constraint satisfied: {al.BestParameterSet.Values[0] + al.BestParameterSet.Values[1]:F6} >= 1");
```

**Constraint Types** (each compares the constraint function $g(\mathbf{x})$ with the constraint's `Value` $b$):
- `ConstraintType.EqualTo`: Equality constraint $g(\mathbf{x}) = b$
- `ConstraintType.LesserThanOrEqualTo`: Inequality constraint $g(\mathbf{x}) \leq b$
- `ConstraintType.GreaterThanOrEqualTo`: Inequality constraint $g(\mathbf{x}) \geq b$

A constraint counts as satisfied when it is violated by no more than its `Tolerance` (default = 1E-8). Each inequality adds the Powell-Hestenes-Rockafellar term $\frac{\rho}{2}\max(0, c + \mu/\rho)^2$ [[8]](#8), which is continuous at the boundary, and the result is the best feasible iterate, or the least-violating one when none is feasible.

**Example: Minimize subject to multiple constraints**:

```cs
// Minimize f(x,y) = (x-3)² + (y-2)²
// Subject to: x + y <= 5
//             x >= 1
//             y >= 1

double ObjectiveFunc(double[] x)
{
    return Math.Pow(x[0] - 3, 2) + Math.Pow(x[1] - 2, 2);
}

var c1 = new Constraint(x => x[0] + x[1], 2, 5.0, ConstraintType.LesserThanOrEqualTo);
var c2 = new Constraint(x => x[0], 2, 1.0, ConstraintType.GreaterThanOrEqualTo);
var c3 = new Constraint(x => x[1], 2, 1.0, ConstraintType.GreaterThanOrEqualTo);

var innerOptimizer = new BFGS(ObjectiveFunc, 2, new[] { 2.0, 2.0 },
                              new[] { 0.0, 0.0 }, new[] { 10.0, 10.0 });
var constrained = new AugmentedLagrange(ObjectiveFunc, innerOptimizer,
                                        new IConstraint[] { c1, c2, c3 });
constrained.Minimize();

Console.WriteLine($"Constrained optimum: [{constrained.BestParameterSet.Values[0]:F4}, " +
                  $"{constrained.BestParameterSet.Values[1]:F4}]");
```

## Pareto Fronts

`ParetoFront` is a static class of pure functions for multi-objective results: dominance, non-dominated sorting, crowding distance, hypervolume, and the generational distance indicators. These are the building blocks of multi-objective optimizers such as NSGA-II, and they also summarize any set of trade-off alternatives, for example a dam safety portfolio scored on cost and risk.

All objectives are **minimized**; negate an objective to maximize it. A set of $n$ points with $m$ objectives is passed as `double[][]` (one row per point) or `double[,]`. NaN is rejected with an `ArgumentException`. Dominance and sorting accept infinite values (for example, a penalized infeasible point); crowding distance, hypervolume and the distance indicators require finite values.

### Dominance and Non-Dominated Sorting

Point $\mathbf{a}$ **dominates** $\mathbf{b}$, written $\mathbf{a} \prec \mathbf{b}$, when it is no worse in every objective and strictly better in at least one:

```math
\mathbf{a} \prec \mathbf{b} \iff a_k \le b_k \;\; \forall k \quad \text{and} \quad a_j < b_j \;\; \text{for some } j
```

`ParetoFront.Dominates(a, b)` tests the relation for two points:

```cs
bool d1 = ParetoFront.Dominates(new[] { 2.0, 3.0 }, new[] { 3.0, 4.0 });   // True
bool d2 = ParetoFront.Dominates(new[] { 1.0, 5.0 }, new[] { 4.0, 1.0 });   // False: a trade-off
bool d3 = ParetoFront.Dominates(new[] { 2.0, 3.0 }, new[] { 2.0, 3.0 });   // False: equal points
```

Equal points do not dominate each other. The **non-dominated fronts** of Deb et al. [[15]](#15) are defined recursively: front $F_0$ holds the points no other point dominates, and front $F_{k+1}$ holds the points dominated only by points of $F_0, \dots, F_k$. Equivalently, the front of a point is one plus the largest front among its dominators. `NonDominatedSort` returns the fronts as lists of point indices, $F_0$ first, each in ascending index order; duplicate points share a front. `NonDominatedIndices` returns $F_0$ alone.

The points are first ordered lexicographically, so only earlier points can dominate a later one, and each point is placed by a binary search over the fronts built so far. With two objectives the sort is $O(n \log n)$, with three it is $O(n \log^2 n)$ comparisons, and with four or more objectives it is the $O(m n^2)$ bound of Deb et al. Every path gives exactly the fronts of the definition.

### Crowding Distance

Within one front, the **crowding distance** of Deb et al. [[15]](#15) measures how isolated a point is. For each objective $k$ the front is ordered by $f_k$ (ties by index); the first and last points get $+\infty$, and each interior point $i$ adds the normalized gap between its neighbors:

```math
d_i = \sum_{k=1}^{m} \frac{f_k^{(i+1)} - f_k^{(i-1)}}{f_k^{\max} - f_k^{\min}}
```

Each objective is normalized by its range over the front. An objective with zero range contributes nothing, boundary points included, and every point of a front with one or two points gets $+\infty$. pymoo's `calc_crowding_distance` uses the same terms but divides the sum by $m$, which does not change the ordering. By default (`filter_out_duplicates=True`) pymoo also removes duplicate points before computing the distances and gives each removed duplicate a distance of zero; `CrowdingDistance` keeps every point, so its values differ from pymoo's on a front with duplicates.

### Hypervolume

The **hypervolume** [[16]](#16) of a set $A$ with respect to a reference point $\mathbf{r}$ is the volume of the region its points dominate inside the box bounded by $\mathbf{r}$:

```math
HV(A, \mathbf{r}) = \lambda\left( \bigcup_{\mathbf{a} \in A,\; \mathbf{a} < \mathbf{r}} [a_1, r_1] \times \cdots \times [a_m, r_m] \right)
```

where $\lambda$ is the Lebesgue measure. Only points strictly better than $\mathbf{r}$ in every objective contribute; dominated and duplicate points change nothing. Larger is better, and the value depends on $\mathbf{r}$, so compare sets with the same reference point. The computation is exact in every dimension:

| Objectives | Method | Cost |
|---|---|---|
| 2 | Sort by $f_1$ and sweep the staircase | $O(n \log n)$ |
| 3 | Sweep in $f_3$, updating the area of the $(f_1, f_2)$ staircase by the exclusive area of each insertion [[17]](#17) | $O(n \log n)$ comparisons, at most $O(n^2)$ moves |
| 4 or more | WFG [[18]](#18): slice on $f_m$ and subtract the hypervolume of the non-dominated limit set from each point's box, down to the 3-objective sweep | exponential in $m$ in the worst case |

### Generational Distance and Inverted Generational Distance

Given an approximation set $A$ and a reference (true) front $Z$, the **generational distance** [[19]](#19) averages how far the approximation is from the front, and the **inverted generational distance** [[20]](#20) averages how far each part of the front is from the approximation:

```math
GD(A, Z) = \frac{1}{|A|} \sum_{\mathbf{a} \in A} \min_{\mathbf{z} \in Z} \lVert \mathbf{a} - \mathbf{z} \rVert_2, \qquad
IGD(A, Z) = \frac{1}{|Z|} \sum_{\mathbf{z} \in Z} \min_{\mathbf{a} \in A} \lVert \mathbf{z} - \mathbf{a} \rVert_2
```

Smaller is better for both. GD measures convergence only; IGD also penalizes gaps in coverage. These are the definitions of pymoo's `GD` and `IGD` without normalization. Van Veldhuizen and Lamont [[19]](#19) define GD with the root of the summed squared distances, $\sqrt{\sum d_i^2}/|A|$, which differs from the mean used here. Objectives are not normalized, so scale them first if their units differ.

```cs
using Numerics.Mathematics.Optimization;

// Six alternatives scored on (cost, risk); both are minimized.
var points = new double[][]
{
    new[] { 1.0, 5.0 },
    new[] { 2.0, 3.0 },
    new[] { 3.0, 4.0 },   // dominated by (2, 3)
    new[] { 4.0, 1.0 },
    new[] { 2.0, 3.0 },   // duplicate of point 1
    new[] { 5.0, 5.0 },   // dominated
};

var fronts = ParetoFront.NonDominatedSort(points);
for (int f = 0; f < fronts.Count; f++)
    Console.WriteLine($"Front {f}: {string.Join(", ", fronts[f])}");

// Crowding distance within the first front.
var distances = ParetoFront.CrowdingDistance(points, fronts[0]);
for (int i = 0; i < fronts[0].Count; i++)
    Console.WriteLine($"Point {fronts[0][i]}: crowding distance {distances[i]:F4}");

// Hypervolume with respect to the reference point (6, 6).
double hv = ParetoFront.Hypervolume(points, new[] { 6.0, 6.0 });
Console.WriteLine($"Hypervolume: {hv:F4}");

// Distance to a reference front on the segment f1 + f2 = 5.
var reference = new double[][]
{
    new[] { 0.0, 5.0 }, new[] { 1.0, 4.0 }, new[] { 2.0, 3.0 },
    new[] { 3.0, 2.0 }, new[] { 4.0, 1.0 }, new[] { 5.0, 0.0 },
};
var first = fronts[0].Select(i => points[i]).ToArray();
Console.WriteLine($"GD:  {ParetoFront.GenerationalDistance(first, reference):F4}");
Console.WriteLine($"IGD: {ParetoFront.InvertedGenerationalDistance(first, reference):F4}");
```

Output:

```
Front 0: 0, 1, 3, 4
Front 1: 2
Front 2: 5
Point 0: crowding distance ∞
Point 1: crowding distance 0.8333
Point 3: crowding distance ∞
Point 4: crowding distance 1.1667
Hypervolume: 17.0000
GD:  0.2500
IGD: 0.8047
```

The duplicate points 1 and 4 share the first front. Ties in each objective are ordered by index, so point 1 lies between points 0 and 4 and point 4 between points 1 and 3, which gives them different crowding distances. The output above is from .NET 10; .NET Framework prints infinity as `Infinity`.

## Shortest Path (Dijkstra)

The `Dijkstra`, `Network`, and `Edge` types solve destination-rooted shortest paths over
directed networks [[9]](#9). They are the routing kernel for agent-based evacuation modeling, where every
agent needs its route to the nearest destination and edge costs (travel times) change as
conditions evolve. The solvers run on an internal indexed binary min-heap sized at the node
count; the public `BinaryHeap<T>` is a separate standalone utility.

### The result table

Every solver returns one `float[,]` routing table with a row per node and three columns:

| Column | Contents |
|--------|----------|
| `[i, 0]` | The next node toward the destination |
| `[i, 1]` | The index of the edge to take (`Edge.Index`) |
| `[i, 2]` | The cumulative cost to the destination |

The destination row is `(itself, -1, 0)`; an unreachable node's row is `(-1, -1, +∞)`.
`Dijkstra.PathExists` tests reachability, and `Dijkstra.GetPath` / `Dijkstra.TryGetPath` walk a
table into the ordered edge-index list with the total cost.

Weight conventions: Dijkstra's algorithm requires **non-negative weights**: every `Dijkstra` solver
and every `Network` solve or path method throws `ArgumentOutOfRangeException` at entry when a
weight it would use is negative, which also rules out negative cycles; a weight of **positive infinity
makes an edge impassable** (the natural encoding for a flooded road segment), and a NaN weight
never relaxes, severing its edge.

Precision: the edge weights are `float`, but every solver accumulates costs in `double` and rounds
each cost to `float` once, in the table. The cost column is therefore the single-precision rounding
of the double-precision distance that `scipy.sparse.csgraph.dijkstra` computes on the same weights;
the unit tests check this bitwise.

```cs
using Numerics.Mathematics.Optimization;

var edges = new List<Edge>
{
    new Edge(0, 1, 2.0f, 0),   // from, to, weight (e.g. travel time), edge index
    new Edge(1, 2, 1.5f, 1),
    new Edge(0, 2, 5.0f, 2),
};

// Route every node to node 2.
float[,] table = Dijkstra.Solve(edges, destinationIndex: 2);

// Route every node to the nearest of several destinations in ONE pass.
float[,] nearest = Dijkstra.SolveNearest(edges, new[] { 0, 2 });

// Walk a route out of the table.
if (Dijkstra.TryGetPath(table, startNodeIndex: 0, out List<int> route, out float cost))
    Console.WriteLine($"Route: {string.Join(" -> ", route)}, cost {cost:F2}");
```

`Dijkstra.Solve(edges, destinationIndices, ...)` (the multi-destination overload) solves each
destination independently and keeps each node's strictly cheapest route, so on an exact cost tie
the earlier destination in the array wins. `SolveNearest` produces the same costs in a single
multi-source pass and is the faster form of the nearest-destination query.

### Multi-source distances on a CSR graph

`Dijkstra.ShortestDistances` takes a directed graph in compressed sparse row (CSR) form, with
double-precision weights, and returns the distance of every node from the nearest of several
sources. It follows `scipy.sparse.csgraph.dijkstra(graph, indices=sources, min_only=True)`: an
optional `limit` leaves farther nodes unreached, and optional arrays receive each node's
predecessor and nearest source (-1 where SciPy reports -9999). Weights must be non-negative
and not NaN.

```cs
// Edges 0->1 (0.1), 0->2 (0.3), 1->2 (0.2), 1->3 (1.2), 2->3 (0.7).
int[] rowStart = { 0, 2, 4, 5, 5 };
int[] toNode = { 1, 2, 2, 3, 3 };
double[] weight = { 0.1, 0.3, 0.2, 1.2, 0.7 };
var predecessors = new int[4];
var nearestSource = new int[4];
double[] distance = Dijkstra.ShortestDistances(4, rowStart, toNode, weight, new[] { 0 },
    predecessors: predecessors, nearestSource: nearestSource);
Console.WriteLine(string.Join(", ", distance));       // 0, 0.1, 0.3, 1
Console.WriteLine(string.Join(", ", predecessors));   // -1, 0, 0, 2
```

**Speed.** Median wall time for ten sources on random directed graphs with five edges per node and
weights uniform on $[0.001, 1.001)$, against `scipy.sparse.csgraph.dijkstra(min_only=True)`. Intel
Core Ultra 9 185H, Windows 11, .NET 10; Python 3.14.3, NumPy 2.5.4, SciPy 1.18.1, pynumerics
`graphs.shortest_distances` on a SciPy CSR matrix. BenchmarkDotNet (ShortRun,
`Benchmarks_Numerics/Optimization/DijkstraBenchmarks.cs`) gives 1.33 ms, 16.8 ms and 479 ms for the
three sizes.

| Nodes | Numerics (.NET) | pynumerics | SciPy |
|---|---|---|---|
| 10,000 | 1.03 ms | 1.17 ms | 1.65 ms |
| 100,000 | 15.4 ms | 17.8 ms | 24.6 ms |
| 1,000,000 | 473 ms | 573 ms | 624 ms |

On graphs larger than the processor cache the time goes to memory latency. The heap keeps each
node's distance next to its slot, so one cache line answers both the comparison and the heap
update, and after each removal the loop reads the edges of the node now at the top of the heap
while the current node's edges relax. Neither changes the settle order, so distances,
predecessors and nearest sources are the same bit for bit as a plain binary heap gives.

### Compiled networks and time-varying weights

`Network` compiles a fixed topology once (node count, incoming and outgoing adjacency, and the
destination set), so a simulation can re-solve every time step against updated edge weights with
no rebuild. Weights are positional with the constructor's edge array, and the table-reuse
overloads run with zero steady-state allocation (they share instance scratch buffers and are not
thread safe; the allocating overloads remain safe for concurrent use):

```cs
var network = new Network(edgeArray, destinationIndices);
var table = new float[network.NodeCount, 3];
var weights = new float[edgeArray.Length];

for (int t = 0; t < timeSteps; t++)
{
    UpdateTravelTimes(weights, t);                  // flooded edges -> float.PositiveInfinity
    network.SolveNearest(weights, table);           // zero-allocation re-solve
    RouteAgents(table);                             // each agent walks its row
}

// Detour routing around blocked segments, splicing onto the precomputed table when possible.
var detour = network.GetPath(blockedEdgeIndices, agentNodeIndex, table);
```

`Network.GetPath` finds the cheapest route to the nearest destination that avoids every edge
bearing an excluded edge index; when the precomputed table's recorded route is untouched by the
exclusions it is returned directly, with no solve. The table-based overload never returns null:
it returns an empty list when the start node is unreachable, is already a destination, or no
detour exists. The two-argument overload (without a table), by contrast, returns null when
every destination is unreachable.

## Discrete Dynamic Programming

`DynamicProgramming` solves sequential decision problems over integer states $s = 0, \dots, S-1$,
decisions $a = 0, \dots, A-1$ and stages $t = 0, \dots, T-1$ by Bellman's backward recursion [[10]](#10).
It is the classic tool for reservoir operation rule studies, in both its deterministic and stochastic
(Markov inflow) forms [[11]](#11), [[12]](#12), and it needs no solver dependency. `Knapsack` gives exact answers for
small budget-constrained project-selection problems.

With stage cost or reward $c_t(s, a)$, transition probabilities $p_t(s' \mid s, a)$, discount factor
$\gamma \in (0, 1]$ and terminal value $V_T(s)$, `Solve()` computes

```math
V_t(s) = \operatorname*{opt}_{a} \left[ c_t(s, a) + \gamma \sum_{s'} p_t(s' \mid s, a)\, V_{t+1}(s') \right], \qquad t = T-1, \dots, 0
```

where opt is min (the default) or max (`Maximize = true`). The work is $O(T \cdot S \cdot A \cdot K)$ for
at most $K$ outcomes per transition.

| Member | Meaning |
|--------|---------|
| `Value[t, s]` | Optimal value for $t = 0, \dots, T$; row $T$ is the terminal value |
| `Policy[t, s]` | Best decision for $t = 0, \dots, T-1$; `-1` where no decision is feasible |
| `DiscountFactor` | $\gamma$ in (0, 1], default 1 |
| `TerminalValue` | $V_T(s)$, one per state; `null` means zero |
| `UseParallel` | Evaluates the states of each stage in parallel; the result is identical to the serial solve |

Conventions:

- A **non-finite stage value** (NaN or $\pm\infty$) marks a decision infeasible; a deterministic transition
  returning **-1** does too. A non-finite terminal value marks that terminal state infeasible.
- A decision that can reach an infeasible state with positive probability is infeasible. A state with no
  feasible decision has value $+\infty$ when minimizing ($-\infty$ when maximizing) and policy `-1`.
- Ties go to the **lowest decision index**.
- The delegates must be deterministic functions of their arguments (and thread safe with `UseParallel`).

### Deterministic reservoir example

A reservoir holds 0 to 3 units, starts with 1, receives inflows of 2, 1 and 1 units over three seasons and
must end with at least 1 unit. Releasing $r$ units is worth $B(r) = (0, 5, 8, 9)$, and storage cannot
exceed capacity (no spill). The concave benefit spreads the 4 releasable units as 2, 1, 1 for a total of 18;
the three orders tie, and the lowest-index rule picks releases 1, 1, 2.

```cs
using Numerics.Mathematics.Optimization;

int[] inflow = { 2, 1, 1 };
double[] benefit = { 0, 5, 8, 9 };

var dp = new DynamicProgramming(stages: 3, states: 4, decisions: 4,
    stageValue: (t, s, r) => benefit[r],
    transition: (t, s, r) =>
    {
        int next = s + inflow[t] - r;
        return next < 0 || next > 3 ? -1 : next;   // -1 marks an infeasible release
    })
{
    Maximize = true,
    TerminalValue = new[] { double.NegativeInfinity, 0, 0, 0 }   // must end with at least 1 unit
};
dp.Solve();

var path = dp.Simulate(initialState: 1);
Console.WriteLine($"Total benefit: {path.TotalValue}");                // 18
Console.WriteLine($"Releases: {string.Join(", ", path.Decisions)}");   // 1, 1, 2
Console.WriteLine($"Storage:  {string.Join(", ", path.States)}");      // 1, 2, 2, 1
```

`Simulate` returns a `DynamicProgrammingPath` with the visited `States` (one per stage plus the terminal
state), the `Decisions`, the undiscounted `StageValues` and the discounted `TotalValue`. For a deterministic
problem `TotalValue` equals `Value[0, initialState]` exactly.

### Stochastic dynamic programming with Markov inflow classes

For stochastic transitions, pass a `TransitionDistribution` delegate and the largest number of outcomes.
The solver hands the delegate two reusable buffers to fill with next states and probabilities, so the inner
loop allocates nothing; the delegate returns the number of outcomes written (0 marks the decision infeasible).
Probabilities must be non-negative and sum to one within $10^{-8}$.

The usual reservoir formulation keeps the current inflow class in the state, $s = \text{storage} \cdot K + j$,
so the release is chosen knowing this season's inflow, and the next class follows a Markov matrix
$P_{j j'}$ [[11]](#11). The stage value is then the cost of the release given storage and class.

```cs
using Numerics.Mathematics.Optimization;

const int Storages = 11, Classes = 3, Releases = 6;
int[] classInflow = { 0, 2, 4 };
double[,] markov = { { 0.6, 0.3, 0.1 }, { 0.25, 0.5, 0.25 }, { 0.1, 0.3, 0.6 } };
int demand = 3;

var sdp = new DynamicProgramming(stages: 12, states: Storages * Classes, decisions: Releases,
    stageValue: (t, s, r) =>
    {
        int x = s / Classes + classInflow[s % Classes] - r;
        if (x < 0) return double.NaN;                        // cannot release more than is stored
        double shortage = Math.Max(demand - r, 0);
        double spill = Math.Max(x - (Storages - 1), 0);
        return shortage * shortage + 0.5 * spill;
    },
    transition: (t, s, r, nextStates, probabilities) =>
    {
        int j = s % Classes;
        int x = s / Classes + classInflow[j] - r;
        if (x < 0) return 0;                                 // infeasible
        int storage = Math.Min(x, Storages - 1);
        for (int jn = 0; jn < Classes; jn++)
        {
            nextStates[jn] = storage * Classes + jn;
            probabilities[jn] = markov[j, jn];
        }
        return Classes;
    },
    maxOutcomes: Classes)
{
    DiscountFactor = 0.98
};
sdp.Solve();

int start = 5 * Classes + 1;   // storage 5, middle inflow class
Console.WriteLine($"Expected cost: {sdp.Value[0, start]:F4}");

// Replay the policy along an observed class sequence, or sample classes with a seeded generator.
var observed = sdp.Simulate(start, new[] { 2, 2, 1, 0, 0, 0, 1, 2, 2, 1, 1, 0 });
var sampled = sdp.Simulate(start, new Numerics.Sampling.MersenneTwister(12345));
Console.WriteLine($"Observed path cost: {observed.TotalValue:F4}, sampled: {sampled.TotalValue:F4}");
```

### Steady-state policies by value iteration

For a stationary infinite-horizon problem, `SolveSteadyState(tolerance, maxIterations)` applies value
iteration with the delegates evaluated once at stage 0 [[13]](#13):

```math
V_{k+1}(s) = \operatorname*{opt}_{a} \left[ c(s, a) + \gamma \sum_{s'} p(s' \mid s, a)\, V_k(s') \right]
```

starting from `TerminalValue` (zero if null) and stopping when $\max_s |V_{k+1}(s) - V_k(s)| \le$ `tolerance`
(default $10^{-8}$, at most 10,000 sweeps). With $\gamma < 1$ the error of the last iterate is at most
$\gamma / (1 - \gamma)$ times the last change. It returns whether it converged and fills `SteadyStateValue`,
`SteadyStatePolicy`, `SteadyStateIterations` and `SteadyStateConverged`. With $\gamma = 1$ it converges only
when values stay bounded, as in a stochastic shortest-path problem with a free absorbing state.

```cs
using Numerics.Mathematics.Optimization;

// Two states, two decisions, discounted cost.
double[,] cost = { { 2.0, 1.0 }, { 3.0, 0.5 } };
double[,,] p = { { { 0.9, 0.1 }, { 0.2, 0.8 } }, { { 0.7, 0.3 }, { 0.05, 0.95 } } };

var mdp = new DynamicProgramming(1, 2, 2, (t, s, a) => cost[s, a],
    (t, s, a, next, prob) =>
    {
        next[0] = 0; next[1] = 1;
        prob[0] = p[s, a, 0]; prob[1] = p[s, a, 1];
        return 2;
    }, maxOutcomes: 2)
{
    DiscountFactor = 0.9
};

bool converged = mdp.SolveSteadyState(tolerance: 1e-10);
Console.WriteLine($"Converged: {converged} after {mdp.SteadyStateIterations} sweeps");
Console.WriteLine($"V = ({mdp.SteadyStateValue[0]:F6}, {mdp.SteadyStateValue[1]:F6}), " +
                  $"policy = ({mdp.SteadyStatePolicy[0]}, {mdp.SteadyStatePolicy[1]})");
// V = (5.838150, 5.260116), policy = (1, 1): the exact solution of the linear Bellman system of policy (1, 1)
```

### Knapsack (project selection)

`Knapsack.Solve` maximizes $\sum_i v_i x_i$ subject to $\sum_i w_i x_i \le W$, $x_i \in \{0, 1\}$, with
integer weights, by the recursion

```math
f_i(c) = \max\left( f_{i-1}(c),\; f_{i-1}(c - w_i) + v_i \right), \qquad c = 0, \dots, W
```

in $O(n \cdot W)$ time and $O(n \cdot W)$ bytes for the decision table [[10]](#10). It returns the optimal value and
the selected item indices in increasing order. An item is selected only when it strictly increases the value,
so items with non-positive value are never selected, and on a tie the later item is left out. The
multi-constraint overload takes `weights[i, k]` and one capacity per constraint (for example one budget per
period); it runs over the product of all capacities, $O(n \cdot m \cdot \prod_k (W_k + 1))$, so it is meant
for small instances and as an exact benchmark for heuristic portfolio optimizers.

```cs
using Numerics.Mathematics.Optimization;

// Risk reduction (value) and cost in $100k of five candidate projects, budget $1M.
double[] riskReduction = { 4.2, 2.9, 6.1, 1.5, 3.4 };
int[] cost = { 4, 3, 6, 1, 3 };

double best = Knapsack.Solve(riskReduction, cost, capacity: 10, out int[] selected);
Console.WriteLine($"Best: {best:F1} from projects {string.Join(", ", selected)}");   // 11.0 from 2, 3, 4

// Two budget periods.
int[,] periodCost = { { 2, 2 }, { 1, 2 }, { 4, 2 }, { 1, 0 }, { 2, 1 } };
double best2 = Knapsack.Solve(riskReduction, periodCost, new[] { 5, 4 }, out int[] selected2);
Console.WriteLine($"Best: {best2:F1} from projects {string.Join(", ", selected2)}");  // 9.1 from 0, 3, 4
```

## Practical Example: Calibrating a Hydrological Model

A complete example of using optimization to calibrate a watershed model:

```cs
using Numerics.Mathematics.Optimization;
using Numerics.Data.Statistics;

// Observed streamflow data
double[] observed = { 12.5, 15.3, 18.7, 16.2, 14.1, 11.8, 10.3 };

// Simple runoff model: Q = C × P^α where Q is flow, P is precipitation, C and α are parameters
double[] precipitation = { 10, 12, 15, 13, 11, 9, 8 };

// Objective function: minimize RMSE
double ObjectiveFunction(double[] parameters)
{
    double C = parameters[0];
    double alpha = parameters[1];
    
    // Simulate streamflow
    double[] simulated = new double[observed.Length];
    for (int i = 0; i < observed.Length; i++)
    {
        simulated[i] = C * Math.Pow(precipitation[i], alpha);
    }
    
    // Compute RMSE
    double rmse = GoodnessOfFit.RMSE(observed, simulated);
    return rmse;
}

// Parameter bounds
var lower = new double[] { 0.1, 0.5 };  // C >= 0.1, α >= 0.5
var upper = new double[] { 5.0, 3.0 };  // C <= 5.0, α <= 3.0

// Use SCE-UA (recommended for hydrological calibration)
var optimizer = new ShuffledComplexEvolution(ObjectiveFunction, 2, lower, upper);
optimizer.Complexes = 5;
optimizer.MaxIterations = 1000;
optimizer.Minimize();

Console.WriteLine($"Calibrated parameters:");
Console.WriteLine($"  C = {optimizer.BestParameterSet.Values[0]:F4}");
Console.WriteLine($"  α = {optimizer.BestParameterSet.Values[1]:F4}");
Console.WriteLine($"  RMSE = {optimizer.BestParameterSet.Fitness:F4}");
Console.WriteLine($"  Function evaluations: {optimizer.FunctionEvaluations}");

// Verify calibration
double[] final_simulated = new double[observed.Length];
for (int i = 0; i < observed.Length; i++)
{
    final_simulated[i] = optimizer.BestParameterSet.Values[0] * 
                         Math.Pow(precipitation[i], optimizer.BestParameterSet.Values[1]);
}

double nse = GoodnessOfFit.NashSutcliffeEfficiency(observed, final_simulated);
Console.WriteLine($"  NSE = {nse:F4}");
```

## Choosing an Optimization Method

| Problem Type | Recommended Method | Notes |
|-------------|-------------------|-------|
| Smooth, unimodal | BFGS or Powell | Fast convergence to local minimum |
| Non-smooth, unimodal | Nelder-Mead | Robust direct search |
| Multimodal, global search | Differential Evolution or SCE-UA | Thorough exploration |
| Quick global search | Particle Swarm | Faster but less thorough |
| Hydrological calibration | SCE-UA | Specifically designed for this |
| Expensive simulation calibration | DDS | Scales to a fixed run budget; polish for precision |
| Correlated or badly scaled parameters | CMA-ES | Learns the covariance; IPOP restarts for multimodal problems |
| Staged decisions (reservoir rules) | Dynamic Programming | Exact on a discretized state space |
| Budget-constrained selection | Knapsack | Exact for integer costs |
| With constraints | Augmented Lagrangian | Handles equality and inequality constraints |
| 1D problem | Brent Search or Golden Section | Specialized efficient methods |
| Machine learning | ADAM | Adaptive learning rates |
| Unknown smoothness | Start with Nelder-Mead | Very robust |
| High dimensions (>20) | Differential Evolution or ADAM | Scale better than others |

## Best Practices

1. **Scale Variables**: Normalize parameters to similar ranges (e.g., [0,1] or [-1,1]) for better convergence.

2. **Set Reasonable Bounds**: Tight bounds improve convergence but shouldn't exclude the optimum.

3. **Multiple Runs**: For global optimization, run multiple times with different random seeds to ensure robustness.

4. **Hybrid Approach**: Use global method followed by local refinement:
   ```cs
   // Global search
   var de = new DifferentialEvolution(ObjectiveFunc, n, lower, upper);
   de.Minimize();
   
   // Local refinement
   var bfgs = new BFGS(ObjectiveFunc, n, de.BestParameterSet.Values, lower, upper);
   bfgs.Minimize();
   ```

5. **Monitor Convergence**: Check the `ParameterSetTrace` to diagnose convergence issues. It holds one entry per function evaluation, each the best point found so far:
   ```cs
   for (int k = 0; k < optimizer.ParameterSetTrace.Count; k++)
   {
       Console.WriteLine($"Evaluation {k + 1}: best fitness so far {optimizer.ParameterSetTrace[k].Fitness}");
   }
   ```

6. **Adjust Tolerances**: Tighter tolerances require more evaluations:
   ```cs
   optimizer.RelativeTolerance = 1e-10; // Very tight
   optimizer.AbsoluteTolerance = 1e-10;
   ```

7. **Population Size**: For global methods, larger populations explore better but cost more:
   ```cs
   de.PopulationSize = 20 * NumberOfParameters; // Rule of thumb
   ```

## Understanding Results

### Parameter Set Structure

The `BestParameterSet` contains:
- `Values`: The optimal parameter vector
- `Fitness`: The objective function value at the optimum
- `Weight`: (Optional) Used internally by some algorithms

### Hessian Matrix

When `ComputeHessian = true`, the Hessian at the solution is computed numerically. The Hessian $\mathbf{H}$ is the matrix of second partial derivatives of the objective function, and its eigenvalues reveal important information about the solution:

- **All eigenvalues positive**: The solution is a local minimum (the Hessian is positive definite)
- **Large eigenvalue**: The objective is highly curved in that direction, so the corresponding parameter is well-determined
- **Small eigenvalue**: The objective is nearly flat, so the parameter is poorly determined or identifiable
- **Near-zero eigenvalue**: Indicates a ridge or valley in the objective surface, suggesting parameter correlation or redundancy

When the objective is a negative log-likelihood, the inverse of the Hessian at the minimum approximates the covariance matrix of the parameters, so standard errors are $SE(\hat{\theta}_i) \approx \sqrt{(\mathbf{H}^{-1})_{ii}}$. Invert the whole matrix: $1/H_{ii}$ is not the same quantity unless the parameters are uncorrelated. When a log-likelihood is maximized instead, the Hessian is that of the log-likelihood itself, so use $-\mathbf{H}^{-1}$.

```cs
optimizer.ComputeHessian = true;
optimizer.Minimize();

// Covariance approximation from the inverse Hessian
var covariance = optimizer.Hessian!.Inverse();

// Parameter standard errors (approximation)
for (int i = 0; i < optimizer.NumberOfParameters; i++)
{
    double se = Math.Sqrt(covariance[i, i]);
    Console.WriteLine($"Parameter {i} ± {se:F4}");
}
```

### Optimization Status

Check the `Status` property to verify success:

```cs
if (optimizer.Status == OptimizationStatus.Success)
{
    Console.WriteLine("Optimization converged successfully");
}
else if (optimizer.Status == OptimizationStatus.MaximumIterationsReached)
{
    Console.WriteLine("Maximum iterations reached - may not have converged");
}
```

## Performance Tips

1. **Function Evaluation Cost**: Most optimization time is spent evaluating the objective function. Optimize your function first.

2. **Parallel Evaluation**: Differential Evolution and Particle Swarm can evaluate each generation in parallel (`UseParallel`), reproducibly for any number of threads. The other optimizers evaluate one point at a time; for an expensive objective, parallelize inside the objective function itself (for example, across the events or reaches of one simulation).

3. **Warm Start**: If solving similar problems repeatedly, use the previous solution as initial guess.

4. **Gradient Information**: `BFGS`, `GradientDescent`, and `ADAM` accept an optional analytic gradient in their constructors; without one they use finite differences, which cost about $2D$ extra evaluations per gradient.

5. **Overhead**: Differential Evolution, Simulated Annealing, Particle Swarm, Powell, Brent, Multi-Start, and MLSL reuse their working arrays rather than allocating per trial point, so the optimizer's own cost is small next to even a cheap objective.

---

## References

<a id="1">[1]</a> Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer.

<a id="2">[2]</a> Nelder, J. A., & Mead, R. (1965). A simplex method for function minimization. *The Computer Journal*, 7(4), 308-313.

<a id="3">[3]</a> Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic optimization. *arXiv preprint arXiv:1412.6980*.

<a id="4">[4]</a> Storn, R., & Price, K. (1997). Differential evolution: A simple and efficient heuristic for global optimization over continuous spaces. *Journal of Global Optimization*, 11(4), 341-359.

<a id="5">[5]</a> Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. *Proceedings of IEEE International Conference on Neural Networks*, 4, 1942-1948.

<a id="6">[6]</a> Duan, Q., Sorooshian, S., & Gupta, V. (1992). Effective and efficient global optimization for conceptual rainfall-runoff models. *Water Resources Research*, 28(4), 1015-1031.

<a id="7">[7]</a> Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. *Science*, 220(4598), 671-680.

<a id="8">[8]</a> Birgin, E. G., & Martínez, J. M. (2014). *Practical Augmented Lagrangian Methods for Constrained Optimization*. SIAM.

<a id="9">[9]</a> Dijkstra, E. W. (1959). A note on two problems in connexion with graphs. *Numerische Mathematik*, 1, 269-271.

<a id="10">[10]</a> Bellman, R. (1957). *Dynamic Programming*. Princeton University Press.

<a id="11">[11]</a> Loucks, D. P., & van Beek, E. (2017). *Water Resource Systems Planning and Management*. Springer. https://doi.org/10.1007/978-3-319-44234-1

<a id="12">[12]</a> Labadie, J. W. (2004). Optimal operation of multireservoir systems: State-of-the-art review. *Journal of Water Resources Planning and Management*, 130(2), 93-111. https://doi.org/10.1061/(ASCE)0733-9496(2004)130:2(93)

<a id="13">[13]</a> Puterman, M. L. (1994). *Markov Decision Processes: Discrete Stochastic Dynamic Programming*. Wiley. https://doi.org/10.1002/9780470316887

<a id="14">[14]</a> Tolson, B. A., & Shoemaker, C. A. (2007). Dynamically dimensioned search algorithm for computationally efficient watershed model calibration. *Water Resources Research*, 43, W01413. https://doi.org/10.1029/2005WR004723

<a id="15">[15]</a> Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on Evolutionary Computation*, 6(2), 182-197. https://doi.org/10.1109/4235.996017

<a id="16">[16]</a> Zitzler, E., & Thiele, L. (1999). Multiobjective evolutionary algorithms: A comparative case study and the strength Pareto approach. *IEEE Transactions on Evolutionary Computation*, 3(4), 257-271. https://doi.org/10.1109/4235.797969

<a id="17">[17]</a> Beume, N., Fonseca, C. M., López-Ibáñez, M., Paquete, L., & Vahrenhold, J. (2009). On the complexity of computing the hypervolume indicator. *IEEE Transactions on Evolutionary Computation*, 13(5), 1075-1082. https://doi.org/10.1109/TEVC.2009.2015575

<a id="18">[18]</a> While, L., Bradstreet, L., & Barone, L. (2012). A fast way of calculating exact hypervolumes. *IEEE Transactions on Evolutionary Computation*, 16(1), 86-95. https://doi.org/10.1109/TEVC.2010.2077298

<a id="19">[19]</a> Van Veldhuizen, D. A., & Lamont, G. B. (2000). On measuring multiobjective evolutionary algorithm performance. *Proceedings of the 2000 Congress on Evolutionary Computation*, 204-211. https://doi.org/10.1109/CEC.2000.870296

<a id="20">[20]</a> Coello Coello, C. A., & Reyes Sierra, M. (2004). A study of the parallelization of a coevolutionary multi-objective evolutionary algorithm. In *MICAI 2004: Advances in Artificial Intelligence*, Lecture Notes in Computer Science (pp. 688-697). Springer. https://doi.org/10.1007/978-3-540-24694-7_71

<a id="21">[21]</a> Rinnooy Kan, A. H. G., & Timmer, G. T. (1987). Stochastic global optimization methods part II: Multi level methods. *Mathematical Programming*, 39, 57-78. https://doi.org/10.1007/BF02592071

<a id="22">[22]</a> Duan, Q., Sorooshian, S., & Gupta, V. K. (1994). Optimal use of the SCE-UA global optimization method for calibrating watershed models. *Journal of Hydrology*, 158(3-4), 265-284. https://doi.org/10.1016/0022-1694(94)90057-4

<a id="23">[23]</a> Hansen, N. (2023). The CMA evolution strategy: A tutorial. *arXiv preprint arXiv:1604.00772v2*.

<a id="24">[24]</a> Hansen, N., & Ostermeier, A. (2001). Completely derandomized self-adaptation in evolution strategies. *Evolutionary Computation*, 9(2), 159-195.

<a id="25">[25]</a> Auger, A., & Hansen, N. (2005). A restart CMA evolution strategy with increasing population size. *2005 IEEE Congress on Evolutionary Computation*, 2, 1769-1776.

<a id="26">[26]</a> Hansen, N., Akimoto, Y., Baudis, P., et al. CMA-ES/pycma (version 4.5.0). Zenodo. https://doi.org/10.5281/zenodo.2559634

<a id="27">[27]</a> Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. *SIAM Journal on Scientific Computing*, 16(5), 1190-1208. https://doi.org/10.1137/0916069

<a id="28">[28]</a> Zhu, C., Byrd, R. H., Lu, P., & Nocedal, J. (1997). Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale bound-constrained optimization. *ACM Transactions on Mathematical Software*, 23(4), 550-560. https://doi.org/10.1145/279232.279236

<a id="29">[29]</a> Morales, J. L., & Nocedal, J. (2011). Remark on "Algorithm 778: L-BFGS-B: Fortran subroutines for large-scale bound constrained optimization". *ACM Transactions on Mathematical Software*, 38(1), 7. https://doi.org/10.1145/2049662.2049669

<a id="30">[30]</a> Moré, J. J., & Thuente, D. J. (1994). Line search algorithms with guaranteed sufficient decrease. *ACM Transactions on Mathematical Software*, 20(3), 286-307. https://doi.org/10.1145/192115.192132

---

[← Previous: Numerical Differentiation](differentiation.md) | [Back to Index](../index.md) | [Next: Root Finding →](root-finding.md)
