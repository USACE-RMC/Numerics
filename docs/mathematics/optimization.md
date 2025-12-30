# Optimization

[← Previous: Numerical Differentiation](differentiation.md) | [Back to Index](../index.md) | [Next: Linear Algebra →](linear-algebra.md)

Optimization is the process of finding the parameter set that minimizes (or maximizes) an objective function. The ***Numerics*** library provides a comprehensive suite of local and global optimization algorithms for both unconstrained and constrained problems. These methods are essential for parameter estimation, model calibration, machine learning, and engineering design optimization.

## Overview

| Method | Type | Best For | Requires Derivatives |
|--------|------|----------|---------------------|
| **Local Methods** | | | |
| BFGS | Local | Smooth, differentiable functions | No (numerical) |
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
| **Constrained** | | | |
| Augmented Lagrangian | Constrained | Equality and inequality constraints | No |

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
- `ReportFailure`: Throw exception on failure (default = true)
- `RecordTraces`: Save optimization trace (default = true)
- `ComputeHessian`: Compute Hessian at solution (default = true)

### Output Properties
- `BestParameterSet`: Optimal solution (Values, Fitness)
- `Iterations`: Number of iterations performed
- `FunctionEvaluations`: Number of function evaluations
- `Status`: Optimization status (Success, Failure, etc.)
- `ParameterSetTrace`: Full trace of parameter evaluations
- `Hessian`: Numerically differentiated Hessian matrix (if computed)

### Methods
- `Minimize()`: Minimize the objective function
- `Maximize()`: Maximize the objective function (minimizes $-f(\mathbf{x})$)

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

// Initial guess
var initial = new double[] { -1.2, 1.0 };
var lower = new double[] { -2, -2 };
var upper = new double[] { 2, 2 };

// Create and configure optimizer
var bfgs = new BFGS(Rosenbrock, 2, initial, lower, upper);
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

**Advantages**: Fast convergence, memory efficient, works well for most smooth problems.

**Disadvantages**: Can fail on non-smooth functions, requires bounded search space.

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

Simple gradient-based optimization with line search:

```cs
var gd = new GradientDescent(Rosenbrock, 2, initial, lower, upper);
gd.LearningRate = 0.001; // Step size
gd.Minimize();
```

### ADAM Optimizer

Adaptive Moment Estimation, popular in machine learning applications [[3]](#3):

```cs
var adam = new ADAM(Rosenbrock, 2, initial, lower, upper);
adam.LearningRate = 0.001;
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

var golden = new GoldenSection(f1d, 1);
golden.LowerBounds = new[] { 0.0 };
golden.UpperBounds = new[] { 5.0 };
golden.Minimize();

Console.WriteLine($"Minimum at x = {golden.BestParameterSet.Values[0]:F6}");
```

#### Brent Search

Combines golden section search with parabolic interpolation for faster convergence:

```cs
var brent = new BrentSearch(f1d, 1);
brent.LowerBounds = new[] { 0.0 };
brent.UpperBounds = new[] { 5.0 };
brent.Minimize();
```

## Global Optimization

Global optimization methods are designed to find the global minimum across the entire search space, avoiding local minima. They typically require more function evaluations but are more robust for multimodal problems.

### Differential Evolution

Differential Evolution (DE) is a population-based evolutionary algorithm that's very robust for continuous optimization [[4]](#4). It creates trial vectors by combining existing population members.

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
de.DifferentialWeight = 0.8;
de.MaxIterations = 1000;
de.Minimize();

Console.WriteLine($"Global minimum: [{de.BestParameterSet.Values[0]:F6}, {de.BestParameterSet.Values[1]:F6}]");
Console.WriteLine($"Function value: {de.BestParameterSet.Fitness:F10}");
```

**Advantages**: Very robust, handles discontinuities, good for difficult problems.

**Parameters**:
- `PopulationSize`: Number of candidate solutions (default = 10 × dimensions)
- `CrossoverProbability`: Probability of crossover (default = 0.9)
- `DifferentialWeight`: Scaling factor for mutation (default = 0.8)

### Particle Swarm Optimization

PSO simulates social behavior of bird flocking or fish schooling [[5]](#5). Particles move through the search space influenced by their own best position and the swarm's best position.

```cs
var pso = new ParticleSwarm(Rastrigin, 2, lower, upper);
pso.PopulationSize = 40;
pso.InertiaWeight = 0.7;
pso.CognitiveWeight = 1.5; // Personal best influence
pso.SocialWeight = 1.5;    // Global best influence
pso.Minimize();

Console.WriteLine($"Solution: [{pso.BestParameterSet.Values[0]:F6}, {pso.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Fast, simple to implement, works well for continuous problems.

**Parameters**:
- `PopulationSize`: Number of particles (default = 10 × dimensions)
- `InertiaWeight`: Momentum (default = 0.7)
- `CognitiveWeight`: Personal best attraction (default = 1.5)
- `SocialWeight`: Global best attraction (default = 1.5)

### Shuffled Complex Evolution (SCE-UA)

SCE-UA was specifically developed for calibrating hydrological models [[6]](#6). It combines complex shuffling with competitive evolution.

```cs
var sce = new ShuffledComplexEvolution(Rastrigin, 2, lower, upper);
sce.NumberOfComplexes = 5;
sce.ComplexSize = 10;
sce.MaxIterations = 1000;
sce.Minimize();

Console.WriteLine($"Calibrated parameters: [{sce.BestParameterSet.Values[0]:F6}, {sce.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Excellent for hydrological model calibration, balances exploration and exploitation.

**Best for**: Calibrating watershed models, water resources applications.

### Simulated Annealing

SA mimics the physical process of annealing in metallurgy [[7]](#7). It accepts uphill moves with decreasing probability, allowing escape from local minima.

```cs
var sa = new SimulatedAnnealing(Rastrigin, 2, lower, upper);
sa.InitialTemperature = 100.0;
sa.CoolingRate = 0.95;
sa.MaxIterations = 10000;
sa.Minimize();

Console.WriteLine($"Solution: [{sa.BestParameterSet.Values[0]:F6}, {sa.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Can escape local minima, works for discrete and continuous problems.

**Parameters**:
- `InitialTemperature`: Starting temperature (default = 100)
- `CoolingRate`: Temperature reduction factor (default = 0.95)

### Multi-Start Optimization

Combines local search with multiple random starting points:

```cs
var ms = new MultiStart(Rastrigin, 2, lower, upper);
ms.LocalMethod = LocalMethod.BFGS; // Choose local optimizer
ms.NumberOfStarts = 20;
ms.Minimize();

Console.WriteLine($"Best solution: [{ms.BestParameterSet.Values[0]:F6}, {ms.BestParameterSet.Values[1]:F6}]");
```

**Advantages**: Simple to implement, leverages fast local search.

**Disadvantages**: May waste evaluations in same basin of attraction.

### MLSL (Multi-Level Single Linkage)

Clustering-based global optimization that avoids redundant local searches:

```cs
var mlsl = new MLSL(Rastrigin, 2, lower, upper);
mlsl.LocalMethod = LocalMethod.BFGS;
mlsl.Minimize();
```

**Advantages**: More efficient than multi-start, avoids redundant searches.

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
    x => x[0] + x[1] - 1,  // g(x) >= 0 form
    ConstraintType.GreaterThanOrEqualTo
);

var lower = new double[] { -5, -5 };
var upper = new double[] { 5, 5 };
var initial = new double[] { 0, 0 };

var al = new AugmentedLagrange(Objective, 2, initial, lower, upper);
al.AddConstraint(constraint);
al.LocalMethod = LocalMethod.BFGS; // Local optimizer for subproblems
al.MaxIterations = 100;
al.Minimize();

Console.WriteLine($"Optimal solution: [{al.BestParameterSet.Values[0]:F6}, {al.BestParameterSet.Values[1]:F6}]");
Console.WriteLine($"Constraint satisfied: {al.BestParameterSet.Values[0] + al.BestParameterSet.Values[1]:F6} >= 1");
```

**Constraint Types**:
- `ConstraintType.EqualTo`: Equality constraint $g(\mathbf{x}) = 0$
- `ConstraintType.LessThanOrEqualTo`: Inequality constraint $g(\mathbf{x}) \leq 0$
- `ConstraintType.GreaterThanOrEqualTo`: Inequality constraint $g(\mathbf{x}) \geq 0$

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

var c1 = new Constraint(x => 5 - x[0] - x[1], ConstraintType.GreaterThanOrEqualTo);
var c2 = new Constraint(x => x[0] - 1, ConstraintType.GreaterThanOrEqualTo);
var c3 = new Constraint(x => x[1] - 1, ConstraintType.GreaterThanOrEqualTo);

var constrained = new AugmentedLagrange(ObjectiveFunc, 2, new[] { 2.0, 2.0 }, 
                                        new[] { 0.0, 0.0 }, new[] { 10.0, 10.0 });
constrained.AddConstraint(c1);
constrained.AddConstraint(c2);
constrained.AddConstraint(c3);
constrained.Minimize();

Console.WriteLine($"Constrained optimum: [{constrained.BestParameterSet.Values[0]:F4}, " +
                  $"{constrained.BestParameterSet.Values[1]:F4}]");
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
    double rmse = Statistics.RMSE(observed, simulated);
    return rmse;
}

// Parameter bounds
var lower = new double[] { 0.1, 0.5 };  // C >= 0.1, α >= 0.5
var upper = new double[] { 5.0, 3.0 };  // C <= 5.0, α <= 3.0

// Use SCE-UA (recommended for hydrological calibration)
var optimizer = new ShuffledComplexEvolution(ObjectiveFunction, 2, lower, upper);
optimizer.NumberOfComplexes = 5;
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

double nse = Statistics.NSE(observed, final_simulated);
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

5. **Monitor Convergence**: Check the `ParameterSetTrace` to diagnose convergence issues:
   ```cs
   foreach (var ps in optimizer.ParameterSetTrace)
   {
       Console.WriteLine($"Iteration {ps.Fitness}");
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

When `ComputeHessian = true`, the Hessian at the solution is computed numerically. This provides information about parameter sensitivity and uncertainty:

```cs
optimizer.ComputeHessian = true;
optimizer.Minimize();

// Access Hessian
var H = optimizer.Hessian;

// Parameter standard errors (approximation)
for (int i = 0; i < optimizer.NumberOfParameters; i++)
{
    double se = Math.Sqrt(Math.Abs(1.0 / H[i, i]));
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
else if (optimizer.Status == OptimizationStatus.MaxIterationsReached)
{
    Console.WriteLine("Maximum iterations reached - may not have converged");
}
```

## Performance Tips

1. **Function Evaluation Cost**: Most optimization time is spent evaluating the objective function. Optimize your function first.

2. **Parallel Evaluation**: For population-based methods (DE, PSO), evaluate population members in parallel if your function is thread-safe.

3. **Warm Start**: If solving similar problems repeatedly, use the previous solution as initial guess.

4. **Gradient Information**: If you can provide analytical gradients, BFGS will converge much faster (though the current implementation uses numerical gradients).

---

## References

<a id="1">[1]</a> Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer.

<a id="2">[2]</a> Nelder, J. A., & Mead, R. (1965). A simplex method for function minimization. *The Computer Journal*, 7(4), 308-313.

<a id="3">[3]</a> Kingma, D. P., & Ba, J. (2014). Adam: A method for stochastic optimization. *arXiv preprint arXiv:1412.6980*.

<a id="4">[4]</a> Storn, R., & Price, K. (1997). Differential evolution–a simple and efficient heuristic for global optimization over continuous spaces. *Journal of Global Optimization*, 11(4), 341-359.

<a id="5">[5]</a> Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. *Proceedings of IEEE International Conference on Neural Networks*, 4, 1942-1948.

<a id="6">[6]</a> Duan, Q., Sorooshian, S., & Gupta, V. (1992). Effective and efficient global optimization for conceptual rainfall-runoff models. *Water Resources Research*, 28(4), 1015-1031.

<a id="7">[7]</a> Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. *Science*, 220(4598), 671-680.

<a id="8">[8]</a> Birgin, E. G., & Martínez, J. M. (2014). *Practical Augmented Lagrangian Methods for Constrained Optimization*. SIAM.

---

[← Previous: Numerical Differentiation](differentiation.md) | [Back to Index](../index.md) | [Next: Linear Algebra →](linear-algebra.md)
