# Optimization

The ***Numerics*** library provides comprehensive optimization algorithms for local and global optimization, including derivative-free methods, quasi-Newton methods, and evolutionary algorithms.

## Overview

### Local Optimization

| Method | Type | Best For |
|--------|------|----------|
| Golden Section | Derivative-free | 1D unimodal functions |
| Brent's Method | Derivative-free | 1D general |
| Nelder-Mead | Derivative-free | Low-dimensional, noisy |
| Powell's Method | Derivative-free | Medium-dimensional |
| BFGS | Quasi-Newton | Smooth functions |

### Global Optimization

| Method | Type | Best For |
|--------|------|----------|
| Differential Evolution | Evolutionary | General purpose |
| SCE-UA | Evolutionary | Hydrologic calibration |
| Particle Swarm | Swarm intelligence | Continuous problems |
| Simulated Annealing | Probabilistic | Combinatorial, continuous |

---

## Local Optimization

### Golden Section Search

Finds the minimum of a unimodal function on an interval [[1]](#ref1):

```cs
using Numerics.Mathematics.Optimization;

Func<double, double> f = x => (x - 2) * (x - 2);

var golden = new GoldenSectionSearch(f, 0, 5);
golden.Minimize();

Console.WriteLine($"Minimum at x = {golden.BestParameterSet[0]:F10}");
Console.WriteLine($"Function value = {golden.BestFitness:F10}");
```

**Convergence**: Linear, reduces interval by factor of 0.618 per iteration.

### Brent's Method

Combines golden section with parabolic interpolation [[1]](#ref1):

```cs
Func<double, double> f = x => Math.Sin(x) - x / 3;

var brent = new BrentSearch(f, 0, 5);
brent.Tolerance = 1e-10;
brent.Minimize();

Console.WriteLine($"Minimum at x = {brent.BestParameterSet[0]:F10}");
```

**Advantages**: Superlinear convergence for smooth functions, guaranteed convergence.

### Nelder-Mead (Simplex)

Derivative-free method using simplex transformations [[2]](#ref2):

```cs
// Rosenbrock function
Func<double[], double> rosenbrock = x =>
{
    double a = 1 - x[0];
    double b = x[1] - x[0] * x[0];
    return a * a + 100 * b * b;
};

double[] initialGuess = { -1, -1 };

var nm = new NelderMead(rosenbrock, initialGuess.Length, initialGuess);
nm.Tolerance = 1e-10;
nm.MaxIterations = 10000;
nm.Minimize();

Console.WriteLine($"Minimum at: ({nm.BestParameterSet[0]:F6}, {nm.BestParameterSet[1]:F6})");
Console.WriteLine($"Function value: {nm.BestFitness:E6}");
Console.WriteLine($"Iterations: {nm.Iterations}");
```

**Advantages**: No derivatives required, handles noisy functions.
**Limitations**: Slow for high dimensions, can stall on non-smooth functions.

### Powell's Method

Conjugate direction method without derivatives [[1]](#ref1):

```cs
var powell = new Powell(rosenbrock, 2, initialGuess);
powell.Tolerance = 1e-10;
powell.Minimize();

Console.WriteLine($"Minimum at: ({powell.BestParameterSet[0]:F6}, {powell.BestParameterSet[1]:F6})");
```

**Advantages**: Faster than Nelder-Mead for smooth functions.

### BFGS (Quasi-Newton)

Uses approximate Hessian for fast convergence [[3]](#ref3):

```cs
var bfgs = new BFGS(rosenbrock, 2, initialGuess);
bfgs.Tolerance = 1e-12;
bfgs.MaxIterations = 1000;
bfgs.Minimize();

Console.WriteLine($"Minimum at: ({bfgs.BestParameterSet[0]:F10}, {bfgs.BestParameterSet[1]:F10})");
Console.WriteLine($"Function value: {bfgs.BestFitness:E10}");
Console.WriteLine($"Gradient evaluations: {bfgs.GradientEvaluations}");
```

**Advantages**: Superlinear convergence, efficient for smooth functions.
**Requirements**: Function must be differentiable (uses numerical gradients).

### Bounded Optimization

All local optimizers support box constraints:

```cs
double[] lower = { -5, -5 };
double[] upper = { 5, 5 };

var bfgs = new BFGS(rosenbrock, 2, initialGuess, lower, upper);
bfgs.Minimize();
```

---

## Global Optimization

### Differential Evolution (DE)

Evolutionary algorithm using vector differences [[4]](#ref4):

```cs
// Rastrigin function (many local minima)
Func<double[], double> rastrigin = x =>
{
    double sum = 10 * x.Length;
    for (int i = 0; i < x.Length; i++)
        sum += x[i] * x[i] - 10 * Math.Cos(2 * Math.PI * x[i]);
    return sum;
};

double[] lower = { -5.12, -5.12, -5.12, -5.12, -5.12 };
double[] upper = { 5.12, 5.12, 5.12, 5.12, 5.12 };

var de = new DifferentialEvolution(rastrigin, 5, lower, upper);
de.PopulationSize = 50;
de.MaxGenerations = 1000;
de.CrossoverRate = 0.9;
de.MutationFactor = 0.8;
de.Minimize();

Console.WriteLine($"Best fitness: {de.BestFitness:E6}");
Console.WriteLine($"Best solution: [{string.Join(", ", de.BestParameterSet.Select(x => x.ToString("F4")))}]");
```

**Parameters**:
- `PopulationSize`: Number of candidate solutions (typically 5-10× dimension)
- `CrossoverRate` (CR): Probability of crossover (0.5-1.0)
- `MutationFactor` (F): Differential weight (0.4-1.0)

### Shuffled Complex Evolution (SCE-UA)

Designed for hydrologic model calibration [[5]](#ref5):

```cs
var sce = new ShuffledComplexEvolution(objectiveFunction, nDimensions, lower, upper);
sce.NumberOfComplexes = 5;
sce.PointsPerComplex = 2 * nDimensions + 1;
sce.MaxIterations = 10000;
sce.Minimize();

Console.WriteLine($"Best fitness: {sce.BestFitness:F6}");
```

**Advantages**: Excellent for hydrologic calibration, balances exploration and exploitation.

### Particle Swarm Optimization (PSO)

Swarm intelligence algorithm [[6]](#ref6):

```cs
var pso = new ParticleSwarm(rastrigin, 5, lower, upper);
pso.SwarmSize = 50;
pso.MaxIterations = 1000;
pso.InertiaWeight = 0.7;
pso.CognitiveWeight = 1.5;
pso.SocialWeight = 1.5;
pso.Minimize();

Console.WriteLine($"Best fitness: {pso.BestFitness:E6}");
```

**Parameters**:
- `InertiaWeight` (ω): Momentum factor (0.4-0.9)
- `CognitiveWeight` (c₁): Personal best attraction
- `SocialWeight` (c₂): Global best attraction

### Simulated Annealing

Probabilistic method inspired by metallurgical annealing [[1]](#ref1):

```cs
var sa = new SimulatedAnnealing(rastrigin, 5, lower, upper);
sa.InitialTemperature = 100;
sa.CoolingRate = 0.95;
sa.MaxIterations = 10000;
sa.Minimize();

Console.WriteLine($"Best fitness: {sa.BestFitness:E6}");
```

**Advantages**: Can escape local minima, simple to implement.

---

## Constrained Optimization

### Augmented Lagrangian

Handles equality and inequality constraints [[3]](#ref3):

```cs
// Minimize f(x) = x₀² + x₁²
// Subject to: x₀ + x₁ = 1 (equality)
//             x₀ ≥ 0, x₁ ≥ 0 (bounds)

Func<double[], double> objective = x => x[0] * x[0] + x[1] * x[1];

// Equality constraint: g(x) = 0
Func<double[], double> equality = x => x[0] + x[1] - 1;

double[] initial = { 0.5, 0.5 };
double[] lower = { 0, 0 };
double[] upper = { double.PositiveInfinity, double.PositiveInfinity };

var al = new AugmentedLagrangian(objective, 2, initial, lower, upper);
al.AddEqualityConstraint(equality);
al.Minimize();

Console.WriteLine($"Solution: ({al.BestParameterSet[0]:F6}, {al.BestParameterSet[1]:F6})");
// Expected: (0.5, 0.5)
```

### Inequality Constraints

```cs
// g(x) ≤ 0 form
Func<double[], double> inequality = x => x[0] + x[1] - 2;  // x₀ + x₁ ≤ 2

al.AddInequalityConstraint(inequality);
```

---

## Multi-Start Optimization

Combine local search with multiple random starts:

```cs
int nStarts = 20;
double bestFitness = double.MaxValue;
double[] bestSolution = null;

var rng = new Random(12345);

for (int i = 0; i < nStarts; i++)
{
    // Random starting point
    double[] start = lower.Zip(upper, (lo, hi) => lo + rng.NextDouble() * (hi - lo)).ToArray();
    
    var bfgs = new BFGS(objective, nDimensions, start, lower, upper);
    bfgs.Minimize();
    
    if (bfgs.BestFitness < bestFitness)
    {
        bestFitness = bfgs.BestFitness;
        bestSolution = bfgs.BestParameterSet.ToArray();
    }
}

Console.WriteLine($"Best found: {bestFitness:E6}");
```

---

## Optimization with Progress Reporting

```cs
var de = new DifferentialEvolution(objective, nDimensions, lower, upper);

de.ProgressChanged += (sender, args) =>
{
    if (args.Generation % 100 == 0)
    {
        Console.WriteLine($"Generation {args.Generation}: Best = {args.BestFitness:E6}");
    }
};

de.Minimize();
```

---

## Parallel Optimization

Enable parallel function evaluations:

```cs
var de = new DifferentialEvolution(objective, nDimensions, lower, upper);
de.ParallelEvaluations = true;
de.MaxDegreeOfParallelism = Environment.ProcessorCount;
de.Minimize();
```

---

## Choosing an Optimizer

| Scenario | Recommended |
|----------|-------------|
| 1D optimization | Brent's Method |
| Smooth, low-dimensional (<10) | BFGS |
| Noisy or non-smooth | Nelder-Mead |
| High-dimensional, smooth | L-BFGS (if available) or DE |
| Many local minima | Differential Evolution |
| Hydrologic calibration | SCE-UA |
| Black-box, expensive | Bayesian optimization |
| Constrained | Augmented Lagrangian |

---

## Example: Distribution Fitting via MLE

```cs
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;

double[] data = { 12.5, 15.2, 11.8, 18.9, 14.2, 16.5, 13.4, 17.8 };

// Negative log-likelihood for Normal distribution
Func<double[], double> negLogLik = parameters =>
{
    double mu = parameters[0];
    double sigma = parameters[1];
    
    if (sigma <= 0) return double.MaxValue;
    
    var dist = new Normal(mu, sigma);
    return -dist.LogLikelihood(data);
};

double[] initial = { Statistics.Mean(data), Statistics.StandardDeviation(data) };
double[] lower = { -100, 0.001 };
double[] upper = { 100, 100 };

var bfgs = new BFGS(negLogLik, 2, initial, lower, upper);
bfgs.Minimize();

Console.WriteLine($"MLE estimates: μ = {bfgs.BestParameterSet[0]:F4}, σ = {bfgs.BestParameterSet[1]:F4}");
```

---

## Example: Hydrologic Model Calibration

```cs
// Objective: Minimize 1 - NSE
Func<double[], double> calibrationObjective = parameters =>
{
    // Run hydrologic model with parameters
    double[] simulated = RunModel(parameters);
    
    // Compute NSE
    double nse = GoodnessOfFit.NSE(observed, simulated);
    
    // Return 1 - NSE (minimization)
    return 1 - nse;
};

// Parameter bounds (e.g., for conceptual rainfall-runoff model)
double[] lower = { 0, 0, 0, 0.1, 1 };
double[] upper = { 500, 5, 1, 10, 100 };

var sce = new ShuffledComplexEvolution(calibrationObjective, 5, lower, upper);
sce.NumberOfComplexes = 5;
sce.MaxIterations = 5000;
sce.Tolerance = 1e-6;
sce.Minimize();

Console.WriteLine($"Calibrated NSE: {1 - sce.BestFitness:F4}");
Console.WriteLine($"Best parameters: [{string.Join(", ", sce.BestParameterSet.Select(p => p.ToString("F3")))}]");
```

---

## Convergence Diagnostics

```cs
var de = new DifferentialEvolution(objective, nDimensions, lower, upper);
de.Minimize();

// Check convergence
Console.WriteLine($"Converged: {de.HasConverged}");
Console.WriteLine($"Iterations: {de.Iterations}");
Console.WriteLine($"Function evaluations: {de.FunctionEvaluations}");
Console.WriteLine($"Final fitness: {de.BestFitness:E10}");

// Check if at bounds (may indicate need for wider bounds)
for (int i = 0; i < nDimensions; i++)
{
    if (Math.Abs(de.BestParameterSet[i] - lower[i]) < 1e-6 ||
        Math.Abs(de.BestParameterSet[i] - upper[i]) < 1e-6)
    {
        Console.WriteLine($"Warning: Parameter {i} at bound");
    }
}
```

---

## References

<a id="ref1">[1]</a> Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). *Numerical Recipes: The Art of Scientific Computing* (3rd ed.). Cambridge University Press.

<a id="ref2">[2]</a> Nelder, J. A., & Mead, R. (1965). A simplex method for function minimization. *The Computer Journal*, 7(4), 308-313.

<a id="ref3">[3]</a> Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.). Springer.

<a id="ref4">[4]</a> Storn, R., & Price, K. (1997). Differential evolution—A simple and efficient heuristic for global optimization over continuous spaces. *Journal of Global Optimization*, 11(4), 341-359.

<a id="ref5">[5]</a> Duan, Q., Sorooshian, S., & Gupta, V. K. (1994). Optimal use of the SCE-UA global optimization method for calibrating watershed models. *Journal of Hydrology*, 158(3-4), 265-284.

<a id="ref6">[6]</a> Kennedy, J., & Eberhart, R. (1995). Particle swarm optimization. *Proceedings of ICNN'95*, 4, 1942-1948.
