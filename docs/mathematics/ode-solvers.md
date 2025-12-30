# ODE Solvers

[← Previous: Special Functions](special-functions.md) | [Back to Index](../index.md) | [Next: Hypothesis Tests →](../statistics/hypothesis-tests.md)

The ***Numerics*** library provides Runge-Kutta methods for solving Ordinary Differential Equations (ODEs). These methods are essential for modeling dynamic systems, population dynamics, chemical reactions, and physical processes.

## Overview

An ordinary differential equation has the form:

```
dy/dt = f(t, y)
```

Given initial condition y(t₀) = y₀, we want to find y(t) for t > t₀.

The `RungeKutta` class provides several methods for numerical integration:
- Second-order Runge-Kutta (RK2)
- Fourth-order Runge-Kutta (RK4)
- Runge-Kutta-Fehlberg (adaptive step size)
- Cash-Karp (adaptive step size)

## Second-Order Runge-Kutta

Simple but less accurate method:

```cs
using Numerics.Mathematics.ODESolvers;

// Solve: dy/dt = y, y(0) = 1
// Analytical solution: y(t) = e^t

Func<double, double, double> f = (t, y) => y;

double y0 = 1.0;        // Initial condition
double t0 = 0.0;        // Start time
double t1 = 1.0;        // End time
int steps = 100;        // Time steps

double[] solution = RungeKutta.SecondOrder(f, y0, t0, t1, steps);

// Check solution at t=1
Console.WriteLine($"Numerical: y(1) = {solution[steps - 1]:F6}");
Console.WriteLine($"Analytical: e^1 = {Math.E:F6}");
Console.WriteLine($"Error: {Math.Abs(solution[steps - 1] - Math.E):E4}");
```

## Fourth-Order Runge-Kutta (RK4)

Most commonly used method - good balance of accuracy and efficiency:

### Array Output (Multiple Time Steps)

```cs
using Numerics.Mathematics.ODESolvers;

// Solve: dy/dt = -2t*y, y(0) = 1
// Analytical solution: y(t) = e^(-t²)

Func<double, double, double> f = (t, y) => -2.0 * t * y;

double y0 = 1.0;
double t0 = 0.0;
double t1 = 2.0;
int steps = 200;

double[] solution = RungeKutta.FourthOrder(f, y0, t0, t1, steps);

// Print solution at several points
double dt = (t1 - t0) / (steps - 1);
for (int i = 0; i < steps; i += 40)
{
    double t = t0 + i * dt;
    double numerical = solution[i];
    double analytical = Math.Exp(-t * t);
    
    Console.WriteLine($"t={t:F2}: y_num={numerical:F6}, y_exact={analytical:F6}, " +
                     $"error={Math.Abs(numerical - analytical):E4}");
}
```

### Single Step

Useful for step-by-step integration:

```cs
// Single RK4 step
double t = 0.0;
double y = 1.0;
double dt = 0.1;

double y_next = RungeKutta.FourthOrder(f, y, t, dt);

Console.WriteLine($"y({t + dt:F1}) = {y_next:F6}");
```

## Adaptive Step Size Methods

Adaptive methods automatically adjust step size to maintain accuracy while minimizing computation.

### Runge-Kutta-Fehlberg (RKF45)

Fifth-order accuracy with embedded fourth-order error estimate:

```cs
using Numerics.Mathematics.ODESolvers;

// Stiff ODE: dy/dt = -15y, y(0) = 1
Func<double, double, double> f = (t, y) => -15.0 * y;

double y0 = 1.0;
double t0 = 0.0;
double dt = 0.1;        // Initial step size
double dtMin = 1e-6;    // Minimum step size
double tolerance = 1e-6; // Error tolerance

double y_final = RungeKutta.Fehlberg(f, y0, t0, dt, dtMin, tolerance);

Console.WriteLine($"y({t0 + dt:F3}) = {y_final:F8}");
Console.WriteLine($"Analytical: {Math.Exp(-15.0 * (t0 + dt)):F8}");
```

### Cash-Karp Method

Alternative adaptive method with different error estimation:

```cs
double y_final_ck = RungeKutta.CashKarp(f, y0, t0, dt, dtMin, tolerance);

Console.WriteLine($"Cash-Karp: y({t0 + dt:F3}) = {y_final_ck:F8}");
```

**When to use adaptive methods:**
- Stiff equations
- Rapidly changing solutions
- Need to guarantee error tolerance
- Variable dynamics (slow then fast)

## Practical Examples

### Example 1: Exponential Decay (Radioactive Decay)

```cs
// Model: dN/dt = -λN, where λ is decay constant
// N(t) = N₀ e^(-λt)

double lambda = 0.693;  // Half-life ≈ 1 year (ln(2))
double N0 = 1000.0;     // Initial amount

Func<double, double, double> decay = (t, N) => -lambda * N;

double[] N = RungeKutta.FourthOrder(decay, N0, 0, 5.0, 100);

Console.WriteLine("Radioactive Decay:");
Console.WriteLine("Time | Amount | Half-Lives");
double dt = 5.0 / 99;
for (int i = 0; i < N.Length; i += 20)
{
    double t = i * dt;
    double halfLives = t / (Math.Log(2) / lambda);
    Console.WriteLine($"{t,4:F2} | {N[i],6:F1} | {halfLives,10:F3}");
}
```

### Example 2: Logistic Growth (Population Dynamics)

```cs
// Model: dP/dt = r*P*(1 - P/K)
// Where: r = growth rate, K = carrying capacity

double r = 0.5;      // Growth rate
double K = 1000.0;   // Carrying capacity
double P0 = 10.0;    // Initial population

Func<double, double, double> logistic = (t, P) => r * P * (1.0 - P / K);

double[] P = RungeKutta.FourthOrder(logistic, P0, 0, 20.0, 200);

Console.WriteLine("Logistic Growth:");
Console.WriteLine("Time | Population | % of Carrying Capacity");
double dt = 20.0 / 199;
for (int i = 0; i < P.Length; i += 40)
{
    double t = i * dt;
    double percent = 100.0 * P[i] / K;
    Console.WriteLine($"{t,4:F1} | {P[i],10:F1} | {percent,22:F1}%");
}

// Find time to reach 95% of carrying capacity
for (int i = 0; i < P.Length; i++)
{
    if (P[i] >= 0.95 * K)
    {
        Console.WriteLine($"\nReaches 95% capacity at t = {i * dt:F2}");
        break;
    }
}
```

### Example 3: Harmonic Oscillator

For second-order ODEs, convert to system of first-order ODEs:

```cs
// Second-order: d²y/dt² + ω²y = 0
// Convert to system:
//   dy₁/dt = y₂
//   dy₂/dt = -ω²y₁
// Where y₁ = position, y₂ = velocity

double omega = 2.0 * Math.PI;  // Angular frequency (1 Hz)

// For systems, use multiple passes
double y1_0 = 1.0;  // Initial position
double y2_0 = 0.0;  // Initial velocity

Func<double, double, double> dy1dt = (t, y1) => y2;  // Needs y2
Func<double, double, double> dy2dt = (t, y2) => -omega * omega * y1;  // Needs y1

// Need to track both variables manually
int steps = 1000;
double t0 = 0.0, t1 = 2.0;
double dt = (t1 - t0) / (steps - 1);

double[] y1 = new double[steps];  // Position
double[] y2 = new double[steps];  // Velocity

y1[0] = y1_0;
y2[0] = y2_0;

for (int i = 1; i < steps; i++)
{
    double t = t0 + (i - 1) * dt;
    
    // RK4 for system
    double k1_y1 = y2[i - 1];
    double k1_y2 = -omega * omega * y1[i - 1];
    
    double k2_y1 = y2[i - 1] + 0.5 * dt * k1_y2;
    double k2_y2 = -omega * omega * (y1[i - 1] + 0.5 * dt * k1_y1);
    
    double k3_y1 = y2[i - 1] + 0.5 * dt * k2_y2;
    double k3_y2 = -omega * omega * (y1[i - 1] + 0.5 * dt * k2_y1);
    
    double k4_y1 = y2[i - 1] + dt * k3_y2;
    double k4_y2 = -omega * omega * (y1[i - 1] + dt * k3_y1);
    
    y1[i] = y1[i - 1] + dt / 6.0 * (k1_y1 + 2 * k2_y1 + 2 * k3_y1 + k4_y1);
    y2[i] = y2[i - 1] + dt / 6.0 * (k1_y2 + 2 * k2_y2 + 2 * k3_y2 + k4_y2);
}

Console.WriteLine("Simple Harmonic Motion:");
Console.WriteLine("Time | Position | Velocity");
for (int i = 0; i < steps; i += 100)
{
    double t = t0 + i * dt;
    Console.WriteLine($"{t,4:F2} | {y1[i],8:F4} | {y2[i],8:F4}");
}
```

### Example 4: Predator-Prey (Lotka-Volterra)

Classic ecology model:

```cs
// dx/dt = αx - βxy  (prey)
// dy/dt = δxy - γy  (predator)

double alpha = 1.0;   // Prey growth rate
double beta = 0.1;    // Predation rate
double delta = 0.075; // Predator growth from predation
double gamma = 1.5;   // Predator death rate

double x0 = 10.0;     // Initial prey
double y0 = 5.0;      // Initial predator

int steps = 1000;
double dt = 0.01;
double t0 = 0.0;

double[] x = new double[steps];  // Prey
double[] y = new double[steps];  // Predator

x[0] = x0;
y[0] = y0;

for (int i = 1; i < steps; i++)
{
    double t = t0 + (i - 1) * dt;
    
    // RK4 for Lotka-Volterra system
    double k1_x = alpha * x[i - 1] - beta * x[i - 1] * y[i - 1];
    double k1_y = delta * x[i - 1] * y[i - 1] - gamma * y[i - 1];
    
    double k2_x = alpha * (x[i - 1] + 0.5 * dt * k1_x) - 
                  beta * (x[i - 1] + 0.5 * dt * k1_x) * (y[i - 1] + 0.5 * dt * k1_y);
    double k2_y = delta * (x[i - 1] + 0.5 * dt * k1_x) * (y[i - 1] + 0.5 * dt * k1_y) - 
                  gamma * (y[i - 1] + 0.5 * dt * k1_y);
    
    double k3_x = alpha * (x[i - 1] + 0.5 * dt * k2_x) - 
                  beta * (x[i - 1] + 0.5 * dt * k2_x) * (y[i - 1] + 0.5 * dt * k2_y);
    double k3_y = delta * (x[i - 1] + 0.5 * dt * k2_x) * (y[i - 1] + 0.5 * dt * k2_y) - 
                  gamma * (y[i - 1] + 0.5 * dt * k2_y);
    
    double k4_x = alpha * (x[i - 1] + dt * k3_x) - 
                  beta * (x[i - 1] + dt * k3_x) * (y[i - 1] + dt * k3_y);
    double k4_y = delta * (x[i - 1] + dt * k3_x) * (y[i - 1] + dt * k3_y) - 
                  gamma * (y[i - 1] + dt * k3_y);
    
    x[i] = x[i - 1] + dt / 6.0 * (k1_x + 2 * k2_x + 2 * k3_x + k4_x);
    y[i] = y[i - 1] + dt / 6.0 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
}

Console.WriteLine("Predator-Prey Dynamics:");
Console.WriteLine("Time  | Prey | Predator");
for (int i = 0; i < steps; i += 100)
{
    double t = t0 + i * dt;
    Console.WriteLine($"{t,5:F2} | {x[i],4:F1} | {y[i],8:F1}");
}
```

## Choosing a Method

| Method | Order | When to Use |
|--------|-------|-------------|
| **Second-Order** | O(h²) | Simple problems, rough estimates |
| **Fourth-Order** | O(h⁴) | General purpose, good accuracy |
| **Fehlberg** | O(h⁵) | Adaptive control, stiff equations |
| **Cash-Karp** | O(h⁵) | Alternative adaptive method |

**Step Size Guidelines:**
- Fixed step: Choose dt such that solution is smooth
- Adaptive: Set tolerance based on required accuracy
- Stiff equations: Use smaller steps or adaptive methods
- For stability: dt < 2/|λ_max| where λ_max is largest eigenvalue

## Best Practices

1. **Verify with known solutions** when possible
2. **Check convergence** by halving step size
3. **Use adaptive methods** for variable dynamics
4. **Monitor conservation** (energy, mass) if applicable
5. **Plot solutions** to detect instabilities
6. **Consider stability** for stiff equations

## Limitations

- Fixed-step methods require careful step size selection
- Stiff equations may require specialized solvers
- Systems require manual coupling of equations
- No built-in event detection
- For very high accuracy, consider specialized ODE libraries

---

[← Previous: Special Functions](special-functions.md) | [Back to Index](../index.md) | [Next: Hypothesis Tests →](../statistics/hypothesis-tests.md)
