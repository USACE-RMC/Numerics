# ODE Solvers

[← Previous: Special Functions](special-functions.md) | [Back to Index](../index.md) | [Next: Interpolation →](../data/interpolation.md)

An ordinary differential equation (ODE) relates a function to its derivatives. The general **initial value problem** has the form:

```math
\frac{dy}{dt} = f(t, y), \quad y(t_0) = y_0
```

Given the initial condition $y(t_0) = y_0$, we seek the solution $y(t)$ for $t > t_0$. Most ODEs arising in practice — population dynamics, chemical kinetics, mechanical systems, heat transfer — cannot be solved analytically and require numerical methods.

The ***Numerics*** library provides Runge-Kutta methods for solving ODEs numerically. These methods advance the solution from $y_n$ at time $t_n$ to $y_{n+1}$ at time $t_{n+1} = t_n + h$ by evaluating $f(t, y)$ at several intermediate points within the step.

## Second-Order Runge-Kutta (RK2)

The simplest Runge-Kutta method evaluates the derivative at the beginning and end of the step, then averages:

```math
k_1 = f(t_n, y_n)
```
```math
k_2 = f(t_n + h, y_n + h \cdot k_1)
```
```math
y_{n+1} = y_n + \frac{h}{2}(k_1 + k_2)
```

This is equivalent to the explicit trapezoidal rule. It is second-order accurate, meaning the local truncation error per step is $O(h^3)$ and the global accumulated error is $O(h^2)$.

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

The classical RK4 method is the workhorse of ODE solving — it provides an excellent balance of accuracy and efficiency. The method evaluates the derivative at four points within each step:

```math
k_1 = f(t_n, y_n)
```
```math
k_2 = f\!\left(t_n + \frac{h}{2}, \; y_n + \frac{h}{2} k_1\right)
```
```math
k_3 = f\!\left(t_n + \frac{h}{2}, \; y_n + \frac{h}{2} k_2\right)
```
```math
k_4 = f(t_n + h, \; y_n + h \cdot k_3)
```
```math
y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
```

The $k$-values can be interpreted as slope estimates: $k_1$ at the beginning, $k_2$ and $k_3$ at the midpoint (using different estimates to get there), and $k_4$ at the end. The weighted average gives $k_2$ and $k_3$ twice the weight of $k_1$ and $k_4$, similar to Simpson's rule in quadrature.

### Butcher Tableau

The coefficients of any Runge-Kutta method can be compactly represented in a **Butcher tableau**. For the classical RK4:

```
  0   |
  1/2 | 1/2
  1/2 |  0   1/2
  1   |  0    0    1
  ----|--------------------
      | 1/6  1/3  1/3  1/6
```

The left column gives the time offsets $c_i$ (where the derivative is evaluated), the body of the table gives the coefficients $a_{ij}$ for combining previous $k$-values, and the bottom row gives the weights $b_i$ for the final combination. The local truncation error is $O(h^5)$ and the global error is $O(h^4)$, making it fourth-order accurate.

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

Useful for step-by-step integration or for building custom solvers:

```cs
// Single RK4 step
double t = 0.0;
double y = 1.0;
double dt = 0.1;

double y_next = RungeKutta.FourthOrder(f, y, t, dt);

Console.WriteLine($"y({t + dt:F1}) = {y_next:F6}");
```

## Adaptive Step Size Methods

Fixed-step methods require the user to choose $h$ in advance, but the optimal step size depends on the local behavior of the solution. If $h$ is too large, the solution may be inaccurate or unstable; if too small, computation is wasted. **Adaptive methods** adjust $h$ automatically by estimating the local error at each step.

### Embedded Runge-Kutta Methods

The key idea behind embedded methods is to compute **two solutions of different orders** using mostly the same function evaluations. The difference between the two solutions provides an error estimate without extra cost.

Both the Fehlberg and Cash-Karp methods compute a fourth-order solution $y^{(4)}$ and a fifth-order solution $y^{(5)}$ using six stages ($k_1$ through $k_6$). The error estimate is:

```math
\text{err} = |y^{(4)} - y^{(5)}|
```

The step size is then adjusted:
- If $\text{err} \leq \text{tolerance}$: **accept** the step. If the error is very small ($< 0.01 \times \text{tolerance}$), double the step size
- If $\text{err} > \text{tolerance}$: **reject** the step, halve $h$, and retry
- If $h$ falls below the minimum step size `dtMin`, accept the step regardless (to avoid infinite loops)

### Runge-Kutta-Fehlberg (RKF45)

The Fehlberg method [[1]](#1) uses specific coefficients optimized to minimize the leading error term. Its Butcher tableau has the structure:

```
    0     |
   1/4    | 1/4
   3/8    | 3/32        9/32
  12/13   | 1932/2197  -7200/2197   7296/2197
    1     | 439/216    -8           3680/513    -845/4104
   1/2    | -8/27       2          -3544/2565    1859/4104   -11/40
  --------|------------------------------------------------------------
  RK4     | 25/216      0           1408/2565    2197/4104   -1/5       0
  RK5     | 16/135      0           6656/12825   28561/56430 -9/50      2/55
```

The fourth-order solution (row labeled RK4) is used for stepping, while the difference from the fifth-order solution provides the error estimate.

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

The Cash-Karp method [[2]](#2) uses a different set of coefficients optimized for a different balance of accuracy and stability:

```cs
double y_final_ck = RungeKutta.CashKarp(f, y0, t0, dt, dtMin, tolerance);

Console.WriteLine($"Cash-Karp: y({t0 + dt:F3}) = {y_final_ck:F8}");
```

**When to use adaptive methods:**
- When the solution dynamics change significantly over the integration interval (e.g., rapid transients followed by slow decay)
- When you need guaranteed error tolerance without manual tuning
- When the problem timescale is unknown in advance

## Stiffness

A differential equation is **stiff** when it contains dynamics on vastly different timescales — some components of the solution change rapidly while others change slowly. The classic example is a system with eigenvalues $\lambda_1 = -1$ and $\lambda_2 = -1000$: the fast component $e^{-1000t}$ decays in milliseconds, but the explicit RK4 step size is constrained by this fast component even long after it has decayed to zero.

For the linear test equation $\frac{dy}{dt} = \lambda y$ (with $\lambda < 0$), the explicit RK4 method is stable only when:

```math
|h \cdot \lambda| \leq 2.785
```

This means that for $\lambda = -1000$, the step size must be smaller than $h < 0.002785$ — regardless of the accuracy requirement. For stiff problems, this stability restriction forces an impractically small step size.

**Mitigation strategies:**
- Use the adaptive methods (Fehlberg, Cash-Karp), which will automatically reduce the step size in stiff regions
- For very stiff problems, consider implicit methods (not currently in the library) which have much larger stability regions
- If the stiff component has decayed, reformulate the problem to remove it

## Error Accumulation

The **local truncation error** is the error introduced in a single step assuming the previous value is exact. For RK4, this is $O(h^5)$ per step. The **global error** after integrating over a fixed interval $[t_0, T]$ with $N = (T-t_0)/h$ steps accumulates to $O(h^4)$, because the $N = O(1/h)$ steps each contribute $O(h^5)$ error.

However, error accumulation depends on the stability of the ODE. For **stable** equations ($\text{Re}(\lambda) < 0$), errors are damped and the global error remains $O(h^4)$. For **unstable** equations ($\text{Re}(\lambda) > 0$), errors grow exponentially, and even small per-step errors can lead to wildly inaccurate solutions over long integration intervals. For such problems, reducing the step size and monitoring the solution for physical plausibility is essential.

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
// Print final time point (the loop step skips the last element)
{
    double t = (N.Length - 1) * dt;
    double halfLives = t / (Math.Log(2) / lambda);
    Console.WriteLine($"{t,4:F2} | {N[N.Length - 1],6:F1} | {halfLives,10:F3}");
}
```

### Example 2: Logistic Growth (Population Dynamics)

The logistic equation models population growth with a carrying capacity:

```math
\frac{dP}{dt} = rP\left(1 - \frac{P}{K}\right)
```

The analytical solution is $P(t) = \frac{K}{1 + \left(\frac{K}{P_0} - 1\right)e^{-rt}}$, which has a sigmoidal shape. The growth rate $r$ controls how quickly the population approaches carrying capacity $K$.

```cs
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

### Example 3: Harmonic Oscillator (Systems of ODEs)

> **Note:** The `RungeKutta` class handles scalar ODEs only (i.e., a single dependent variable). For systems of ODEs (multiple coupled equations), you must implement the RK4 stepping logic manually, as shown below.

A second-order ODE can always be converted to a system of first-order ODEs. For the harmonic oscillator $\frac{d^2y}{dt^2} + \omega^2 y = 0$, we introduce $y_1 = y$ (position) and $y_2 = \frac{dy}{dt}$ (velocity):

```math
\frac{dy_1}{dt} = y_2, \qquad \frac{dy_2}{dt} = -\omega^2 y_1
```

The RK4 method is applied to both equations simultaneously, using the same $k$-values structure but evaluated for each component:

```cs
double omega = 2.0 * Math.PI;  // Angular frequency (1 Hz)

double y1_0 = 1.0;  // Initial position
double y2_0 = 0.0;  // Initial velocity

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

The Lotka-Volterra equations model the dynamics of two interacting species:

```math
\frac{dx}{dt} = \alpha x - \beta xy \qquad \text{(prey)}
```
```math
\frac{dy}{dt} = \delta xy - \gamma y \qquad \text{(predator)}
```

The prey population grows exponentially in the absence of predators (rate $\alpha$) but is reduced by predation (rate $\beta xy$). Predators grow when prey is abundant (rate $\delta xy$) but die in its absence (rate $\gamma y$). The system exhibits periodic oscillations — predator peaks follow prey peaks with a phase lag.

```cs
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

| Method | Order | Error per Step | When to Use |
|--------|-------|----------------|-------------|
| **Second-Order** | $O(h^2)$ | $O(h^3)$ | Simple problems, rough estimates |
| **Fourth-Order** | $O(h^4)$ | $O(h^5)$ | General purpose, good accuracy-to-cost ratio |
| **Fehlberg** | $O(h^5)$ | Embedded error | Adaptive control, variable dynamics |
| **Cash-Karp** | $O(h^5)$ | Embedded error | Alternative adaptive method |

**Step Size Guidelines:**
- Fixed step: Choose $h$ such that the solution is smooth and verify by halving $h$
- Adaptive: Set tolerance based on required accuracy — the method handles the rest
- For stability: $h < 2.785 / |\lambda_{\max}|$ where $\lambda_{\max}$ is the largest eigenvalue magnitude of the Jacobian

## Best Practices

1. **Verify with known solutions** when possible — compare against analytical results for test problems
2. **Check convergence** by halving the step size and confirming the solution doesn't change significantly (the error should decrease by $2^4 = 16 \times$ for RK4)
3. **Use adaptive methods** for problems with variable dynamics — they automatically concentrate effort where the solution changes rapidly
4. **Monitor conservation** — for Hamiltonian systems (e.g., harmonic oscillator), check that energy is conserved over long integrations
5. **Plot solutions** to detect instabilities — oscillations growing in amplitude usually indicate the step size is too large
6. **Consider reformulation** for stiff equations — if possible, analytically eliminate the fast component or use a change of variables

## Limitations

- Fixed-step methods require careful step size selection
- Very stiff equations may require specialized implicit solvers (not currently in the library)
- Systems of ODEs require manual coupling of equations — the `RungeKutta` class handles scalar ODEs only
- No built-in event detection (e.g., finding when the solution crosses zero)
- For very high accuracy or long-time integrations, consider symplectic integrators for Hamiltonian systems

---

## References

<a id="1">[1]</a> E. Fehlberg, "Low-order classical Runge-Kutta formulas with stepsize control and their application to some heat transfer problems," NASA Technical Report 315, 1969.

<a id="2">[2]</a> J. R. Cash and A. H. Karp, "A variable order Runge-Kutta method for initial value problems with rapidly varying right-hand sides," *ACM Transactions on Mathematical Software*, vol. 16, no. 3, pp. 201-222, 1990.

<a id="3">[3]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="4">[4]</a> J. C. Butcher, *Numerical Methods for Ordinary Differential Equations*, 3rd ed., Chichester: Wiley, 2016.

---

[← Previous: Special Functions](special-functions.md) | [Back to Index](../index.md) | [Next: Interpolation →](../data/interpolation.md)
