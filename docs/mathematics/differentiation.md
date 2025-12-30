# Numerical Differentiation

[← Previous: Numerical Integration](integration.md) | [Back to Index](../index.md) | [Next: Optimization →](optimization.md)

Numerical differentiation is a fundamental technique in various scientific and engineering fields. Many optimization algorithms, such as gradient descent and Newton's method, rely on calculating the gradient of a function to find its minimum or maximum points. Some optimization methods, like quasi-Newton methods, use the Hessian matrix (the matrix of second-order derivatives) to determine the curvature of the function. Numerical differentiation provides a way to approximate both the gradient and the Hessian matrix efficiently.

## Derivatives

In ***Numerics***, the derivative is evaluated using the two-point (central difference) formula by default:

```math
\frac{df}{dx} = \frac{f(x + h) - f(x - h)}{2h}
```

where $x$ is the input point and $h$ represents a small change in $x$. In ***Numerics***, the step size $h$ is automatically determined according to the magnitude of $x$:

```math
\begin{equation}
  h =
    \begin{cases}
      \mid x \mid \cdot \epsilon^\frac{1}{2}  & x \neq 0\\
      \epsilon^\frac{1}{2}  & x = 0\\
    \end{cases}       
\end{equation}
```

where $\epsilon$ is double precision machine epsilon. The step size $h$ can also be user-defined.

For example, consider the simple function:

```math
f(x)=x^3
```

Differentiating with respect to $x$ gives:

```math
\frac{df}{dx} =3x^2
```

Evaluating the function at $x=2$ yields a derivative equal to 12:

```math
\frac{df}{dx} =3\cdot2^2 = 12
```

Now, let's implement this in ***Numerics***. First, we need to reference _Numerics_ and the _Mathematics_ namespace:

```cs
using Numerics.Mathematics;
```

Next, create the test function:

```cs
/// <summary>
/// Test function: f(x) = x^3
/// </summary>
public double FX(double x)
{
    return Math.Pow(x, 3);
}
```

And then compute the derivative using the _NumericalDerivative_ class and the _Derivative_ method, which uses the two-point (central difference) formula:

```cs
double dFdx = NumericalDerivative.Derivative(FX, 2); // 11.999999949167176
```

The accuracy of numerical differentiation depends on the smoothness of the function and the step size. For functions with discontinuities or sharp gradients, the two-point formula might not accurately capture the derivative at those points. If the function changes rapidly, a larger step size might miss important details.

### Alternative Finite Difference Formulas

The ***Numerics*** library provides additional finite difference methods:

#### Forward Difference

The forward difference approximation uses:

```math
\frac{df}{dx} \approx \frac{f(x + h) - f(x)}{h}
```

```cs
double dFdx = NumericalDerivative.ForwardDifference(FX, 2); // 12.000001490104875
```

This method is useful when you can only evaluate the function ahead of the point, but it's less accurate than the central difference method.

#### Backward Difference

The backward difference approximation uses:

```math
\frac{df}{dx} \approx \frac{f(x) - f(x - h)}{h}
```

```cs
double dFdx = NumericalDerivative.BackwardDifference(FX, 2); // 11.999998508229537
```

This method is useful when you can only evaluate the function behind the point.

#### Central Difference

The central difference method (which is also used by the `Derivative` method) provides better accuracy:

```cs
double dFdx = NumericalDerivative.CentralDifference(FX, 2); // 11.999999949167176
```

### Ridder's Method

The ***Numerics*** library also provides Ridder's method for computing the numerical derivative, which can be more accurate than the two-point method in some cases. This method also outputs an estimate of the error in the derivative:

```cs
double dFdx = NumericalDerivative.RiddersMethod(FX, 2, out var err); // 11.99999994392223
Console.WriteLine($"Derivative: {dFdx}");
Console.WriteLine($"Estimated error: {err}");
```

Ridder's method uses Richardson extrapolation to refine the estimate by evaluating the derivative at multiple step sizes and extrapolating to zero step size.

### Custom Step Size

You can specify a custom step size if the automatic determination is not suitable for your problem:

```cs
double h = 0.001; // Custom step size
double dFdx = NumericalDerivative.Derivative(FX, 2, h);
```

## Second Derivatives

The ***Numerics*** library provides methods for computing second derivatives. The second derivative measures the rate of change of the first derivative, or the curvature of the function.

### Central Second Derivative

The central difference approximation for the second derivative is:

```math
\frac{d^2f}{dx^2} \approx \frac{f(x + h) - 2f(x) + f(x - h)}{h^2}
```

```cs
double d2Fdx2 = NumericalDerivative.SecondDerivative(FX, 2); // 11.999997613071552
```

For our test function $f(x) = x^3$, the second derivative is $\frac{d^2f}{dx^2} = 6x$, so at $x=2$, we expect 12.

### Forward Second Derivative

```cs
double d2Fdx2 = NumericalDerivative.SecondDerivativeForward(FX, 2);
```

### Backward Second Derivative

```cs
double d2Fdx2 = NumericalDerivative.SecondDerivativeBackward(FX, 2);
```

## Gradient

The gradient is a vector of first-order partial derivatives of a scalar-valued function:

```math
\nabla f = \left(\frac{\partial f}{\partial x_1}, \frac{\partial f}{\partial x_2}, \ldots , \frac{\partial f}{\partial x_n}\right)
```

For example, consider a function with three variables:

```math
f(x,y,z)=x^3+y^4+z^5
```

Differentiating with respect to each variable gives:

```math
\frac{\partial f}{\partial x} =3x^2  \quad \frac{\partial f}{\partial y} =4y^3 \quad \frac{\partial f}{\partial z} =5z^4
```

Evaluating the function at $x=2$, $y=2$, and $z=2$ yields partial derivatives equal to $12$, $32$, and $80$, respectively:

```math
\frac{\partial f}{\partial x} =3\cdot2^2=12  \quad \frac{\partial f}{\partial y} =4\cdot2^3=32 \quad \frac{\partial f}{\partial z} =5\cdot2^4=80
```

The gradient is the vector of these first-order partial derivatives:

```math
\nabla f = \{ 12, 32, 80\}
```

In ***Numerics***, the gradient is computed using the two-point formula described earlier. Let's create the test function:

```cs
/// <summary>
/// Test function: f(x, y, z) = x^3 + y^4 + z^5
/// </summary>
public double FXYZ(double[] x)
{
    return Math.Pow(x[0], 3) + Math.Pow(x[1], 4) + Math.Pow(x[2], 5);
}
```

And then compute the gradient:

```cs
double[] gradient = NumericalDerivative.Gradient(FXYZ, new double[] {2, 2, 2}); 
// {12.000000019411923, 32.000000014301264, 79.999999754774166}

Console.WriteLine($"∂f/∂x = {gradient[0]:F6}");
Console.WriteLine($"∂f/∂y = {gradient[1]:F6}");
Console.WriteLine($"∂f/∂z = {gradient[2]:F6}");
```

## Hessian

The Hessian is a square matrix of second-order partial derivatives of a scalar-valued function:

```math
\mathbf{H}_{i,j} = \frac{\partial^2 f}{\partial x_i\partial x_j}
```

Using the same 3-variable function as before $f(x,y,z)=x^3+y^4+z^5$, the Hessian matrix becomes:

```math
\mathbf{H} =
  \left[ {\begin{array}{ccc}
   \frac{\partial^2 f}{\partial x\partial x} & \frac{\partial^2 f}{\partial x\partial y} & \frac{\partial^2 f}{\partial x\partial z} \\
   \frac{\partial^2 f}{\partial y\partial x} & \frac{\partial^2 f}{\partial y\partial y} & \frac{\partial^2 f}{\partial y\partial z} \\
   \frac{\partial^2 f}{\partial z\partial x} & \frac{\partial^2 f}{\partial z\partial y} & \frac{\partial^2 f}{\partial z\partial z} \\
  \end{array} } \right]
```

Since none of the variables interact with each other, the Hessian reduces to:

```math
\mathbf{H} =
  \left[ {\begin{array}{ccc}
   \frac{\partial^2 f}{\partial x^2} & 0 & 0 \\
   0 & \frac{\partial^2 f}{\partial y^2} & 0 \\
   0 & 0 & \frac{\partial^2 f}{\partial z^2} \\
  \end{array} } \right]
```

Taking the second-order derivatives with respect to each variable gives:

```math
\frac{\partial^2 f}{\partial x^2} =6x  \quad \frac{\partial^2 f}{\partial y^2} =12y^2 \quad \frac{\partial^2 f}{\partial z^2} =20z^3
```

Evaluating the function at $x=2$, $y=2$, and $z=2$ yields:

```math
\frac{\partial^2 f}{\partial x^2} =6\cdot2=12  \quad \frac{\partial^2 f}{\partial y^2} =12\cdot2^2=48 \quad \frac{\partial^2 f}{\partial z^2} =20\cdot2^3=160
```

Now, to compute the Hessian in ***Numerics***, simply do the following:

```cs
double[,] hessian = NumericalDerivative.Hessian(FXYZ, new double[] {2, 2, 2}); 
// [0,0] = 12.000009765258449
// [0,1] = 8.5443864603330376E-06
// [0,2] = 0
// [1,0] = 8.5443864603330376E-06
// [1,1] = 48.000004883487954
// [1,2] = 0
// [2,0] = 0
// [2,1] = 0
// [2,2] = 159.99999349326262

Console.WriteLine("Hessian matrix:");
for (int i = 0; i < 3; i++)
{
    for (int j = 0; j < 3; j++)
    {
        Console.Write($"{hessian[i, j],12:F6}  ");
    }
    Console.WriteLine();
}
```

The small off-diagonal elements (on the order of $10^{-6}$) are numerical errors and should be zero for this particular function.

### Example: Function with Interacting Variables

Now, consider another example where the function variables interact:

```math
f(x,y)=x^3-2xy-y^6
```

The Hessian matrix is:

```math
\mathbf{H} =
  \left[ {\begin{array}{cc}
   \frac{\partial^2 f}{\partial x\partial x} & \frac{\partial^2 f}{\partial x \partial y}  \\
   \frac{\partial^2 f}{\partial y\partial x} & \frac{\partial^2 f}{\partial y \partial y} \\
  \end{array} } \right]
```

Taking the second-order derivatives with respect to each variable gives:

```math
\frac{\partial^2 f}{\partial x^2} =6x  \quad \frac{\partial^2 f}{\partial x \partial y} = \frac{\partial^2 f}{\partial y \partial x} = -2 \quad \frac{\partial^2 f}{\partial y^2} = -30y^4
```

Evaluating the function at $x=1$ and $y=2$ yields:

```math
\frac{\partial^2 f}{\partial x^2} =6\cdot1=6  \quad \frac{\partial^2 f}{\partial x \partial y} = -2 \quad \frac{\partial^2 f}{\partial y^2} =-30\cdot2^4=-480
```

Create the test function:

```cs
/// <summary>
/// Test function: f(x, y) = x^3 - 2xy - y^6
/// </summary>
public double FXY(double[] x)
{
    return Math.Pow(x[0], 3) - 2 * x[0] * x[1] - Math.Pow(x[1], 6);
}
```

Now, compute the Hessian using ***Numerics***:

```cs
double[,] hessian = NumericalDerivative.Hessian(FXY, new double[] {1, 2}); 
// [0,0] = 5.9999414101667655
// [0,1] = -2.000001627543075
// [1,0] = -2.000001627543075
// [1,1] = -480.00004883487958

Console.WriteLine("Hessian matrix:");
Console.WriteLine($"[{hessian[0,0],10:F4}  {hessian[0,1],10:F4}]");
Console.WriteLine($"[{hessian[1,0],10:F4}  {hessian[1,1],10:F4}]");
```

Note that the Hessian is symmetric, as expected from the equality of mixed partial derivatives (Schwarz's theorem), and the off-diagonal elements correctly capture the interaction between variables.

## Jacobian

The Jacobian matrix represents the first-order partial derivatives of a vector-valued function. For a function $\mathbf{f}: \mathbb{R}^n \rightarrow \mathbb{R}^m$, the Jacobian is an $m \times n$ matrix:

```math
\mathbf{J}_{i,j} = \frac{\partial f_i}{\partial x_j}
```

The ***Numerics*** library provides two overloads for computing the Jacobian:

```cs
// For a single output function f(x_i, x[])
double[,] jacobian1 = NumericalDerivative.Jacobian(
    (xi, x) => /* function of xi and x[] */,
    xValues,
    point
);

// For a vector-valued function f(x[]) -> y[]
double[,] jacobian2 = NumericalDerivative.Jacobian(
    x => /* returns double[] */,
    point
);
```

Example of computing a Jacobian for a system of equations:

```cs
// System: f1(x,y) = x^2 + y^2, f2(x,y) = xy
double[] F(double[] vars)
{
    double x = vars[0];
    double y = vars[1];
    return new double[] 
    { 
        x * x + y * y,  // f1
        x * y            // f2
    };
}

var point = new double[] { 2, 3 };
double[,] jacobian = NumericalDerivative.Jacobian(F, point);

// Jacobian at (2,3):
// [∂f1/∂x  ∂f1/∂y]   [2x    2y]   [4   6]
// [∂f2/∂x  ∂f2/∂y] = [ y     x] = [3   2]
```

## Calculating Step Size

The `CalculateStepSize` method computes an appropriate step size for numerical differentiation based on the magnitude of the point and the order of the derivative:

```cs
double h = NumericalDerivative.CalculateStepSize(x: 2.0, order: 1); 
// Returns approximately 1.49e-08 for first derivative

double h2 = NumericalDerivative.CalculateStepSize(x: 2.0, order: 2); 
// Returns approximately 3.45e-06 for second derivative
```

The step size is calculated as:

```math
h = |x| \cdot \epsilon^{1/(1+\text{order})}
```

where $\epsilon$ is machine epsilon. For $x=0$, the formula simplifies to $h = \epsilon^{1/(1+\text{order})}$.

## Best Practices

1. **Use Central Differences**: When possible, use the central difference method (default `Derivative` method) as it provides better accuracy than forward or backward differences.

2. **Ridder's Method for Critical Applications**: When you need both a derivative and an error estimate, use Ridder's method.

3. **Automatic Step Sizing**: The default automatic step sizing works well for most problems. Only specify a custom step size if you have specific numerical issues.

4. **Beware of Noise**: Numerical differentiation amplifies noise in function evaluations. If your function has numerical noise (e.g., from Monte Carlo simulations), consider smoothing or using a larger step size.

5. **Check for Symmetry**: For Hessian matrices, check that the result is symmetric (within numerical tolerance). Significant asymmetry indicates numerical issues.

6. **Scale Considerations**: For problems with variables at very different scales, consider normalizing variables before computing derivatives.

## Accuracy Considerations

The central difference formula has truncation error $O(h^2)$ and roundoff error $O(\epsilon/h)$, where $\epsilon$ is machine epsilon. The optimal step size balances these errors at approximately $h \approx \epsilon^{1/3}$ for first derivatives and $h \approx \epsilon^{1/4}$ for second derivatives. The automatic step sizing in ***Numerics*** uses $h \approx \epsilon^{1/2}$, which is a conservative choice that works well in practice.

For the second derivative, the truncation error is $O(h^2)$ and roundoff error is $O(\epsilon/h^2)$, making it more sensitive to numerical noise than first derivatives.

---

## References

The numerical differentiation methods implemented in ***Numerics*** are based on standard finite difference formulas well-documented in numerical analysis literature [[1]](#1). Ridder's method for derivative estimation with error bounds was introduced by Ridders [[2]](#2).

<a id="1">[1]</a> W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery, *Numerical Recipes: The Art of Scientific Computing*, 3rd ed., Cambridge, UK: Cambridge University Press, 2007.

<a id="2">[2]</a> C. J. F. Ridders, "Accurate computation of F'(x) and F'(x) F''(x)," *Advances in Engineering Software*, vol. 4, no. 2, pp. 75-76, 1982.

---

[← Previous: Numerical Integration](integration.md) | [Back to Index](../index.md) | [Next: Optimization →](optimization.md)
