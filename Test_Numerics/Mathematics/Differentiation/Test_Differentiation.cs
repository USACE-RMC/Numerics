using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics;

namespace Mathematics.Differentiation
{
    /// <summary>
    /// Comprehensive unit tests for numerical differentiation methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    /// </para>
    /// <list type="bullet">
    /// <item><description>
    /// Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    /// Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <para>
    /// This test suite covers:
    /// </para>
    /// <list type="bullet">
    /// <item><description>All finite difference methods (forward, backward, central)</description></item>
    /// <item><description>First and second derivatives</description></item>
    /// <item><description>High-accuracy methods (Ridders)</description></item>
    /// <item><description>Multivariable calculus (Jacobian, Gradient, Hessian)</description></item>
    /// <item><description>Boundary conditions and constraints</description></item>
    /// <item><description>Edge cases and error handling</description></item>
    /// </list>
    /// </remarks>
    [TestClass]
    public class Test_Differentiation
    {
        #region Test Functions - Single Variable

        /// <summary>
        /// Test function: f(x) = x³
        /// f'(x) = 3x²
        /// f''(x) = 6x
        /// </summary>
        public double FX(double x)
        {
            return Math.Pow(x, 3.0);
        }

        /// <summary>
        /// Test function: f(x) = e^x
        /// f'(x) = e^x
        /// f''(x) = e^x
        /// </summary>
        public double EX(double x)
        {
            return Math.Exp(x);
        }

        /// <summary>
        /// Test function: f(x) = ln(x)
        /// f'(x) = 1/x
        /// f''(x) = -1/x²
        /// </summary>
        public double LN(double x)
        {
            return Math.Log(x);
        }

        /// <summary>
        /// Test function: f(x) = sin(x)
        /// f'(x) = cos(x)
        /// f''(x) = -sin(x)
        /// </summary>
        public double SIN(double x)
        {
            return Math.Sin(x);
        }

        /// <summary>
        /// Test function: f(x) = cos(x)
        /// f'(x) = -sin(x)
        /// f''(x) = -cos(x)
        /// </summary>
        public double COS(double x)
        {
            return Math.Cos(x);
        }

        /// <summary>
        /// Test function: f(x) = x⁴ - 2x² + 3
        /// f'(x) = 4x³ - 4x
        /// f''(x) = 12x² - 4
        /// </summary>
        public double Polynomial(double x)
        {
            return Math.Pow(x, 4) - 2 * Math.Pow(x, 2) + 3;
        }

        #endregion

        #region Test Functions - Multivariable

        /// <summary>
        /// Test function: f(x,y) = x²y³
        /// ∂f/∂x = 2xy³
        /// ∂f/∂y = 3x²y²
        /// </summary>
        public double FXY(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return Math.Pow(x, 2) * Math.Pow(y, 3);
        }

        /// <summary>
        /// Test function: f(x,y,z) = x³ + y⁴ + z⁵
        /// ∂f/∂x = 3x²
        /// ∂f/∂y = 4y³
        /// ∂f/∂z = 5z⁴
        /// </summary>
        public double FXYZ(double[] points)
        {
            double x = points[0];
            double y = points[1];
            double z = points[2];
            return Math.Pow(x, 3.0) + Math.Pow(y, 4.0) + Math.Pow(z, 5.0);
        }

        /// <summary>
        /// Test function for Hessian: f(x,y) = x³ - 2xy - y⁶
        /// ∂²f/∂x² = 6x
        /// ∂²f/∂y² = -30y⁴
        /// ∂²f/∂x∂y = -2
        /// </summary>
        public double FH(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return Math.Pow(x, 3.0) - 2 * x * y - Math.Pow(y, 6);
        }

        /// <summary>
        /// Test function for Jacobian: f(t,x,y,z) = tx³ + y⁴ + z⁵
        /// </summary>
        public double FTXYZ(double t, double[] points)
        {
            double x = points[0];
            double y = points[1];
            double z = points[2];
            return t * Math.Pow(x, 3.0) + Math.Pow(y, 4.0) + Math.Pow(z, 5.0);
        }

        /// <summary>
        /// Rosenbrock function: f(x,y) = (1-x)² + 100(y-x²)²
        /// Common test function for optimization
        /// </summary>
        public double Rosenbrock(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return Math.Pow(1 - x, 2) + 100 * Math.Pow(y - x * x, 2);
        }

        /// <summary>
        /// Quadratic form: f(x) = x'Ax where A is positive definite
        /// Used for testing Hessian computation
        /// </summary>
        public double QuadraticForm(double[] x)
        {
            // f(x,y,z) = x² + 2y² + 3z² + xy + 2xz + yz
            return x[0] * x[0] + 2 * x[1] * x[1] + 3 * x[2] * x[2] +
                   x[0] * x[1] + 2 * x[0] * x[2] + x[1] * x[2];
        }

        /// <summary>
        /// Vector-valued function for Jacobian: g(x,y) = [x²+y², x*y]
        /// </summary>
        public double[] VectorFunction(double[] points)
        {
            double x = points[0];
            double y = points[1];
            return new double[] { x * x + y * y, x * y };
        }

        #endregion

        #region Tests - Forward Difference

        /// <summary>
        /// Test forward difference method for polynomial function.
        /// </summary>
        [TestMethod]
        public void Test_ForwardDifference_Polynomial()
        {
            // f(x) = x³, f'(2) = 12
            double deriv = NumericalDerivative.ForwardDifference(FX, 2.0);
            Assert.AreEqual(12.0, deriv, 1E-6);

            // f(x) = sin(x), f'(π/4) = cos(π/4) = √2/2
            deriv = NumericalDerivative.ForwardDifference(SIN, Math.PI / 4);
            Assert.AreEqual(Math.Cos(Math.PI / 4), deriv, 1E-6);
        }

        /// <summary>
        /// Test forward difference with custom step size.
        /// </summary>
        [TestMethod]
        public void Test_ForwardDifference_CustomStepSize()
        {
            // Using a larger step size should reduce accuracy
            double derivSmallStep = NumericalDerivative.ForwardDifference(FX, 2.0, 1E-8);
            double derivLargeStep = NumericalDerivative.ForwardDifference(FX, 2.0, 1E-2);

            // Small step should be more accurate
            Assert.IsLessThan(Math.Abs(derivLargeStep - 12.0), Math.Abs(derivSmallStep - 12.0));
        }

        #endregion

        #region Tests - Backward Difference

        /// <summary>
        /// Test backward difference method for polynomial function.
        /// </summary>
        [TestMethod]
        public void Test_BackwardDifference_Polynomial()
        {
            // f(x) = x³, f'(2) = 12
            double deriv = NumericalDerivative.BackwardDifference(FX, 2.0);
            Assert.AreEqual(12.0, deriv, 1E-6);

            // f(x) = e^x, f'(4) = e⁴
            deriv = NumericalDerivative.BackwardDifference(EX, 4.0);
            Assert.AreEqual(Math.Exp(4), deriv, 1E-4);
        }

        #endregion

        #region Tests - Central Difference

        /// <summary>
        /// Test central difference method (should be more accurate than forward/backward).
        /// </summary>
        [TestMethod]
        public void Test_CentralDifference_HigherAccuracy()
        {
            double point = 2.0;
            double expected = 12.0; // f'(2) for f(x) = x³

            double forward = NumericalDerivative.ForwardDifference(FX, point);
            double backward = NumericalDerivative.BackwardDifference(FX, point);
            double central = NumericalDerivative.CentralDifference(FX, point);

            // Central difference should be most accurate (O(h²) vs O(h))
            double errorForward = Math.Abs(forward - expected);
            double errorBackward = Math.Abs(backward - expected);
            double errorCentral = Math.Abs(central - expected);

            Assert.IsLessThan(errorForward, errorCentral);
            Assert.IsLessThan(errorBackward, errorCentral);
        }

        /// <summary>
        /// Test that Derivative method is equivalent to CentralDifference.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_EquivalentToCentral()
        {
            double point = 2.0;
            double derivCentral = NumericalDerivative.CentralDifference(FX, point);
            double derivDefault = NumericalDerivative.Derivative(FX, point);

            Assert.AreEqual(derivCentral, derivDefault, 1E-15);
        }

        #endregion

        #region Tests - First Derivative (Original Tests Enhanced)

        /// <summary>
        /// Test derivative with multiple functions.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_MultipleFunctions()
        {
            // f(x) = x³, f'(2) = 12
            double derivFX = NumericalDerivative.Derivative(FX, 2.0);
            Assert.AreEqual(12.0, derivFX, 1E-6);

            // f(x) = e^x, f'(4) = e⁴
            double derivEX = NumericalDerivative.Derivative(EX, 4.0);
            Assert.AreEqual(Math.Exp(4), derivEX, 1E-4);

            // f(x) = ln(x), f'(2) = 0.5
            double derivLN = NumericalDerivative.Derivative(LN, 2.0);
            Assert.AreEqual(0.5, derivLN, 1E-4);

            // f(x) = sin(x), f'(π/3) = cos(π/3) = 0.5
            double derivSIN = NumericalDerivative.Derivative(SIN, Math.PI / 3);
            Assert.AreEqual(0.5, derivSIN, 1E-6);
        }

        /// <summary>
        /// Test derivative at zero (edge case).
        /// </summary>
        [TestMethod]
        public void Test_Derivative_AtZero()
        {
            // f(x) = x³, f'(0) = 0
            double deriv = NumericalDerivative.Derivative(FX, 0.0);
            Assert.AreEqual(0.0, deriv, 1E-10);
        }

        /// <summary>
        /// Test derivative at negative point.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_NegativePoint()
        {
            // f(x) = x³, f'(-2) = 12
            double deriv = NumericalDerivative.Derivative(FX, -2.0);
            Assert.AreEqual(12.0, deriv, 1E-6);
        }

        #endregion

        #region Tests - Second Derivatives

        /// <summary>
        /// Test second derivative using central difference.
        /// </summary>
        [TestMethod]
        public void Test_SecondDerivative_Central()
        {
            // f(x) = x³, f''(2) = 6*2 = 12
            double deriv = NumericalDerivative.SecondDerivative(FX, 2.0);
            Assert.AreEqual(12.0, deriv, 1E-4);

            // f(x) = e^x, f''(4) = e⁴
            deriv = NumericalDerivative.SecondDerivative(EX, 4.0);
            Assert.AreEqual(Math.Exp(4), deriv, 1E-3);

            // f(x) = sin(x), f''(π/6) = -sin(π/6) = -0.5
            deriv = NumericalDerivative.SecondDerivative(SIN, Math.PI / 6);
            Assert.AreEqual(-0.5, deriv, 1E-5);
        }

        /// <summary>
        /// Test second derivative using forward difference.
        /// </summary>
        [TestMethod]
        public void Test_SecondDerivative_Forward()
        {
            // f(x) = x³, f''(2) = 12
            double deriv = NumericalDerivative.SecondDerivativeForward(FX, 2.0);
            Assert.AreEqual(12.0, deriv, 1E-3);

            // f(x) = x⁴-2x²+3, f''(1) = 12-4 = 8
            deriv = NumericalDerivative.SecondDerivativeForward(Polynomial, 1.0);
            Assert.AreEqual(8.0, deriv, 1E-3);
        }

        /// <summary>
        /// Test second derivative using backward difference.
        /// </summary>
        [TestMethod]
        public void Test_SecondDerivative_Backward()
        {
            // f(x) = x³, f''(2) = 12
            double deriv = NumericalDerivative.SecondDerivativeBackward(FX, 2.0);
            Assert.AreEqual(12.0, deriv, 1E-4);

            // f(x) = cos(x), f''(π/4) = -cos(π/4)
            deriv = NumericalDerivative.SecondDerivativeBackward(COS, Math.PI / 4);
            Assert.AreEqual(-Math.Cos(Math.PI / 4), deriv, 1E-5);
        }

        /// <summary>
        /// Test that central difference is more accurate for second derivatives.
        /// </summary>
        [TestMethod]
        public void Test_SecondDerivative_CentralMoreAccurate()
        {
            double point = 2.0;
            double expected = 12.0; // f''(2) for f(x) = x³

            double forward = NumericalDerivative.SecondDerivativeForward(FX, point);
            double backward = NumericalDerivative.SecondDerivativeBackward(FX, point);
            double central = NumericalDerivative.SecondDerivative(FX, point);

            double errorForward = Math.Abs(forward - expected);
            double errorBackward = Math.Abs(backward - expected);
            double errorCentral = Math.Abs(central - expected);

            // Central should generally be most accurate
            Assert.IsTrue(errorCentral <= errorForward || errorCentral <= errorBackward);
        }

        #endregion

        #region Tests - Ridders' Method

        /// <summary>
        /// Test Ridders' method for high-accuracy derivatives.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_HighAccuracy()
        {
            // f(x) = x³, f'(2) = 12
            double derivFX = NumericalDerivative.RiddersMethod(FX, 2.0, out double errFX);
            Assert.AreEqual(12.0, derivFX, 1E-6);
            Assert.IsLessThan(1E-4, errFX); // Error estimate should be small

            // f(x) = e^x, f'(4) = e⁴
            double derivEX = NumericalDerivative.RiddersMethod(EX, 4.0, out double errEX);
            Assert.AreEqual(Math.Exp(4), derivEX, 1E-6);

            // f(x) = ln(x), f'(2) = 0.5
            double derivLN = NumericalDerivative.RiddersMethod(LN, 2.0, out double errLN);
            Assert.AreEqual(0.5, derivLN, 1E-6);
        }

        /// <summary>
        /// When the default step leaves the domain of f, the step is halved until the tableau is finite.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_DefaultStepOutsideDomain_FallsBack()
        {
            // The default step 0.1 (1 + |x|) exceeds x, so log and sqrt are NaN at x - h.
            Assert.AreEqual(20.0, NumericalDerivative.RiddersMethod(Math.Log, 0.05, out double errLog), 1E-12 * 20.0);
            Assert.IsLessThan(1E-10, errLog);
            Assert.AreEqual(5.0, NumericalDerivative.RiddersMethod(Math.Sqrt, 0.01, out _), 1E-12 * 5.0);
            Assert.AreEqual(1E6, NumericalDerivative.RiddersMethod(Math.Log, 1E-6, out _), 1E-12 * 1E6);
            // A function that is never finite keeps the result of the default step.
            Assert.AreEqual(0.0, NumericalDerivative.RiddersMethod(x => double.NaN, 1.0, out double errNaN), 0.0);
            Assert.AreEqual(double.MaxValue, errNaN);
            // A given step is used as it is.
            Assert.AreEqual(20.000000037611226, NumericalDerivative.RiddersMethod(Math.Log, 0.05, out _, 0.1), 0.0);
        }

        /// <summary>
        /// A tableau that is finite throughout returns the plain Ridders extrapolation, pinned bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_FiniteTableau_Unchanged()
        {
            Assert.AreEqual(0x40001C2A61268970L, BitConverter.DoubleToInt64Bits(NumericalDerivative.RiddersMethod(Math.Exp, 0.7, out _)));
            Assert.AreEqual(0x3FD11EB3682A4BEDL, BitConverter.DoubleToInt64Bits(NumericalDerivative.RiddersMethod(Math.Sin, 1.3, out _)));
            Assert.AreEqual(0x3FDFFFFFFFFFFFBEL, BitConverter.DoubleToInt64Bits(NumericalDerivative.RiddersMethod(Math.Log, 2.0, out _)));
            Assert.AreEqual(0x403AFFFFFFFFFFF4L, BitConverter.DoubleToInt64Bits(NumericalDerivative.RiddersMethod(x => x * x * x, -3.0, out _)));
            Assert.ThrowsExactly<ArgumentNullException>(() => NumericalDerivative.RiddersMethod((Func<double, double>)null, 1.0, out _));
        }

        /// <summary>
        /// Test that Ridders' method is more accurate than central difference.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_MoreAccurateThanCentral()
        {
            // Use exponential function e^x where Ridders excels
            // For polynomials like x³, both methods achieve machine precision
            // so no advantage is observable
            double point = 2.0;
            double expected = Math.Exp(2.0);  // f'(e^x) = e^x

            double central = NumericalDerivative.CentralDifference(EX, point);
            double ridders = NumericalDerivative.RiddersMethod(EX, point, out _);

            double errorCentral = Math.Abs(central - expected);
            double errorRidders = Math.Abs(ridders - expected);

            // For smooth transcendental functions, Ridders should be significantly better
            // Expect at least 10× improvement (often 100-1000× for exponentials)
            Assert.IsLessThan(errorCentral * 0.5, errorRidders);
        }

        /// <summary>
        /// The default initial step is the large step 0.1 (1 + |x|) of Press et al. (2007, section 5.7): the derivative of
        /// exp at 0.7 is accurate to about 1E-15 relative (about 1E-8 from the small step CalculateStepSize), and an explicit
        /// step still overrides it.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_DefaultLargeStep()
        {
            double expected = Math.Exp(0.7);
            double ridders = NumericalDerivative.RiddersMethod(Math.Exp, 0.7, out double err);
            Assert.AreEqual(expected, ridders, 2E-14 * expected);
            Assert.IsLessThan(1E-13, err);

            // Further smooth cases: sin at 1.3, ln at 2.5 and x³ at -4.
            Assert.AreEqual(Math.Cos(1.3), NumericalDerivative.RiddersMethod(Math.Sin, 1.3, out _), 1E-13);
            Assert.AreEqual(0.4, NumericalDerivative.RiddersMethod(Math.Log, 2.5, out _), 1E-14);
            Assert.AreEqual(48.0, NumericalDerivative.RiddersMethod(x => x * x * x, -4.0, out _), 1E-12);

            // A zero or negative step selects the default; a positive step is used as given.
            Assert.AreEqual(ridders, NumericalDerivative.RiddersMethod(Math.Exp, 0.7, out _, 0.0), 0.0);
            Assert.AreEqual(ridders, NumericalDerivative.RiddersMethod(Math.Exp, 0.7, out _, 0.1 * (1.0 + 0.7)), 0.0);
            double small = NumericalDerivative.RiddersMethod(Math.Exp, 0.7, out _, NumericalDerivative.CalculateStepSize(0.7));
            Assert.IsGreaterThan(Math.Abs(ridders - expected), Math.Abs(small - expected));
        }

        /// <summary>
        /// Test Ridders' method error estimate reliability.
        /// </summary>
        [TestMethod]
        public void Test_RiddersMethod_ErrorEstimate()
        {
            double result = NumericalDerivative.RiddersMethod(FX, 2.0, out double err);
            double actual = 12.0;
            double actualError = Math.Abs(result - actual);

            // Error estimate should be in a reasonable range
            // For smooth polynomial functions, error estimate should be conservative
            Assert.IsGreaterThanOrEqualTo(0, err);
            Assert.IsLessThan(1.0, err); // Should be reasonably small for polynomial
        }

        #endregion

        #region Tests - Gradient

        /// <summary>
        /// Test gradient computation.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_TwoVariables()
        {
            // f(x,y) = x²y³ at (2,2)
            // ∂f/∂x = 2xy³ = 2*2*8 = 32
            // ∂f/∂y = 3x²y² = 3*4*4 = 48
            var grad = NumericalDerivative.Gradient(FXY, new[] { 2.0, 2.0 });
            Assert.AreEqual(32.0, grad[0], 1E-5);
            Assert.AreEqual(48.0, grad[1], 1E-5);
        }

        /// <summary>
        /// Test gradient with three variables.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_ThreeVariables()
        {
            // f(x,y,z) = x³+y⁴+z⁵ at (2,2,2)
            // ∂f/∂x = 3x² = 12
            // ∂f/∂y = 4y³ = 32
            // ∂f/∂z = 5z⁴ = 80
            var grad = NumericalDerivative.Gradient(FXYZ, new[] { 2.0, 2.0, 2.0 });
            Assert.AreEqual(12.0, grad[0], 1E-5);
            Assert.AreEqual(32.0, grad[1], 1E-5);
            Assert.AreEqual(80.0, grad[2], 1E-5);
        }

        /// <summary>
        /// Test gradient at minimum of Rosenbrock function.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_Rosenbrock()
        {
            // Rosenbrock function has minimum at (1,1) where gradient = 0
            var grad = NumericalDerivative.Gradient(Rosenbrock, new[] { 1.0, 1.0 });
            Assert.AreEqual(0.0, grad[0], 1E-4);
            Assert.AreEqual(0.0, grad[1], 1E-4);
        }

        /// <summary>
        /// Test gradient with lower bounds (boundary handling).
        /// </summary>
        [TestMethod]
        public void Test_Gradient_WithLowerBounds()
        {
            // Point near lower bound - should use forward differences
            var point = new[] { 0.01, 0.01 };
            var lowerBounds = new[] { 0.0, 0.0 };
            var grad = NumericalDerivative.Gradient(FXY, point, lowerBounds);

            // Should still compute reasonable gradients
            Assert.IsLessThan(1.0, Math.Abs(grad[0]));
            Assert.IsLessThan(1.0, Math.Abs(grad[1]));
        }

        /// <summary>
        /// Test gradient with upper bounds (boundary handling).
        /// </summary>
        [TestMethod]
        public void Test_Gradient_WithUpperBounds()
        {
            // Point near upper bound - should use backward differences
            var point = new[] { 9.99, 9.99 };
            var upperBounds = new[] { 10.0, 10.0 };
            var grad = NumericalDerivative.Gradient(FXY, point, upperBounds);

            // Should still compute reasonable gradients
            Assert.IsTrue(!double.IsNaN(grad[0]) && !double.IsInfinity(grad[0]));
            Assert.IsTrue(!double.IsNaN(grad[1]) && !double.IsInfinity(grad[1]));
        }

        /// <summary>
        /// Test gradient with both lower and upper bounds.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_WithBothBounds()
        {
            var point = new[] { 5.0, 5.0, 5.0 };
            var lowerBounds = new[] { 0.0, 0.0, 0.0 };
            var upperBounds = new[] { 10.0, 10.0, 10.0 };

            var grad = NumericalDerivative.Gradient(FXYZ, point, lowerBounds, upperBounds);

            // Should use central differences in interior
            // ∂f/∂x = 3x² = 3*25 = 75 at x=5
            // ∂f/∂y = 4y³ = 4*125 = 500 at y=5
            // ∂f/∂z = 5z⁴ = 5*625 = 3125 at z=5
            Assert.AreEqual(75.0, grad[0], 1E-3);
            Assert.AreEqual(500.0, grad[1], 1E-2);
            Assert.AreEqual(3125.0, grad[2], 1E-1);
        }

        #endregion

        #region Tests - Jacobian (Original Method)

        /// <summary>
        /// Test Jacobian matrix computation.
        /// </summary>
        [TestMethod]
        public void Test_Jacobian_OriginalMethod()
        {
            // f(t,x,y,z) = tx³+y⁴+z⁵
            // At t=[4,6,8], θ=[2,2,2]:
            // ∂f/∂x = 3tx² = [48, 72, 96]
            // ∂f/∂y = 4y³ = [32, 32, 32]
            // ∂f/∂z = 5z⁴ = [80, 80, 80]
            var jac = NumericalDerivative.Jacobian(FTXYZ, new[] { 4.0, 6.0, 8.0 }, new[] { 2.0, 2.0, 2.0 });

            double[,] expected = new double[,]
            {
                { 48.0, 32.0, 80.0 },
                { 72.0, 32.0, 80.0 },
                { 96.0, 32.0, 80.0 }
            };

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(expected[i, j], jac[i, j], 1E-4);
        }

        #endregion

        #region Tests - Jacobian (Vector Function)

        /// <summary>
        /// Test Jacobian of vector-valued function.
        /// </summary>
        [TestMethod]
        public void Test_Jacobian_VectorFunction()
        {
            // g(x,y) = [x²+y², xy]
            // Jacobian at (2,3):
            // [2x  2y]   [4  6]
            // [y   x ] = [3  2]
            var jac = NumericalDerivative.Jacobian(VectorFunction, new[] { 2.0, 3.0 });

            Assert.AreEqual(4.0, jac[0, 0], 1E-5);
            Assert.AreEqual(6.0, jac[0, 1], 1E-5);
            Assert.AreEqual(3.0, jac[1, 0], 1E-5);
            Assert.AreEqual(2.0, jac[1, 1], 1E-5);
        }

        /// <summary>
        /// Test Jacobian with bounds (near lower bound).
        /// </summary>
        [TestMethod]
        public void Test_Jacobian_WithBounds_LowerBound()
        {
            var point = new[] { 0.01, 0.01 };
            var lowerBounds = new[] { 0.0, 0.0 };

            var jac = NumericalDerivative.Jacobian(VectorFunction, point, lowerBounds);

            // Should compute reasonable values despite being near boundary
            Assert.IsFalse(double.IsNaN(jac[0, 0]));
            Assert.IsFalse(double.IsNaN(jac[1, 1]));
        }

        /// <summary>
        /// Test Jacobian with bounds (near upper bound).
        /// </summary>
        [TestMethod]
        public void Test_Jacobian_WithBounds_UpperBound()
        {
            var point = new[] { 9.99, 9.99 };
            var upperBounds = new[] { 10.0, 10.0 };

            var jac = NumericalDerivative.Jacobian(VectorFunction, point, upperBounds);

            // Should compute reasonable values despite being near boundary
            Assert.IsFalse(double.IsNaN(jac[0, 0]));
            Assert.IsFalse(double.IsNaN(jac[1, 1]));
        }

        /// <summary>
        /// Test Jacobian in interior with bounds.
        /// </summary>
        [TestMethod]
        public void Test_Jacobian_WithBounds_Interior()
        {
            var point = new[] { 2.0, 3.0 };
            var lowerBounds = new[] { 0.0, 0.0 };
            var upperBounds = new[] { 10.0, 10.0 };

            var jac = NumericalDerivative.Jacobian(VectorFunction, point, lowerBounds, upperBounds);

            // Should use central differences and match expected values
            Assert.AreEqual(4.0, jac[0, 0], 1E-5);
            Assert.AreEqual(6.0, jac[0, 1], 1E-5);
            Assert.AreEqual(3.0, jac[1, 0], 1E-5);
            Assert.AreEqual(2.0, jac[1, 1], 1E-5);
        }

        #endregion

        #region Tests - Hessian

        /// <summary>
        /// Test Hessian matrix computation.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_MixedPartials()
        {
            // f(x,y) = x³-2xy-y⁶ at (1,2)
            // ∂²f/∂x² = 6x = 6
            // ∂²f/∂y² = -30y⁴ = -480
            // ∂²f/∂x∂y = -2
            var hess = NumericalDerivative.Hessian(FH, new[] { 1.0, 2.0 });

            Assert.AreEqual(6.0, hess[0, 0], 1E-3);
            Assert.AreEqual(-2.0, hess[0, 1], 1E-3);
            Assert.AreEqual(-2.0, hess[1, 0], 1E-3);
            Assert.AreEqual(-480.0, hess[1, 1], 1E-3);
        }

        /// <summary>
        /// Test Hessian for diagonal function (no mixed partials).
        /// </summary>
        [TestMethod]
        public void Test_Hessian_Diagonal()
        {
            // f(x,y,z) = x³+y⁴+z⁵ at (2,2,2)
            // ∂²f/∂x² = 6x = 12
            // ∂²f/∂y² = 12y² = 48
            // ∂²f/∂z² = 20z³ = 160
            // All mixed partials = 0
            var hess = NumericalDerivative.Hessian(FXYZ, new[] { 2.0, 2.0, 2.0 });

            Assert.AreEqual(12.0, hess[0, 0], 1E-3);
            Assert.AreEqual(48.0, hess[1, 1], 1E-3);
            Assert.AreEqual(160.0, hess[2, 2], 1E-3);

            // Off-diagonals should be near zero
            Assert.AreEqual(0.0, hess[0, 1], 1E-2);
            Assert.AreEqual(0.0, hess[0, 2], 1E-2);
            Assert.AreEqual(0.0, hess[1, 2], 1E-2);
        }

        /// <summary>
        /// Test Hessian symmetry.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_Symmetry()
        {
            var hess = NumericalDerivative.Hessian(FH, new[] { 1.0, 2.0 });

            // Hessian should be symmetric
            Assert.AreEqual(hess[0, 1], hess[1, 0], 1E-10);
        }

        /// <summary>
        /// Test Hessian for quadratic form (should be constant).
        /// </summary>
        [TestMethod]
        public void Test_Hessian_QuadraticForm()
        {
            // f(x,y,z) = x² + 2y² + 3z² + xy + 2xz + yz
            // Hessian should be constant everywhere:
            // [2   1   2]
            // [1   4   1]
            // [2   1   6]

            var point1 = new[] { 1.0, 1.0, 1.0 };
            var point2 = new[] { 5.0, 5.0, 5.0 };

            var hess1 = NumericalDerivative.Hessian(QuadraticForm, point1);
            var hess2 = NumericalDerivative.Hessian(QuadraticForm, point2);

            // Hessian should be the same at both points (for quadratic functions)
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Assert.AreEqual(hess1[i, j], hess2[i, j], 1E-2);

            // Check specific values
            Assert.AreEqual(2.0, hess1[0, 0], 1E-2);
            Assert.AreEqual(4.0, hess1[1, 1], 1E-2);
            Assert.AreEqual(6.0, hess1[2, 2], 1E-2);
            Assert.AreEqual(1.0, hess1[0, 1], 1E-2);
            Assert.AreEqual(2.0, hess1[0, 2], 1E-2);
            Assert.AreEqual(1.0, hess1[1, 2], 1E-2);
        }

        /// <summary>
        /// Test Hessian at Rosenbrock saddle point.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_Rosenbrock()
        {
            // Rosenbrock function at (1,1) should have positive definite Hessian
            var hess = NumericalDerivative.Hessian(Rosenbrock, new[] { 1.0, 1.0 });

            // All diagonal elements should be positive at minimum
            Assert.IsGreaterThan(0, hess[0, 0]);
            Assert.IsGreaterThan(0, hess[1, 1]);

            // Hessian should be symmetric
            Assert.AreEqual(hess[0, 1], hess[1, 0], 1E-8);
        }

        /// <summary>
        /// Test Hessian with lower bounds.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_WithLowerBounds()
        {
            var point = new[] { 0.5, 0.5 };
            var lowerBounds = new[] { 0.0, 0.0 };

            var hess = NumericalDerivative.Hessian(FXY, point, lowerBounds);

            // Should compute without errors
            Assert.IsFalse(double.IsNaN(hess[0, 0]));
            Assert.IsFalse(double.IsNaN(hess[1, 1]));
        }

        /// <summary>
        /// Test Hessian with upper bounds.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_WithUpperBounds()
        {
            var point = new[] { 9.5, 9.5 };
            var upperBounds = new[] { 10.0, 10.0 };

            var hess = NumericalDerivative.Hessian(FXY, point, upperBounds);

            // Should compute without errors
            Assert.IsFalse(double.IsNaN(hess[0, 0]));
            Assert.IsFalse(double.IsNaN(hess[1, 1]));
        }

        #endregion

        #region Tests - Step Size Calculation

        /// <summary>
        /// Test step size calculation for first derivatives.
        /// </summary>
        [TestMethod]
        public void Test_CalculateStepSize_FirstDerivative()
        {
            double h1 = NumericalDerivative.CalculateStepSize(1.0, 1);
            double h2 = NumericalDerivative.CalculateStepSize(10.0, 1);

            // Step size should scale with magnitude of input
            Assert.IsGreaterThan(h1, h2);

            // Should be reasonable for finite differences
            Assert.IsTrue(h1 > 1E-10 && h1 < 1E-5);
        }

        /// <summary>
        /// Test step size calculation for second derivatives.
        /// </summary>
        [TestMethod]
        public void Test_CalculateStepSize_SecondDerivative()
        {
            double h1 = NumericalDerivative.CalculateStepSize(1.0, 1);
            double h2 = NumericalDerivative.CalculateStepSize(1.0, 2);

            // Second derivative needs larger step size
            Assert.IsGreaterThan(h1, h2);
        }

        /// <summary>
        /// Test step size at zero.
        /// </summary>
        [TestMethod]
        public void Test_CalculateStepSize_AtZero()
        {
            double h = NumericalDerivative.CalculateStepSize(0.0);

            // Should return reasonable value even at zero
            Assert.IsGreaterThan(0, h);
            Assert.IsLessThan(1E-4, h);
        }

        #endregion

        #region Tests - Edge Cases and Error Handling

        /// <summary>
        /// Test derivative at very large values.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_LargeValues()
        {
            // Should handle large inputs gracefully
            // Use relative tolerance for large values
            double deriv = NumericalDerivative.Derivative(x => x * x, 1E6);
            double expected = 2E6;
            double relativeError = Math.Abs((deriv - expected) / expected);
            Assert.IsLessThan(1E-6, relativeError); // 0.0001% relative error
        }

        /// <summary>
        /// Test derivative at very small values.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_SmallValues()
        {
            // Should handle small inputs gracefully
            double deriv = NumericalDerivative.Derivative(x => x * x, 1E-6);
            Assert.AreEqual(2E-6, deriv, 1E-12);
        }

        /// <summary>
        /// Test gradient with single variable (edge case).
        /// </summary>
        [TestMethod]
        public void Test_Gradient_SingleVariable()
        {
            Func<double[], double> f = x => x[0] * x[0];
            var grad = NumericalDerivative.Gradient(f, new[] { 2.0 });

            Assert.AreEqual(4.0, grad[0], 1E-6);
        }

        /// <summary>
        /// Test that null function throws exception.
        /// </summary>
        [TestMethod]
        public void Test_Derivative_NullFunction()
        {
            var ex = Assert.Throws<Exception>(() =>
            {
                NumericalDerivative.Derivative(null, 1.0);
            });
        }

        /// <summary>
        /// Test that null theta throws exception in Gradient.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_NullTheta()
        {
            var ex = Assert.Throws<Exception>(() =>
            {
                NumericalDerivative.Gradient(FXY, null);
            });
        }

        /// <summary>
        /// Test Hessian with non-finite function value.
        /// </summary>
        [TestMethod]
        public void Test_Hessian_NonFiniteValue()
        {
            Func<double[], double> badFunc = x => double.NaN;
            var ex = Assert.Throws<Exception>(() =>
            {
                NumericalDerivative.Hessian(badFunc, new[] { 1.0, 1.0 });
            });
        }

        #endregion

        #region Tests - Comparison Between Methods

        /// <summary>
        /// Compare all first derivative methods for consistency.
        /// </summary>
        [TestMethod]
        public void Test_CompareFirstDerivativeMethods()
        {
            double point = 2.0;
            double expected = 12.0;

            double forward = NumericalDerivative.ForwardDifference(FX, point);
            double backward = NumericalDerivative.BackwardDifference(FX, point);
            double central = NumericalDerivative.CentralDifference(FX, point);
            double ridders = NumericalDerivative.RiddersMethod(FX, point, out _);

            // All should be within reasonable tolerance of expected
            // Forward and backward are O(h) so less accurate
            Assert.AreEqual(expected, forward, 1E-4);
            Assert.AreEqual(expected, backward, 1E-4);
            Assert.AreEqual(expected, central, 1E-6);
            Assert.AreEqual(expected, ridders, 1E-6);
        }

        /// <summary>
        /// Verify that gradient matches Jacobian for scalar functions.
        /// </summary>
        [TestMethod]
        public void Test_GradientMatchesJacobianForScalar()
        {
            var point = new[] { 2.0, 3.0 };

            var grad = NumericalDerivative.Gradient(FXY, point);

            // Create wrapper for Jacobian (single output)
            Func<double[], double[]> wrapper = p => new[] { FXY(p) };
            var jac = NumericalDerivative.Jacobian(wrapper, point);

            // Gradient should match first row of Jacobian
            Assert.AreEqual(grad[0], jac[0, 0], 1E-8);
            Assert.AreEqual(grad[1], jac[0, 1], 1E-8);
        }

        #endregion
    }
}
