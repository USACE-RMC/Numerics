using System;

namespace Mathematics.RootFinding
{
    /// <summary>
    /// Functions designed to test root finding algorithms.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    public class TestFunctions
    {
        /// <summary>
        /// A quadratic test function. [0, 4] x = sqrt(2)
        /// </summary>
        public static double Quadratic(double x)
        {
            return Math.Pow(x, 2) - 2;
        }

        /// <summary>
        /// First derivative of quadratic function. 
        /// </summary>
        public static double Quadratic_Deriv(double x)
        {
            return 2 * x;
        }

        /// <summary>
        /// A cubic test function. [-1, 5] x = 1.32472
        /// </summary>
        public static double Cubic(double x)
        {
            return x * x * x - x - 1d;
        }

        /// <summary>
        /// First derivative of cubic function.
        /// </summary>
        public static double Cubic_Deriv(double x)
        {
            return 3d * (x * x) - 1d;
        }

        /// <summary>
        /// A trigonometric test function. [0, 3.14] x = 1.12191713 
        /// </summary>
        public static double Trigonometric(double x)
        {
            return 2 * Math.Sin(x) - 3 * Math.Cos(x) - 0.5;
        }

        /// <summary>
        /// First derivative of the trigonometric function.
        /// </summary>
        public static double Trigonometric_Deriv(double x)
        {
            return 2 * Math.Cos(x) + 3 * Math.Sin(x);
        }

        /// <summary>
        /// An exponential test function. [-2, 2] x = 0.567143290
        /// </summary>
        public static double Exponential(double x)
        {
            return Math.Exp(-x) - x;
        }

        /// <summary>
        /// First derivative of the exponential function.
        /// </summary>
        public static double Exponential_Deriv(double x)
        {
            return -Math.Exp(-x) - 1;
        }

        /// <summary>
        /// A power test function. [0, 2] x = 1.0
        /// </summary>
        public static double Power(double x)
        {
            return Math.Pow(x, 10) - 1;
        }

        /// <summary>
        /// First derivative of the power function.
        /// </summary>
        public static double Power_Deriv(double x)
        {
            return 10 * Math.Pow(x, 9);
        }
    }
}
