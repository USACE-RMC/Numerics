using Numerics.Distributions;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Test functions used to test the integration algorithms. There are both multidimensional and single dimensional functions
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    public class Integrands
    {
        /// <summary>
        /// Test function. The integral of x^3, should equal 0.25
        /// </summary>
        public static double FX3(double x)
        {
            return Math.Pow(x, 3d);
        }

        /// <summary>
        /// Test function. The integral of Cos(x), should equal ~ 1.6829419
        /// </summary>
        public static double Cosine(double x)
        {
            return Math.Cos(x);
        }

        /// <summary>
        /// Test function. The integral of Sine(x), should equal ~ 0.459697694131
        /// </summary>
        public static double Sine(double x)
        {
            return Math.Sin(x);
        }

        /// <summary>
        /// Test function. The integral 2rd order polynomial, should equal 57.
        /// </summary>
        public static double FXX(double x)
        {
            return 0.5 + 24 * x + 3 * x * x;
        }

        /// <summary>
        /// Test function. The integral 3rd order polynomial, should equal 89.
        /// </summary>
        public static double FXXX(double x)
        {
            return 0.5 + 24 * x + 3 * x * x + 8 * x * x * x;
        }

        /// <summary>
        /// Test function. The integral of Pi. Should equal ~3.14
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        public static double PI2D(double x, double y)
        {
            return (x * x + y * y < 1) ? 1 : 0;
        }

        /// <summary>
        /// Test function. The integral of Pi. Should equal ~3.14
        /// </summary>
        /// <param name="vals">Array of values.</param>
        public static double PI(double[] vals)
        {
            var x = vals[0];
            var y = vals[1];
            return (x * x + y * y < 1) ? 1 : 0;
        }

        /// <summary>
        /// Test function from GNU Scientific Library, GSL. Result should equal 1.3932039296856768591842462603255
        /// </summary>
        public static double GSL(double[] x)
        {
            double A = 1.0 / (Math.PI * Math.PI * Math.PI);
            return A / (1.0 - Math.Cos(x[0]) * Math.Cos(x[1]) * Math.Cos(x[2]));
        }

        /// <summary>
        /// Mean values for the 20-dimensional sum-of-normals integration fixture.
        /// </summary>
        public static double[] mu20 = new double[] { 10, 30, 17, 99, 68, 26, 35, 55, 13, 59, 12, 28, 49, 54, 20, 47, 12, 76, 70, 57 };

        /// <summary>
        /// Standard deviation values for the 20-dimensional sum-of-normals integration fixture.
        /// </summary>
        public static double[] sigma20 = new double[] { 2, 15, 5, 14, 7, 24, 29, 22, 22, 1, 3, 28, 19, 18, 4, 24, 23, 26, 26, 19 };

        /// <summary>
        /// Test function for the sum of normal distributions. Max D=20.
        /// </summary>
        /// <param name="p">Array of probability values.</param>
        public static double SumOfNormals(double[] p)
        {
            double result = 0;
            for (int i = 0; i < p.Length; i++)
            {
                result += mu20[i] + sigma20[i] * Normal.StandardZ(p[i]);
            }
            return  result;
        }

        /// <summary>
        /// Test function for the sum of 2 normal distributions.
        /// </summary>
        /// <param name="p1">The probability of distribution 1.</param>
        /// <param name="p2">The probability of distribution 2.</param>
        public static double SumOfNormals2D(double p1, double p2)
        {
            double result = 0;
            result += mu20[0] + sigma20[0] * Normal.StandardZ(p1);
            result += mu20[1] + sigma20[1] * Normal.StandardZ(p2);
            return result;
        }

        // The six Genz test families in two dimensions over the unit square [0,1]^2, each with its
        // exact integral in closed form. Reference: Genz, A. (1984). "Testing multidimensional
        // integration routines." Tools, Methods and Languages for Scientific and Engineering
        // Computation, 81-94.

        /// <summary>
        /// Genz oscillatory family: cos(2*pi*w1 + c1*x + c2*y).
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The phase shift parameter.</param>
        public static double GenzOscillatory(double x, double y, double c1, double c2, double w1)
        {
            return Math.Cos(2d * Math.PI * w1 + c1 * x + c2 * y);
        }

        /// <summary>
        /// Exact integral of the Genz oscillatory family over the unit square:
        /// (4/(c1*c2)) * sin(c1/2) * sin(c2/2) * cos(2*pi*w1 + c1/2 + c2/2).
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The phase shift parameter.</param>
        public static double GenzOscillatoryExact(double c1, double c2, double w1)
        {
            return 4d / (c1 * c2) * Math.Sin(c1 / 2d) * Math.Sin(c2 / 2d) * Math.Cos(2d * Math.PI * w1 + c1 / 2d + c2 / 2d);
        }

        /// <summary>
        /// Genz product-peak family: 1 / ((c1^-2 + (x-w1)^2) * (c2^-2 + (y-w2)^2)).
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the peak.</param>
        /// <param name="w2">The y-location of the peak.</param>
        public static double GenzProductPeak(double x, double y, double c1, double c2, double w1, double w2)
        {
            return 1d / ((1d / (c1 * c1) + (x - w1) * (x - w1)) * (1d / (c2 * c2) + (y - w2) * (y - w2)));
        }

        /// <summary>
        /// Exact integral of the Genz product-peak family over the unit square: the product of
        /// c * (atan(c*(1-w)) + atan(c*w)) over both axes.
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the peak.</param>
        /// <param name="w2">The y-location of the peak.</param>
        public static double GenzProductPeakExact(double c1, double c2, double w1, double w2)
        {
            double ix = c1 * (Math.Atan(c1 * (1d - w1)) + Math.Atan(c1 * w1));
            double iy = c2 * (Math.Atan(c2 * (1d - w2)) + Math.Atan(c2 * w2));
            return ix * iy;
        }

        /// <summary>
        /// Genz corner-peak family: (1 + c1*x + c2*y)^-3.
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        public static double GenzCornerPeak(double x, double y, double c1, double c2)
        {
            double b = 1d + c1 * x + c2 * y;
            return 1d / (b * b * b);
        }

        /// <summary>
        /// Exact integral of the Genz corner-peak family over the unit square:
        /// (1/(2*c1*c2)) * (1 - 1/(1+c1) - 1/(1+c2) + 1/(1+c1+c2)).
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        public static double GenzCornerPeakExact(double c1, double c2)
        {
            return 1d / (2d * c1 * c2) * (1d - 1d / (1d + c1) - 1d / (1d + c2) + 1d / (1d + c1 + c2));
        }

        /// <summary>
        /// Genz Gaussian family: exp(-c1^2*(x-w1)^2 - c2^2*(y-w2)^2).
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the peak.</param>
        /// <param name="w2">The y-location of the peak.</param>
        public static double GenzGaussian(double x, double y, double c1, double c2, double w1, double w2)
        {
            return Math.Exp(-c1 * c1 * (x - w1) * (x - w1) - c2 * c2 * (y - w2) * (y - w2));
        }

        /// <summary>
        /// Exact integral of the Genz Gaussian family over the unit square: the product of
        /// sqrt(pi)/(2c) * (erf(c*(1-w)) + erf(c*w)) over both axes.
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the peak.</param>
        /// <param name="w2">The y-location of the peak.</param>
        public static double GenzGaussianExact(double c1, double c2, double w1, double w2)
        {
            double ix = Math.Sqrt(Math.PI) / (2d * c1) * (Numerics.Mathematics.SpecialFunctions.Erf.Function(c1 * (1d - w1)) + Numerics.Mathematics.SpecialFunctions.Erf.Function(c1 * w1));
            double iy = Math.Sqrt(Math.PI) / (2d * c2) * (Numerics.Mathematics.SpecialFunctions.Erf.Function(c2 * (1d - w2)) + Numerics.Mathematics.SpecialFunctions.Erf.Function(c2 * w2));
            return ix * iy;
        }

        /// <summary>
        /// Genz continuous (C0) family: exp(-c1*|x-w1| - c2*|y-w2|). Continuous with a gradient
        /// discontinuity along both peak lines.
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the ridge.</param>
        /// <param name="w2">The y-location of the ridge.</param>
        public static double GenzContinuous(double x, double y, double c1, double c2, double w1, double w2)
        {
            return Math.Exp(-c1 * Math.Abs(x - w1) - c2 * Math.Abs(y - w2));
        }

        /// <summary>
        /// Exact integral of the Genz continuous family over the unit square: the product of
        /// (2 - exp(-c*w) - exp(-c*(1-w)))/c over both axes.
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the ridge.</param>
        /// <param name="w2">The y-location of the ridge.</param>
        public static double GenzContinuousExact(double c1, double c2, double w1, double w2)
        {
            double ix = (2d - Math.Exp(-c1 * w1) - Math.Exp(-c1 * (1d - w1))) / c1;
            double iy = (2d - Math.Exp(-c2 * w2) - Math.Exp(-c2 * (1d - w2))) / c2;
            return ix * iy;
        }

        /// <summary>
        /// Genz discontinuous family: exp(c1*x + c2*y) inside the corner box [0,w1]x[0,w2] and zero
        /// outside it.
        /// </summary>
        /// <param name="x">The x-coordinate.</param>
        /// <param name="y">The y-coordinate.</param>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the discontinuity.</param>
        /// <param name="w2">The y-location of the discontinuity.</param>
        public static double GenzDiscontinuous(double x, double y, double c1, double c2, double w1, double w2)
        {
            return x > w1 || y > w2 ? 0d : Math.Exp(c1 * x + c2 * y);
        }

        /// <summary>
        /// Exact integral of the Genz discontinuous family over the unit square:
        /// ((exp(c1*w1)-1)/c1) * ((exp(c2*w2)-1)/c2).
        /// </summary>
        /// <param name="c1">The x-direction difficulty parameter.</param>
        /// <param name="c2">The y-direction difficulty parameter.</param>
        /// <param name="w1">The x-location of the discontinuity.</param>
        /// <param name="w2">The y-location of the discontinuity.</param>
        public static double GenzDiscontinuousExact(double c1, double c2, double w1, double w2)
        {
            return (Math.Exp(c1 * w1) - 1d) / c1 * ((Math.Exp(c2 * w2) - 1d) / c2);
        }

    }
}
