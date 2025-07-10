/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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

        public static double[] mu20 = new double[] { 10, 30, 17, 99, 68, 26, 35, 55, 13, 59, 12, 28, 49, 54, 20, 47, 12, 76, 70, 57 };
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

    }
}
