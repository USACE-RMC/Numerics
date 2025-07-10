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
