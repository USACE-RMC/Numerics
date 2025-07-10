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

namespace Mathematics.Optimization
{
    /// <summary>
    /// Functions designed to test optimization algorithms.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    public class TestFunctions
    {

        /// <summary>
        /// A simple one-dimensional test function
        /// </summary>
        /// <remarks>
        /// <para>
        /// Test area is restricted to [-3, 3]
        /// </para>
        /// <para>
        /// <b> Local minimum:</b>
        /// f(1) = 0
        /// </para>
        /// <para>
        /// <b> Local maximum:</b>
        /// f(-1.6667) = 9.4815
        /// </para>
        /// </remarks>
        public static double FX(double x)
        {
            double F = (x + 3d) * Math.Pow(x - 1d, 2d);
            return F;
        }

        /// <summary>
        /// A simple 3-dimensional test function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0.125, 0.2, 0.35) = 0
        /// </para>
        /// </remarks>
        public static double FXYZ(double[] parms)
        {
            double x = parms[0];
            double y = parms[1];
            double z = parms[2];
            double F = Math.Pow(4d * x - 0.5d, 2d) + Math.Pow(3d * y - 0.6d, 2d) + Math.Pow(2d * z - 0.7d, 2d);
            return F;
        }

        /// <summary>
        /// The De Jong (or sphere) function. 
        /// </summary>
        /// <remarks>
        /// <para>
        /// Test area is usually restricted to hypercube [-5.12, 5.12]
        /// </para>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// </remarks>
        public static double DeJong(double[] x)
        {
            int n = x.Length;
            double F = 0;
            for (int i = 0; i < n; i++)
                F += Math.Pow(x[i], 2);
            return F;
        }

        /// <summary>
        /// The Sum of Different Power Functions. 
        /// </summary>
        /// <remarks>
        /// <para>
        /// Test area is usually restricted to hypercube [-1, 1]
        /// </para>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// </remarks>
        public static double SumOfPowerFunctions(double[] x)
        {
            int n = x.Length;
            double F = 0;
            for (int i = 0; i < n; i++)
                F += Math.Pow(Math.Abs(x[i]), i + 2);
            return F;
        }

        /// <summary>
        /// The Rosenberg Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Test area is usually restricted to hypercube [-2.048, 2.048]
        /// </para>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(1, 1) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Rosenbrock(double[] parms)
        {
            int n = parms.Length;
            double F = 0;
            for (int i = 0; i < n - 1; i++)
                F += 100 * Math.Pow(parms[i + 1] - parms[i] * parms[i], 2) + Math.Pow(1 - parms[i], 2);
            return F;
        }

        /// <summary>
        /// The Booth Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(1, 3) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Booth(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = Math.Pow(x + 2 * y - 7, 2) + Math.Pow(2 * x + y - 5, 2);
            return F;
        }

        /// <summary>
        /// The Matyas Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Matyas(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = 0.26 * (x * x + y * y) - 0.48 * x * y;
            return F;
        }

        /// <summary>
        /// The McCormick Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(-0.54719, -1.54719) = -1.9133
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double McCormick(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = Math.Sin(x + y) + Math.Pow(x - y, 2) - 1.5 * x + 2.5 * y + 1.0;
            return F;
        }

        /// <summary>
        /// The Rastrigin Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Rastrigin(double[] parms)
        {
            double A = 10;
            int n = parms.Length;
            double F = A * n;
            for (int i = 0; i < n; i++)
                F += parms[i] * parms[i] - A * Math.Cos(2 * Math.PI * parms[i]);
            return F;
        }

        /// <summary>
        /// The Ackley Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Ackley(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = -20 * Math.Exp(-0.2 * Math.Sqrt(0.5 * (x * x + y * y))) - Math.Exp(0.5 * (Math.Cos(2 * Math.PI * x) + Math.Cos(2 * Math.PI * y))) + Math.E + 20;
            return F;
        }

        /// <summary>
        /// The Beale Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(3, 0.5) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Beale(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = Math.Pow(1.5 - x + x * y, 2) + Math.Pow(2.25 - x + x * y * y, 2) + Math.Pow(2.625 - x + x * y * y * y, 2);
            return F;
        }

        /// <summary>
        /// The Goldstein-Price Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, -1) = 3
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double GoldsteinPrice(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = (1 + Math.Pow(x + y + 1, 2) * (19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 3 * y * y)) *
                       (30 + Math.Pow(2 * x - 3 * y, 2) * (18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y));
            return F;
        }

        /// <summary>
        /// The Bukin Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(-10, 1) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Bukin(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = 100 * Math.Sqrt(Math.Abs(y - 0.01 * x * x)) + 0.01 * Math.Abs(x + 10);
            return F;
        }

        /// <summary>
        /// The three hump camel Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(0, 0) = 0
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double ThreeHumpCamel(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = 2 * Math.Pow(x, 2) - 1.05 * Math.Pow(x, 4) + Math.Pow(x, 6) / 6 + x * y + Math.Pow(y, 2);
            return F;
        }

        /// <summary>
        /// The Eggholder Function.
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(512, 404.2319) = -959.6407
        /// </para>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization" />
        /// </remarks>
        public static double Eggholder(double[] parms)
        {
            var x = parms[0];
            var y = parms[1];
            double F = -(y + 47) * Math.Sin(Math.Sqrt(Math.Abs((x / 2) + (y + 47)))) - x * Math.Sin(Math.Sqrt(Math.Abs(x - (y + 47))));
            return F;
        }

        /// <summary>
        /// tp2 function - Multiple local optima and 2 global minimizers, by virtue of symmetry
        /// </summary>
        /// <remarks>
        /// <para>
        /// <b> Global minimum:</b>
        /// f(2/3, 1) = f(1, 2/3) = 0
        /// </para>
        /// </remarks>
        public static double tp2(double[] parms)
        {
            var p1 = parms[0];
            var p2 = parms[1];
            int i;
            double x, F = 0;

            for (i = 0; i <= 86; i++)
            {
                x = 3.1 + 0.15 * i;
                F = F + Math.Pow(((Math.Sin(p1 * x) + Math.Sin(p2 * x)) - (Math.Sin((2.0 / 3.0) * x) + Math.Sin(1 * x))), 2);
            }
            return F;
        }


    }
}
