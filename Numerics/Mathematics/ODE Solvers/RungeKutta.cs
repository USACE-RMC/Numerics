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

namespace Numerics.Mathematics.ODESolvers
{

    /// <summary>
    /// Runge-Kutta method for solving Ordinary Differential Equations (ODE).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// In numerical analysis, the Runge–Kutta methods are a family of implicit and explicit iterative methods,
    /// which include the well-known routine called the Euler Method,
    /// used in temporal discretization for the approximate solutions of ordinary differential equations.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods"/>
    /// </para>
    /// </remarks>
    public class RungeKutta
    {

        /// <summary>
        /// Second Order Runge-Kutta method.
        /// </summary>
        /// <param name="f">ODE function F(t,Y).</param>
        /// <param name="initialValue">Initial value of Y.</param>
        /// <param name="startTime">Start time.</param>
        /// <param name="endTime">End time.</param>
        /// <param name="timeSteps">The number of time steps between the start and end time.</param>
        /// <returns>
        /// An array of the ODE function solved at the number of time steps specified, equally spaced between the start and end time.
        /// </returns>
        public static double[] SecondOrder(Func<double, double, double> f, double initialValue, double startTime, double endTime, int timeSteps)
        {
            double dt = (endTime - startTime) / (timeSteps - 1);
            double k1 = 0d;
            double k2 = 0d;
            double t = startTime;
            var y = new double[timeSteps];
            double y0 = initialValue;
            y[0] = y0;
            for (int i = 1; i < timeSteps; i++)
            {
                k1 = f(t, y0);
                k2 = f(t + dt, y0 + k1 * dt);
                y[i] = y0 + dt * 0.5d * (k1 + k2);
                t += dt;
                y0 = y[i];
            }

            return y;
        }

        /// <summary>
        /// The Fourth Order Runge-Kutta method.
        /// </summary>
        /// <param name="f">ODE function F(t,Y).</param>
        /// <param name="initialValue">Initial value of Y.</param>
        /// <param name="startTime">Start time.</param>
        /// <param name="endTime">End time.</param>
        /// <param name="timeSteps">The number of time steps between the start and end time.</param>
        /// <returns>
        /// An array of the ODE function solved at the number of time steps specified, equally space between the start and end time.
        /// </returns>
        public static double[] FourthOrder(Func<double, double, double> f, double initialValue, double startTime, double endTime, int timeSteps)
        {
            double dt = (endTime - startTime) / (timeSteps - 1);
            double k1 = 0d;
            double k2 = 0d;
            double k3 = 0d;
            double k4 = 0d;
            double t = startTime;
            var y = new double[timeSteps];
            double y0 = initialValue;
            y[0] = y0;
            for (int i = 1; i < timeSteps; i++)
            {
                k1 = f(t, y0);
                k2 = f(t + dt / 2d, y0 + k1 * dt / 2d);
                k3 = f(t + dt / 2d, y0 + k2 * dt / 2d);
                k4 = f(t + dt, y0 + k3 * dt);
                y[i] = y0 + dt / 6d * (k1 + 2d * k2 + 2d * k3 + k4);
                t += dt;
                y0 = y[i];
            }
            return y;
        }

        /// <summary>
        /// The Fourth Order Runge-Kutta method.
        /// </summary>
        /// <param name="f">The ODE function f(t,Y).</param>
        /// <param name="initialValue">Initial value of Y.</param>
        /// <param name="startTime">Start time.</param>
        /// <param name="dt">The time step size.</param>
        /// <returns>
        /// The ODE solved at the given time step (i.e. the ODE solved at startTime + dt)
        /// </returns>
        public static double FourthOrder(Func<double, double, double> f, double initialValue, double startTime, double dt)
        {
            double t = startTime;
            double y = initialValue;
            double k1 = f(t, y);
            double k2 = f(t + dt / 2d, y + k1 * dt / 2d);
            double k3 = f(t + dt / 2d, y + k2 * dt / 2d);
            double k4 = f(t + dt, y + k3 * dt);          
            return y + dt / 6d * (k1 + 2d * k2 + 2d * k3 + k4);
        }

        /// <summary>
        /// The Runge-Kutta-Fehlberg method.
        /// </summary>
        /// <param name="f">The ODE function dY/dt = f(t,Y).</param>
        /// <param name="initialValue">Initial value of Y.</param>
        /// <param name="startTime">Start time.</param>
        /// <param name="dt">The user defined step size, which is also the maximum step size.</param>
        /// <param name="dtMin">The minimum step size.</param>
        /// <param name="tolerance">The absolute tolerance. Default = 1E-3.</param>
        /// <remarks>
        /// References:
        /// <list type="bullet">
        /// <item>Runge–Kutta–Fehlberg: <see href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method"/></item>
        /// <item>Adaptive step size: <see href="https://en.wikipedia.org/wiki/Adaptive_step_size"/></item>
        /// </list>
        /// </remarks>
        /// <returns>
        /// The ODE solved at the given time step (i.e. the ODE solved at startTime + dt)
        /// </returns>
        public static double Fehlberg(Func<double, double, double> f, double initialValue, double startTime, double dt, double dtMin, double tolerance = 1E-3)
        {
            double t = startTime;
            double tf = startTime + dt;
            double y = initialValue;
            double k1, k2, k3, k4, k5, k6, rk4, rk5, error;
            double h = dt;
            double minTolerance = 1E-2 * tolerance;

            while (t < tf)
            {
                // Make sure the time step is constrained
                if (h < dtMin) h = dtMin;
                if (h > dt) h = dt;
                if (t + h > tf) h = tf - t;

                // Compute RK4, RK5, and error = |RK4 - RK5|
                k1 = h * f(t, y);
                k2 = h * f(t + 1d / 4d * h, y + 1d / 4d * k1);
                k3 = h * f(t + 3d / 8d * h, y + 3d / 32d * k1 + 9d / 32d * k2);
                k4 = h * f(t + 12d / 13d * h, y + 1932d / 2197d * k1 - 7200d / 2197d * k2 + 7296d / 2197d * k3);
                k5 = h * f(t + 1d / 1d * h, y + 439d / 216d * k1 - 8d * k2 + 3680d / 513d * k3 - 845d / 4104d * k4);
                k6 = h * f(t + 1d / 2d * h, y - 8d / 27d * k1 + 2d * k2 - 3544d / 2565d * k3 + 1859d / 4104d * k4 - 11d / 40d * k5);
                rk4 = y + 25d / 216d * k1 + 1408d / 2565d * k3 + 2197d / 4104d * k4 - 1d / 5d * k6;
                rk5 = y + 16d / 135d * k1 + 6656d / 12825d * k3 + 28561d / 56430d * k4 - 9d / 50d * k5 + 2d / 55d * k6;
                error = Math.Max(minTolerance, Math.Abs(rk4 - rk5));

                // Check convergence
                if (error <= tolerance || h <= dtMin)
                {
                    // accept the step
                    t += h;
                    y = rk5;

                    // Determine if the step size should be increased
                    if (error < minTolerance)
                        h *= 2d;
                }
                else
                {
                    // Reject the step and decrease the step size
                    h /= 2d;
                }
            }

            return y;
        }

        /// <summary>
        /// The adaptive Runge-Kutta-Cash-Karp method.
        /// </summary>
        /// <param name="f">The ODE function dY/dt = f(t,Y).</param>
        /// <param name="initialValue">Initial value of Y.</param>
        /// <param name="startTime">Start time.</param>
        /// <param name="dt">The user defined step size, which is also the maximum step size.</param>
        /// <param name="dtMin">The minimum step size.</param>
        /// <param name="tolerance">The absolute tolerance. Default = 1E-3.</param>
        /// <remarks>
        /// References:
        /// <list type="bullet">
        /// <item>Runge–Kutta–Cash-Karp: <see href="https://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method"/></item>
        /// <item>Adaptive step size: <see href="https://en.wikipedia.org/wiki/Adaptive_step_size"/></item>
        /// </list>
        /// </remarks>
        /// <returns>
        /// The ODE solved at the given time step (i.e. the ODE solved at startTime + dt)
        /// </returns>
        public static double CashKarp(Func<double, double, double> f, double initialValue, double startTime, double dt, double dtMin, double tolerance = 1E-3)
        {
            double t = startTime;
            double tf = startTime + dt;
            double y = initialValue;
            double k1, k2, k3, k4, k5, k6, rk4, rk5, error;
            double h = dt;
            double minTolerance = 1E-2 * tolerance;

            while (t < tf)
            {
                // Make sure the time step is constrained
                if (h < dtMin) h = dtMin;
                if (h > dt) h = dt;
                if (t + h > tf) h = tf - t;

                // Compute RK4, RK5, and error = |RK4 - RK5|
                k1 = h * f(t, y);
                k2 = h * f(t + 1d / 5d * h, y + 1d / 5d * k1);
                k3 = h * f(t + 3d / 10d * h, y + 3d / 40d * k1 + 9d / 40d * k2);
                k4 = h * f(t + 3d / 5d * h, y + 3d / 10d * k1 - 9d / 10d * k2 + 6d / 5d * k3);
                k5 = h * f(t + 1d / 1d * h, y - 11d / 54d * k1 + 5d / 2d * k2 - 70d / 27d * k3 + 35d / 27d * k4);
                k6 = h * f(t + 7d / 8d * h, y + 1631d / 55296d * k1 + 175d / 512d * k2 + 575d / 13824d * k3 + 44275d / 110592d * k4 + 253d / 4096d * k5);
                rk4 = y + 37d / 378d * k1 + 250d / 621d * k3 + 125d / 594d * k4 + 512d / 1771d * k6;
                rk5 = y + 2825d / 27648d * k1 + 18575d / 48384d * k3 + 13525d / 55296d * k4 + 277d / 14336d * k5 + 1d / 4d * k6;
                error = Math.Max(minTolerance, Math.Abs(rk4 - rk5));

                // Check convergence
                if (error <= tolerance || h <= dtMin)
                {
                    // accept the step
                    t += h;
                    y = rk5;

                    // Determine if the step size should be increased
                    if (error < minTolerance)
                        h *= 2d;
                }
                else
                {
                    // Reject the step and decrease the step size
                    h /= 2d;
                }

            }

            return y;
        }

    }
}