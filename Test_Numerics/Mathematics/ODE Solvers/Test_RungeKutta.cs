using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.ODESolvers;
using System;

namespace Mathematics.ODESolvers
{
    /// <summary>
    /// Unit tests for the RungeKutta method
    /// </summary>
    /// <remarks>
    /// <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    /// Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> References: </b>
    /// <para>
    /// The example values were taken from OK State's MATH4513 course
    /// </para>
    /// <see href="https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf"/>
    /// </remarks>
    [TestClass]
    public class Test_RungeKutta
    {
        /// <summary>
        /// Test the second order Runge-Kutta method
        /// </summary>
        [TestMethod]
        public void Test_2ndRK()
        {
            var testValid = new double[] { 0.5d, 1.425639364649936d, 2.640859085770477d, 4.009155464830968d, 5.305471950534675d };

            Func<double, double, double> ode = (double t, double y) =>
            {
                return y - (t * t) + 1;
            };

            double initial = 0.5d;
            double start = 0d;
            double end = 2d;
            int steps = 5;
            var testResults = RungeKutta.SecondOrder(ode, initial, start, end, steps);

            for (int i = 0; i < testValid.Length; i++)
            {
                Assert.AreEqual(testValid[i], testResults[i], 1);
            }
        }

        /// <summary>
        /// Test the fourth order Runge-Kutta method
        /// </summary>
        [TestMethod]
        public void Test_4thRK()
        {
            var testValid = new double[] { 0.5d, 1.425639364649936d, 2.640859085770477d, 4.009155464830968d, 5.305471950534675d };

            Func<double, double, double> ode = (double t, double y) =>
            {
                return y - (t * t) + 1;
            };

            double initial = 0.5d;
            double start = 0d;
            double end = 2d;
            int steps = 5;
            var testResults = RungeKutta.FourthOrder(ode, initial, start, end, steps);

            for (int i = 0; i < testValid.Length; i++)
            {
                Assert.AreEqual(testValid[i], testResults[i], 1E-2);
            }
        }


        /// <summary>
        /// Test the second fourth order Runge-Kutta method that allows you to specify the time step size (dt) instead of the number of time steps
        [TestMethod]
        public void Test_4thRK_2()
        {
            var testValid = new double[] { 0.5d, 1.425639364649936d, 2.640859085770477d, 4.009155464830968d, 5.305471950534675d };

            Func<double, double, double> ode = (double t, double y) =>
            {
                return y - (t * t) + 1;
            };

            double initial = 0.5d;
            double start = 0d;
            double dt = 0.5d;
            var testResults = new double[testValid.Length];
            testResults[0] = initial;

            // loop to define the results, as this method only returns one value at a time
            for (int i = 1; i < testResults.Length; i++)
            {
                testResults[i] = RungeKutta.FourthOrder(ode, initial, start, dt);
                initial = testResults[i];
                start += dt;
            }

            for (int i = 0; i < testValid.Length; i++)
            {
                Assert.AreEqual(testValid[i], testResults[i], 1E-2);
            }
        }

        /// <summary>
        /// Test the Runge-Kutta-Fehlberg method
        /// </summary>
        [TestMethod]
        public void Test_RKF()
        {
            var testValid = new double[] { 0.5d, 1.425639364649936d, 2.640859085770477d, 4.009155464830968d, 5.305471950534675d };

            Func<double, double, double> ode = (double t, double y) =>
            {
                return y - (t * t) + 1;
            };

            double initial = 0.5d;
            double start = 0d;
            double dt = 0.5d;
            double dtMin = 0.001;
            var testResults = new double[testValid.Length];
            testResults[0] = initial;

            // loop to define the results, as this method only returns one value at a time
            for (int i = 1; i < testResults.Length; i++)
            {
                testResults[i] = RungeKutta.Fehlberg(ode, initial, start, dt, dtMin);
                initial = testResults[i];
                start += dt;
            }

            for (int i = 0; i < testValid.Length; i++)
            {
                Assert.AreEqual(testValid[i], testResults[i], 1E-4);
            }
        }


        /// <summary>
        /// Test the adaptive Runge-Kutta-Cash-Karp method
        /// </summary>
        [TestMethod]
        public void Test_RKCK()
        {
            var testValid = new double[] { 0.5d, 1.425639364649936d, 2.640859085770477d, 4.009155464830968d, 5.305471950534675d };

            Func<double, double, double> ode = (double t, double y) =>
            {
                return y - (t * t) + 1;
            };

            double initial = 0.5d;
            double start = 0d;
            double dt = 0.5d;
            double dtMin = 0.001;
            var testResults = new double[testValid.Length];
            testResults[0] = initial;

            // loop to define the results, as this method only returns one value at a time
            for (int i = 1; i < testResults.Length; i++)
            {
                testResults[i] = RungeKutta.CashKarp(ode, initial, start, dt, dtMin);
                initial = testResults[i];
                start += dt;
            }

            for (int i = 0; i < testValid.Length; i++)
            {
                Assert.AreEqual(testValid[i], testResults[i], 1E-3);
            }
        }
    }
}
