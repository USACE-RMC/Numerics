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

using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using Numerics.Mathematics.Optimization;
using Numerics;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Augmented Lagrange optimization algorithm
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet"> 
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_AugmentedLagrange
    {
        // direct reference for these 2?
        /// <summary>
        /// Solution is from Water Economics graduate course
        /// </summary>
        [TestMethod]
        public void Test_1()
        {
            var constraint = new Constraint((x) => { return Tools.Sum(x); }, 3, 22, ConstraintType.EqualTo);
            Func<double[], double> func = (double[] x) =>
            {
                var NB = new double[3];
                for (int i = 0; i < 3; i++)
                {
                    NB[i] = (20 * x[i] - x[i] * x[i] - 24) / Math.Pow(1.10, i);
                }
                return -Tools.Sum(NB);
            };
            var innerSolver = new BFGS(func, 3, new double[] { 7, 7, 8 }, new double[] { double.MinValue, double.MinValue, double.MinValue }, new double[] { double.MaxValue, double.MaxValue, double.MaxValue });
            var solver = new AugmentedLagrange(func, innerSolver, new IConstraint[] { constraint });
            solver.Minimize();

            // Point
            Assert.AreEqual(7.583082, solver.BestParameterSet.Values[0], 1E-3);
            Assert.AreEqual(7.341390, solver.BestParameterSet.Values[1], 1E-3);
            Assert.AreEqual(7.075529, solver.BestParameterSet.Values[2], 1E-3);
            // Function
            Assert.AreEqual(188.5655, -solver.BestParameterSet.Fitness, 1E-3);
            // Multiplier
            Assert.AreEqual(4.833835, solver.Lambda[0], 1E-3);
        }

        /// <summary>
        /// Solution is from Water Economics graduate course
        /// </summary>
        [TestMethod]
        public void Test_2()
        {
            var constraint = new Constraint((x) => { return Tools.Sum(x); }, 2, 100, ConstraintType.EqualTo);
            Func<double[], double> func = (double[] x) =>
            {
                var NB = new double[2];
                NB[0] = 60 * x[0] - 0.5 * x[0] * x[0];
                NB[1] = (64 * x[1] - 0.5 * x[1] * x[1]) / 1.5;
                return -Tools.Sum(NB);
            };
            var innerSolver = new BFGS(func, 2, new double[] { 50, 50 }, new double[] { double.MinValue, double.MinValue }, new double[] { double.MaxValue, double.MaxValue });
            var solver = new AugmentedLagrange(func, innerSolver, new IConstraint[] { constraint });
            solver.Minimize();

            // Point
            Assert.AreEqual(50.4, solver.BestParameterSet.Values[0], 1E-3);
            Assert.AreEqual(49.6, solver.BestParameterSet.Values[1], 1E-3);
            // Function
            Assert.AreEqual(3050.133, -solver.BestParameterSet.Fitness, 1E-3);
            // Multiplier
            Assert.AreEqual(9.6, solver.Lambda[0], 1E-3);
        }

        /// <summary>
        /// Test the Lagrange Augmented Algorithm with example problem 5.2 from "Risk Modeling, Assessment, and Management"
        /// </summary>
        [TestMethod]
        public void Test_Haimes_5_2()
        {
            Func<double[], double> primaryFunc = (double[] x) =>
            {
                return Math.Pow(x[0] - 2, 2) + Math.Pow(x[1] - 4, 2) + 5;
            };
            Func<double[], double> secondaryFunc = (double[] x) =>
            {
                return Math.Pow(x[0] - 6, 2) + Math.Pow(x[1] - 10, 2) + 6;
            };

            // Set up inner solver
            var initial = new double[] { 5, 5 };
            var lower = new double[] { 0, 0 };
            var upper = new double[] { 10, 10 };
            var innerSolver = new BFGS(primaryFunc, 2, initial, lower, upper);
            // Set up constraint
            var constraint = new Constraint(secondaryFunc, 2, 13.31, ConstraintType.LesserThanOrEqualTo);
            // Solve
            var solver = new AugmentedLagrange(primaryFunc, innerSolver, new IConstraint[] { constraint });
            solver.Minimize();

            // Point
            Assert.AreEqual(4.5, solver.BestParameterSet.Values[0], 1E-2);
            Assert.AreEqual(7.75, solver.BestParameterSet.Values[1], 1E-2);
            // Function
            Assert.AreEqual(25.31, solver.BestParameterSet.Fitness, 1E-2);
            // Multiplier
            Assert.AreEqual(1.67, solver.Mu[0], 1E-2);

        }

        /// <summary>
        /// Test the Augmented Lagrange algorithm on the Rosenbrock Function constrained to a disk
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <see href="https://en.wikipedia.org/wiki/Test_functions_for_optimization"/>
        /// </remarks>
        [TestMethod]
        public void Test_RosenbrockDisk()
        {
            Func<double[], double> cfunc = (double[] x) =>
            {
                return (x[0] * x[0]) + (x[1] * x[1]);
            };


            Func<double[], double> func = (double[] x) =>
            {
                return Math.Pow(1 - x[0], 2) + 100 * Math.Pow(x[1] - x[0] * x[0], 2);
            };

            // Set up inner solver
            var initial = new double[] { 0, 0 };
            var lower = new double[] { -1.5, -1.5 };
            var upper = new double[] { 1.5, 1.5};
            var innerSolver = new BFGS(func, 2, initial, lower, upper);
            // Set up constraint
            var constraint = new Constraint(cfunc, 2, 2, ConstraintType.LesserThanOrEqualTo);
            // Solve
            var solver = new AugmentedLagrange(func, innerSolver, new IConstraint[] { constraint });
            solver.Minimize();

            // Point
            Assert.AreEqual(1d, solver.BestParameterSet.Values[0], 1E-4);
            Assert.AreEqual(1d, solver.BestParameterSet.Values[1], 1E-4);
            // Function
            Assert.AreEqual(0d, solver.BestParameterSet.Fitness, 1E-4);
            // Multiplier
            Assert.AreEqual(0d, solver.Mu[0]);
        }
    }
}
