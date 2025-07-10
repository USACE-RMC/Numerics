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
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the MLSL optimization algorithm
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
    public class Test_MLSL
    {
        /// <summary>
        /// Test the MLSL algorithm with a simple 3-dimensional test function.
        /// </summary>
        [TestMethod]
        public void Test_FXYZ()
        {
            var initial = new double[] { 0.2d, 0.5d, 0.5d };
            var lower = new double[] { 0d, 0d, 0d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MLSL(TestFunctions.FXYZ, 3, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            double x = solution[0];
            double y = solution[1];
            double z = solution[2];
            double validX = 0.125d;
            double validY = 0.2d;
            double validZ = 0.35d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
            Assert.AreEqual(z, validZ, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the De Jong Function in 5-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            var initial = new double[] { 1.0d, -1.0d, 2.0d, -2.0d, 1.0d };
            var lower = new double[] { -5.12d, -5.12d, -5.12d, -5.12d, -5.12d };
            var upper = new double[] { 5.12d, 5.12d, 5.12d, 5.12d, 5.12d };
            var solver = new MLSL(TestFunctions.DeJong, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Sum of Power functions in 3-D.
        /// </summary>
        [TestMethod]
        public void Test_SumOfPowerFunctions()
        {
            var initial = new double[] { 0.5d, -0.5d, 0.5d };
            var lower = new double[] { -1d, -1d, -1d };
            var upper = new double[] { 1d, 1d, 1d };
            var solver = new MLSL(TestFunctions.SumOfPowerFunctions, 3, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0d, 0.0d, 0.0d };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Rosenbrock Function
        /// </summary>
        [TestMethod]
        public void Test_Rosenbrock()
        {
            var initial = new double[] { 0, 0, 0, 0, 0 };
            var lower = new double[] { -2.048, -2.048, -2.048, -2.048, -2.048 };
            var upper = new double[] { 2.048, 2.048, 2.048, 2.048, 2.048 };
            var solver = new MLSL(TestFunctions.Rosenbrock, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 1.0, 1.0, 1.0, 1.0, 1.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Booth Function
        /// </summary>
        [TestMethod]
        public void Test_Booth()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MLSL(TestFunctions.Booth, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 1.0d;
            var validY = 3.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Matyas Function
        /// </summary>
        [TestMethod]
        public void Test_Matyas()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -10d, -10d };
            var upper = new double[] { 10d, 10d };
            var solver = new MLSL(TestFunctions.Matyas, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the McCormick Function
        /// </summary>
        [TestMethod]
        public void Test_McCormick()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -1.5d, -3d };
            var upper = new double[] { 4d, 4d };
            var solver = new MLSL(TestFunctions.McCormick, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -1.9133;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = -0.54719d;
            var validY = -1.54719d;
            Assert.AreEqual(x, validX, 1E-3);
            Assert.AreEqual(y, validY, 1E-3);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Beale Function
        /// </summary
        [TestMethod]
        public void Test_Beale()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -4.5d, -4.5d };
            var upper = new double[] { 4.5d, 4.5d };
            var solver = new MLSL(TestFunctions.Beale, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 3.0d;
            var validY = 0.5d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Goldstein-Price Function
        /// </summary>
        [TestMethod]
        public void Test_GoldsteinPrice()
        {
            var initial = new double[] { -1d, 1d };
            var lower = new double[] { -2d, -2d };
            var upper = new double[] { 2d, 2d };
            var solver = new MLSL(TestFunctions.GoldsteinPrice, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 3.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = -1.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Rastrigin Function
        /// </summary>
        [TestMethod]
        public void Test_Rastrigin()
        {
            var initial = new double[] { 1, 1, 1, 1, 1 };
            var lower = new double[] { -5.12, -5.12, -5.12, -5.12, -5.12 };
            var upper = new double[] { 5.12, 5.12, 5.12, 5.12, 5.12 };
            // Need to run a lot of starts
            var solver = new MLSL(TestFunctions.Rastrigin, 5, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var valid = new double[] { 0.0, 0.0, 0.0, 0.0, 0.0 };
            for (int i = 0; i < valid.Length; i++)
                Assert.AreEqual(solution[i], valid[i], 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Ackley Function
        /// </summary>
        [TestMethod]
        public void Test_Ackley()
        {
            var initial = new double[] { 1d, 1d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5d, 5d };
            var solver = new MLSL(TestFunctions.Ackley, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the three hump camel Function
        /// </summary>
        [TestMethod]
        public void Test_ThreeHumpCamel()
        {
            var initial = new double[] { 2d, -2d };
            var lower = new double[] { -5d, -5d };
            var upper = new double[] { 5, 5d };
            var solver = new MLSL(TestFunctions.ThreeHumpCamel, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 0.0d;
            var validY = 0.0d;
            Assert.AreEqual(x, validX, 1E-4);
            Assert.AreEqual(y, validY, 1E-4);
        }

        /// <summary>
        /// Test the MLSL algorithm with the Eggholder Function
        /// </summary>
        [TestMethod]
        public void Test_Eggholder()
        {
            var initial = new double[] { 0d, 0d };
            var lower = new double[] { -512d, -512d };
            var upper = new double[] { 512d, 512d };
            // Need to run a lot of starts
            var solver = new MLSL(TestFunctions.Eggholder, 2, initial, lower, upper);
            solver.SampleSize = 200;
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = -959.6407;
            Assert.AreEqual(F, trueF, 0.5);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 512d;
            var validY = 404.2319d;
            Assert.AreEqual(x, validX, 1E-1);
            Assert.AreEqual(y, validY, 1E-1);
        }

        /// <summary>
        /// Test the MLSL algorithm with the tp2 Function
        /// </summary>
        [TestMethod]
        public void Test_TP2()
        {
            var initial = new double[] { 2d, 2d };
            var lower = new double[] { 0d, 0d };
            var upper = new double[] { 2d, 2d };
            var solver = new MLSL(TestFunctions.tp2, 2, initial, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            var solution = solver.BestParameterSet.Values;
            var x = solution[0];
            var y = solution[1];
            var validX = 1d;
            var validY = 0.666667d;

            bool match1 = Math.Abs(x - validX) < 1E-4 && Math.Abs(y - validY) < 1E-4;
            bool match2 = Math.Abs(x - validY) < 1E-4 && Math.Abs(y - validX) < 1E-4;
            Assert.IsTrue(match1 || match2);
        }
    }
}
