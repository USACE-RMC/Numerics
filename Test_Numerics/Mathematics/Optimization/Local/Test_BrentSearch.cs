using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.Optimization;

namespace Mathematics.Optimization
{
    /// <summary>
    /// Unit tests for the Brent optimization algorithm
    /// </summary>
    [TestClass]
    public class Test_BrentSearch
    {
        /// <summary>
        /// Test to find the minimum of a one dimensional function using Brent's method
        /// </summary>
        [TestMethod]
        public void Test_Minimize()
        {
            double lower = -3d;
            double upper = 3d;
            var solver = new BrentSearch(TestFunctions.FX, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = 1.0d;
            Assert.AreEqual(X, trueX, 1E-4);
        }

        /// <summary>
        /// Test to find the maximum of a one dimensional function using Brent's method
        /// </summary>
        [TestMethod]
        public void Test_Maximize()
        {
            double lower = -3;
            double upper = 3d;
            var solver = new BrentSearch(TestFunctions.FX, lower, upper);
            solver.Maximize();
            double F = -1*solver.BestParameterSet.Fitness;
            double trueF = 9.4815;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = -1.6667d;
            Assert.AreEqual(X, trueX, 1E-4);
        }

        /// <summary>
        /// Test the Brent algorithm with De Jong's function in 1-D.
        /// </summary>
        [TestMethod]
        public void Test_DeJong()
        {
            double lower = -5.12d;
            double upper = 5.12d;
            var solver = new BrentSearch((x) => { return TestFunctions.DeJong(new double[] { x }); }, lower, upper);
            solver.Minimize();
            double F = solver.BestParameterSet.Fitness;
            double trueF = 0.0;
            Assert.AreEqual(F, trueF, 1E-4);
            double X = solver.BestParameterSet.Values[0];
            double trueX = 0.0;
            Assert.AreEqual(X, trueX, 1E-4);
        }

    }
}
