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
