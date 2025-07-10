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
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.SpecialFunctions;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Unit tests for the Beta class
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> Description: </b>
    /// Some of these special functions are tested directly with MS Excel and with values provided by:
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// Others must be tested indirectly through the continuous distribution classes.
    /// </remarks>
    [TestClass]
    public class Test_Beta
    {
        /// <summary>
        /// Test the Beta Function
        /// </summary>
        [TestMethod]
        public void Test_Function()
        {
            var testX = new double[] { 1d, 0.2d, 1d, 0.4d, 1d, 0.6d, 0.8d, 6d, 6d, 6d, 6d, 6d, 7d, 5d, 4d, 3d, 2d };
            var testY = new double[] { 1d, 1d, 0.2d, 1d, 0.4d, 1d, 1d, 6d, 5d, 4d, 3d, 2d, 7d, 5d, 4d, 3d, 2d };
            var testValid = new double[] { 1d, 5d, 5d, 2.5d, 2.5d, 1.666667d, 1.25d, 0.0003607504d, 0.0007936508d, 0.001984127d, 0.005952381d, 0.0238095d, 0.00008325008d, 0.001587302d, 0.007142857d, 0.03333333d, 0.1666667d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Beta.Function(testX[i], testY[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the Incomplete Beta function Ix(a,b)
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete()
        {
            double a = 5;
            double b = 4;
            double x = 0.5;

            double actual = Beta.Incomplete(a, b, x);
            double expected = 0.36328125;

            Assert.AreEqual(expected, actual, 1e-6);
        }

        /// <summary>
        /// Test the continued fraction expansion #1 for the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_Incbf()
        {
            double actual = Beta.Incbcf(4, 2, 4.2);
            Assert.AreEqual(-0.23046874999999992, actual);
        }

        /// <summary>
        /// Test the continue fraction expansion #2 for incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_Incbd()
        {
            double actual = Beta.Incbd(4, 2, 4.2);
            Assert.AreEqual(0.7375, actual);
        }

        /// <summary>
        /// Test the inverse of the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_IncompleteInverse()
        {
            double actual = Beta.IncompleteInverse(0.5, 0.6, 0.1);
            Assert.AreEqual(0.019145979066925722, actual);

        }

        /// <summary>
        /// Test the power series for the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_PowerSeries()
        {
            double actual = Beta.PowerSeries(4, 2, 4.2);
            Assert.AreEqual(-3671.801280000001, actual);
        }

        /// <summary>
        /// Test the incomplete beta function ratio
        /// </summary>
        /// <remarks>
        /// This function was test against values calculated using the "Ibeta( )" and "Cbeta( )" functions the R "zipfR" package
        /// <para></para>
        /// References:
        /// <para></para>
        /// Evert, Stefan and Baroni, Marco (2007). "zipfR: Word frequency distributions in
        /// R." In Proceedings of the 45th Annual Meeting of the Association for Computational
        /// Linguistics, Posters and Demonstrations Sessions, pages 29-32, Prague, Czech
        /// Republic. (R package version 0.6-70 of 2020-10-10)
        /// </remarks>
        [TestMethod()]
        public void Test_IncompleteRatio()
        {
            var testX = new double[] { 0.2d, 0.2d, 0.2d, 0.5d, 0.5d, 0.5d, 0.8d, 0.8d, 0.8d };
            var testP = new double[] { 2d, 3d, 4d, 2d, 3d, 4d, 2d, 3d, 4d};
            var testQ = new double[] { 4d, 6d, 9d, 4d, 6d, 9d, 4d, 6d, 9d};
            var testBeta = new double[testX.Length];
            for (int i = 0; i < testBeta.Length; i++)
            {
                testBeta[i] = Math.Log(Gamma.Function(testP[i]) * Gamma.Function(testQ[i]) / Gamma.Function(testP[i] + testQ[i]) ); 
            }
            var testValid = new double[] { 0.26272d, 0.2030822d, 0.2054311d, 0.8125d, 0.8554688d, 0.927002d, 0.99328d, 0.9987686d, 0.9999378d};
            var testResults = new double[testValid.Length];

            for(int i = 0; i < testResults.Length; i++)
            {
                testResults[i] = Beta.IncompleteRatio(testX[i], testP[i], testQ[i], testBeta[i]);
                Assert.AreEqual(testValid[i], testResults[i], 1e-4);
            }
        }
    }
}
