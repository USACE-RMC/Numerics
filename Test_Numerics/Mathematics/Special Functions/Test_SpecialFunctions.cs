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
    /// Unit tests for Special Functions, outside of the Beta and Gamma Class.
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
    public class Test_SpecialFunctions
    {
        /// <summary>
        /// Test the error function
        /// </summary>
        [TestMethod]
        public void Test_Erf()
        {
            var testX = new double[] { 0d, 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d, 1d, 1.1d, 1.2d, 1.3d, 1.4d, 1.5d, 1.6d, 1.7d, 1.8d, 1.9d };
            var testValid = new double[] { 0d, 0.112462916018285d, 0.222702589210478d, 0.328626759459127d, 0.428392355046668d, 0.520499877813047d, 0.603856090847926d, 0.677801193837418d, 0.742100964707661d, 0.796908212422832d, 0.842700792949715d, 0.880205069574082d, 0.910313978229635d, 0.934007944940652d, 0.952285119762649d, 0.966105146475311d, 0.976348383344644d, 0.983790458590775d, 0.989090501635731d, 0.992790429235258d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Erf.Function(testX[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the error complement function
        /// </summary>
        [TestMethod]
        public void Test_Erfc()
        {
            var testX = new double[] { 0d, 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d, 1d, 1.1d, 1.2d, 1.3d, 1.4d, 1.5d, 1.6d, 1.7d, 1.8d, 1.9d };
            var testValid = new double[] { 1d, 0.887537083981715d, 0.777297410789522d, 0.671373240540873d, 0.571607644953332d, 0.479500122186953d, 0.396143909152074d, 0.322198806162582d, 0.257899035292339d, 0.203091787577168d, 0.157299207050285d, 0.119794930425918d, 0.0896860217703646d, 0.0659920550593475d, 0.0477148802373512d, 0.0338948535246893d, 0.023651616655356d, 0.0162095414092254d, 0.0109094983642693d, 0.00720957076474253d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Erf.Erfc(testX[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the inverse of the error function
        /// </summary>
        [TestMethod]
        public void Test_InverseErf()
        {
            var testY = new double[] { 0d, 0.112462916018285d, 0.222702589210478d, 0.328626759459127d, 0.428392355046668d, 0.520499877813047d, 0.603856090847926d, 0.677801193837418d, 0.742100964707661d, 0.796908212422832d, 0.842700792949715d, 0.880205069574082d, 0.910313978229635d, 0.934007944940652d, 0.952285119762649d, 0.966105146475311d, 0.976348383344644d, 0.983790458590775d, 0.989090501635731d, 0.992790429235258d };
            var testValid = new double[] { 0d, 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d, 1d, 1.1d, 1.2d, 1.3d, 1.4d, 1.5d, 1.6d, 1.7d, 1.8d, 1.9d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Erf.InverseErf(testY[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the inverse of the complemented error function
        /// </summary>
        [TestMethod]
        public void Test_InverseErfc()
        {
            var testY = new double[] { 1d, 0.887537083981715d, 0.777297410789522d, 0.671373240540873d, 0.571607644953332d, 0.479500122186953d, 0.396143909152074d, 0.322198806162582d, 0.257899035292339d, 0.203091787577168d, 0.157299207050285d, 0.119794930425918d, 0.0896860217703646d, 0.0659920550593475d, 0.0477148802373512d, 0.0338948535246893d, 0.023651616655356d, 0.0162095414092254d, 0.0109094983642693d, 0.00720957076474253d };
            var testValid = new double[] { 0d, 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d, 1d, 1.1d, 1.2d, 1.3d, 1.4d, 1.5d, 1.6d, 1.7d, 1.8d, 1.9d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Erf.InverseErfc(testY[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the Debye function
        /// </summary>
        [TestMethod]
        public void Test_Debye()
        {
            var testX = new double[] { 0.1, 1.0, 2.8, 9.5, 10, 15, 25, 100 };
            var testValid = new double[] { 0.9629999, 0.6744156, 0.3099952, 0.02241066, 0.01929577, 0.005771263, 0.001246836, 1.948182e-05 };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Debye.Function(testX[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// Test the factorial function
        /// </summary>
        [TestMethod]
        public void Test_Factorial()
        {
            var testX = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            var testValid = new double[] { 1d, 2d, 6d, 24d, 120d, 720d, 5040d, 40320d, 362880d, 3628800d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testResults.Length; i++)
            {
                testResults[i] = Factorial.Function(testX[i]);
                Assert.AreEqual(testValid[i], testResults[i]);
            }
        }

        /// <summary>
        /// Test the log of the factorial functions
        /// </summary>
        [TestMethod]
        public void Test_LogFactorial()
        {
            var testX = new int[] { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
            var testFact = new double[] { 1d, 2d, 6d, 24d, 120d, 720d, 5040d, 40320d, 362880d, 3628800d };
            var testValid = new double[testX.Length];
            for (int i = 0; i < testValid.Length; i++)
            {
                testValid[i] = Math.Log(testFact[i]);
            }
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testResults.Length; i++)
            {
                testResults[i] = Factorial.LogFactorial(testX[i]);
                Assert.AreEqual(testValid[i], testResults[i]);
            }
        }

        /// <summary>
        /// Test the binomial coefficient calculation
        /// </summary>
        [TestMethod]
        public void Test_BinomialCoefficient()
        {
            var testN = new int[] { 4, 4, 4, 6, 6, 6, 8, 8, 8, 10, 10, 10 };
            var testK = new int[] { 1, 2, 3, 2, 3, 4, 2, 3, 4, 2, 3, 4 };
            var testValid = new double[] { 4d, 6d, 4d, 15d, 20d, 15d, 28d, 56d, 70d, 45d, 120d, 210d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testResults.Length; i++)
            {
                testResults[i] = Factorial.BinomialCoefficient(testN[i], testK[i]);
                Assert.AreEqual(testValid[i], testResults[i]);
            }
        }

        /// <summary>
        /// Test the function that finds all combinations of m within n without replacement by validating the size of the outputs
        /// </summary>
        [TestMethod]
        public void Test_CombinationsNum()
        {
            var cc = Factorial.AllCombinations(5);
            // Total number of possible combinations
            double possible = Math.Pow(2, 5) - 1;

            // Length of cc should be the possible number of combinations * the number of elements in each combination (5)
            Assert.AreEqual(possible * 5, cc.Length);

            // How many of subsets of combinations there should be
            // For example, there are 5 ways to have only one #1 in the array, with the other 4 elements being #0
            // This can be calculated from the binomial coefficient (i.e 5 choose 1)
            var valid = new double[5];
            for (int i = 0; i < 5; i++)
            {
                valid[i] = Factorial.BinomialCoefficient(5, i + 1);
            }
            // All of the subsets should add up to the total number of possible combinations
            Assert.AreEqual(valid[0] + valid[1] + valid[2] + valid[3] + valid[4], possible);

            // Check cc for the subsets
            var counts = new double[5];
            for (int i = 0; i < possible; i++)
            {
                int counter = 0;
                for (int j = 0; j < 5; j++)
                {
                    if (cc[i, j] == 1)
                    {
                        counter++;
                    }
                }
                counts[counter - 1]++;
            }

            // Validate that the number subsets found in cc match the calculated number
            for (int i = 0; i < 5; i++)
            {
                Assert.AreEqual(valid[i], counts[i]);
            }

        }

        /// <summary>
        /// Test the function that finds all combinations of m within n without replacement by finding matches between the true int[,] and the output int[,]
        /// </summary>
        [TestMethod]
        public void Test_Combinations()
        {
            int[,] valid = new int[7,3]
                { 
                    { 1, 0 ,0 },
                    { 0, 1, 0 },
                    { 0, 0, 1 },
                    { 1, 1, 0 },
                    { 1, 0, 1 },
                    { 0, 1, 1 },
                    { 1, 1, 1 }
                };

            var actual = Factorial.AllCombinations(3);

            int found = 0;
            // Checking every row in valid against rows in actual to find a match
            for(int i = 0; i < 7; i++)
            {
                for(int j = 0; j < 7; j++)
                {
                    if (valid[i,0] == actual[j, 0])
                    {
                        if (valid[i,1] == actual[j, 1])
                        {
                            if (valid[i,2] == actual[j, 2])
                            {
                                // Match found!
                                found++;
                            }
                        }
                    }
                }
            }

            // There are 7 possible outcomes (2^n - 1)
            Assert.AreEqual(7, found);
        }

        /// <summary>
        /// Test the polynomial function
        /// </summary>
        [TestMethod]
        public void Test_Polynomial()
        {
            var coeffs = new double[] { 3, 5, 7 };
            double x = 4;

            double valid = 3 + (5 * x) + (7 * x * x);
            double actual = Evaluate.Polynomial(coeffs, x);
            Assert.AreEqual(valid, actual);
        }

        /// <summary>
        /// Test the polynomial function that takes the coefficients in the reverse order
        /// </summary>
        [TestMethod]
        public void Test_PolynomialRev()
        {
            var coeffs = new double[] { 3, 5, 7 };
            double x = 4;

            // test without a given n
            double valid = (3 * x * x) + (5 * x) + 7;
            double actual = Evaluate.PolynomialRev(coeffs, x);
            Assert.AreEqual(valid, actual);

            // test with a given n
            int n = 1;
            double validN = (3 * x) + 5;
            double actualN = Evaluate.PolynomialRev(coeffs, x, n);
            Assert.AreEqual(validN, actualN);
        }

        /// <summary>
        /// Test the polynomial function that takes the coefficients in reverse order, with the first coefficient automatically always 1
        /// </summary>
        [TestMethod]
        public void Test_PolynomialRev_1()
        {
            var coeffs = new double[] { 3, 5, 7 };
            double x = 4;

            double valid = (x * x * x) + (3 * x * x) + (5 * x) + 7;
            double actual = Evaluate.PolynomialRev_1(coeffs, x);
            Assert.AreEqual(valid, actual);
        }
    }
}