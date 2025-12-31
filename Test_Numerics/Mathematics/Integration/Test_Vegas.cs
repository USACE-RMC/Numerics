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
using Numerics.Mathematics.Integration;
using System;

namespace Mathematics.Integration
{
    /// <summary>
    /// Unit tests for the Vegas method of adaptive multidimensional integration
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </remarks>
    [TestClass]
    public class Test_Vegas
    {
        /// <summary>
        /// Test the Vegas method with the Pi function
        /// </summary>
        [TestMethod()]
        public void Test_PI()
        {
            var vegas = new Vegas((x, y) => { return Integrands.PI(x); }, 2, new double[] { -1, -1 }, new double[] { 1, 1 });
            vegas.Integrate();
            var result = vegas.Result;
            double trueResult = 3.1416;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Vegas method with the GSL function
        /// </summary>
        [TestMethod()]
        public void Test_GSL()
        {
            var vegas = new Vegas((x, y) => { return Integrands.GSL(x); }, 3, new double[] { 0, 0, 0 }, new double[] { Math.PI, Math.PI, Math.PI });
            vegas.Integrate();
            var result = vegas.Result;
            double trueResult = 1.3932039296856768591842462603255;
            Assert.AreEqual(trueResult, result, 1E-2 * trueResult);
        }

        /// <summary>
        /// Test the Vegas method with the sum of three normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfThreeNormals()
        {
            var min = new double[3];
            var max = new double[3];
            for (int i = 0; i < 3; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var vegas = new Vegas((x, y) => { return Integrands.SumOfNormals(x); }, 3, min, max);
            vegas.Integrate();
            var result = vegas.Result;
            double trueResult = 57;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Vegas method with the sum of five normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfFiveNormals()
        {
            var min = new double[5];
            var max = new double[5];
            for (int i = 0; i < 5; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var vegas = new Vegas((x, y) => { return Integrands.SumOfNormals(x); }, 5, min, max);
            vegas.Integrate();
            var result = vegas.Result;
            double trueResult = 224;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }

        /// <summary>
        /// Test the Vegas method with the sum of 20 normal distributions
        /// </summary>
        [TestMethod()]
        public void Test_SumOfTwentyNormals()
        {
            var min = new double[20];
            var max = new double[20];
            for (int i = 0; i < 20; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            var vegas = new Vegas((x, y) => { return Integrands.SumOfNormals(x); }, 20, min, max);
            vegas.Integrate();
            var result = vegas.Result;
            double trueResult = 837;
            Assert.AreEqual(trueResult, result, 1E-3 * trueResult);
        }


        /// <summary>
        /// Test 1: Verify backward compatibility - Power Transform with γ=1 should match standard Vegas
        /// </summary>
        [TestMethod()]
        public void Test_PowerTransform_BackwardCompatibility()
        {
            var min = new double[5];
            var max = new double[5];
            for (int i = 0; i < 5; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            // Standard Vegas (no transform)
            var vegasStandard = new Vegas((x, y) => { return Integrands.SumOfNormals(x); }, 5, min, max);
            vegasStandard.Integrate();
            var resultStandard = vegasStandard.Result;

            // Power Transform with γ=1.0 (should be identical)
            var vegasPower = new Vegas((x, y) => { return Integrands.SumOfNormals(x); }, 5, min, max);
            vegasPower.TailFocusParameter = 1.0;  // No transformation
            vegasPower.Integrate();
            var resultPower = vegasPower.Result;

            // Results should be essentially the same
            double trueResult = 224;
            Assert.AreEqual(trueResult, resultStandard, 1E-2 * trueResult, "Standard Vegas failed");
            Assert.AreEqual(trueResult, resultPower, 1E-2 * trueResult, "Power Transform with γ=1 failed");
            Assert.AreEqual(resultStandard, resultPower, 1E-2 * trueResult, "Results should match");
        }

        /// <summary>
        /// Test 2: Find probability of rare upper tail event P(Sum > threshold)
        /// This tests the core rare event capability of Power Transform
        /// </summary>
        [TestMethod()]
        public void Test_PowerTransform_RareUpperTailEvent()
        {
            var min = new double[5];
            var max = new double[5];
            for (int i = 0; i < 5; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            // Define a high threshold for rare event
            // For sum of 5 normals: mean = 224, std dev = sqrt(2² + 15² + 5² + 14² + 7²) = sqrt(499) ≈ 22.34
            // threshold = mean + 3*std ≈ 224 + 3*22.34 ≈ 291.0
            // P(Sum > 291.0) = P(Z > 3.0) ≈ 0.00135 (about 1.35e-3)
            double threshold = 291.0;

            var vegas = new Vegas(
                (x, y) =>
                {
                    double sum = Integrands.SumOfNormals(x);
                    // Return weight if in failure region, 0 otherwise
                    return (sum > threshold) ? 1.0 : 0.0;
                },
                5, min, max
            );

            // Configure for rare event with moderate tail focus
            vegas.TailFocusParameter = 2.0;  // Moderate focus
            vegas.NumberOfBins = 100;
            vegas.FunctionCalls = 50000;
            vegas.IndependentEvaluations = 5;
            vegas.Alpha = 1.8;

            vegas.Integrate();
            var failureProbability = vegas.Result;

            // Analytical approximation: P(Sum > mean + 3σ) = P(Z > 3) ≈ 0.00135
            // We expect something in the ballpark of 1e-3 to 2e-3
            Assert.IsGreaterThan(5E-4, failureProbability, $"Probability too small: {failureProbability:E6}");
            Assert.IsLessThan(5E-3, failureProbability, $"Probability too large: {failureProbability:E6}");

            // Standard error should be reasonable (less than 50% of estimate)
            double relativeError = vegas.StandardError / Math.Abs(failureProbability);
            Assert.IsLessThan(0.5, relativeError, $"Relative error too large: {relativeError:P1}");
        }

        /// <summary>
        /// Test 3: Adaptive configuration for very rare events (p ~ 1e-6)
        /// Tests the ConfigureForRareEvents() helper method
        /// </summary>
        [TestMethod()]
        public void Test_PowerTransform_VeryRareEvent()
        {
            var min = new double[3];
            var max = new double[3];
            for (int i = 0; i < 3; i++)
            {
                min[i] = 1E-16;
                max[i] = 1 - 1E-16;
            }

            // For sum of 3 normals: mean = 57, std dev ≈ sqrt(2² + 15² + 5²) ≈ 15.9
            // threshold = mean + 4.5*std ≈ 57 + 4.5*15.9 ≈ 128.5
            // P(Sum > 128.5) ≈ Normal.CDF(-4.5) ≈ 3.4e-6
            double threshold = 128.5;

            var vegas = new Vegas(
                (x, y) =>
                {
                    double sum = Integrands.SumOfNormals(x);
                    return (sum > threshold) ? 1.0 : 0.0;
                },
                3, min, max
            );

            // Use automatic configuration for very rare event
            vegas.ConfigureForRareEvents(targetProbability: 1e-5);
            vegas.FunctionCalls = 100000;  // More samples for very rare event
            vegas.IndependentEvaluations = 5;

            vegas.Integrate();
            var failureProbability = vegas.Result;

            // Should be in ballpark of 1e-6 to 1e-5
            Assert.IsGreaterThan(1E-7, failureProbability, $"Probability too small: {failureProbability:E2}");
            Assert.IsLessThan(1E-4,failureProbability, $"Probability too large: {failureProbability:E2}");

            // For very rare events, relative error can be higher but should be finite
            double relativeError = vegas.StandardError / Math.Abs(failureProbability);
            Assert.IsLessThan(1.0, relativeError, $"Relative error too large: {relativeError:P1}");
            Assert.IsFalse(double.IsNaN(failureProbability), "Result should not be NaN");
            Assert.IsFalse(double.IsInfinity(failureProbability), "Result should not be infinite");
        }

        /// <summary>
        /// Test 4: BONUS - Verify that integrand receives correct probabilities
        /// This ensures the power transform doesn't break the probability input range
        /// </summary>
        [TestMethod()]
        public void Test_PowerTransform_ProbabilityRange()
        {
            var min = new double[3];
            var max = new double[3];
            for (int i = 0; i < 3; i++)
            {
                min[i] = 0.0;
                max[i] = 1.0;
            }

            double minObserved = 1.0;
            double maxObserved = 0.0;
            int sampleCount = 0;

            var vegas = new Vegas(
                (x, y) =>
                {
                    sampleCount++;
                    // Track min/max observed probabilities
                    foreach (var p in x)
                    {
                        minObserved = Math.Min(minObserved, p);
                        maxObserved = Math.Max(maxObserved, p);

                        // Verify probabilities are in valid range
                        Assert.IsTrue(p >= 0.0 && p <= 1.0,
                            $"Probability out of range: {p}");
                    }
                    return Integrands.SumOfNormals(x);
                },
                3, min, max
            );

            // Use strong tail focus
            vegas.TailFocusParameter = 4.0;
            vegas.FunctionCalls = 10000;
            vegas.IndependentEvaluations = 2;

            vegas.Integrate();

            // Verify samples were generated
            Assert.IsGreaterThan(0,sampleCount , "No samples generated");

            // With γ=4, should see strong tail focus (many samples near 1.0)
            Assert.IsGreaterThan(0.99,maxObserved, $"Max probability too low: {maxObserved}");

            // Should still have some diversity (not all at 1.0)
            Assert.IsLessThan(0.5,minObserved, $"Min probability too high: {minObserved}");

        }

    }

}
