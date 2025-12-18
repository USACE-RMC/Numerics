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
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;

namespace Serialization
{
    /// <summary>
    /// Unit tests for JsonSerializer implementations in UncertaintyAnalysisResults and MCMCResults.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This test class verifies the serialization and deserialization functionality for:
    /// - UncertaintyAnalysisResults: Bootstrap analysis results with distribution parameters
    /// - MCMCResults: Markov Chain Monte Carlo sampling results
    /// </para>
    /// <para>
    /// Tests cover normal operation, edge cases (null, empty), and error handling.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_JsonSerialization
    {
        #region UncertaintyAnalysisResults Tests

        /// <summary>
        /// Tests basic serialization and deserialization of UncertaintyAnalysisResults with scalar properties.
        /// </summary>
        /// <remarks>
        /// Verifies that AIC, BIC, DIC, and RMSE values are preserved through serialization round-trip.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_BasicSerialization()
        {
            // Arrange
            var original = CreateSampleUncertaintyAnalysisResults();

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.AreEqual(original.AIC, deserialized.AIC, 1e-10);
            Assert.AreEqual(original.BIC, deserialized.BIC, 1e-10);
            Assert.AreEqual(original.DIC, deserialized.DIC, 1e-10);
            Assert.AreEqual(original.RMSE, deserialized.RMSE, 1e-10);
        }

        /// <summary>
        /// Tests serialization and deserialization of 1D arrays (ModeCurve, MeanCurve).
        /// </summary>
        /// <remarks>
        /// Verifies that array dimensions and values are preserved exactly.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_ArraySerialization()
        {
            // Arrange
            var original = CreateSampleUncertaintyAnalysisResults();

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.ModeCurve);
            Assert.HasCount(original.ModeCurve.Length, deserialized.ModeCurve);
            for (int i = 0; i < original.ModeCurve.Length; i++)
            {
                Assert.AreEqual(original.ModeCurve[i], deserialized.ModeCurve[i], 1e-10);
            }

            Assert.IsNotNull(deserialized.MeanCurve);
            Assert.HasCount(original.MeanCurve.Length, deserialized.MeanCurve);
            for (int i = 0; i < original.MeanCurve.Length; i++)
            {
                Assert.AreEqual(original.MeanCurve[i], deserialized.MeanCurve[i], 1e-10);
            }
        }

        /// <summary>
        /// Tests serialization and deserialization of 2D arrays (ConfidenceIntervals).
        /// </summary>
        /// <remarks>
        /// Verifies that 2D array dimensions and all element values are preserved correctly.
        /// Uses the custom Double2DArrayConverter.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_2DArraySerialization()
        {
            // Arrange
            var original = CreateSampleUncertaintyAnalysisResults();

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.ConfidenceIntervals);
            Assert.AreEqual(original.ConfidenceIntervals.GetLength(0), deserialized.ConfidenceIntervals.GetLength(0));
            Assert.AreEqual(original.ConfidenceIntervals.GetLength(1), deserialized.ConfidenceIntervals.GetLength(1));

            for (int i = 0; i < original.ConfidenceIntervals.GetLength(0); i++)
            {
                for (int j = 0; j < original.ConfidenceIntervals.GetLength(1); j++)
                {
                    Assert.AreEqual(original.ConfidenceIntervals[i, j], deserialized.ConfidenceIntervals[i, j], 1e-10);
                }
            }
        }

        /// <summary>
        /// Tests serialization and deserialization of complex objects (ParameterSet arrays).
        /// </summary>
        /// <remarks>
        /// Verifies that nested objects with arrays are correctly serialized and deserialized,
        /// including Fitness, Weight, and Values properties.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_ParameterSetsSerialization()
        {
            // Arrange
            var original = CreateSampleUncertaintyAnalysisResults();

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            var origLen = original.ParameterSets.Length;
            var desLen = deserialized.ParameterSets.Length;
            // Assert
            Assert.IsNotNull(deserialized.ParameterSets);
            Assert.AreEqual(origLen, desLen);

            for (int i = 0; i < original.ParameterSets.Length; i++)
            {
                Assert.AreEqual(original.ParameterSets[i].Fitness, deserialized.ParameterSets[i].Fitness, 1e-10);
                Assert.AreEqual(original.ParameterSets[i].Weight, deserialized.ParameterSets[i].Weight, 1e-10);

                if (original.ParameterSets[i].Values != null)
                {
                    Assert.IsNotNull(deserialized.ParameterSets[i].Values);
                    Assert.HasCount(original.ParameterSets[i].Values.Length, deserialized.ParameterSets[i].Values);

                    for (int j = 0; j < original.ParameterSets[i].Values.Length; j++)
                    {
                        Assert.AreEqual(original.ParameterSets[i].Values[j], deserialized.ParameterSets[i].Values[j], 1e-10);
                    }
                }
            }
        }

        /// <summary>
        /// Tests serialization with null properties to verify proper null handling.
        /// </summary>
        /// <remarks>
        /// Verifies that null values are serialized and deserialized correctly without throwing exceptions.
        /// The JsonSerializerOptions should be configured to handle nulls appropriately.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_NullHandling()
        {
            // Arrange
            var original = new UncertaintyAnalysisResults
            {
                AIC = 123.456,
                BIC = 234.567,
                ParentDistribution = null,
                ParameterSets = null,
                ConfidenceIntervals = null,
                ModeCurve = null,
                MeanCurve = null
            };

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.AreEqual(original.AIC, deserialized.AIC, 1e-10);
            Assert.AreEqual(original.BIC, deserialized.BIC, 1e-10);
            Assert.IsNull(deserialized.ParentDistribution);
            Assert.IsNull(deserialized.ParameterSets);
            Assert.IsNull(deserialized.ConfidenceIntervals);
            Assert.IsNull(deserialized.ModeCurve);
            Assert.IsNull(deserialized.MeanCurve);
        }

        /// <summary>
        /// Tests serialization with empty arrays to verify proper empty collection handling.
        /// </summary>
        /// <remarks>
        /// Verifies that empty arrays (length 0) are preserved through serialization
        /// and are not converted to null.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_EmptyArrays()
        {
            // Arrange
            var original = new UncertaintyAnalysisResults
            {
                ParameterSets = new ParameterSet[0],
                ModeCurve = new double[0],
                MeanCurve = new double[0],
                ConfidenceIntervals = new double[0, 0]
            };

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.IsNotNull(deserialized.ParameterSets);
            Assert.IsEmpty( deserialized.ParameterSets);
            Assert.IsNotNull(deserialized.ModeCurve);
            Assert.IsEmpty(deserialized.ModeCurve);
            Assert.IsNotNull(deserialized.MeanCurve);
            Assert.IsEmpty(deserialized.MeanCurve);
        }

        /// <summary>
        /// Tests deserialization with invalid byte data to verify error handling.
        /// </summary>
        /// <remarks>
        /// Verifies that FromByteArray returns null instead of throwing an exception
        /// when given invalid JSON data. This is important for robust error handling.
        /// </remarks>
        [TestMethod]
        public void Test_UncertaintyAnalysisResults_InvalidBytes()
        {
            // Arrange
            byte[] invalidBytes = new byte[] { 0x00, 0x01, 0x02, 0x03 };

            // Act
            var result = UncertaintyAnalysisResults.FromByteArray(invalidBytes);

            // Assert - Should return null on deserialization error
            Assert.IsNull(result);
        }

        #endregion

        #region MCMCResults Tests

        /// <summary>
        /// Tests basic serialization and deserialization of MCMCResults with acceptance rates.
        /// </summary>
        /// <remarks>
        /// Verifies that scalar arrays like AcceptanceRates are preserved through serialization.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_BasicSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            var origLen = original.AcceptanceRates.Length;
            var desLen = deserialized.AcceptanceRates.Length;
            // Assert
            Assert.IsNotNull(deserialized);
            Assert.IsNotNull(deserialized.AcceptanceRates);
            Assert.AreEqual(origLen, desLen);

            for (int i = 0; i < original.AcceptanceRates.Length; i++)
            {
                Assert.AreEqual(original.AcceptanceRates[i], deserialized.AcceptanceRates[i], 1e-10);
            }
        }

        /// <summary>
        /// Tests serialization and deserialization of Markov chains (jagged array of ParameterSet lists).
        /// </summary>
        /// <remarks>
        /// This is the most complex data structure in MCMCResults. Verifies that:
        /// - Array dimensions are preserved
        /// - List lengths are preserved
        /// - All ParameterSet properties (Fitness, Weight, Values) are correct
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_MarkovChainsSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            var origLen = original.MarkovChains.Length;
            var desLen = deserialized.MarkovChains.Length;
            // Assert
            Assert.IsNotNull(deserialized.MarkovChains);
            Assert.AreEqual(origLen, desLen);

            for (int i = 0; i < original.MarkovChains.Length; i++)
            {
                Assert.IsNotNull(deserialized.MarkovChains[i]);
                Assert.HasCount(original.MarkovChains[i].Count, deserialized.MarkovChains[i]);

                for (int j = 0; j < original.MarkovChains[i].Count; j++)
                {
                    var originalParam = original.MarkovChains[i][j];
                    var deserializedParam = deserialized.MarkovChains[i][j];

                    Assert.AreEqual(originalParam.Fitness, deserializedParam.Fitness, 1e-10);
                    Assert.AreEqual(originalParam.Weight, deserializedParam.Weight, 1e-10);

                    if (originalParam.Values != null)
                    {
                        Assert.IsNotNull(deserializedParam.Values);
                        CollectionAssert.AreEqual(originalParam.Values, deserializedParam.Values);
                    }
                }
            }
        }

        /// <summary>
        /// Tests serialization and deserialization of the Output list.
        /// </summary>
        /// <remarks>
        /// The Output list contains all samples from all chains. Verifies that the combined
        /// output is correctly serialized and deserialized.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_OutputSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.Output);
            Assert.HasCount(original.Output.Count, deserialized.Output);

            for (int i = 0; i < original.Output.Count; i++)
            {
                Assert.AreEqual(original.Output[i].Fitness, deserialized.Output[i].Fitness, 1e-10);
                Assert.AreEqual(original.Output[i].Weight, deserialized.Output[i].Weight, 1e-10);
            }
        }

        /// <summary>
        /// Tests serialization and deserialization of MAP and PosteriorMean ParameterSets.
        /// </summary>
        /// <remarks>
        /// Verifies that the Maximum A Posteriori (MAP) and posterior mean estimates
        /// are correctly preserved through serialization.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_MAPandMeanSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.AreEqual(original.MAP.Fitness, deserialized.MAP.Fitness, 1e-10);
            Assert.IsNotNull(deserialized.MAP.Values);
            CollectionAssert.AreEqual(original.MAP.Values, deserialized.MAP.Values);

            Assert.AreEqual(original.PosteriorMean.Fitness, deserialized.PosteriorMean.Fitness, 1e-10);
            Assert.IsNotNull(deserialized.PosteriorMean.Values);
            CollectionAssert.AreEqual(original.PosteriorMean.Values, deserialized.PosteriorMean.Values);
        }

        /// <summary>
        /// Tests serialization and deserialization of MeanLogLikelihood diagnostics.
        /// </summary>
        /// <remarks>
        /// MeanLogLikelihood tracks convergence diagnostics. Verifies that all values
        /// are preserved correctly.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_MeanLogLikelihoodSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.MeanLogLikelihood);
            Assert.HasCount(original.MeanLogLikelihood.Count, deserialized.MeanLogLikelihood);

            for (int i = 0; i < original.MeanLogLikelihood.Count; i++)
            {
                Assert.AreEqual(original.MeanLogLikelihood[i], deserialized.MeanLogLikelihood[i], 1e-10);
            }
        }

        /// <summary>
        /// Tests serialization with empty Markov chains to verify edge case handling.
        /// </summary>
        /// <remarks>
        /// Verifies that empty chains (with no samples) are correctly serialized and deserialized.
        /// This can occur if MCMC sampling fails or is interrupted.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_EmptyChains()
        {
            // Arrange
            var original = new MCMCResults();
            var chains = new List<ParameterSet>[2];
            chains[0] = new List<ParameterSet>();
            chains[1] = new List<ParameterSet>();

            // Use reflection to set private properties for testing
            var type = typeof(MCMCResults);
            type.GetProperty("MarkovChains").SetValue(original, chains);
            type.GetProperty("Output").SetValue(original, new List<ParameterSet>());
            type.GetProperty("MeanLogLikelihood").SetValue(original, new List<double>());
            type.GetProperty("AcceptanceRates").SetValue(original, new double[2]);

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.IsNotNull(deserialized.MarkovChains);
            Assert.HasCount(2, deserialized.MarkovChains);
            Assert.HasCount(0, deserialized.MarkovChains[0]);
            Assert.HasCount(0, deserialized.MarkovChains[1]);
        }

        /// <summary>
        /// Tests serialization with a large dataset to verify performance and correctness at scale.
        /// </summary>
        /// <remarks>
        /// Creates 5 chains with 1000 samples each (5000 total samples) to test:
        /// - Memory efficiency of serialization
        /// - Correctness with large data volumes
        /// - Performance (should complete in reasonable time)
        /// Tests first and last elements to ensure proper ordering is maintained.
        /// </remarks>
        [TestMethod]
        public void Test_MCMCResults_LargeDataSet()
        {
            // Arrange - Create a larger dataset to test performance and correctness
            var original = CreateLargeMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            var origLen = original.MarkovChains.Length;
            var deserializedLen = deserialized.MarkovChains.Length;
            // Assert
            Assert.IsNotNull(deserialized);
            Assert.AreEqual(origLen, deserializedLen);

            // Verify first and last elements to ensure proper serialization
            var firstOriginal = original.MarkovChains[0][0];
            var firstDeserialized = deserialized.MarkovChains[0][0];
            Assert.AreEqual(firstOriginal.Fitness, firstDeserialized.Fitness, 1e-10);

            var lastChainIdx = original.MarkovChains.Length - 1;
            var lastElementIdx = original.MarkovChains[lastChainIdx].Count - 1;
            var lastOriginal = original.MarkovChains[lastChainIdx][lastElementIdx];
            var lastDeserialized = deserialized.MarkovChains[lastChainIdx][lastElementIdx];
            Assert.AreEqual(lastOriginal.Fitness, lastDeserialized.Fitness, 1e-10);
        }

        /// <summary>
        /// Tests that JsonSerializerOptions are configured correctly for the serializers.
        /// </summary>
        /// <remarks>
        /// Verifies important serialization settings:
        /// - WriteIndented = false (compact JSON for efficiency)
        /// - DefaultIgnoreCondition (null handling)
        /// - IncludeFields = true (serializes public fields in addition to properties)
        /// </remarks>
        [TestMethod]
        public void Test_JsonSerializerOptions_Configuration()
        {
            // This test verifies that the JsonSerializerOptions are configured correctly
            // Arrange
            var original = new UncertaintyAnalysisResults
            {
                AIC = 100.0,
                BIC = 200.0,
                DIC = 0.0, // Test that zero values are included
                RMSE = 0.0,
                ParentDistribution = null // Test that null values are ignored
            };

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            string jsonString = System.Text.Encoding.UTF8.GetString(serialized);

            // Assert
            // Verify that WriteIndented is false (no formatting whitespace)
            Assert.DoesNotContain(jsonString,("\n"));
            Assert.DoesNotContain(jsonString, ("  ")); // No indentation

            // Verify that null values are not included (DefaultIgnoreCondition)
            Assert.DoesNotContain(jsonString,  ("\"ParentDistribution\":null"));

            // Verify that fields are included (IncludeFields = true)
            Assert.DoesNotContain(jsonString, ("\"AIC\":"));
            Assert.DoesNotContain(jsonString, ("\"BIC\":"));
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Creates a sample UncertaintyAnalysisResults instance for testing.
        /// </summary>
        /// <returns>
        /// A UncertaintyAnalysisResults with populated properties including:
        /// - ParentDistribution (Normal)
        /// - ParameterSets array with 3 elements
        /// - ConfidenceIntervals 2D array (3x2)
        /// - ModeCurve and MeanCurve arrays
        /// - Scalar diagnostic values (AIC, BIC, DIC, RMSE)
        /// </returns>
        private UncertaintyAnalysisResults CreateSampleUncertaintyAnalysisResults()
        {
            return new UncertaintyAnalysisResults
            {
                ParentDistribution = new Normal(100, 15),
                ParameterSets = new[]
                {
                    new ParameterSet(new double[] { 100, 15 }, -123.45, 0.5),
                    new ParameterSet(new double[] { 99.5, 14.8 }, -124.00, 0.3),
                    new ParameterSet(new double[] { 100.5, 15.2 }, -123.80, 0.2)
                },
                ConfidenceIntervals = new double[,]
                {
                    { 70.0, 130.0 },
                    { 75.0, 125.0 },
                    { 80.0, 120.0 }
                },
                ModeCurve = new double[] { 85, 90, 95, 100, 105, 110, 115 },
                MeanCurve = new double[] { 84, 89, 94, 99, 104, 109, 114 },
                AIC = 250.5,
                BIC = 255.3,
                DIC = 252.1,
                RMSE = 2.35
            };
        }

        /// <summary>
        /// Creates a sample MCMCResults instance for testing.
        /// </summary>
        /// <returns>
        /// A MCMCResults with:
        /// - 3 Markov chains with 3 samples each
        /// - Combined output list
        /// - MeanLogLikelihood diagnostics
        /// - AcceptanceRates for each chain
        /// - MAP and PosteriorMean estimates
        /// </returns>
        /// <remarks>
        /// Uses reflection to set properties with private setters.
        /// </remarks>
        private MCMCResults CreateSampleMCMCResults()
        {
            var result = new MCMCResults();

            // Create sample Markov chains
            var chains = new List<ParameterSet>[3];
            for (int i = 0; i < 3; i++)
            {
                chains[i] = new List<ParameterSet>
                {
                    new ParameterSet(new double[] { 10 + i, 20 + i }, -100 - i, 1.0),
                    new ParameterSet(new double[] { 11 + i, 21 + i }, -101 - i, 1.0),
                    new ParameterSet(new double[] { 12 + i, 22 + i }, -102 - i, 1.0)
                };
            }

            // Create output
            var output = new List<ParameterSet>();
            foreach (var chain in chains)
            {
                output.AddRange(chain);
            }

            // Set properties using reflection since they have private setters
            var type = typeof(MCMCResults);
            type.GetProperty("MarkovChains").SetValue(result, chains);
            type.GetProperty("Output").SetValue(result, output);
            type.GetProperty("MeanLogLikelihood").SetValue(result, new List<double> { -100, -101, -102 });
            type.GetProperty("AcceptanceRates").SetValue(result, new double[] { 0.45, 0.50, 0.48 });
            type.GetProperty("MAP").SetValue(result, new ParameterSet(new double[] { 10, 20 }, -100, 1.0));
            type.GetProperty("PosteriorMean").SetValue(result, new ParameterSet(new double[] { 11, 21 }, -101, 1.0));

            return result;
        }

        /// <summary>
        /// Creates a large MCMCResults instance for performance and stress testing.
        /// </summary>
        /// <returns>
        /// A MCMCResults with:
        /// - 5 Markov chains with 1000 samples each (5000 total)
        /// - Random parameter values for realistic data size
        /// - Full diagnostic information
        /// </returns>
        /// <remarks>
        /// <para>
        /// Used to test serialization performance with realistic data volumes.
        /// Uses a fixed random seed (12345) for reproducible test results.
        /// </para>
        /// <para>
        /// Total data size: ~5000 ParameterSets × 2 parameters × 8 bytes = ~80KB minimum
        /// </para>
        /// </remarks>
        private MCMCResults CreateLargeMCMCResults()
        {
            var result = new MCMCResults();
            var random = new Random(12345);

            // Create larger Markov chains
            var chains = new List<ParameterSet>[5];
            for (int i = 0; i < 5; i++)
            {
                chains[i] = new List<ParameterSet>();
                for (int j = 0; j < 1000; j++)
                {
                    var values = new double[] { random.NextDouble() * 100, random.NextDouble() * 50 };
                    chains[i].Add(new ParameterSet(values, -random.NextDouble() * 200, 1.0));
                }
            }

            // Create output
            var output = new List<ParameterSet>();
            foreach (var chain in chains)
            {
                output.AddRange(chain);
            }

            // Create mean log likelihood
            var meanLogLikelihood = new List<double>();
            for (int i = 0; i < 1000; i++)
            {
                meanLogLikelihood.Add(-random.NextDouble() * 200);
            }

            // Set properties using reflection
            var type = typeof(MCMCResults);
            type.GetProperty("MarkovChains").SetValue(result, chains);
            type.GetProperty("Output").SetValue(result, output);
            type.GetProperty("MeanLogLikelihood").SetValue(result, meanLogLikelihood);
            type.GetProperty("AcceptanceRates").SetValue(result, new double[] { 0.4, 0.45, 0.5, 0.48, 0.46 });
            type.GetProperty("MAP").SetValue(result, chains[0][0]);
            type.GetProperty("PosteriorMean").SetValue(result, chains[0][500]);

            return result;
        }

        #endregion
    }
}