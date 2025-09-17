using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling.MCMC;

namespace Test_Numerics.Serialization
{
    /// <summary>
    /// Unit tests for JsonSerializer implementations in UncertaintyAnalysisResults and MCMCResults.
    /// </summary>
    [TestClass]
    public class Test_JsonSerialization
    {
        #region UncertaintyAnalysisResults Tests

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
            Assert.AreEqual(original.ModeCurve.Length, deserialized.ModeCurve.Length);
            for (int i = 0; i < original.ModeCurve.Length; i++)
            {
                Assert.AreEqual(original.ModeCurve[i], deserialized.ModeCurve[i], 1e-10);
            }

            Assert.IsNotNull(deserialized.MeanCurve);
            Assert.AreEqual(original.MeanCurve.Length, deserialized.MeanCurve.Length);
            for (int i = 0; i < original.MeanCurve.Length; i++)
            {
                Assert.AreEqual(original.MeanCurve[i], deserialized.MeanCurve[i], 1e-10);
            }
        }

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

        [TestMethod]
        public void Test_UncertaintyAnalysisResults_ParameterSetsSerialization()
        {
            // Arrange
            var original = CreateSampleUncertaintyAnalysisResults();

            // Act
            byte[] serialized = UncertaintyAnalysisResults.ToByteArray(original);
            var deserialized = UncertaintyAnalysisResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.ParameterSets);
            Assert.AreEqual(original.ParameterSets.Length, deserialized.ParameterSets.Length);

            for (int i = 0; i < original.ParameterSets.Length; i++)
            {
                Assert.AreEqual(original.ParameterSets[i].Fitness, deserialized.ParameterSets[i].Fitness, 1e-10);
                Assert.AreEqual(original.ParameterSets[i].Weight, deserialized.ParameterSets[i].Weight, 1e-10);

                if (original.ParameterSets[i].Values != null)
                {
                    Assert.IsNotNull(deserialized.ParameterSets[i].Values);
                    Assert.AreEqual(original.ParameterSets[i].Values.Length, deserialized.ParameterSets[i].Values.Length);

                    for (int j = 0; j < original.ParameterSets[i].Values.Length; j++)
                    {
                        Assert.AreEqual(original.ParameterSets[i].Values[j], deserialized.ParameterSets[i].Values[j], 1e-10);
                    }
                }
            }
        }

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
            Assert.AreEqual(0, deserialized.ParameterSets.Length);
            Assert.IsNotNull(deserialized.ModeCurve);
            Assert.AreEqual(0, deserialized.ModeCurve.Length);
            Assert.IsNotNull(deserialized.MeanCurve);
            Assert.AreEqual(0, deserialized.MeanCurve.Length);
        }

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

        [TestMethod]
        public void Test_MCMCResults_BasicSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.IsNotNull(deserialized.AcceptanceRates);
            Assert.AreEqual(original.AcceptanceRates.Length, deserialized.AcceptanceRates.Length);

            for (int i = 0; i < original.AcceptanceRates.Length; i++)
            {
                Assert.AreEqual(original.AcceptanceRates[i], deserialized.AcceptanceRates[i], 1e-10);
            }
        }

        [TestMethod]
        public void Test_MCMCResults_MarkovChainsSerialization()
        {
            // Arrange
            var original = CreateSampleMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized.MarkovChains);
            Assert.AreEqual(original.MarkovChains.Length, deserialized.MarkovChains.Length);

            for (int i = 0; i < original.MarkovChains.Length; i++)
            {
                Assert.IsNotNull(deserialized.MarkovChains[i]);
                Assert.AreEqual(original.MarkovChains[i].Count, deserialized.MarkovChains[i].Count);

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
            Assert.AreEqual(original.Output.Count, deserialized.Output.Count);

            for (int i = 0; i < original.Output.Count; i++)
            {
                Assert.AreEqual(original.Output[i].Fitness, deserialized.Output[i].Fitness, 1e-10);
                Assert.AreEqual(original.Output[i].Weight, deserialized.Output[i].Weight, 1e-10);
            }
        }

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
            Assert.AreEqual(original.MeanLogLikelihood.Count, deserialized.MeanLogLikelihood.Count);

            for (int i = 0; i < original.MeanLogLikelihood.Count; i++)
            {
                Assert.AreEqual(original.MeanLogLikelihood[i], deserialized.MeanLogLikelihood[i], 1e-10);
            }
        }

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
            Assert.AreEqual(2, deserialized.MarkovChains.Length);
            Assert.AreEqual(0, deserialized.MarkovChains[0].Count);
            Assert.AreEqual(0, deserialized.MarkovChains[1].Count);
        }

        [TestMethod]
        public void Test_MCMCResults_LargeDataSet()
        {
            // Arrange - Create a larger dataset to test performance and correctness
            var original = CreateLargeMCMCResults();

            // Act
            byte[] serialized = MCMCResults.ToByteArray(original);
            var deserialized = MCMCResults.FromByteArray(serialized);

            // Assert
            Assert.IsNotNull(deserialized);
            Assert.AreEqual(original.MarkovChains.Length, deserialized.MarkovChains.Length);

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
            Assert.IsFalse(jsonString.Contains("\n"));
            Assert.IsFalse(jsonString.Contains("  ")); // No indentation

            // Verify that null values are not included (DefaultIgnoreCondition)
            Assert.IsFalse(jsonString.Contains("\"ParentDistribution\":null"));

            // Verify that fields are included (IncludeFields = true)
            Assert.IsTrue(jsonString.Contains("\"AIC\":"));
            Assert.IsTrue(jsonString.Contains("\"BIC\":"));
        }

        #endregion

        #region Helper Methods

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