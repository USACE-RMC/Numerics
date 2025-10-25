using System;
using System.Text.Json;
using Numerics.Distributions;
using Numerics.Utilities;

namespace Test_Numerics.Serialization
{
    /// <summary>
    /// Demonstration of custom JSON converters for complex types.
    /// </summary>
    [TestClass]
    public static class JsonConverterDemo
    {
        /// <summary>
        /// Demonstrates how the custom converters handle 2D arrays and complex distribution objects.
        /// </summary>
        /// 
        [TestMethod]
        public static void RunDemo()
        {
            Console.WriteLine("=== JSON Converter Demo ===\n");

            // Demo 1: 2D Double Array
            Console.WriteLine("1. Double[,] Array Serialization:");
            double[,] matrix = new double[,]
            {
                { 1.1, 2.2, 3.3 },
                { 4.4, 5.5, 6.6 }
            };

            var options = new JsonSerializerOptions
            {
                WriteIndented = true,
                Converters = { new Double2DArrayConverter() }
            };
            options.Converters.Add(new UnivariateDistributionConverter());

            string json = JsonSerializer.Serialize(matrix, options);
            Console.WriteLine("Serialized 2D array:");
            Console.WriteLine(json);

            double[,] deserializedMatrix = JsonSerializer.Deserialize<double[,]>(json, options);
            Console.WriteLine($"\nDeserialized successfully: {deserializedMatrix[0, 0]} == {matrix[0, 0]}");

            // Demo 2: UnivariateDistribution
            Console.WriteLine("\n2. UnivariateDistribution Serialization:");
            var normalDist = new Normal(100, 15);

            json = JsonSerializer.Serialize<UnivariateDistributionBase>(normalDist, options);
            Console.WriteLine("Serialized Normal distribution:");
            Console.WriteLine(json);

            var deserializedDist = JsonSerializer.Deserialize<UnivariateDistributionBase>(json, options);
            Console.WriteLine($"\nDeserialized successfully: Type = {deserializedDist?.Type}");
            Console.WriteLine($"Parameters match: Mean={deserializedDist?.Mean:F2}, StdDev={deserializedDist?.StandardDeviation:F2}");

            // Demo 3: Full UncertaintyAnalysisResults
            Console.WriteLine("\n3. UncertaintyAnalysisResults with all complex types:");
            var results = new UncertaintyAnalysisResults
            {
                ParentDistribution = new Normal(50, 10),
                ConfidenceIntervals = new double[,]
                {
                    { 30.0, 70.0 },
                    { 35.0, 65.0 },
                    { 40.0, 60.0 }
                },
                AIC = 123.45,
                BIC = 234.56
            };

            // Using the built-in serialization methods that now include converters
            byte[] bytes = UncertaintyAnalysisResults.ToByteArray(results);
            var deserializedResults = UncertaintyAnalysisResults.FromByteArray(bytes);

            Console.WriteLine($"Full object serialized and deserialized successfully!");
            Console.WriteLine($"AIC: {deserializedResults?.AIC}");
            Console.WriteLine($"BIC: {deserializedResults?.BIC}");
            Console.WriteLine($"Parent Distribution Type: {deserializedResults?.ParentDistribution?.Type}");
            Console.WriteLine($"Confidence Intervals shape: [{deserializedResults?.ConfidenceIntervals?.GetLength(0)}, {deserializedResults?.ConfidenceIntervals?.GetLength(1)}]");

            Console.WriteLine("\n=== Demo Complete ===");
        }

        /// <summary>
        /// Shows how the converters handle edge cases.
        /// </summary>
        public static void RunEdgeCasesDemo()
        {
            Console.WriteLine("=== Edge Cases Demo ===\n");

            var options = new JsonSerializerOptions
            {
                WriteIndented = true,
                Converters =
                {
                    new Double2DArrayConverter(),
                    new String2DArrayConverter(),
                    new UnivariateDistributionConverter()
                }
            };

            // Empty 2D array
            double[,] emptyMatrix = new double[0, 0];
            string json = JsonSerializer.Serialize(emptyMatrix, options);
            Console.WriteLine("Empty 2D array:");
            Console.WriteLine(json);

            // Null distribution
            UnivariateDistributionBase nullDist = null;
            json = JsonSerializer.Serialize(nullDist, options);
            Console.WriteLine("\nNull distribution:");
            Console.WriteLine(json);

            // Large 2D array
            double[,] largeMatrix = new double[100, 50];
            Random rand = new Random(42);
            for (int i = 0; i < 100; i++)
            {
                for (int j = 0; j < 50; j++)
                {
                    largeMatrix[i, j] = rand.NextDouble() * 1000;
                }
            }

            json = JsonSerializer.Serialize(largeMatrix, options);
            var deserializedLarge = JsonSerializer.Deserialize<double[,]>(json, options);
            Console.WriteLine($"\nLarge matrix [{largeMatrix.GetLength(0)}x{largeMatrix.GetLength(1)}] serialized/deserialized");
            Console.WriteLine($"Data integrity check: First={largeMatrix[0, 0]:F2}, Last={largeMatrix[99, 49]:F2}");
            Console.WriteLine($"Deserialized matches: {deserializedLarge[0, 0] == largeMatrix[0, 0]}");

            Console.WriteLine("\n=== Edge Cases Complete ===");
        }
    }
}