using System;
using System.Text.Json;
using System.Text.Json.Serialization;
using Numerics.Distributions;

namespace Numerics.Utilities
{
    /// <summary>
    /// Custom JSON converter for 2D double arrays.
    /// Serializes 2D arrays as an object with dimensions and flattened data.
    /// </summary>
    public class Double2DArrayConverter : JsonConverter<double[,]>
    {
        public override double[,] Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            int rows = 0;
            int cols = 0;
            double[] data = null;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string propertyName = reader.GetString();
                    reader.Read();

                    switch (propertyName)
                    {
                        case "rows":
                            rows = reader.GetInt32();
                            break;
                        case "cols":
                            cols = reader.GetInt32();
                            break;
                        case "data":
                            data = JsonSerializer.Deserialize<double[]>(ref reader, options);
                            break;
                    }
                }
            }

            if (data == null || rows == 0 || cols == 0)
                return new double[0, 0];

            double[,] result = new double[rows, cols];
            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (index < data.Length)
                        result[i, j] = data[index++];
                }
            }

            return result;
        }

        public override void Write(Utf8JsonWriter writer, double[,] value, JsonSerializerOptions options)
        {
            if (value == null)
            {
                writer.WriteNullValue();
                return;
            }

            int rows = value.GetLength(0);
            int cols = value.GetLength(1);
            double[] data = new double[rows * cols];

            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    data[index++] = value[i, j];
                }
            }

            writer.WriteStartObject();
            writer.WriteNumber("rows", rows);
            writer.WriteNumber("cols", cols);
            writer.WritePropertyName("data");
            JsonSerializer.Serialize(writer, data, options);
            writer.WriteEndObject();
        }
    }

    /// <summary>
    /// Custom JSON converter for 2D string arrays.
    /// Serializes 2D arrays as an object with dimensions and flattened data.
    /// </summary>
    public class String2DArrayConverter : JsonConverter<string[,]>
    {
        public override string[,] Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            int rows = 0;
            int cols = 0;
            string[] data = null;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string propertyName = reader.GetString();
                    reader.Read();

                    switch (propertyName)
                    {
                        case "rows":
                            rows = reader.GetInt32();
                            break;
                        case "cols":
                            cols = reader.GetInt32();
                            break;
                        case "data":
                            data = JsonSerializer.Deserialize<string[]>(ref reader, options);
                            break;
                    }
                }
            }

            if (data == null || rows == 0 || cols == 0)
                return new string[0, 0];

            string[,] result = new string[rows, cols];
            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    if (index < data.Length)
                        result[i, j] = data[index++];
                }
            }

            return result;
        }

        public override void Write(Utf8JsonWriter writer, string[,] value, JsonSerializerOptions options)
        {
            if (value == null)
            {
                writer.WriteNullValue();
                return;
            }

            int rows = value.GetLength(0);
            int cols = value.GetLength(1);
            string[] data = new string[rows * cols];

            int index = 0;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    data[index++] = value[i, j];
                }
            }

            writer.WriteStartObject();
            writer.WriteNumber("rows", rows);
            writer.WriteNumber("cols", cols);
            writer.WritePropertyName("data");
            JsonSerializer.Serialize(writer, data, options);
            writer.WriteEndObject();
        }
    }

    /// <summary>
    /// Custom JSON converter for UnivariateDistributionBase.
    /// Serializes only essential properties needed for reconstruction.
    /// </summary>
    public class UnivariateDistributionConverter : JsonConverter<UnivariateDistributionBase>
    {
        public override UnivariateDistributionBase Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            UnivariateDistributionType? distributionType = null;
            double[] parameters = null;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string propertyName = reader.GetString();
                    reader.Read();

                    switch (propertyName)
                    {
                        case "Type":
                            distributionType = JsonSerializer.Deserialize<UnivariateDistributionType>(ref reader, options);
                            break;
                        case "Parameters":
                            parameters = JsonSerializer.Deserialize<double[]>(ref reader, options);
                            break;
                    }
                }
            }

            if (!distributionType.HasValue || parameters == null)
                return null;

            // Use the factory to create a default distribution, then set its parameters
            try
            {
                var distribution = UnivariateDistributionFactory.CreateDistribution(distributionType.Value);
                if (distribution != null && parameters != null && parameters.Length > 0)
                {
                    distribution.SetParameters(parameters);
                }
                return distribution;
            }
            catch
            {
                // If we can't recreate it, return null
                return null;
            }
        }

        public override void Write(Utf8JsonWriter writer, UnivariateDistributionBase value, JsonSerializerOptions options)
        {
            if (value == null)
            {
                writer.WriteNullValue();
                return;
            }

            writer.WriteStartObject();

            // Write the distribution type
            writer.WritePropertyName("Type");
            JsonSerializer.Serialize(writer, value.Type, options);

            // Write the parameters
            writer.WritePropertyName("Parameters");
            try
            {
                var parameters = value.GetParameters;
                JsonSerializer.Serialize(writer, parameters, options);
            }
            catch
            {
                // If we can't get parameters, write an empty array
                JsonSerializer.Serialize(writer, new double[0], options);
            }

            writer.WriteEndObject();
        }
    }
}