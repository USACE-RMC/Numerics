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
using System.Text.Json;
using System.Text.Json.Serialization;
using Numerics.Distributions;

namespace Numerics.Utilities
{
    /// <summary>
    /// Custom JSON converter for 2D double arrays.
    /// Serializes 2D arrays as an object with dimensions and flattened data.
    /// </summary>
    /// <remarks>
    /// <para>
    /// This converter serializes a double[,] array into a JSON object with three properties:
    /// "rows", "cols", and "data" (a flattened 1D array). This format is more efficient than
    /// serializing nested arrays and preserves the exact dimensions.
    /// </para>
    /// <para>
    /// <b>JSON Format:</b>
    /// <code>
    /// {
    ///   "rows": 2,
    ///   "cols": 3,
    ///   "data": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    /// }
    /// </code>
    /// </para>
    /// </remarks>
    public class Double2DArrayConverter : JsonConverter<double[,]>
    {

        /// <summary>
        /// Reads and converts JSON to a 2D double array.
        /// </summary>
        /// <param name="reader">The JSON reader.</param>
        /// <param name="typeToConvert">The type to convert.</param>
        /// <param name="options">Serialization options.</param>
        /// <returns>
        /// A 2D double array reconstructed from the JSON, or null if the JSON is null.
        /// Returns an empty array (0x0) if the data is invalid.
        /// </returns>
        /// <exception cref="JsonException">
        /// Thrown when the JSON format is invalid (e.g., missing StartObject token).
        /// </exception>
        /// <remarks>
        /// Expects JSON in the format: { "rows": int, "cols": int, "data": double[] }
        /// </remarks>
        public override double[,] Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null!;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            int rows = 0;
            int cols = 0;
            double[]? data = null!;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string? propertyName = reader.GetString();
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

        /// <summary>
        /// Writes a 2D double array as JSON.
        /// </summary>
        /// <param name="writer">The JSON writer.</param>
        /// <param name="value">The 2D double array to serialize.</param>
        /// <param name="options">Serialization options.</param>
        /// <remarks>
        /// <para>
        /// Serializes the array as: { "rows": int, "cols": int, "data": double[] }
        /// where data contains row-major flattened values.
        /// </para>
        /// <para>
        /// Null arrays are serialized as JSON null.
        /// </para>
        /// </remarks>
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
    /// <remarks>
    /// <para>
    /// This converter serializes a string[,] array into a JSON object with three properties:
    /// "rows", "cols", and "data" (a flattened 1D array). This format is more efficient than
    /// serializing nested arrays and preserves the exact dimensions.
    /// </para>
    /// <para>
    /// <b>JSON Format:</b>
    /// <code>
    /// {
    ///   "rows": 2,
    ///   "cols": 2,
    ///   "data": ["A", "B", "C", "D"]
    /// }
    /// </code>
    /// </para>
    /// </remarks>
    public class String2DArrayConverter : JsonConverter<string[,]>
    {

        /// <summary>
        /// Reads and converts JSON to a 2D string array.
        /// </summary>
        /// <param name="reader">The JSON reader.</param>
        /// <param name="typeToConvert">The type to convert.</param>
        /// <param name="options">Serialization options.</param>
        /// <returns>
        /// A 2D string array reconstructed from the JSON, or null if the JSON is null.
        /// Returns an empty array (0x0) if the data is invalid.
        /// </returns>
        /// <exception cref="JsonException">
        /// Thrown when the JSON format is invalid (e.g., missing StartObject token).
        /// </exception>
        /// <remarks>
        /// Expects JSON in the format: { "rows": int, "cols": int, "data": string[] }
        /// </remarks>
        public override string[,] Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null!;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            int rows = 0;
            int cols = 0;
            string[]? data = null!;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string? propertyName = reader.GetString();
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

        /// <summary>
        /// Writes a 2D string array as JSON.
        /// </summary>
        /// <param name="writer">The JSON writer.</param>
        /// <param name="value">The 2D string array to serialize.</param>
        /// <param name="options">Serialization options.</param>
        /// <remarks>
        /// <para>
        /// Serializes the array as: { "rows": int, "cols": int, "data": string[] }
        /// where data contains row-major flattened values.
        /// </para>
        /// <para>
        /// Null arrays are serialized as JSON null.
        /// </para>
        /// </remarks>
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
    /// <remarks>
    /// <para>
    /// This converter enables serialization of distribution objects by storing
    /// their type and parameters. Upon deserialization, it uses the 
    /// UnivariateDistributionFactory to recreate the distribution.
    /// </para>
    /// <para>
    /// <b>JSON Format:</b>
    /// <code>
    /// {
    ///   "Type": "Normal",
    ///   "Parameters": [0.0, 1.0]
    /// }
    /// </code>
    /// </para>
    /// <para>
    /// <b>Supported Distributions:</b> All distributions in UnivariateDistributionType enum.
    /// </para>
    /// <para>
    /// <b>Limitations:</b> Only the distribution type and parameters are preserved.
    /// Other properties (like Name, Description) are not serialized.
    /// </para>
    /// </remarks>
    public class UnivariateDistributionConverter : JsonConverter<UnivariateDistributionBase>
    {
        /// <summary>
        /// Reads and converts JSON to a UnivariateDistributionBase instance.
        /// </summary>
        /// <param name="reader">The JSON reader.</param>
        /// <param name="typeToConvert">The type to convert.</param>
        /// <param name="options">Serialization options.</param>
        /// <returns>
        /// A UnivariateDistributionBase instance reconstructed from the JSON, or null if:
        /// - The JSON is null
        /// - The distribution type is missing or invalid
        /// - The parameters are missing or invalid
        /// - The distribution cannot be created
        /// </returns>
        /// <exception cref="JsonException">
        /// Thrown when the JSON format is invalid (e.g., missing StartObject token).
        /// </exception>
        /// <remarks>
        /// <para>
        /// Uses UnivariateDistributionFactory to create a default distribution of the
        /// specified type, then applies the parameters from JSON.
        /// </para>
        /// <para>
        /// If the distribution cannot be reconstructed (e.g., invalid parameters),
        /// returns null rather than throwing an exception.
        /// </para>
        /// </remarks>
        public override UnivariateDistributionBase Read(ref Utf8JsonReader reader, Type typeToConvert, JsonSerializerOptions options)
        {
            if (reader.TokenType == JsonTokenType.Null)
                return null!;

            if (reader.TokenType != JsonTokenType.StartObject)
                throw new JsonException("Expected StartObject token");

            UnivariateDistributionType? distributionType = null;
            double[]? parameters = null!;

            while (reader.Read())
            {
                if (reader.TokenType == JsonTokenType.EndObject)
                    break;

                if (reader.TokenType == JsonTokenType.PropertyName)
                {
                    string? propertyName = reader.GetString();
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
                return null!;

            // Use the factory to create a default distribution, then set its parameters
            try
            {
                var distribution = UnivariateDistributionFactory.CreateDistribution(distributionType.Value);
                if (distribution != null! && parameters != null && parameters.Length > 0)
                {
                    distribution.SetParameters(parameters);
                }
                return distribution!;
            }
            catch
            {
                // If we can't recreate it, return null
                return null!;
            }
        }

        /// <summary>
        /// Writes a UnivariateDistributionBase instance as JSON.
        /// </summary>
        /// <param name="writer">The JSON writer.</param>
        /// <param name="value">The distribution to serialize.</param>
        /// <param name="options">Serialization options.</param>
        /// <remarks>
        /// <para>
        /// Serializes the distribution as: { "Type": string, "Parameters": double[] }
        /// </para>
        /// <para>
        /// The Type property contains the distribution type (e.g., "Normal", "Exponential").
        /// The Parameters property contains the distribution's parameter values in the order
        /// defined by the distribution's GetParameters property.
        /// </para>
        /// <para>
        /// If the distribution is null, writes JSON null.
        /// If parameters cannot be retrieved, writes an empty array.
        /// </para>
        /// </remarks>
        public override void Write(Utf8JsonWriter writer, UnivariateDistributionBase value, JsonSerializerOptions options)
        {
            if (value == null!)
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