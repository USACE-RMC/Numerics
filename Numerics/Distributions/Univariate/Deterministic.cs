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

using Numerics.Data.Statistics;
using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// Deterministic point value estimate.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// This is used in cases, such as event tree analysis, to represent a value or event with a single probability estimate.
    /// </para>
    /// </remarks>
    [Serializable]
    public class Deterministic : UnivariateDistributionBase, IEstimation
    {
       
        /// <summary>
        /// Constructs a deterministic distribution with default value of 0.5.
        /// </summary>
        /// <remarks>
        /// This is used in cases, such as event tree analysis, to represent a value or event with a single probability estimate.
        /// </remarks>
        public Deterministic()
        {
            SetParameters(0.5d);
        }

        /// <summary>
        /// Constructs a deterministic distribution.
        /// </summary>
        /// <param name="deterministicValue">The constant value.</param>
        /// <remarks>
        /// This is used in cases, such as event tree analysis, to represent a value or event with a single probability estimate.
        /// </remarks>
        public Deterministic(double deterministicValue)
        {
            SetParameters(deterministicValue);
        }
       
        /// <summary>
        /// Gets and sets the point value estimate.
        /// </summary>
        public double Value { get; set; }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 1; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Deterministic; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Deterministic"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "D"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                parmString[0, 0] = "Value";
                parmString[0, 1] = Value.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["Value"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Value)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Value]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Value; }
        }

        /// <summary>
        /// Gets the median of the distribution. Not supported.
        /// </summary>
        public override double Median
        {
            get { return Value; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Value; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return double.NaN; }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return double.NaN; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return double.NaN; }
        }
 
        /// <summary>
        /// Gets the minimum of the distribution. Not supported.
        /// </summary>
        public override double Minimum
        {
            get { return Value; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Value; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(Statistics.Mean(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="value">The point value estimate.</param>
        public void SetParameters(double value)
        {
            Value = value;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0]);
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            // Validate probability
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Probability), "The point value must be a number.");
                return new ArgumentOutOfRangeException(nameof(Probability), "The point value must be a number.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override double[,] CreatePDFGraph()
        {
            var PDFgraph = new double[5, 2];
            PDFgraph[0, 0] = Value - 1d;
            PDFgraph[0, 1] = 0d;
            PDFgraph[1, 0] = Value;
            PDFgraph[1, 1] = 0d;
            PDFgraph[2, 0] = Value;
            PDFgraph[2, 1] = 1d;
            PDFgraph[3, 0] = Value;
            PDFgraph[3, 1] = 0d;
            PDFgraph[4, 0] = Value + 1d;
            PDFgraph[4, 1] = 0d;
            return PDFgraph;
        }

        /// <inheritdoc/>
        public override double PDF(double X)
        {
            if (X != Value) return 0.0;
            return 1d;
        }

        /// <inheritdoc/>
        public override double CDF(double X)
        {
            if (X < Value) return 0;
            return 1d;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            return Value;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Deterministic(Value);
        }
   
    }
}