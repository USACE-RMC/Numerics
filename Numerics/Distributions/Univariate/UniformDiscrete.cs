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

namespace Numerics.Distributions
{

    /// <summary>
    /// The uniform (discrete) distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Discrete_uniform_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class UniformDiscrete : UnivariateDistributionBase
    {
   
        /// <summary>
        /// Constructs a Uniform (discrete) distribution with min = 0 and max = 1.
        /// </summary>
        public UniformDiscrete()
        {
            SetParameters(0d, 1d);
        }

        /// <summary>
        /// Constructs a Uniform (discrete) distribution with specified min and max.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        public UniformDiscrete(double min, double max)
        {
            SetParameters(min, max);
        }

        private double _min;
        private double _max;

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set
            {
                _parametersValid = ValidateParameters(value, Max, false) is null;
                _min = value;
            }
        }

        /// <summary>
        /// Get and set the max of the distribution.
        /// </summary>
        public double Max
        {
            get { return _max; }
            set
            {
                _parametersValid = ValidateParameters(Min, value, false) is null;
                _max = value;
            }
        }

        /// <summary>
        /// Returns the sample size.
        /// </summary>
        public double N
        {
            get { return Max - Min + 1d; }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.UniformDiscrete; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Uniform (Discrete)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "UD"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Min";
                parmString[1, 0] = "Max";
                parmString[0, 1] = Min.ToString();
                parmString[1, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["Min", "Max"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Min), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Min, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return (Min + Max) / 2.0d; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return (Min + Max) / 2.0d; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return double.NaN; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt((Max - Min) * (Max - Min) / 12.0d); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3.0d - 6d * (N * N + 1d) / (5d * (N * N - 1d)); }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Min; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Max; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, Min]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Max, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        public void SetParameters(double min, double max)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(min, max, false) is null;
            // Set parameters
            _min = min;
            _max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double min, double max, bool throwException)
        {
            if (double.IsNaN(min) || double.IsInfinity(min) ||
                double.IsNaN(max) || double.IsInfinity(max) || min > max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (k < Minimum || k > Maximum) return 0.0d;
            return 1d / N;
        }

        /// <inheritdoc/>
        public override double CDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (k < Minimum)
                return 0d;
            if (k > Maximum)
                return 1d;
            return (k - Min + 1.0d) / N;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            return Math.Floor(Min + probability * N);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new UniformDiscrete(Min, Max);
        }
   
    }
}