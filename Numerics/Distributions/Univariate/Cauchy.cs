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
    /// The Cauchy distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Cauchy_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Cauchy : UnivariateDistributionBase
    {
      
        /// <summary>
        /// Constructs a Cauchy distribution with location of 0 and scale of 1.
        /// </summary>
        public Cauchy()
        {
            SetParameters([0d, 1d]);
        }

        /// <summary>
        /// Constructs a Cauchy distribution with a given X0 and γ.
        /// </summary>
        /// <param name="location">The location parameter (X0).</param>
        /// <param name="scale">The scale parameter γ (gamma).</param>
        public Cauchy(double location, double scale)
        {
            SetParameters([location, scale]);
        }
     
        private double _gamma;

        /// <summary>
        /// Gets and sets the location parameter (X0).
        /// </summary>
        public double X0 { get; set; }

        /// <summary>
        /// Gets and sets the scale parameter γ (gamma).
        /// </summary>
        public double Gamma
        {
            get { return _gamma; }
            set
            {
                _parametersValid = ValidateParameters([X0, value], false) is null;
                _gamma = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Cauchy; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Cauchy"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "C"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Location (X0)";
                parmString[1, 0] = "Scale (γ)";
                parmString[0, 1] = X0.ToString();
                parmString[1, 1] = Gamma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["X0", "γ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(X0), nameof(Gamma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [X0, Gamma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return double.NaN; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return X0; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return X0; }
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

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            X0 = parameters[0];
            Gamma = parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(X0), "The location parameter X0 must be valid.");
                return new ArgumentOutOfRangeException(nameof(X0), "The location parameter X0 must be valid.");
            }

            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "The scale parameter γ (gamma) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "The scale parameter γ (gamma) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([X0, Gamma], true);
            double z = (x - X0) / Gamma;
            return 1.0d / (Math.PI * Gamma * (1.0d + z * z));
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([X0, Gamma], true);
            return 1.0d / Math.PI * Math.Atan2(x - X0, Gamma) + 0.5d;
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
                ValidateParameters([X0, Gamma], true);
            return X0 + Gamma * Math.Tan(Math.PI * (probability - 0.5d));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Cauchy(X0, Gamma);
        }
    
    }
}