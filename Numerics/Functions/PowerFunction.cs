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

using Numerics.Distributions;
using System;
using System.Collections.Generic;

namespace Numerics.Functions
{
    /// <summary>
    /// A class for a power function with normally distributed noise.
    /// Y = [α * (X - ξ)^β] * ϵ, where ϵ ~ N(0,σ) 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class PowerFunction : IUnivariateFunction
    {

        /// <summary>
        /// Construct a new deterministic power function with α=1 and β=1.5 and ξ=0. 
        /// </summary>
        public PowerFunction()
        {
            Alpha = 1;
            Beta = 1.5;
            Xi = 0;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new power function with a given α, β, and ξ. 
        /// </summary>
        /// <param name="alpha">The coefficient parameter α.</param>
        /// <param name="beta">The exponent parameter β.</param>
        /// <param name="xi">The location parameter ξ.</param>
        public PowerFunction(double alpha, double beta, double xi = 0)
        {
            Alpha = alpha;
            Beta = beta;
            Xi = xi;
            IsDeterministic = true;
        }

        /// <summary>
        /// Construct a new power function with a given α, β, ξ and standard error σ.
        /// </summary>
        /// <param name="alpha">The coefficient parameter α.</param>
        /// <param name="beta">The exponent parameter β.</param>
        /// <param name="xi">The location parameter ξ.</param>
        /// <param name="sigma">The log-space standard error σ.</param>
        public PowerFunction(double alpha, double beta, double xi, double sigma)
        {
            Alpha = alpha;
            Beta = beta;
            Sigma = sigma;
            Xi = xi;
            IsDeterministic = false;
        }

        private bool _parametersValid = true;
        private double _alpha, _beta, _xi, _sigma;
        private Normal _normal = new Normal();

        /// <summary>
        /// The coefficient α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, Beta, Xi, Sigma }, false) is null;
                _alpha = value;
            }
        }

        /// <summary>
        /// The exponent β (beta).
        /// </summary>
        public double Beta
        {
            get { return _beta; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, value, Xi, Sigma }, false) is null;
                _beta = value;
            }
        }

        /// <summary>
        /// The location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, Beta, value, Sigma }, false) is null;
                _xi = value;
            }
        }

        /// <summary>
        /// The standard error parameter σ (sigma) in log-space.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, Beta, Xi, value }, false) is null;
                _sigma = value;
                _normal.SetParameters(0, _sigma);
            }
        }

        /// <inheritdoc/>
        public int NumberOfParameters => 4;

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public double Minimum
        {
            get { return Xi; }
            set { return; }
        }

        /// <inheritdoc/>
        public double Maximum { get; set; } = double.MaxValue;

        /// <inheritdoc/>
        public double[] MinimumOfParameters => new double[] { 0, -10, 0 };

        /// <inheritdoc/>
        public double[] MaximumOfParameters => new double[] { double.MaxValue, 10, double.MaxValue };

        /// <inheritdoc/>
        public bool IsDeterministic { get; set; }

        /// <summary>
        /// Determines if the power function should be inverted.
        /// </summary>
        public bool IsInverse { get; set; }

        /// <inheritdoc/>
        public double ConfidenceLevel { get; set; } = -1;

        /// <inheritdoc/>
        public void SetParameters(IList<double> parameters)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(parameters, false) is null;
            // Set parameters
            _alpha = parameters[0];
            _beta = parameters[1];
            _xi = parameters[2];
            _sigma = parameters[3];
        }

        /// <inheritdoc/>
        public ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (IsDeterministic == false && parameters[3] <= 0)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard error must be greater than zero.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard error must be greater than zero.");
            }
            return null;
        }

        /// <inheritdoc/>
        public double Function(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Alpha, Beta,  Xi, Sigma }, true);

            // Check support
            if (x <= Xi)
                x = Xi + Tools.DoubleMachineEpsilon;
            else if (x >= Maximum)
                x = Maximum;

            double y = 0;
            if (IsInverse)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    y = Math.Exp((Math.Log(x) - Math.Log(Alpha)) / Beta) + Xi;
                }
                else
                {
                    y = Math.Exp((Math.Log(x) - Math.Log(Alpha) - _normal.InverseCDF(ConfidenceLevel)) / Beta) + Xi;
                }
            }
            else
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    y = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(x - Xi));
                }
                else
                {
                    y = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(x - Xi) + _normal.InverseCDF(ConfidenceLevel));
                }
            }
            return y;
        }

        /// <inheritdoc/>
        public double InverseFunction(double y)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(new[] { Alpha, Beta, Xi, Sigma }, true);

            double x = 0;
            if (IsInverse)
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    x = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(y - Xi));
                }
                else
                {
                    x = Math.Exp(Math.Log(Alpha) + Beta * Math.Log(y - Xi) + _normal.InverseCDF(ConfidenceLevel));
                }
            }
            else
            {
                if (IsDeterministic == true || ConfidenceLevel < 0 || ConfidenceLevel > 1)
                {
                    x = Math.Exp((Math.Log(y) - Math.Log(Alpha)) / Beta) + Xi;
                }
                else
                {
                    x = Math.Exp((Math.Log(y) - Math.Log(Alpha) - _normal.InverseCDF(ConfidenceLevel)) / Beta) + Xi;
                }
            }
            if (x < Minimum) return Minimum;
            if (x > Maximum) return Maximum;
            return x;
        }

    }
}
