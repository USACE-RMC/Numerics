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
    /// The Beta distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Beta_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class BetaDistribution : UnivariateDistributionBase
    {

        
        /// <summary>
        /// Constructs a Beta distribution with α and β = 2, defined in the interval (0,1).
        /// </summary>
        public BetaDistribution()
        {
            SetParameters([2d, 2d]);
        }

        /// <summary>
        /// Constructs a Beta distribution with the given parameters α and β, defined in the interval (0,1).
        /// </summary>
        /// <param name="alpha">The shape parameter α (alpha).</param>
        /// <param name="beta">The shape parameter β (beta).</param>
        public BetaDistribution(double alpha, double beta)
        {
            SetParameters([alpha, beta]);
        }
    
        private double _alpha;
        private double _beta;

        /// <summary>
        /// Gets and sets the shape parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, Beta }, false) is null;
                _alpha = value;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return _beta; }
            set
            {
                _parametersValid = ValidateParameters(new[] { Alpha, value }, false) is null;
                _beta = value;
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
            get { return UnivariateDistributionType.Beta; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Beta"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Beta"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Shape (α)";
                parmString[1, 0] = "Shape (β)";
                parmString[0, 1] = Alpha.ToString();
                parmString[1, 1] = Beta.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["α", "β"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Alpha), nameof(Beta)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Alpha, Beta]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Alpha / (Alpha + Beta); }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Alpha == 1d && Beta == 1d)
                    return 0.5d;
                if (Alpha <= 1d && Beta > 1d)
                    return 0.0d;
                if (Alpha > 1d && Beta <= 1d)
                    return 1.0d;
                return (Alpha - 1d) / (Alpha + Beta - 2d);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(Alpha * Beta / ((Alpha + Beta) * (Alpha + Beta) * (Alpha + Beta + 1.0d))); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 2.0d * (Beta - Alpha) * Math.Sqrt(Alpha + Beta + 1.0d) / ((Alpha + Beta + 2.0d) * Math.Sqrt(Alpha * Beta)); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double num = 6d * ((Alpha - Beta) * (Alpha - Beta) * (Alpha + Beta + 1d) - Alpha * Beta * (Alpha + Beta + 2d));
                double den = Alpha * Beta * (Alpha + Beta + 2d) * (Alpha + Beta + 3d);
                return 3.0d + num / den;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return 1.0d; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            Alpha = parameters[0];
            Beta = parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
            }
            // 
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Beta), "The shape parameter β (beta) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Beta), "The shape parameter β (beta) must be positive.");
            }
            // 
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Alpha, Beta], true);
            // 
            if (x < Minimum || x > Maximum) return 0.0d;
            double constant = 1.0d / Mathematics.SpecialFunctions.Beta.Function(Alpha, Beta);
            double a = Math.Pow(x, Alpha - 1d);
            double b = Math.Pow(1d - x, Beta - 1d);
            return constant * a * b;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Alpha, Beta], true);
            // 
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            return Mathematics.SpecialFunctions.Beta.Incomplete(Alpha, Beta, x);
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
                ValidateParameters([Alpha, Beta], true);
            // 
            return Mathematics.SpecialFunctions.Beta.IncompleteInverse(Alpha, Beta, probability);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new BetaDistribution(Alpha, Beta);
        }
       
    }
}