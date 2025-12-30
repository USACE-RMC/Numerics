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
    /// The four-parameter Beta distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Beta_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class GeneralizedBeta : UnivariateDistributionBase
    {
     
        /// <summary>
        /// Constructs a Beta distribution with α and β = 2, defined in the interval (0,1).
        /// </summary>
        public GeneralizedBeta()
        {
            SetParameters(2d, 2d, 0d, 1d);
        }

        /// <summary>
        /// Constructs a Beta distribution with the given parameters α and β, defined in the interval (0,1).
        /// </summary>
        /// <param name="alpha">The shape parameter α (alpha).</param>
        /// <param name="beta">The shape parameter β (beta).</param>
        public GeneralizedBeta(double alpha, double beta)
        {
            SetParameters(alpha, beta, 0d, 1d);
        }

        /// <summary>
        /// Constructs a Beta distribution with the given parameters α, β, a and b, defined in the interval (a,b).
        /// </summary>
        /// <param name="alpha">The shape parameter α (alpha).</param>
        /// <param name="beta">The shape parameter β (beta).</param>
        /// <param name="min">The minimum possible value.</param>
        /// <param name="max">The maximum possible value.</param>
        public GeneralizedBeta(double alpha, double beta, double min, double max)
        {
            SetParameters(alpha, beta, min, max);
        }

        /// <summary>
        /// Constructs a Beta distribution using the PERT estimation method with the parameters a, b, and mode, defined in the interval (a,b).
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <returns>
        /// A Beta distribution parameterized using the PERT method.
        /// </returns>
        /// <remarks>
        /// <see href = "https://en.wikipedia.org/wiki/PERT_distribution" />
        /// </remarks>
        public static GeneralizedBeta PERT(double min, double mode, double max)
        {
            return PERT(min, mode, max, 4d);
        }

        /// <summary>
        /// Constructs a Beta distribution using the PERT estimation method with the parameters a, b, mode and λ, defined in the interval (a,b).
        /// </summary>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="mode">The mode, or most likely, value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <param name="scale">The scale parameter λ (lambda). Default is 4.</param>
        /// <returns>
        /// A Beta distribution parameterized using the PERT method.
        /// </returns>
        /// <remarks>
        /// <see href = "https://en.wikipedia.org/wiki/PERT_distribution" />
        /// </remarks>
        public static GeneralizedBeta PERT(double min, double mode, double max, double scale)
        {
            // validate parameters
            if (min > max)
            {
                throw new ArgumentOutOfRangeException(nameof(max), "The maximum value must be greater than the minimum value.");
            }
            else if (mode < min || mode > max)
            {
                throw new ArgumentOutOfRangeException(nameof(mode), "The mode must be between the minimum and maximum values.");
            }
            // get alpha and beta
            double mean = (min + scale * mode + max) / (scale + 2d);
            double alpha = 1d + scale / 2d;
            if (!mean.AlmostEquals(mode))
            {
                alpha = (mean - min) * (2d * mode - min - max) / ((mode - mean) * (max - min));
            }
            double beta = alpha * (max - mean) / (mean - min);
            // return beta
            return new GeneralizedBeta(alpha, beta, min, max);
        }
      
        private double _alpha;
        private double _beta;
        private double _min = double.MinValue;
        private double _max = double.MaxValue;
    
        /// <summary>
        /// Gets and sets the shape parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters(value, Beta, Min, Max, false) is null;
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
                _parametersValid = ValidateParameters(Alpha, value, Min, Max, false) is null;
                _beta = value;
            }
        }

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set
            {
                _parametersValid = ValidateParameters(Alpha, Beta, value, Max, false) is null;
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
                _parametersValid = ValidateParameters(Alpha, Beta, Min, value, false) is null;
                _max = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 4; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.GeneralizedBeta; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Beta"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GB"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[4, 2];
                parmString[0, 0] = "Shape (α)";
                parmString[1, 0] = "Shape (β)";
                parmString[2, 0] = "Min";
                parmString[3, 0] = "Max";
                parmString[0, 1] = Alpha.ToString();
                parmString[1, 1] = Beta.ToString();
                parmString[2, 1] = Min.ToString();
                parmString[3, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["α", "β", "Min", "Max"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Alpha), nameof(Beta), nameof(Min), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Alpha, Beta, Min, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                double _mean = Alpha / (Alpha + Beta);
                return _mean * (Max - Min) + Min;
            }
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
                double _mode = (Alpha - 1.0d) / (Alpha + Beta - 2.0d);
                return _mode * (Max - Min) + Min;
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                double var = Alpha * Beta / ((Alpha + Beta) * (Alpha + Beta) * (Alpha + Beta + 1d));
                return Math.Sqrt(var) * (Max - Min);
            }
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
                double num = 6d * ((Alpha - Beta) * (Alpha - Beta) * (Alpha + Beta + 1.0d) - Alpha * Beta * (Alpha + Beta + 2.0d));
                double den = Alpha * Beta * (Alpha + Beta + 2.0d) * (Alpha + Beta + 3.0d);
                return 3d + num / den;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return _min; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return _max; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d, double.NegativeInfinity, Min]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, Max, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name ="alpha"> The shape parameter α (alpha).</param>
        /// <param name="beta">The shape parameter β (beta).</param>
        /// <param name="min">The minimum possible value.</param>
        /// <param name="max">The maximum possible value.</param>
        public void SetParameters(double alpha, double beta, double min, double max)
        {
            // validate parameters
            _parametersValid = ValidateParameters(alpha, beta, min, max, false) is null;
            // Set parameters
            _alpha = alpha;
            _beta = beta;
            _min = min;
            _max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2], parameters[3]);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mu"></param>
        /// <param name="sigma"></param>
        /// <param name="min"></param>
        /// <param name="max"></param>
        public void SetParametersFromMoments(double mu, double sigma, double min, double max)
        {
            var s2 = sigma * sigma;
            _alpha = (min - mu) * (min * max - min * mu - max * mu + mu * mu + s2) / (s2 * (max - min));
            _beta = -(max - mu) * (min * max - min * mu - max * mu + mu * mu + s2) / (s2 * (max - min));
            _min = min;
            _max = max;
            // validate parameters
            _parametersValid = ValidateParameters(_alpha, _beta, _min, _max, false) is null;
        }


        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name ="alpha"> The shape parameter α (alpha).</param>
        /// <param name="beta">The shape parameter β (beta).</param>
        /// <param name="min">The minimum possible value.</param>
        /// <param name="max">The maximum possible value.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double alpha, double beta, double min, double max, bool throwException)
        {
            if (double.IsNaN(alpha) || double.IsInfinity(alpha) || alpha <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The shape parameter α (alpha) must be positive.");
            }
            if (double.IsNaN(beta) || double.IsInfinity(beta) || beta <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Beta), "The shape parameter β (beta) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Beta), "The shape parameter β (beta) must be positive.");
            }
            if (double.IsNaN(min) || double.IsInfinity(min) ||
                double.IsNaN(max) || double.IsInfinity(max) || min >= max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than or equal to the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than or equal to the max.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], parameters[3], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Alpha, Beta, Min, Max, true);
            if (X < Minimum || X > Maximum) return 0.0d;

            if (Min == Max) return 0d;
            double constant = 1.0d / Mathematics.SpecialFunctions.Beta.Function(Alpha, Beta);
            double z = (X - Min) / (Max - Min);
            double a = Math.Pow(z, Alpha - 1d);
            double b = Math.Pow(1d - z, Beta - 1d);
            return constant * a * b / (Max - Min);
        }

        /// <inheritdoc/>
        public override double CDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Alpha, Beta, Min, Max, true);

            if (Min == Max) return 1d;
            if (X <= Minimum) return 0d;
            if (X >= Maximum) return 1d;        
            double z = (X - Min) / (Max - Min);
            return Mathematics.SpecialFunctions.Beta.Incomplete(Alpha, Beta, z);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Alpha, Beta, Min, Max, true);
            double z = Mathematics.SpecialFunctions.Beta.IncompleteInverse(Alpha, Beta, probability);
            return z * (Max - Min) + Min;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedBeta(Alpha, Beta, Min, Max);
        }
    
    }
}