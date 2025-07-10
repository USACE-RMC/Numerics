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
    /// The Truncated Normal (Gaussian) probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Wikipedia contributors, "Truncated Normal distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Truncated_normal_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class TruncatedNormal : UnivariateDistributionBase, IEstimation
    {
   
        /// <summary>
        /// Constructs a Truncated Normal (Gaussian) distribution with a mean of 0.5 and standard deviation of 0.1.
        /// Min is set to 0.0 and max is set to 1.0.
        /// </summary>
        public TruncatedNormal()
        {
            SetParameters(0.5d, 0.2d, 0.0d, 1.0d);
        }

        /// <summary>
        /// Constructs a Normal (Gaussian) distribution with given mean, standard deviation, min and max.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public TruncatedNormal(double mean, double standardDeviation, double min, double max)
        {
            SetParameters(mean, standardDeviation, min, max);
        }

        private double _mu;
        private double _sigma;
        private double _min;
        private double _max;
        private double Xi;
        private double alpha;
        private double beta;
        private double Z;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, Min, Max, false) is null;
                _mu = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter σ (sigma).
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, Min, Max, false) is null;
                _sigma = value;
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
                _parametersValid = ValidateParameters(Mu, Sigma, value, Max, false) is null;
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
                _parametersValid = ValidateParameters(Mu, Sigma, Min, value, false) is null;
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
            get { return UnivariateDistributionType.TruncatedNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Truncated Normal"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Trunc. N"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[4, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[2, 0] = "Min";
                parmString[3, 0] = "Max";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                parmString[2, 1] = Min.ToString();
                parmString[3, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ", "Min", "Max"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma), nameof(Min), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma, Min, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                alpha = (Min - Mu) / Sigma;
                beta = (Max - Mu) / Sigma;
                Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
                if (Z == 0) return Mu;
                return Mu + Sigma * ((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / Z);
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get
            {
                alpha = (Min - Mu) / Sigma;
                beta = (Max - Mu) / Sigma;
                Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
                if (Z == 0) return Mu;
                return Mu + Sigma * Normal.StandardZ((Normal.StandardCDF(alpha) + Normal.StandardCDF(beta)) / 2d);
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Mu < Min)
                {
                    return Min;
                }
                else if (Mu >= Min && Mu <= Max)
                {
                    return Mu;
                }
                else if (Mu > Max)
                {
                    return Max;
                }
                return default;
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                alpha = (Min - Mu) / Sigma;
                beta = (Max - Mu) / Sigma;
                Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
                if (Z == 0) return Sigma;
                return Math.Sqrt(Math.Pow(Sigma, 2d) * (1d + (alpha * Normal.StandardPDF(alpha) - beta * Normal.StandardPDF(beta)) / Z - Math.Pow((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / Z, 2d)));
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Reference: Pearson diagrams for truncated normal and truncated Weibull distributions. Sugiura, Gomi. 1985.
        /// </remarks>
        public override double Skewness
        {
            get
            {
                double K0 = (Min - Mu) / Sigma;
                double K1 = (Max - Mu) / Sigma;
                Z = Normal.StandardCDF(K1) - Normal.StandardCDF(K0);
                if (Z == 0) return 0;
                double Z0 = Normal.StandardPDF(K0) / Z;
                double Z1 = Normal.StandardPDF(K1) / Z;
                double V = 1d - (K1 * Z1 - K0 * Z0) - Math.Pow(Z1 - Z0, 2d);
                return -(1d / Math.Pow(V, 3d / 2d)) * (2d * Math.Pow(Z1 - Z0, 3d) + (3d * K1 * Z1 - 3d * K0 * Z0 - 1d) * (Z1 - Z0) + Math.Pow(K1, 2d) * Z1 - Math.Pow(K0, 2d) * Z0);
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double K0 = (Min - Mu) / Sigma;
                double K1 = (Max - Mu) / Sigma;
                Z = Normal.StandardCDF(K1) - Normal.StandardCDF(K0);
                if (Z == 0) return 3d;
                double Z0 = Normal.StandardPDF(K0) / Z;
                double Z1 = Normal.StandardPDF(K1) / Z;
                double V = 1d - (K1 * Z1 - K0 * Z0) - Math.Pow(Z1 - Z0, 2d);
                return (1d / Math.Pow(V, 2d)) * (-3 * Math.Pow(Z1 - Z0, 4d) - 6d * (K1 * Z1 - K0 * Z0) * Math.Pow(Z1 - Z0, 2d) - 2d * Math.Pow(Z1 - Z0, 2d) - 4d * (Math.Pow(K1, 2d) * Z1 - Math.Pow(K0, 2d) * Z0) * (Z1 - Z0) - 3d * (K1 * Z1 - K0 * Z0) - (Math.Pow(K1, 3d) * Z1 - Math.Pow(K0, 3d) * Z0) + 3d);
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
            get { return [double.NegativeInfinity, 0.0d, double.NegativeInfinity, Min]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, Max, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                var mom = Statistics.ProductMoments(sample);
                var min = Statistics.Minimum(sample);
                var max = Statistics.Maximum(sample);
                SetParameters(mom[0], mom[1], min, max);
            }
            else
            {
                throw new NotImplementedException();
            }
        }


        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public void SetParameters(double mean, double standardDeviation, double min, double max)
        {
            // Set parameters
            Mu = mean;
            Sigma = standardDeviation;
            Min = min;
            Max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2], parameters[3]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double mean, double standardDeviation, double min, double max, bool throwException)
        {
            if (double.IsNaN(mean) || double.IsInfinity(Mean))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mean must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mean must be a number.");
            }
            if (double.IsNaN(standardDeviation) || double.IsInfinity(standardDeviation) || standardDeviation <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
            }
            if (double.IsNaN(min) || double.IsNaN(max) || double.IsInfinity(min) || double.IsInfinity(max) || min > max)
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
            return ValidateParameters(parameters[0], parameters[1], parameters[2], parameters[3], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Min, Max, true);
            // 
            if (x < Minimum || x > Maximum) return 0.0d;
            Xi = (x - Mu) / Sigma;
            alpha = (Min - Mu) / Sigma;
            beta = (Max - Mu) / Sigma;
            Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
            if (Z == 0) return Normal.StandardPDF(Xi);
            return Normal.StandardPDF(Xi) / (Sigma * Z);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Min, Max, true);
            
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return 1d;
            Xi = (x - Mu) / Sigma;
            alpha = (Min - Mu) / Sigma;
            beta = (Max - Mu) / Sigma;
            Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
            if (Z == 0) return Normal.StandardCDF(Xi);
            return (Normal.StandardCDF(Xi) - Normal.StandardCDF(alpha)) / Z;
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
                ValidateParameters(Mu, Sigma, Min, Max, true);
            alpha = (Min - Mu) / Sigma;
            beta = (Max - Mu) / Sigma;
            Z = Normal.StandardCDF(beta) - Normal.StandardCDF(alpha);
            if (Z == 0) return Mu + Sigma * Normal.StandardZ(probability);
            return Mu + Sigma * Normal.StandardZ(Normal.StandardCDF(alpha) + probability * Z);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new TruncatedNormal(Mu, Sigma, Min, Max);
        }
          
    }
}