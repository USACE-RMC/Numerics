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
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Inverse Chi-Squared (Inv-χ²) probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Inverse-chi-squared_distribution" />
    /// <see href = "https://en.wikipedia.org/wiki/Scaled_inverse_chi-squared_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class InverseChiSquared : UnivariateDistributionBase
    {
   
        /// <summary>
        /// Constructs an Inverse Chi-Squared distribution with 10 degrees of freedom.
        /// </summary>
        public InverseChiSquared()
        {
            SetParameters([10d, 1d]);
        }

        /// <summary>
        /// Constructs an Inverse Chi-Squared distribution with given degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom for the distribution.</param>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        public InverseChiSquared(int degreesOfFreedom, double scale)
        {
            SetParameters([degreesOfFreedom, scale]);
        }

        private int _degreesOfFreedom;
        private double _sigma;

        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public int DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters([value, Sigma], false) is null;
                _degreesOfFreedom = value;
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
                _parametersValid = ValidateParameters([DegreesOfFreedom, value], false) is null;
                _sigma = value;
            }
        }

        /// <summary>
        /// Gets the precision, or inverse variance 1/σ^2 = τ^2
        /// </summary>
        public double Precision
        {
            get { return 1d / (Sigma * Sigma); }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.InverseChiSquared; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Inverse Chi-Squared (χ²)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Inv-χ²"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[1, 0] = "Scale (σ)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                parmString[1, 1] = Sigma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν", "σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames => [nameof(DegreesOfFreedom), nameof(Sigma)];

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [DegreesOfFreedom, Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (DegreesOfFreedom <= 2)
                    return double.NaN;
                return DegreesOfFreedom  * Sigma / (DegreesOfFreedom - 2);
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
            get { return DegreesOfFreedom * Sigma / (DegreesOfFreedom + 2); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (DegreesOfFreedom <= 4)
                    return double.NaN;
                int v = DegreesOfFreedom;
                double t2 = Sigma;
                return Math.Sqrt(2.0d * v * v * t2 * t2 / ((v - 2) * (v - 2) * (v - 4)));
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (DegreesOfFreedom <= 6)
                    return double.NaN;
                int v = DegreesOfFreedom;
                return 4.0d / (v - 6) * Math.Sqrt(2.0d * (v - 4));
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (DegreesOfFreedom <= 8)
                    return double.NaN;
                int v = DegreesOfFreedom;
                return 3.0d + 12.0d * (5.0d * v - 22.0d) / ((v - 6) * (v - 8));
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
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [1.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            DegreesOfFreedom = (int)parameters[0];
            Sigma = parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters[0] < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "The scale parameter σ (sigma) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "The scale parameter σ (sigma) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x < Minimum) return 0.0d;
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            double a = Math.Pow(t2 * v / 2d, v / 2d);
            double b = Gamma.Function(v / 2d);
            double c = Math.Exp(-v * t2 / (2d * x));
            double d = Math.Pow(x, 1d + v / 2d);
            return a / b * (c / d);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            if (x <= Minimum)
                return 0d;
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            return Gamma.UpperIncomplete(v / 2.0d, v * t2 / 2.0d * x);
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
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters([DegreesOfFreedom, Sigma], true);
            double v = DegreesOfFreedom;
            double t2 = Sigma;
            return  Gamma.InverseUpperIncomplete(v / 2.0d, probability) / (t2 * v / 2.0d);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new InverseChiSquared(DegreesOfFreedom, Sigma);
        }

    }
}