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
    /// The Chi-Squared (χ²) probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para> 
    /// <b> Description:</b>
    /// </para>
    /// <para>
    /// In probability theory and statistics, the chi-square distribution (also chi-squared
    /// or χ²-distribution) with k degrees of freedom Is the distribution of a sum of the
    /// squares of k independent standard normal random variables. It Is one of the most
    /// widely used probability distributions in inferential statistics, e.g. in hypothesis
    /// testing, or in construction of confidence intervals.
    /// </para>
    /// <para>
    /// <b>References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Wikipedia contributors, "Chi-squared distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Chi-squared_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class ChiSquared : UnivariateDistributionBase
    {
       
        /// <summary>
        /// Constructs a Chi-Squared distribution with 10 degrees of freedom.
        /// </summary>
        public ChiSquared()
        {
            SetParameters(10d);
        }

        /// <summary>
        /// Constructs a Chi-Squared distribution with given degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom for the distribution.</param>
        public ChiSquared(int degreesOfFreedom)
        {
            SetParameters(degreesOfFreedom);
        }
       
        private int _degreesOfFreedom;
       
        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public int DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(value, false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 1; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.ChiSquared; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Chi-Squared (χ²)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "χ²"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(DegreesOfFreedom)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Convert.ToDouble(DegreesOfFreedom)]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return DegreesOfFreedom; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Math.Max(DegreesOfFreedom - 2.0d, 0.0d); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(2.0d * DegreesOfFreedom); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return Math.Sqrt(8.0d / DegreesOfFreedom); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + 12.0d / DegreesOfFreedom; }
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
            get { return [1.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity]; }
        }
     
        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        public void SetParameters(double v)
        {
            DegreesOfFreedom = (int)v;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom for the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(int degreesOfFreedom, bool throwException)
        {
            if (degreesOfFreedom < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }

            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters((int)parameters[0], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double X)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(DegreesOfFreedom, true);
            if (X < Minimum) return 0.0d;
            double v = DegreesOfFreedom;
            double a = Math.Pow(X, (v - 2.0d) / 2.0d);
            double b = Math.Exp(-X / 2.0d);
            double c = Math.Pow(2d, v / 2.0d) * Gamma.Function(v / 2.0d);
            return a * b / c;
        }

        /// <inheritdoc/>
        public override double CDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(DegreesOfFreedom, true);
            if (X <= Minimum)
                return 0d;
            return Gamma.LowerIncomplete(DegreesOfFreedom / 2.0d, X / 2.0d);
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
                ValidateParameters(DegreesOfFreedom, true);
            return Gamma.InverseLowerIncomplete(DegreesOfFreedom / 2.0d, probability) * 2.0d;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new ChiSquared(DegreesOfFreedom);
        }
      
    }
}