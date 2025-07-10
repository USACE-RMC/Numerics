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
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Poisson distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Poisson_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Poisson : UnivariateDistributionBase
    {
 
        /// <summary>
        /// Constructs a Poisson distribution with λ = 1.
        /// </summary>
        public Poisson()
        {
            SetParameters([1.0d]);
        }

        /// <summary>
        /// Constructs a Poisson distribution with a given rate λ (lambda).
        /// </summary>
        /// <param name="rate">The rate parameter λ (lambda). Range: λ > 0.</param>
        public Poisson(double rate)
        {
            SetParameters([rate]);
        }

        private double _lambda;

        /// <summary>
        /// Gets and sets the Poisson's rate parameter λ (lambda). Range: λ > 0.
        /// </summary>
        public double Lambda
        {
            get { return _lambda; }

            set
            {
                _parametersValid = ValidateParameters([value], false) is null;
                _lambda = value;
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
            get { return UnivariateDistributionType.Poisson; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Poisson"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "P"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                parmString[0, 0] = "Rate (λ)";
                parmString[0, 1] = Lambda.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["λ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Lambda)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Lambda]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Lambda; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Math.Floor(Lambda); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(Lambda); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 1.0d / Math.Sqrt(Lambda); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3.0d + 1.0d / Lambda; }
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
            get { return [0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            Lambda = parameters[0];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Lambda), "The rate (λ) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Lambda), "The rate (λ) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            k = Math.Floor(k);
            if (k < Minimum || k > Maximum) return 0.0d;
            return Math.Exp(-Lambda + k * Math.Log(Lambda) - Factorial.LogFactorial((int)k));
        }

        /// <inheritdoc/>
        public override double CDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            k = Math.Floor(k);
            if (k < Minimum)
                return 0.0d;
            if (k > Maximum)
                return 1.0d;
            return Gamma.UpperIncomplete(k + 1d, Lambda);
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
                ValidateParameters([Lambda], true);

            int k = (int)Math.Max(0, Math.Floor(Lambda));  // start near the mean

            // Move downward if needed
            while (k > 0 && CDF(k - 1) >= probability)
                k--;

            // Move upward to find first k such that CDF(k) >= probability
            while (CDF(k) < probability)
                k++;

            return k;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Poisson(Lambda);
        }
    
    }
}