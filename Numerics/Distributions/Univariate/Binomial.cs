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
    /// The Binomial distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Binomial_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Binomial : UnivariateDistributionBase
    {

        
        /// <summary>
        /// Constructs a Binomial distribution with p=0.5 and n=10.
        /// </summary>
        public Binomial()
        {
            SetParameters([0.5d, 10d]);
        }

        /// <summary>
        /// Constructs a Binomial distribution with a given probability and sample size.
        /// </summary>
        /// <param name="probability">The success probability (p) in each trial. Range: 0 ≤ p ≤ 1.</param>
        /// <param name="numberOfTrials">The number of trials (n). Range: n ≥ 0.</param>
        public Binomial(double probability, int numberOfTrials)
        {
            SetParameters([probability, numberOfTrials]);
        }
     
        private double _probabilityOfSuccess;
        private int _numberOfTrials;

        /// <summary>
        /// Gets and sets the success probability in each trial. Range: 0 ≤ p ≤ 1.
        /// </summary>
        public double ProbabilityOfSuccess
        {
            get { return _probabilityOfSuccess; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, NumberOfTrials }, false) is null;
                _probabilityOfSuccess = value;
            }
        }

        /// <summary>
        /// Gets the complement of the probability.
        /// </summary>
        public double Complement
        {
            get { return 1d - ProbabilityOfSuccess; }
        }

        /// <summary>
        /// Gets and sets the number of trials. Range: n ≥ 0.
        /// </summary>
        public int NumberOfTrials
        {
            get { return _numberOfTrials; }
            set
            {
                _parametersValid = ValidateParameters(new[] { ProbabilityOfSuccess, value }, false) is null;
                _numberOfTrials = value;
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
            get { return UnivariateDistributionType.Binomial; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Binomial"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Bin"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Probability of Success (p)";
                parmString[1, 0] = "Number of Trials (n)";
                parmString[0, 1] = ProbabilityOfSuccess.ToString();
                parmString[1, 1] = NumberOfTrials.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["p", "n"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(ProbabilityOfSuccess), nameof(NumberOfTrials)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [ProbabilityOfSuccess, NumberOfTrials]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return NumberOfTrials * ProbabilityOfSuccess; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Math.Ceiling(NumberOfTrials * ProbabilityOfSuccess); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (ProbabilityOfSuccess == 1.0d)
                    return NumberOfTrials;
                if (ProbabilityOfSuccess == 0.0d)
                    return 0.0d;
                return (int)Math.Floor((NumberOfTrials + 1) * ProbabilityOfSuccess);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return (1.0d - 2.0d * ProbabilityOfSuccess) / Math.Sqrt(NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + (1.0d - 6d * Complement * ProbabilityOfSuccess) / (NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return NumberOfTrials; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [1.0d, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            ProbabilityOfSuccess = parameters[0];
            NumberOfTrials = (int)parameters[1];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            // Validate probability
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] < 0.0d || parameters[0] > 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "Probability must be between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "Probability must be between 0 and 1.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "The number of trials (n) must be positive.");
                return new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "The number of trials (n) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
            k = Math.Floor(k);
            if (k < Minimum || k > Maximum) return 0.0d;
            if (ProbabilityOfSuccess == 0.0d)
                return k == 0d ? 1.0d : 0.0d;
            if (ProbabilityOfSuccess == 1.0d)
                return k == NumberOfTrials ? 1.0d : 0.0d;
            return Factorial.BinomialCoefficient(NumberOfTrials, (int)k) * Math.Pow(ProbabilityOfSuccess, k) * Math.Pow(Complement, NumberOfTrials - k);
        }

        /// <inheritdoc/>
        public override double CDF(double k)
        {
            if (_parametersValid == false)
                ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
            k = Math.Floor(k);
            if (k < Minimum)
                return 0.0d;
            if (k > Maximum)
                return 1.0d;
            return Beta.Incomplete(NumberOfTrials - k, k + 1d, Complement);
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
                ValidateParameters([probability, NumberOfTrials], true);
            double k = 0d;
            for (int i = 0; i < NumberOfTrials; i++)
            {
                if (CDF(i) >= probability)
                {
                    k = i;
                    break;
                }
            }
            return k;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Binomial(ProbabilityOfSuccess, NumberOfTrials);
        }
     
    }
}