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

using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// A general truncated probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class TruncatedDistribution : UnivariateDistributionBase
    {

        /// <summary>
        /// Constructs a truncated distribution. 
        /// </summary>
        /// <param name="basDistribution">The base distribution to truncate.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public TruncatedDistribution(UnivariateDistributionBase basDistribution, double min, double max)
        {
            _baseDist = basDistribution;
            _min = min;
            _max = max;
            _Fmin = _baseDist.CDF(_min);
            _Fmax = _baseDist.CDF(_max);
            _momentsComputed = false;
        }

        private readonly UnivariateDistributionBase _baseDist;
        private double _min = double.NegativeInfinity;
        private double _max = double.PositiveInfinity;
        private double _Fmin, _Fmax;
        private bool _momentsComputed = false;
        private double[] u = [double.NaN, double.NaN, double.NaN, double.NaN];

        /// <summary>
        /// Gets the base distribution. 
        /// </summary>
        public UnivariateDistributionBase BaseDistribution => _baseDist;

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min => _min;

        /// <summary>
        /// Get and set the max of the distribution.
        /// </summary>
        public double Max => _max;

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return _baseDist.NumberOfParameters + 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => _baseDist.Type;

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Truncated " + _baseDist.DisplayName; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Trunc. " +_baseDist.ShortDisplayName; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[_baseDist.NumberOfParameters + 2, 2];
                for (int i = 0; i < _baseDist.NumberOfParameters; i++)
                {
                    parmString[i, 0] = _baseDist.ParametersToString[i, 0];
                    parmString[i, 1] = _baseDist.ParametersToString[i, 1];
                }
                parmString[NumberOfParameters - 2, 0] = "Min";
                parmString[NumberOfParameters - 1, 0] = "Max";
                parmString[NumberOfParameters - 2, 1] = Min.ToString();
                parmString[NumberOfParameters - 1, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get
            {
                var parms = _baseDist.ParameterNamesShortForm.ToList();
                parms.AddRange(new[] { "Min", "Max" });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get
            {
                var parms = _baseDist.GetParameterPropertyNames.ToList();
                parms.AddRange(new[] { nameof(Min), nameof(Max) });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get
            {
                var parms = _baseDist.GetParameters.ToList();
                parms.AddRange(new[] { Min, Max });                
                return parms.ToArray(); 
            }
        }


        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[0];
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
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[1];
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[2];
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[3];
            }
        }

        /// <inheritdoc/>
        public override double Minimum => Math.Max(_baseDist.Minimum, Min);

        /// <inheritdoc/>
        public override double Maximum => Math.Min(_baseDist.Maximum, Max);

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get 
            {
                var parms = _baseDist.MinimumOfParameters.ToList();
                parms.AddRange(new[] { double.NegativeInfinity, Min });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get
            {
                var parms = _baseDist.MaximumOfParameters.ToList();
                parms.AddRange(new[] { Max, double.PositiveInfinity });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            _baseDist.SetParameters(parameters.ToArray().Subset(parameters.Count - 2));
            _min = parameters[parameters.Count - 2];
            _max = parameters[parameters.Count - 1];
            _Fmin = _baseDist.CDF(_min);
            _Fmax = _baseDist.CDF(_max);
            _momentsComputed = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (_baseDist != null!) _baseDist.ValidateParameters(parameters.ToArray().Subset(0, parameters.Count - 2), throwException);
            if (double.IsNaN(Min) || double.IsNaN(Max) || double.IsInfinity(Min) || double.IsInfinity(Max) || Min >= Max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min must be less than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min must be less than the max.");
            }
            if (_Fmin == _Fmax)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "Truncation interval has zero probability mass.");
                return new ArgumentOutOfRangeException(nameof(Min), "Truncation interval has zero probability mass.");
            }              
            return null!;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);
            if (x < Min || x > Max) return 0.0;
            return _baseDist.PDF(x) / (_Fmax - _Fmin);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);
            if (x <= Min) return 0.0;
            if (x >= Max) return 1.0;
            return (_baseDist.CDF(x) - _Fmin) / (_Fmax - _Fmin);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(GetParameters, true);
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            return _baseDist.InverseCDF(probability * (_Fmax - _Fmin) + _Fmin);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new TruncatedDistribution(_baseDist, Min, Max);
        }

    }
}
