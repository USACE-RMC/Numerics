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
            var parameters = _baseDist.GetParameters.ToList();
            parameters.Add(min);
            parameters.Add(max);
            _parametersValid = ValidateParametersCore(parameters, false, out _Fmin, out _Fmax) is null;
            _momentsComputed = false;
        }

        private readonly UnivariateDistributionBase _baseDist;
        private double _min;
        private double _max;
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
        /// <exception cref="ArgumentException">Thrown when the flattened parameter count does not match the base distribution and two truncation bounds.</exception>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }
            var validation = ValidateParametersCore(parameters, false, out double lowerProbability, out double upperProbability);
            _baseDist.SetParameters(parameters.ToArray().Subset(0, parameters.Count - 3));
            _min = parameters[parameters.Count - 2];
            _max = parameters[parameters.Count - 1];
            _Fmin = lowerProbability;
            _Fmax = upperProbability;
            _parametersValid = validation is null;
            _momentsComputed = false;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParametersCore(parameters, throwException, out _, out _);
        }

        /// <summary>
        /// Validates a flattened base-distribution parameter vector and its truncation bounds.
        /// </summary>
        /// <param name="parameters">The base parameters followed by the lower and upper truncation bounds.</param>
        /// <param name="throwException">Whether to throw the first validation exception.</param>
        /// <param name="lowerProbability">The base-distribution CDF at the lower bound when valid; otherwise <see cref="double.NaN"/>.</param>
        /// <param name="upperProbability">The base-distribution CDF at the upper bound when valid; otherwise <see cref="double.NaN"/>.</param>
        /// <returns><see langword="null"/> when all parameters are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="throwException"/> is true and validation fails.</exception>
        private ArgumentOutOfRangeException? ValidateParametersCore(IList<double> parameters, bool throwException, out double lowerProbability, out double upperProbability)
        {
            lowerProbability = double.NaN;
            upperProbability = double.NaN;
            if (parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "The flattened parameter count is invalid.");
                if (throwException) throw exception;
                return exception;
            }
            var baseParameters = parameters.ToArray().Subset(0, parameters.Count - 3);
            var baseValidation = _baseDist.ValidateParameters(baseParameters, false);
            if (baseValidation is not null)
            {
                if (throwException) _baseDist.ValidateParameters(baseParameters, true);
                return baseValidation;
            }
            double min = parameters[parameters.Count - 2];
            double max = parameters[parameters.Count - 1];
            if (double.IsNaN(min) || double.IsNaN(max) || double.IsInfinity(min) || double.IsInfinity(max) || min >= max)
            {
                var exception = new ArgumentOutOfRangeException(nameof(Min), "The min must be finite and less than the max.");
                if (throwException) throw exception;
                return exception;
            }
            var candidate = _baseDist.Clone();
            candidate.SetParameters(baseParameters);
            lowerProbability = candidate.CDF(min);
            upperProbability = candidate.CDF(max);
            if (Math.Abs(lowerProbability - upperProbability) < 1e-15)
            {
                var exception = new ArgumentOutOfRangeException(nameof(Min), "Truncation interval has zero probability mass.");
                if (throwException) throw exception;
                return exception;
            }
            return null;
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
