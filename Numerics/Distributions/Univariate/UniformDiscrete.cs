using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;

namespace Numerics.Distributions
{

    /// <summary>
    /// The uniform (discrete) distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Discrete_uniform_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class UniformDiscrete : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
   
        /// <summary>
        /// Constructs a Uniform (discrete) distribution with min = 0 and max = 1.
        /// </summary>
        public UniformDiscrete()
        {
            SetParameters(0d, 1d);
        }

        /// <summary>
        /// Constructs a Uniform (discrete) distribution with specified min and max.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        public UniformDiscrete(double min, double max)
        {
            SetParameters(min, max);
        }

        private double _min;
        private double _max;

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set
            {
                _parametersValid = ValidateParameters(value, Max, false) is null;
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
                _parametersValid = ValidateParameters(Min, value, false) is null;
                _max = value;
            }
        }

        /// <summary>
        /// Gets the number of support points, Max - Min + 1.
        /// </summary>
        public double N
        {
            get { return Max - Min + 1d; }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.UniformDiscrete; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Uniform (Discrete)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "UD"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Min";
                parmString[1, 0] = "Max";
                parmString[0, 1] = Min.ToString();
                parmString[1, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["Min", "Max"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Min), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Min, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return (Min + Max) / 2.0d; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return (Min + Max) / 2.0d; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return double.NaN; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                double n = Max - Min + 1;
                return Math.Sqrt((n * n - 1.0d) / 12.0d);
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Undefined, and returned as not-a-number, for the point mass Min = Max, whose variance is zero.
        /// </remarks>
        public override double Kurtosis
        {
            get
            {
                if (Min == Max) return double.NaN;
                return 3.0d - 6d * (N * N + 1d) / (5d * (N * N - 1d));
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Min; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Max; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, Min]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [Max, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty or contains a value that is not a
        /// finite integer.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not
        /// <see cref="ParameterEstimationMethod.MaximumLikelihood"/>.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            ValidateSample(sample);
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new UniformDiscrete(Min, Max);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The initial values are the closed-form maximum likelihood estimates. Min is bounded above by the sample
        /// minimum and Max is bounded below by the sample maximum, because every observation must lie in
        /// [Min, Max]; the outer bounds are scale-aware finite limits.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty or contains a value that is not a
        /// finite integer.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = MLE(sample);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // A constant sample is a valid point mass; one support step stands in for its zero range.
            double scale = Math.Max(1d, initialVals[1] - initialVals[0]);
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], scale, initialVals[0], initialVals[1],
                true, out lowerVals[0], out upperVals[0]);
            DistributionNumerics.LocationParameterBounds(ref initialVals[1], scale, initialVals[0], initialVals[1],
                false, out _, out upperVals[1]);
            lowerVals[1] = initialVals[1];
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The maximum likelihood estimates have the closed form Min = min(x) and Max = max(x), the narrowest support
        /// that contains every observation.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty or contains a value that is not a
        /// finite integer.</exception>
        public double[] MLE(IList<double> sample)
        {
            ValidateSample(sample);
            return [Statistics.Minimum(sample), Statistics.Maximum(sample)];
        }

        /// <summary>
        /// Checks that a sample can be used to estimate the support.
        /// </summary>
        /// <param name="sample">The observations to validate.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty or contains a value that is not a
        /// finite integer.</exception>
        /// <remarks>
        /// The support is a range of integers, so every observation must be a finite integer. A constant sample
        /// estimates a valid point mass and is accepted.
        /// </remarks>
        private static void ValidateSample(IList<double> sample)
        {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new ArgumentOutOfRangeException(nameof(sample), "Insufficient observations to initialize the distribution.");
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]) || sample[i] != Math.Floor(sample[i]))
                    throw new ArgumentOutOfRangeException(nameof(sample), "Observations must be finite integers.");
            }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        public void SetParameters(double min, double max)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(min, max, false) is null;
            // Set parameters
            _min = min;
            _max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        /// <returns><see langword="null"/> when the parameters are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="throwException"/> is <see langword="true"/>
        /// and the min or max is not a finite integer, or the min is greater than the max.</exception>
        /// <remarks>
        /// The support is the integers from min through max, so both must be integers. Min = Max is a valid point mass.
        /// </remarks>
        public ArgumentOutOfRangeException? ValidateParameters(double min, double max, bool throwException)
        {
            if (double.IsNaN(min) || double.IsInfinity(min) ||
                double.IsNaN(max) || double.IsInfinity(max) || min > max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
            }
            if (min != Math.Floor(min) || max != Math.Floor(max))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min and max must be integers.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min and max must be integers.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The mass is 1/N at each integer from Min through Max and zero elsewhere, including between the integers.
        /// A not-a-number argument returns not-a-number.
        /// </remarks>
        public override double PDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (double.IsNaN(k)) return double.NaN;
            if (k < Minimum || k > Maximum) return 0.0d;
            // The probability mass is zero away from the integers.
            if (k != Math.Floor(k)) return 0.0d;
            return 1d / N;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Non-integer k is floored, so the CDF is the step function (⌊k⌋ − Min + 1)/N on [Min, Max].
        /// </remarks>
        public override double CDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (k < Minimum)
                return 0d;
            if (k > Maximum)
                return 1d;
            return (Math.Floor(k) - Min + 1.0d) / N;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Returns the smallest integer k with CDF(k) ≥ p, which is Min + ⌈p·N⌉ − 1, checked against the CDF so
        /// that rounding in p·N cannot move the result by one.
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // j is the number of support points at or below the quantile, so CDF(Min + j - 1) = j/N.
            double j = Tools.Clamp(Math.Ceiling(probability * N), 1d, N);
            while (j > 1d && (j - 1d) / N >= probability)
                j--;
            while (j < N && j / N < probability)
                j++;
            return Min + j - 1d;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new UniformDiscrete(Min, Max);
        }
   
    }
}