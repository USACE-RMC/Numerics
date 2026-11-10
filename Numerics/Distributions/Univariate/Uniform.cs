using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Uniform probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)#Probability_density_function" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Uniform : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a Uniform distribution with min = 0 and max = 1.
        /// </summary>
        public Uniform()
        {
            SetParameters(0d, 1d);
        }

        /// <summary>
        /// Constructs a Uniform distribution with specified min and max.
        /// </summary>
        /// <param name="min">The min of the distribution.</param>
        /// <param name="max">The max of the distribution.</param>
        public Uniform(double min, double max)
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

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Uniform; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Uniform"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "U"; }
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
            get { return (Max - Min) / Math.Sqrt(12.0d); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 9d / 5d; }
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
        /// <remarks>
        /// Maximum likelihood gives Min = min(x) and Max = max(x). The method of moments gives
        /// mean(x) ∓ √3·s, where s is the sample standard deviation with divisor n − 1.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite value, or is constant.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not supported.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 2);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(ParametersFromMoments([Statistics.Mean(sample), Statistics.StandardDeviation(sample)]));
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
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
            var newDistribution = new Uniform(Min, Max);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the mean and standard deviation, the first two entries of <paramref name="moments"/>: the
        /// parameters are mean ∓ √3·standard deviation.
        /// </remarks>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            double halfWidth = Math.Sqrt(3d) * moments[1];
            return [moments[0] - halfWidth, moments[0] + halfWidth];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new Uniform(parameters[0], parameters[1]);
            return [dist.Mean, dist.StandardDeviation, dist.Skewness, dist.Kurtosis];
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The initial values are the closed-form maximum likelihood estimates. Min is bounded above by the sample
        /// minimum and Max is bounded below by the sample maximum, because every observation must lie in
        /// [Min, Max]; the outer bounds are scale-aware finite limits.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite value, or is constant.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = MLE(sample);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            double range = initialVals[1] - initialVals[0];
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], range, initialVals[0], initialVals[1],
                true, out lowerVals[0], out upperVals[0]);
            DistributionNumerics.LocationParameterBounds(ref initialVals[1], range, initialVals[0], initialVals[1],
                false, out _, out upperVals[1]);
            lowerVals[1] = initialVals[1];
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The maximum likelihood estimates have the closed form Min = min(x) and Max = max(x); they equal
        /// scipy.stats.uniform.fit, whose location and scale are Min and Max − Min.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than two observations, contains a
        /// nonfinite value, or is constant.</exception>
        public double[] MLE(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2);
            return [Statistics.Minimum(sample), Statistics.Maximum(sample)];
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
        public ArgumentOutOfRangeException? ValidateParameters(double min, double max, bool throwException)
        {
            if (double.IsNaN(min) || double.IsInfinity(min) ||
                double.IsNaN(max) || double.IsInfinity(max) || min > max)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
                return new ArgumentOutOfRangeException(nameof(Min), "The min cannot be greater than the max.");
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
        /// When Min = Max the distribution is a point mass, which has no density, and the PDF is 0 everywhere.
        /// A not-a-number argument returns not-a-number.
        /// </remarks>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (double.IsNaN(x)) return double.NaN;
            // A point mass has no density.
            if (Min == Max) return 0d;
            if (x < Minimum || x > Maximum) return 0.0d;
            return 1.0d / (Max - Min);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// When Min = Max the distribution is a point mass at Min, so the CDF is 0 below Min and 1 from Min upward.
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (Min == Max) return double.IsNaN(x) ? double.NaN : x < Min ? 0d : 1d;
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return 1d;
            return (x - Min) / (Max - Min);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// When Min = Max every probability maps to Min.
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Min, Max, true);
            if (Min == Max) return Min;
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            return Min + probability * (Max - Min);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Uniform(Min, Max);
        }

    }
}