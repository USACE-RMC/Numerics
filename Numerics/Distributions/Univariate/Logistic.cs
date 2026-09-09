using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The logistic distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Logistic_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Logistic : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, IStandardError, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a Logistic distribution with a location of 0 and scale of 0.1.
        /// </summary>
        public Logistic()
        {
            SetParameters(0d, 0.1d);
        }

        /// <summary>
        /// Constructs an Logistic distribution with a given ξ and α.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public Logistic(double location, double scale)
        {
            SetParameters(location, scale);
        }

        private double _xi; // location
        private double _alpha; // scale

        /// <summary>
        /// Gets and sets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters([value, Alpha], false) is null;
                _xi = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters([Xi, value], false) is null;
                _alpha = value;
            }
        }

        /// <inheritdoc/>
        /// <remarks>The symmetric logit form avoids exponential overflow in either tail.</remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            double magnitude = Math.Abs(DistributionNumerics.Standardize(x, Xi, Alpha));
            return -Math.Log(Alpha) - magnitude - 2d * Tools.Log1p(Math.Exp(-magnitude));
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            double z = DistributionNumerics.Standardize(x, Xi, Alpha);
            return z >= 0d ? -Tools.Log1p(Math.Exp(-z)) : z - Tools.Log1p(Math.Exp(z));
        }

        /// <inheritdoc/>
        public override double CCDF(double x)
        {
            return Math.Exp(LogCCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            double z = DistributionNumerics.Standardize(x, Xi, Alpha);
            return z >= 0d ? -z - Tools.Log1p(Math.Exp(-z)) : -Tools.Log1p(Math.Exp(z));
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Logistic; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Logistic"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LO"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Location (ξ)";
                parmString[1, 0] = "Scale (α)";
                parmString[0, 1] = Xi.ToString();
                parmString[1, 1] = Alpha.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ξ", "α"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xi), nameof(Alpha)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xi, Alpha]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Xi; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Xi; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Alpha * Math.PI / Math.Sqrt(3d); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + 6d / 5d; }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(ParametersFromMoments(Statistics.ProductMoments(sample)));
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
            var newDistribution = new Logistic(Xi, Alpha);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public void SetParameters(double location, double scale)
        {
            _xi = location;
            _alpha = scale;
            _parametersValid = ValidateParameters(new[] { location, scale }, false) is null;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var parms = new double[NumberOfParameters];
            parms[0] = moments[0];
            parms[1] = moments[1] * Math.Sqrt(3d) / Math.PI;
            return parms;
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new Logistic();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            return DistributionNumerics.PreferLegacyConstraints(
                () => GetLegacyParameterConstraints(sample), () => GetRobustParameterConstraints(sample));
        }

        /// <summary>Preserves the established initialization and family-specific prior envelope.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>The legacy initial values and lower and upper bounds.</returns>
        private Tuple<double[], double[], double[]> GetLegacyParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            initialVals = ParametersFromMoments(Statistics.ProductMoments(sample));
            // Get bounds of location
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            var normal = new Normal().GetRobustParameterConstraints(sample);
            double correction = Math.Sqrt(3d) / Math.PI;
            normal.Item1[1] *= correction;
            normal.Item2[1] *= correction;
            normal.Item3[1] *= correction;
            return normal;
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var LO = new Logistic();
                LO.SetParameters(x);
                return LO.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            double magnitude = Math.Abs(DistributionNumerics.Standardize(x, Xi, Alpha));
            if (magnitude > 36d) return Math.Exp(LogPDF(x));
            double tail = Math.Exp(-magnitude);
            return (tail / (1d + tail)) / (1d + tail) / Alpha;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return Math.Exp(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (double.IsNaN(probability) || probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, Alpha], true);
            return Xi + Alpha * (Math.Log(probability) - Tools.Log1p(-probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Logistic(Xi, Alpha);
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments &&
                estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, _alpha], true);

            // Compute covariance
            double a = Alpha / Math.Sqrt(sampleSize);
            var covar = new double[2, 2];
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                covar[0, 0] = Math.PI * Math.PI / 3d * (a * a); // location
                covar[1, 1] = 4d / 5d * (a * a); // scale
                covar[0, 1] = 0.0; // location & scale
                covar[1, 0] = covar[0, 1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                covar[0, 0] = 3d * (a * a); // location
                covar[1, 1] = 9d / (3d + Math.PI * Math.PI) * (a * a); // scale
                covar[0, 1] = 0.0; // location & scale
                covar[1, 0] = covar[0, 1];
            }
            return covar;
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            var grad = QuantileGradient(probability);
            var covariance = new Logistic(0, 1).ParameterCovariance(sampleSize, estimationMethod);
            return DistributionNumerics.ScaledQuantileVariance(covariance, grad, Alpha);
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, _alpha], true);
            var gradient = new double[]
            {
                1.0d, // location
                (Math.Log(probability) - Tools.Log1p(-probability)) // scale
            };
            return gradient;
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

    }
}
