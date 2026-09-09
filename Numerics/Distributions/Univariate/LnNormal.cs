using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Ln-Normal (Galton) probability distribution.
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
    /// Wikipedia contributors, "Log-normal distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Log-normal_distribution"/>
    /// </para>
    /// <para>
    /// <see cref="SetParameters(double, double)"/> and <see cref="GetParameters"/> use physical
    /// mean and standard deviation; <see cref="Mu"/> and <see cref="Sigma"/> store the natural-log
    /// mean and standard deviation. Quantile derivatives and covariance use the physical coordinates.
    /// The method-of-moments estimator is the existing indirect estimator of log observations,
    /// so its leading covariance is the transformed Normal covariance, as for maximum likelihood.
    /// <see cref="ParametersFromMoments"/> instead converts explicitly supplied physical moments.
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class LnNormal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
    
        /// <summary>
        /// Constructs a Ln-Normal distribution with a mean of 10 and standard deviation of 10.
        /// </summary>
        public LnNormal()
        {
            SetParameters(10d, 10d);
        }

        /// <summary>
        /// Constructs a Ln-Normal (Galton) distribution with given mean and standard deviation.
        /// </summary>
        /// <param name="mean">The mean of the distribution.</param>
        /// <param name="standardDeviation">The standard deviation of the distribution.</param>
        /// <remarks>
        /// Enter the real-space mean and standard deviation of the distribution. The two parameters μ and σ are not
        /// location and scale parameters for a log-normally distributed random variable X, but they are respectively
        /// location and scale parameters for the normally distributed logarithm ln(X).
        /// </remarks>
        public LnNormal(double mean, double standardDeviation)
        {
            SetParameters(mean, standardDeviation);
        }

        private double _mu;
        private double _sigma;
        private bool _hasPhysicalMoments;
        private double _physicalMean;
        private double _physicalStandardDeviation;

        /// <summary>Whether the stored physical-moment surface is active, for snapshot capture without computing moments.</summary>
        internal bool PhysicalMomentModeForSnapshot => _hasPhysicalMoments;

        /// <summary>The stored physical mean field, for snapshot capture without computing moments.</summary>
        internal double PhysicalMeanForSnapshot => _physicalMean;

        /// <summary>The stored physical standard-deviation field, for snapshot capture without computing moments.</summary>
        internal double PhysicalStandardDeviationForSnapshot => _physicalStandardDeviation;

        /// <summary>
        /// Gets and sets the mean µ (Mu) of the natural logarithm of the observation.
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateLogParameters(value, Sigma, false) is null;
                _mu = value;
                _hasPhysicalMoments = false;
            }
        }

        /// <summary>
        /// Gets and sets the standard deviation σ (sigma) of the natural logarithm of the observation.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateLogParameters(Mu, value, false) is null;
                _sigma = value;
                _hasPhysicalMoments = false;
            }
        }

        /// <summary>Validates the internal natural-log coordinates without applying physical-mean constraints.</summary>
        /// <param name="mean">The proposed finite natural-log mean.</param>
        /// <param name="standardDeviation">The proposed finite positive natural-log standard deviation.</param>
        /// <param name="throwException"><see langword="true"/> to throw the validation error; <see langword="false"/> to return it.</param>
        /// <returns><see langword="null"/> when both log-coordinate parameters are valid; otherwise, the corresponding validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="throwException"/> is <see langword="true"/> and either log-coordinate parameter is invalid.</exception>
        private static ArgumentOutOfRangeException? ValidateLogParameters(double mean, double standardDeviation, bool throwException)
        {
            ArgumentOutOfRangeException? error = null;
            if (double.IsNaN(mean) || double.IsInfinity(mean))
                error = new ArgumentOutOfRangeException(nameof(Mu), "The logarithmic mean must be finite.");
            else if (!(standardDeviation > 0d) || double.IsInfinity(standardDeviation))
                error = new ArgumentOutOfRangeException(nameof(Sigma), "The logarithmic standard deviation must be finite and positive.");
            if (throwException && error != null) throw error;
            return error;
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            return x <= 0d ? double.NegativeInfinity : DistributionNumerics.NormalLogCDF(DistributionNumerics.Standardize(Math.Log(x), Mu, Sigma));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            return x <= 0d ? 0d : DistributionNumerics.NormalLogSurvival(DistributionNumerics.Standardize(Math.Log(x), Mu, Sigma));
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.LnNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Log-Normal (base e)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LN"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[0, 1] = Mean.ToString();
                parmString[1, 1] = StandardDeviation.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mean), nameof(StandardDeviation)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mean, StandardDeviation]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return _hasPhysicalMoments ? _physicalMean : Math.Exp(Mu + Sigma * Sigma / 2.0d); }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Math.Exp(Mu); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Math.Exp(Mu - Sigma * Sigma); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (_hasPhysicalMoments) return _physicalStandardDeviation;
                double variance = Sigma * Sigma;
                double logExcess = variance > 0.5d ? variance + Tools.Log1p(-Math.Exp(-variance)) : Math.Log(Tools.Expm1(variance));
                return Math.Exp(Mu + 0.5d * variance + 0.5d * logExcess);
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return (Math.Exp(Sigma * Sigma) + 2d) * Math.Sqrt(Tools.Expm1(Sigma * Sigma)); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double variance = Sigma * Sigma;
                return 3d + Tools.Expm1(4d * variance) + 2d * Tools.Expm1(3d * variance) + 3d * Tools.Expm1(2d * variance);
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
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                // Estimate using the method of moments(a.k.a product moments).
                var parms = IndirectMethodOfMoments(sample);
                Mu = parms[0];
                Sigma = parms[1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                var parms = ParametersFromLinearMoments(IndirectMethodOfLinearMoments(sample));
                Mu = parms[0];
                Sigma = parms[1];
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
            // Create a new distribution and estimate parameters from the bootstrap sample 
            var newDistribution = new LnNormal() { Mu = Mu, Sigma = Sigma };
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters using the "direct method."
        /// </summary>
        /// <param name="mean">The mean of the distribution.</param>
        /// <param name="standardDeviation">The standard deviation of the distribution.</param>
        /// <remarks>
        /// The direct method for setting parameters is used so that users can set the parameters directly
        /// from real-space data, which is more intuitive.
        /// Physical inputs are retained exactly for parameter-vector and XML round trips until a log
        /// parameter changes. The legacy minimum log-scale rule still applies when conversion reaches it.
        /// </remarks>
        public void SetParameters(double mean, double standardDeviation)
        {
            var parms = DirectMethodOfMoments(mean, standardDeviation);
            // Validate parameters
            Mu = parms[0];
            Sigma = parms[1];
            if (_parametersValid && Sigma == parms[1])
            {
                _physicalMean = mean;
                _physicalStandardDeviation = standardDeviation;
                _hasPhysicalMoments = true;
            }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validates physical mean and standard deviation supplied through the parameter-vector API.
        /// </summary>
        /// <param name="mean">Mean.</param>
        /// <param name="standardDeviation">Standard deviation.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        /// <returns>The validation error, or null when both physical moments are finite and positive.</returns>
        /// <exception cref="ArgumentOutOfRangeException">A physical moment is invalid and <paramref name="throwException"/> is true.</exception>
        public ArgumentOutOfRangeException? ValidateParameters(double mean, double standardDeviation, bool throwException)
        {
            if (!(mean > 0d) || double.IsInfinity(mean))
            {
                var error = new ArgumentOutOfRangeException(nameof(mean), "The physical mean must be finite and positive.");
                if (throwException) throw error;
                return error;
            }
            if (!(standardDeviation > 0d) || double.IsInfinity(standardDeviation))
            {
                var error = new ArgumentOutOfRangeException(nameof(standardDeviation), "The physical standard deviation must be finite and positive.");
                if (throwException) throw error;
                return error;
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }
       
        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <returns>The product moments of the natural-log-transformed observations.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is insufficient, constant, nonfinite, or contains a nonpositive observation.</exception>
        public static double[] IndirectMethodOfMoments(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformedSample[i] = Math.Log(sample[i]);
            return Statistics.ProductMoments(transformedSample);
        }

        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <returns>The linear moments of the natural-log-transformed observations.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is insufficient, constant, nonfinite, or contains a nonpositive observation.</exception>
        public double[] IndirectMethodOfLinearMoments(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformedSample[i] = Math.Log(sample[i]);
            return Statistics.LinearMoments(transformedSample);
        }

        /// <summary>
        /// Sets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="mean">The real-space mean of the data.</param>
        /// <param name="standardDeviation">The real-space standard deviation of the data.</param>
        /// <returns>The natural-log mean and standard deviation, or two not-a-number values when the physical moments are invalid.</returns>
        public static double[] DirectMethodOfMoments(double mean, double standardDeviation)
        {
            if (!(mean > 0d) || !(standardDeviation > 0d) || double.IsInfinity(mean) || double.IsInfinity(standardDeviation))
                return [double.NaN, double.NaN];
            double logRatio = Math.Log(standardDeviation) - Math.Log(mean);
            double variance = logRatio > 0d
                ? 2d * logRatio + Tools.Log1p(Math.Exp(-2d * logRatio))
                : Tools.Log1p(Math.Exp(2d * logRatio));
            return [Math.Log(mean) - 0.5d * variance, Math.Sqrt(variance)];
        }

        /// <inheritdoc/>
        /// <remarks>Returns the supplied physical mean and standard deviation in the same coordinates used by <see cref="SetParameters(double, double)"/>.</remarks>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            ValidateParameters(moments[0], moments[1], true);
            return [moments[0], moments[1]];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new LnNormal();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            double mu = moments[0];
            double sigma = moments[1] * Math.Sqrt(Math.PI);
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double L1 = parameters[0];
            double L2 = parameters[1] * Math.Pow(Math.PI, -0.5);
            double T3 = 0d;
            double T4 = 30d * Math.Pow(Math.PI, -1d) * Math.Atan(Tools.Sqrt2) - 9d;
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        /// <remarks>Preserves usable legacy physical-moment initialization and rounded prior bounds,
        /// including samples containing nonpositive observations. Density support is unchanged.</remarks>
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
            var moments = Statistics.ProductMoments(sample);
            initialVals[0] = moments[0];
            initialVals[1] = moments[1];
            // Get bounds of mean
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of standard deviation
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var constraints = new Normal().GetRobustParameterConstraints(sample);
            constraints.Item2[0] = Math.Min(Tools.DoubleMachineEpsilon, constraints.Item1[0] / 10d);
            return constraints;
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
                var LN = new LnNormal();
                LN.SetParameters(x);
                return LN.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            return Math.Exp(LogPDF(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated in log space, so far-tail densities that underflow <see cref="PDF(double)"/>
        /// keep a finite log density.
        /// </remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            if (x <= 0d || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            double logX = Math.Log(x);
            double z = DistributionNumerics.Standardize(logX, Mu, Sigma);
            return -0.5d * z * z - Math.Log(Sigma) - Math.Log(Tools.Sqrt2PI) - logX;
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
                ValidateLogParameters(Mu, Sigma, true);
            return Math.Exp(Mu + Sigma * Normal.StandardZ(probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var clone = new LnNormal() { Mu = Mu, Sigma = Sigma };
            clone._hasPhysicalMoments = _hasPhysicalMoments;
            clone._physicalMean = _physicalMean;
            clone._physicalStandardDeviation = _physicalStandardDeviation;
            return clone;
        }

        /// <inheritdoc/>
        /// <remarks>The returned covariance is for physical mean and standard deviation. Both supported estimators use log observations; this method applies the complete two-coordinate delta transformation.</remarks>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments && estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException();
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            // Both supported estimators operate on log observations. Transform the leading
            // Normal covariance diag(sigma^2, sigma^2/2)/n into physical (mean, SD) coordinates.
            double mean = Mean, sd = StandardDeviation, variance = Sigma * Sigma;
            double excess = Tools.Expm1(variance);
            double sdDerivative = excess < 0.5d
                ? mean * (Sigma / Math.Sqrt(excess)) * (1d + 2d * excess)
                : sd * Sigma * (2d + 1d / excess);
            double seMean = Sigma / Math.Sqrt(sampleSize), seSd = Sigma / Math.Sqrt(2d * sampleSize);
            double m0 = mean * seMean, m1 = mean * Sigma * seSd;
            double s0 = sd * seMean, s1 = sdDerivative * seSd;
            return new[,] { { m0 * m0 + m1 * m1, m0 * s0 + m1 * s1 }, { m0 * s0 + m1 * s1, s0 * s0 + s1 * s1 } };
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments && estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException();
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            // The two physical-coordinate chain rules cancel to the original log-estimator delta method.
            double z = Normal.StandardZ(probability);
            double logStandardError = Mu + Sigma * z + Math.Log(Sigma) - .5 * Math.Log(sampleSize);
            return Math.Exp(2 * logStandardError + Tools.Log1p(.5 * z * z));
        }

        /// <inheritdoc/>
        /// <remarks>Returns derivatives of the physical quantile with respect to physical mean and standard deviation, including both log-mean and log-variance dependencies.</remarks>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateLogParameters(Mu, Sigma, true);
            double z = Normal.StandardZ(probability), variance = Sigma * Sigma;
            double weight = -Tools.Expm1(-variance);
            double relativeQuantile = Math.Exp(z * Sigma - 0.5d * variance);
            double meanDerivative = relativeQuantile * (1d + weight - z * weight / Sigma);
            double sdDerivative = relativeQuantile * Math.Sqrt(weight) * Math.Exp(-0.5d * variance) * (z / Sigma - 1d);
            return [meanDerivative, sdDerivative];
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

        /// <inheritdoc/>
        public override double[] ConditionalMoments(double a, double b)
        {
            if (a >= b)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Log-space Normal parameters: Y = ln X ~ N(mu, sigma^2)
            double mu = Mu;
            double sigma = Sigma;
            if (!(sigma > 0.0))
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Map bounds to log-space [A, B], allowing a <= 0 (A = -inf) and b = +inf
            double A = (a > 0.0) ? Math.Log(a) : double.NegativeInfinity;
            double B = double.IsPositiveInfinity(b) ? double.PositiveInfinity
                                                    : (b > 0.0 ? Math.Log(b) : double.NegativeInfinity);

            // Standardized limits for Y
            double alpha = (A - mu) / sigma;   // can be -inf
            double beta = (B - mu) / sigma;   // can be +inf

            // Standard normal CDF helper (consistent with your Normal code)
            static double Phi(double x) => 0.5 * (1.0 + Mathematics.SpecialFunctions.Erf.Function(x / Math.Sqrt(2.0)));

            double PhiA = double.IsNegativeInfinity(alpha) ? 0.0 : Phi(alpha);
            double PhiB = double.IsPositiveInfinity(beta) ? 1.0 : Phi(beta);

            double Z = PhiB - PhiA;    // normalization (P(a < X < b))
            if (Z <= 1e-15)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Raw truncated moments: E[X^k | a < X < b] for k = 1..4
            double[] Eraw = new double[5]; // we'll fill 1..4

            for (int k = 1; k <= 4; k++)
            {
                double ks = k * sigma;
                // shifted limits: Φ(β - kσ) − Φ(α - kσ)
                double PhiBk = double.IsPositiveInfinity(beta) ? 1.0 : Phi(beta - ks);
                double PhiAk = double.IsNegativeInfinity(alpha) ? 0.0 : Phi(alpha - ks);

                double scale = Math.Exp(k * mu + 0.5 * k * k * sigma * sigma);
                Eraw[k] = scale * (PhiBk - PhiAk) / Z;
            }

            // Unconditional mean of X (about which we take central moments)
            double muX = Math.Exp(mu + 0.5 * sigma * sigma);

            // Central moments about μ_X (unconditional)
            double m1 = Eraw[1];
            double m2 = Eraw[2] - 2.0 * muX * Eraw[1] + muX * muX;
            double m3 = Eraw[3] - 3.0 * muX * Eraw[2] + 3.0 * muX * muX * Eraw[1] - muX * muX * muX;
            double m4 = Eraw[4]
                      - 4.0 * muX * Eraw[3]
                      + 6.0 * muX * muX * Eraw[2]
                      - 4.0 * muX * muX * muX * Eraw[1]
                      + muX * muX * muX * muX;

            return new[] { m1, m2, m3, m4 };
        }

    }
}
