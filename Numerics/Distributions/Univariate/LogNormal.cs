using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Log-Normal probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class LogNormal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Log-Normal distribution with a mean (of log) of 3 and standard deviation (of log) of 0.5
        /// </summary>
        public LogNormal()
        {
            SetParameters(3d, 0.5d);
        }

        /// <summary>
        /// Constructs a Log-Normal distribution with given mean (of log) and standard deviation (of log).
        /// </summary>
        /// <param name="meanOfLog">The mean of the log transformed data.</param>
        /// <param name="standardDeviationOfLog">The standard deviation of the log transformed data.</param>
        public LogNormal(double meanOfLog, double standardDeviationOfLog)
        {
            SetParameters(meanOfLog, standardDeviationOfLog);
        }

        // Private variables
        private double _mu;
        private double _sigma;
        private double _base = 10d;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, false) is null;
                _mu = value;
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
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, false) is null;
                _sigma = value;
            }
        }

        /// <summary>
        /// Gets and sets the finite logarithm base, which must be greater than one.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">The value is nonfinite or no greater than one.</exception>
        public double Base
        {
            get { return _base; }
            set
            {
                if (!(value > 1d) || double.IsInfinity(value))
                    throw new ArgumentOutOfRangeException(nameof(Base), "The logarithm base must be finite and greater than one.");
                _base = value;
            }
        }

        /// <summary>
        /// Gets the log correction factor.
        /// </summary>
        private double K
        {
            get { return 1d / Math.Log(Base); }
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Mu, Sigma, true);
            return x <= 0d ? double.NegativeInfinity : DistributionNumerics.NormalLogCDF(DistributionNumerics.Standardize(Math.Log(x, Base), Mu, Sigma));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Mu, Sigma, true);
            return x <= 0d ? 0d : DistributionNumerics.NormalLogSurvival(DistributionNumerics.Standardize(Math.Log(x, Base), Mu, Sigma));
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.LogNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Log-Normal (base 10)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LogN"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean (of log) (µ)";
                parmString[1, 0] = "Std Dev (of log) (σ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
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
            get { return [nameof(Mu), nameof(Sigma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                double lnB = Math.Log(Base);
                return Math.Exp((Mu + 0.5 * Sigma * Sigma * lnB) * lnB);
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <summary>
        /// Gets the mode of the distribution.
        /// </summary>
        public override double Mode
        {
            get { double logBase = Math.Log(Base); return Math.Exp(Mu * logBase - Math.Pow(Sigma * logBase, 2d)); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                double logBase = Math.Log(Base), variance = Math.Pow(Sigma * logBase, 2d);
                double logExcess = variance > 0.5d ? variance + Tools.Log1p(-Math.Exp(-variance)) : Math.Log(Tools.Expm1(variance));
                return Math.Exp(Mu * logBase + 0.5d * variance + 0.5d * logExcess);
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                double variance = Math.Pow(Sigma * Math.Log(Base), 2d);
                return (Math.Exp(variance) + 2d) * Math.Sqrt(Tools.Expm1(variance));
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double variance = Math.Pow(Sigma * Math.Log(Base), 2d);
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
            // The mean of the base-log observations is a location parameter and can be any
            // finite value; only the log-space standard deviation is bounded below by zero.
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
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(IndirectMethodOfMoments(sample));
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(IndirectMethodOfLinearMoments(sample)));
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
        /// <remarks>The bootstrap distribution retains the configured <see cref="Base"/>.</remarks>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new LogNormal(Mu, Sigma) { Base = Base };
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters based on the moments of the log transformed data.
        /// </summary>
        /// <param name="meanOfLog">The mean of the log transformed data.</param>
        /// <param name="standardDeviationOfLog">The standard deviation of the log transformed data.</param>
        public void SetParameters(double meanOfLog, double standardDeviationOfLog)
        {
            // Set parameters
            Mu = meanOfLog;
            Sigma = standardDeviationOfLog;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mu">The mean (of log).</param>
        /// <param name="sigma">The standard deviation (of log).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double mu, double sigma, bool throwException)
        {
            if (double.IsNaN(mu) || double.IsInfinity(mu))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
            }
            if (double.IsNaN(sigma) || double.IsInfinity(sigma) || sigma <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Sigma must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Sigma must be positive.");
            }
            return null!;
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
        public double[] IndirectMethodOfMoments(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformedSample[i] = Math.Log(sample[i], Base);
            return Statistics.ProductMoments(transformedSample);
        }

        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        public double[] IndirectMethodOfLinearMoments(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformedSample[i] = Math.Log(sample[i], Base);
            return Statistics.LinearMoments(transformedSample);
        }

        /// <summary>
        /// Sets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="mean">The real-space mean of the data.</param>
        /// <param name="standardDeviation">The real-space standard deviation of the data.</param>
        public double[] DirectMethodOfMoments(double mean, double standardDeviation)
        {
            double[] natural = LnNormal.DirectMethodOfMoments(mean, standardDeviation);
            double logBase = Math.Log(Base);
            return [natural[0] / logBase, natural[1] / logBase];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            return DirectMethodOfMoments(moments[0], moments[1]);
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new LogNormal() { Base = Base };
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
        /// <remarks>Preserves the legacy 0.1 substitution for nonpositive observations when constructing
        /// initial values and rounded prior bounds. This does not modify the sample or density support.</remarks>
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
            // Estimate initial values using the method of moments (a.k.a product moments).
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++)
                transformedSample[i] = Math.Log(sample[i] > 0d ? sample[i] : 0.1d, Base);
            var mom = Statistics.ProductMoments(transformedSample);
            initialVals = new double[] { mom[0], mom[1] };
            // Get bounds of mean. The mean is a location parameter on the log scale and is
            // legitimately negative whenever the data are mostly below 1, so the bounds are
            // symmetric about zero from the magnitude of the initial value, matching Normal's
            // location bounds. A machine-epsilon floor here would reject any sub-unity sample
            // before a fit could start.
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of standard deviation
            double real = Math.Exp(initialVals[1] / K);
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Ceiling(Math.Log(Math.Pow(10d, Math.Ceiling(Math.Log10(real) + 1d)), Base));
            upperVals[1] = double.IsNaN(upperVals[1]) ? 4 : upperVals[1];
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformed = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformed[i] = Math.Log(sample[i], Base);
            return new Normal().GetRobustParameterConstraints(transformed);
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
                var LogN = new LogNormal() { Base = Base };
                LogN.SetParameters(x);
                return LogN.LogLikelihood(sample);
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
            if (!_parametersValid) ValidateParameters(Mu, Sigma, true);
            if (x <= 0d || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            double logX = Math.Log(x), logBase = Math.Log(Base);
            double z = DistributionNumerics.Standardize(logX / logBase, Mu, Sigma);
            return -0.5d * z * z - Math.Log(Sigma) - Math.Log(Tools.Sqrt2PI) - Math.Log(logBase) - logX;
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
                ValidateParameters(Mu, Sigma, true);
            return Math.Exp((Mu + Sigma * Normal.StandardZ(probability)) / K);
        }

        /// <summary>
        /// Get confidence intervals using Monte Carlo simulation.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="realizations">The number of Monte Carlo realizations.</param>
        /// <param name="quantiles">List of exceedance probabilities for output frequency curves.</param>
        /// <param name="percentiles">List of confidence percentiles for confidence interval output.</param>
        /// <remarks>
        /// This is the same sampling approach as used in HEC-FDA.
        /// Each simulated distribution retains the configured <see cref="Base"/>.
        /// </remarks>
        public double[,] MonteCarloConfidenceIntervals(int sampleSize, int realizations, IList<double> quantiles, IList<double> percentiles)
        {
            DistributionNumerics.ValidateConfidenceInputs(sampleSize, quantiles, percentiles, 2);
            if (realizations <= 0) throw new ArgumentOutOfRangeException(nameof(realizations), "At least one realization is required.");
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Dimension output array
            int q = quantiles.Count;
            int p = percentiles.Count;
            var Output = new double[q, p];

            // Variables
            double OriginalMean = Mu;
            double OriginalStdDev = Sigma;

            // Create random numbers for mean and standard deviation
            var r = new MersenneTwister(12345);
            var rndMean = r.NextDoubles(realizations);
            var rndStdDev = r.NextDoubles(realizations);


            // Create list of Monte Carlo distributions
            var MonteCarloDistributions = new UnivariateDistributionBase[realizations];
            Parallel.For(0, realizations, idx =>
            {

                // Generate new mean
                var Normal = new Normal(OriginalMean, OriginalStdDev / Math.Sqrt(sampleSize));
                double NewMu = Normal.InverseCDF(rndMean[idx]);
                // Generate new standard deviation
                var Chi = new ChiSquared(sampleSize - 1);
                double NewSigma = Math.Sqrt((sampleSize - 1) * Math.Pow(OriginalStdDev, 2d) / Chi.InverseCDF(rndStdDev[idx]));
                // Create a new distribution with the new parameters
                MonteCarloDistributions[idx] = new LogNormal(NewMu, NewSigma) { Base = Base };
            });

            // Create confidence intervals
            for (int i = 0; i < q; i++)
            {
                // Create array of X values across user-defined probabilities
                var XValues = new double[realizations];
                // Record X values
                Parallel.For(0, realizations, idx => XValues[idx] = MonteCarloDistributions[idx].InverseCDF(quantiles[i]));
                // Record percentiles for user-defined probabilities
                for (int j = 0; j < p;  j++)
                    Output[i, j] = Statistics.Percentile(XValues, percentiles[j]);
            }

            // Return confidence percentile output
            return Output;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new LogNormal(Mu, Sigma) { Base = Base };
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
                ValidateParameters(Mu, _sigma, true);
            // Compute covariance in (μ, σ) parameterization of log-space.
            // Var(μ̂) = σ²/n, Var(σ̂) = σ²/(2n), Cov = 0.
            // Log-moment and maximum-likelihood estimators share this leading asymptotic covariance.
            double scaled = Sigma / Math.Sqrt(sampleSize);
            double s2 = scaled * scaled;
            var covar = new double[2, 2];
            covar[0, 0] = s2; // Var(μ̂)
            covar[1, 1] = s2 / 2d; // Var(σ̂)
            covar[0, 1] = 0.0;
            covar[1, 0] = covar[0, 1];
            return covar;
        }

        /// <inheritdoc/>
        /// <remarks>Combines the physical quantile and parameter standard errors in log space before squaring, so representable variance is retained when the quantile square would overflow.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments && estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException();
            if (!_parametersValid) ValidateParameters(Mu, Sigma, true);
            double z = Normal.StandardZ(probability), logBase = Math.Log(Base);
            double logQuantile = (Mu + Sigma * z) * logBase;
            double logStandardError = logQuantile + Math.Log(logBase) + Math.Log(Sigma) - 0.5d * Math.Log(sampleSize);
            // For the mean and SD of base-log observations, covariance is zero and the SD variance is half the mean variance.
            return Math.Exp(2d * logStandardError + Tools.Log1p(0.5d * z * z));
        }

        /// <inheritdoc/>
        /// <remarks>Returns derivatives of the physical quantile with respect to the configured base-log <see cref="Mu"/> and <see cref="Sigma"/>. The change-of-base and exponential Jacobian are applied exactly once.</remarks>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Mu, Sigma, true);
            double z = Normal.StandardZ(probability);
            double factor = InverseCDF(probability) * Math.Log(Base);
            return [factor, factor * z];
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }


    }
}
