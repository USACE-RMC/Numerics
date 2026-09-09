using System;
using System.Collections.Generic;
using System.Linq;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Weibull probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Weibull_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Weibull : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IStandardError, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a Weibull distribution with scale = 10 and shape = 2.
        /// </summary>
        public Weibull()
        {
            SetParameters(10d, 2d);
        }

        /// <summary>
        /// Constructs a Weibull distribution with the given parameters λ and k.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        public Weibull(double scale, double shape)
        {
            SetParameters(scale, shape);
        }

        private double _lambda;
        private double _kappa;

        /// <summary>
        /// Gets and sets the scale parameter λ (lambda).
        /// </summary>
        public double Lambda
        {
            get { return _lambda; }
            set
            {
                _parametersValid = ValidateParameters(value, Kappa, false) is null;
                _lambda = value;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter κ (kappa).
        /// </summary>
        public double Kappa
        {
            get { return _kappa; }
            set
            {
                _parametersValid = ValidateParameters(Lambda, value, false) is null;
                _kappa = value;
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
            get { return UnivariateDistributionType.Weibull; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Weibull"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "W"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Scale (λ)";
                parmString[1, 0] = "Shape (κ)";
                parmString[0, 1] = Lambda.ToString();
                parmString[1, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["λ", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Lambda), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Lambda, Kappa]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (Kappa == 1) return Lambda;
                return Math.Exp(Math.Log(Lambda) + Gamma.LogGamma(1 + 1 / Kappa));
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(.5); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Kappa <= 1.0d)
                {
                    return 0.0d;
                }
                else
                {
                    return Lambda * Math.Pow((Kappa - 1.0d) / Kappa, 1.0d / Kappa);
                }
            }
        }

        /// <inheritdoc/>
        /// <remarks>Uses the exact exponential-power relationship to GEV and restores scale in
        /// logarithms, without forming lambda squared or overflowing raw Gamma moments.</remarks>
        public override double StandardDeviation
        {
            get
            {
                if (Kappa == 1) return Lambda;
                double power = 1 / Kappa;
                if (double.IsPositiveInfinity(power)) return double.PositiveInfinity;
                double logarithm = power <= .05
                    ? Math.Log(power) + Math.Log(new GeneralizedExtremeValue(0, 1, power).StandardDeviation)
                    : .5 * GeneralizedExtremeValue.LogPowerVariance(power);
                return Math.Exp(Math.Log(Lambda) + logarithm);
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                double power = 1 / Kappa;
                // For T unit exponential, T^power = 1 - power*GEV(0,1,power).
                return double.IsPositiveInfinity(power) ? double.PositiveInfinity
                    : -new GeneralizedExtremeValue(0, 1, power).Skewness;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double power = 1 / Kappa;
                return double.IsPositiveInfinity(power) ? double.PositiveInfinity
                    : new GeneralizedExtremeValue(0, 1, power).Kurtosis;
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
            var newDistribution = new Weibull(Lambda, Kappa);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        public void SetParameters(double scale, double shape)
        {
            Lambda = scale;
            Kappa = shape;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(parameters), "Exactly two parameters are required.");
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double scale, double shape, bool throwException)
        {
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Lambda), "The scale parameter λ (lambda) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Lambda), "The scale parameter λ (lambda) must be positive.");
            }
            if (double.IsNaN(shape) || double.IsInfinity(shape) || shape <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "Exactly two parameters are required.");
                if (throwException) throw exception;
                return exception;
            }
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        /// <remarks>Requires finite, nonconstant observations. Preserves the legacy initializer's treatment
        /// of nonpositive observations whenever its initial values and rounded prior bounds are usable.
        /// The exceptional-input fallback requires positive observations and retains representable small scales.
        /// Density support is unchanged.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample or a finite feasible initialization is invalid.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2);
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
            initialVals = LegacyConstraintSolveMLE(sample);
            // Get bounds of scale
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of shape
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Retains the established constraint initializer arithmetic for ordinary samples.</summary>
        /// <param name="samples">Observations.</param>
        /// <returns>The legacy initial parameter values.</returns>
        /// <exception cref="Exception">Fewer than two observations are supplied.</exception>
        private double[] LegacyConstraintSolveMLE(IList<double> samples)
        {
            double n = samples.Count;
            if (n <= 1d)
            {
                throw new Exception("Observations not sufficient. There must be more than 1 data point.");
            }

            double s1 = 0d;
            double s2 = 0d;
            double s3 = 0d;
            double previousC = int.MinValue;
            double QofC = 0d;
            double c = 10d; // shape
            double b = 0d; // scale

            // solve for the shape parameter
            while (Math.Abs(c - previousC) >= 0.0001d)
            {
                s1 = 0d;
                s2 = 0d;
                s3 = 0d;
                foreach (double x in samples)
                {
                    if (x > 0d)
                    {
                        s1 += Math.Log(x);
                        s2 += Math.Pow(x, c);
                        s3 += Math.Pow(x, c) * Math.Log(x);
                    }
                }

                QofC = n * s2 / (n * s3 - s1 * s2);
                previousC = c;
                c = (c + QofC) / 2d;
            }

            // solve for scale
            foreach (double x in samples)
            {
                if (x > 0d)
                {
                    b += Math.Pow(x, c);
                }
            }

            b = Math.Pow(b / n, 1d / c);

            // return parameters
            return [b, c];
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The sample is invalid, insufficient, constant, or produces an invalid positive parameter.</exception>
        /// <exception cref="InvalidOperationException">The initialization iteration is numerically unresolved.</exception>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2, true);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            var initialVals = SolveMLE(sample);
            DistributionNumerics.PositiveParameterBounds(initialVals[0], out lowerVals[0], out upperVals[0]);
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
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
                var W = new Weibull();
                W.SetParameters(x);
                return W.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <summary>
        /// The Maximum Likelihood Estimation method for the Weibull distribution.
        /// </summary>
        /// <param name="samples">The array of sample data.</param>
        /// <remarks>
        /// Implemented according to: Parameter estimation of the Weibull probability distribution, 1994, Hongzhu Qiao, Chris P. Tsokos
        /// <para>All observations must be finite and strictly positive, with at least two distinct
        /// observations. The same fixed-point iteration is evaluated with log-relative bounded
        /// weights, retaining every observation and the existing convergence threshold.</para>
        /// <para>
        /// References:
        /// This code was copied and modified from the Math.NET Library.
        /// <list type="bullet">
        /// <item><description>
        /// Math.NET Numerics Library, http://numerics.mathdotnet.com
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <returns>The initial scale and shape estimates.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The observations are invalid, insufficient or constant.</exception>
        /// <exception cref="InvalidOperationException">The initialization iteration is numerically unresolved.</exception>
        public double[] SolveMLE(IList<double> samples)
        {
            DistributionNumerics.ValidateSample(samples, 2, true);
            double n = samples.Count;
            double scale = DistributionNumerics.InitializationScale(samples);
            var logRatios = new double[samples.Count];
            for (int i = 0; i < samples.Count; i++)
            {
                double ratio = samples[i] / scale;
                logRatios[i] = ratio > 0 ? Math.Log(ratio) : Math.Log(samples[i]) - Math.Log(scale);
            }

            double s1 = 0d;
            double s2 = 0d;
            double s3 = 0d;
            double previousC = int.MinValue;
            double QofC = 0d;
            double c = 10d; // shape
            double b = 0d; // scale

            // solve for the shape parameter
            while (Math.Abs(c - previousC) >= 0.0001d)
            {
                s1 = 0d;
                s2 = 0d;
                s3 = 0d;
                foreach (double logarithm in logRatios)
                {
                    double weight = Math.Exp(c * logarithm);
                    s1 += logarithm;
                    s2 += weight;
                    s3 += weight * logarithm;
                }

                QofC = n * s2 / (n * s3 - s1 * s2);
                previousC = c;
                c = (c + QofC) / 2d;
                if (!(c > 0) || !Tools.IsFinite(c))
                    throw new InvalidOperationException("The Weibull initialization iteration did not produce a finite positive shape.");
            }

            // solve for scale
            foreach (double logarithm in logRatios)
            {
                b += Math.Exp(c * logarithm);
            }

            b = scale * Math.Pow(b / n, 1d / c);

            // return parameters
            return [b, c];
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
        /// When <c>x = 0</c> and the shape <c>κ &lt; 1</c>, the Weibull density has a genuine
        /// integrable singularity and this method intentionally returns positive infinity.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Lambda, Kappa, true);
            if (x < Minimum || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            if (x == 0) return Kappa == 1 ? -Math.Log(Lambda) : Kappa < 1 ? double.PositiveInfinity : double.NegativeInfinity;
            double logarithm = LogStandardizedValue(x);
            double power = Math.Exp(Kappa * logarithm);
            if (double.IsPositiveInfinity(power)) return double.NegativeInfinity;
            double lf = Math.Log(Kappa) - Math.Log(Lambda) + (Kappa - 1) * logarithm - power;
            return double.IsNaN(lf) ? double.NegativeInfinity : lf;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Lambda, Kappa, true);
            if (x <= Minimum)
                return 0d;
            return -Tools.Expm1(-Math.Exp(Kappa * LogStandardizedValue(x)));
        }

        /// <inheritdoc/>
        /// <remarks>Retains the lower-tail logarithm when the positive power underflows.</remarks>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Lambda, Kappa, true);
            if (x <= 0) return double.NegativeInfinity;
            double logarithm = Kappa * LogStandardizedValue(x);
            double power = Math.Exp(logarithm);
            return power == 0 ? logarithm : DistributionNumerics.Log1mExp(-power);
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Lambda, Kappa, true);
            return x <= 0 ? 0 : -Math.Exp(Kappa * LogStandardizedValue(x));
        }

        /// <summary>Forms log(x/lambda) without an overflowing or underflowing intermediate ratio.</summary>
        /// <param name="x">The positive observation in physical coordinates.</param>
        /// <returns>The logarithm of <c>x/lambda</c>.</returns>
        private double LogStandardizedValue(double x)
        {
            double ratio = x / Lambda;
            return ratio > 0 && !double.IsInfinity(ratio) ? Math.Log(ratio) : Math.Log(x) - Math.Log(Lambda);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Lambda, Kappa, true);
            // Compute the inverse CDF
            // The exact exponential identity preserves subnormal probabilities on .NET Framework.
            if (Kappa == 1) return -Lambda * Tools.Log1p(-probability);
            double logarithm = Math.Log(-Tools.Log1p(-probability)) / Kappa;
            double unitQuantile = Math.Exp(logarithm);
            return unitQuantile < 1E-200 || double.IsInfinity(unitQuantile)
                ? Math.Exp(Math.Log(Lambda) + logarithm) : Lambda * unitQuantile;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Weibull(Lambda, Kappa);
        }

        /// <inheritdoc/>
        /// <remarks>Retains the published rounded MLE covariance constants in scale and shape coordinates.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or sample size is not positive.</exception>
        /// <exception cref="NotImplementedException">The method is not maximum likelihood.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_lambda, _kappa, true);

            // Compute covariance
            double a = Lambda;
            double b = Kappa;
            var covar = new double[2, 2];
            double scaledA = (a / b) / Math.Sqrt(sampleSize);
            double scaledB = b / Math.Sqrt(sampleSize);
            covar[0, 0] = 1.108665d * (scaledA * scaledA); // scale
            covar[1, 1] = 0.607927d * (scaledB * scaledB); // shape
            covar[0, 1] = 0.257022d * a / sampleSize;
            covar[1, 0] = covar[0, 1];
            return covar;
        }

        /// <inheritdoc/>
        /// <remarks>Uses the actual inverse-CDF gradient and restores scale after the normalized covariance contraction.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Lambda, Kappa, true);
            var covar = new Weibull(1, Kappa).ParameterCovariance(sampleSize, estimationMethod);
            var grad = QuantileGradient(probability);
            return DistributionNumerics.ScaledQuantileVariance(covar, [InverseCDF(probability), grad[1]]);
        }

        /// <inheritdoc/>
        /// <remarks>For t=-log(1-p), returns [Q/lambda, -Q*log(t)/kappa squared].</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_lambda, _kappa, true);
            double logT = Math.Log(-Tools.Log1p(-probability));
            double quantile = InverseCDF(probability);
            var gradient = new double[]
            {
                Math.Exp(logT / Kappa), // scale
                -(quantile / Kappa) * (logT / Kappa) // shape
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
