using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Rayleigh probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References:</b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Rayleigh_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Rayleigh : UnivariateDistributionBase, IEstimation, IBootstrappable
    {
       
        /// <summary>
        /// Constructs a Rayleigh distribution with a scale (σ) = 10.
        /// </summary>
        public Rayleigh()
        {
            SetParameters(10d);
        }

        /// <summary>
        /// Constructs a Rayleigh distribution with a specified scale (σ).
        /// </summary>
        public Rayleigh(double scale)
        {
            SetParameters(scale);
        }
   
        private double _sigma;

        /// <summary>
        /// Gets and sets the scale parameter σ (sigma).
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters(value, false) is null;
                _sigma = value;
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
            get { return UnivariateDistributionType.Rayleigh; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Rayleigh"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "RAY"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                parmString[0, 0] = "Scale (σ)";
                parmString[0, 1] = Sigma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Sigma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Sigma * Math.Sqrt(Math.PI / 2d); }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Sigma * Math.Sqrt(Math.Log(4.0d)); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Sigma; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt((4.0d - Math.PI) / 2.0d * Sigma * Sigma); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 2.0d * Math.Sqrt(Math.PI) * (Math.PI - 3.0d) / Math.Pow(4.0d - Math.PI, 1.5d); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double num = -(6d * Math.Pow(Math.PI, 2d) - 24d * Math.PI + 16d);
                double den = Math.Pow(4d - Math.PI, 2d);
                return 3d + num / den;
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
            get { return [0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The method of moments sets σ = mean(x) / √(π/2).
        /// </para>
        /// <para>
        /// Maximum likelihood first forms the closed-form estimator σ̂ = √(Σxᵢ² / (2n)) and then multiplies it by
        /// the bias correction Γ(n)·√n / Γ(n + 1/2), which makes the estimate of σ unbiased for a sample of size n
        /// (see the Rayleigh distribution reference above). The corrected value is therefore slightly larger than the
        /// plain maximum likelihood estimate returned by, for example, scipy.stats.rayleigh.fit with the location
        /// fixed at zero; the two agree as n grows.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a nonfinite or negative value,
        /// or has no positive value.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not supported.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            ValidateSample(sample);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(Statistics.Mean(sample) / Math.Sqrt(Math.PI / 2.0d));
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                // https://en.wikipedia.org/wiki/Rayleigh_distribution
                // compute the biased estimator first
                double sum = 0d;
                for (int i = 0; i < sample.Count; i++)
                    sum += Math.Pow(sample[i], 2d);
                double n = sample.Count;
                double c = 1d / (2d * n);
                double biased = Math.Sqrt(c * sum);
                // then compute the biased corrected estimator
                double num = Gamma.LogGamma(n) + 0.5d * Math.Log(n);
                double den = Gamma.LogGamma(n + 0.5d);
                double correction = Math.Exp(num - den);
                if (double.IsPositiveInfinity(num) && double.IsPositiveInfinity(den))
                {
                    correction = 1d;
                }

                SetParameters(biased * correction);
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Checks that a sample can be used to estimate the scale parameter.
        /// </summary>
        /// <param name="sample">The observations to validate.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a nonfinite or negative value,
        /// or has no positive value.</exception>
        /// <remarks>
        /// The support is [0, ∞). A single observation, or a constant sample, still determines σ, so neither is
        /// rejected; an all-zero sample would give σ = 0 and is rejected.
        /// </remarks>
        private static void ValidateSample(IList<double> sample)
        {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new ArgumentOutOfRangeException(nameof(sample), "Insufficient observations to initialize the distribution.");
            bool positive = false;
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]) || sample[i] < 0d)
                    throw new ArgumentOutOfRangeException(nameof(sample), "Observations must be finite and nonnegative.");
                positive |= sample[i] > 0d;
            }
            if (!positive) throw new ArgumentOutOfRangeException(nameof(sample), "At least one observation must be positive.");
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new Rayleigh(Sigma);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        public void SetParameters(double scale)
        {
            Sigma = scale;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double scale, bool throwException)
        {
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be greater than zero.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be greater than zero.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], throwException);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated through z = x/σ as (z/σ)·exp(−z²/2), so σ² is never formed and cannot overflow.
        /// </remarks>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Sigma, true);
            if (x < Minimum || double.IsPositiveInfinity(x)) return 0.0d;
            double z = x / Sigma;
            return z / Sigma * Math.Exp(-0.5d * z * z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as ln(x) − 2·ln(σ) − (x/σ)²/2, which stays finite after the density underflows.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Sigma, true);
            if (x < Minimum || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            double z = x / Sigma;
            return Math.Log(x) - 2.0d * Math.Log(Sigma) - 0.5d * z * z;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as −expm1(−(x/σ)²/2), which keeps full relative precision for small x.
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Sigma, true);
            if (x <= Minimum)
                return 0d;
            double z = x / Sigma;
            return -Tools.Expm1(-0.5d * z * z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated directly as exp(−(x/σ)²/2), so upper-tail probabilities are not lost to the subtraction 1 − CDF.
        /// </remarks>
        public override double CCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Sigma, true);
            if (x <= Minimum)
                return 1d;
            double z = x / Sigma;
            return Math.Exp(-0.5d * z * z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as −(x/σ)²/2, which stays finite after the survival probability underflows.
        /// </remarks>
        public override double LogCCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Sigma, true);
            if (x <= Minimum)
                return 0d;
            double z = x / Sigma;
            return -0.5d * z * z;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as σ·√(−2·log1p(−p)), which keeps full relative precision for small p.
        /// </remarks>
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
                ValidateParameters(Sigma, true);
            return Sigma * Math.Sqrt(-2d * Tools.Log1p(-probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Rayleigh(Sigma);
        }

    }
}