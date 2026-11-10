using System;
using System.Buffers.Text;
using System.Collections.Generic;
using System.Linq;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The log-Pearson Type III distribution.
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
    /// <see href = "http://mathworld.wolfram.com/PearsonTypeIIIDistribution.html" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class LogPearsonTypeIII : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a log-Pearson Type III distribution with a mean (of log) of 3, std dev (of log) of 0.5, and skew (of log) of 0.
        /// </summary>
        public LogPearsonTypeIII()
        {
            SetParameters(3d, 0.5d, 0d);
        }

        /// <summary>
        /// Constructs a log-Pearson Type III distribution with the given moments (of log) µ, σ, and γ.
        /// </summary>
        /// <param name="meanOfLog">The mean of the log transformed data.</param>
        /// <param name="standardDeviationOfLog">The standard deviation of the log transformed data.</param>
        /// <param name="skewOfLog">The skew of the log transformed data.</param>
        public LogPearsonTypeIII(double meanOfLog, double standardDeviationOfLog, double skewOfLog)
        {
            SetParameters(meanOfLog, standardDeviationOfLog, skewOfLog);
        }

        private double _mu;
        private double _sigma;
        private double _gamma;
        private double _base = 10d;
        private double _logBase = Math.Log(10d); // Math.Log(Base)
        private double _logLogBase = Math.Log(Math.Log(10d)); // Math.Log(Math.Log(Base))
        private double _boundary = double.NaN; // Math.Exp(Xi * Math.Log(Base)), the finite support endpoint when Gamma is nonzero

        /// <summary>
        /// Gets and sets the Mean (of log) of the distribution.
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, Gamma, false) is null;
                SetParameters(value, Sigma, Gamma);
            }
        }

        /// <summary>
        /// Gets and sets the Standard Deviation (of log) of the distribution.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, Gamma, false) is null;
                SetParameters(Mu, value, Gamma);
            }
        }

        /// <summary>
        /// Gets and sets the Skew (of log) of the distribution.
        /// </summary>
        public double Gamma
        {
            get { return _gamma; }
            set
            {
                _parametersValid = ValidateParameters(Mu, Sigma, value, false) is null;
                SetParameters(Mu, Sigma, value);
            }
        }

        /// <summary>
        /// Gets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return Mu - Sigma * (2d / Gamma); }
        }

        /// <summary>
        /// Gets the scale parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return 0.5d * Sigma * Gamma; }
        }

        /// <summary>
        /// Gets the shape parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return Math.Pow(2d / Gamma, 2d); }
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
                UpdateLogBaseCache();
            }
        }

        /// <summary>
        /// Recomputes the cached logarithms of the base and the support endpoint from the current parameters and base.
        /// </summary>
        private void UpdateLogBaseCache()
        {
            _logBase = Math.Log(_base);
            _logLogBase = Math.Log(_logBase);
            _boundary = Math.Exp(Xi * _logBase);
        }

        /// <summary>
        /// Gets the log correction factor.
        /// </summary>
        private double K
        {
            get { return 1d / _logBase; }
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0d;
            return PearsonTypeIII.LogTail(Mu, Sigma, Gamma, Math.Log(x, Base), false);
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return double.NegativeInfinity;
            return PearsonTypeIII.LogTail(Mu, Sigma, Gamma, Math.Log(x, Base), true);
        }

        /// <summary>Log raw-moment contribution after removing the location, with each moment's existence checked separately.</summary>
        /// <param name="order">The positive raw-moment order.</param>
        /// <returns>The logarithm of the shape-dependent raw-moment factor, or positive infinity when that moment does not exist.</returns>
        private double LogMomentShape(int order)
        {
            double scale = order * Sigma * Math.Log(Base), argument = Gamma * scale / 2d;
            if (!(argument < 1d)) return double.PositiveInfinity;
            if (Math.Abs(argument) < 0.01d)
            {
                // [-log(1-u)-u]/u^2 = sum u^k/(k+2), analytic at u=0.
                double sum = 0.5d, term = 1d;
                for (int k = 1; k <= 12; k++) { term *= argument; sum += term / (k + 2d); }
                return scale * scale * sum;
            }
            // -ln(1-u)-u from the cancellation-free log(1+x)-x.
            return Alpha * -Tools.Log1pMinusX(-argument);
        }

        /// <summary>
        /// Value of <c>σ·ln(Base)·(1 + |γ|)</c> below which the skewness and kurtosis are evaluated in cancellation-free form.
        /// </summary>
        /// <remarks>
        /// The closed-form exponential moments lose about <c>ε/x²</c> relative to cancellation as the log scale
        /// <c>x = σ·ln(Base)·(1 + |γ|)</c> decreases. Below this value <see cref="SmallScaleStandardizedMoments"/> removes the
        /// cancelling terms analytically. On a reference grid of log scales from 1E-5 to 0.5 and skews up to 6 in magnitude,
        /// both forms agree with 60-digit references to about 1.1E-12 relative in the skewness and 5E-13 in the kurtosis.
        /// </remarks>
        private const double SmallScaleLimit = 0.15d;

        /// <summary>Evaluates the standardized third and fourth central moments without cancellation for a small log scale.</summary>
        /// <param name="skewness">The resulting standardized third central moment.</param>
        /// <param name="kurtosis">The resulting standardized fourth central moment.</param>
        /// <remarks>
        /// With <c>Kᵣ</c> the location-free log moment of order <c>r</c> and <c>dᵣ = Kᵣ - rK₁</c>, the normalized central
        /// moments are <c>e^d₃ - 3e^d₂ + 2</c> and <c>e^d₄ - 4e^d₃ + 6e^d₂ - 3</c>, whose leading terms cancel for a small log
        /// scale <c>t = σ·ln(Base)</c>. They are formed as <c>F(d₃) - 3F(d₂) + Δ₃</c> and <c>F(d₄) - 4F(d₃) + 6F(d₂) + Δ₄</c>,
        /// with <c>F(d) = e^d - 1 - d</c> and the finite differences <c>Δ₃ = K₃ - 3K₂ + 3K₁</c> and
        /// <c>Δ₄ = K₄ - 4K₃ + 6K₂ - 4K₁</c>. Because <c>Kᵣ = α[-ln(1 - ru) - ru]</c> with <c>u = tγ/2</c> and <c>αu² = t²</c>,
        /// each difference is the power series <c>t² Σₖ (cₖ/k) u^(k-2)</c> whose coefficients <c>cₖ</c> are the same finite
        /// differences of <c>rᵏ</c>; these vanish below orders 3 and 4, so the series carry no cancellation and are analytic
        /// at zero skew. The <c>dᵣ</c> are summed the same way, from the coefficients <c>rᵏ - r</c>. For <c>4|u| ≥ 0.5</c> the
        /// differences are formed directly from the log moments.
        /// </remarks>
        private void SmallScaleStandardizedMoments(out double skewness, out double kurtosis)
        {
            double scale = Sigma * Math.Log(Base), u = Gamma * scale / 2d, scale2 = scale * scale;
            double d2, d3, d4, delta3, delta4;
            if (4d * Math.Abs(u) < 0.5d)
            {
                // dᵣ = t² Σₖ (rᵏ - r)/k u^(k-2) and the differences Δ₃, Δ₄ share the powers of u.
                double sum2 = 0d, sum3 = 0d, sum4 = 0d, diff3 = 0d, diff4 = 0d, power = 1d, two = 4d, three = 9d, four = 16d;
                for (int k = 2; k < 200; k++)
                {
                    double t2 = (two - 2d) / k * power;
                    double t3 = (three - 3d) / k * power;
                    double t4 = (four - 4d) / k * power;
                    double s3 = (three - 3d * two + 3d) / k * power;
                    double s4 = (four - 4d * three + 6d * two - 4d) / k * power;
                    sum2 += t2; sum3 += t3; sum4 += t4; diff3 += s3; diff4 += s4;
                    if (k >= 4 && Math.Abs(t4) <= 1E-17 * Math.Abs(sum4) && Math.Abs(s3) <= 1E-17 * Math.Abs(diff3)
                        && Math.Abs(s4) <= 1E-17 * Math.Abs(diff4)) break;
                    power *= u; two *= 2d; three *= 3d; four *= 4d;
                }
                d2 = scale2 * sum2;
                d3 = scale2 * sum3;
                d4 = scale2 * sum4;
                delta3 = scale2 * diff3;
                delta4 = scale2 * diff4;
            }
            else
            {
                double k1 = LogMomentShape(1), k2 = LogMomentShape(2), k3 = LogMomentShape(3), k4 = LogMomentShape(4);
                d2 = k2 - 2d * k1;
                d3 = k3 - 3d * k1;
                d4 = k4 - 4d * k1;
                delta3 = k3 - 3d * k2 + 3d * k1;
                delta4 = k4 - 4d * k3 + 6d * k2 - 4d * k1;
            }
            double variance = Tools.Expm1(d2);
            double f2 = Tools.Expm1MinusX(d2), f3 = Tools.Expm1MinusX(d3);
            skewness = (f3 - 3d * f2 + delta3) / Math.Pow(variance, 1.5d);
            kurtosis = (Tools.Expm1MinusX(d4) - 4d * f3 + 6d * f2 + delta4) / (variance * variance);
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.LogPearsonTypeIII; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Log-Pearson Type III"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LPIII"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Mean (of log) (µ)";
                parmString[1, 0] = "Std Dev (of log) (σ)";
                parmString[2, 0] = "Skew (of log) (γ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                parmString[2, 1] = Gamma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ", "γ"]; }
        }
        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma), nameof(Gamma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma, Gamma]; }
        }

        /// <inheritdoc/>
        /// <remarks>Moment order r exists only when 1-r*Beta*log(Base) is positive. Divergent positive raw moments return positive infinity.</remarks>
        public override double Mean
        {
            get { return Math.Exp(Mu * Math.Log(Base) + LogMomentShape(1)); }
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
                double logBase = Math.Log(Base), scale = Sigma * logBase;
                if (Gamma == 0d) return Math.Exp(Mu * logBase - scale * scale);
                double beta = Beta * logBase, shape = Alpha;
                if (Gamma > 0d && shape <= 1d) return Minimum;
                if (Gamma < 0d)
                {
                    if (shape < 1d) return Maximum;
                    if (shape == 1d) return beta < -1d ? 0d : Maximum;
                    if (beta <= -1d) return 0d;
                }
                // The transformation density includes 1/x, moving the gamma mode.
                return Math.Exp(Mu * logBase - (scale * scale + beta) / (1d + beta));
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                double first = LogMomentShape(1), second = LogMomentShape(2);
                if (double.IsPositiveInfinity(first) || double.IsPositiveInfinity(second)) return double.PositiveInfinity;
                double delta = second - 2d * first;
                double logExcess = delta > 0.5d ? delta + Tools.Log1p(-Math.Exp(-delta)) : Math.Log(Tools.Expm1(delta));
                return Math.Exp(Mu * Math.Log(Base) + first + 0.5d * logExcess);
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (double.IsPositiveInfinity(LogMomentShape(3)))
                    return double.IsPositiveInfinity(LogMomentShape(2)) ? double.NaN : double.PositiveInfinity;
                double scale = Sigma * Math.Log(Base);
                if (Math.Abs(scale) * (1d + Math.Abs(Gamma)) < SmallScaleLimit)
                {
                    SmallScaleStandardizedMoments(out double skewness, out _);
                    return skewness;
                }
                double first = LogMomentShape(1), d2 = LogMomentShape(2) - 2d * first, d3 = LogMomentShape(3) - 3d * first;
                double logVariance = d2 > 0.5d ? d2 + Tools.Log1p(-Math.Exp(-d2)) : Math.Log(Tools.Expm1(d2));
                // Factor exp(d3) out of the signed numerator, then combine it with
                // the variance denominator before exponentiating the final ratio.
                double normalizedThird = -Tools.Expm1(-d3)
                    - 3d * Math.Exp(d2 - d3) * -Tools.Expm1(-d2);
                if (normalizedThird == 0d) return 0d;
                return (normalizedThird < 0d ? -1d : 1d) * Math.Exp(d3 - 1.5d * logVariance + Math.Log(Math.Abs(normalizedThird)));
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (double.IsPositiveInfinity(LogMomentShape(4)))
                    return double.IsPositiveInfinity(LogMomentShape(2)) ? double.NaN : double.PositiveInfinity;
                double scale = Sigma * Math.Log(Base);
                if (Math.Abs(scale) * (1d + Math.Abs(Gamma)) < SmallScaleLimit)
                {
                    SmallScaleStandardizedMoments(out _, out double kurtosis);
                    return kurtosis;
                }
                double first = LogMomentShape(1), d2 = LogMomentShape(2) - 2d * first;
                double d3 = LogMomentShape(3) - 3d * first, d4 = LogMomentShape(4) - 4d * first;
                double logVariance = d2 > 0.5d ? d2 + Tools.Log1p(-Math.Exp(-d2)) : Math.Log(Tools.Expm1(d2));
                // All exponential arguments in the normalized numerator are nonpositive.
                // The last exponential alone determines genuine result overflow.
                double normalizedFourth = -Tools.Expm1(-d4)
                    - 4d * Math.Exp(d3 - d4) * -Tools.Expm1(-d3)
                    + 6d * Math.Exp(d2 - d4) * -Tools.Expm1(-d2);
                return Math.Exp(d4 - 2d * logVariance + Math.Log(normalizedFourth));
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Gamma > 0d ? _boundary : 0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Gamma < 0d ? _boundary : double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            // The mean of the log10-transformed variable is a location parameter and can be any
            // finite value; only the log-space standard deviation is bounded below by zero.
            get { return [double.NegativeInfinity, 0.0d, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity]; }
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
            var newDistribution = new LogPearsonTypeIII(Mu, Sigma, Gamma) { Base = Base };
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
        /// <param name="skewOfLog">The skew of the log transformed data.</param>
        public void SetParameters(double meanOfLog, double standardDeviationOfLog, double skewOfLog)
        {
            _parametersValid = ValidateParameters(meanOfLog, standardDeviationOfLog, skewOfLog, false) is null;
            _mu = meanOfLog;
            _sigma = standardDeviationOfLog;
            _gamma = skewOfLog;
            UpdateLogBaseCache();
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mu">The mean of the log transformed data.</param>
        /// <param name="sigma">The standard deviation of the log transformed data.</param>
        /// <param name="gamma">The skew of the log transformed data.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double mu, double sigma, double gamma, bool throwException)
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
            if (double.IsNaN(gamma) || double.IsInfinity(gamma))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma must be a number.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma must be a number.");
            }
            if (gamma > 6)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be less than 6.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be less than 6.");
            }
            if (gamma < -6)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be greater than -6.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be greater than -6.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }
      
        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        /// <returns>The product moments of the base-log-transformed observations.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is insufficient, constant, nonfinite, or contains a nonpositive observation.</exception>
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
        /// <returns>The linear moments of the base-log-transformed observations.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is insufficient, constant, nonfinite, or contains a nonpositive observation.</exception>
        public double[] IndirectMethodOfLinearMoments(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, positive: true);
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) transformedSample[i] = Math.Log(sample[i], Base);
            return Statistics.LinearMoments(transformedSample);
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            return moments.ToArray().Subset(0, 2);
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new LogPearsonTypeIII() { Base = Base };
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <inheritdoc/>
        /// <remarks>The base-log Pearson type III parameters are those of <see cref="PearsonTypeIII"/> fitted to the
        /// L-moments of the base-log observations, so this method returns
        /// <see cref="PearsonTypeIII.ParametersFromLinearMoments(IList{double})"/>.</remarks>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            return new PearsonTypeIII().ParametersFromLinearMoments(moments);
        }

        /// <inheritdoc/>
        /// <remarks>Returns the L-moments of the base-log variable, which are those of <see cref="PearsonTypeIII"/>
        /// from <see cref="PearsonTypeIII.LinearMomentsFromParameters(IList{double})"/>.</remarks>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            return new PearsonTypeIII().LinearMomentsFromParameters(parameters);
        }

        /// <inheritdoc/>
        /// <remarks>Preserves the legacy 0.01 substitution for nonpositive observations when constructing
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
            //
            // Estimate initial values using the method of moments.
            var transformedSample = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++)
                transformedSample[i] = Math.Log(sample[i] > 0d ? sample[i] : 0.01d, Base);
            var mom = Statistics.ProductMoments(transformedSample);
            initialVals = [mom[0], mom[1], mom[2]];
            // Get bounds of mean. The mean is a location parameter on the log scale and is
            // legitimately negative whenever the data are mostly below 1, so the bounds are
            // symmetric about zero from the magnitude of the initial value, matching Normal's
            // location bounds. A machine-epsilon floor here would reject any sub-unity sample
            // before a fit could start.
            double real = Math.Exp(initialVals[0] / K);
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = Math.Floor(Math.Log(Math.Pow(10d, Math.Floor(Math.Log10(real)) - 1d), Base));
            upperVals[0] = Math.Ceiling(Math.Log(Math.Pow(10d, Math.Ceiling(Math.Log10(real)) + 1d), Base));
            // Get bounds of standard deviation
            real = Math.Exp(initialVals[1] / K);
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Ceiling(Math.Log(Math.Pow(10d, Math.Ceiling(Math.Log10(real) + 1d)), Base));
            upperVals[1] = double.IsNaN(upperVals[1]) ? 4 : upperVals[1];

            // Get bounds of skew
            lowerVals[2] = -6d;
            upperVals[2] = 6d;
            // Correct initial value of skew if necessary
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0.01;
            }
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
            return new PearsonTypeIII().GetRobustParameterConstraints(transformed);
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
                var LP3 = new LogPearsonTypeIII() { Base = Base };
                LP3.SetParameters(x);
                return LP3.LogLikelihood(sample);
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
        /// When the shape <c>α &lt; 1</c>, the density has a genuine integrable singularity at
        /// the transformed support boundary <c>x = Base^ξ</c> and this method intentionally
        /// returns positive infinity for either skew direction.
        /// </remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            double logBase = _logBase;
            double boundary = Gamma == 0d ? double.NaN : _boundary;
            double minimum = Gamma > 0d ? boundary : 0d;
            double maximum = Gamma < 0d ? boundary : double.PositiveInfinity;
            if (x < minimum || x > maximum || double.IsPositiveInfinity(x)) return double.NegativeInfinity;
            if (x == 0d)
            {
                if (Gamma >= 0d) return double.NegativeInfinity;
                double rate = -1d / (Beta * logBase);
                if (rate > 1d) return double.NegativeInfinity;
                if (rate < 1d || Alpha > 1d) return double.PositiveInfinity;
                return Alpha == 1d ? -Xi * logBase : double.NegativeInfinity;
            }
            double logX = Math.Log(x);
            // Preserve exact transformed endpoint density limits despite logarithm roundoff.
            double transformed = (Gamma > 0d && x == boundary) || (Gamma < 0d && x == boundary) ? Xi : logX / logBase;
            return PearsonTypeIII.LogPDF(Mu, Sigma, Gamma, transformed) - _logLogBase - logX;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return Math.Exp(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            if (double.IsNaN(probability) || probability < 0d || probability > 1d)
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between zero and one.");
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            if (probability == 0d) return Minimum;
            if (probability == 1d) return Maximum;
            return Math.Exp(PearsonTypeIII.InverseCDF(Mu, Sigma, Gamma, probability) * Math.Log(Base));
        }

        /// <summary>
        /// Returns the inverse CDF using the modified Wilson-Hilferty transformation.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>The approximate log-Pearson type III quantile in physical coordinates.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The probability is outside the closed unit interval or the distribution parameters are invalid.</exception>
        /// <remarks>
        /// Cornish-Fisher transformation (Fisher and Cornish, 1960) for abs(skew) less than or equal to 2. If abs(skew) > 2 then use Modified Wilson-Hilferty transformation (Kirby,1972).
        /// </remarks>
        public double WilsonHilfertyInverseCDF(double probability)
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
                ValidateParameters(Mu, Sigma, Gamma, true);
            // 
            return Math.Exp((Mu + Sigma * GammaDistribution.FrequencyFactorKp(Gamma, probability)) / K);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            var clone = new LogPearsonTypeIII(Mu, Sigma, Gamma);
            clone.Base = Base;
            return clone;
        }

       
        /// <inheritdoc/>
        /// <remarks>Uses public base-log mean, standard deviation and skew coordinates, including the full zero-skew limit. Maximum-likelihood covariance requires absolute skew below sqrt(2).</remarks>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            return new PearsonTypeIII(Mu, Sigma, Gamma).ParameterCovariance(sampleSize, estimationMethod);
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments && estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                return double.NaN;
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            var standardized = new PearsonTypeIII(0, 1, Gamma);
            double variance = standardized.QuantileVariance(probability, sampleSize, estimationMethod);
            double logBase = Math.Log(Base);
            double logQuantile = new PearsonTypeIII(Mu, Sigma, Gamma).InverseCDF(probability) * logBase;
            return Math.Exp(2 * (logQuantile + Math.Log(logBase) + Math.Log(Sigma)) + Math.Log(variance));
        }

        /// <inheritdoc/>
        /// <remarks>Transforms the actual Pearson quantile derivatives into physical-quantile derivatives, applying log(Base)*quantile once in each public base-log parameter coordinate.</remarks>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Mu, Sigma, Gamma, true);
            double[] gradient = new PearsonTypeIII(Mu, Sigma, Gamma).QuantileGradient(probability);
            double quantile = Math.Exp((Mu + Sigma * gradient[1]) * Math.Log(Base));
            double factor = quantile * Math.Log(Base);
            for (int i = 0; i < gradient.Length; i++) gradient[i] *= factor;
            return gradient;
        }

        /// <summary>
        /// Returns a list of partial derivatives of X given probability with respect to each moment.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>The physical-quantile gradient in public moment coordinates.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The probability is not finite and strictly interior or the distribution parameters are invalid.</exception>
        public IList<double> QuantileGradientForMoments(double probability)
        {
            return QuantileGradient(probability);
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

    }
}
