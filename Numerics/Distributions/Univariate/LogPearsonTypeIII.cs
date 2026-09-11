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
        /// Gets and sets the scale parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return 0.5d * Sigma * Gamma; }
        }

        /// <summary>
        /// Gets and sets the shape parameter α (alpha).
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
            return Alpha * (-Tools.Log1p(-argument) - argument);
        }

        /// <summary>Expands centered exponential moments before evaluating them, avoiding cancellation for tiny log scale.</summary>
        /// <param name="skewness">The resulting standardized third central moment.</param>
        /// <param name="kurtosis">The resulting standardized fourth central moment.</param>
        private void SmallScaleStandardizedMoments(out double skewness, out double kurtosis)
        {
            const int order = 12;
            var raw = new double[4][];
            for (int r = 1; r <= 4; r++)
            {
                raw[r - 1] = new double[order + 1];
                raw[r - 1][0] = 1d;
                for (int n = 2; n <= order; n++)
                    for (int k = 2; k <= n; k++)
                        raw[r - 1][n] += Math.Pow(r, k) * Math.Pow(Gamma / 2d, k - 2) * raw[r - 1][n - k] / n;
            }
            double[] square = MultiplySeries(raw[0], raw[0]), cube = MultiplySeries(square, raw[0]);
            double[] fourth = MultiplySeries(square, square), secondFirst = MultiplySeries(raw[1], raw[0]);
            double[] thirdFirst = MultiplySeries(raw[2], raw[0]), secondSquare = MultiplySeries(raw[1], square);
            double scale = Sigma * Math.Log(Base), variance = 0d, third = 0d, fourthCentral = 0d;
            for (int n = order; n >= 2; n--) variance = variance * scale + raw[1][n] - square[n];
            for (int n = order; n >= 3; n--) third = third * scale + raw[2][n] - 3d * secondFirst[n] + 2d * cube[n];
            for (int n = order; n >= 4; n--) fourthCentral = fourthCentral * scale + raw[3][n] - 4d * thirdFirst[n] + 6d * secondSquare[n] - 3d * fourth[n];
            skewness = third / Math.Pow(variance, 1.5d);
            kurtosis = fourthCentral / (variance * variance);
        }

        /// <summary>Multiplies equal-length truncated power series used for centered moment evaluation.</summary>
        /// <param name="left">The first coefficient vector.</param>
        /// <param name="right">The second coefficient vector of the same length.</param>
        /// <returns>The product truncated to the input vector length.</returns>
        private static double[] MultiplySeries(double[] left, double[] right)
        {
            var result = new double[left.Length];
            for (int n = 0; n < result.Length; n++)
                for (int k = 0; k <= n; k++) result[n] += left[k] * right[n - k];
            return result;
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
                if (Math.Abs(scale) * (1d + Math.Abs(Gamma)) < 0.01d)
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
                if (Math.Abs(scale) * (1d + Math.Abs(Gamma)) < 0.01d)
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
            get { return Gamma > 0d ? Math.Exp(Xi * Math.Log(Base)) : 0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return Gamma < 0d ? Math.Exp(Xi * Math.Log(Base)) : double.PositiveInfinity; }
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
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            if (T3 == 0.0d)
            {
                return [L1, L2 * Math.Sqrt(Math.PI), 0.0d];
            }
            var alpha = default(double);
            double z;
            // The following approximation has relative accuracy better than 5x10-5 for all values of alpha.
            if (Math.Abs(T3) > 0.0d && Math.Abs(T3) < 1d / 3d)
            {
                z = 3.0d * Math.PI * Math.Pow(T3, 2d);
                alpha = (1.0d + 0.2906d * z) / (z + 0.1882d * Math.Pow(z, 2d) + 0.0442d * Math.Pow(z, 3d));
            }
            else if (Math.Abs(T3) >= 1d / 3d && Math.Abs(T3) < 1.0d)
            {
                z = 1.0d - Math.Abs(T3);
                alpha = (0.36067d * z - 0.59567d * Math.Pow(z, 2d) + 0.25361d * Math.Pow(z, 3d)) / (1.0d - 2.78861d * z + 2.56096d * Math.Pow(z, 2d) - 0.77045d * Math.Pow(z, 3d));
            }

            double mu = L1;
            double gamma = 2.0d * Math.Pow(alpha, -0.5d) * Math.Sign(T3);
            double sigma;
            if (alpha < 100d)
            {
                sigma = L2 * Math.Sqrt(Math.PI) * Math.Sqrt(alpha) * Mathematics.SpecialFunctions.Gamma.Function(alpha) / Mathematics.SpecialFunctions.Gamma.Function(alpha + 0.5d);
            }
            else
            {
                double inverseAlpha = 1.0d / alpha;
                double correction = 1.0d - inverseAlpha / 8.0d + inverseAlpha * inverseAlpha / 128.0d;
                sigma = Math.Sqrt(Math.PI) * L2 / correction;
            }

            return [mu, sigma, gamma];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double mu = parameters[0];
            double sigma = parameters[1];
            double gamma = parameters[2];
            if (gamma == 0.0d)
            {
                return [mu, sigma / Math.Sqrt(Math.PI), 0.0d, 0.12260172d];
            }
            double alpha = 4.0d / Math.Pow(gamma, 2d);
            double beta = 0.5d * sigma * gamma;
            double L1 = mu;
            double L2;
            if (alpha < 100.0d)
            {
                L2 = Math.Abs(beta * Mathematics.SpecialFunctions.Gamma.Function(alpha + 0.5d) / (Math.Sqrt(Math.PI) * Mathematics.SpecialFunctions.Gamma.Function(alpha)));
            }
            else
            {
                double inverseAlpha = 1.0d / alpha;
                double correction = 1.0d - inverseAlpha / 8.0d + inverseAlpha * inverseAlpha / 128.0d;
                L2 = sigma / Math.Sqrt(Math.PI) * correction;
            }
            // The following approximations are accurate to 10-6. 
            double A0 = 0.32573501d;
            double A1 = 0.1686915d;
            double A2 = 0.078327243d;
            double A3 = -0.0029120539d;
            double B1 = 0.46697102d;
            double B2 = 0.24255406d;
            double C0 = 0.12260172d;
            double C1 = 0.05373013d;
            double C2 = 0.043384378d;
            double C3 = 0.011101277d;
            double D1 = 0.18324466d;
            double D2 = 0.20166036d;
            double E1 = 2.3807576d;
            double E2 = 1.5931792d;
            double E3 = 0.11618371d;
            double F1 = 5.1533299d;
            double F2 = 7.142526d;
            double F3 = 1.9745056d;
            double G1 = 2.1235833d;
            double G2 = 4.1670213d;
            double G3 = 3.1925299d;
            double H1 = 9.0551443d;
            double H2 = 26.649995d;
            double H3 = 26.193668d;
            double T3;
            double T4;
            if (alpha >= 1d)
            {
                T3 = Math.Pow(alpha, -0.5d) * (A0 + A1 * Math.Pow(alpha, -1) + A2 * Math.Pow(alpha, -2) + A3 * Math.Pow(alpha, -3)) / (1d + B1 * Math.Pow(alpha, -1) + B2 * Math.Pow(alpha, -2));
                T4 = (C0 + C1 * Math.Pow(alpha, -1) + C2 * Math.Pow(alpha, -2) + C3 * Math.Pow(alpha, -3)) / (1d + D1 * Math.Pow(alpha, -1) + D2 * Math.Pow(alpha, -2));
            }
            else
            {
                T3 = (1d + E1 * alpha + E2 * Math.Pow(alpha, 2d) + E3 * Math.Pow(alpha, 3d)) / (1d + F1 * alpha + F2 * Math.Pow(alpha, 2d) + F3 * Math.Pow(alpha, 3d));
                T4 = (1d + G1 * alpha + G2 * Math.Pow(alpha, 2d) + G3 * Math.Pow(alpha, 3d)) / (1d + H1 * alpha + H2 * Math.Pow(alpha, 2d) + H3 * Math.Pow(alpha, 3d));
            }
            T3 *= Math.Sign(gamma);

            return [L1, L2, T3, T4];
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
            double logBase = Math.Log(Base);
            double boundary = Gamma == 0d ? double.NaN : Math.Exp(Xi * logBase);
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
            return PearsonTypeIII.LogPDF(Mu, Sigma, Gamma, transformed) - Math.Log(logBase) - logX;
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
