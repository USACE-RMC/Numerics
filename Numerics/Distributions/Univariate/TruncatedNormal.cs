using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using System;
using System.Collections.Generic;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Truncated Normal (Gaussian) probability distribution.
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
    /// Wikipedia contributors, "Truncated Normal distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Truncated_normal_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// <para>
    /// When the interval contains the parent mean (α ≤ 0 ≤ β, with α = (Min − µ)/σ and β = (Max − µ)/σ) the
    /// truncated mass is Z = Φ(β) − Φ(α). When the interval lies wholly in one tail of the parent, the mass,
    /// density, distribution, quantile and moment functions are evaluated from logarithms of the lower-tail
    /// probabilities of the reflected interval (Φ(−α) − Φ(−β) above the mean, Φ(β) − Φ(α) below it), so the
    /// results stay accurate where Φ(β) − Φ(α) cancels or underflows.
    /// </para>
    /// <para>
    /// Min = Max is a valid configuration that represents a point mass at Min, using the conventions of
    /// <see cref="Deterministic"/>.
    /// </para>
    /// </remarks>
    [Serializable]
    public class TruncatedNormal : UnivariateDistributionBase, IEstimation, IBootstrappable
    {
   
        /// <summary>
        /// Constructs a Truncated Normal (Gaussian) distribution with a parent Normal mean of 0.5 and standard deviation of 0.2.
        /// Min is set to 0.0 and max is set to 1.0.
        /// </summary>
        public TruncatedNormal()
        {
            SetParameters(0.5d, 0.2d, 0.0d, 1.0d);
        }

        /// <summary>
        /// Constructs a Truncated Normal (Gaussian) distribution with the given parent Normal mean and standard deviation, min and max.
        /// </summary>
        /// <param name="mean">Mean of the parent (untruncated) Normal distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the parent (untruncated) Normal distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public TruncatedNormal(double mean, double standardDeviation, double min, double max)
        {
            SetParameters(mean, standardDeviation, min, max);
        }

        private double _mu;
        private double _sigma;
        private double _min;
        private double _max;
        private double _alpha;
        private double _beta;
        private double _cdfAlpha;
        private double _z;
        private int _tail;
        private double _tailA;
        private double _tailB;
        private double _logTailA;
        private double _logTailB;
        private double _logZ;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, Min, Max, false) is null;
                _mu = value;
                UpdateStandardizedBounds();
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter σ (sigma).
        /// </summary>
        /// <remarks>
        /// Values in [0, 1E-16) are raised to 1E-16; negative values are stored as given and fail validation.
        /// </remarks>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, Min, Max, false) is null;
                _sigma = value;
                UpdateStandardizedBounds();
            }
        }

        /// <summary>
        /// Get and set the min of the distribution.
        /// </summary>
        public double Min
        {
            get { return _min; }
            set
            {
                _parametersValid = ValidateParameters(Mu, Sigma, value, Max, false) is null;
                _min = value;
                UpdateStandardizedBounds();
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
                _parametersValid = ValidateParameters(Mu, Sigma, Min, value, false) is null;
                _max = value;
                UpdateStandardizedBounds();
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 4; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.TruncatedNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Truncated Normal"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Trunc. N"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[4, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[2, 0] = "Min";
                parmString[3, 0] = "Max";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                parmString[2, 1] = Min.ToString();
                parmString[3, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ", "Min", "Max"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma), nameof(Min), nameof(Max)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma, Min, Max]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (_min == _max) return _min;
                if (_tail != 0)
                {
                    TailDensityRatios(out double ratioAlpha, out double ratioBeta);
                    return Mu + Sigma * (ratioAlpha - ratioBeta);
                }
                double alpha = _alpha, beta = _beta, Z = _z;
                if (Z == 0) return Mu;
                return Mu + Sigma * ((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / Z);
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get
            {
                if (_min == _max) return _min;
                if (_tail != 0) return TailQuantile(0.5d);
                if (_z == 0) return Mu;
                return Mu + Sigma * Normal.StandardZ((_cdfAlpha + Normal.StandardCDF(_beta)) / 2d);
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Mu < Min)
                {
                    return Min;
                }
                else if (Mu >= Min && Mu <= Max)
                {
                    return Mu;
                }
                else if (Mu > Max)
                {
                    return Max;
                }
                return default;
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (_min == _max) return 0.0d;
                if (_tail != 0)
                {
                    TailDensityRatios(out double ratioAlpha, out double ratioBeta);
                    return Math.Sqrt(Math.Pow(Sigma, 2d) * (1d + (_alpha * ratioAlpha - _beta * ratioBeta) - Math.Pow(ratioAlpha - ratioBeta, 2d)));
                }
                double alpha = _alpha, beta = _beta, Z = _z;
                if (Z == 0) return Sigma;
                return Math.Sqrt(Math.Pow(Sigma, 2d) * (1d + (alpha * Normal.StandardPDF(alpha) - beta * Normal.StandardPDF(beta)) / Z - Math.Pow((Normal.StandardPDF(alpha) - Normal.StandardPDF(beta)) / Z, 2d)));
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Reference: Pearson diagrams for truncated normal and truncated Weibull distributions. Sugiura, Gomi. 1985.
        /// A point mass (Min = Max) has an undefined skewness and returns <see cref="double.NaN"/>.
        /// </remarks>
        public override double Skewness
        {
            get
            {
                if (_min == _max) return double.NaN;
                double K0 = (Min - Mu) / Sigma;
                double K1 = (Max - Mu) / Sigma;
                double Z0, Z1;
                if (_tail != 0)
                {
                    TailDensityRatios(out Z0, out Z1);
                }
                else
                {
                    double Z = _z;
                    if (Z == 0) return 0;
                    Z0 = Normal.StandardPDF(K0) / Z;
                    Z1 = Normal.StandardPDF(K1) / Z;
                }
                double V = 1d - (K1 * Z1 - K0 * Z0) - Math.Pow(Z1 - Z0, 2d);
                return -(1d / Math.Pow(V, 3d / 2d)) * (2d * Math.Pow(Z1 - Z0, 3d) + (3d * K1 * Z1 - 3d * K0 * Z0 - 1d) * (Z1 - Z0) + Math.Pow(K1, 2d) * Z1 - Math.Pow(K0, 2d) * Z0);
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Returns the kurtosis (3 for the Normal), not the excess kurtosis. A point mass (Min = Max) has an
        /// undefined kurtosis and returns <see cref="double.NaN"/>.
        /// </remarks>
        public override double Kurtosis
        {
            get
            {
                if (_min == _max) return double.NaN;
                double K0 = (Min - Mu) / Sigma;
                double K1 = (Max - Mu) / Sigma;
                double Z0, Z1;
                if (_tail != 0)
                {
                    TailDensityRatios(out Z0, out Z1);
                }
                else
                {
                    double Z = _z;
                    if (Z == 0) return 3d;
                    Z0 = Normal.StandardPDF(K0) / Z;
                    Z1 = Normal.StandardPDF(K1) / Z;
                }
                double V = 1d - (K1 * Z1 - K0 * Z0) - Math.Pow(Z1 - Z0, 2d);
                return (1d / Math.Pow(V, 2d)) * (-3 * Math.Pow(Z1 - Z0, 4d) - 6d * (K1 * Z1 - K0 * Z0) * Math.Pow(Z1 - Z0, 2d) - 2d * Math.Pow(Z1 - Z0, 2d) - 4d * (Math.Pow(K1, 2d) * Z1 - Math.Pow(K0, 2d) * Z0) * (Z1 - Z0) - 3d * (K1 * Z1 - K0 * Z0) - (Math.Pow(K1, 3d) * Z1 - Math.Pow(K0, 3d) * Z0) + 3d);
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return _min; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return _max; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d, double.NegativeInfinity, Min]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, Max, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="sample"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the sample has fewer than four observations or contains a value that is not finite.</exception>
        /// <exception cref="NotImplementedException">Thrown for estimation methods other than the method of moments.</exception>
        /// <remarks>
        /// <para>
        /// The method of moments sets Min and Max to the sample minimum and maximum and solves for the parent µ and σ
        /// such that the truncated <see cref="Mean"/> and <see cref="StandardDeviation"/> equal the sample mean and the
        /// sample standard deviation of <see cref="Statistics.ProductMoments(IList{double})"/>.
        /// </para>
        /// <para>
        /// For a fixed σ the truncated mean increases monotonically in µ, and along the curve of matched means the
        /// truncated standard deviation increases monotonically in σ, from zero toward the standard deviation of the
        /// exponentially tilted uniform distribution on [Min, Max] with the same mean (σ → ∞). The two equations are
        /// therefore solved by nested Brent root finding, µ for each σ and then σ, starting from the sample mean and
        /// sample standard deviation. σ is limited to 5 times the sample range, beyond which the closed-form moments of
        /// the strongly tilted interval lose accuracy; at that limit the attainable standard deviation is within about
        /// 0.1% of its supremum. When the sample standard deviation is not attainable, which is always the case when it
        /// exceeds the uniform value (Max − Min)/√12, no truncated normal reproduces both moments; σ is then set to
        /// the limit and µ minimizes the sum of the squared differences between the fitted and sample means and
        /// standard deviations.
        /// </para>
        /// <para>
        /// A constant sample yields the point mass Min = Max = µ at the common value.
        /// </para>
        /// </remarks>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                if (sample == null) throw new ArgumentNullException(nameof(sample));
                if (sample.Count >= 4 && Tools.IsFinite(sample[0]))
                {
                    bool constant = true;
                    for (int i = 1; i < sample.Count && constant; i++)
                        constant = sample[i] == sample[0];
                    if (constant)
                    {
                        SetParameters(sample[0], 0d, sample[0], sample[0]);
                        return;
                    }
                }
                DistributionNumerics.ValidateSample(sample, 4);
                var mom = Statistics.ProductMoments(sample);
                var min = Statistics.Minimum(sample);
                var max = Statistics.Maximum(sample);
                SetParameters(MomentParameters(mom[0], mom[1], min, max));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Solves for the parent parameters whose truncation to [min, max] has the given mean and standard deviation.
        /// </summary>
        /// <param name="mean">The target mean, strictly between <paramref name="min"/> and <paramref name="max"/>.</param>
        /// <param name="standardDeviation">The target standard deviation, which must be positive.</param>
        /// <param name="min">The lower truncation bound.</param>
        /// <param name="max">The upper truncation bound, greater than <paramref name="min"/>.</param>
        /// <returns>The parameters {µ, σ, min, max}.</returns>
        /// <remarks>
        /// The equations are solved on the unit interval [-1, 1], to which the truncated normal family is closed under
        /// the affine map x = c + h·y with c the interval midpoint and h its half-width, and the solution is mapped back.
        /// See <see cref="Estimate(IList{double}, ParameterEstimationMethod)"/> for the method and the fallback when the
        /// standard deviation cannot be matched.
        /// </remarks>
        private static double[] MomentParameters(double mean, double standardDeviation, double min, double max)
        {
            double center = min / 2d + max / 2d, half = max / 2d - min / 2d;
            double targetMean = (mean - center) / half, targetSD = standardDeviation / half;
            var unit = new TruncatedNormal(targetMean, targetSD, -1d, 1d);
            // Largest σ, in half-widths of the interval (5 times its width). Beyond it the closed-form truncated
            // moments of the strongly tilted interval lose accuracy to cancellation.
            const double maximumSigma = 10d;

            double MeanResidual(double mu)
            {
                unit.Mu = mu;
                return unit.Mean - targetMean;
            }

            double MuForMean(double sigma)
            {
                unit.Sigma = sigma;
                double step = Math.Max(sigma, 1d), lower = targetMean - step, upper = targetMean + step;
                for (int i = 0; i < 2000 && MeanResidual(lower) > 0d; i++) { lower -= step; step *= 2d; }
                step = Math.Max(sigma, 1d);
                for (int i = 0; i < 2000 && MeanResidual(upper) < 0d; i++) { upper += step; step *= 2d; }
                return Brent.Solve(MeanResidual, lower, upper, 1E-14 * Math.Max(1d, Math.Abs(targetMean)), 1000, false);
            }

            double SDResidual(double logSigma)
            {
                double sigma = Math.Exp(logSigma);
                double mu = MuForMean(sigma);
                unit.SetParameters(mu, sigma, -1d, 1d);
                return unit.StandardDeviation - targetSD;
            }

            // The truncated variance never exceeds the parent variance, so σ = target SD is a lower bracket.
            double logLower = Math.Log(targetSD), logUpper = logLower;
            for (int i = 0; i < 100 && SDResidual(logLower) > 0d; i++) logLower -= Math.Log(2d);
            double upperResidual = SDResidual(logUpper);
            while (upperResidual < 0d && logUpper < Math.Log(maximumSigma))
            {
                logUpper = Math.Min(logUpper + Math.Log(4d), Math.Log(maximumSigma));
                upperResidual = SDResidual(logUpper);
            }

            double muUnit, sigmaUnit;
            if (upperResidual >= 0d)
            {
                sigmaUnit = Math.Exp(Brent.Solve(SDResidual, logLower, logUpper, 1E-14, 1000, false));
                muUnit = MuForMean(sigmaUnit);
            }
            else
            {
                // No solution: least squares on the mean and standard deviation at the largest σ.
                sigmaUnit = maximumSigma;
                double matched = MuForMean(sigmaUnit);
                double LeastSquares(double mu)
                {
                    unit.SetParameters(mu, sigmaUnit, -1d, 1d);
                    return Math.Pow(unit.Mean - targetMean, 2d) + Math.Pow(unit.StandardDeviation - targetSD, 2d);
                }
                if (matched == 0d)
                {
                    muUnit = 0d;
                }
                else
                {
                    var search = new BrentSearch(LeastSquares, Math.Min(matched, 0d), Math.Max(matched, 0d))
                    {
                        RelativeTolerance = 1E-12
                    };
                    search.Minimize();
                    muUnit = search.BestParameterSet.Values[0];
                }
            }
            return [center + half * muUnit, half * sigmaUnit, min, max];
        }

        /// <inheritdoc/>
        /// <exception cref="Exception">Thrown when the parameters estimated from the bootstrap sample are invalid.</exception>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new TruncatedNormal(Mu, Sigma, Min, Max);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }


        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="mean">Mean of the parent (untruncated) Normal distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the parent (untruncated) Normal distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public void SetParameters(double mean, double standardDeviation, double min, double max)
        {
            // Set parameters
            Mu = mean;
            Sigma = standardDeviation;
            Min = min;
            Max = max;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2], parameters[3]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mean">Mean of the parent (untruncated) Normal distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the parent (untruncated) Normal distribution.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double mean, double standardDeviation, double min, double max, bool throwException)
        {
            if (double.IsNaN(mean) || double.IsInfinity(mean))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mean must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mean must be a number.");
            }
            if (double.IsNaN(standardDeviation) || double.IsInfinity(standardDeviation) || standardDeviation <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
            }
            if (double.IsNaN(min) || double.IsNaN(max) || double.IsInfinity(min) || double.IsInfinity(max) || min > max)
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
            return ValidateParameters(parameters[0], parameters[1], parameters[2], parameters[3], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Min, Max, true);
            // 
            if (x < Minimum || x > Maximum) return 0.0d;
            if (_min == _max) return 1.0d;
            double Xi = (x - Mu) / Sigma;
            if (_tail != 0) return Math.Exp(Normal.StandardLogPDF(Xi) - _logZ) / Sigma;
            double Z = _z;
            if (Z == 0) return Normal.StandardPDF(Xi);
            return Normal.StandardPDF(Xi) / (Sigma * Z);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Min, Max, true);
            if (_min == _max) return x < _min ? 0d : 1d;
            if (x <= Minimum) return 0d;
            if (x >= Maximum) return 1d;
            double Xi = (x - Mu) / Sigma;
            if (_tail != 0)
            {
                // Lower-tail probability of the (possibly reflected) standardized value, kept inside the bounds.
                double logTail = Tools.Clamp(DistributionNumerics.NormalLogCDF(_tail > 0 ? -Xi : Xi), _logTailA, _logTailB);
                double logMass = _tail > 0 ? DistributionNumerics.LogDifference(_logTailB, logTail) : DistributionNumerics.LogDifference(logTail, _logTailA);
                return Math.Min(1d, Math.Exp(logMass - _logZ));
            }
            double Z = _z;
            if (Z == 0) return Normal.StandardCDF(Xi);
            return (Normal.StandardCDF(Xi) - _cdfAlpha) / Z;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Min, Max, true);
            if (_min == _max) return _min;
            if (_tail != 0) return TailQuantile(probability);
            double Z = _z;
            if (Z == 0) return Mu + Sigma * Normal.StandardZ(probability);
            return Mu + Sigma * Normal.StandardZ(_cdfAlpha + probability * Z);
        }

        /// <summary>
        /// Caches the standardized bounds, Φ(α) and the truncated mass Z = Φ(β) - Φ(α) shared by the density,
        /// distribution, quantile and moment functions.
        /// </summary>
        /// <remarks>
        /// <para>
        /// When the interval lies wholly above the parent mean (α &gt; 0) or wholly below it (β &lt; 0), the cache also
        /// holds the bounds [a, b] of the equivalent lower-tail interval ([−β, −α] above the mean, [α, β] below it),
        /// log Φ(a), log Φ(b) and log Z = log(Φ(b) − Φ(a)). An interval too narrow for the logarithmic difference to
        /// resolve uses the midpoint rule log φ((a + b)/2) + log(b − a).
        /// </para>
        /// <para>
        /// The cache is written only when a parameter changes, so concurrent evaluations on one instance do not
        /// share intermediate state.
        /// </para>
        /// </remarks>
        private void UpdateStandardizedBounds()
        {
            _alpha = (_min - _mu) / _sigma;
            _beta = (_max - _mu) / _sigma;
            _cdfAlpha = Normal.StandardCDF(_alpha);
            _z = Normal.StandardCDF(_beta) - _cdfAlpha;
            _tail = _alpha > 0d ? 1 : _beta < 0d ? -1 : 0;
            if (_tail == 0) return;
            _tailA = _tail > 0 ? -_beta : _alpha;
            _tailB = _tail > 0 ? -_alpha : _beta;
            _logTailA = DistributionNumerics.NormalLogCDF(_tailA);
            _logTailB = DistributionNumerics.NormalLogCDF(_tailB);
            _logZ = DistributionNumerics.LogDifference(_logTailB, _logTailA);
            if (!(_logZ > double.NegativeInfinity))
                _logZ = Normal.StandardLogPDF(_tailA / 2d + _tailB / 2d) + Math.Log(_tailB - _tailA);
        }

        /// <summary>
        /// Returns φ(α)/Z and φ(β)/Z for an interval in one tail of the parent, evaluated from logarithms.
        /// </summary>
        /// <param name="ratioAlpha">The standard Normal density at α divided by the truncated mass.</param>
        /// <param name="ratioBeta">The standard Normal density at β divided by the truncated mass.</param>
        private void TailDensityRatios(out double ratioAlpha, out double ratioBeta)
        {
            ratioAlpha = Math.Exp(Normal.StandardLogPDF(_alpha) - _logZ);
            ratioBeta = Math.Exp(Normal.StandardLogPDF(_beta) - _logZ);
        }

        /// <summary>
        /// Returns the quantile of an interval in one tail of the parent.
        /// </summary>
        /// <param name="probability">The non-exceedance probability, strictly between 0 and 1.</param>
        /// <returns>The quantile, kept inside [Min, Max].</returns>
        /// <remarks>
        /// Works on the equivalent lower-tail interval [a, b]. The target log Φ(u) is log(Φ(a) + qZ) when the lower-tail
        /// probability q is at most one half and log(Φ(b) − (1 − q)Z) otherwise, so neither form subtracts nearly equal
        /// probabilities. The inversion uses <see cref="Normal.StandardZ(double)"/> while Φ(u) is a normal double, and
        /// otherwise solves log Φ(u) equal to the target with Brent's method on [a, b].
        /// </remarks>
        private double TailQuantile(double probability)
        {
            double q = _tail > 0 ? 1d - probability : probability;
            double complement = _tail > 0 ? probability : 1d - probability;
            double logTarget = q <= 0.5d
                ? DistributionNumerics.LogSum(_logTailA, Math.Log(q) + _logZ)
                : DistributionNumerics.LogDifference(_logTailB, Math.Min(Math.Log(complement) + _logZ, _logTailB));
            logTarget = Tools.Clamp(logTarget, _logTailA, _logTailB);
            double u = logTarget > -700d
                ? Normal.StandardZ(Math.Exp(logTarget))
                : Brent.Solve(t => DistributionNumerics.NormalLogCDF(t) - logTarget, _tailA, _tailB, 1E-15 * Math.Abs(_tailB), 1000, false);
            u = Tools.Clamp(u, _tailA, _tailB);
            double x = Mu + Sigma * (_tail > 0 ? -u : u);
            return Tools.Clamp(x, Min, Max);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new TruncatedNormal(Mu, Sigma, Min, Max);
        }
          
    }
}