using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Noncentral t probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Noncentral_t-distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class NoncentralT : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
    {
 
        /// <summary>
        /// Constructs a Noncentral t distribution with 10 degrees of freedom and noncentrality = 0.
        /// </summary>
        public NoncentralT()
        {
            SetParameters(10d, 0d);
        }

        /// <summary>
        /// Constructs a Noncentral t distribution with given degrees of freedom and noncentrality.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        /// <param name="noncentrality">The noncentrality parameter.</param>
        public NoncentralT(double degreesOfFreedom, double noncentrality)
        {
            SetParameters(degreesOfFreedom, noncentrality);
        }

        private double _degreesOfFreedom;
        private double _noncentrality;
       
        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public double DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(value, Noncentrality, false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <summary>
        /// Gets and sets the noncentrality parameter μ (mu) of the distribution.
        /// </summary>
        public double Noncentrality
        {
            get { return _noncentrality; }
            set
            {
                _parametersValid = ValidateParameters(DegreesOfFreedom, value, false) is null;
                _noncentrality = value;
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
            get { return UnivariateDistributionType.NoncentralT; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Noncentral t"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "NCT"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Degrees of Freedom (ν)";
                parmString[1, 0] = "Noncentrality (μ)";
                parmString[0, 1] = DegreesOfFreedom.ToString();
                parmString[1, 1] = Noncentrality.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ν", "μ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(DegreesOfFreedom), nameof(Noncentrality)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [DegreesOfFreedom, Noncentrality]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (_degreesOfFreedom > 1)
                {
                    return Noncentrality * Math.Sqrt(DegreesOfFreedom / 2d) * Gamma.RatioHalf((DegreesOfFreedom - 1d) / 2d);
                }
                else
                {
                    return double.NaN;
                }
            }
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
                var brent = new BrentSearch(PDF, InverseCDF(0.001), InverseCDF(0.999));
                brent.Maximize();
                return brent.BestParameterSet.Values[0];
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (DegreesOfFreedom > 2)
                {
                    double a = DegreesOfFreedom * (1d + Math.Pow(Noncentrality, 2d)) / (DegreesOfFreedom - 2);
                    double b = DegreesOfFreedom * Math.Pow(Noncentrality, 2d) / 2d;
                    double c = Gamma.RatioHalf((DegreesOfFreedom - 1d) / 2d);
                    return Math.Sqrt(a - b * c * c);
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Computed in closed form from the raw moments E[T^k] = (ν/2)^(k/2) Γ((ν−k)/2) / Γ(ν/2) E[Z^k], with
        /// Z ~ N(μ, 1), where μ is the noncentrality parameter, E[Z^2] = μ² + 1 and E[Z^3] = μ³ + 3μ. With m = E[T] and
        /// s² = Var[T], the third central moment is m [ν(2ν − 3 + μ²) / ((ν − 2)(ν − 3)) − 2s²].
        /// </para>
        /// <para>
        /// The skewness exists only for ν &gt; 3; otherwise NaN is returned.
        /// </para>
        /// </remarks>
        public override double Skewness
        {
            get
            {
                double v = DegreesOfFreedom;
                if (!(v > 3)) return double.NaN;
                double d2 = Noncentrality * Noncentrality;
                double mu = Noncentrality * Math.Sqrt(v / 2d) * Gamma.RatioHalf((v - 1d) / 2d);
                double variance = v * (1d + d2) / (v - 2d) - mu * mu;
                double m3 = mu * (v * (2d * v - 3d + d2) / ((v - 2d) * (v - 3d)) - 2d * variance);
                return m3 / Math.Pow(variance, 1.5);
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Returns the kurtosis (not the excess kurtosis), computed in closed form from the raw moments
        /// E[T^k] = (ν/2)^(k/2) Γ((ν−k)/2) / Γ(ν/2) E[Z^k], with Z ~ N(μ, 1), where μ is the noncentrality parameter, and
        /// E[Z^4] = μ⁴ + 6μ² + 3. With m = E[T] and s² = Var[T], the fourth central moment is
        /// ν²(μ⁴ + 6μ² + 3) / ((ν − 2)(ν − 4)) − m² [ν((ν + 1)μ² + 3(3ν − 5)) / ((ν − 2)(ν − 3)) − 3s²].
        /// </para>
        /// <para>
        /// The kurtosis exists only for ν &gt; 4; otherwise NaN is returned.
        /// </para>
        /// </remarks>
        public override double Kurtosis
        {
            get
            {
                double v = DegreesOfFreedom;
                if (!(v > 4)) return double.NaN;
                double d2 = Noncentrality * Noncentrality;
                double mu = Noncentrality * Math.Sqrt(v / 2d) * Gamma.RatioHalf((v - 1d) / 2d);
                double variance = v * (1d + d2) / (v - 2d) - mu * mu;
                double e4 = v * v * (d2 * d2 + 6d * d2 + 3d) / ((v - 2d) * (v - 4d));
                double m4 = e4 - mu * mu * (v * ((v + 1d) * d2 + 3d * (3d * v - 5d)) / ((v - 2d) * (v - 3d)) - 3d * variance);
                return m4 / (variance * variance);
            }
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
            get { return [1.0d, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        public void SetParameters(double v, double mu)
        {
            _degreesOfFreedom = v;
            _noncentrality = mu;
            _parametersValid = ValidateParameters(v, mu, false) is null;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        /// <param name="mu">The noncentrality parameter μ (mu).</param>
        /// <param name="throwException">Whether to throw the validation exception.</param>
        /// <returns><see langword="null"/> when the parameters are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="throwException"/> is true and either parameter is invalid.</exception>
        public ArgumentOutOfRangeException? ValidateParameters(double v, double mu, bool throwException)
        {
            if (double.IsNaN(v) || double.IsInfinity(v) || v < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            if (double.IsNaN(mu) || double.IsInfinity(mu))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Noncentrality), "The noncentrality parameter μ (mu) must be a number.");
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
        /// Only <see cref="ParameterEstimationMethod.MaximumLikelihood"/> is supported, with the degrees of freedom and
        /// the noncentrality estimated jointly. The distribution has no location or scale parameter.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than four values, a nonfinite value, or no variation.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not maximum likelihood.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 4);
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
            var newDistribution = new NoncentralT(DegreesOfFreedom, Noncentrality);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The initial values match the sample mean m and variance s² to the moments of the distribution,
        /// E[T] = μ c(ν) and Var[T] = ν(1 + μ²)/(ν − 2) − E[T]², with c(ν) = √(ν/2) Γ((ν − 1)/2) / Γ(ν/2), by a few
        /// fixed-point updates μ = m/c(ν), ν = 2(s² + m²)/(s² + m² − 1 − μ²) from ν = 30. The degrees of freedom are
        /// limited to [2.5, 1000] and set to 1000 when the update has no positive solution.
        /// </para>
        /// <para>
        /// The degrees of freedom are bounded by [1, 1000] and the noncentrality by ±10^(⌈log10(max(|μ|, 1))⌉ + 1).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than four values, a nonfinite value, or no variation.</exception>
        /// <exception cref="ArgumentException">The sample moments are not finite.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            var moments = Statistics.ProductMoments(sample);
            double mean = moments[0];
            double secondMoment = moments[1] * moments[1] + mean * mean;
            if (!Tools.IsFinite(mean) || !Tools.IsFinite(secondMoment))
                throw new ArgumentException("Sample moments must be finite.", nameof(sample));
            double v = 30d;
            double noncentrality = mean;
            for (int i = 0; i < 20; i++)
            {
                noncentrality = mean / (Math.Sqrt(v / 2d) * Gamma.RatioHalf((v - 1d) / 2d));
                double denominator = secondMoment - 1d - noncentrality * noncentrality;
                v = denominator > 0d ? Tools.Clamp(2d * secondMoment / denominator, 2.5d, 1000d) : 1000d;
            }
            noncentrality = mean / (Math.Sqrt(v / 2d) * Gamma.RatioHalf((v - 1d) / 2d));
            double noncentralityBound = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Max(Math.Abs(noncentrality), 1d)) + 1d));
            var initialVals = new double[] { v, noncentrality };
            var lowerVals = new double[] { 1d, -noncentralityBound };
            var upperVals = new double[] { 1000d, noncentralityBound };
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The log-likelihood is maximized over the degrees of freedom and noncentrality with the Nelder-Mead method,
        /// within the bounds from <see cref="GetParameterConstraints(IList{double})"/>.
        /// </remarks>
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
                var T = new NoncentralT();
                T.SetParameters(x);
                return T.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// At x = 0 the density is Γ((ν + 1)/2) e^(−μ²/2) / (√(πν) Γ(ν/2)), with the gamma ratio from
        /// exp(ln B(ν/2, 1/2) − ln Γ(1/2)) for ν ≥ 2, which does not overflow for large ν.
        /// </para>
        /// <para>
        /// Elsewhere the density is evaluated from its integral representation over S = √(χ²_ν/ν),
        /// f(x) = ∫ s φ(xs − μ) f_S(s) ds, whose integrand is positive and log-concave; see
        /// <see cref="LogPDF(double)"/>. This keeps full relative precision in both tails and near x = 0.
        /// </para>
        /// </remarks>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            double u = Noncentrality;
            double v = DegreesOfFreedom;
            if (x != 0)
            {
                return Math.Exp(LogDensity(x, v, u));
            }
            else
            {
                // Γ((ν + 1)/2) / Γ(ν/2) is the reciprocal of Γ(ν/2) / Γ(ν/2 + 1/2).
                return Math.Exp(-(u * u) / 2) / (Gamma.RatioHalf(v / 2) * Math.Sqrt(Math.PI * v));
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// With S = √(χ²_ν/ν), the density is f(x) = (2/√(2π)) ∫₀^∞ exp(g(s)) ds, where
        /// g(s) = ln(w^(ν/2) e^(−w) / Γ(ν/2)) − (xs − μ)²/2 and w = νs²/2. The integrand is maximized at
        /// s* = (xμ + √(x²μ² + 4ν(x² + ν))) / (2(x² + ν)) and g''(s) = −ν/s² − x² − ν, so the integral is evaluated
        /// by adaptive Gauss-Kronrod quadrature of exp(g(s) − g(s*)) on the window where g(s) &gt; g(s*) − 45,
        /// scaled by the local curvature. The log density is therefore finite far beyond the range where the density
        /// underflows.
        /// </para>
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            if (x == 0) return Math.Log(PDF(0d));
            return LogDensity(x, DegreesOfFreedom, Noncentrality);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Let λ = μ²/2, p_j = e^(−λ) λ^j / j! and q_j = μ e^(−λ) λ^j / (√2 Γ(j + 3/2)). For x &gt; 0 the CDF is
        /// Φ(−μ) + ½ Σ_j [p_j I_y(j + ½, ν/2) + q_j I_y(j + 1, ν/2)], y = x²/(x² + ν) (Guenther 1978; Lenth 1989).
        /// The Poisson-weighted series is started at its modal term j = ⌊λ⌋ with weights evaluated in log space and
        /// summed in both directions with the incomplete beta recurrences, so it keeps its mass for any noncentrality.
        /// Summation stops when an upper bound on the remaining terms falls below machine precision relative to the sum.
        /// </para>
        /// <para>
        /// For x &lt; 0 and μ ≤ 0 the CDF equals the upper tail at −x for noncentrality −μ, which is summed directly as
        /// ½ Σ_j [p_j I_{1−y}(ν/2, j + ½) + q_j I_{1−y}(ν/2, j + 1)] with all terms positive. For x &lt; 0 and μ &gt; 0
        /// the series would subtract from Φ(−μ), so the CDF is evaluated from the positive integral
        /// F(x) = ∫₀^∞ φ(μ + u) P(ν/2, νu²/(2x²)) du, with P the regularized lower incomplete gamma function, by
        /// adaptive Gauss-Kronrod quadrature around the maximum of the log-concave integrand.
        /// </para>
        /// <para>
        /// CDF(+∞) = 1, CDF(−∞) = 0 and CDF(NaN) = NaN.
        /// </para>
        /// <para>
        /// References: W. C. Guenther (1978), "Evaluation of probabilities for the noncentral distributions and the
        /// difference of two t-variables with a desk calculator", Journal of Statistical Computation and Simulation,
        /// 6(3-4), 199-206, https://doi.org/10.1080/00949657808810188; R. V. Lenth (1989), "Algorithm AS 243:
        /// Cumulative distribution function of the non-central t distribution", Journal of the Royal Statistical
        /// Society, Series C (Applied Statistics), 38(1), 185-189, https://doi.org/10.2307/2347693; D. Benton and
        /// K. Krishnamoorthy (2003), "Computing discrete mixtures of continuous distributions: noncentral chisquare,
        /// noncentral t and the distribution of the square of the sample multiple correlation coefficient",
        /// Computational Statistics &amp; Data Analysis, 43(2), 249-267, https://doi.org/10.1016/S0167-9473(02)00283-9.
        /// </para>
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            return LowerTail(x, DegreesOfFreedom, Noncentrality);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated directly by reflection, P(T &gt; x) = F(−x; ν, −μ), so the upper tail keeps full relative
        /// precision instead of being formed as 1 − CDF.
        /// </remarks>
        public override double CCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            return LowerTail(-x, DegreesOfFreedom, -Noncentrality);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// A probability outside [0, 1], including NaN, throws <see cref="ArgumentOutOfRangeException"/>. The quantile is
        /// bracketed from the Jennett and Welch starting estimate and refined with Brent's method to near machine
        /// precision in x.
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_degreesOfFreedom, Noncentrality, true);
            // 
            return InverseCdf(probability, DegreesOfFreedom, Noncentrality);
        }

        /// <summary>
        /// Returns the log density of the noncentral t distribution at a nonzero point.
        /// </summary>
        /// <param name="t">The point at which to evaluate the density, not zero.</param>
        /// <param name="df">The degrees of freedom ν.</param>
        /// <param name="delta">The noncentrality parameter μ.</param>
        /// <returns>The natural logarithm of the density at <paramref name="t"/>.</returns>
        /// <remarks>
        /// Integrates f(t) = (2/√(2π)) ∫₀^∞ exp(g(s)) ds, g(s) = ln(w^(ν/2) e^(−w) / Γ(ν/2)) − (ts − μ)²/2,
        /// w = νs²/2, around its exact maximizer s*. Writing s = s* + σd with σ = (t² + ν + ν/s*²)^(−1/2),
        /// g(s) − g(s*) = ν[ln(1 + σd/s*) − σd(2s* + σd)/2] − tσd(t(2s* + σd) − 2μ)/2 has no large cancelling terms.
        /// Since g'' ≤ −(ν/s² + t² + ν), g falls below g(s*) − 45 for d &lt; −√90 and for
        /// d &gt; √90 / (σ √(t² + ν)), so the integral is taken over that window.
        /// </remarks>
        private static double LogDensity(double t, double df, double delta)
        {
            if (double.IsNaN(t)) return double.NaN;
            if (double.IsInfinity(t)) return double.NegativeInfinity;
            double absT = Math.Abs(t);
            double sign = Math.Sign(t);
            // Exact maximizer s* of g, from (t² + ν)s² − tμs − ν = 0, scaled by |t| when |t| ≥ 1 to avoid overflow,
            // and rationalized when tμ < 0 to avoid cancellation.
            double sStar;
            double sqrtCurvatureTerm;
            if (absT >= 1d)
            {
                double ratio = df / (absT * absT);
                double root = Math.Sqrt(delta * delta + 4d * df * (1d + ratio));
                sStar = sign * delta >= 0d
                    ? (sign * delta + root) / (2d * absT * (1d + ratio))
                    : 2d * df / (absT * (root - sign * delta));
                sqrtCurvatureTerm = absT * Math.Sqrt(1d + ratio);
            }
            else
            {
                double sum = t * t + df;
                double root = Math.Sqrt(t * t * delta * delta + 4d * df * sum);
                sStar = t * delta >= 0d
                    ? (t * delta + root) / (2d * sum)
                    : 2d * df / (root - t * delta);
                sqrtCurvatureTerm = Math.Sqrt(sum);
            }
            // Local curvature scale σ = (t² + ν + ν/s*²)^(−1/2), formed without overflow.
            double c = Math.Sqrt(df) / sStar;
            double m = Math.Max(sqrtCurvatureTerm, c);
            double sigma = 1d / (m * Math.Sqrt((sqrtCurvatureTerm / m) * (sqrtCurvatureTerm / m) + (c / m) * (c / m)));
            // Log integrand at the maximum.
            double a = df / 2d;
            double w = a * sStar * sStar;
            double logKernel = w >= 1E-300
                ? DistributionNumerics.GammaLogDensity(a, w) + Math.Log(w)
                : a * (Math.Log(a) + 2d * Math.Log(sStar)) - w - Gamma.LogGamma(a);
            double residual = t * sStar - delta;
            double gStar = logKernel - 0.5d * residual * residual;
            // Integrate exp(g(s* + σd) − g(s*)) over d.
            double lower = Math.Max(-sStar / sigma, -Math.Sqrt(90d));
            double upper = Math.Sqrt(90d) / (sigma * sqrtCurvatureTerm);
            double integrand(double d)
            {
                double step = sigma * d;
                double g = df * (Tools.Log1p(step / sStar) - 0.5d * step * (2d * sStar + step))
                    - 0.5d * t * step * (t * (2d * sStar + step) - 2d * delta);
                return Math.Exp(g);
            }
            double integral = Integrate(integrand, lower, upper);
            return Math.Log(2d) - Tools.LogSqrt2PI + gStar + Math.Log(sigma) + Math.Log(integral);
        }

        /// <summary>
        /// Integrates a smooth, positive, unimodal function over a finite interval with adaptive Gauss-Kronrod
        /// quadrature at near machine precision.
        /// </summary>
        /// <param name="function">The integrand, scaled so that its maximum is about one.</param>
        /// <param name="min">The lower limit.</param>
        /// <param name="max">The upper limit.</param>
        /// <returns>The integral.</returns>
        private static double Integrate(Func<double, double> function, double min, double max)
        {
            if (!(max > min)) return 0d;
            var quadrature = new AdaptiveGaussKronrod(function, min, max)
            {
                RelativeTolerance = 1E-14,
                AbsoluteTolerance = 1E-15,
                ReportFailure = false
            };
            quadrature.Integrate();
            return quadrature.Result;
        }

        /// <summary>
        /// Returns the cumulative probability P(T ≤ t) of the noncentral t distribution.
        /// </summary>
        /// <param name="t">The point at which to evaluate the probability.</param>
        /// <param name="df">The degrees of freedom ν.</param>
        /// <param name="delta">The noncentrality parameter μ.</param>
        /// <returns>The lower-tail probability at <paramref name="t"/>.</returns>
        private static double LowerTail(double t, double df, double delta)
        {
            if (double.IsNaN(t)) return double.NaN;
            if (double.IsPositiveInfinity(t)) return 1d;
            if (double.IsNegativeInfinity(t)) return 0d;
            if (t == 0d) return Normal.StandardCDF(-delta);
            if (t > 0d)
            {
                double normalPart = Normal.StandardCDF(-delta);
                return Math.Min(1d, normalPart + PoissonMixtureSeries(t, df, delta, false, normalPart));
            }
            if (delta <= 0d) return PoissonMixtureSeries(-t, df, -delta, true, 0d);
            return NegativeTailIntegral(t, df, delta);
        }

        /// <summary>
        /// Sums the Poisson-weighted incomplete beta series of the noncentral t distribution from its modal term.
        /// </summary>
        /// <param name="t">The positive, finite point.</param>
        /// <param name="df">The degrees of freedom ν.</param>
        /// <param name="delta">The noncentrality parameter μ.</param>
        /// <param name="upper"><see langword="false"/> to sum ½ Σ_j [p_j I_y(j + ½, ν/2) + q_j I_y(j + 1, ν/2)];
        /// <see langword="true"/> to sum the complementary series ½ Σ_j [p_j I_{1−y}(ν/2, j + ½) + q_j I_{1−y}(ν/2, j + 1)],
        /// with y = t²/(t² + ν).</param>
        /// <param name="offset">The magnitude of any value the sum is added to, used in the relative stopping rule.</param>
        /// <returns>The series sum.</returns>
        /// <exception cref="InvalidOperationException">The series does not converge within its iteration limit.</exception>
        /// <remarks>
        /// With λ = μ²/2, the weights p_j = e^(−λ) λ^j / j! and q_j = μ e^(−λ) λ^j / (√2 Γ(j + 3/2)) are evaluated in log
        /// space at the modal index k = ⌊λ⌋, as are I_y(k + ½, ν/2) and I_y(k + 1, ν/2). The other terms follow from
        /// I_y(a + 1, b) = I_y(a, b) − Γ(a + b) / (Γ(a + 1) Γ(b)) y^a (1 − y)^b, applied forward and backward from k
        /// (Benton and Krishnamoorthy 2003). Each direction stops when a geometric bound on its remaining terms is below
        /// machine precision relative to the magnitude of the result.
        /// </remarks>
        private static double PoissonMixtureSeries(double t, double df, double delta, bool upper, double offset)
        {
            // y = t²/(t² + ν) and 1 − y = ν/(t² + ν), without overflow.
            double y, yc;
            if (t * t > df)
            {
                double ratio = df / (t * t);
                y = 1d / (1d + ratio);
                yc = ratio / (1d + ratio);
            }
            else
            {
                double ratio = t * t / df;
                y = ratio / (1d + ratio);
                yc = 1d / (1d + ratio);
            }
            double b = df / 2d;
            double lambda = delta * delta / 2d;
            // Limits where the incomplete beta functions are constant: Σ p_j = 1 and ½ Σ q_j = Φ(μ) − ½.
            if (y == 0d || yc == 0d)
            {
                bool full = (y == 0d) == upper;
                return full ? Normal.StandardCDF(delta) : 0d;
            }
            const double epsilon = 1E-16;
            double logLambda = Math.Log(lambda);
            double logY = Math.Log(y);
            double logYc = Math.Log(yc);
            double logGammaB = Gamma.LogGamma(b);
            double sign = delta < 0d ? -1d : 1d;
            double logAbsDelta = Math.Log(Math.Abs(delta));
            // Log of the recurrence term Γ(a + b) / (Γ(a + 1) Γ(b)) y^a (1 − y)^b.
            double logRecurrenceTerm(double a) => Gamma.LogGamma(a + b) - Gamma.LogGamma(a + 1d) - logGammaB + a * logY + b * logYc;
            // Each tail is evaluated directly on the side of the beta mean where it is the smaller value.
            double incompleteBeta(double a)
            {
                bool lowerIsSmall = y <= a / (a + b);
                if (upper) return lowerIsSmall ? 1d - Beta.Incomplete(a, b, y) : Beta.Incomplete(b, a, yc);
                return lowerIsSmall ? Beta.Incomplete(a, b, y) : 1d - Beta.Incomplete(b, a, yc);
            }
            // Applies one recurrence step; a subtraction that would cancel more than one bit is replaced by direct
            // evaluation at the new index, which keeps the decreasing direction of the recurrence stable.
            double recurrenceStep(double current, double change, double a)
            {
                if (change < 0d && -change > 0.5d * current) return incompleteBeta(a);
                return current + change;
            }
            if (lambda == 0d) return 0.5d * incompleteBeta(0.5d);

            // Start at the modal Poisson index. When the lower-tail incomplete beta values underflow there, start at
            // the largest index whose recurrence term is representable, so the terms below it are not lost.
            int k = (int)Math.Min(Math.Floor(lambda), int.MaxValue - 1);
            if (!upper && incompleteBeta(k + 0.5d) < 1E-280)
            {
                int low = 0, high = k;
                while (high - low > 1)
                {
                    int middle = low + (high - low) / 2;
                    if (logRecurrenceTerm(middle + 0.5d) >= -600d) low = middle; else high = middle;
                }
                k = low;
            }
            int maxTerms = 1000 + (int)Math.Min(100d * Math.Sqrt(lambda), 1E8);

            // Weights and incomplete beta values at the starting index.
            double p0 = Math.Exp(-lambda + k * logLambda - Gamma.LogGamma(k + 1d));
            double q0 = sign * Math.Exp(-lambda + k * logLambda + logAbsDelta - 0.5d * Math.Log(2d) - Gamma.LogGamma(k + 1.5d));
            double aP = k + 0.5d, aQ = k + 1d;
            double iP0 = incompleteBeta(aP);
            double iQ0 = incompleteBeta(aQ);
            double gP0 = Math.Exp(logRecurrenceTerm(aP));
            double gQ0 = Math.Exp(logRecurrenceTerm(aQ));
            double direction = upper ? 1d : -1d;
            double sum = p0 * iP0 + q0 * iQ0;

            // Forward from the starting index.
            double p = p0, q = q0, iP = iP0, iQ = iQ0, gP = gP0, gQ = gQ0;
            for (int j = k; ; j++)
            {
                if (j - k > maxTerms) throw new InvalidOperationException("The noncentral t series did not converge.");
                double nextP = p * lambda / (j + 1d);
                double nextQ = q * lambda / (j + 1.5d);
                iP = recurrenceStep(iP, direction * gP, j + 1.5d);
                iQ = recurrenceStep(iQ, direction * gQ, j + 2d);
                gP *= y * (j + 0.5d + b) / (j + 1.5d);
                gQ *= y * (j + 1d + b) / (j + 2d);
                p = nextP;
                q = nextQ;
                double term = p * iP + q * iQ;
                sum += term;
                // Beyond the mode the weights fall at least geometrically with ratio λ/(j + 2). The lower-tail
                // incomplete beta values decrease with j, and Σ p_j ≤ 1 and Σ |q_j| ≤ 1, so the remaining lower-tail
                // terms are also bounded by the current incomplete beta values. The upper-tail values are bounded by one.
                double rate = lambda / (j + 2d);
                double bound = upper ? double.PositiveInfinity : iP + iQ;
                if (rate < 1d) bound = Math.Min(bound, (p + Math.Abs(q)) * (upper ? 1d : Math.Max(iP, iQ)) * rate / (1d - rate));
                if (bound <= epsilon * (Math.Abs(sum) + offset)) break;
            }

            // Backward from the starting index.
            p = p0; q = q0; iP = iP0; iQ = iQ0;
            gP = gP0 * aP / ((aP - 1d + b) * y);
            gQ = gQ0 * aQ / ((aQ - 1d + b) * y);
            for (int j = k - 1; j >= 0; j--)
            {
                if (k - j > maxTerms) throw new InvalidOperationException("The noncentral t series did not converge.");
                p *= (j + 1d) / lambda;
                q *= (j + 1.5d) / lambda;
                iP = recurrenceStep(iP, -direction * gP, j + 0.5d);
                iQ = recurrenceStep(iQ, -direction * gQ, j + 1d);
                gP *= (j + 0.5d) / ((j - 0.5d + b) * y);
                gQ *= (j + 1d) / ((j + b) * y);
                double term = p * iP + q * iQ;
                sum += term;
                // Below the mode the weights fall at least geometrically with ratio j/λ; the upper-tail incomplete
                // beta values decrease as j decreases and the lower-tail values are bounded by one.
                double rate = j / lambda;
                double bound = (p + Math.Abs(q)) * (upper ? Math.Max(iP, iQ) : 1d) * rate / (1d - rate);
                if (bound <= epsilon * (Math.Abs(sum) + offset)) break;
            }
            return 0.5d * sum;
        }

        /// <summary>
        /// Returns P(T ≤ t) for t &lt; 0 and a positive noncentrality from the positive integral representation.
        /// </summary>
        /// <param name="t">The negative, finite point.</param>
        /// <param name="df">The degrees of freedom ν.</param>
        /// <param name="delta">The positive noncentrality parameter μ.</param>
        /// <returns>The lower-tail probability at <paramref name="t"/>.</returns>
        /// <remarks>
        /// With T = (Z + μ)/S, Z standard Normal and S = √(χ²_ν/ν), P(T ≤ t) = P(Z + μ ≤ 0, S ≤ (Z + μ)/t), so
        /// F(t) = ∫₀^∞ φ(μ + u) P(ν/2, νu²/(2t²)) du. The log integrand h(u) is concave with h'' ≤ −1, and its
        /// maximizer lies in [0, u₀] with u₀ = 2ν/(μ + √(μ² + 4ν)), because d ln P/du ≤ ν/u. The integral of
        /// exp(h(u) − h(û)) is taken over û ± 10, beyond which h is below its maximum by more than 45.
        /// </remarks>
        private static double NegativeTailIntegral(double t, double df, double delta)
        {
            double a = df / 2d;
            double scale = df / (2d * t * t);
            if (double.IsPositiveInfinity(scale)) return Normal.StandardCDF(-delta);
            double logIntegrand(double u)
            {
                double value = -0.5d * (delta + u) * (delta + u) + DistributionNumerics.GammaLogCDF(a, scale * u * u);
                // Keep the optimizer arithmetic finite where the incomplete gamma function underflows.
                return Math.Max(value, -1E300);
            }
            double upperBound = 2d * df / (delta + Math.Sqrt(delta * delta + 4d * df));
            var search = new BrentSearch(logIntegrand, 0d, upperBound)
            {
                ReportFailure = false,
                ComputeHessian = false,
                RecordTraces = false
            };
            search.Maximize();
            double uHat = search.BestParameterSet.Values[0];
            double hHat = logIntegrand(uHat);
            if (hHat <= -1E300) return 0d;
            // P(ν/2, νu²/(2t²)) rises from zero to one around u = |t| with width about |t|/√(2ν), which can be far
            // narrower than the window, so the window is split where the rise is complete.
            double lower = Math.Max(0d, uHat - 10d), upper = uHat + 10d;
            double split = Tools.Clamp(Math.Abs(t) * (1d + 12d / Math.Sqrt(df)), lower, upper);
            double integral = Integrate(u => Math.Exp(logIntegrand(u) - hHat), lower, split)
                + Integrate(u => Math.Exp(logIntegrand(u) - hHat), split, upper);
            return Math.Exp(hHat - Tools.LogSqrt2PI + Math.Log(integral));
        }

        /// <summary>
        /// The inverse of the non-central t distribution
        /// </summary>
        /// <param name="p">Probability strictly between 0 and 1.</param>
        /// <param name="df">The degrees of freedom.</param>
        /// <param name="delta">The noncentrality parameter.</param>
        /// <returns>The quantile associated with <paramref name="p"/>.</returns>
        /// <remarks>
        /// Starting from the Jennett and Welch estimate, the bracket is widened by doubling steps until the CDF
        /// straddles <paramref name="p"/>, and the root is then refined with Brent's method to a relative tolerance of
        /// about 2ε in x. Above the median the root of P(T &gt; x) = 1 − p is found instead.
        /// </remarks>
        private static double InverseCdf(double p, double df, double delta)
        {
            double guess = InitialInverseCdfEstimate(p, df, delta);
            if (!Tools.IsFinite(guess)) guess = delta;
            double lower = guess, upper = guess;
            double step = Math.Max(1d, 0.1d * Math.Abs(guess));
            if (LowerTail(guess, df, delta) < p)
            {
                upper = guess + step;
                while (LowerTail(upper, df, delta) < p && !double.IsInfinity(upper))
                {
                    lower = upper;
                    step *= 2d;
                    upper += step;
                }
            }
            else
            {
                lower = guess - step;
                while (LowerTail(lower, df, delta) > p && !double.IsInfinity(lower))
                {
                    upper = lower;
                    step *= 2d;
                    lower -= step;
                }
            }
            if (double.IsInfinity(lower) || double.IsInfinity(upper)) return double.IsInfinity(lower) ? lower : upper;
            // Above the median the upper tail is matched instead, since 1 − p is exact there and the upper tail keeps
            // its relative precision.
            if (p > 0.5d) return Brent.Solve(x => (1d - p) - LowerTail(-x, df, -delta), lower, upper, 1E-300, 1000, false);
            return Brent.Solve(x => LowerTail(x, df, delta) - p, lower, upper, 1E-300, 1000, false);
        }

        /// <summary>
        /// Computes the Johnson and Kotz starting estimate for noncentral-t inversion.
        /// </summary>
        /// <param name="probability">The cumulative probability.</param>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        /// <param name="noncentrality">The noncentrality parameter.</param>
        /// <returns>An initial quantile estimate for the root search.</returns>
        /// <remarks>Uses the Jennett and Welch approximation and falls back to an offset central-t quantile when necessary.</remarks>
        private static double InitialInverseCdfEstimate(double probability, double degreesOfFreedom, double noncentrality)
        {
            // Approximates percentage points of the non-central t distribution.
            // Source: Johnson & Kotz, Continuous Univariate Distributions, Volume 2.
            double z = Normal.StandardZ(probability);

            // Jennett & Welch approximation, formula (14.1).
            // Intended for the large noncentralities common in tolerance-interval calculations.
            double gammaRatio = Math.Sqrt(2d / degreesOfFreedom) / Gamma.RatioHalf(degreesOfFreedom / 2d);
            double zSquared = z * z;
            double gammaRatioSquared = gammaRatio * gammaRatio;
            double denominator = gammaRatioSquared - zSquared * (1d - gammaRatioSquared);
            double discriminant = gammaRatioSquared + (1d - gammaRatioSquared) * (noncentrality * noncentrality - zSquared);
            if (discriminant > 0d && Math.Abs(denominator) > 1e-12)
            {
                return (noncentrality * gammaRatio + z * Math.Sqrt(discriminant)) / denominator;
            }

            // Fall back to an offset central-t quantile.
            var studentT = new StudentT(degreesOfFreedom);
            return studentT.InverseCDF(probability) + noncentrality;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new NoncentralT(DegreesOfFreedom, Noncentrality);
        }

    }
}