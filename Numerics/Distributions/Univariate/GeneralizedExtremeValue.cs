using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Generalized Extreme Value distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The shape κ follows the Hosking sign convention: a positive κ gives an upper bound at ξ + α/κ.
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class GeneralizedExtremeValue : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
    
        /// <summary>
        /// Constructs a Generalized Extreme Value with a location of 100, scale of 10, and shape of 0.
        /// </summary>
        public GeneralizedExtremeValue()
        {
            SetParameters(100d, 10d, 0d);
        }

        /// <summary>
        /// Constructs a Generalized Extreme Value (GEV) distribution with the given parameters ξ, α, and κ.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GeneralizedExtremeValue(double location, double scale, double shape)
        {
            SetParameters(location, scale, shape);
        }

        /// <summary>
        /// Creates the GEV distribution of the annual maximum of a Poisson number of independent
        /// generalized Pareto exceedances of a threshold (the peaks-over-threshold model).
        /// </summary>
        /// <param name="rate">The mean annual number of exceedances λ of the Poisson count. Must be positive and finite.</param>
        /// <param name="location">The GPD location ξ, the threshold. Must be finite.</param>
        /// <param name="scale">The GPD scale α. Must be positive and finite.</param>
        /// <param name="shape">The GPD shape κ, in the Hosking sign convention of <see cref="GeneralizedPareto"/>. Must be finite.</param>
        /// <returns>The GEV distribution with location ξ + α (1 − λ^(−κ))/κ, scale α λ^(−κ), and the same shape κ.</returns>
        /// <remarks>
        /// <para>
        /// If the number of exceedances in a year is Poisson with mean λ and each exceedance is
        /// GPD(ξ, α, κ), the annual maximum M satisfies P(M ≤ x) = exp(−λ [1 − κ (x − ξ)/α]^(1/κ)) for x
        /// above the threshold, which is exactly a GEV CDF with the same shape κ, scale α λ^(−κ), and
        /// location ξ + α (1 − λ^(−κ))/κ. At κ = 0 the exceedances are exponential and the annual
        /// maximum is Gumbel with location ξ + α ln λ and scale α. The relation holds above the
        /// threshold only; below it the annual maximum has an atom P(M ≤ ξ) = exp(−λ), which the GEV
        /// extends smoothly.
        /// </para>
        /// <para>
        /// The shape is used exactly: a nonzero κ, however small, keeps its value and is not collapsed
        /// to zero. With t = ln λ and h = −κ, the location is evaluated as ξ + α (expm1(h t)/h) through
        /// <see cref="Tools.Expm1(double)"/>, which stays accurate as κ approaches zero where
        /// (λ^h − 1)/h loses digits by cancellation, and the scale as α λ^h through
        /// <see cref="Math.Pow(double, double)"/>. At κ = 0 the location is ξ + α t and the scale is α.
        /// For λ ≥ 1 the location equals the GPD quantile at exceedance probability 1/λ, the level
        /// exceeded once per year on average.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// <list type="bullet">
        /// <item><description>Madsen, H., Rasmussen, P. F., and Rosbjerg, D. (1997). Comparison of annual maximum series and partial duration series methods for modeling extreme hydrologic events: 1. At-site modeling. Water Resources Research, 33(4), 747-757.</description></item>
        /// <item><description>Coles, S. (2001). An Introduction to Statistical Modeling of Extreme Values, Section 4.2. Springer, London.</description></item>
        /// </list>
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the rate or the scale is not positive and finite,
        /// the location or the shape is not finite, or the mapped GEV location or scale is outside the range of a double.</exception>
        public static GeneralizedExtremeValue FromPoissonGeneralizedPareto(double rate, double location, double scale, double shape)
        {
            if (!Tools.IsFinite(rate) || rate <= 0d)
                throw new ArgumentOutOfRangeException(nameof(rate), "The Poisson rate must be positive and finite.");
            if (!Tools.IsFinite(location))
                throw new ArgumentOutOfRangeException(nameof(location), "The location parameter ξ (Xi) must be a number.");
            if (!Tools.IsFinite(scale) || scale <= 0d)
                throw new ArgumentOutOfRangeException(nameof(scale), "The scale parameter α (alpha) must be positive.");
            if (!Tools.IsFinite(shape))
                throw new ArgumentOutOfRangeException(nameof(shape), "The shape parameter κ (kappa) must be a number.");

            double logRate = Math.Log(rate);
            double gevLocation, gevScale;
            if (shape == 0d)
            {
                gevLocation = location + scale * logRate;
                gevScale = scale;
            }
            else
            {
                double tail = -shape;
                gevLocation = location + scale * (Tools.Expm1(tail * logRate) / tail);
                gevScale = scale * Math.Pow(rate, tail);
            }
            if (!Tools.IsFinite(gevLocation) || !Tools.IsFinite(gevScale) || gevScale <= 0d)
                throw new ArgumentOutOfRangeException(nameof(rate), "The rate and shape give a GEV location or scale outside the range of a double.");
            return new GeneralizedExtremeValue(gevLocation, gevScale, shape);
        }

        private double _xi; // location
        private double _alpha; // scale
        private double _kappa; // shape

        /// <summary>
        /// Gets and sets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters(value, Alpha, Kappa, false) is null;
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
                _parametersValid = ValidateParameters(Xi, value, Kappa, false) is null;
                _alpha = value;
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
                _parametersValid = ValidateParameters(Xi, Alpha, value, false) is null;
                _kappa = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.GeneralizedExtremeValue; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Extreme Value"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GEV"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Location (ξ)";
                parmString[1, 0] = "Scale (α)";
                parmString[2, 0] = "Shape (κ)";
                parmString[0, 1] = Xi.ToString();
                parmString[1, 1] = Alpha.ToString();
                parmString[2, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ξ", "α", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xi), nameof(Alpha), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xi, Alpha, Kappa]; }
        }

        /// <inheritdoc/>
        /// <remarks>The mean exists for kappa &gt; -1. Log-Gamma arithmetic preserves finite scaled
        /// values even when a raw Gamma function overflows.</remarks>
        public override double Mean
        {
            get
            {
                if (Kappa <= -1) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return Xi + Alpha * SmallShapeStandardizedMoments()[0];
                double logarithm = Gamma.LogGamma(1 + Kappa);
                if (logarithm == 0) return Xi;
                double magnitude = logarithm > 0 ? logarithm + Tools.Log1mExp(-logarithm)
                    : Tools.Log1mExp(logarithm);
                return Xi - Math.Sign(Kappa) * Math.Sign(logarithm)
                    * Math.Exp(Math.Log(Alpha) + magnitude - Math.Log(Math.Abs(Kappa)));
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(.5); }
        }

        /// <inheritdoc/>
        /// <remarks>The mode is the finite upper endpoint when kappa is at least one; otherwise it is the interior stationary point.</remarks>
        public override double Mode
        {
            get
            {
                if (Kappa >= 1) return Maximum;
                double logarithm = Tools.Log1p(-Kappa);
                return Xi - Alpha * logarithm * DistributionNumerics.Exprel(Kappa * logarithm);
            }
        }

        /// <inheritdoc/>
        /// <remarks>The variance exists for kappa &gt; -1/2. Nonexistent moments return NaN.</remarks>
        public override double StandardDeviation
        {
            get
            {
                if (Kappa <= -.5) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return Alpha * SmallShapeStandardizedMoments()[1];
                if (Kappa == 1) return Alpha;
                double logVariance = LogPowerVariance(Kappa);
                return Math.Exp(Math.Log(Alpha) + .5 * logVariance - Math.Log(Math.Abs(Kappa)));
            }
        }

        /// <inheritdoc/>
        /// <remarks>The third moment exists for kappa &gt; -1/3; positive bounded shapes are not excluded.</remarks>
        public override double Skewness
        {
            get
            {
                if (Kappa <= -1d / 3d) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return SmallShapeStandardizedMoments()[2];
                if (Kappa == 1) return -2;
                double l1 = Gamma.LogGamma(1 + Kappa), l2 = Gamma.LogGamma(1 + 2 * Kappa), l3 = Gamma.LogGamma(1 + 3 * Kappa);
                if (double.IsPositiveInfinity(l3)) return double.NegativeInfinity;
                double largest = Math.Max(l3, Math.Max(l1 + l2, 3 * l1));
                double centered = Math.Exp(l3 - largest) - 3 * Math.Exp(l1 + l2 - largest) + 2 * Math.Exp(3 * l1 - largest);
                return centered == 0 ? 0 : -Math.Sign(Kappa) * Math.Sign(centered)
                    * Math.Exp(largest + Math.Log(Math.Abs(centered)) - 1.5 * LogPowerVariance(Kappa));
            }
        }

        /// <inheritdoc/>
        /// <remarks>Ordinary kurtosis requires kappa &gt; -1/4. Analytical normalized central moments are used without raw-moment overflow.</remarks>
        public override double Kurtosis
        {
            get
            {
                if (Kappa <= -.25) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return SmallShapeStandardizedMoments()[3];
                if (Kappa == 1) return 9;
                double l1 = Gamma.LogGamma(1 + Kappa), l2 = Gamma.LogGamma(1 + 2 * Kappa);
                double l3 = Gamma.LogGamma(1 + 3 * Kappa), l4 = Gamma.LogGamma(1 + 4 * Kappa);
                if (double.IsPositiveInfinity(l4)) return double.PositiveInfinity;
                double largest = Math.Max(Math.Max(l4, l1 + l3), Math.Max(2 * l1 + l2, 4 * l1));
                double centered = Math.Exp(l4 - largest) - 4 * Math.Exp(l1 + l3 - largest)
                    + 6 * Math.Exp(2 * l1 + l2 - largest) - 3 * Math.Exp(4 * l1 - largest);
                return Math.Exp(largest + Math.Log(centered) - 2 * LogPowerVariance(Kappa));
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                if (Kappa >= 0)
                {
                    return double.NegativeInfinity;
                }
                else
                {
                    return Xi + Alpha / Kappa;
                }
            }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                if (Kappa <= 0)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return Xi + Alpha / Kappa;
                }
            }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
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
            DistributionNumerics.ValidateSample(sample, 4);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(DirectMethodOfMoments(Statistics.ProductMoments(sample)));
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(Statistics.LinearMoments(sample)));
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
            var newDistribution = new GeneralizedExtremeValue(Xi, Alpha, Kappa);
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
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public void SetParameters(double location, double scale, double shape)
        {
            _parametersValid = ValidateParameters(location, scale, shape, false) is null;
            _xi = location;
            _alpha = scale;
            _kappa = shape;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(parameters), "Exactly three parameters are required.");
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double location, double scale, double shape, bool throwException)
        {
            if (double.IsNaN(location) || double.IsInfinity(location))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            if (double.IsNaN(shape) || double.IsInfinity(shape))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "Exactly three parameters are required.");
                if (throwException) throw exception;
                return exception;
            }
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <summary>
        /// Gets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="moments">The array of sample moments.</param>
        public double[] DirectMethodOfMoments(IList<double> moments)
        {
            // Solve for kappa
            double k = SolveForKappa(moments[2]);
            double a;
            double x;
            if (Math.Abs(k) <= NearZero)
            {
                a = Math.Sqrt(6d) / Math.PI * moments[1];
                x = moments[0] - a * Tools.Euler;
            }
            else
            {
                double U1 = Gamma.Function(1d + k);
                double U2 = Gamma.Function(1d + 2d * k);
                a = Math.Sqrt(moments[1] * moments[1] * k * k / (U2 - Math.Pow(U1, 2d)));
                x = moments[0] - a / k * (1d - U1);
            }

            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            // Solve for kappa
            double k = SolveForKappa(moments[2]);
            double a;
            double x;
            if (Math.Abs(k) <= NearZero)
            {
                a = Math.Sqrt(6d) / Math.PI * moments[1];
                x = moments[0] - a * Tools.Euler;
            }
            else
            {
                double U1 = Gamma.Function(1d + k);
                double U2 = Gamma.Function(1d + 2d * k);
                a = Math.Sqrt(moments[1] * moments[1] * k * k / (U2 - Math.Pow(U1, 2d)));
                x = moments[0] - a / k * (1d - U1);
            }

            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new GeneralizedExtremeValue();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <summary>
        /// Solve for the shape parameter κ (kappa) given the skewness coefficient.
        /// </summary>
        /// <param name="skew">The skewness coefficient</param>
        /// <returns>
        /// Kappa
        /// </returns>
        /// <remarks>
        /// <para>
        /// For 0 &lt;= skew &lt;= 3.7, regression polynomials give kappa (zero at skew = 1.14); on 1.14 &lt; skew &lt;= 3.7
        /// the polynomial is within 2E-4 of the exact root. Elsewhere the skewness equation is solved by the Brent
        /// method, using that skewness decreases monotonically in kappa for kappa &gt; -1/3: kappa in (-1/3, 0) for
        /// 3.7 &lt; skew &lt; 10, (-1/3, 1] for -2 &lt;= skew &lt; 0, and [1, 100] for skew &lt; -2. A skew of 10 or more
        /// returns NaN.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentException">The skew is below the skewness at kappa = 100 (about -1.4E52).</exception>
        public double SolveForKappa(double skew)
        {
            if (skew > 3.7d && skew < 10d)
            {
                // Extreme Value Type II, beyond the range where the regression polynomial is accurate.
                // The third moment exists only for kappa > -1/3.
                return Brent.Solve((x) => new GeneralizedExtremeValue(0d, 1d, x).Skewness - skew, -1d / 3d + 1E-12, 0d);
            }
            else if (skew > 1.14d && skew <= 3.7d)
            {
                // Extreme Value Type II
                return 0.2858221d - 0.357983d * skew + 0.116659d * Math.Pow(skew, 2d) - 0.022725d * Math.Pow(skew, 3d) + 0.002604d * Math.Pow(skew, 4d) - 0.000161d * Math.Pow(skew, 5d) + 0.000004d * Math.Pow(skew, 6d);
            }
            else if (skew == 1.14d)
            {
                // Extreme Value Type I
                return 0d;
            }
            else if (skew >= 0d && skew < 1.14d)
            {
                // Extreme Value Type III
                // This regression equation works for -2 < Cs < 1.14
                return 0.277648d - 0.322016d * skew + 0.060278d * Math.Pow(skew, 2d) + 0.016759d * Math.Pow(skew, 3d) - 0.005873d * Math.Pow(skew, 4d) - 0.00244d * Math.Pow(skew, 5d) - 0.00005d * Math.Pow(skew, 6d);
            }
            else if (skew < 0d && skew >= -2)
            {
                // For negative values of skew, there exists two possible values of kappa.
                // Either extreme value type II or III
                // Therefor it must be solved. The Brent method is used here. 
                return Brent.Solve((x) =>
                {
                    double U1 = Gamma.Lanczos(1d + x);
                    double U2 = Gamma.Lanczos(1d + 2d * x);
                    double U3 = Gamma.Lanczos(1d + 3d * x);
                    double k = Math.Sign(x) * (-U3 + 3d * U1 * U2 - 2d * Math.Pow(U1, 3d)) / Math.Pow(U2 - Math.Pow(U1, 2d), 3d / 2d);
                    return k - skew;
                }, -(1d / 3d), 1d);
            }
            else if (skew < -2)
            {
                // Extreme Value Type III with kappa > 1, where the skewness falls below -2 at kappa = 1
                // and decreases monotonically without bound.
                return Brent.Solve((x) => new GeneralizedExtremeValue(0d, 1d, x).Skewness - skew, 1d, 100d);
            }
            else
            {
                return double.NaN;
            }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Hosking's routine PELGEV from the LMOMENTS Fortran library, version 3.04, followed by one Newton refinement.
        /// The shape κ solves τ3 = 2(1 − 3^(−κ))/(1 − 2^(−κ)) − 3, with κ in the convention of this class (κ &gt; 0
        /// bounds the upper tail). PELGEV gives the rational-function approximation of Donaldson (1996) for
        /// −0.8 ≤ τ3 &lt; 1 and Newton-Raphson iteration on τ3 below −0.8 (at most 20 steps, stopping when the step is
        /// at most 1E-6 κ, started at 1 − ln(1 + τ3)/ln 2 when τ3 ≤ −0.97). One Newton step on the exact equation, in
        /// the expm1 form that stays accurate as κ tends to zero, then refines it. PELGEV's snap of |κ| &lt; 1E-5 to
        /// the Gumbel limit is not applied, since the refined shape is accurate there. The scale and location are
        /// α = λ2 κ/(Γ(1 + κ)(1 − 2^(−κ))) and ξ = λ1 − α(1 − Γ(1 + κ))/κ; for |κ| ≤ 0.05 they are evaluated from
        /// ln Γ(1 + κ) and the divided exponential, which are continuous through the Gumbel values α = λ2/ln 2 and
        /// ξ = λ1 − γα at κ = 0.
        /// </para>
        /// <para>
        /// Against the exact root of the τ3 equation (mpmath, 50 digits) on a grid of τ3 in [−0.999, 0.999], the
        /// absolute error in κ is at most 3E-14 for −0.8 ≤ τ3 &lt; 1 and the relative error 3E-14 for τ3 &lt; −0.8, from
        /// 3.2E-7 and 1.1E-11 for PELGEV alone. The constants ln 2, ln 3 and Euler's γ are used at full double
        /// precision instead of the eight-digit values of the Fortran source. The fourth L-moment ratio, if present,
        /// is not used.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// </para>
        /// <para>
        /// Hosking, J. R. M. (1996). Fortran routines for use with the method of L-moments, Version 3. Research
        /// Report RC20525, IBM Research Division, Yorktown Heights, N.Y.
        /// </para>
        /// <para>
        /// Donaldson, R. W. (1996). Calculating inverse CV, skew and PWM functions for Pearson-3, log-normal,
        /// extreme-value, and log-logistic distributions. Communications in Statistics - Simulation and
        /// Computation, 25(3), 741-747.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="moments"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Fewer than three L-moments are given, λ1 or λ2 is not finite,
        /// λ2 is not positive, or τ3 is not finite or |τ3| ≥ 1.</exception>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            if (moments == null) throw new ArgumentNullException(nameof(moments));
            if (moments.Count < 3 || !Tools.IsFinite(moments[0]) || !Tools.IsFinite(moments[1])
                || moments[1] <= 0 || !Tools.IsFinite(moments[2]) || Math.Abs(moments[2]) >= 1)
                throw new ArgumentOutOfRangeException(nameof(moments), "Finite L-moments require positive L-scale and absolute L-skewness below one.");
            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            double ln2 = Math.Log(2d), ln3 = Math.Log(3d);
            double g;
            if (T3 > 0d)
            {
                // Rational-function approximation for 0 < tau3 < 1.
                double z = 1d - T3;
                g = (-1d + z * (1.59921491d + z * (-0.48832213d + z * 0.01573152d))) / (1d + z * (-0.64363929d + z * 0.08985247d));
            }
            else
            {
                // Rational-function approximation for -0.8 <= tau3 <= 0.
                g = (0.28377530d + T3 * (-1.21096399d + T3 * (-2.50728214d + T3 * (-1.13455566d + T3 * -0.07138022d))))
                    / (1d + T3 * (2.06189696d + T3 * (1.31912239d + T3 * 0.25077104d)));
                if (T3 < -0.8d)
                {
                    // Newton-Raphson iteration on (1 - 3^-k)/(1 - 2^-k) = (tau3 + 3)/2 for tau3 < -0.8.
                    if (T3 <= -0.97d) g = 1d - Math.Log(1d + T3) / ln2;
                    double t0 = (T3 + 3d) * 0.5d;
                    for (int iteration = 0; iteration < 20; iteration++)
                    {
                        double x2 = Math.Pow(2d, -g), x3 = Math.Pow(3d, -g);
                        double xx2 = 1d - x2, xx3 = 1d - x3;
                        double derivative = (xx2 * x3 * ln3 - xx3 * x2 * ln2) / (xx2 * xx2);
                        double previous = g;
                        g -= (xx3 / xx2 - t0) / derivative;
                        if (Math.Abs(g - previous) <= 1E-6 * g) break;
                    }
                }
            }
            g = RefineLinearMomentShape(g, T3, ln2, ln3);
            if (Math.Abs(g) <= 0.05d)
            {
                // Small shapes: Gamma(1 + k) = exp(l) with l = ln Gamma(1 + k), (1 - Gamma(1 + k))/k = -(l/k) exprel(l)
                // and (1 - 2^-k)/k = ln 2 exprel(-k ln 2), so the Gumbel limit k = 0 is continuous.
                double logarithm = Gamma.LogGammaOnePlus(g);
                double standardMean = g == 0d ? Tools.Euler : -(logarithm / g) * DistributionNumerics.Exprel(logarithm);
                double smallAlpha = L2 / (Math.Exp(logarithm) * ln2 * DistributionNumerics.Exprel(-g * ln2));
                return [L1 - smallAlpha * standardMean, smallAlpha, g];
            }
            double gamma = Gamma.Function(1d + g);
            double alpha = L2 * g / (gamma * (1d - Math.Pow(2d, -g)));
            double xi = L1 - alpha * (1d - gamma) / g;
            return [xi, alpha, g];
        }

        /// <summary>
        /// One Newton step on the exact L-skewness equation of the GEV shape.
        /// </summary>
        /// <param name="k">The PELGEV estimate of the shape.</param>
        /// <param name="tau3">The L-skewness τ3, with |τ3| &lt; 1.</param>
        /// <param name="ln2">ln 2.</param>
        /// <param name="ln3">ln 3.</param>
        /// <returns>The refined shape.</returns>
        /// <remarks>
        /// The step solves f(k) = 2 a/b − 3 − τ3 = 0 with a = 1 − 3^(−k) = −expm1(−k ln 3) and b = 1 − 2^(−k) =
        /// −expm1(−k ln 2), so a/b keeps full relative precision as k tends to zero, and
        /// f′(k) = 2 (ln 3 · 3^(−k) b − ln 2 · 2^(−k) a)/b². At k = 0 the limits a/b = ln 3/ln 2 and
        /// f′(0) = 2 ln 3 (ln 2 − ln 3)/ln 2 are used. From the PELGEV error of at most 3.2E-7 the step is quadratically
        /// convergent, which leaves an absolute error of about 1E-14 in k.
        /// </remarks>
        private static double RefineLinearMomentShape(double k, double tau3, double ln2, double ln3)
        {
            double ratio, derivative;
            if (k == 0d)
            {
                ratio = ln3 / ln2;
                derivative = 2d * ln3 * (ln2 - ln3) / ln2;
            }
            else
            {
                double a = -Tools.Expm1(-k * ln3), b = -Tools.Expm1(-k * ln2);
                ratio = a / b;
                derivative = 2d * (ln3 * Math.Exp(-k * ln3) * b - ln2 * Math.Exp(-k * ln2) * a) / (b * b);
            }
            double refined = k - (2d * ratio - 3d - tau3) / derivative;
            return Tools.IsFinite(refined) ? refined : k;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The L-moments are λ1 = ξ + α (1 - Γ(1 + κ))/κ, λ2 = α Γ(1 + κ)(1 - 2^(-κ))/κ,
        /// τ3 = 2(1 - 3^(-κ))/(1 - 2^(-κ)) - 3 and τ4 = (5(1 - 4^(-κ)) - 10(1 - 3^(-κ)) + 6(1 - 2^(-κ)))/(1 - 2^(-κ))
        /// (Hosking, 1996, routine LMRGEV). Each divided difference is evaluated without cancellation: with
        /// l = ln Γ(1 + κ) from its rational approximation about 1, (1 - Γ(1 + κ))/κ = -(l/κ) (e^l - 1)/l, and
        /// (1 - j^(-κ))/κ = ln j (e^x - 1)/x with x = -κ ln j. The expressions are continuous through κ = 0, where they
        /// give the Gumbel values λ1 = ξ + γα (γ the Euler-Mascheroni constant), λ2 = α ln 2,
        /// τ3 = 2 ln 3/ln 2 - 3 and τ4 = 16 - 10 ln 3/ln 2.
        /// </para>
        /// <para>
        /// Against mpmath (50 digits) at ξ = 10 and α = 2, for κ = 0, ±1E-12 to ±0.9, 2 and 5, the relative errors of
        /// λ1 and λ2 and the absolute errors of τ3 and τ4 are below 1E-14.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">The parameters are invalid or κ ≤ -1.</exception>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            if (kappa <= -1.0d)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "L-moments can only be defined for kappa > -1.");
            double ln2 = Math.Log(2d), ln3 = Math.Log(3d), ln4 = 2d * ln2;
            double logarithm = Gamma.LogGammaOnePlus(kappa);
            double standardMean = kappa == 0d ? Tools.Euler : -(logarithm / kappa) * DistributionNumerics.Exprel(logarithm);
            double b2 = ln2 * DistributionNumerics.Exprel(-kappa * ln2);
            double b3 = ln3 * DistributionNumerics.Exprel(-kappa * ln3);
            double b4 = ln4 * DistributionNumerics.Exprel(-kappa * ln4);
            double L1 = xi + alpha * standardMean;
            double L2 = alpha * Math.Exp(logarithm) * b2;
            double T3 = 2d * b3 / b2 - 3d;
            double T4 = (5d * b4 - 10d * b3 + 6d * b2) / b2;
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        /// <remarks>Requires at least four finite, nonconstant observations. Initialization uses
        /// the existing linear-moment estimator and preserves the fitting shape bounds.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample or a finite feasible initialization is invalid.</exception>
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
            // initialVals = DirectMethodOfMoments(Statistics.ComputeProductMoments(sample))
            initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(sample));
            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[1])) + 1d));
            // Get bounds of shape
            lowerVals[2] = -10;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
            //
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            double normalization = DistributionNumerics.InitializationScale(sample);
            var normalized = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) normalized[i] = sample[i] / normalization;
            var initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(normalized));
            initialVals[0] *= normalization;
            initialVals[1] *= normalization;
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], initialVals[1], Statistics.Minimum(sample),
                Statistics.Maximum(sample), false, out lowerVals[0], out upperVals[0]);
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
            // Get bounds of shape
            lowerVals[2] = -10;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
            if (!Tools.IsFinite(initialVals[2]) || initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
            // 
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
                var GEV = new GeneralizedExtremeValue();
                GEV.SetParameters(x);
                return GEV.LogLikelihood(sample);
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
        /// The exact shape controls support. At the finite upper endpoint, shape one has density
        /// 1/alpha and shapes above one have a genuine integrable infinite density.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x < Minimum || x > Maximum || double.IsInfinity(x)) return double.NegativeInfinity;
            if (Kappa > 0 && x == Maximum)
                return Kappa < 1 ? double.NegativeInfinity : Kappa == 1 ? -Math.Log(Alpha) : double.PositiveInfinity;
            if (Kappa < 0 && x == Minimum) return double.NegativeInfinity;
            double y = TransformedValue(x);
            double lf = -(1d - Kappa) * y - Math.Exp(-y) - Math.Log(Alpha);
            return double.IsNaN(lf) ? double.NegativeInfinity : lf;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return Math.Exp(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0;
            return -Math.Exp(-TransformedValue(x));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => -Tools.Expm1(LogCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return 0;
            if (x >= Maximum) return double.NegativeInfinity;
            double y = TransformedValue(x);
            double exponential = Math.Exp(-y);
            return exponential == 0 ? -y : Tools.Log1mExp(-exponential);
        }

        /// <summary>Maps an interior observation to its Gumbel coordinate using the exact nonzero shape.</summary>
        /// <param name="x">The observation in physical coordinates.</param>
        /// <returns>The corresponding Gumbel coordinate.</returns>
        private double TransformedValue(double x)
        {
            return DistributionNumerics.HoskingShapeTransform(x, Xi, Alpha, Kappa);
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
                ValidateParameters(Xi, Alpha, Kappa, true);
            return InteriorQuantile(Xi, Alpha, Kappa, Math.Log(-Math.Log(probability)));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Closed form in survival space: the quantile at survival probability q is that at the
        /// Gumbel logarithm g = log(−log(1 − q)), evaluated as log(−log1p(−q)) through
        /// <see cref="Tools.Log1p(double)"/>, so 1 − q is never formed and a small q keeps its full
        /// relative precision. The quantile ξ + α (1 − e^(κ g))/κ is then evaluated by the same
        /// private kernel as <see cref="InverseCDF(double)"/>, including its extended-range handling.
        /// Survival probabilities below about 1.1E-16, where 1 − q rounds to one, still give finite
        /// heavy-tail quantiles instead of the infinite upper limit. A survival probability of zero returns
        /// <see cref="Maximum"/> and one returns <see cref="Minimum"/>.
        /// </para>
        /// <para>
        /// Against 200-bit reference values over nine shapes in [−0.5, 0.5] and survival
        /// probabilities from 0.89 down to 1E-300, the error is within 3 + |κ g| units in the last
        /// place (the rounding of g is amplified by the condition number |κ g| of a heavy tail; the
        /// largest is 65 units at κ = −0.5 and q = 1E-300), whereas <see cref="InverseCDF(double)"/>
        /// at 1 − q loses up to about 1.3E11 units for q down to 1E-12.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Hosking, J. R. M., and Wallis, J. R. (1997). Regional Frequency Analysis: An Approach Based on
        /// L-Moments, Appendix A. Cambridge University Press.
        /// </para>
        /// </remarks>
        public override double InverseCCDF(double probability)
        {
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Maximum;
            if (probability == 1.0d)
                return Minimum;
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            return InteriorQuantile(Xi, Alpha, Kappa, Math.Log(-Tools.Log1p(-probability)));
        }

        /// <summary>
        /// The quantile at an interior probability p, given the Gumbel logarithm log(−log p).
        /// </summary>
        /// <param name="xi">The location parameter ξ.</param>
        /// <param name="alpha">The scale parameter α.</param>
        /// <param name="kappa">The shape parameter κ.</param>
        /// <param name="logarithm">log(−log p) for a probability p strictly between zero and one.</param>
        /// <returns>The quantile ξ + α (1 − e^(κ log(−log p)))/κ, evaluated through the scaled divided exponential.</returns>
        /// <remarks>
        /// Shared by <see cref="InverseCDF(double)"/>, <see cref="InverseCCDF(double)"/> and
        /// <see cref="ParameterSetQuantiles(double[,], IList{double}, out bool[])"/>, so all evaluate exactly the same operations.
        /// </remarks>
        private static double InteriorQuantile(double xi, double alpha, double kappa, double logarithm)
        {
            double product = kappa * logarithm;
            double unitQuantile = double.IsNegativeInfinity(product) ? 1 / kappa : DistributionNumerics.ScaledExprelProduct(1, -logarithm, product);
            double displacement = double.IsNegativeInfinity(product) ? alpha / kappa : DistributionNumerics.ScaledExprelProduct(alpha, -logarithm, product);
            return double.IsInfinity(displacement) && Tools.IsFinite(unitQuantile)
                ? alpha * (xi / alpha + unitQuantile) : xi + displacement;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The values are bitwise those of the base member (the per-row loop <c>ValidateParameters</c> →
        /// <c>SetParameters</c> → <c>InverseCDF</c> on a clone): a row is valid exactly when
        /// <see cref="ValidateParameters(double, double, double, bool)"/> accepts it, and each interior probability
        /// is evaluated by the same private kernel <see cref="InverseCDF(double)"/> uses. The Gumbel logarithm
        /// log(−log p) of each probability is computed once instead of once per row, and the parameters are read
        /// from the row without being applied to the clone. Probabilities of exactly zero or one, and invalid
        /// probabilities, go through <see cref="InverseCDF(double)"/> on the clone with the row applied, so the
        /// support limits and the exception are those of the base member.
        /// </para>
        /// <para>
        /// The rows are evaluated serially in order, so the result does not depend on the machine or the thread
        /// count.
        /// </para>
        /// </remarks>
        public override double[,] ParameterSetQuantiles(double[,] parameterSets, IList<double> probabilities, out bool[] valid)
        {
            if (parameterSets == null) throw new ArgumentNullException(nameof(parameterSets));
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            int count = parameterSets.GetLength(0);
            int arity = NumberOfParameters;
            if (parameterSets.GetLength(1) != arity)
                throw new ArgumentException("Each parameter set must have exactly " + arity + " parameters.", nameof(parameterSets));

            int m = probabilities.Count;
            var logarithms = new double[m];
            var interior = new bool[m];
            for (int k = 0; k < m; k++)
            {
                double probability = probabilities[k];
                interior[k] = probability > 0.0d && probability < 1.0d;
                if (interior[k]) logarithms[k] = Math.Log(-Math.Log(probability));
            }

            var result = new double[count, m];
            valid = new bool[count];
            var probe = (GeneralizedExtremeValue)Clone();
            for (int i = 0; i < count; i++)
            {
                double xi = parameterSets[i, 0], alpha = parameterSets[i, 1], kappa = parameterSets[i, 2];
                valid[i] = ValidateParameters(xi, alpha, kappa, false) is null;
                if (!valid[i])
                {
                    for (int k = 0; k < m; k++)
                        result[i, k] = double.NaN;
                    continue;
                }
                bool applied = false;
                for (int k = 0; k < m; k++)
                {
                    if (interior[k])
                    {
                        result[i, k] = InteriorQuantile(xi, alpha, kappa, logarithms[k]);
                        continue;
                    }
                    if (!applied)
                    {
                        probe.SetParameters(xi, alpha, kappa);
                        applied = true;
                    }
                    result[i, k] = probe.InverseCDF(probabilities[k]);
                }
            }
            return result;
        }

        /// <summary>
        /// Gets the expected Fisher information matrix.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <returns>The full three-parameter expected information in location, scale and shape coordinates.</returns>
        /// <remarks>Regular information requires kappa &lt; 1/2. At zero shape the shape parameter
        /// remains estimated, so this is not the two-parameter Gumbel information.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Sample size, distribution parameters or regularity are invalid.</exception>
        /// <exception cref="InvalidOperationException">Information cannot be resolved or represented numerically.</exception>
        public Matrix ExpectedInformationMatrix(int sampleSize)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            var information = KappaExpectedInformation.ExpectedInformation(Kappa, 0, 3, out _, out _, out _);
            for (int i = 0; i < 3; i++)
            for (int j = i; j < 3; j++)
            {
                double value = information[i, j];
                if (value != 0)
                {
                    double logarithm = Math.Log(Math.Abs(value)) + Math.Log(sampleSize)
                        - (i < 2 ? Math.Log(Alpha) : 0) - (j < 2 ? Math.Log(Alpha) : 0);
                    value = Math.Sign(value) * Math.Exp(logarithm);
                    if (!Tools.IsFinite(value))
                        throw new InvalidOperationException("Expected information is outside the finite floating-point range.");
                }
                information[i, j] = information[j, i] = value;
            }
            return new Matrix(information);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedExtremeValue(Xi, Alpha, Kappa);
        }

        /// <inheritdoc/>
        /// <remarks>Full three-parameter MLE covariance requires positive sample size and kappa &lt; 1/2.
        /// The zero-shape limit still estimates shape; uncertainty regularity does not restrict the distribution's valid shape domain.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters, sample size or information regularity are invalid.</exception>
        /// <exception cref="NotImplementedException">The method is not maximum likelihood.</exception>
        /// <exception cref="InvalidOperationException">The information or covariance is numerically unresolved.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            return KappaExpectedInformation.ParameterCovariance(Alpha, Kappa, 0, sampleSize, 3);
        }

        /// <inheritdoc/>
        /// <remarks>Analytical exponential divided differences retain the exact nonzero shape and
        /// the full [1,-L,-alpha*L squared/2] gradient at zero shape, where L=log(-log(p)).</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            double logarithm = Math.Log(-Math.Log(probability));
            double product = Kappa * logarithm;
            var gradient = new double[]
            {
                1.0d, // location
                double.IsNegativeInfinity(product) ? 1 / Kappa : DistributionNumerics.ScaledExprelProduct(1, -logarithm, product), // scale
                double.IsNegativeInfinity(product) ? -(Alpha / Kappa) / Kappa : -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, logarithm, product) // shape
            };
            return gradient;
        }

        /// <inheritdoc/>
        /// <remarks>Uses finite physical common-coordinate gradients and unit-scale covariance,
        /// permitting scalar variance evaluation independently of the physical covariance matrix's range.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            var unit = new GeneralizedExtremeValue(0, 1, Kappa);
            double logarithm = Math.Log(-Math.Log(probability)), product = Kappa * logarithm;
            double scaleGradient = DistributionNumerics.ScaledExprelProduct(Alpha, -logarithm, product);
            double shapeGradient = -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, logarithm, product);
            return DistributionNumerics.ScaledQuantileVariance(unit.ParameterCovariance(sampleSize, estimationMethod),
                [Alpha, scaleGradient, shapeGradient]);
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

        /// <summary>Returns log Var(T^power), T unit exponential, without forming raw gamma moments.</summary>
        /// <param name="power">The real exponent, established by the caller to be greater than negative one-half.</param>
        /// <returns>The logarithm of the variance of a unit-exponential variate raised to <paramref name="power"/>, including positive infinity when the second moment overflows.</returns>
        /// <remarks>The caller establishes power &gt; -1/2 and uses a divided series near zero.</remarks>
        internal static double LogPowerVariance(double power)
        {
            double first = Gamma.LogGamma(1 + power), second = Gamma.LogGamma(1 + 2 * power);
            return double.IsPositiveInfinity(second) ? double.PositiveInfinity
                : second + Tools.Log1mExp(2 * first - second);
        }

        /// <summary>Evaluates analytical standardized central moments with log-Gamma divided differences near zero shape.</summary>
        /// <returns>The standardized mean shift, standard deviation, skewness, and kurtosis in that order.</returns>
        /// <remarks>Finite differences of log Gamma remove the cancelling powers before evaluation.
        /// Exponential polynomial identities then retain the second through fourth centered moments.
        /// The exact zero is the analytical Gumbel limit; every nonzero kappa remains in the series.</remarks>
        private double[] SmallShapeStandardizedMoments()
        {
            double k = Kappa;
            if (k == 0) return [Tools.Euler, Math.PI / Math.Sqrt(6), 1.1395470994046487, 5.4];
            double logarithm = Gamma.LogGammaOnePlus(k);
            double mean = -(logarithm / k) * DistributionNumerics.Exprel(logarithm);
            double a = NormalizedLogGammaDifference(k, 2);
            double b = NormalizedLogGammaDifference(k, 3);
            double c = NormalizedLogGammaDifference(k, 4);
            double a2 = k * k * a, b3 = k * k * k * b, c4 = k * k * k * k * c;
            double u = Tools.Expm1(a2), v = Tools.Expm1(b3);
            double variance = a * DistributionNumerics.Exprel(a2);
            double third = b * DistributionNumerics.Exprel(b3);
            double fourth = c * DistributionNumerics.Exprel(c4);
            double skew = -(k * variance * variance * (3 + u) + Math.Exp(3 * a2) * third)
                / (variance * Math.Sqrt(variance));
            double kurtosis = (variance * variance * (3 + u * (16 + u * (15 + u * (6 + u))))
                + 12 * k * third * Math.Exp(3 * a2) * a * DistributionNumerics.Exprel(3 * a2)
                + Math.Exp(6 * a2) * (k * k * third * third * (6 + v * (4 + v)) + Math.Exp(4 * b3) * fourth))
                / (variance * variance);
            return [mean, Math.Exp(logarithm) * Math.Sqrt(variance), skew, kurtosis];
        }

        /// <summary>Returns the order-r forward difference of log Gamma(1+j*kappa), divided by kappa to power r.</summary>
        /// <param name="kappa">The nonzero shape increment.</param>
        /// <param name="order">The forward-difference order, expected to be two, three, or four.</param>
        /// <returns>The normalized log-gamma forward difference.</returns>
        private static double NormalizedLogGammaDifference(double kappa, int order)
        {
            double sum = 0, power = 1;
            for (int n = order; n <= 32; n++)
            {
                double factor = order == 2 ? Math.Pow(2, n) - 2
                    : order == 3 ? Math.Pow(3, n) - 3 * Math.Pow(2, n) + 3
                    : Math.Pow(4, n) - 4 * Math.Pow(3, n) + 6 * Math.Pow(2, n) - 4;
                sum += (n % 2 == 0 ? 1 : -1) * Gamma.ZetaInteger(n) * factor * power / n;
                power *= kappa;
            }
            return sum;
        }

    }
}
