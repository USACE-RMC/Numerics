using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Frank copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class FrankCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Frank copula with a dependency θ = 2.
        /// </summary>
        public FrankCopula()
        {
            Theta = 2d;
        }

        /// <summary>
        /// Constructs a Frank copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public FrankCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Frank copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public FrankCopula(double theta, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Frank; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Frank"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "F"; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException)
        {
            if (double.IsNaN(parameter) || double.IsInfinity(parameter))
            {
                var exception = new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter must be finite.");
                if (throwException) throw exception;
                return exception;
            }
            if (parameter < ThetaMinimum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be greater than or equal to " + ThetaMinimum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be greater than or equal to " + ThetaMinimum.ToString() + ".");
            }
            if (parameter > ThetaMaximum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be less than or equal to " + ThetaMaximum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The dependency parameter θ (theta) must be less than or equal to " + ThetaMaximum.ToString() + ".");
            }
            return null;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ(t) = −ln(R) with R = Expm1(−θt)/Expm1(−θ) where R ≤ 1/2, that is for t ≤ −Log1p(Expm1(−θ)/2)/θ, and
        /// otherwise −Log1p(r) with r = R − 1 = (e^(−θt) − e^(−θ))/Expm1(−θ). Where θ(1 − t) is within 1/2 of 0 the
        /// difference is formed as −e^(−θt)·Expm1(−θ(1 − t)), so it keeps its relative accuracy near t = 1 and for small
        /// |θ|; elsewhere the two exponentials differ by a factor of at least e^(1/2) and, while |θt| ≤ 64, are subtracted
        /// directly. The θ-dependent constants are computed once per value of θ. θ = 0 returns the independence generator
        /// −ln t, its limit.
        /// </remarks>
        public override double Generator(double t)
        {
            double theta = Theta;
            if (theta == 0d) return -Math.Log(t);
            var terms = GetThetaTerms(theta);
            double x = -theta * t;
            if (t <= terms.GeneratorCut) return -Math.Log(Tools.Expm1(x) * terms.InverseExpMinusOne);
            // Near t = 1 the ratio is close to 1: take ln(1 + r) from its difference r.
            double y = -theta * (1d - t);
            double exp = Math.Exp(x);
            // The direct difference amplifies the rounding of θt by at most e^(1/2)/(e^(1/2) − 1), about 2.5, so it is
            // used while |θt| ≤ 64, where that rounding is below 1E-14 of e^(−θt).
            if (Math.Abs(y) > 0.5d && Math.Abs(x) <= 64d) return -Tools.Log1p((exp - terms.Exp) * terms.InverseExpMinusOne);
            return -Tools.Log1p(-exp * Tools.Expm1(y) * terms.InverseExpMinusOne);
        }

        /// <summary>
        /// The θ-dependent constants of the generator members and the θ they belong to, replaced as one object when θ
        /// changes.
        /// </summary>
        private sealed class ThetaTerms
        {
            /// <summary>Creates the cache entry.</summary>
            /// <param name="theta">The dependency parameter, not zero.</param>
            internal ThetaTerms(double theta)
            {
                Theta = theta;
                InverseTheta = 1d / theta;
                double expMinusOne = Tools.Expm1(-theta);
                ExpMinusOne = expMinusOne;
                InverseExpMinusOne = 1d / expMinusOne;
                Exp = Math.Exp(-theta);
                TwoMinusExp = 1d - expMinusOne;
                // R = Expm1(−θt)/Expm1(−θ) ≤ 1/2 exactly where t ≤ −Log1p(Expm1(−θ)/2)/θ.
                GeneratorCut = -Tools.Log1p(0.5d * expMinusOne) / theta;
                // r = e^(−s)·Expm1(−θ) < −1/2 exactly where s < ln(−2·Expm1(−θ)).
                InverseCut = -expMinusOne > 0.5d ? Math.Log(-2d * expMinusOne) : double.NegativeInfinity;
            }

            /// <summary>The dependency parameter.</summary>
            internal readonly double Theta;

            /// <summary>1/θ.</summary>
            internal readonly double InverseTheta;

            /// <summary>e^(−θ) − 1.</summary>
            internal readonly double ExpMinusOne;

            /// <summary>1/(e^(−θ) − 1).</summary>
            internal readonly double InverseExpMinusOne;

            /// <summary>e^(−θ).</summary>
            internal readonly double Exp;

            /// <summary>2 − e^(−θ).</summary>
            internal readonly double TwoMinusExp;

            /// <summary>The argument at and below which the generator takes the logarithm of the ratio directly.</summary>
            internal readonly double GeneratorCut;

            /// <summary>The generator value below which the inverse generator uses its half-angle form.</summary>
            internal readonly double InverseCut;
        }

        /// <summary>The cached θ-dependent constants.</summary>
        private ThetaTerms? _thetaTerms;

        /// <summary>
        /// Returns the θ-dependent constants, computed once per value of θ; the generator members use them on every call.
        /// </summary>
        /// <param name="theta">The current dependency parameter, not zero.</param>
        /// <returns>The constants for <paramref name="theta"/>.</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private ThetaTerms GetThetaTerms(double theta)
        {
            var terms = _thetaTerms;
            return terms != null && terms.Theta == theta ? terms : UpdateThetaTerms(theta);
        }

        /// <summary>Computes and caches the θ-dependent constants for a new value of θ.</summary>
        /// <param name="theta">The current dependency parameter, not zero.</param>
        /// <returns>The constants for <paramref name="theta"/>.</returns>
        private ThetaTerms UpdateThetaTerms(double theta)
        {
            var terms = new ThetaTerms(theta);
            _thetaTerms = terms;
            return terms;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ⁻¹(s) = −Log1p(r)/θ with r = e^(−s)·Expm1(−θ), which is e^(−s) at θ = 0. Where r &lt; −1/2, that is for
        /// s &lt; ln(−2·Expm1(−θ)) with θ > 0, 1 + r is formed from τ = tanh(s/2), with e^(−s) = (1 − τ)/(1 + τ), as
        /// (e^(−θ) + (2 − e^(−θ))·τ)/(1 + τ), a ratio of positive terms.
        /// </remarks>
        public override double GeneratorInverse(double t)
        {
            double theta = Theta;
            if (theta == 0d) return Math.Exp(-t);
            var terms = GetThetaTerms(theta);
            if (t < terms.InverseCut)
            {
                double tau = Math.Tanh(0.5d * t);
                return -Math.Log((terms.Exp + terms.TwoMinusExp * tau) / (1d + tau)) * terms.InverseTheta;
            }
            return -Tools.Log1p(Math.Exp(-t) * terms.ExpMinusOne) * terms.InverseTheta;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ′(t) = −θ/Expm1(θt), which is −1/t at θ = 0. Where |θt| ≤ 1/2 it is evaluated as θ/2 − (θ/2)·coth(θt/2),
        /// the same value, in which coth(θt/2) is at least 4 in magnitude, so the difference loses no accuracy; below
        /// |θt| = 1E-8 it is θ/2 − 1/t, its expansion to double precision.
        /// </remarks>
        public override double GeneratorPrime(double t)
        {
            double theta = Theta;
            double x = theta * t;
            if (Math.Abs(x) > 0.5d) return -theta / (Math.Exp(x) - 1d);
            // Below 1E-8, x/Expm1(x) = 1 − x/2 to double precision, and θt may be subnormal; θ = 0 returns −1/t here.
            if (Math.Abs(x) < 1E-8d) return 0.5d * theta - 1d / t;
            // 1/Expm1(x) = (coth(x/2) − 1)/2.
            double half = 0.5d * theta;
            return half - half / Math.Tanh(half * t);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ″(t) = θ²·e^(θt)/Expm1(θt)², evaluated for θt > 0 as θ²·e^(−θt)/Expm1(−θt)² so it does not overflow; it is
        /// 1/t² at θ = 0.
        /// </remarks>
        public override double GeneratorPrime2(double t)
        {
            double theta = Theta;
            if (theta == 0d) return 1d / (t * t);
            double x = theta * t;
            double ratio = x > 0d ? theta / Tools.Expm1(-x) : theta / Tools.Expm1(x);
            return ratio * ratio * Math.Exp(-Math.Abs(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The closed form Log1p(−θ/s)/θ of φ′(t) = s, which is −1/s at θ = 0. Where |θ/s| > 1/2, that is |s| &lt; 2|θ|,
        /// it is ln((s − θ)/s)/θ, whose rounding is small against the logarithm there.
        /// </remarks>
        public override double GeneratorPrimeInverse(double t)
        {
            double theta = Theta;
            if (Math.Abs(t) < 2d * Math.Abs(theta)) return Math.Log((t - theta) / t) / theta;
            if (theta == 0d) return -1d / t;
            return Tools.Log1p(-theta / t) / theta;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The density is <c>Math.Exp(LogPDF(u, v))</c>; see <see cref="LogPDF(double, double)"/>. At θ = 0 it is the
        /// independence density 1.
        /// </remarks>
        public override double PDF(double u, double v)
        {
            return Math.Exp(LogPDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The Frank density c(u,v) = θ(1 − e^(−θ))·e^(−θ(u+v)) / D², with
        /// D = (1 − e^(−θ)) − (1 − e^(−θu))(1 − e^(−θv)) (Nelsen, 2006), is evaluated in log space. For θ &lt; 0 it
        /// uses the symmetry c_θ(u, v) = c_(−θ)(u, 1 − v). For θ > 0, with m ≤ M the two arguments,
        /// D = e^(−θm)·[(1 − e^(−θ(1−m))) + e^(−θ(M−m))·(1 − e^(−θm))], a sum of two non-negative terms, so ln D has
        /// no cancellation and no exponential overflows for any θ.
        /// </para>
        /// <para>
        /// Each term keeps its relative accuracy as θ → 0, where ln θ, ln(1 − e^(−θ)) and −2·ln D cancel smoothly to 0;
        /// θ = 0 returns the independence limit 0.
        /// </para>
        /// </remarks>
        public override double LogPDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double theta = Theta, w = v;
            if (theta == 0d) return 0d;
            if (theta < 0d)
            {
                theta = -theta;
                w = 1d - v;
            }
            return Math.Log(theta) + Math.Log(-Tools.Expm1(-theta)) - theta * (u + w) - 2d * LogD(theta, u, w);
        }

        /// <summary>
        /// Returns ln D for D = (1 − e^(−θ)) − (1 − e^(−θu))(1 − e^(−θw)) and θ > 0.
        /// </summary>
        /// <param name="theta">The positive dependency parameter.</param>
        /// <param name="u">The first argument.</param>
        /// <param name="w">The second argument.</param>
        /// <returns>The log of D.</returns>
        /// <remarks>
        /// With m ≤ M the two arguments, D = e^(−θm)·[(1 − e^(−θ(1−m))) + e^(−θ(M−m))·(1 − e^(−θm))], a sum of two
        /// non-negative terms, so ln D has no cancellation and no exponential overflows.
        /// </remarks>
        private static double LogD(double theta, double u, double w)
        {
            double low = Math.Min(u, w), high = Math.Max(u, w);
            return -theta * low + Math.Log(-Tools.Expm1(-theta * (1d - low)) + Math.Exp(-theta * (high - low)) * -Tools.Expm1(-theta * low));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// C(u,v) = −Log1p(r)/θ with r = Expm1(−θu)·Expm1(−θv)/Expm1(−θ). For θ > 0, where r lies in (−1, 0], an r below
        /// −1/2 takes ln(1 + r) as ln D − ln(1 − e^(−θ)) with ln D from <see cref="LogD(double, double, double)"/>, since
        /// 1 + r = D/(1 − e^(−θ)). Both forms are smooth through θ = 0, and θ = 0 returns the independence limit uv.
        /// </remarks>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double theta = Theta;
            if (theta == 0d) return u * v;
            double r = Tools.Expm1(-theta * u) * Tools.Expm1(-theta * v) / Tools.Expm1(-theta);
            if (theta < 0d || r > -0.5d) return -Tools.Log1p(r) / theta;
            // 1 + r = D / (1 − e^(−θ)) is below one half: take its log from ln D rather than from 1 + r.
            return -(LogD(theta, u, v) - Math.Log(-Tools.Expm1(-theta))) / theta;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// h(v|u) = e^(−θu)·(e^(−θv) − 1) / ((e^(−θ) − 1) + (e^(−θu) − 1)(e^(−θv) − 1)). For θ > 0 the denominator is −D,
        /// so h = exp(−θu + ln(1 − e^(−θv)) − ln D) with ln D from <see cref="LogD(double, double, double)"/>; for θ &lt; 0
        /// every term is positive and is evaluated with Expm1. Both forms are smooth through θ = 0, and θ = 0 returns the
        /// independence limit v.
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (v == 0d) return 0d;
            if (v == 1d) return 1d;
            double theta = Theta;
            if (theta == 0d) return v;
            if (theta > 0d) return Math.Exp(-theta * u + Math.Log(-Tools.Expm1(-theta * v)) - LogD(theta, u, v));
            double eu = Tools.Expm1(-theta * u), ev = Tools.Expm1(-theta * v);
            return (1d + eu) * ev / (Tools.Expm1(-theta) + eu * ev);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the Frank closed-form conditional inversion v = −ln(1 + t·(e^(−θ) − 1)/(t + (1 − t)·e^(−θu)))/θ. For
        /// θ &lt; 0 it is evaluated as Log1p(t·Expm1(−θ)/(1 + (1 − t)·Expm1(−θu)))/(−θ), where every term is
        /// non-negative. For θ > 0 the Log1p argument r lies in (−1, 0]; below −1/2, ln(1 + r) is taken as the log of
        /// ((1 − t)·e^(−θu) + t·e^(−θ))/(t + (1 − t)·e^(−θu)), a ratio of positive terms. Both forms keep their relative
        /// accuracy for small v and are smooth through θ = 0; θ = 0 returns the independence limit t.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double theta = Theta;
            if (theta == 0d) return t;
            if (theta < 0d)
                return Tools.Log1p(t * Tools.Expm1(-theta) / (1d + (1d - t) * Tools.Expm1(-theta * u))) / -theta;
            // θ > 0: v = −ln(1 + t·(e^(−θ) − 1)/(t + (1 − t)·e^(−θu)))/θ, which keeps relative
            // accuracy for small v.
            double denominator = t + (1d - t) * Math.Exp(-theta * u);
            double r = t * Tools.Expm1(-theta) / denominator;
            if (r > -0.5d) return -Tools.Log1p(r) / theta;
            // 1 + r = ((1 − t)·e^(−θu) + t·e^(−θ))/denominator, a ratio of positive terms.
            return -(-theta * u + Math.Log((1d - t) + t * Math.Exp(-theta * (1d - u))) - Math.Log(denominator)) / theta;
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 0.
        /// The Frank copula has no tail dependence.
        /// </summary>
        public override double UpperTailDependence => 0.0;

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 0.
        /// The Frank copula has no tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new FrankCopula(Theta, CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

        /// <summary>
        /// Returns Kendall's τ (tau) implied by the Frank copula dependency parameter θ (theta).
        /// </summary>
        /// <param name="theta">The dependency parameter, θ. Must be finite and non-zero.</param>
        /// <returns>Kendall's τ for the given θ. The value is negative for θ &lt; 0 and positive for θ &gt; 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when θ is zero or is not finite. θ = 0 is the independence limit, at which the relation is
        /// indeterminate; it is rejected rather than returned as a silent NaN.
        /// </exception>
        /// <remarks>
        /// <para>
        /// The Frank copula relates θ to Kendall's τ through the order-1 Debye function D₁:
        /// </para>
        /// <code>
        ///     τ(θ) = 1 - (4 / θ) * (1 - D₁(θ))
        /// </code>
        /// <para>
        /// The relation is odd in θ, spans τ in (-1, 1) as θ ranges over the real line, and has the removable
        /// independence limit τ = 0 at θ = 0, where this expression is indeterminate. The negative branch
        /// depends on the reflection D₁(-y) = D₁(y) + y/2; omitting it makes τ diverge instead of approaching
        /// -1.
        /// </para>
        /// <para>
        /// <b> Accuracy. </b> Agrees with pyvinecopulib 0.7.6 reference values to 3.4E-16.
        /// Measured against mpmath at 60 decimal digits on a grid over 0.1 &#8804; |θ| &#8804; 100, the worst
        /// absolute error is 1.8E-15, at θ = 0.2. Below |θ| of about 0.01 the leading terms of the expression
        /// cancel against each other and the absolute error grows, reaching 2.0E-13 at θ = 0.001, which is the
        /// floor of the fitting bracket and where τ itself is only 1.1E-4. The cancellation is a property of
        /// this form of the relation rather than of the Debye evaluation.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// </para>
        /// <list type="bullet">
        /// <item><description>
        /// Genest, C. (1987). Frank's family of bivariate distributions. Biometrika, 74(3), 549-555.
        /// </description></item>
        /// <item><description>
        /// Nelsen, R. B. (2006). An Introduction to Copulas, 2nd ed., Table 4.1. Springer, New York.
        /// </description></item>
        /// </list>
        /// </remarks>
        internal static double KendallsTauFromTheta(double theta)
        {
            if (theta == 0d || double.IsNaN(theta) || double.IsInfinity(theta))
                throw new ArgumentOutOfRangeException(nameof(theta), "The dependency parameter θ (theta) must be finite and non-zero. θ = 0 is the independence limit, at which the Kendall's tau relation is indeterminate.");

            return 1d - 4d / theta * (1d - Debye.FunctionOrderOne(theta));
        }

        /// <summary>
        /// Estimates the dependency parameter using the method of moments.
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        /// <exception cref="ArgumentException">
        /// Thrown when Kendall's τ for the sample data lies outside the range the fitting bracket can reach.
        /// </exception>
        /// <remarks>
        /// Kendall's τ is estimated from the sample data and the closed-form Kendall's τ(θ) relation is
        /// inverted with Brent's method over the bracket returned by
        /// <see cref="ParameterConstraints(IList{double}, IList{double})"/>, which is [0.001, 100] for a
        /// positive τ and [-100, -0.001] for a non-positive one. That bracket reaches |τ| in
        /// about [1.1E-4, 0.96065797], so a |τ| above the upper end means the dependence is too strong to fit
        /// within the bracket and throws. The independence limit θ = 0 leaves the Frank generator and
        /// distribution functions indeterminate and is excluded from the bracket, so a τ of 0, or any τ too
        /// small for the bracket to straddle a root, is assigned the endpoint of its own bracket nearest
        /// independence rather than handed to the solver. A τ of exactly 0 takes the negative bracket, matching
        /// <see cref="ParameterConstraints(IList{double}, IList{double})"/>, and so returns θ = -0.001; at that
        /// magnitude the copula is indistinguishable from independence on either branch.
        /// </remarks>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);

            if (Math.Abs(tau) > KendallsTauFromTheta(100d))
                throw new ArgumentException("For the Frank copula, tau must be in ~= [-0.96065, 0.96065], the range attainable over the fitting bracket θ (theta) of [-100, 100]. The dependency in the data is too strong to use the Frank copula.");

            // θ = 0 is the independence limit and is excluded from the bracket, so a τ smaller than the bracket
            // can reach is assigned the endpoint nearest independence. The sign test matches the one used by
            // ParameterConstraints and by the bracket below, so a τ of exactly 0 resolves the same way in all
            // three places.
            double nearIndependence = tau > 0d ? 0.001d : -0.001d;
            if (Math.Abs(tau) <= Math.Abs(KendallsTauFromTheta(nearIndependence)))
            {
                Theta = nearIndependence;
                return;
            }

            double L = tau > 0 ? 0.001d : -100d;
            double U = tau > 0 ? 100d : -0.001d;
            Theta = Brent.Solve(t => KendallsTauFromTheta(t) - tau, L, U);
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);
            double L = tau > 0 ? 0.001d : -100d;
            double U = tau > 0 ? 100d : -0.001d;
            return new double[,] { { L, U } };
        }

    }
}