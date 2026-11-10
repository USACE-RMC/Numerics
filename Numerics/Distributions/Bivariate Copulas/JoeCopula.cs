using System;
using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;
using System.Collections.Generic;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Joe copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class JoeCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Joe copula with a dependency θ = 2.
        /// </summary>
        public JoeCopula()
        {
            Theta = 2d;
        }

        /// <summary>
        /// Constructs a Joe copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public JoeCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Joe copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public JoeCopula(double theta, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Joe; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Joe"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "J"; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return 1.0d; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ(t) = −ln(1 − (1 − t)^θ), evaluated as −Log1mExp(θ·Log1p(−t)), which keeps its relative accuracy near both
        /// ends of [0, 1]; it is the independence generator −ln t at θ = 1.
        /// </remarks>
        public override double Generator(double t)
        {
            return -Tools.Log1mExp(Theta * Tools.Log1p(-t));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ⁻¹(s) = 1 − (1 − e^(−s))^(1/θ) = −Expm1(y) with y = ln(1 − e^(−s))/θ. Below s = ln 2 the logarithm is
        /// ln(2·sinh(s/2)) − s/2, a sum of two terms of one sign, and above it Log1p(−e^(−s)). Where y ≥ −1/2 the
        /// result is −2τ/(1 − τ) with τ = tanh(y/2), which keeps its relative accuracy for small y (and −y below 1E-300 in
        /// magnitude), and otherwise 1 − e^y.
        /// </remarks>
        public override double GeneratorInverse(double t)
        {
            double logW = t < Ln2 ? Math.Log(2d * Math.Sinh(0.5d * t)) - 0.5d * t : Tools.Log1p(-Math.Exp(-t));
            double y = logW / Theta;
            if (y < -0.5d) return 1d - Math.Exp(y);
            // Below 1E-300 in magnitude 1 − e^y = −y, which keeps a subnormal y exact.
            if (y > -1E-300d) return -y;
            // 1 − e^y = −2·tanh(y/2)/(1 − tanh(y/2)).
            double tau = Math.Tanh(0.5d * y);
            return -2d * tau / (1d - tau);
        }

        /// <summary>The natural logarithm of 2.</summary>
        private const double Ln2 = 0.69314718055994530942d;

        /// <inheritdoc/>
        /// <remarks>
        /// φ′(t) = −θ(1 − t)^(θ−1)/(1 − (1 − t)^θ), with the power p = (1 − t)^(θ−1) formed from Log1p(−t). With
        /// x = θ·Log1p(−t), the denominator is 1 − p·(1 − t) where x &lt; −1/2, and −2τ/(1 − τ) with τ = tanh(x/2)
        /// otherwise, which keeps its relative accuracy for small t.
        /// </remarks>
        public override double GeneratorPrime(double t)
        {
            double theta = Theta;
            double logBar = Tools.Log1p(-t);
            double x = theta * logBar;
            double power = Math.Exp((theta - 1d) * logBar);
            double denominator;
            if (x < -0.5d)
            {
                denominator = 1d - power * (1d - t);
            }
            else
            {
                // 1 − e^x = −2·tanh(x/2)/(1 − tanh(x/2)).
                double tau = Math.Tanh(0.5d * x);
                denominator = -2d * tau / (1d - tau);
            }
            return -theta * power / denominator;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ″(t) = θ(1 − t)^(θ−2)·((θ − 1) + (1 − t)^θ)/(1 − (1 − t)^θ)², with the powers formed from Log1p(−t).
        /// </remarks>
        public override double GeneratorPrime2(double t)
        {
            double theta = Theta;
            double logBar = Tools.Log1p(-t);
            double denominator = -Tools.Expm1(theta * logBar);
            return theta * Math.Exp((theta - 2d) * logBar) * ((theta - 1d) + Math.Exp(theta * logBar)) / (denominator * denominator);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Solves φ′(t) = s for t. With x = −ln(1 − t) the equation is (θ − 1)·x + ln(1 − e^(−θx)) = ln(−θ/s); in
        /// z = ln x the left side is increasing, and a safeguarded Newton iteration from the small-x and large-x
        /// approximations converges to full precision in a few steps; t = −Expm1(−x). Each step evaluates the function
        /// and its slope from one pair of exponentials: with a = θx and q = e^(−a), ln(1 − q) is Log1p(−q) above
        /// a = ln 2 and, below it, q and 1 − q come from τ = tanh(a/2) as (1 − τ)/(1 + τ) and 2τ/(1 + τ). At θ = 1 the
        /// closed form −1/s is used. Slopes at or above φ′(1) = 0 return 1.
        /// </remarks>
        public override double GeneratorPrimeInverse(double t)
        {
            double theta = Theta;
            if (!(t < 0d)) return t >= 0d ? 1d : double.NaN;
            if (theta == 1d) return Math.Min(1d, -1d / t);
            double c = Math.Log(theta) - Math.Log(-t);
            double small = -1d / t, large = c / (theta - 1d);
            // The small-x root −1/s overshoots when x is large, and the large-x root c/(θ − 1) when x is small, so the
            // smaller positive one starts nearer.
            double start = Math.Log(large > 0d ? Math.Min(small, large) : small);
            double z = SolveIncreasing(w =>
            {
                double x = Math.Exp(w), a = theta * x, q, oneMinusQ, logOneMinusQ;
                if (a > Ln2)
                {
                    q = Math.Exp(-a);
                    oneMinusQ = 1d - q;
                    logOneMinusQ = Tools.Log1p(-q);
                }
                else if (a < 1E-300d)
                {
                    // 1 − e^(−a) = a to double precision, which keeps a subnormal a exact.
                    q = 1d;
                    oneMinusQ = a;
                    logOneMinusQ = Math.Log(a);
                }
                else
                {
                    // e^(−a) = (1 − τ)/(1 + τ) with τ = tanh(a/2), so 1 − e^(−a) = 2τ/(1 + τ) keeps its relative accuracy.
                    double tau = Math.Tanh(0.5d * a), scale = 1d / (1d + tau);
                    q = (1d - tau) * scale;
                    oneMinusQ = 2d * tau * scale;
                    logOneMinusQ = Math.Log(oneMinusQ);
                }
                // The slope in z is x·((θ − 1) + θ/Expm1(θx)) = (θ − 1)·x + a·q/(1 − q).
                return ((theta - 1d) * x + logOneMinusQ - c, (theta - 1d) * x + a * q / oneMinusQ);
            }, start);
            return -Tools.Expm1(-Math.Exp(z));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Inside the open unit square the density is <c>Math.Exp(LogPDF(u, v))</c>, so it shares the log-space
        /// evaluation of <see cref="LogPDF(double, double)"/>. On the boundary and outside the square the generic
        /// Archimedean density of the base class is returned.
        /// </remarks>
        public override double PDF(double u, double v)
        {
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.PDF(u, v);
            return Math.Exp(LogPDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The Joe density c(u,v) = A^(1/θ−2)·(1−u)^(θ−1)·(1−v)^(θ−1)·(θ − 1 + A), with
        /// A = (1−u)^θ + (1−v)^θ − (1−u)^θ(1−v)^θ (Joe, 1997), is evaluated in log space:
        /// ln(1−u) and ln(1−v) come from Log1p, and ln A = h + Log1p(e^(l−h)·(1 − e^h)), where h ≥ l are
        /// θ·ln(1−u) and θ·ln(1−v). The powers of 1−u and 1−v never form, so the log density stays finite
        /// for every θ ≥ 1 and every u and v strictly inside the unit square.
        /// </para>
        /// <para>
        /// On the boundary and outside the square the log of the generic Archimedean density is returned.
        /// </para>
        /// </remarks>
        public override double LogPDF(double u, double v)
        {
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.LogPDF(u, v);
            double theta = Theta;
            double logUBar = Tools.Log1p(-u), logVBar = Tools.Log1p(-v);
            double logA = LogA(theta, logUBar, logVBar);
            return (1d / theta - 2d) * logA + (theta - 1d) * (logUBar + logVBar) + Math.Log(theta - 1d + Math.Exp(logA));
        }

        /// <summary>
        /// Returns ln A for A = (1−u)^θ + (1−v)^θ − (1−u)^θ(1−v)^θ, given ln(1−u) and ln(1−v).
        /// </summary>
        /// <param name="theta">The dependency parameter.</param>
        /// <param name="logUBar">ln(1 − u).</param>
        /// <param name="logVBar">ln(1 − v).</param>
        /// <returns>The log of A.</returns>
        /// <remarks>
        /// With h ≥ l the two values θ·ln(1−u) and θ·ln(1−v), ln A = h + Log1p(e^(l−h)·(1 − e^h)), so the powers never form.
        /// </remarks>
        private static double LogA(double theta, double logUBar, double logVBar)
        {
            double high = Math.Max(theta * logUBar, theta * logVBar);
            double low = Math.Min(theta * logUBar, theta * logVBar);
            return high + Tools.Log1p(Math.Exp(low - high) * -Tools.Expm1(high));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The closed form h(v|u) = (1−u)^(θ−1)·[1 − (1−v)^θ]·A^(1/θ−1), evaluated in log space inside the unit square as
        /// exp(<see cref="LogConditionalCDF(double, double)"/>).
        /// Unlike the generic ratio φ′(u)/φ′(C(u,v)), it does not pass through the rounded C, which loses relative accuracy
        /// where C is close to 1. On the boundary the generic conditional is returned.
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (v == 0d) return 0d;
            if (v == 1d) return 1d;
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.ConditionalCDF(u, v);
            return Math.Exp(LogConditionalCDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Inside the unit square the conditional and the density share ln(1−u), ln(1−v) and ln(1 + r), computed once:
        /// ln A = ln a + ln(1 + r) gives the density (1/θ − 2)·ln A + (θ − 1)·(ln(1−u) + ln(1−v)) + ln(θ − 1 + A).
        /// </remarks>
        protected override double LogConditionalAndDensity(double u, double v, out double density)
        {
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.LogConditionalAndDensity(u, v, out density);
            double theta = Theta;
            double logUBar = Tools.Log1p(-u), logVBar = Tools.Log1p(-v);
            double logA = theta * logUBar, logB = theta * logVBar;
            double logR = logB - logA + Tools.Log1mExp(logA);
            double log1pR = logR > 30d ? logR + Tools.Log1p(Math.Exp(-logR)) : Tools.Log1p(Math.Exp(logR));
            double logBigA = logA + log1pR;
            density = Math.Exp((1d / theta - 2d) * logBigA + (theta - 1d) * (logUBar + logVBar) + Math.Log(theta - 1d + Math.Exp(logBigA)));
            return Tools.Log1mExp(logB) - (1d - 1d / theta) * log1pR;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Inside the unit square, with a = (1−u)^θ, b = (1−v)^θ and r = b(1 − a)/a, the conditional is
        /// h = (1 − b)·(1 + r)^(−(1 − 1/θ)), so ln h = Log1mExp(ln b) − (1 − 1/θ)·Log1p(r), evaluated from ln a and ln b;
        /// both terms are small where h is close to 1, so 1 − h keeps its relative accuracy. Elsewhere it is the log of
        /// the conditional.
        /// </remarks>
        protected override double LogConditionalCDF(double u, double v)
        {
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return Math.Log(ConditionalCDF(u, v));
            double theta = Theta;
            double logA = theta * Tools.Log1p(-u), logB = theta * Tools.Log1p(-v);
            // h = (1 − b)·(1 + r)^(−(1 − 1/θ)) with a = (1−u)^θ, b = (1−v)^θ and r = b(1 − a)/a, so ln h is a sum of two
            // terms that are both small where h is close to 1.
            double logR = logB - logA + Tools.Log1mExp(logA);
            double log1pR = logR > 30d ? logR + Tools.Log1p(Math.Exp(-logR)) : Tools.Log1p(Math.Exp(logR));
            return Tools.Log1mExp(logB) - (1d - 1d / theta) * log1pR;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The Joe conditional has no closed-form inverse, so the closed-form conditional
        /// <see cref="ConditionalCDF(double, double)"/> is inverted with a safeguarded Newton iteration whose derivative is
        /// the density, started at v = t (see <see cref="ArchimedeanCopula.InvertConditionalCDF(double, double, double)"/>),
        /// to full double precision. For 0 &lt; u &lt; 1, t ≤ 0 returns 0 and t ≥ 1 returns 1. The iteration evaluates this class's closed forms
        /// through <see cref="LogConditionalAndDensity(double, double, out double)"/>, so a subclass that overrides
        /// <see cref="ConditionalCDF(double, double)"/> or <see cref="PDF(double, double)"/> must override it as well.
        /// At u = 0 the conditional is h(x|0) = 1 − (1−x)^θ, so the inverse is 1 − (1−t)^(1/θ), evaluated as
        /// −Expm1(Log1p(−t)/θ) to keep its relative accuracy for small t. At u = 1, where h
        /// is not defined, the inverse is its limit: for θ > 1 and 0 &lt; t &lt; 1, h(x|u) tends to 0 for every
        /// 0 &lt; x &lt; 1, so the inverse tends to 1; at θ = 1, and for t = 0 or t = 1, it is t.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (u == 0d) return -Tools.Expm1(Tools.Log1p(-t) / Theta);
            if (u == 1d)
            {
                if (Theta == 1d || !(t > 0d && t < 1d)) return t;
                return 1d;
            }

            if (!(t > 0d && t < 1d)) return t >= 1d ? 1d : t <= 0d ? 0d : double.NaN;
            if (Theta == 1d) return t;
            return InvertConditionalCDF(u, t, t);
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 2 - 2^(1/θ).
        /// </summary>
        public override double UpperTailDependence
        {
            get
            {
                return 2.0 - Math.Pow(2.0, 1.0 / Theta);
            }
        }

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 0.
        /// The Joe copula has no lower tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new JoeCopula(Theta, CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

        /// <summary>
        /// The number of terms of the Kendall's τ series that are summed directly before the closed-form tail is added.
        /// </summary>
        private const int TauSeriesTerms = 1000;

        /// <summary>
        /// The highest inverse power of the summation index retained in the closed-form tail of the Kendall's τ series.
        /// </summary>
        private const int TauSeriesTailOrder = 6;

        /// <summary>
        /// Returns Kendall's τ (tau) implied by the Joe copula dependency parameter θ (theta).
        /// </summary>
        /// <param name="theta">The dependency parameter, θ. Must be finite and greater than or equal to 1.</param>
        /// <returns>Kendall's τ for the given θ. The value is 0 at θ = 1 and increases to 1 as θ → ∞.</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when θ is below <see cref="ThetaMinimum"/> or is not finite. The series is defined outside
        /// that domain but the value it returns has no meaning for this copula, so an inadmissible θ is
        /// rejected rather than returned as a plausible-looking number.
        /// </exception>
        /// <remarks>
        /// <para>
        /// The Joe copula has no closed-form relation between θ and Kendall's τ. The relation is the series
        /// </para>
        /// <code>
        ///     τ(θ) = 1 - 4 * Σ[k = 1 to ∞] 1 / ( k * (θ*k + 2) * (θ*(k - 1) + 2) )
        /// </code>
        /// <para>
        /// The factor θ*(k - 1) + 2 equals 2 at k = 1, so no admissible θ divides by zero. The Joe copula
        /// models positive dependence only, so τ is confined to [0, 1).
        /// </para>
        /// <para>
        /// <b> Truncation. </b> The terms decay like 1 / (θ² k³), so simply truncating the series after K terms
        /// leaves an error of about 1 / (2 θ² K²) in the sum, which is 5E-7 at θ = 1 and K = 1,000 and would
        /// still be near 1E-13 after two million terms. Truncating on the size of the last term is therefore not
        /// sufficient. The first K = 1,000 terms are instead summed directly, from the smallest term to the
        /// largest to limit accumulated rounding, and the remainder is added as a closed form: the term is
        /// expanded in inverse powers of k, giving 1 / (θ² k³) * Σ[j ≥ 0] (-1)^j h_j / k^j with
        /// h_j = Σ[i = 0 to j] p^i q^(j - i), p = 2/θ and q = (2 - θ)/θ, and each power is summed over k &gt; K
        /// with the Euler-Maclaurin form of the Hurwitz zeta function. Retaining terms through j = 6 leaves a
        /// residual below 1E-25, so the accuracy of the result is limited only by double rounding. Agrees
        /// with pyvinecopulib 0.7.6 reference values to 4.5E-16, and against mpmath at 60 decimal digits
        /// on a grid over θ in [1, 100] the worst absolute error is 1.7E-16.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// </para>
        /// <list type="bullet">
        /// <item><description>
        /// Joe, H. (1997). Multivariate Models and Dependence Concepts. Chapman and Hall, London.
        /// </description></item>
        /// <item><description>
        /// Nelsen, R. B. (2006). An Introduction to Copulas, 2nd ed. Springer, New York.
        /// </description></item>
        /// </list>
        /// </remarks>
        internal static double KendallsTauFromTheta(double theta)
        {
            if (theta < 1d || double.IsNaN(theta) || double.IsInfinity(theta))
                throw new ArgumentOutOfRangeException(nameof(theta), "The dependency parameter θ (theta) must be finite and greater than or equal to 1.");

            // The terms are summed from the smallest to the largest to limit accumulated rounding error.
            double sum = 0d;
            for (int k = TauSeriesTerms; k >= 1; k--)
                sum += 1d / (k * (theta * k + 2d) * (theta * (k - 1d) + 2d));

            return 1d - 4d * (sum + TauSeriesTail(theta));
        }

        /// <summary>
        /// Returns the closed-form remainder of the Kendall's τ series beyond the directly summed terms.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        /// <returns>The sum of the series terms for k greater than <see cref="TauSeriesTerms"/>.</returns>
        /// <remarks>
        /// The term 1 / (k * (θ*k + 2) * (θ*(k - 1) + 2)) is written as 1 / (θ² k³ (1 + p/k)(1 + q/k)) with
        /// p = 2/θ and q = (2 - θ)/θ, and the product is expanded as Σ[j ≥ 0] (-1)^j h_j / k^j, where the
        /// complete homogeneous symmetric polynomials h_j satisfy h_(j+1) = p * h_j + q^(j+1). Both |p| and |q|
        /// are far smaller than the truncation point, so the expansion converges geometrically. Each inverse
        /// power is summed over the tail with the Euler-Maclaurin expansion of the Hurwitz zeta function
        /// ζ(m, a) = a^(1-m)/(m-1) + a^(-m)/2 + m/12 * a^(-m-1) - m(m+1)(m+2)/720 * a^(-m-3) + ..., which is
        /// exact to about 1E-25 at a = 1,001.
        /// </remarks>
        private static double TauSeriesTail(double theta)
        {
            double a = TauSeriesTerms + 1d;
            double inverse = 1d / a;
            double p = 2d / theta;
            double q = (2d - theta) / theta;
            double h = 1d;
            double qPower = 1d;
            double aPower = inverse * inverse * inverse;
            double sign = 1d;
            double tail = 0d;

            for (int j = 0; j <= TauSeriesTailOrder; j++)
            {
                double m = 3d + j;
                double zeta = aPower * a / (m - 1d) + 0.5d * aPower + m / 12d * aPower * inverse
                              - m * (m + 1d) * (m + 2d) / 720d * aPower * inverse * inverse * inverse;
                tail += sign * h * zeta;

                sign = -sign;
                qPower *= q;
                h = p * h + qPower;
                aPower *= inverse;
            }

            return tail / (theta * theta);
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
        /// <see cref="ParameterConstraints(IList{double}, IList{double})"/>, θ in [1, 100]. That bracket reaches
        /// τ in [0, 0.98025359]. The Joe copula models positive dependence only, so a negative τ is not attainable
        /// at any θ, and a τ above the upper end means the dependence is too strong to fit within the bracket;
        /// both throw. A τ of 0 is the independence limit, reached at θ = 1, and is assigned directly rather
        /// than handed to a bracket that does not straddle a root.
        /// </remarks>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);

            double L = 1d;
            double U = 100d;

            if (tau < 0d || tau > KendallsTauFromTheta(U))
                throw new ArgumentException("For the Joe copula, tau must be in ~= [0, 0.98025], the range attainable over the fitting bracket θ (theta) of [1, 100]. The Joe copula models positive dependence only, and the dependency in the data is too strong to use the Joe copula.");

            // τ(1) is 0 up to rounding, so any τ at or below it is the independence limit at θ = 1.
            if (tau <= KendallsTauFromTheta(L))
            {
                Theta = L;
                return;
            }

            Theta = Brent.Solve(t => KendallsTauFromTheta(t) - tau, L, U);
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { 1, 100 } };
        }
    }
}