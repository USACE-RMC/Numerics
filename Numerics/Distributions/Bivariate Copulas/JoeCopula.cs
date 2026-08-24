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
        public override double Generator(double t)
        {
            double a = 1.0d - t;
            return -Math.Log(1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), Theta));
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            double a = 1.0d - Math.Exp(-t);
            return 1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), 1.0d / Theta);
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            double a = 1.0d - t;
            return -(Theta * Math.Sign(a) * Math.Pow(Math.Abs(a), Theta - 1.0d)) / (1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), Theta));
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            double a = 1.0d - t;
            double num = Theta * (Theta + Math.Sign(a) * Math.Pow(Math.Abs(a), Theta) - 1.0d) * Math.Sign(a) * Math.Pow(Math.Abs(a), Theta - 2.0d);
            double aa = 1.0d - Math.Sign(a) * Math.Pow(Math.Abs(a), Theta);
            double den = Math.Sign(aa) * Math.Pow(Math.Abs(aa), 2d);
            return num / den;
        }

        /// <inheritdoc/>
        public override double GeneratorPrimeInverse(double t)
        {
            return Brent.Solve(x => GeneratorPrime(x) - t, 0d, 1d);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The Joe conditional has no closed-form inverse, so the conditional probability
        /// function h(x|u) = (1−u)^(θ−1)·[1 − (1−x)^θ]·A^(1/θ−1),
        /// A = (1−u)^θ + (1−x)^θ − (1−u)^θ(1−x)^θ, is inverted numerically with Brent's
        /// method on x ∈ [0, 1]. In exact arithmetic h(1|u) = 1, but the floating-point
        /// evaluation can round below 1, so a conditional probability within rounding distance
        /// of 1 would otherwise leave the bracket without a sign change; the inverse saturates
        /// at the boundary there.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);

            // Use conditional probability function
            double p = t;
            Func<double, double> f = x =>
            {
                double vu = -(Math.Pow(1d - x, Theta) - 1d) * Math.Pow(Math.Pow(1d - u, Theta) - Math.Pow(1d - u, Theta) * Math.Pow(1d - x, Theta) + Math.Pow(1d - x, Theta), (-Theta + 1d) / Theta) * Math.Pow(1d - u, Theta - 1d);
                return vu - p;
            };
            if (f(1d) <= 0d) return 1d;
            return Brent.Solve(f, 0d, 1d);
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
        /// <param name="theta">The dependency parameter, θ. Must be greater than or equal to 1.</param>
        /// <returns>Kendall's τ for the given θ. The value is 0 at θ = 1 and increases to 1 as θ → ∞.</returns>
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
        /// residual below 1E-25, so the accuracy of the result is limited only by double rounding. The
        /// implementation agrees with pyvinecopulib 0.7.6 to 5E-16 or better over θ in [1, 100].
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
        public static double KendallsTauFromTheta(double theta)
        {
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
        /// Kendall's τ is estimated from the sample data and <see cref="KendallsTauFromTheta(double)"/> is
        /// inverted with Brent's method over the bracket returned by
        /// <see cref="ParameterConstraints(IList{double}, IList{double})"/>, θ in [1, 100]. That bracket reaches
        /// τ in [0, 0.9803]. The Joe copula models positive dependence only, so a negative τ is not attainable
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
                throw new ArgumentException("For the Joe copula, tau must be in [0, 0.9803], the range attainable over the fitting bracket θ (theta) of [1, 100]. The Joe copula models positive dependence only, and the dependency in the data is too strong to use the Joe copula.");

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