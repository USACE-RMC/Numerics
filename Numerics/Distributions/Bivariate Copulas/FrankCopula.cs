using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;
using System;
using System.Collections.Generic;

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
        public override double Generator(double t)
        {
            if (Math.Abs(Theta) < 1e-10) return t;
            return -Math.Log((Math.Exp(-Theta * t) - 1d) / (Math.Exp(-Theta) - 1d));
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            return -Math.Log(Math.Exp(-Theta - t) - Math.Exp(-t) + 1d) / Theta;
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            return Theta / (1d - Math.Exp(Theta * t));
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            return Theta * Theta * Math.Exp(Theta * t) / Math.Pow(1d - Math.Exp(Theta * t), 2d);
        }

        /// <inheritdoc/>
        public override double GeneratorPrimeInverse(double t)
        {
            return Math.Log((Theta - t) / -t) / Theta;
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double num = Theta * (Math.Exp(Theta) - 1d) * Math.Exp(Theta * (1d + u + v));
            double den = Math.Pow(Math.Exp(Theta * (u + v)) - Math.Exp(Theta * (1d + u)) - Math.Exp(Theta * (1d + v)) + Math.Exp(Theta), 2d);
            return num / den;
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return -(1d / Theta) * Math.Log(1d + (Math.Exp(-Theta * u) - 1d) * (Math.Exp(-Theta * v) - 1d) / (Math.Exp(-Theta) - 1d));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the Frank closed-form conditional inversion. The evaluation always runs at the
        /// negated dependency −|θ|, where the logarithm's argument is well conditioned, and maps
        /// a positive-θ request through the Frank reflection identity
        /// h_θ(v|u) = 1 − h_{−θ}(1−v|u): solving h_θ(v|u) = t therefore requires handing the
        /// negative-θ evaluation the complemented conditional 1 − t and reflecting its result,
        /// so that the round trip ConditionalCDF(u, InverseConditionalCDF(u, t)) = t holds on
        /// both dependency branches.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double a = -Math.Abs(Theta);
            double s = Theta > 0d ? 1d - t : t;
            double v = -1d / a * Math.Log((-s * (Math.Exp(-a) - 1d) / (Math.Exp(-a * u) * (s - 1d) - s)) + 1d);
            return Theta > 0d ? 1d - v : v;
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
        /// The Maclaurin coefficients of the order-1 Debye function for the even powers x², x⁴, ... x²⁴.
        /// </summary>
        /// <remarks>
        /// The k-th entry is B(2k) / ((2k + 1) * (2k)!), where B(2k) is the 2k-th Bernoulli number. The entries
        /// were evaluated as exact rationals and rounded once to double.
        /// </remarks>
        private static readonly double[] DebyeSeriesCoefficients =
        {
            2.7777777777777776E-02, -2.7777777777777778E-04, 4.7241118669690098E-06,
            -9.1857730746619641E-08, 1.8978869988971000E-09, -4.0647616451442256E-11,
            8.9216910204564523E-13, -1.9939295860721074E-14, 4.5189800296199183E-16,
            -1.0356517612181247E-17, 2.3952186210261870E-19, -5.5817858743250090E-21
        };

        /// <summary>
        /// The largest number of exponential terms summed in the large-argument branch of the order-1 Debye function.
        /// </summary>
        private const int DebyeMaximumTerms = 1000;

        /// <summary>
        /// The size at which an exponential term is small enough to end the large-argument Debye summation.
        /// </summary>
        private const double DebyeTermTolerance = 1E-20;

        /// <summary>
        /// Returns the order-1 Debye function D₁(x) for any real argument.
        /// </summary>
        /// <param name="x">The argument to evaluate.</param>
        /// <returns>The order-1 Debye function evaluated at the given argument.</returns>
        /// <remarks>
        /// <para>
        /// The order-1 Debye function is D₁(x) = (1/x) ∫[0 to x] t / (e^t - 1) dt, with D₁(0) = 1. Note that
        /// <see cref="Numerics.Mathematics.SpecialFunctions.Debye.Function(double)"/> is the order-3 Debye
        /// function, not this one, so it cannot be used here.
        /// </para>
        /// <para>
        /// Three branches are used. A negative argument is reduced by the reflection D₁(-y) = D₁(y) + y/2 for
        /// y &gt; 0, which is what keeps the negative dependency branch of the Frank copula finite. For
        /// 0 &lt; x ≤ 1 the Maclaurin series D₁(x) = 1 - x/4 + Σ[k ≥ 1] B(2k) x^(2k) / ((2k + 1) (2k)!) is used,
        /// truncated after x²⁴, where the next term is below 1E-22. For x &gt; 1 the integral is written against
        /// its limit π²/6, giving D₁(x) = π²/(6x) - Σ[k ≥ 1] e^(-k x) (1/k + 1/(k² x)), which converges
        /// geometrically and is summed until the term falls below 1E-20.
        /// </para>
        /// </remarks>
        private static double DebyeOrderOne(double x)
        {
            if (x == 0d) return 1d;

            // D₁(-y) = D₁(y) + y/2 for y > 0.
            if (x < 0d) return DebyeOrderOne(-x) - 0.5d * x;

            if (x <= 1d)
            {
                double squared = x * x;
                double power = 1d;
                double sum = 0d;
                for (int k = 0; k < DebyeSeriesCoefficients.Length; k++)
                {
                    power *= squared;
                    sum += DebyeSeriesCoefficients[k] * power;
                }
                return 1d - 0.25d * x + sum;
            }

            double remainder = 0d;
            for (int k = 1; k <= DebyeMaximumTerms; k++)
            {
                double term = Math.Exp(-k * x) * (1d / k + 1d / ((double)k * k * x));
                remainder += term;
                if (term < DebyeTermTolerance) break;
            }
            return Math.PI * Math.PI / (6d * x) - remainder;
        }

        /// <summary>
        /// Returns Kendall's τ (tau) implied by the Frank copula dependency parameter θ (theta).
        /// </summary>
        /// <param name="theta">The dependency parameter, θ. Must be non-zero.</param>
        /// <returns>Kendall's τ for the given θ. The value is negative for θ &lt; 0 and positive for θ &gt; 0.</returns>
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
        /// -1. The implementation agrees with pyvinecopulib 0.7.6 to 5E-16 or better for |θ| in [0.1, 100].
        /// For |θ| below about 0.01 the leading terms of the expression cancel and the absolute accuracy
        /// degrades to about 2E-14, which is far smaller than the τ values involved there.
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
        public static double KendallsTauFromTheta(double theta)
        {
            return 1d - 4d / theta * (1d - DebyeOrderOne(theta));
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
        /// <see cref="ParameterConstraints(IList{double}, IList{double})"/>, which is [0.001, 100] for a
        /// positive τ and [-100, -0.001] for a negative one. That bracket reaches |τ| in [0.000111, 0.9607], so
        /// a |τ| above the upper end means the dependence is too strong to fit within the bracket and throws.
        /// The independence limit θ = 0 leaves the Frank generator and distribution functions indeterminate and
        /// is excluded from the bracket, so a τ of 0, or any τ too small for the bracket to straddle a root, is
        /// assigned the bracket endpoint of matching sign rather than handed to the solver.
        /// </remarks>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);

            if (Math.Abs(tau) > KendallsTauFromTheta(100d))
                throw new ArgumentException("For the Frank copula, tau must be in [-0.9607, 0.9607], the range attainable over the fitting bracket θ (theta) of [-100, 100]. The dependency in the data is too strong to use the Frank copula.");

            // θ = 0 is the independence limit and is excluded from the bracket, so a τ smaller than the bracket
            // can reach is assigned the endpoint nearest independence.
            double nearIndependence = tau < 0d ? -0.001d : 0.001d;
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