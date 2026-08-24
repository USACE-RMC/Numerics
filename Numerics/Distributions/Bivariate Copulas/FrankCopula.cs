using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;
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
        /// <b> Accuracy. </b> The eight pinned pyvinecopulib 0.7.6 oracle points are matched to 3.4E-16.
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
        /// Kendall's τ is estimated from the sample data and <see cref="KendallsTauFromTheta(double)"/> is
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