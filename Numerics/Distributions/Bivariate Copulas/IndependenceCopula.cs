using System;
using System.Collections.Generic;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Independence (product) copula, Π(u,v) = u·v.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Independence copula couples two marginal distributions with no dependence at all:
    /// C(u,v) = u·v, the copula of any pair of independent random variables per Sklar's theorem.
    /// It is the natural default when no dependence structure has been asserted, and the
    /// identity against which dependence modeling is compared — every copula family that
    /// admits independence in its parameter interior reduces to it there.
    /// </para>
    /// <para>
    /// The copula has no parameters. It is therefore permanently valid: <see cref="ValidateParameter"/>
    /// always returns null, and no assignment to <see cref="BivariateCopula.Theta"/> (including
    /// non-finite values, which the parameterized families reject) can invalidate it, because the
    /// dependency parameter carries no meaning here. <see cref="SetCopulaParameters"/> is a no-op
    /// and <see cref="GetCopulaParameters"/> is empty, so parameter estimation has nothing to do.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Nelsen, R.B. (2006). "An Introduction to Copulas." 2nd ed. Springer. Section 2.4.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Copula_(probability_theory)"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class IndependenceCopula : BivariateCopula
    {

        /// <summary>
        /// Constructs an Independence copula.
        /// </summary>
        public IndependenceCopula()
        {
        }

        /// <summary>
        /// Constructs an Independence copula with marginal distributions.
        /// </summary>
        /// <param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        /// <param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public IndependenceCopula(IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Independence; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Independence"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Π"; }
        }

        /// <inheritdoc/>
        public override string[,] ParameterToString
        {
            get { return new string[0, 2]; }
        }

        /// <inheritdoc/>
        public override string ParameterNameShortForm
        {
            get { return ""; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override int NumberOfCopulaParameters => 0;

        /// <inheritdoc/>
        public override double[] GetCopulaParameters => Array.Empty<double>();

        /// <inheritdoc/>
        /// <remarks>
        /// The Independence copula has no parameters, so this method is a no-op and the
        /// input is ignored.
        /// </remarks>
        public override void SetCopulaParameters(double[] parameters)
        {
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The Independence copula has no parameters, so every input is valid and this
        /// method always returns null. The copula can never enter an invalid state.
        /// </remarks>
        public override ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException)
        {
            return null;
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[0, 2];
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            return 1.0d;
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            return u * v;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Under independence the conditional distribution of V given U = u is the
        /// unconditional uniform, so h(v|u) = v.
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            return v;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Under independence the conditional inversion is the identity in t.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            return t;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Under independence conditional simulation is the identity: the pair (u, v) is
        /// returned unchanged.
        /// </remarks>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, v];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 0.
        /// Independent variables have no tail dependence.
        /// </summary>
        public override double UpperTailDependence => 0.0;

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 0.
        /// Independent variables have no tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new IndependenceCopula(CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

    }
}
