using System;
using System.Collections.Generic;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Gaussian (Normal) copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class NormalCopula : BivariateCopula
    {

        /// <summary>
        /// Constructs a bivariate Gaussian copula with a correlation rho ρ = 0.0.
        /// </summary>
        public NormalCopula()
        {
            Theta = 0.0d;
        }

        /// <summary>
        /// Constructs a bivariate Gaussian copula with a specified correlation rho ρ.
        /// </summary>
        public NormalCopula(double rho)
        {
            Theta = rho;
        }

        /// <summary>
        /// Constructs a bivariate Gaussian copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="rho">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public NormalCopula(double rho, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = rho;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Normal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Normal"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "N"; }
        }

        /// <inheritdoc/>
        public override string[,] ParameterToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Correlation (ρ)";
                parmString[0, 1] = Theta.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string ParameterNameShortForm
        {
            get { return "ρ"; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return -1.0d; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return 1.0d; }
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException)
        {
            if (parameter < ThetaMinimum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be greater than " + ThetaMinimum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be greater than " + ThetaMinimum.ToString() + ".");
            }
            if (parameter > ThetaMaximum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be less than " + ThetaMaximum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be less than " + ThetaMaximum.ToString() + ".");
            }
            return null;
        }

        /// <inheritdoc/>
        public override int NumberOfCopulaParameters => 1;

        /// <inheritdoc/>
        public override double[] GetCopulaParameters => new double[] { Theta };

        /// <inheritdoc/>
        public override void SetCopulaParameters(double[] parameters)
        {
            Theta = parameters[0];
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { -1 + Tools.DoubleMachineEpsilon, 1 - Tools.DoubleMachineEpsilon } };
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            var r = _theta;
            var s = Normal.StandardZ(u);
            var t = Normal.StandardZ(v);
            return 1d / Math.Sqrt(1d - r * r) * Math.Exp(-(r * r * s * s + r * r * t * t - 2 * r * s * t) / (2 * (1d - r * r)));
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            // BivariateCDF implements Genz's BVND which computes Phi2(-h,-k;r).
            // To get the copula C(u,v) = Phi2(Phi^-1(u), Phi^-1(v); r),
            // we pass -Phi^-1(u) = Phi^-1(1-u) as arguments.
            return MultivariateNormal.BivariateCDF(-Normal.StandardZ(u), -Normal.StandardZ(v), _theta);
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double z1 = Normal.StandardZ(u);
            double z2 = Normal.StandardZ(v);
            double r = _theta;
            double w2 = r * z1 + Math.Sqrt(1d - r * r) * z2;
            v = Normal.StandardCDF(w2);
            return [u, v];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 0.
        /// The Normal copula has no upper tail dependence.
        /// </summary>
        public override double UpperTailDependence => 0.0;

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 0.
        /// The Normal copula has no lower tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new NormalCopula(Theta, MarginalDistributionX, MarginalDistributionY);
        }

    }
}