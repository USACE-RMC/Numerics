using System;
using System.Collections.Generic;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The bivariate Student's t-copula.
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
    /// The Student's t-copula is an elliptical copula derived from the bivariate Student's t-distribution.
    /// Unlike the Normal (Gaussian) copula which has zero tail dependence, the t-copula exhibits symmetric
    /// upper and lower tail dependence controlled by the degrees of freedom parameter. As the degrees of
    /// freedom increase to infinity, the t-copula converges to the Normal copula.
    /// </para>
    /// <para>
    /// The copula is parameterized by the correlation coefficient rho (stored in Theta) and the degrees
    /// of freedom nu. The copula density is:
    /// <code>
    ///     c(u,v) = f₂(t⁻¹_ν(u), t⁻¹_ν(v); ρ, ν) / [f_ν(t⁻¹_ν(u)) · f_ν(t⁻¹_ν(v))]
    /// </code>
    /// where f₂ is the bivariate t density with correlation ρ and ν degrees of freedom,
    /// f_ν is the univariate standard t density, and t⁻¹_ν is the univariate t quantile function.
    /// </para>
    /// <para>
    /// Key applications include:
    /// <list type="bullet">
    /// <item><description>Joint flood risk analysis in TotalRisk where tail dependence between flood variables is critical.</description></item>
    /// <item><description>Bivariate frequency analysis in RMC-BestFit for modeling dependent extremes.</description></item>
    /// <item><description>Hydrologic risk assessment where normal copulas underestimate joint extreme event probabilities.</description></item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Demarta, S. and McNeil, A.J. (2005). "The t Copula and Related Copulas."
    /// International Statistical Review, 73(1), 111-129.
    /// </description></item>
    /// <item><description>
    /// Nelsen, R.B. (2006). "An Introduction to Copulas." 2nd ed. Springer. Chapter 4.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Copula_(probability_theory)#Student%27s_t-copula"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class StudentTCopula : BivariateCopula
    {

        /// <summary>
        /// Constructs a bivariate Student's t-copula with default parameters: ρ = 0 and ν = 5.
        /// </summary>
        public StudentTCopula()
        {
            _nu = 5.0;
            Theta = 0.0d;
        }

        /// <summary>
        /// Constructs a bivariate Student's t-copula with specified correlation and degrees of freedom.
        /// </summary>
        /// <param name="rho">The correlation parameter ρ (rho). Must be in [-1, 1].</param>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than 2.</param>
        public StudentTCopula(double rho, double degreesOfFreedom)
        {
            _nu = degreesOfFreedom;
            Theta = rho;
        }

        /// <summary>
        /// Constructs a bivariate Student's t-copula with specified parameters and marginal distributions.
        /// </summary>
        /// <param name="rho">The correlation parameter ρ (rho). Must be in [-1, 1].</param>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than 2.</param>
        /// <param name="marginalDistributionX">The X marginal distribution.</param>
        /// <param name="marginalDistributionY">The Y marginal distribution.</param>
        public StudentTCopula(double rho, double degreesOfFreedom, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            _nu = degreesOfFreedom;
            Theta = rho;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        private double _nu;

        /// <summary>
        /// Gets or sets the degrees of freedom ν (nu). Continuous; must be greater than 2.
        /// </summary>
        /// <remarks>
        /// ν is represented as a <see cref="double"/> so that gradient-free MCMC samplers can
        /// explore the parameter space smoothly; previously the value was rounded to an integer
        /// on every <see cref="SetCopulaParameters"/> call, which produced a step-function
        /// likelihood surface and unnecessary plateaus in the posterior. Both
        /// <see cref="StudentT"/> and <see cref="MultivariateStudentT"/> accept non-integer
        /// degrees of freedom.
        /// </remarks>
        public double DegreesOfFreedom
        {
            get { return _nu; }
            set
            {
                _parametersValid = ValidateParameters(Theta, value, false) is null;
                _nu = value;
            }
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.StudentT; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Student's t"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "t"; }
        }

        /// <inheritdoc/>
        public override string[,] ParameterToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Correlation (ρ)";
                parmString[0, 1] = Theta.ToString();
                parmString[1, 0] = "Degrees of Freedom (ν)";
                parmString[1, 1] = _nu.ToString("G4");
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
        /// <remarks>
        /// The base <see cref="BivariateCopula.Theta"/> setter calls this method on the new
        /// rho value; we delegate to <see cref="ValidateParameters(double, double, bool)"/>
        /// so that <c>_parametersValid</c> reflects the validity of both ρ and the current ν.
        /// </remarks>
        public override ArgumentOutOfRangeException? ValidateParameter(double parameter, bool throwException)
        {
            return ValidateParameters(parameter, _nu, throwException);
        }

        /// <summary>
        /// Validates the correlation ρ and degrees of freedom ν together. Returns null when
        /// both parameters are in their valid domains.
        /// </summary>
        /// <param name="rho">The correlation parameter ρ. Must be in [<see cref="ThetaMinimum"/>, <see cref="ThetaMaximum"/>].</param>
        /// <param name="degreesOfFreedom">The degrees of freedom ν. Must be finite and greater than 2.</param>
        /// <param name="throwException">When true, throws on the first invalid parameter; when false, returns the exception instead.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double rho, double degreesOfFreedom, bool throwException)
        {
            if (rho < ThetaMinimum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be greater than " + ThetaMinimum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be greater than " + ThetaMinimum.ToString() + ".");
            }
            if (rho > ThetaMaximum)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be less than " + ThetaMaximum.ToString() + ".");
                return new ArgumentOutOfRangeException(nameof(Theta), "The correlation parameter ρ (rho) must be less than " + ThetaMaximum.ToString() + ".");
            }
            if (double.IsNaN(degreesOfFreedom) || double.IsInfinity(degreesOfFreedom) || degreesOfFreedom <= 2.0)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom must be greater than 2.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom must be greater than 2.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override int NumberOfCopulaParameters => 2;

        /// <inheritdoc/>
        public override double[] GetCopulaParameters => new double[] { Theta, DegreesOfFreedom };

        /// <inheritdoc/>
        /// <remarks>
        /// ν is clamped to the valid domain (ν &gt; 2) but otherwise kept continuous — no
        /// integer rounding. MCMC samplers that propose continuous ν values see a smooth
        /// likelihood, avoiding the step-function plateaus produced by previous rounding
        /// behavior. The clamp uses 2 + 1e-10 because <c>2.0 + Tools.DoubleMachineEpsilon</c>
        /// rounds back to 2.0 in IEEE 754 (the ULP at 2.0 is 2^-51, larger than ε = 2^-53).
        /// </remarks>
        public override void SetCopulaParameters(double[] parameters)
        {
            Theta = parameters[0];
            DegreesOfFreedom = Math.Max(2.0 + 1E-10, parameters[1]);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The upper bound on ν is 30 rather than ∞ for a statistical reason: for ν ≳ 30 the
        /// Student's t copula is empirically indistinguishable from the Gaussian copula at
        /// typical hydrologic sample sizes (n ≤ a few hundred), so the likelihood is nearly
        /// flat on [30, ∞). A Uniform prior across that flat region would drag the posterior
        /// mean toward mid-range regardless of the data. Users who genuinely need very high ν
        /// should prefer the Normal copula directly.
        /// </remarks>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,]
            {
                { -1 + Tools.DoubleMachineEpsilon, 1 - Tools.DoubleMachineEpsilon },
                { 2.0 + 1E-10, 30 }
            };
        }

        /// <summary>
        /// Returns the copula probability density function (PDF) evaluated at (u, v).
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        /// <returns>The copula density c(u, v).</returns>
        /// <remarks>
        /// <para>
        /// The t-copula density is computed in log-space for numerical stability:
        /// <code>
        ///     log c = LogΓ((ν+2)/2) + LogΓ(ν/2) - 2·LogΓ((ν+1)/2)
        ///           - ½·log(1-ρ²)
        ///           - ((ν+2)/2)·log(1 + Q/(ν(1-ρ²)))
        ///           + ((ν+1)/2)·log(1 + x₁²/ν) + ((ν+1)/2)·log(1 + x₂²/ν)
        /// </code>
        /// where x₁ = t⁻¹_ν(u), x₂ = t⁻¹_ν(v), and Q = x₁² - 2ρx₁x₂ + x₂².
        /// </para>
        /// </remarks>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);

            double r = _theta;
            double nu = _nu;

            // Transform to t-quantiles
            var tDist = new StudentT(0, 1, _nu);
            double x1 = tDist.InverseCDF(u);
            double x2 = tDist.InverseCDF(v);

            // Compute log copula density in log-space for numerical stability
            double r2 = r * r;
            double logC = Gamma.LogGamma((nu + 2.0) / 2.0)
                        + Gamma.LogGamma(nu / 2.0)
                        - 2.0 * Gamma.LogGamma((nu + 1.0) / 2.0)
                        - 0.5 * Math.Log(1.0 - r2);

            // Quadratic form: Q = (x1^2 - 2*rho*x1*x2 + x2^2)
            double Q = x1 * x1 - 2.0 * r * x1 * x2 + x2 * x2;
            logC -= ((nu + 2.0) / 2.0) * Math.Log(1.0 + Q / (nu * (1.0 - r2)));
            logC += ((nu + 1.0) / 2.0) * Math.Log(1.0 + x1 * x1 / nu);
            logC += ((nu + 1.0) / 2.0) * Math.Log(1.0 + x2 * x2 / nu);

            return Math.Exp(logC);
        }

        /// <summary>
        /// Returns the copula cumulative distribution function (CDF) evaluated at (u, v).
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        /// <returns>The copula CDF C(u, v).</returns>
        /// <remarks>
        /// <para>
        /// The t-copula CDF is computed using the bivariate Student's t CDF:
        /// <code>
        ///     C(u, v) = F₂(t⁻¹_ν(u), t⁻¹_ν(v); ρ, ν)
        /// </code>
        /// where F₂ is the bivariate Student's t CDF evaluated via the <see cref="MultivariateStudentT"/> class.
        /// </para>
        /// </remarks>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);

            double r = _theta;

            // Transform to t-quantiles
            var tDist = new StudentT(0, 1, _nu);
            double x1 = tDist.InverseCDF(u);
            double x2 = tDist.InverseCDF(v);

            // Evaluate bivariate t CDF
            var scaleMatrix = new double[,] { { 1.0, r }, { r, 1.0 } };
            var mvt = new MultivariateStudentT(_nu, new double[] { 0.0, 0.0 }, scaleMatrix);
            return mvt.CDF(new double[] { x1, x2 });
        }

        /// <summary>
        /// Returns the inverse CDF (conditional sampling) for the copula.
        /// </summary>
        /// <param name="u">The first uniform variate in (0, 1).</param>
        /// <param name="v">The second uniform variate in (0, 1), used as the conditional probability.</param>
        /// <returns>
        /// A 2-element array [u, v'] where v' is the conditionally sampled variate.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Uses the conditional distribution of the bivariate Student's t:
        /// X₂ | X₁ = x₁ ~ t_{ν+1}(ρ·x₁, √((1-ρ²)(ν+x₁²)/(ν+1)))
        /// </para>
        /// <para>
        /// The algorithm is:
        /// <list type="number">
        /// <item><description>Transform u to t-quantile: x₁ = t⁻¹_ν(u)</description></item>
        /// <item><description>Sample from the conditional t_{ν+1} distribution using v</description></item>
        /// <item><description>Transform the conditional sample back to uniform: v' = t_ν(x₂)</description></item>
        /// </list>
        /// </para>
        /// </remarks>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);

            double r = _theta;
            double nu = _nu;

            // Transform u to t-quantile
            var tNu = new StudentT(0, 1, _nu);
            double x1 = tNu.InverseCDF(u);

            // Conditional distribution: X2|X1=x1 ~ t_{ν+1} with location = ρ·x1, scale = √((1-ρ²)(ν+x1²)/(ν+1))
            var tNu1 = new StudentT(0, 1, _nu + 1);
            double z2 = tNu1.InverseCDF(v);

            double conditionalScale = Math.Sqrt((1.0 - r * r) * (nu + x1 * x1) / (nu + 1.0));
            double x2 = r * x1 + conditionalScale * z2;

            // Transform back to uniform
            v = tNu.CDF(x2);
            return [u, v];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U.
        /// </summary>
        /// <remarks>
        /// <para>
        /// The t-copula has symmetric upper and lower tail dependence:
        /// <code>
        ///     λ_U = λ_L = 2 · t_{ν+1}(-√((ν+1)(1-ρ)/(1+ρ)))
        /// </code>
        /// where t_{ν+1} is the CDF of the univariate Student's t with ν+1 degrees of freedom.
        /// For ρ = -1, λ = 0. For ρ = 1, λ = 1. As ν → ∞, λ → 0 (Normal copula limit).
        /// </para>
        /// </remarks>
        public override double UpperTailDependence
        {
            get
            {
                double r = _theta;
                double nu = _nu;

                if (r >= 1.0) return 1.0;
                if (r <= -1.0) return 0.0;

                double arg = -Math.Sqrt((nu + 1.0) * (1.0 - r) / (1.0 + r));
                var tDist = new StudentT(0, 1, _nu + 1);
                return 2.0 * tDist.CDF(arg);
            }
        }

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L.
        /// The t-copula has symmetric tail dependence, so λ_L = λ_U.
        /// </summary>
        public override double LowerTailDependence => UpperTailDependence;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new StudentTCopula(Theta, _nu, MarginalDistributionX, MarginalDistributionY);
        }

    }
}
