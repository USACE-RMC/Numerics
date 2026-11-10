using Numerics.Data.Statistics;
using System;
using System.Collections.Generic;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Clayton copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class ClaytonCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Clayton copula with a dependency θ = 2.
        /// </summary>
        public ClaytonCopula()
        {
            Theta = 2d;
        }

        /// <summary>
        /// Constructs a Clayton copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public ClaytonCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Clayton copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public ClaytonCopula(double theta, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Clayton; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Clayton"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "C"; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return -1.0d; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The generator of Nelsen (2006, Table 4.1), φ(t) = (t^(−θ) − 1)/θ, evaluated as Expm1(−θ·ln t)/θ where
        /// |θ·ln t| ≤ 1/2 and as (e^(−θ·ln t) − 1)/θ otherwise, so it has no cancellation for small |θ|. θ = 0 returns
        /// the independence generator −ln t, its limit.
        /// </remarks>
        public override double Generator(double t)
        {
            double theta = Theta;
            double logT = Math.Log(t);
            if (theta == 0d) return -logT;
            double x = -theta * logT;
            return Tools.Expm1(x) / theta;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ⁻¹(s) = (1 + θs)^(−1/θ), evaluated as exp(−Log1p(θs)/θ); where θs &lt; −1/2 the sum 1 + θs is formed with the
        /// exact residual of the product (Dekker's two-product), so it keeps its relative accuracy near zero. For θ &lt; 0
        /// the pseudo-inverse is 0 where 1 + θs ≤ 0. θ = 0 returns e^(−s), its limit.
        /// </remarks>
        public override double GeneratorInverse(double t)
        {
            double theta = Theta;
            if (theta == 0d) return Math.Exp(-t);
            double y = theta * t;
            if (y < -0.5d)
            {
                // 1 + θs near zero: add the exact residual of the product so the small sum keeps its relative accuracy.
                double sum = (1d + y) + Mathematics.SpecialFunctions.Beta.ProductError(theta, t, y);
                return sum > 0d ? Math.Exp(-Math.Log(sum) / theta) : 0d;
            }
            if (!(y > -1d)) return double.NaN;
            return Math.Exp(-Tools.Log1p(y) / theta);
        }

        /// <inheritdoc/>
        /// <remarks>φ′(t) = −t^(−θ−1), which is −1/t at θ = 0.</remarks>
        public override double GeneratorPrime(double t)
        {
            return -Math.Pow(t, -Theta - 1d);
        }

        /// <inheritdoc/>
        /// <remarks>φ″(t) = (θ + 1)·t^(−θ−2), which is 1/t² at θ = 0.</remarks>
        public override double GeneratorPrime2(double t)
        {
            return (Theta + 1d) * Math.Pow(t, -Theta - 2d);
        }

        /// <inheritdoc/>
        /// <remarks>The closed form (−s)^(−1/(θ+1)) of φ′(t) = s, which is −1/s at θ = 0.</remarks>
        public override double GeneratorPrimeInverse(double t)
        {
            return Math.Pow(-t, -1d / (Theta + 1d));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// C(u,v) = max(S, 0)^(−1/θ) with S = u^(−θ) + v^(−θ) − 1. Inside the unit square it is evaluated as
        /// exp(−ln S / θ), with ln S from <see cref="LogS(double, double)"/>, which is smooth through θ = 0: as θ → 0,
        /// ln S / θ → ln u + ln v, and θ = 0 returns the independence limit uv. A zero argument returns 0.
        /// </remarks>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double theta = Theta;
            if (u >= 0d && u <= 1d && v >= 0d && v <= 1d)
            {
                if (theta == 0d) return u * v;
                if (u == 0d || v == 0d) return 0d;
                return Math.Exp(-LogS(-theta * Math.Log(u), -theta * Math.Log(v)) / theta);
            }
            return Math.Pow(Math.Max(Math.Pow(u, -theta) + Math.Pow(v, -theta) - 1d, 0d), -1d / theta);
        }

        /// <summary>
        /// Returns ln S for S = u^(−θ) + v^(−θ) − 1, given a = −θ·ln u and b = −θ·ln v; negative infinity where S ≤ 0.
        /// </summary>
        /// <param name="a">−θ·ln u.</param>
        /// <param name="b">−θ·ln v.</param>
        /// <returns>The log of S, or negative infinity when S is not positive.</returns>
        /// <remarks>
        /// With m ≥ n the larger and smaller of a and b: when m > 1, ln S = m + ln(e^(a−m) + e^(b−m) − e^(−m)); when
        /// n ≥ −1/2, ln S = Log1p(Expm1(a) + Expm1(b)); otherwise (θ &lt; 0) ln S = ln(e^n + Expm1(m)). Each form adds
        /// terms of one sign or a small term to a larger one, the powers u^(−θ) and v^(−θ) never form, and for small
        /// |θ| the result keeps its relative accuracy, so ln S / θ is smooth through θ = 0.
        /// </remarks>
        private static double LogS(double a, double b)
        {
            double m = Math.Max(a, b), n = Math.Min(a, b);
            if (m > 1d)
                return m + Math.Log(Math.Exp(a - m) + Math.Exp(b - m) - Math.Exp(-m));
            if (n >= -0.5d)
            {
                // Both exponential terms are within a factor of e^0.5 of 1, and Expm1 of each has one sign when θ ≠ 0.
                return Tools.Log1p(Tools.Expm1(a) + Tools.Expm1(b));
            }
            // θ < 0 with e^n small: S = e^n + Expm1(m) adds a small positive term to a term in (−1, 0].
            double s = Math.Exp(n) + Tools.Expm1(m);
            return s > 0d ? Math.Log(s) : double.NegativeInfinity;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Inside the open unit square the density is <c>Math.Exp(LogPDF(u, v))</c>, so it shares the log-space
        /// evaluation of <see cref="LogPDF(double, double)"/>. At θ = 0 it is the independence density 1 on the closed
        /// unit square. Otherwise the generic Archimedean density of the base class is returned.
        /// </remarks>
        public override double PDF(double u, double v)
        {
            if (Theta == 0d && u >= 0d && u <= 1d && v >= 0d && v <= 1d)
            {
                if (_parametersValid == false) ValidateParameter(Theta, true);
                return 1d;
            }
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.PDF(u, v);
            return Math.Exp(LogPDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The Clayton density c(u,v) = (1 + θ)·(uv)^(−θ−1)·S^(−2−1/θ), with S = u^(−θ) + v^(−θ) − 1 (Nelsen,
        /// 2006), is evaluated in log space from a = −θ·ln u and b = −θ·ln v, with ln S from
        /// <see cref="LogS(double, double)"/>, so the powers u^(−θ) and v^(−θ) never form and the log density stays finite
        /// for large θ. For θ &lt; 0 the density is zero, and the log density negative infinity, where S ≤ 0.
        /// </para>
        /// <para>
        /// The form is smooth through θ = 0, where ln S / θ tends to ln u + ln v, and θ = 0 returns the independence
        /// limit 0 on the closed unit square. On the boundary and outside the square for θ ≠ 0 the log of the generic
        /// Archimedean density is returned.
        /// </para>
        /// </remarks>
        public override double LogPDF(double u, double v)
        {
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (Theta == 0d && u >= 0d && u <= 1d && v >= 0d && v <= 1d) return 0d;
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.LogPDF(u, v);
            double theta = Theta;
            double logU = Math.Log(u), logV = Math.Log(v);
            double logS = LogS(-theta * logU, -theta * logV);
            if (double.IsNegativeInfinity(logS)) return double.NegativeInfinity;
            return Tools.Log1p(theta) - (theta + 1d) * (logU + logV) - 2d * logS - logS / theta;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// h(v|u) = u^(−θ−1)·S^(−1−1/θ), evaluated inside the unit square as exp(−(θ+1)·ln u − ln S − ln S / θ) with ln S
        /// from <see cref="LogS(double, double)"/>, which is smooth through θ = 0; θ = 0 returns the independence limit v.
        /// Where S ≤ 0 (θ &lt; 0) the conditional is 0. On the boundary the generic Archimedean conditional is returned.
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (v == 0d) return 0d;
            if (v == 1d) return 1d;
            double theta = Theta;
            if (theta == 0d) return v;
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.ConditionalCDF(u, v);
            double logU = Math.Log(u);
            double logS = LogS(-theta * logU, -theta * Math.Log(v));
            if (double.IsNegativeInfinity(logS)) return 0d;
            return Math.Exp(-(theta + 1d) * logU - logS - logS / theta);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the Clayton closed-form conditional inversion
        /// v = (u^(−θ)·(t^(−θ/(θ+1)) − 1) + 1)^(−1/θ), evaluated as exp(−L/θ) with
        /// L = Log1p(u^(−θ)·e), e = Expm1(−θ/(θ+1)·ln t), and, where u^(−θ)·e exceeds 1 or u^(−θ) would overflow, from
        /// z = −θ·ln u + ln e as L = z + Log1p(e^(−z)) or Log1p(e^z). For θ &lt; 0, where u^(−θ)·e is below −1/2,
        /// L = ln((1 − u^(−θ)) + u^(−θ)·t^(−θ/(θ+1))), a sum of two non-negative terms. The form is smooth through
        /// θ = 0, where L/θ tends to −ln t, and
        /// θ = 0 returns the independence limit t. At u = 0 with θ > 0 it gives the limit 0 for t &lt; 1; t = 1 returns 1
        /// for every u.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (t == 1d) return 1d;
            double theta = Theta;
            if (theta == 0d) return t;
            if (!(u >= 0d && u <= 1d && t >= 0d && t <= 1d))
                return Math.Pow(Math.Pow(u, -theta) * (Math.Pow(t, -theta / (theta + 1d)) - 1d) + 1d, -1d / theta);
            double a = -theta * Math.Log(u);
            double c = -theta / (theta + 1d) * Math.Log(t);
            double e = Tools.Expm1(c);
            double log;
            double z = e > 0d ? a + Math.Log(e) : double.NegativeInfinity;
            if (z > 0d)
                log = z + Tools.Log1p(Math.Exp(-z));
            else if (a < 700d && !(theta < 0d && Math.Exp(a) * e < -0.5d))
                log = Tools.Log1p(Math.Exp(a) * e);
            else if (theta < 0d)
                // 1 + u^(−θ)·e = (1 − u^(−θ)) + u^(−θ)·t^(−θ/(θ+1)), a sum of two non-negative terms, when the product is near −1.
                log = Math.Log(-Tools.Expm1(a) + Math.Exp(a + c));
            else
                log = Tools.Log1p(Math.Exp(z));
            return Math.Exp(-log / theta);
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 0.
        /// The Clayton copula has no upper tail dependence.
        /// </summary>
        public override double UpperTailDependence => 0.0;

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 2^(-1/θ).
        /// </summary>
        public override double LowerTailDependence
        {
            get
            {
                if (Theta <= 0.0) return 0.0;
                return Math.Pow(2.0, -1.0 / Theta);
            }
        }

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new ClaytonCopula(Theta, CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

        /// <summary>
        /// Estimates the dependency parameter using the method of moments.
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);
            Theta = (2d * tau) / (1 - tau);
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { -1, 100 } };
        }

    }
}