using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Gumbel copula. Sometimes referred to as Gumbel-Hougaard copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class GumbelCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Gumbel copula with a dependency θ = 2.
        /// </summary>
        public GumbelCopula()
        {
            Theta = 2d;
        }

        /// <summary>
        /// Constructs a Gumbel copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public GumbelCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Gumbel copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public GumbelCopula(double theta, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Gumbel; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Gumbel"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "G"; }
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
        /// <remarks>φ(t) = (−ln t)^θ, which is the independence generator −ln t at θ = 1.</remarks>
        public override double Generator(double t)
        {
            return Math.Pow(-Math.Log(t), Theta);
        }

        /// <inheritdoc/>
        /// <remarks>φ⁻¹(s) = exp(−s^(1/θ)).</remarks>
        public override double GeneratorInverse(double t)
        {
            return Math.Exp(-Math.Pow(t, 1d / Theta));
        }

        /// <inheritdoc/>
        /// <remarks>φ′(t) = −θ(−ln t)^(θ−1)/t.</remarks>
        public override double GeneratorPrime(double t)
        {
            return -Theta * Math.Pow(-Math.Log(t), Theta - 1d) / t;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// φ″(t) = θ(−ln t)^(θ−2)·((θ − 1) − ln t)/t², which is positive on 0 &lt; t &lt; 1 because φ is convex. The factor
        /// (θ − 1) − ln t is formed in that order, so it keeps its relative accuracy at θ = 1 and t near 1.
        /// </remarks>
        public override double GeneratorPrime2(double t)
        {
            double y = -Math.Log(t);
            return Theta * Math.Pow(y, Theta - 2d) * ((Theta - 1d) + y) / (t * t);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Solves φ′(t) = s for t. With y = −ln t the equation is (θ − 1)·ln y + y = ln(−s/θ); in z = ln y the left side
        /// (θ − 1)·z + e^z is convex and increasing, so a safeguarded Newton iteration converges to full precision in a
        /// few steps. At θ = 1 the closed form y = ln(−s) is used. Slopes at or above φ′(1), the largest slope, return 1.
        /// </remarks>
        public override double GeneratorPrimeInverse(double t)
        {
            double theta = Theta;
            if (!(t < 0d)) return t >= 0d ? 1d : double.NaN;
            double c = Math.Log(-t) - Math.Log(theta);
            if (theta == 1d) return c > 0d ? Math.Exp(-c) : 1d;
            double start = c > 0d ? Math.Log(c) : Math.Max(c / (theta - 1d), -745d);
            double z = SolveIncreasing(w =>
            {
                double e = Math.Exp(w);
                return ((theta - 1d) * w + e - c, (theta - 1d) + e);
            }, start);
            return Math.Exp(-Math.Exp(z));
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
        /// With x = −ln u, y = −ln v and w = (x^θ + y^θ)^(1/θ), the Gumbel density is
        /// c(u,v) = e^(x + y − w)·(xy)^(θ−1)·w^(1−2θ)·(w + θ − 1) (Nelsen, 2006). It is evaluated in log
        /// space with w = M·e^g, g = Log1p((m/M)^θ)/θ, where M ≥ m are x and y, and x + y − w = m − M·Expm1(g),
        /// which keeps the exponent accurate where w is close to x + y. The log density stays finite for every
        /// θ ≥ 1 and every u and v strictly inside the unit square.
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
            double x = -Math.Log(u), y = -Math.Log(v);
            double high = Math.Max(x, y), low = Math.Min(x, y);
            double g = Tools.Log1p(Math.Pow(low / high, theta)) / theta;
            double w = high * Math.Exp(g);
            double exponent = low - high * Tools.Expm1(g);
            return exponent + (theta - 1d) * (Math.Log(x) + Math.Log(y)) + (1d - 2d * theta) * (Math.Log(high) + g) + Math.Log(w + (theta - 1d));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The closed form h(v|u) = C(u,v)·A^(1/θ−1)·x^(θ−1)/u with x = −ln u, y = −ln v and A = x^θ + y^θ, evaluated in log
        /// space inside the unit square. With M ≥ m the two of x and y and g = Log1p((m/M)^θ)/θ, A^(1/θ) = M·e^g and
        /// ln h = (x − M) − M·Expm1(g) + (1 − θ)·(ln M − ln x + g). Unlike the generic ratio φ′(u)/φ′(C(u,v)), it does not
        /// pass through the rounded C, which loses relative accuracy where C is close to 1. On the boundary the generic
        /// conditional is returned.
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
        /// Inside the unit square the conditional and the density share x, y, M, g and their logarithms, which are
        /// computed once; the values equal <see cref="LogConditionalCDF(double, double)"/> and
        /// <see cref="PDF(double, double)"/> bit for bit.
        /// </remarks>
        protected override double LogConditionalAndDensity(double u, double v, out double density)
        {
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return base.LogConditionalAndDensity(u, v, out density);
            double theta = Theta;
            double x = -Math.Log(u), y = -Math.Log(v);
            double high = Math.Max(x, y), low = Math.Min(x, y);
            double g = Tools.Log1p(Math.Pow(low / high, theta)) / theta;
            double expm1G = Tools.Expm1(g);
            double logX = Math.Log(x), logY = Math.Log(y), logHigh = Math.Log(high);
            double logC = (low - high * expm1G) + (theta - 1d) * (logX + logY) + (1d - 2d * theta) * (logHigh + g) + Math.Log(high * Math.Exp(g) + (theta - 1d));
            density = Math.Exp(logC);
            return (x - high) - high * expm1G + (1d - theta) * (logHigh - logX + g);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Inside the unit square, ln h = (x − M) − M·Expm1(g) + (1 − θ)·(ln M − ln x + g), the exponent of
        /// <see cref="ConditionalCDF(double, double)"/>; elsewhere the log of the conditional.
        /// </remarks>
        protected override double LogConditionalCDF(double u, double v)
        {
            if (!(u > 0d && u < 1d && v > 0d && v < 1d)) return Math.Log(ConditionalCDF(u, v));
            double theta = Theta;
            double x = -Math.Log(u), y = -Math.Log(v);
            double high = Math.Max(x, y), low = Math.Min(x, y);
            double g = Tools.Log1p(Math.Pow(low / high, theta)) / theta;
            return (x - high) - high * Tools.Expm1(g) + (1d - theta) * (Math.Log(high) - Math.Log(x) + g);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The Gumbel conditional has no closed-form inverse, so the closed-form conditional
        /// <see cref="ConditionalCDF(double, double)"/> is inverted with a safeguarded Newton iteration whose derivative is
        /// the density, started at v = t (see <see cref="ArchimedeanCopula.InvertConditionalCDF(double, double, double)"/>),
        /// to full double precision. For 0 &lt; u &lt; 1, t ≤ 0 returns 0 and t ≥ 1 returns 1. The iteration evaluates this class's closed forms
        /// through <see cref="LogConditionalAndDensity(double, double, out double)"/>, so a subclass that overrides
        /// <see cref="ConditionalCDF(double, double)"/> or <see cref="PDF(double, double)"/> must override it as well.
        /// At u = 0 and u = 1, where h is not defined, the inverse is its limit: for θ > 1 and 0 &lt; t &lt; 1,
        /// h(x|u) tends to 1 as u → 0 and to 0 as u → 1 for every 0 &lt; x &lt; 1, so the inverse tends to 0 and 1;
        /// at θ = 1, and for t = 0 or t = 1, it is t.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (u == 0d || u == 1d)
            {
                if (Theta == 1d || !(t > 0d && t < 1d)) return t;
                return u == 0d ? 0d : 1d;
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
        /// The Gumbel copula has no lower tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new GumbelCopula(Theta, CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

        /// <summary>
        /// Estimates the dependency parameter using the method of moments.
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);
            Theta = 1d / (1d - tau);
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { 1, 100 } };
        }

    }
}