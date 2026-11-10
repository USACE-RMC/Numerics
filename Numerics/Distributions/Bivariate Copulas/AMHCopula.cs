using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// The Ali-Mikhail-Haq (AHM) copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class AMHCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Ali-Mikhail-Haq (AHM) copula with a dependency θ = 2.
        /// </summary>
        public AMHCopula()
        {
            Theta = 0.0;
        }

        /// <summary>
        /// Constructs a Ali-Mikhail-Haq (AHM) copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public AMHCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Ali-Mikhail-Haq (AHM) copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public AMHCopula(double theta, IUnivariateDistribution? marginalDistributionX, IUnivariateDistribution? marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.AliMikhailHaq; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Ali-Mikhail-Haq"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "AHM"; }
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
        /// φ(t) = ln((1 − θ(1 − t))/t) = Log1p(r) with r = (1 − θ)(1 − t)/t, the exact difference of the ratio from 1,
        /// which is non-negative for θ &lt; 1. Above r = 1/4 it is evaluated as ln(1 + r), whose rounding is then at most
        /// about 2 ulps of the logarithm. It is the independence generator −ln t at θ = 0.
        /// </remarks>
        public override double Generator(double t)
        {
            double r = (1d - Theta) * (1d - t) / t;
            return r > 0.25d ? Math.Log(1d + r) : Tools.Log1p(r);
        }

        /// <inheritdoc/>
        /// <remarks>φ⁻¹(s) = (1 − θ)/(e^s − θ).</remarks>
        public override double GeneratorInverse(double t)
        {
            return (1d - Theta) / (Math.Exp(t) - Theta);
        }

        /// <inheritdoc/>
        /// <remarks>φ′(t) = (θ − 1)/(t·((1 − θ) + θt)).</remarks>
        public override double GeneratorPrime(double t)
        {
            double theta = Theta;
            return (theta - 1d) / (t * ((1d - theta) + theta * t));
        }

        /// <inheritdoc/>
        /// <remarks>φ″(t) = (1 − θ)·((1 − θ) + 2θt)/(t·((1 − θ) + θt))².</remarks>
        public override double GeneratorPrime2(double t)
        {
            double theta = Theta;
            double a = t * ((1d - theta) + theta * t);
            return (1d - theta) * ((1d - theta) + 2d * theta * t) / (a * a);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The closed form of φ′(t) = s: the root in (0, 1] of sθ·t² + s(1 − θ)·t + (1 − θ) = 0, taken as
        /// t = (1 − θ)/q with q = (√D − s(1 − θ))/2 and D = s(1 − θ)·(s(1 − θ) − 4θ), which adds two positive terms and
        /// is −1/s at θ = 0. Slopes at or above φ′(1) = θ − 1 return 1.
        /// </remarks>
        public override double GeneratorPrimeInverse(double t)
        {
            double theta = Theta;
            if (!(t <= theta - 1d)) return t > theta - 1d ? 1d : double.NaN;
            double b = t * (1d - theta);
            double q = 0.5d * (Math.Sqrt(b * (b - 4d * theta)) - b);
            return (1d - theta) / q;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// At θ = 1 the generator is identically zero, so the closed form C(u, v) = uv / (1 − θ(1 − u)(1 − v)) is used
        /// there; below θ = 1 the generator form of <see cref="ArchimedeanCopula"/> is used.
        /// </remarks>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (Theta < 1d) return base.CDF(u, v);
            if (u <= 0d || v <= 0d) return 0d;
            return u * v / (1d - (1d - u) * (1d - v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// At θ = 1 the generator is identically zero, so the closed form h(v|u) = v(1 − θ(1 − v)) / D², with
        /// D = 1 − θ(1 − u)(1 − v), is used there; below θ = 1 the generator ratio of <see cref="ArchimedeanCopula"/>
        /// is used.
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (Theta < 1d) return base.ConditionalCDF(u, v);
            if (v <= 0d) return 0d;
            if (v >= 1d) return 1d;
            double d = 1d - (1d - u) * (1d - v);
            return v * v / (d * d);
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return (-1 + Tools.Sqr(Theta) * (-1 + u + v - u * v) - Theta * (-2 + u + v + u * v)) / Math.Pow(-1 + Theta * (-1 + u) * (-1 + v), 3);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Uses the Ali-Mikhail-Haq closed-form conditional inversion: the quadratic of Johnson (1987, p.362), written in v
        /// itself rather than in 1 − v so a small v keeps its relative accuracy, and solved with the cancellation-free pair
        /// of root formulas. It is defined at θ = 0 (independence, where the quadratic term vanishes and v = t). t = 0
        /// returns 0 and t = 1 returns 1, where at θ = 1 and u = 0 the quadratic degenerates.
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            // h(v|u) = t is the quadratic a·v² + b·v + c = 0 in v itself (Johnson, 1987, p. 362, solves it in 1 − v),
            // with k = 1 − θ(1 − u): a = θ(1 − tθ(1 − u)²), b = (1 − θ) − 2tθ(1 − u)k, c = −t·k².
            if (t == 0d) return 0d;
            if (t == 1d) return 1d;
            double theta = Theta;
            double bar = 1d - u;
            double k = 1d - theta * bar;
            double a = theta * (1d - t * theta * bar * bar);
            double b = (1d - theta) - 2d * t * theta * bar * k;
            double c = -t * k * k;
            // The root in [0, 1], from the cancellation-free pair q/a and c/q with q = −(b + sign(b)·√D)/2.
            double q = -0.5d * (b + (b >= 0d ? 1d : -1d) * Math.Sqrt(b * b - 4d * a * c));
            return b >= 0d ? c / q : q / a;
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

        /// <summary>
        /// Gets the upper tail dependence coefficient λ_U = 0.
        /// The AMH copula has no tail dependence.
        /// </summary>
        public override double UpperTailDependence => 0.0;

        /// <summary>
        /// Gets the lower tail dependence coefficient λ_L = 0.
        /// The AMH copula has no tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new AMHCopula(Theta, CloneMarginal(MarginalDistributionX), CloneMarginal(MarginalDistributionY));
        }

        /// <summary>
        /// Estimates the dependency parameter using the method of moments.
        /// </summary>
        /// <param name="sampleDataX">The sample data for the X variable.</param>
        /// <param name="sampleDataY">The sample data for the Y variable.</param>
        public void SetThetaFromTau(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);

            if ((tau < (5d - 8d * Math.Log(2d)) / 3d) || (tau > 1d / 3d))
                throw new Exception("For the AMH copula, tau must be in [(5 - 8 log 2) / 3, 1 / 3] ~= [-0.1817, 0.3333]. The dependency in the data is too strong to use the AMH copula.");

            double L = tau > 0 ? 0.001d : -1d + Tools.DoubleMachineEpsilon;
            double U = tau > 0 ? 1d - Tools.DoubleMachineEpsilon : -0.001d;

            Theta = Brent.Solve(t => 
            {
                var x = 1d - 2d * (Math.Pow(1d - t, 2d) * Math.Log(-t + 1d) + t) / (3d * t * t);
                return x - tau;
            }, L, U);

        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { -1 + Tools.DoubleMachineEpsilon, 1 - Tools.DoubleMachineEpsilon } };
        }

    }
}