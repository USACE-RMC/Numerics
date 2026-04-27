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
            return Math.Log((1.0d - Theta * (1.0d - t)) / t);
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            return (1.0d - Theta) / (Math.Exp(t) - Theta);
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            return (Theta - 1.0d) / (t * (Theta * (t - 1.0d) + 1.0d));
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            double num = (Theta - 1.0d) * (Theta * (2.0d * t - 1.0d) + 1.0d);
            double a = Theta * (t - 1.0d) * t + t;
            double den = Math.Sign(a) * Math.Pow(Math.Abs(a), 2d);
            return -num / den;
        }

        /// <inheritdoc/>
        public override double GeneratorPrimeInverse(double t)
        {
            return Brent.Solve(x => GeneratorPrime(x) - t, 0d, 1d);
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return (-1 + Tools.Sqr(Theta) * (-1 + u + v - u * v) - Theta * (-2 + u + v + u * v)) / Math.Pow(-1 + Theta * (-1 + u) * (-1 + v), 3);
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            //Johnson (1987, p.362).
            double w = v;
            double b = 1d - u;
            double A = w * Math.Pow(Theta * b, 2) - Theta;
            double B = Theta + 1d - 2d * Theta * b * w;
            double C = w - 1d;
            v = (-B + Math.Sqrt(B * B - 4d * A * C)) / 2d / A;
            v = 1d - v;
            return [u, v];
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
            return new AMHCopula(Theta, MarginalDistributionX, MarginalDistributionY);
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