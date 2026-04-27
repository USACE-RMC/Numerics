using System;
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
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);

            // Use conditional probability function 
            double p = v;
            v = Brent.Solve(x =>
            {
                double vu = -(Math.Pow(1d - x, Theta) - 1d) * Math.Pow(Math.Pow(1d - u, Theta) - Math.Pow(1d - u, Theta) * Math.Pow(1d - x, Theta) + Math.Pow(1d - x, Theta), (-Theta + 1d) / Theta) * Math.Pow(1d - u, Theta - 1d);
                return vu - p;
            }, 0d, 1d);
            return [u, v];
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
            return new JoeCopula(Theta, MarginalDistributionX, MarginalDistributionY);
        }

        /// <inheritdoc/>
        public override double[,] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return new double[,] { { 1, 100 } };
        }
    }
}