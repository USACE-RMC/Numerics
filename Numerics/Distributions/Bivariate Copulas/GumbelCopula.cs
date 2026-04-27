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
        public override double Generator(double t)
        {
            double a = -Math.Log(t);
            return Math.Sign(a) * Math.Pow(Math.Abs(a), Theta);
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            double a = -t;
            return Math.Exp(Math.Sign(a) * Math.Pow(Math.Abs(a), 1.0d / Theta));
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            double a = Math.Log(t);
            return -Theta * Math.Sign(a) * Math.Pow(Math.Abs(a), Theta - 1.0d) / t;
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            double a = -Math.Log(t);
            return Theta * Math.Sign(a) * Math.Pow(Math.Abs(a), Theta - 2.0d) * (-Theta + Math.Log(t) + 1.0d) / (t * t);
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
                double vu = Math.Pow(-Math.Log(u), Theta - 1d) * Math.Exp(-Math.Pow(Math.Pow(-Math.Log(u), Theta) + Math.Pow(-Math.Log(x), Theta), 1d / Theta)) * Math.Pow(Math.Pow(-Math.Log(u), Theta) + Math.Pow(-Math.Log(x), Theta), 1d / Theta - 1d) / u;
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
        /// The Gumbel copula has no lower tail dependence.
        /// </summary>
        public override double LowerTailDependence => 0.0;

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new GumbelCopula(Theta, MarginalDistributionX, MarginalDistributionY);
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