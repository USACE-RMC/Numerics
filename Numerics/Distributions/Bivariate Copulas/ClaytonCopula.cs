/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
        public ClaytonCopula(double theta, IUnivariateDistribution marginalDistributionX, IUnivariateDistribution marginalDistributionY)
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
        public override double Generator(double t)
        {
            double a = t;
            return Math.Sign(a) * Math.Pow(Math.Abs(a), -Theta) - 1.0d;
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            double a = 1.0d + t;
            return Math.Sign(a) * Math.Pow(Math.Abs(a), -1.0d / Theta);
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            double a = t;
            return -Theta * Math.Sign(a) * Math.Pow(Math.Abs(a), -Theta - 1.0d);
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            double a = t;
            return -Theta * Math.Sign(a) * Math.Pow(Math.Abs(a), -Theta - 2.0d) * (-Theta - 1.0d);
        }

        /// <inheritdoc/>
        public override double GeneratorPrimeInverse(double t)
        {
            double a = t / -Theta;
            return Math.Sign(a) * Math.Pow(Math.Abs(a), -1.0d / (Theta + 1.0d));
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return Math.Pow(Math.Max(Math.Pow(u,-Theta) + Math.Pow(v, -Theta) - 1d,0d), -1d / Theta);
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            v = Math.Pow(Math.Pow(u, -Theta) * (Math.Pow(v, -Theta / (Theta + 1d)) - 1d) + 1d, -1d / Theta);
            return [u, v];
        }

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new ClaytonCopula(Theta, MarginalDistributionX, MarginalDistributionY);
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
        public override double[] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            return [-1, 100];
        }

    }
}