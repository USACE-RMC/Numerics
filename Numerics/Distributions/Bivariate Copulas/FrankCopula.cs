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
    /// The Frank copula.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class FrankCopula : ArchimedeanCopula
    {

        /// <summary>
        /// Constructs a Frank copula with a dependency θ = 2.
        /// </summary>
        public FrankCopula()
        {
            Theta = 2d;
        }

        /// <summary>
        /// Constructs a Frank copula with a specified θ.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        public FrankCopula(double theta)
        {
            Theta = theta;
        }

        /// <summary>
        /// Constructs a Frank copula with a specified θ and marginal distributions.
        /// </summary>
        /// <param name="theta">The dependency parameter, θ.</param>
        ///<param name="marginalDistributionX">The X marginal distribution for the copula.</param>
        ///<param name="marginalDistributionY">The Y marginal distribution for the copula.</param>
        public FrankCopula(double theta, IUnivariateDistribution marginalDistributionX, IUnivariateDistribution marginalDistributionY)
        {
            Theta = theta;
            MarginalDistributionX = marginalDistributionX;
            MarginalDistributionY = marginalDistributionY;
        }

        /// <inheritdoc/>
        public override CopulaType Type
        {
            get { return CopulaType.Frank; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Frank"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "F"; }
        }

        /// <inheritdoc/>
        public override double ThetaMinimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double ThetaMaximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameter(double parameter, bool throwException)
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
            return -Math.Log((Math.Exp(-Theta * t) - 1d) / (Math.Exp(-Theta) - 1d));
        }

        /// <inheritdoc/>
        public override double GeneratorInverse(double t)
        {
            return -Math.Log(Math.Exp(-Theta - t) - Math.Exp(-t) + 1d) / Theta;
        }

        /// <inheritdoc/>
        public override double GeneratorPrime(double t)
        {
            return Theta / (1d - Math.Exp(Theta * t));
        }

        /// <inheritdoc/>
        public override double GeneratorPrime2(double t)
        {
            return Theta * Theta * Math.Exp(Theta * t) / Math.Pow(1d - Math.Exp(Theta * t), 2d);
        }

        /// <inheritdoc/>
        public override double GeneratorPrimeInverse(double t)
        {
            return Math.Log((Theta - t) / -t) / Theta;
        }

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double num = Theta * (Math.Exp(Theta) - 1d) * Math.Exp(Theta * (1d + u + v));
            double den = Math.Pow(Math.Exp(Theta * (u + v)) - Math.Exp(Theta * (1d + u)) - Math.Exp(Theta * (1d + v)) + Math.Exp(Theta), 2d);
            return num / den;
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return -(1d / Theta) * Math.Log(1d + (Math.Exp(-Theta * u) - 1d) * (Math.Exp(-Theta * v) - 1d) / (Math.Exp(-Theta) - 1d));
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double a = -Math.Abs(Theta);
            v = -1d / a * Math.Log((-v * (Math.Exp(-a) - 1d) / (Math.Exp(-a * u) * (v - 1d) - v)) + 1d);
            v = Theta > 0d ? 1d - v : v;
            return [u, v];
        }

        /// <inheritdoc/>
        public override BivariateCopula Clone()
        {
            return new FrankCopula(Theta, MarginalDistributionX, MarginalDistributionY);
        }

        /// <inheritdoc/>
        public override double[] ParameterConstraints(IList<double> sampleDataX, IList<double> sampleDataY)
        {
            var tau = Correlation.KendallsTau(sampleDataX, sampleDataY);
            double L = tau > 0 ? 0.001d : -100d;
            double U = tau > 0 ? 100d : -0.001d;
            return [L, U];
        }

    }
}