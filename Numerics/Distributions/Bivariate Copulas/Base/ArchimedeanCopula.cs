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

using System;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// Declares common functionality of all Archimedean Copulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class ArchimedeanCopula : BivariateCopula, IArchimedeanCopula
    {

        /// <inheritdoc/>
        public override string[,] ParameterToString
        {
            get
            {
              var parmString = new string[2, 2];
              parmString[0, 0] = "Dependency (θ)";
              parmString[0, 1] = Theta.ToString();
              return parmString;
            }
        }

        /// <inheritdoc/>
        public override string ParameterNameShortForm
        {
            get { return "θ"; }
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
        public abstract double Generator(double t);

        /// <inheritdoc/>
        public abstract double GeneratorInverse(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrime(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrime2(double t);

        /// <inheritdoc/>
        public abstract double GeneratorPrimeInverse(double t);

        /// <inheritdoc/>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double num = -GeneratorPrime2(CDF(u, v)) * GeneratorPrime(u) * GeneratorPrime(v);
            double den = Math.Pow(GeneratorPrime(CDF(u, v)), 3d);
            return num / den;
        }

        /// <inheritdoc/>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return GeneratorInverse(Generator(u) + Generator(v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This method is based on Genest et al. 1986
        /// 1) Two independent uniformly distributed U(0,1) random variates, u and v, are generated.
        /// 2) Two new variables, s and w, are obtained as s = GeneratorPrime(u) / v and w = GeneratorPrimeInverse(s).
        /// 3) Another variable v is obtained as v = GeneratorInverse(Generator(w) - Generator(u))
        /// 4) The pairs u and v are the simulated pair, preserving the dependence structure.
        /// 5) Both these u and v in the range [0,1]. These simulated pairs of u and v are then
        /// back-transformed through their corresponding marginal distributions.
        /// </remarks>
        public override double[] InverseCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double s = GeneratorPrime(u) / v;
            double w = GeneratorPrimeInverse(s);
            v = GeneratorInverse(Generator(w) - Generator(u));
            return [u, v];
        }

    }
}