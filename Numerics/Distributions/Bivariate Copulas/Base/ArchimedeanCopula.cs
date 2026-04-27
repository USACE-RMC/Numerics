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
        public override int NumberOfCopulaParameters => 1;

        /// <inheritdoc/>
        public override double[] GetCopulaParameters => new double[] { Theta };

        /// <inheritdoc/>
        public override void SetCopulaParameters(double[] parameters)
        {
            Theta = parameters[0];
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
            return new ArgumentOutOfRangeException(nameof(Theta),"Parameter is valid");
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