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
            // A valid parameter must return null: the Theta setter derives ParametersValid
            // from 'ValidateParameter(value, false) is null', so returning a sentinel
            // exception here left ParametersValid permanently false for every Archimedean
            // family that did not override this method (Clayton, Gumbel, Joe).
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
        /// <para>
        /// For an Archimedean copula with generator φ, differentiating C(u,v) = φ⁻¹(φ(u) + φ(v))
        /// with respect to u gives the exact generic conditional h(v|u) = φ′(u) / φ′(C(u,v)).
        /// The ratio is insensitive to a family's internal sign convention for φ′.
        /// </para>
        /// <para>
        /// Per family the ratio reduces to the closed forms:
        /// Clayton h = u^(−θ−1)·(u^(−θ) + v^(−θ) − 1)^(−1−1/θ);
        /// Frank h = e^(−θu)(e^(−θv) − 1) / [(e^(−θ) − 1) + (e^(−θu) − 1)(e^(−θv) − 1)];
        /// Gumbel h = C(u,v)·A^(1/θ−1)·(−ln u)^(θ−1)/u with A = (−ln u)^θ + (−ln v)^θ;
        /// Joe h = (1−u)^(θ−1)·[1 − (1−v)^θ]·A^(1/θ−1) with A = (1−u)^θ + (1−v)^θ − (1−u)^θ(1−v)^θ;
        /// Ali-Mikhail-Haq h = v(1 − θ(1−v))/D² with D = 1 − θ(1−u)(1−v).
        /// </para>
        /// <para>
        /// The generic generator ratio has conditioning domain 0 &lt; u &lt; 1. Within that domain,
        /// the dependent-variable boundaries are exact: h(0|u) = 0 and h(1|u) = 1. Conditioning
        /// endpoint limits are family-specific, so this base implementation does not impose a
        /// family-independent value at u = 0 or u = 1.
        /// </para>
        /// </remarks>
        public override double ConditionalCDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            if (v == 0d) return 0d;
            if (v == 1d) return 1d;
            return GeneratorPrime(u) / GeneratorPrime(CDF(u, v));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This method is based on Genest et al. 1986
        /// 1) Two independent uniformly distributed U(0,1) random variates, u and t, are generated.
        /// 2) Two new variables, s and w, are obtained as s = GeneratorPrime(u) / t and w = GeneratorPrimeInverse(s).
        /// 3) The dependent variate is obtained as v = GeneratorInverse(Generator(w) - Generator(u))
        /// 4) The pairs u and v are the simulated pair, preserving the dependence structure.
        /// 5) Both these u and v are in the range [0,1]. Simulated pairs of u and v are then
        /// back-transformed through their corresponding marginal distributions.
        /// This is the exact inverse of the generic conditional h(v|u) = φ′(u)/φ′(C(u,v)):
        /// setting h = t and solving gives C = φ′⁻¹(φ′(u)/t) and v = φ⁻¹(φ(C) − φ(u)).
        /// </remarks>
        public override double InverseConditionalCDF(double u, double t)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double s = GeneratorPrime(u) / t;
            double w = GeneratorPrimeInverse(s);
            return GeneratorInverse(Generator(w) - Generator(u));
        }

        /// <inheritdoc/>
        public override double[] InverseCDF(double u, double v)
        {
            return [u, InverseConditionalCDF(u, v)];
        }

    }
}
