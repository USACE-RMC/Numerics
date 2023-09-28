﻿using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System;
using System.Collections.Generic;

namespace Numerics.Distributions.Copulas
{

    /// <summary>
    /// Declares common functionality of all Archimedean Copulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     Authors:
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public abstract class ArchimedeanCopula : BivariateCopula, IArchimedeanCopula
    {

        /// <summary>
        /// Returns the name and value of the theta parameter in 2-column array of string.
        /// </summary>
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

        /// <summary>
        /// Returns the distribution parameter name in short form (e.g. θ).
        /// </summary>
        public override string ParameterNameShortForm
        {
            get { return "θ"; }
        }

        /// <summary>
        /// Test to see if distribution parameters are valid.
        /// </summary>
        /// <param name="parameter">Dependency parameter.</param>
        /// <param name="throwException">Boolean indicating whether to throw the exception or not.</param>
        /// <returns>Nothing if the parameters are valid and the exception if invalid parameters were found.</returns>
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

        /// <summary>
        /// The generator function of the copula.
        /// </summary>
        /// <param name="t">The reduced variate.</param>
        public abstract double Generator(double t);

        /// <summary>
        /// The inverse of the generator function.
        /// </summary>
        /// <param name="t">The reduced variate.</param>
        public abstract double GeneratorInverse(double t);

        /// <summary>
        /// The first derivative of the generator function.
        /// </summary>
        /// <param name="t">The reduced variate.</param>
        public abstract double GeneratorPrime(double t);

        /// <summary>
        /// The second derivative of the generator function.
        /// </summary>
        /// <param name="t">The reduced variate.</param>
        public abstract double GeneratorPrime2(double t);

        /// <summary>
        /// The inverse of the first derivative of the generator function.
        /// </summary>
        /// <param name="t">The reduced variate.</param>
        public abstract double GeneratorPrimeInverse(double t);

        /// <summary>
        /// The probability density function (PDF) of the copula evaluated at reduced variates u and v.
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        public override double PDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            double num = -GeneratorPrime2(CDF(u, v)) * GeneratorPrime(u) * GeneratorPrime(v);
            double den = Math.Pow(GeneratorPrime(CDF(u, v)), 3d);
            return num / den;
        }

        /// <summary>
        /// The cumulative distribution function (CDF) of the copula evaluated at reduced variates u and v.
        /// </summary>
        /// <param name="u">The reduced variate between 0 and 1.</param>
        /// <param name="v">The reduced variate between 0 and 1.</param>
        public override double CDF(double u, double v)
        {
            // Validate parameters
            if (_parametersValid == false) ValidateParameter(Theta, true);
            return GeneratorInverse(Generator(u) + Generator(v));
        }

        /// <summary>
        /// The inverse cumulative distribution function (InverseCDF) of the copula evaluated at probabilities u and v.
        /// </summary>
        /// <param name="u">Probability between 0 and 1.</param>
        /// <param name="v">Probability between 0 and 1.</param>
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
            return new[] { u, v };
        }

    }
}