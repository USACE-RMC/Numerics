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

// Since I use functions from the Accord Math Library, here is the required license header:
// Haden Smith (January 2019)
// 
// Accord Math Library
// The Accord.NET Framework
// http://accord-framework.net
// 
// Copyright © César Souza, 2009-2017
// cesarsouza at gmail.com
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

using System;
using System.Collections.Generic;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Student's t probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Student%27s_t-distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class StudentT : UnivariateDistributionBase
    {
       
        /// <summary>
        /// Constructs a standard Student's t distribution with 10 degrees of freedom.
        /// </summary>
        public StudentT()
        {
            SetParameters(0d, 1d, 10d);
        }

        /// <summary>
        /// Constructs a standard Student's t distribution with a given degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        public StudentT(double degreesOfFreedom)
        {
            SetParameters(0d, 1d, degreesOfFreedom);
        }

        /// <summary>
        /// Constructs a Student's t distribution with a given location, scale, and degrees of freedom.
        /// </summary>
        /// <param name="location">The location parameter.</param>
        /// <param name="scale">The scale parameter.</param>
        /// <param name="degreesOfFreedom">The degrees of freedom.</param>
        public StudentT(double location, double scale, double degreesOfFreedom)
        {
            SetParameters(location, scale, degreesOfFreedom);
        }

        private double _mu;
        private double _sigma;
        private int _degreesOfFreedom;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, DegreesOfFreedom, false) is null;
                _mu = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter σ (sigma).
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                _parametersValid = ValidateParameters(Mu, value, DegreesOfFreedom, false) is null;
                _sigma = value;
            }
        }

        /// <summary>
        /// Gets and sets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public int DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(Mu, Sigma, value, false) is null;
                _degreesOfFreedom = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.StudentT; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Student's t"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "T"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Location (µ)";
                parmString[1, 0] = "Scale (σ)";
                parmString[2, 0] = "Degrees of Freedom (ν)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                parmString[2, 1] = DegreesOfFreedom.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ", "ν"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma), nameof(DegreesOfFreedom)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma, DegreesOfFreedom]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (DegreesOfFreedom > 1)
                {
                    return Mu;
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (DegreesOfFreedom > 2)
                {
                    return Math.Sqrt(Sigma * Sigma * DegreesOfFreedom / (DegreesOfFreedom - 2));
                }
                else if (DegreesOfFreedom > 1 && DegreesOfFreedom <= 2)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (DegreesOfFreedom > 3)
                {
                    return 0d;
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (DegreesOfFreedom > 4)
                {
                    return 3d + 6d / (DegreesOfFreedom - 4);
                }
                else if (DegreesOfFreedom > 2 && DegreesOfFreedom <= 4)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return double.NaN;
                }
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d, 1.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="mu">The location µ (Mu).</param>
        /// <param name="sigma">The scale σ (sigma). Range: σ > 0.</param>
        /// <param name="v">The degrees of freedom ν (nu). Range: ν > 0.</param>
        public void SetParameters(double mu, double sigma, double v)
        {
            Mu = mu;
            Sigma = sigma;
            DegreesOfFreedom = (int)v;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter µ (Mu).</param>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        /// <param name="v">The degrees of freedom ν (nu).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double location, double scale, double v, bool throwException)
        {
            if (double.IsNaN(location) || double.IsInfinity(location))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
            }
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
            }
            if (v < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double Z = (x - Mu) / Sigma;
            if (DegreesOfFreedom >= 100000000.0d)
            {
                return Normal.StandardPDF(Z);
            }
            return Math.Exp(Gamma.LogGamma((DegreesOfFreedom + 1.0d) / 2.0d) - Gamma.LogGamma(DegreesOfFreedom / 2.0d)) * Math.Pow(1.0d + Z * Z / DegreesOfFreedom, -0.5d * (DegreesOfFreedom + 1.0d)) / Math.Sqrt(DegreesOfFreedom * Math.PI);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            // 
            double Z = (x - Mu) / Sigma;
            if (DegreesOfFreedom >= 100000000.0d)
            {
                return Normal.StandardCDF(Z);
            }
            double betaX = (Z + Math.Sqrt(Z * Z + DegreesOfFreedom)) / (2d * Math.Sqrt(Z * Z + DegreesOfFreedom));
            return Beta.Incomplete(DegreesOfFreedom / 2.0d, DegreesOfFreedom / 2.0d, betaX);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This function is also know as the Quantile Function.
        /// <para>
        /// References:
        /// This code was copied from the Accord Math Library.
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, http://accord-framework.net
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            // 
            if (probability > 0.25d && probability < 0.75d)
            {
                if (probability == 0.5d)
                {
                    return Mu;
                }

                double z = Beta.IncompleteInverse(0.5d, 0.5d * DegreesOfFreedom, Math.Abs(1.0d - 2.0d * probability));
                double t = Math.Sqrt(DegreesOfFreedom * z / (1.0d - z));
                if (probability < 0.5d)
                {
                    t = -t;
                }

                return Mu + Sigma * t;
            }
            else
            {
                int rflg = -1;
                if (probability >= 0.5d)
                {
                    probability = 1.0d - probability;
                    rflg = 1;
                }

                double z = Beta.IncompleteInverse(0.5d * DegreesOfFreedom, 0.5d, 2.0d * probability);
                double t = Math.Sqrt(DegreesOfFreedom / z - DegreesOfFreedom);
                if (double.MaxValue * z < DegreesOfFreedom)
                {
                    return rflg * double.MaxValue;
                }

                return Mu + Sigma * (rflg * t);
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new StudentT(Mu, Sigma, DegreesOfFreedom);
        }
              
    }
}