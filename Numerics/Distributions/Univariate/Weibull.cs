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
using System.Collections.Generic;
using System.Linq;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Weibull probability distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Weibull_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Weibull : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IStandardError, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a Weibull distribution with scale = 10 and shape = 2.
        /// </summary>
        public Weibull()
        {
            SetParameters(10d, 2d);
        }

        /// <summary>
        /// Constructs a Weibull distribution with the given parameters λ and k.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        public Weibull(double scale, double shape)
        {
            SetParameters(scale, shape);
        }

        private double _lambda;
        private double _kappa;

        /// <summary>
        /// Gets and sets the scale parameter λ (lambda).
        /// </summary>
        public double Lambda
        {
            get { return _lambda; }
            set
            {
                _parametersValid = ValidateParameters(value, Kappa, false) is null;
                _lambda = value;
            }
        }

        /// <summary>
        /// Gets and sets the shape parameter κ (kappa).
        /// </summary>
        public double Kappa
        {
            get { return _kappa; }
            set
            {
                _parametersValid = ValidateParameters(Lambda, value, false) is null;
                _kappa = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Weibull; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Weibull"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "W"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Scale (λ)";
                parmString[1, 0] = "Shape (κ)";
                parmString[0, 1] = Lambda.ToString();
                parmString[1, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["λ", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Lambda), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Lambda, Kappa]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Lambda * Gamma.Function(1.0d + 1.0d / Kappa); }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Lambda * Math.Pow(Math.Log(2.0d), 1.0d / Kappa); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Kappa <= 1.0d)
                {
                    return 0.0d;
                }
                else
                {
                    return Lambda * Math.Pow((Kappa - 1.0d) / Kappa, 1.0d / Kappa);
                }
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(Lambda * Lambda * Gamma.Function(1.0d + 2.0d / Kappa) - Mean * Mean); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                double mu = Mean;
                double sigma = StandardDeviation;
                return (Gamma.Function(1.0d + 3.0d / Kappa) * Math.Pow(Lambda, 3.0d) - 3.0d * mu * sigma * sigma - Math.Pow(mu, 3.0d)) / Math.Pow(sigma, 3.0d);
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double g1 = Gamma.Function(1d + 1d / Kappa);
                double g2 = Gamma.Function(1d + 2d / Kappa);
                double g3 = Gamma.Function(1d + 3d / Kappa);
                double g4 = Gamma.Function(1d + 4d / Kappa);
                double num = -6 * Math.Pow(g1, 4d) + 12d * g2 * Math.Pow(g1, 2d) - 3d * Math.Pow(g2, 2d) - 4d * g1 * g3 + g4;
                double den = Math.Pow(g2 - Math.Pow(g1, 2d), 2d);
                return 3d + num / den;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                SetParameters(MLE(sample));
            }
            else
            {
                throw new NotImplementedException();
            }
        }

        /// <inheritdoc/>
        public IUnivariateDistribution Bootstrap(ParameterEstimationMethod estimationMethod, int sampleSize, int seed = -1)
        {
            var newDistribution = new Weibull(Lambda, Kappa);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        public void SetParameters(double scale, double shape)
        {
            Lambda = scale;
            Kappa = shape;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="scale">The scale parameter λ (lambda). Range: λ > 0.</param>
        /// <param name="shape">The shape parameter κ (kappa). Range: k > 0.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double scale, double shape, bool throwException)
        {
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Lambda), "The scale parameter λ (lambda) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Lambda), "The scale parameter λ (lambda) must be positive.");
            }
            if (double.IsNaN(shape) || double.IsInfinity(shape) || shape <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            initialVals = SolveMLE(sample);
            // Get bounds of scale
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of shape
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        public double[] MLE(IList<double> sample)
        {
            // Set constraints
            var tuple = GetParameterConstraints(sample);
            var Initials = tuple.Item1;
            var Lowers = tuple.Item2;
            var Uppers = tuple.Item3;

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var W = new Weibull();
                W.SetParameters(x);
                return W.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <summary>
        /// The Maximum Likelihood Estimation method for the Weibull distribution.
        /// </summary>
        /// <param name="samples">The array of sample data.</param>
        /// <remarks>
        /// Implemented according to: Parameter estimation of the Weibull probability distribution, 1994, Hongzhu Qiao, Chris P. Tsokos
        /// <para>
        /// References:
        /// This code was copied and modified from the Math.NET Library.
        /// <list type="bullet">
        /// <item><description>
        /// Math.NET Numerics Library, http://numerics.mathdotnet.com
        /// </description></item>
        /// </list>
        /// </para>
        /// </remarks>
        public double[] SolveMLE(IList<double> samples)
        {
            double n = samples.Count;
            if (n <= 1d)
            {
                throw new Exception("Observations not sufficient. There must be more than 1 data point.");
            }

            double s1 = 0d;
            double s2 = 0d;
            double s3 = 0d;
            double previousC = int.MinValue;
            double QofC = 0d;
            double c = 10d; // shape
            double b = 0d; // scale

            // solve for the shape parameter
            while (Math.Abs(c - previousC) >= 0.0001d)
            {
                s1 = 0d;
                s2 = 0d;
                s3 = 0d;
                foreach (double x in samples)
                {
                    if (x > 0d)
                    {
                        s1 += Math.Log(x);
                        s2 += Math.Pow(x, c);
                        s3 += Math.Pow(x, c) * Math.Log(x);
                    }
                }

                QofC = n * s2 / (n * s3 - s1 * s2);
                previousC = c;
                c = (c + QofC) / 2d;
            }

            // solve for scale
            foreach (double x in samples)
            {
                if (x > 0d)
                {
                    b += Math.Pow(x, c);
                }
            }

            b = Math.Pow(b / n, 1d / c);

            // return parameters
            return [b, c];
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Lambda, Kappa, true);
            if (x < Minimum) return 0.0d;
            if (x == 0.0d && Kappa == 1.0d)
            {
                return Kappa / Lambda;
            }
            else
            {
                return Kappa / Lambda * Math.Pow(x / Lambda, Kappa - 1.0d) * Math.Exp(-Math.Pow(x / Lambda, Kappa));
            }
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Lambda, Kappa, true);
            if (x < Minimum)
                return 0d;
            return 1d - Math.Exp(-Math.Pow(x / Lambda, Kappa));
        }

        /// <inheritdoc/>
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
                ValidateParameters(Lambda, Kappa, true);
            // Compute the inverse CDF
            return Lambda * Math.Pow(Math.Log(1d / (1d - probability)), 1d / Kappa);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Weibull(Lambda, Kappa);
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_lambda, _kappa, true);

            // Compute covariance
            double a = Lambda;
            double b = Kappa;
            var covar = new double[2, 2];
            covar[0, 0] = 1.108665d * a * a / (sampleSize * b * b); // scale
            covar[1, 1] = 0.607927d * b * b / sampleSize; // shape
            covar[0, 1] = 0.257022d * a / sampleSize;
            covar[1, 0] = covar[0, 1];
            return covar;
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            var covar = ParameterCovariance(sampleSize, estimationMethod);
            var grad = QuantileGradient(probability);
            double varA = covar[0, 0];
            double varB = covar[1, 1];
            double covAB = covar[1, 0];
            double dQx1 = grad[0];
            double dQx2 = grad[1];
            return Math.Pow(dQx1, 2d) * varA + Math.Pow(dQx2, 2d) * varB + 2d * dQx1 * dQx2 * covAB;
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_lambda, _kappa, true);
            double a = Lambda;
            double b = Kappa;
            var gradient = new double[]
            {
                Math.Log(Math.Pow(1d / (1d - probability), 1d / Kappa)), // scale
                a * Math.Log(1d - probability) / (b * b) // shape
            };
            return gradient;
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            if (probabilities.Count != NumberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(probabilities), "The number of probabilities must be the same length as the number of distribution parameters.");
            }
            // Get gradients
            var dQp1 = QuantileGradient(probabilities[0]);
            var dQp2 = QuantileGradient(probabilities[1]);
            // Compute determinant
            // |a b|
            // |c d|
            // |A| = ad − bc
            double a = dQp1[0];
            double b = dQp1[1];
            double c = dQp2[0];
            double d = dQp2[1];
            determinant = a * d - b * c;
            // Return Jacobian
            var jacobian = new double[,] { { a, b }, { c, d } };
            return jacobian;
        }

    }
}