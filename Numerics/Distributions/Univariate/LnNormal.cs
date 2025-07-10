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
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Ln-Normal (Galton) probability distribution.
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
    /// Wikipedia contributors, "Log-normal distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Log-normal_distribution"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class LnNormal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
    
        /// <summary>
        /// Constructs a Ln-Normal distribution with a mean of 10 and standard deviation of 10.
        /// </summary>
        public LnNormal()
        {
            SetParameters(10d, 10d);
        }

        /// <summary>
        /// Constructs a Ln-Normal (Galton) distribution with given mean and standard deviation.
        /// </summary>
        /// <param name="mean">The mean of the distribution.</param>
        /// <param name="standardDeviation">The standard deviation of the distribution.</param>
        /// <remarks>
        /// Enter the real-space mean and standard deviation of the distribution. The two parameters μ and σ are not
        /// location and scale parameters for a log-normally distributed random variable X, but they are respectively
        /// location and scale parameters for the normally distributed logarithm ln(X).
        /// </remarks>
        public LnNormal(double mean, double standardDeviation)
        {
            SetParameters(mean, standardDeviation);
        }

        private double _mu;
        private double _sigma;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, false) is null;
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
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, false) is null;
                _sigma = value;
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
            get { return UnivariateDistributionType.LnNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Log-Normal (base e)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LN"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[0, 1] = Mean.ToString();
                parmString[1, 1] = StandardDeviation.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mean), nameof(StandardDeviation)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mean, StandardDeviation]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Math.Exp(Mu + Sigma * Sigma / 2.0d); }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Math.Exp(Mu); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Math.Exp(Mu - Sigma * Sigma); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt((Math.Exp(Sigma * Sigma) - 1.0d) * Math.Exp(2d * Mu + Sigma * Sigma)); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return (Math.Exp(Sigma * Sigma) + 2.0d) * Math.Sqrt(Math.Exp(Sigma * Sigma) - 1d); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                double siqma2 = Math.Pow(Sigma, 2d);
                return 3d + (Math.Exp(4d * siqma2) + 2d * Math.Exp(3d * siqma2) + 3d * Math.Exp(2d * siqma2) - 6d)  ;
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
            get { return [double.NegativeInfinity, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                // Estimate using the method of moments(a.k.a product moments).
                var parms = IndirectMethodOfMoments(sample);
                Mu = parms[0];
                Sigma = parms[1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                var parms = ParametersFromLinearMoments(IndirectMethodOfLinearMoments(sample));
                Mu = parms[0];
                Sigma = parms[1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
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
            // Create a new distribution and estimate parameters from the bootstrap sample 
            var newDistribution = new LnNormal() { Mu = Mu, Sigma = Sigma };
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters using the "direct method."
        /// </summary>
        /// <param name="mean">The mean of the distribution.</param>
        /// <param name="standardDeviation">The standard deviation of the distribution.</param>
        /// <remarks>
        /// The direct method for setting parameters is used so that users can set the parameters directly
        /// from real-space data, which is more intuitive.
        /// </remarks>
        public void SetParameters(double mean, double standardDeviation)
        {
            var parms = DirectMethodOfMoments(mean, standardDeviation);
            // Validate parameters
            Mu = parms[0];
            Sigma = parms[1];
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mean">Mean.</param>
        /// <param name="standardDeviation">Standard deviation.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double mean, double standardDeviation, bool throwException)
        {
            if (double.IsNaN(mean) || double.IsInfinity(mean))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
            }
            if (double.IsNaN(standardDeviation) || double.IsInfinity(standardDeviation) || standardDeviation <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Sigma must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Sigma must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }
       
        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        public static double[] IndirectMethodOfMoments(IList<double> sample)
        {
            // Transform the sample
            var transformedSample = new List<double>();
            for (int i = 0; i < sample.Count; i++)
            {
                if (sample[i] > 0d)
                {
                    transformedSample.Add(Math.Log(sample[i]));
                }
                else
                {
                    transformedSample.Add(Math.Log(0.1d));
                }
            }
            return Statistics.ProductMoments(transformedSample);
        }

        /// <summary>
        /// The indirect method of moments derives the moments from the log transformed data.
        /// This method was proposed by the U.S. Water Resources Council (WRC, 1967).
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        public double[] IndirectMethodOfLinearMoments(IList<double> sample)
        {
            // Transform the sample
            var transformedSample = new List<double>();
            for (int i = 0; i < sample.Count; i++)
            {
                if (sample[i] > 0d)
                {
                    transformedSample.Add(Math.Log(sample[i]));
                }
                else
                {
                    transformedSample.Add(Math.Log(0.1d));
                }
            }
            return Statistics.LinearMoments(transformedSample);
        }

        /// <summary>
        /// Sets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="mean">The real-space mean of the data.</param>
        /// <param name="standardDeviation">The real-space standard deviation of the data.</param>
        public static double[] DirectMethodOfMoments(double mean, double standardDeviation)
        {
            if (standardDeviation <= 0)
                return [double.NaN, double.NaN];
            double variance = Math.Pow(standardDeviation, 2d);
            double mu = Math.Log(Math.Pow(mean, 2d) / Math.Sqrt(variance + Math.Pow(mean, 2d)));
            double sigma = Math.Sqrt(Math.Log(1.0d + variance / Math.Pow(mean, 2d)));
            if (sigma < 1E-16 && Math.Sign(sigma) != -1) sigma = Tools.DoubleMachineEpsilon;
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var mean = moments[0];
            var standardDeviation = moments[1];
            if (standardDeviation <= 0)
                return [double.NaN, double.NaN];
            double variance = Math.Pow(standardDeviation, 2d);
            double mu = Math.Log(Math.Pow(mean, 2d) / Math.Sqrt(variance + Math.Pow(mean, 2d)));
            double sigma = Math.Sqrt(Math.Log(1.0d + variance / Math.Pow(mean, 2d)));
            if (sigma < 1E-16 && Math.Sign(sigma) != -1) sigma = Tools.DoubleMachineEpsilon;
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new LnNormal();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            double mu = moments[0];
            double sigma = moments[1] * Math.Sqrt(Math.PI);
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double L1 = parameters[0];
            double L2 = parameters[1] * Math.Pow(Math.PI, -0.5);
            double T3 = 0d;
            double T4 = 30d * Math.Pow(Math.PI, -1d) * Math.Atan(Tools.Sqrt2) - 9d;
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            var moments = Statistics.ProductMoments(sample);
            initialVals[0] = moments[0];
            initialVals[1] = moments[1];
            // Get bounds of mean
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of standard deviation
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
                var LN = new LnNormal();
                LN.SetParameters(x);
                return LN.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            if (x <= Minimum) return 0.0d;
            double d = (Math.Log(x) - Mu) / Sigma;
            return Math.Exp(-0.5d * d * d) / (Tools.Sqrt2PI * Sigma * x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            if (x <= Minimum)
                return 0d;
            return 0.5d * (1.0d + Erf.Function((Math.Log(x) - Mu) / (Sigma * Math.Sqrt(2.0d))));
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
                ValidateParameters(Mu, Sigma, true);
            return Math.Exp(Mu - Sigma * Math.Sqrt(2.0d) * Erf.InverseErfc(2.0d * probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new LnNormal() { Mu = Mu, Sigma = Sigma };
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments &&
                estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Compute covariance
            double u2 = Sigma;
            var covar = new double[2, 2];
            covar[0, 0] = Math.Pow(u2, 2d) / sampleSize; // location
            covar[1, 1] = 2d * Math.Pow(u2, 4d) / sampleSize; // scale
            covar[0, 1] = 0.0;
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
                ValidateParameters(Mu, _sigma, true);
            double u1 = Mu;
            double u2 = Sigma;
            double z = Normal.StandardZ(probability);
            var gradient = new double[]
            {
                Math.Exp(u1 + z * u2), // location
                z * Math.Exp(u1 + z * u2) / (2d * u2) // scale
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