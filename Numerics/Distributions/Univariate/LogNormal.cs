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
using System.Threading.Tasks;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Log-Normal probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class LogNormal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Log-Normal distribution with a mean (of log) of 3 and standard deviation (of log) of 0.5
        /// </summary>
        public LogNormal()
        {
            SetParameters(3d, 0.5d);
        }

        /// <summary>
        /// Constructs a Log-Normal distribution with given mean (of log) and standard deviation (of log).
        /// </summary>
        /// <param name="meanOfLog">The mean of the log transformed data.</param>
        /// <param name="standardDeviationOfLog">The standard deviation of the log transformed data.</param>
        public LogNormal(double meanOfLog, double standardDeviationOfLog)
        {
            SetParameters(meanOfLog, standardDeviationOfLog);
        }

        // Private variables
        private double _mu;
        private double _sigma;
        private double _base = 10d;
        private bool _momentsComputed = false;
        private double[] u = [double.NaN, double.NaN, double.NaN, double.NaN];

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
                _momentsComputed = false;
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
                _momentsComputed = false;
            }
        }

        /// <summary>
        /// Gets and sets the base of the logarithm
        /// </summary>
        public double Base
        {
            get { return _base; }
            set
            {
                if (value < 1d)
                {
                    _base = 1d;
                }
                else
                {
                    _base = value;
                }
                _momentsComputed = false;
            }
        }

        /// <summary>
        /// Gets the log correction factor.
        /// </summary>
        private double K
        {
            get { return 1d / Math.Log(Base); }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.LogNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Log-Normal (base 10)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "LogN"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean (of log) (µ)";
                parmString[1, 0] = "Std Dev (of log) (σ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
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
            get { return [nameof(Mu), nameof(Sigma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[0];
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <summary>
        /// Gets the mode of the distribution.
        /// </summary>
        public override double Mode
        {
            get
            {    
                return Math.Exp(Mu / K);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[1];
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[2];
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = CentralMoments(1000);
                    _momentsComputed = true;
                }
                return u[3];
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
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(IndirectMethodOfMoments(sample));
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(IndirectMethodOfLinearMoments(sample)));
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
            var newDistribution = new LogNormal(Mu, Sigma);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters based on the moments of the log transformed data.
        /// </summary>
        /// <param name="meanOfLog">The mean of the log transformed data.</param>
        /// <param name="standardDeviationOfLog">The standard deviation of the log transformed data.</param>
        public void SetParameters(double meanOfLog, double standardDeviationOfLog)
        {
            // Set parameters
            Mu = meanOfLog;
            Sigma = standardDeviationOfLog;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mu">The mean (of log).</param>
        /// <param name="sigma">The standard deviation (of log).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double mu, double sigma, bool throwException)
        {
            if (double.IsNaN(mu) || double.IsInfinity(mu))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "Mu must be a number.");
            }
            if (double.IsNaN(sigma) || double.IsInfinity(sigma) || sigma <= 0.0d)
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
        public double[] IndirectMethodOfMoments(IList<double> sample)
        {
            // Transform the sample
            var transformedSample = new List<double>();
            for (int i = 0; i < sample.Count; i++)
            {
                if (sample[i] > 0d)
                {
                    transformedSample.Add(Math.Log(sample[i], Base));
                }
                else
                {
                    transformedSample.Add(Math.Log(0.1d, Base));
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
                    transformedSample.Add(Math.Log(sample[i], Base));
                }
                else
                {
                    transformedSample.Add(Math.Log(0.1d, Base));
                }
            }
            return Statistics.LinearMoments(transformedSample);
        }

        /// <summary>
        /// Sets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="mean">The real-space mean of the data.</param>
        /// <param name="standardDeviation">The real-space standard deviation of the data.</param>
        public double[] DirectMethodOfMoments(double mean, double standardDeviation)
        {
            double variance = Math.Pow(standardDeviation, 2d);
            double mu = Math.Log(Math.Pow(mean, 2d) / Math.Sqrt(variance + Math.Pow(mean, 2d)), Base);
            double sigma = Math.Sqrt(Math.Log(1.0d + variance / Math.Pow(mean, 2d), Base));
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var mean = moments[0];
            var standardDeviation = moments[1];
            double variance = Math.Pow(standardDeviation, 2d);
            double mu = Math.Log(Math.Pow(mean, 2d) / Math.Sqrt(variance + Math.Pow(mean, 2d)), Base);
            double sigma = Math.Sqrt(Math.Log(1.0d + variance / Math.Pow(mean, 2d), Base));
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new LogNormal();
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
            // Estimate initial values using the method of moments (a.k.a product moments).
            var mom = IndirectMethodOfMoments(sample);
            initialVals = new double[] { mom[0], mom[1] };
            // Get bounds of mean
            double real = Math.Exp(initialVals[0] / K);
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Ceiling(Math.Log(Math.Pow(10d, Math.Ceiling(Math.Log10(real) + 1d)), Base));
            upperVals[0] = double.IsNaN(upperVals[0]) ? 5 : upperVals[0];
            // Get bounds of standard deviation
            real = Math.Exp(initialVals[1] / K);
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Ceiling(Math.Log(Math.Pow(10d, Math.Ceiling(Math.Log10(real) + 1d)), Base));
            upperVals[1] = double.IsNaN(upperVals[1]) ? 4 : upperVals[1];
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
                var LogN = new LogNormal() { Base = Base };
                LogN.SetParameters(x);
                return LogN.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;

        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            if (x <= Minimum) return 0.0d;
            double d = (Math.Log(x, Base) - Mu) / Sigma;
            return Math.Exp(-0.5d * d * d) / (Tools.Sqrt2PI * Sigma) * (K / x);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            if (x <= Minimum)
                return 0d;
            return 0.5d * (1.0d + Erf.Function((Math.Log(x, Base) - Mu) / (Sigma * Math.Sqrt(2.0d))));
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
            return Math.Exp((Mu - Sigma * Math.Sqrt(2.0d) * Erf.InverseErfc(2.0d * probability)) / K);
        }

        /// <summary>
        /// Get confidence intervals using Monte Carlo simulation.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="realizations">The number of Monte Carlo realizations.</param>
        /// <param name="quantiles">List of exceedance probabilities for output frequency curves.</param>
        /// <param name="percentiles">List of confidence percentiles for confidence interval output.</param>
        /// <remarks>
        /// This is the same sampling approach as used in HEC-FDA.
        /// </remarks>
        public double[,] MonteCarloConfidenceIntervals(int sampleSize, int realizations, IList<double> quantiles, IList<double> percentiles)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Dimension output array
            int q = quantiles.Count;
            int p = percentiles.Count;
            var Output = new double[q, p];

            // Variables
            double OriginalMean = Mu;
            double OriginalStdDev = Sigma;

            // Create random numbers for mean and standard deviation
            var r = new MersenneTwister(12345);
            var rndMean = r.NextDoubles(realizations);
            r = new MersenneTwister(45678);
            var rndStdDev = r.NextDoubles(realizations);


            // Create list of Monte Carlo distributions
            var MonteCarloDistributions = new UnivariateDistributionBase[realizations];
            Parallel.For(0, realizations, idx =>
            {

                // Generate new mean
                var Normal = new Normal(OriginalMean, OriginalStdDev / Math.Sqrt(sampleSize));
                double NewMu = Normal.InverseCDF(rndMean[idx]);
                // Generate new standard deviation
                var Chi = new ChiSquared(sampleSize - 1);
                double NewSigma = Math.Sqrt((sampleSize - 1) * Math.Pow(OriginalStdDev, 2d) / Chi.InverseCDF(rndStdDev[idx]));
                // Create a new distribution with the new parameters
                MonteCarloDistributions[idx] = new LogNormal(NewMu, NewSigma);
            });

            // Create confidence intervals
            for (int i = 0; i < q; i++)
            {
                // Create array of X values across user-defined probabilities
                var XValues = new double[realizations];
                // Record X values
                Parallel.For(0, realizations, idx => XValues[idx] = MonteCarloDistributions[idx].InverseCDF(1d - quantiles[i]));
                // Record percentiles for user-defined probabilities
                for (int j = 0; j < p;  j++)
                    Output[i, j] = Statistics.Percentile(XValues, percentiles[j]);
            }

            // Return confidence percentile output
            return Output;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new LogNormal(Mu, Sigma);
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
            double varQ = Math.Pow(dQx1, 2d) * varA + Math.Pow(dQx2, 2d) * varB + 2d * dQx1 * dQx2 * covAB;
            return varQ * Math.Pow(InverseCDF(probability) / K, 2d);
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            double z = Normal.StandardZ(probability);
            var gradient = new double[]
            {
                1.0d, // location
                z / (2d * Sigma) // scale
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
            double p0 = InverseCDF(probabilities[0]) / K;
            double p1 = InverseCDF(probabilities[1]) / K;
            double a = dQp1[0] * p0;
            double b = dQp1[1] * p0;
            double c = dQp2[0] * p1;
            double d = dQp2[1] * p1;
            determinant = a * d - b * c;
            // Return Jacobian
            var jacobian = new double[,] { { a, b }, { c, d } };
            return jacobian;
        }
    }
}