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
using Numerics.Mathematics.RootFinding;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Gumbel (Extreme Value Type I) probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href="https://en.wikipedia.org/wiki/Gumbel_distribution"/>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Gumbel : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Gumbel (Extreme Value Type I) distribution with a location of 100 and scale of 10.
        /// </summary>
        public Gumbel()
        {
            SetParameters(100d, 10d);
        }

        /// <summary>
        /// Constructs a Gumbel (Extreme Value Type I) distribution with a given ξ and α.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public Gumbel(double location, double scale)
        {
            SetParameters(location, scale);
        }

        private double _xi; // location
        private double _alpha; // scale

        /// <summary>
        /// Gets and sets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters([value, Alpha], false) is null;
                _xi = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return _alpha; }
            set
            {
                _parametersValid = ValidateParameters([Xi, value], false) is null;
                _alpha = value;
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
            get { return UnivariateDistributionType.Gumbel; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Gumbel (EVI)"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "EVI"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Location (ξ)";
                parmString[1, 0] = "Scale (α)";
                parmString[0, 1] = Xi.ToString();
                parmString[1, 1] = Alpha.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ξ", "α"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xi), nameof(Alpha)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xi, Alpha]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Xi + Alpha * Tools.Euler; }
        }

        /// <summary>
        /// Gets the median of the distribution.
        /// </summary>
        public override double Median
        {
            get { return Xi - Alpha * Math.Log(Math.Log(2.0d)); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Xi; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(Math.Pow(Math.PI, 2.0d) / 6.0d * Math.Pow(Alpha, 2.0d)); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 1.1396d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + 12.0d / 5.0d; }
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
                SetParameters(ParametersFromMoments(Statistics.ProductMoments(sample)));
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(Statistics.LinearMoments(sample)));
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
            var newDistribution = new Gumbel(Xi, Alpha);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public void SetParameters(double location, double scale)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(location, scale, false) is null;
            // Set parameters
            Xi = location;
            _alpha = scale;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double location, double scale, bool throwException)
        {
            if (double.IsNaN(location) || double.IsInfinity(location))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            // Solve for alpha
            double a = Math.Sqrt(6d) / Math.PI * moments[1];
            // Solve for Xi
            double x = moments[0] - a * Tools.Euler;
            // return parameters
            return [x, a];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new Gumbel();
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
            double L1 = moments[0];
            double L2 = moments[1];
            double alpha = L2 / Tools.Log2;
            double xi = L1 - alpha * Tools.Euler;
            return [xi, alpha];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0];
            double alpha = parameters[1];
            double L1 = xi + alpha * Tools.Euler;
            double L2 = alpha * Math.Log(2.0d);
            double T3 = Math.Log(9d / 8d) / Math.Log(2.0d);
            double T4 = (16.0d * Math.Log(2.0d) - 10.0d * Math.Log(3.0d)) / Math.Log(2.0d);
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // 
            // Get initial values
            // initialVals = DirectMethodOfMoments(Statistics.ComputeProductMoments(sample))
            initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(sample));
            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of scale
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
                var GUM = new Gumbel();
                GUM.SetParameters(x);
                return GUM.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <summary>
        /// This function is used to calculate the maximum likelihood estimates of location and scale parameters.
        /// </summary>
        /// <param name="sample"></param>
        /// <references>
        /// Handbook of Statistical Distributions with Application
        /// </references>
        public void SetParametersFromMLE(IList<double> sample)
        {
            // compute moments from sample
            var moments = Statistics.ProductMoments(sample);
            try
            {

                // get lower and upper bounds
                double LB = 1.0d / (Math.Sqrt(6d) / Math.PI * moments[1]) / 4d;
                double UB = 1.0d / (Math.Sqrt(6d) / Math.PI * moments[1]) * 4d;
                // get first guess
                double theta_hat = 1.0d / (Math.Sqrt(6d) / Math.PI * moments[1]);

                // solve using robust newton-raphson
                double thetan = NewtonRaphson.RobustSolve((t) =>
                {
                    double xb = moments[0];
                    double SumXY = 0d;
                    double SumYY = 0d;
                    for (int i = 0; i < sample.Count; i++)
                    {
                        SumXY += sample[i] * Math.Exp(-sample[i] * t);
                        SumYY += Math.Exp(-sample[i] * t);
                    }

                    return 1.0d / t - xb + SumXY / SumYY;
                }, (t) =>
                {
                    double SumXY = 0d;
                    double SumX2Y = 0d;
                    double SumYY = 0d;
                    for (int i = 0; i < sample.Count; i++)
                    {
                        SumXY += sample[i] * Math.Exp(-sample[i] * t);
                        SumX2Y += sample[i] * sample[i] * Math.Exp(-sample[i] * t);
                        SumYY += Math.Exp(-sample[i] * t);
                    }

                    return Math.Pow(SumXY, 2.0d) / Math.Pow(SumYY, 2.0d) - SumX2Y / SumYY - 1.0d / Math.Pow(t, 2.0d);
                }, theta_hat, LB, UB, 0.000001d, 100, true);
                double SumY = 0d;
                for (int i = 0; i < sample.Count; i++)
                    SumY += Math.Exp(-sample[i] * thetan);
                double AvgY = SumY / sample.Count;

                // set parameters
                double x = -Math.Log(AvgY) / thetan;
                double a = 1.0d / thetan;
                SetParameters(x, a);
            }
            catch (ArgumentException)
            {
                // If the newton method fails to converge, fall back to sample moments
                SetParameters(ParametersFromMoments(moments));
            }
            catch (Exception ex)
            {
                throw ex;
            }
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, true);
            double z = (x - Xi) / Alpha;
            return 1d / Alpha * Math.Exp(-(z + Math.Exp(-z)));
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, true);
            double z = (x - Xi) / Alpha;
            return Math.Exp(-Math.Exp(-z));
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
                ValidateParameters(Xi, _alpha, true);
            return Xi - Alpha * Math.Log(-Math.Log(probability));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Gumbel(Xi, Alpha);
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
                ValidateParameters(Xi, _alpha, true);
            // Compute covariance
            double a = Alpha;
            var covar = new double[2, 2];
            covar[0, 0] = 1.1087d * a * a / sampleSize; // location
            covar[1, 1] = 0.6079d * a * a / sampleSize; // scale
            covar[0, 1] = 0.257d * a * a / sampleSize;
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
                ValidateParameters(Xi, _alpha, true);
            var gradient = new double[]
            {
                1.0d, // location
                -Math.Log(-Math.Log(probability)) // scale
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