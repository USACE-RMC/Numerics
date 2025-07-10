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
    /// The generalized Pareto distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Generalized_Pareto_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class GeneralizedPareto : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs an Generalized Pareto distribution with a location of 100, scale of 10, and shape of 0.
        /// </summary>
        public GeneralizedPareto()
        {
            SetParameters(new[] { 100d, 10d, 0d });
        }

        /// <summary>
        /// Constructs a Generalized Pareto (GPA) distribution with the given parameters ξ, α, and κ.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GeneralizedPareto(double location, double scale, double shape)
        {
            SetParameters(location, scale, shape);
        }

        private double _xi; // location
        private double _alpha; // scale
        private double _kappa; // shape
        private double _lambda = 1d;

        /// <summary>
        /// Gets and sets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return _xi; }
            set
            {
                _parametersValid = ValidateParameters([value, Alpha, Kappa], false) is null;
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
                _parametersValid = ValidateParameters(Xi, value, Kappa, false) is null;
                _alpha = value;
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
                _parametersValid = ValidateParameters([Xi, Alpha, value], false) is null;
                _kappa = value;
            }
        }

        /// <summary>
        /// Gets and sets the average number of peaks per block.
        /// </summary>
        public double Lambda
        {
            get { return _lambda; }
            set { _lambda = value; }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.GeneralizedPareto; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Pareto"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GPA"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Location (ξ)";
                parmString[1, 0] = "Scale (α)";
                parmString[2, 0] = "Shape (κ)";
                parmString[0, 1] = Xi.ToString();
                parmString[1, 1] = Alpha.ToString();
                parmString[2, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["ξ", "α", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Xi), nameof(Alpha), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Xi, Alpha, Kappa]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return Xi + Alpha;
                }
                else if (Math.Abs(Kappa) < 1d)
                {
                    return Xi + Alpha / (1d + Kappa);
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
            get
            {
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return Xi - Math.Log(0.5d) * Alpha;
                }
                else
                {
                    return Xi + Alpha * (Math.Pow(2.0d, -Kappa) - 1d) / Kappa;
                }
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return Xi;
                }
                else
                {
                    return Xi + Alpha * (Math.Pow(1d + Kappa, -Kappa) - 1d) / Kappa;
                }
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return Alpha;
                }
                else if (Math.Abs(Kappa) < 0.5d)
                {
                    return Math.Sqrt(Math.Pow(Alpha, 2d) / ((1d + 2d * Kappa) * Math.Pow(1d + Kappa, 2d)));
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
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return 2.0d;
                }
                else if (Math.Abs(Kappa) < 1d / 3d)
                {
                    double num = 2d * (1d - Kappa) * Math.Sqrt(1d + 2d * Kappa);
                    double den = 1d + 3d * Kappa;
                    return num / den;
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
                if (Math.Abs(Kappa) <= NearZero)
                {
                    return 9.0d;
                }
                else if (Math.Abs(Kappa) < 0.25d)
                {
                    double num = 3d * (1d + 2d * Kappa) * (3d - Kappa + 2d * Math.Pow(Kappa, 2d));
                    double den = (1d + 3d * Kappa) * (1d + 4d * Kappa);
                    return num / den;
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
            get { return Xi; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                if (Kappa <= NearZero)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return Xi + Alpha / Kappa;
                }
            }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d, double.NegativeInfinity]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity]; }
        }


        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(DirectMethodOfMoments(Statistics.ProductMoments(sample)));
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
            var newDistribution = new GeneralizedPareto(Xi, Alpha, Kappa);
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
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public void SetParameters(double location, double scale, double shape)
        {
            _parametersValid = ValidateParameters(location, scale, shape, false) is null;
            Xi = location;
            _alpha = scale;
            Kappa = shape;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double location, double scale, double shape, bool throwException)
        {
            if (double.IsNaN(location) || double.IsInfinity(location))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException) throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            if (double.IsNaN(shape) || double.IsInfinity(shape))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be a number.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        
        /// <summary>
        /// Gets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="moments">The array of sample moments.</param>
        public double[] DirectMethodOfMoments(IList<double> moments)
        {

            // Solve for kappa
            double k = SolveForKappa(moments[2]);
            double a;
            double x;
            if (Math.Abs(k) <= NearZero)
            {
                x = moments[0] - moments[1];
                a = moments[1];
            }
            else
            {
                a = Math.Sqrt(moments[1] * moments[1] * Math.Pow(1d + k, 2d) * (1d + 2d * k));
                x = moments[0] - a / (1d + k);
            }

            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var parms = new double[NumberOfParameters];
            parms[0] = 1d / (moments[0] / Math.Pow(moments[1], 2d));
            parms[1] = Math.Pow(moments[0], 2d) / Math.Pow(moments[1], 2d);
            return parms;
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new GeneralizedPareto();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <summary>
        /// Estimate parameters using the modified method of moments.
        /// </summary>
        /// <param name="sample">The array of sample data.</param>
        public double[] ModifiedMethodOfMoments(IList<double> sample)
        {
            int N = sample.Count;
            var moments = Statistics.ProductMoments(sample);
            double min = Statistics.Minimum(sample);
            double m1 = moments[0];
            double m2 = Math.Pow(moments[1], 2d);
            double b = (N - 1) * m2 / (m1 - min) - m1;
            double c = m1 * m1 - m2 + 2d * m2 * (m1 - N * min) / (m1 - min);
            double x = -b + Math.Sqrt(b * b - c);
            double k = 0.5d * (Math.Pow(m1 - x, 2d) / m2 - 1d);
            double a = 0.5d * ((m1 - x) * (Math.Pow(m1 - x, 2d) / m2 + 1d));
            // return parameters
            return [x, a, k];
        }


        /// <summary>
        /// Solve for the shape parameter κ (kappa) given the skewness coefficient.
        /// </summary>
        /// <param name="skew">The skewness coefficient</param>
        /// <returns>
        /// Kappa
        /// </returns>
        public double SolveForKappa(double skew)
        {
            if (Math.Abs(skew) < 10d)
            {
                // Kappa must be solved for. The Brent method is used here. 
                return Brent.Solve((x) =>
                {
                    double k = 2d * (1d - x) * Math.Sqrt(1d + 2d * x) / (1d + 3d * x);
                    return k - skew;
                }, -(1d / 3d), 1d / 3d);
            }
            else
            {
                return double.NaN;
            }
        }

        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            double kappa = (1.0d - 3.0d * T3) / (1.0d + T3);
            double alpha = (1.0d + kappa) * (2.0d + kappa) * L2;
            double xi = L1 - (2.0d + kappa) * L2;
            return [xi, alpha, kappa];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            if (kappa <= -1.0d)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "L-moments can only be defined for kappa > -1.");
            double L1 = xi + alpha / (1.0d + kappa);
            double L2 = alpha / ((1.0d + kappa) * (2.0d + kappa));
            double T3 = (1.0d - kappa) / (3.0d + kappa);
            double T4 = (1.0d - kappa) * (2.0d - kappa) / ((3.0d + kappa) * (4.0d + kappa));
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
            initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(sample));
            double minData = Statistics.Minimum(sample);
            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = initialVals[0] - Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0]))));
            upperVals[0] = minData + Tools.DoubleMachineEpsilon;

            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[1])) + 1d));
            // Get bounds of shape
            lowerVals[2] = -10d;
            upperVals[2] = 10d;
            // Correct initial values if necessary
            if (initialVals[0] <= lowerVals[0] || initialVals[0] >= upperVals[0])
            {
                initialVals[0] = Statistics.Mean(new[] { lowerVals[0], upperVals[0] });
            }
            if (initialVals[1] <= lowerVals[1] || initialVals[1] >= upperVals[1])
            {
                initialVals[1] = Statistics.Mean(new[] { lowerVals[1], upperVals[1] });
            }
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
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

            var xi = Statistics.Minimum(sample);

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var GPA = new GeneralizedPareto();
                GPA.SetParameters(new[] { xi, x[0], x[1] });
                //GPA.SetParameters(x);
                return GPA.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, 2, Initials.Subset(1), Lowers.Subset(1), Uppers.Subset(1));
            solver.ReportFailure = true;
            solver.Maximize();
            return [xi, solver.BestParameterSet.Values[0], solver.BestParameterSet.Values[1]];
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x < Minimum || x > Maximum) return 0.0d;
            double y = (x - Xi) / Alpha;
            if (Math.Abs(Kappa) > NearZero)
                y = -Math.Log(1d - Kappa * y) / Kappa;
            return Math.Exp(-(1d - Kappa) * y) / Alpha;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            double y = (x - Xi) / Alpha;
            if (Math.Abs(Kappa) > NearZero)
                y = -Math.Log(1d - Kappa * y) / Kappa;
            return 1d - Math.Exp(-y);
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
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (Math.Abs(Kappa) <= NearZero)
            {
                return Xi - Alpha * Math.Log(1d - probability);
            }
            else
            {
                return Xi + Alpha / Kappa * (1d - Math.Pow(1d - probability, Kappa));
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedPareto(Xi, Alpha, Kappa) { Lambda = Lambda };
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
                ValidateParameters(new[] { Xi, _alpha }, true);
            double a = Alpha;
            double k = Kappa;
            int N = sampleSize;
            var covar = new double[3, 3];
            covar[0, 0] = N * a * a / ((N + 2d * k) * Math.Pow(N + k, 2d)); // location
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                double num = Math.Pow(1d + k, 2d) * (1d + 6d * k + 12d * Math.Pow(k, 2d));
                double den = (1d + 2d * k) * (1d + 3d * k) * (1d + 4d * k);
                covar[1, 1] = 2d * a * a / N * num / den; // scale
                // 
                num = Math.Pow(1d + k, 2d) * Math.Pow(1d + 2d * k, 2d) * (1d + k + 6d * Math.Pow(k, 2d));
                covar[2, 2] = 1d / N * num / den; // shape
                //
                covar[0, 1] = 0.0;
                covar[0, 2] = 0.0;
                covar[1, 0] = 0.0;
                covar[2, 0] = 0.0;
                //
                num = Math.Pow(1d + k, 2d) * (1d + 2d * k) * (1d + 4d * k + 12d * Math.Pow(k, 2d));
                den = (1d + 2d * k) * (1d + 3d * k) * (1d + 4d * k);
                covar[2, 1] = a / N * num / den; // scale & shape
                covar[1, 2] = covar[2, 1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                covar[1, 1] = (1d - k) * (2d * a * a) / N; // scale
                covar[2, 2] = 1d / N * Math.Pow(1d - k, 2d); // shape
                //
                covar[0, 1] = 0.0;
                covar[0, 2] = 0.0;
                covar[1, 0] = 0.0;
                covar[2, 0] = 0.0;
                //
                covar[2, 1] = a / N * (1d - k); // scale & shape
                covar[1, 2] = covar[2, 1];
            }

            return covar;
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            // The variance in the location parameter is of order N-2 and thus
            // can be neglected relative to the variances of scale and shape (order N-1)
            // to obtain approximate confidence intervals. 
            var covar = ParameterCovariance(sampleSize, estimationMethod);
            var grad = QuantileGradient(probability);
            double varA = covar[1, 1];
            double varB = covar[2, 2];
            double covAB = covar[2, 1];
            double dQx1 = grad[1];
            double dQx2 = grad[2];
            return Math.Pow(dQx1, 2d) * varA + Math.Pow(dQx2, 2d) * varB + 2d * dQx1 * dQx2 * covAB;
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            double p = probability;
            double a = Alpha;
            double k = Kappa;
            var gradient = new double[]
            {
                1.0d, // location
                1d / k * (1d - Math.Pow(1d - p, Kappa)), // scale
                -a / (k * k) * (1d - Math.Pow(1d - p, k)) - a / k * Math.Log(1d - p) * Math.Pow(1d - p, k) // shape
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
            var dQp3 = QuantileGradient(probabilities[2]);
            // Compute determinant
            // |a b c|
            // |d e f|
            // |g h i|
            // |A| = a(ei − fh) − b(di − fg) + c(dh − eg)
            double a = dQp1[0];
            double b = dQp1[1];
            double c = dQp1[2];
            double d = dQp2[0];
            double e = dQp2[1];
            double f = dQp2[2];
            double g = dQp3[0];
            double h = dQp3[1];
            double i = dQp3[2];
            determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
            // Return Jacobian
            var jacobian = new double[,] { { a, b, c }, { d, e, f }, { g, h, i } };
            return jacobian;
        }

    }
}