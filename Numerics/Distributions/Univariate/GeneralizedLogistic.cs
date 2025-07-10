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
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The generalized logistic distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Generalized_logistic_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class GeneralizedLogistic : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Generalized Logistic distribution with a location of 100, scale of 10, and shape of 0.
        /// </summary>
        public GeneralizedLogistic()
        {
            SetParameters(100d, 10d, 0d);
        }

        /// <summary>
        /// Constructs a Generalized Logistic distribution with the given parameters ξ, α, and κ.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GeneralizedLogistic(double location, double scale, double shape)
        {
            SetParameters(location, scale, shape);
        }

        private double _xi; // location
        private double _alpha; // scale
        private double _kappa; // shape

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

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.GeneralizedLogistic; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Logistic"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GLO"; }
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
                    return Xi;
                }
                else if (Math.Abs(Kappa) < 1d)
                {
                    double U1 = Gamma.Function(1d + Kappa) * Gamma.Function(1d - Kappa);
                    return Xi + Alpha / Kappa * (1d - U1);
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
            get { return InverseCDF(0.5d); }
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
                    return Alpha * Math.PI / Math.Sqrt(3d);
                }
                else if (Math.Abs(Kappa) < 0.5d)
                {
                    double U1 = Gamma.Function(1d + Kappa) * Gamma.Function(1d - Kappa);
                    double U2 = Gamma.Function(1d + 2d * Kappa) * Gamma.Function(1d - 2d * Kappa);
                    return Math.Sqrt(Math.Pow(Alpha, 2d) / Math.Pow(Kappa, 2d) * (U2 - Math.Pow(U1, 2d)));
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
                    return 0.0d;
                }
                else if (Math.Abs(Kappa) < 1d / 3d)
                {
                    double U1 = Gamma.Function(1d + Kappa) * Gamma.Function(1d - Kappa);
                    double U2 = Gamma.Function(1d + 2d * Kappa) * Gamma.Function(1d - 2d * Kappa);
                    double U3 = Gamma.Function(1d + 3d * Kappa) * Gamma.Function(1d - 3d * Kappa);
                    return Math.Sign(Kappa) * (-U3 + 3d * U1 * U2 - 2d * Math.Pow(U1, 3d)) / Math.Pow(U2 - Math.Pow(U1, 2d), 3d / 2d);
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
                    return 3d + 6d / 5d;
                }
                else if (Math.Abs(Kappa) < 0.25d)
                {
                    double U1 = Gamma.Function(1d + Kappa) * Gamma.Function(1d - Kappa);
                    double U2 = Gamma.Function(1d + 2d * Kappa) * Gamma.Function(1d - 2d * Kappa);
                    double U3 = Gamma.Function(1d + 3d * Kappa) * Gamma.Function(1d - 3d * Kappa);
                    double U4 = Gamma.Function(1d + 4d * Kappa) * Gamma.Function(1d - 4d * Kappa);
                    double kNum = U4 - 4d * U3 * U1 - 3d * Math.Pow(U2, 2d) + 12d * U2 * Math.Pow(U1, 2d) - 6d * Math.Pow(U1, 4d);
                    double kDen = Math.Pow(U2 - Math.Pow(U1, 2d), 2d);
                    return 3 + kNum / kDen;
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
            get
            {
                if (Kappa >= -NearZero)
                {
                    return double.NegativeInfinity;
                }
                else
                {
                    return Xi + Alpha / Kappa;
                }
            }
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
            var newDistribution = new GeneralizedLogistic(Xi, Alpha, Kappa);
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
            // Validate parameters
            _parametersValid = ValidateParameters(location, scale, shape, false) is null;
            // Set parameters
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
                x = moments[0];
                a = moments[1] * Math.Sqrt(3d) / Math.PI;
            }
            else
            {
                double U1 = Gamma.Function(1d + k) * Gamma.Function(1d - k);
                double U2 = Gamma.Function(1d + 2d * k) * Gamma.Function(1d - 2d * k);
                a = Math.Sqrt(moments[1] * moments[1] * k * k / (U2 - Math.Pow(U1, 2d)));
                x = moments[0] - a / k * (1d - U1);
            }
            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new GeneralizedLogistic();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
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
                    double U1 = Gamma.Lanczos(1d + x) * Gamma.Lanczos(1d - x);
                    double U2 = Gamma.Lanczos(1d + 2d * x) * Gamma.Lanczos(1d - 2d * x);
                    double U3 = Gamma.Lanczos(1d + 3d * x) * Gamma.Lanczos(1d - 3d * x);
                    double k = Math.Sign(x) * (-U3 + 3d * U1 * U2 - 2d * Math.Pow(U1, 3d)) / Math.Pow(U2 - Math.Pow(U1, 2d), 3d / 2d);
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
            double T4 = moments[3];
            double kappa = -T3;
            double alpha = L2 * Math.Sin(kappa * Math.PI) / (kappa * Math.PI);
            double xi = L1 - alpha * (1.0d / kappa - Math.PI / Math.Sin(kappa * Math.PI));
            return [xi, alpha, kappa];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            if (Math.Abs(kappa) >= 1.0d)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "L-moments can only be defined for -1 < kappa < 1.");
            double L1 = xi + alpha * (1.0d / kappa - Math.PI / Math.Sin(kappa * Math.PI));
            double L2 = alpha * kappa * Math.PI / Math.Sin(kappa * Math.PI);
            double T3 = -kappa;
            double T4 = (1.0d + 5.0d * Math.Pow(kappa, 2.0d)) / 6.0d;
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            // Estimate initial values using the method of moments (a.k.a product moments).
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            // initialVals = DirectMethodOfMoments(Statistics.ComputeProductMoments(sample))
            initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(sample));
            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[1]))) + 1d);
            // Get bounds of shape
            lowerVals[2] = -10;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
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

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var GLO = new GeneralizedLogistic();
                GLO.SetParameters(x);
                return GLO.LogLikelihood(sample);
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
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x < Minimum || x > Maximum) return 0.0d;
            double y = (x - Xi) / Alpha;
            if (Math.Abs(Kappa) > NearZero)
                y = -Math.Log(1d - Kappa * y) / Kappa;
            return 1d / Alpha * Math.Exp(-(1d - Kappa) * y) / Math.Pow(1d + Math.Exp(-y), 2d);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            double y = (x - Xi) / Alpha;
            if (Math.Abs(Kappa) > NearZero)
                y = -Math.Log(1d - Kappa * y) / Kappa;
            return 1d / (1d + Math.Exp(-y));
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
                return Xi - Alpha * Math.Log((1d - probability) / probability);
            }
            else
            {
                return Xi + Alpha / Kappa * (1d - Math.Pow((1d - probability) / probability, Kappa));
            }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedLogistic(Xi, Alpha, Kappa);
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            throw new NotImplementedException();
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            double a = Alpha;
            double k = Kappa;
            var gradient = new double[3];
            gradient[0] = 1.0d; // location
            gradient[1] = 1d / k * (1d - Math.Pow((1d - probability) / probability, k)); // scale
            gradient[2] = -(a / (k * k)) * (1d - Math.Pow((1d - probability) / probability, k)) - a / k * Math.Pow((1d - probability) / probability, k) * Math.Log((1d - probability) / probability); // shape
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