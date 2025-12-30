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
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Pearson Type III distribution.
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
    /// <see href = "http://mathworld.wolfram.com/PearsonTypeIIIDistribution.html" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class PearsonTypeIII : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Pearson Type III distribution with a mean of 100, standard deviation of 10, and skew = 0.
        /// </summary>
        public PearsonTypeIII()
        {
            SetParameters(100d, 10d, 0d);
        }

        /// <summary>
        /// Constructs a Pearson Type III distribution with the given moments (of log) µ, σ, and γ.
        /// </summary>
        /// <param name="mean">The mean of the data.</param>
        /// <param name="standardDeviation">The standard deviation of the data.</param>
        /// <param name="skew">The skew of the data.</param>
        public PearsonTypeIII(double mean, double standardDeviation, double skew)
        {
            SetParameters(mean, standardDeviation, skew);
        }

        private double _mu;
        private double _sigma;
        private double _gamma;
        
        /// <summary>
        /// Gets and sets the Mean of the distribution.
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, Gamma, false) is null;
                SetParameters(value, Sigma, Gamma);
            }
        }

        /// <summary>
        /// Gets and sets the Standard Deviation of the distribution.
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, Gamma, false) is null;
                SetParameters(Mu, value, Gamma);
            }
        }

        /// <summary>
        /// Gets and sets the Skew of the distribution.
        /// </summary>
        public double Gamma
        {
            get { return _gamma; }
            set
            {
                _parametersValid = ValidateParameters(Mu, Sigma, value, false) is null;
                SetParameters(Mu, Sigma, value);
            }
        }

        /// <summary>
        /// Gets the location parameter ξ (Xi).
        /// </summary>
        public double Xi
        {
            get { return Mu - 2.0d * Sigma / Gamma; }
        }

        /// <summary>
        /// Gets and sets the scale parameter β (beta).
        /// </summary>
        public double Beta
        {
            get { return 0.5d * Sigma * Gamma; }
        }

        /// <summary>
        /// Gets and sets the shape parameter α (alpha).
        /// </summary>
        public double Alpha
        {
            get { return 4.0d / Math.Pow(Gamma, 2d); }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 3; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.PearsonTypeIII; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Pearson Type III"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "PIII"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[3, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[2, 0] = "Skew (γ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                parmString[2, 1] = Gamma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ", "γ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma), nameof(Gamma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get {  return [Mu, Sigma, Gamma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                return Mu;
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
                if (Math.Abs(Gamma) <= NearZero)
                {
                    // Use Normal
                    return Mu;
                }
                else
                {
                    return Xi + (Alpha - 1d) * Beta;
                }
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                return Sigma;
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                return Gamma;
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                return 3d + 6d / Alpha;
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                if (Math.Abs(Gamma) <= NearZero)
                {
                    // Use Normal
                    return double.NegativeInfinity;
                }
                else if (Beta > 0d)
                {
                    return Xi;
                }
                else
                {
                    return double.NegativeInfinity;
                }
            }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get
            {
                if (Math.Abs(Skewness) <= NearZero)
                {
                    // Use Normal
                    return double.PositiveInfinity;
                }
                else if (Beta > 0d)
                {
                    return double.PositiveInfinity;
                }
                else
                {
                    return Xi;
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
                SetParameters(Statistics.ProductMoments(sample));
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
            // create a new distribution with the sample 
            var newDistribution = new PearsonTypeIII(Mu, Sigma, Gamma);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters based on the moments of the data.
        /// </summary>
        /// <param name="mean">The mean of the data.</param>
        /// <param name="standardDeviation">The standard deviation of the data.</param>
        /// <param name="skew">The skew of the data.</param>
        public void SetParameters(double mean, double standardDeviation, double skew)
        {
            _parametersValid = ValidateParameters(mean, standardDeviation, skew, false) is null;
            _mu = mean;
            _sigma = standardDeviation;
            _gamma = skew;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mu">The mean of the distribution.</param>
        /// <param name="sigma">The standard deviation of the distribution.</param>
        /// <param name="gamma">The skew of the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double mu, double sigma, double gamma, bool throwException)
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
            if (double.IsNaN(gamma) || double.IsInfinity(gamma))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma must be a number.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma must be a number.");
            }
            if (gamma > 5)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be less than 5.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be less than 5.");
            }
            if (gamma < -5)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be greater than -5.");
                return new ArgumentOutOfRangeException(nameof(Gamma), "Gamma = " + gamma + ". Gamma must be greater than -5.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            return moments.ToArray().Subset(0, 2);
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new PearsonTypeIII();
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
            double T3 = moments[2];
            double alpha = double.NaN;
            double z;
            // The following approximation has relative accuracy better than 5x10-5 for all values of alpha.
            if (Math.Abs(T3) > 0.0d && Math.Abs(T3) < 1d / 3d)
            {
                z = 3.0d * Math.PI * Math.Pow(T3, 2d);
                alpha = (1.0d + 0.2906d * z) / (z + 0.1882d * Math.Pow(z, 2d) + 0.0442d * Math.Pow(z, 3d));
            }
            else if (Math.Abs(T3) >= 1d / 3d && Math.Abs(T3) < 1.0d)
            {
                z = 1.0d - Math.Abs(T3);
                alpha = (0.36067d * z - 0.59567d * Math.Pow(z, 2d) + 0.25361d * Math.Pow(z, 3d)) / (1.0d - 2.78861d * z + 2.56096d * Math.Pow(z, 2d) - 0.77045d * Math.Pow(z, 3d));
            }

            double gamma = 2.0d * Math.Pow(alpha, -0.5d) * Math.Sign(T3);
            double sigma = L2 * Math.Pow(Math.PI, 0.5d) * Math.Pow(alpha, 0.5d) * Mathematics.SpecialFunctions.Gamma.Function(alpha) / Mathematics.SpecialFunctions.Gamma.Function(alpha + 0.5d);
            double mu = L1;
            return [mu, sigma, gamma];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double mu = parameters[0];
            double sigma = parameters[1];
            double gamma = parameters[2];
            double xi = mu - 2.0d * sigma / gamma;
            double alpha = 4.0d / Math.Pow(gamma, 2d);
            double beta = 0.5d * sigma * gamma;
            double L1 = xi + alpha * beta;
            double L2 = Math.Abs(Math.Pow(Math.PI, -0.5d) * beta * Mathematics.SpecialFunctions.Gamma.Function(alpha + 0.5d) / Mathematics.SpecialFunctions.Gamma.Function(alpha));
            // The following approximations are accurate to 10-6. 
            double A0 = 0.32573501d;
            double A1 = 0.1686915d;
            double A2 = 0.078327243d;
            double A3 = -0.0029120539d;
            double B1 = 0.46697102d;
            double B2 = 0.24255406d;
            double C0 = 0.12260172d;
            double C1 = 0.05373013d;
            double C2 = 0.043384378d;
            double C3 = 0.011101277d;
            double D1 = 0.18324466d;
            double D2 = 0.20166036d;
            double E1 = 2.3807576d;
            double E2 = 1.5931792d;
            double E3 = 0.11618371d;
            double F1 = 5.1533299d;
            double F2 = 7.142526d;
            double F3 = 1.9745056d;
            double G1 = 2.1235833d;
            double G2 = 4.1670213d;
            double G3 = 3.1925299d;
            double H1 = 9.0551443d;
            double H2 = 26.649995d;
            double H3 = 26.193668d;
            double T3;
            double T4;
            if (alpha >= 1d)
            {
                T3 = Math.Pow(alpha, -0.5d) * (A0 + A1 * Math.Pow(alpha, -1) + A2 * Math.Pow(alpha, -2) + A3 * Math.Pow(alpha, -3)) / (1d + B1 * Math.Pow(alpha, -1) + B2 * Math.Pow(alpha, -2));
                T4 = (C0 + C1 * Math.Pow(alpha, -1) + C2 * Math.Pow(alpha, -2) + C3 * Math.Pow(alpha, -3)) / (1d + D1 * Math.Pow(alpha, -1) + D2 * Math.Pow(alpha, -2));
            }
            else
            {
                T3 = (1d + E1 * alpha + E2 * Math.Pow(alpha, 2d) + E3 * Math.Pow(alpha, 3d)) / (1d + F1 * alpha + F2 * Math.Pow(alpha, 2d) + F3 * Math.Pow(alpha, 3d));
                T4 = (1d + G1 * alpha + G2 * Math.Pow(alpha, 2d) + G3 * Math.Pow(alpha, 3d)) / (1d + H1 * alpha + H2 * Math.Pow(alpha, 2d) + H3 * Math.Pow(alpha, 3d));
            }

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
            initialVals = moments.Subset(0, 2);
            // Get bounds of mean
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of standard deviation
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            // Get bounds of skew
            lowerVals[2] = -5d;
            upperVals[2] = 5d;

            // Correct initial value of skew if necessary
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0.01;
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
                var P3 = new PearsonTypeIII();
                P3.SetParameters(x);
                return P3.LogLikelihood(sample);
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
                ValidateParameters(Mu, Sigma, Gamma, true);
            if (x < Minimum || x > Maximum) return 0.0d;

            if (Math.Abs(Gamma) <= NearZero)
            {
                // Use Normal distribution
                double z = (x - Mu) / Sigma;
                return Math.Exp(-0.5d * z * z) / (Tools.Sqrt2PI * Sigma);
            }
            // Use Gamma distribution
            if (Beta > 0d)
            {
                double shiftedX = x - Xi;
                return Math.Exp(-shiftedX / Math.Abs(Beta) + (Alpha - 1.0d) * Math.Log(shiftedX) - Alpha * Math.Log(Math.Abs(Beta)) - Mathematics.SpecialFunctions.Gamma.LogGamma(Alpha));
            }
            else
            {
                double shiftedX = Xi - x;
                return Math.Exp(-shiftedX / Math.Abs(Beta) + (Alpha - 1.0d) * Math.Log(shiftedX) - Alpha * Math.Log(Math.Abs(Beta)) - Mathematics.SpecialFunctions.Gamma.LogGamma(Alpha));
            }
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Gamma, true);
            if (x <= Minimum)
                return 0d;
            if (x >= Maximum)
                return 1d;
            if (Math.Abs(Gamma) <= NearZero)
            {
                return 0.5d * (1.0d + Mathematics.SpecialFunctions.Erf.Function((x - Mu) / (Sigma * Math.Sqrt(2.0d))));
            }
            else if (Beta > 0d)
            {
                double shiftedX = x - Xi;
                return Mathematics.SpecialFunctions.Gamma.LowerIncomplete(Alpha, shiftedX / Math.Abs(Beta));
            }
            else
            {
                double shiftedX = Xi - x;
                return 1.0d - Mathematics.SpecialFunctions.Gamma.LowerIncomplete(Alpha, shiftedX / Math.Abs(Beta));
            }
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
                ValidateParameters(Mu, Sigma, Gamma, true);
            if (Math.Abs(Gamma) <= NearZero)
            {
                return Mu + Sigma * Normal.StandardZ(probability);
            }
            else if (Beta > 0d)
            {
                return Xi + Mathematics.SpecialFunctions.Gamma.InverseLowerIncomplete(Alpha, probability) * Math.Abs(Beta);
            }
            else
            {
                return Xi - Mathematics.SpecialFunctions.Gamma.InverseLowerIncomplete(Alpha, 1d - probability) * Math.Abs(Beta);
            }
        }

        /// <summary>
        /// Returns the inverse CDF using the modified Wilson-Hilferty transformation.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <remarks>
        /// Cornish-Fisher transformation (Fisher and Cornish, 1960) for abs(skew) less than or equal to 2. If abs(skew) > 2 then use Modified Wilson-Hilferty transformation (Kirby,1972).
        /// </remarks>
        public double WilsonHilfertyInverseCDF(double probability)
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
                ValidateParameters(Mu, Sigma, Gamma, true);
            // 
            return Mu + Sigma * GammaDistribution.FrequencyFactorKp(Gamma, probability);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new PearsonTypeIII(Mu, Sigma, Gamma);
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
                ValidateParameters(_mu, _sigma, _gamma, true);
            // Compute covariance
            double alpha = 1d / Beta;
            double lambda = Alpha;
            double A = 2d * Mathematics.SpecialFunctions.Gamma.Trigamma(lambda) - 2d / (lambda - 1d) + 1d / Math.Pow(lambda - 1d, 2d);
            var covar = new double[3, 3];
            covar[0, 0] = (lambda - 2d) / (sampleSize * A) * (1d / Math.Pow(alpha, 2d)) * (Mathematics.SpecialFunctions.Gamma.Trigamma(lambda) * lambda - 1d); // location
            covar[1, 1] = (lambda - 2d) * Math.Pow(alpha, 2d) / (sampleSize * A) * (Mathematics.SpecialFunctions.Gamma.Trigamma(lambda) / (lambda - 2d) - 1d / Math.Pow(lambda - 1d, 2d)); // scale
            covar[2, 2] = 2d / (sampleSize * A); // shape
            covar[0, 1] = 1d / sampleSize * ((lambda - 2d) / A) * (Mathematics.SpecialFunctions.Gamma.Trigamma(lambda) - 1d / (lambda - 1d)); // location & scale
            covar[1, 0] = covar[0, 1];
            covar[0, 2] = (2d - lambda) / (sampleSize * alpha * A * (lambda - 1d)); // location & shape
            covar[2, 0] = covar[0, 2];
            covar[1, 2] = alpha / (sampleSize * A * (lambda - 1d)); // scale & shape
            covar[2, 1] = covar[1, 2];
            return covar;
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                double u2 = _sigma;
                double c = _gamma;
                double c2 = Math.Pow(c, 2d);
                int N = sampleSize;
                return Math.Pow(u2, 2d) / N * (1d + Math.Pow(GammaDistribution.FrequencyFactorKp(c, probability), 2d) / 2d * (1d + 0.75d * c2) + GammaDistribution.FrequencyFactorKp(c, probability) * c + 6d * (1d + 0.25d * c2) * GammaDistribution.PartialKp(c, probability) * (GammaDistribution.PartialKp(c, probability) * (1d + 5d * c2 / 4d) + GammaDistribution.FrequencyFactorKp(c, probability) / 2d * c));
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                var covar = ParameterCovariance(sampleSize, estimationMethod);
                var grad = QuantileGradient(probability);
                double varA = covar[0, 0];
                double varB = covar[1, 1];
                double varG = covar[2, 2];
                double covAB = covar[1, 0];
                double covAG = covar[2, 0];
                double covBG = covar[2, 1];
                double dQx1 = grad[0];
                double dQx2 = grad[1];
                double dQx3 = grad[2];
                return Math.Pow(dQx1, 2d) * varA + Math.Pow(dQx2, 2d) * varB + Math.Pow(dQx3, 2d) * varG + 2d * dQx1 * dQx2 * covAB + 2d * dQx1 * dQx3 * covAG + 2d * dQx2 * dQx3 * covBG;
            }

            return double.NaN;
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Gamma, true);
            double alpha = 1d / Beta;
            double lambda = Alpha;
            double eps = Math.Sign(alpha);
            var gradient = new double[]
            {
                1.0d, // location
                -lambda / Math.Pow(alpha, 2d) * (1.0d + eps / Math.Sqrt(lambda) * GammaDistribution.FrequencyFactorKp(Gamma, probability)), // scale
                1.0d / alpha * (1.0d + eps / Math.Sqrt(lambda) * (GammaDistribution.FrequencyFactorKp(Gamma, probability) / 2.0d) - 1.0d / lambda * GammaDistribution.PartialKp(Gamma, probability)) // shape
            };
            return gradient;
        }

        /// <summary>
        /// Returns a list of partial derivatives of X given probability with respect to each moment.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        public double[] QuantileGradientForMoments(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, Gamma, true);
            var gradient = new double[]
            {
                1.0d,
                GammaDistribution.FrequencyFactorKp(Gamma, probability),
                Sigma * GammaDistribution.PartialKp(Gamma, probability)
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
            var dQp1 = QuantileGradientForMoments(probabilities[0]);
            var dQp2 = QuantileGradientForMoments(probabilities[1]);
            var dQp3 = QuantileGradientForMoments(probabilities[2]);
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

        /// <inheritdoc/>
        public override double[] ConditionalMoments(double a, double b)
        {
            if (a >= b)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            var pearson = TryPearsonConditionalMoments(a, b, out var mPearson, 1e-15);
            if (!pearson.stable || !pearson.divisable) 
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };
            return mPearson;
        }

        /// <summary>
        /// Pearson Type III conditional moments with a smooth Normal ↔ WH ↔ Exact blend.
        /// Returns [E[X], m2, m3, m4] where mk are central moments about the **unconditional** mean μ0.
        /// </summary>
        /// <remarks>Returns (stable=false,false) on hard numerical failure; when stable=true but divisable=false, probability mass was too small.</remarks>
        private (bool stable, bool divisable) TryPearsonConditionalMoments(
            double a, double b, out double[] moments, double pMin)
        {
            moments = null!;

            // ---- small math helpers -------------------------------------------------
            static double Phi(double x) => 0.5 * (1.0 + Mathematics.SpecialFunctions.Erf.Function(x / Math.Sqrt(2.0)));
            static double phi(double x) => Math.Exp(-0.5 * x * x) / Math.Sqrt(2.0 * Math.PI);

            // Standardized truncated-normal integrals J_0..J_N with recursion
            static void TruncStdNormalJ(double aStar, double bStar, int N, double[] J)
            {
                double phiA = double.IsNegativeInfinity(aStar) ? 0.0 : phi(aStar);
                double phiB = double.IsPositiveInfinity(bStar) ? 0.0 : phi(bStar);
                double PhiA = double.IsNegativeInfinity(aStar) ? 0.0 : Phi(aStar);
                double PhiB = double.IsPositiveInfinity(bStar) ? 1.0 : Phi(bStar);

                J[0] = PhiB - PhiA;                   // D
                if (N == 0) return;
                J[1] = phiA - phiB;                   // E[T] * D
                for (int n = 2; n <= N; n++)
                {
                    double termA = double.IsNegativeInfinity(aStar) ? 0.0 : Math.Pow(aStar, n - 1) * phiA;
                    double termB = double.IsPositiveInfinity(bStar) ? 0.0 : Math.Pow(bStar, n - 1) * phiB;
                    J[n] = termA - termB + (n - 1) * J[n - 2]; // E[T^n] * D
                }
            }

            // E[(μ+σT)^n | a*<=T<=b*], given standardized J's
            static double EZPow(int n, double mu, double sigma, double[] J)
            {
                double D = J[0];
                if (!(D > 0)) return double.NaN;

                double sum = 0.0;
                for (int j = 0; j <= n; j++)
                {
                    double bin = Mathematics.SpecialFunctions.Factorial.BinomialCoefficient(n, j);
                    double Ej = J[j] / D; // E[T^j | trunc]
                    sum += bin * Math.Pow(mu, n - j) * Math.Pow(sigma, j) * Ej;
                }
                return sum;
            }

            // smooth logistic weights (differentiable)
            static double Sigmoid(double x, double x0, double s) => 1.0 / (1.0 + Math.Exp((x - x0) / s));

            // ---- guard Gamma(alpha) -------------------------------------------------
            double lgAlpha = Mathematics.SpecialFunctions.Gamma.LogGamma(Alpha);
            if (double.IsNaN(lgAlpha) || double.IsInfinity(lgAlpha))
                return (false, false);

            try
            {
                // Map [a,b] in X to [L,U] in Y for unit-scale Gamma Y~Gamma(Alpha,1)
                bool pos = Beta > 0.0;
                double absβ = Math.Abs(Beta);

                double L, U;
                if (pos)
                {
                    L = Math.Max(0.0, (a - Xi) / absβ);
                    U = double.IsPositiveInfinity(b) ? double.PositiveInfinity
                        : Math.Max(0.0, (b - Xi) / absβ);
                }
                else
                {
                    L = Math.Max(0.0, (Xi - b) / absβ);
                    U = double.IsPositiveInfinity(a) ? double.PositiveInfinity
                        : Math.Max(0.0, (Xi - a) / absβ);
                }

                // =====================================================================
                // (A) Exact truncated-Gamma path (numerically stable centralization)
                // =====================================================================
                bool exactOk = true;
                double EX_exact = double.NaN, m2_exact = double.NaN, m3_exact = double.NaN, m4_exact = double.NaN;

                double CdfGamma(int r, double x)
                {
                    double y = (pos ? (x - Xi) : (Xi - x)) / absβ;
                    double P = Mathematics.SpecialFunctions.Gamma.LowerIncomplete(Alpha + r, y);
                    return pos ? P : (1.0 - P);
                }

                double[] EY = new double[5]; // E[Y^r | L<=Y<=U], r=0..4
                {
                    double[] dP = new double[5];
                    for (int r = 0; r <= 4; r++) dP[r] = CdfGamma(r, b) - CdfGamma(r, a);

                    double P0 = dP[0];
                    if (!(P0 > pMin)) exactOk = false;
                    else
                    {
                        double[] gr = new double[5]; // Γ(α+r)/Γ(α)
                        for (int r = 0; r <= 4; r++)
                            gr[r] = Math.Exp(Mathematics.SpecialFunctions.Gamma.LogGamma(Alpha + r) - lgAlpha);

                        for (int r = 0; r <= 4; r++)
                            EY[r] = gr[r] * (dP[r] / P0); // truncated raw moments of Y
                    }
                }

                if (exactOk)
                {
                    double EY1 = EY[1], EY2 = EY[2], EY3 = EY[3], EY4 = EY[4];
                    EX_exact = Xi + Beta * EY1;

                    // conditional central moments of Y about its conditional mean
                    double c2Y = Math.Max(0.0, EY2 - EY1 * EY1);
                    double c3Y = EY3 - 3.0 * EY1 * EY2 + 2.0 * EY1 * EY1 * EY1;
                    double c4Y = EY4 - 4.0 * EY1 * EY3 + 6.0 * EY1 * EY1 * EY2 - 3.0 * Math.Pow(EY1, 4);

                    // shift to moments about α (unconditional mean of Y)
                    double dY = EY1 - Alpha;

                    double b1 = Beta; double b2 = b1 * b1; double b3 = b2 * b1; double b4 = b3 * b1;
                    m2_exact = b2 * (c2Y + dY * dY);
                    m3_exact = b3 * (c3Y + 3.0 * dY * c2Y + dY * dY * dY);
                    m4_exact = b4 * (c4Y + 4.0 * dY * c3Y + 6.0 * dY * dY * c2Y + Math.Pow(dY, 4));
                }

                // =====================================================================
                // (B) Wilson–Hilferty (Gamma ≈ Z^3) truncated-Normal path
                // =====================================================================
                bool whOk = true;
                double EX_wh = double.NaN, m2_wh = double.NaN, m3_wh = double.NaN, m4_wh = double.NaN;

                if (!(Alpha > 0.0) || double.IsInfinity(Alpha) || double.IsNaN(Alpha)) whOk = false;
                else
                {
                    // Z ~ N(muZ, sigZ^2); Y ≈ Alpha * Z^3
                    double muZ = 1.0 - 1.0 / (9.0 * Alpha);
                    double sigZ = Math.Sqrt(1.0 / (9.0 * Alpha));

                    double ell = (L <= 0.0) ? 0.0 : Math.Pow(L / Alpha, 1.0 / 3.0);
                    double uuu = double.IsPositiveInfinity(U) ? double.PositiveInfinity
                               : (U <= 0.0 ? 0.0 : Math.Pow(U / Alpha, 1.0 / 3.0));

                    double aStarZ = double.IsNegativeInfinity(ell) ? double.NegativeInfinity : (ell - muZ) / sigZ;
                    double bStarZ = double.IsPositiveInfinity(uuu) ? double.PositiveInfinity : (uuu - muZ) / sigZ;

                    double[] Jz = new double[13];
                    TruncStdNormalJ(aStarZ, bStarZ, 12, Jz);
                    double Dz = Jz[0];
                    if (!(Dz > pMin) || double.IsNaN(Dz)) whOk = false;
                    else
                    {
                        // Moments of Z we need
                        double EZ3 = EZPow(3, muZ, sigZ, Jz);
                        double EZ6 = EZPow(6, muZ, sigZ, Jz);
                        double EZ9 = EZPow(9, muZ, sigZ, Jz);
                        double EZ12 = EZPow(12, muZ, sigZ, Jz);

                        if (!Tools.IsFinite(EZ3) || !Tools.IsFinite(EZ6) || !Tools.IsFinite(EZ9) || !Tools.IsFinite(EZ12))
                            whOk = false;
                        else
                        {
                            // Mean
                            EX_wh = Xi + Beta * (Alpha * EZ3);

                            // Central moments about μ0 using W = Z^3 - 1
                            double EW2 = (EZ6 - 2.0 * EZ3 + 1.0);
                            double EW3 = (EZ9 - 3.0 * EZ6 + 3.0 * EZ3 - 1.0);
                            double EW4 = (EZ12 - 4.0 * EZ9 + 6.0 * EZ6 - 4.0 * EZ3 + 1.0);

                            double BA = Beta * Alpha;
                            double BA2 = BA * BA, BA3 = BA2 * BA, BA4 = BA3 * BA;
                            m2_wh = BA2 * EW2;
                            m3_wh = BA3 * EW3;
                            m4_wh = BA4 * EW4;
                        }
                    }
                }

                // =====================================================================
                // (C) Normal small-skew path: Y ~ N(α, α) truncated to [L,U]
                // =====================================================================
                bool normalOk = true;
                double EX_n = double.NaN, m2_n = double.NaN, m3_n = double.NaN, m4_n = double.NaN;

                {
                    if (!(Alpha > 0.0) || double.IsNaN(Alpha) || double.IsInfinity(Alpha)) normalOk = false;
                    else
                    {
                        double s = Math.Sqrt(Alpha);
                        double aStar = (L - Alpha) / s;
                        double bStar = double.IsPositiveInfinity(U) ? double.PositiveInfinity : (U - Alpha) / s;

                        double[] J = new double[5];
                        TruncStdNormalJ(aStar, bStar, 4, J);
                        double D = J[0];
                        if (!(D > pMin)) normalOk = false;
                        else
                        {
                            // E[T^m]
                            double[] ET = new double[5];
                            for (int m = 0; m <= 4; m++) ET[m] = J[m] / D;

                            // E[Y^j], j=0..4 with Y = α + s T
                            double[] EYnorm = new double[5];
                            EYnorm[0] = 1.0;
                            for (int j = 1; j <= 4; j++)
                            {
                                double sum = 0.0;
                                for (int m = 0; m <= j; m++)
                                {
                                    double bin = Mathematics.SpecialFunctions.Factorial.BinomialCoefficient(j, m);
                                    sum += bin * Math.Pow(Alpha, j - m) * Math.Pow(s, m) * ET[m];
                                }
                                EYnorm[j] = sum;
                            }

                            EX_n = Xi + Beta * EYnorm[1];

                            // central moments about α for Y
                            double[] CY = new double[5];
                            CY[0] = 1.0; CY[1] = 0.0; // centered at α
                            for (int k = 2; k <= 4; k++)
                            {
                                double sum = 0.0;
                                for (int m = 0; m <= k; m++)
                                {
                                    double bin = Mathematics.SpecialFunctions.Factorial.BinomialCoefficient(k, m);
                                    sum += bin * Math.Pow(-Alpha, k - m) * EYnorm[m];
                                }
                                CY[k] = sum;
                            }

                            double b1 = Beta, b2 = b1 * b1, b3 = b2 * b1, b4 = b3 * b1;
                            m2_n = b2 * CY[2];
                            m3_n = b3 * CY[3];
                            m4_n = b4 * CY[4];

                            if (!(Tools.IsFinite(EX_n) && Tools.IsFinite(m2_n) && Tools.IsFinite(m3_n) && Tools.IsFinite(m4_n)))
                                normalOk = false;
                        }
                    }
                }

                // If no path has usable mass, return divisable=false
                if (!exactOk && !whOk && !normalOk) return (true, false);

                // =====================================================================
                // (D) Smooth weights and blend
                // =====================================================================
                // Unconditional skew magnitude for Gamma(alpha,·): |γ1| = 2 / sqrt(α)
                double absGamma1 = 2.0 / Math.Sqrt(Alpha);

                // WH ↔ Exact around 0.10 (scale 0.03)
                double wWH = Sigmoid(absGamma1, 0.10, 0.03);   // near 1 at very small skew (favor WH)
                if (!whOk) wWH = 0.0;
                if (!exactOk) wWH = 1.0;

                double EX_we = wWH * EX_wh + (1.0 - wWH) * EX_exact;
                double m2_we = wWH * m2_wh + (1.0 - wWH) * m2_exact;
                double m3_we = wWH * m3_wh + (1.0 - wWH) * m3_exact;
                double m4_we = wWH * m4_wh + (1.0 - wWH) * m4_exact;

                // Normal small-skew vs (WH/Exact) around 1e-2 (narrow scale)
                double wN = Sigmoid(absGamma1, 1e-3, 3e-4);    // near 1 below ~1e-3 (favor Normal)
                if (!normalOk) wN = 0.0;

                double EX = wN * EX_n + (1.0 - wN) * EX_we;
                double m2 = wN * m2_n + (1.0 - wN) * m2_we;
                double m3 = wN * m3_n + (1.0 - wN) * m3_we;
                double m4 = wN * m4_n + (1.0 - wN) * m4_we;

                moments = new[] { EX, m2, m3, m4 };
                return (true, true);
            }
            catch
            {
                return (false, true);
            }
        }

    }
}