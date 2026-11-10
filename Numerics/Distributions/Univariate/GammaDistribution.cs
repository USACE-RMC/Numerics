// This file uses functions from the Accord Math Library; the required license header follows:
// Haden Smith (November 2017)
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
// 

using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// Gamma distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b>Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description:</b>
    /// </para>
    /// <para>
    /// The gamma distribution is a two-parameter family of continuous probability
    /// distributions. There are three different parameterizations in common use:</para>
    /// <list type="bullet">
    /// <item><description>
    /// With a shape parameter k and a
    /// scale parameter θ.</description></item>
    /// <item><description>
    /// With a shape parameter α = k and an inverse scale parameter
    /// β = 1/θ, called a rate parameter.</description></item>
    /// <item><description>
    /// With a shape parameter k and a mean
    /// parameter μ = k/β.</description></item>
    /// </list>
    /// <para>
    /// In each of these three forms, both parameters are positive real numbers. The
    /// parameterization with k and θ appears to be more common in econometrics and
    /// certain other applied fields, where e.g. the gamma distribution is frequently
    /// used to model waiting times. For instance, in life testing, the waiting time
    /// until death is a random variable that is frequently modeled with a gamma
    /// distribution. This is the default construction method for this class.
    /// </para>
    /// <para>
    /// This class only uses the parameterization described in the first bullet.
    /// </para>
    /// <para>
    /// If k is an integer, then the distribution represents an Erlang distribution; i.e.,
    /// the sum of k independent exponentially distributed random variables, each of which
    /// has a mean of θ (which is equivalent to a rate parameter of 1/θ).
    /// </para>
    /// <para>
    /// The gamma distribution is the maximum entropy probability distribution for a random
    /// variable X for which E[X] = kθ = α/β is fixed and greater than zero, and <c>E[ln(X)] =
    /// ψ(k) + ln(θ) = ψ(α) − ln(β)</c> is fixed (ψ is the digamma function).
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// This code was developed using two primary sources: 1) Wikipedia; and 2) the Accord Math Library.
    /// <list type="bullet">
    /// <item><description>
    /// Wikipedia, The Free Encyclopedia. Gamma distribution. Available on:
    /// <see href="http://en.wikipedia.org/wiki/Gamma_distribution"/>
    /// </description></item>
    /// <item><description>
    /// Accord Math Library, <see href="http://accord-framework.net"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class GammaDistribution : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
        /*
        * There are three different parameterizations in common use:
        * 1. With a shape parameter κ and a scale parameter θ.
        * 2. With a shape parameter α = κ and an inverse scale parameter β = 1/θ, called a rate parameter.
        * 3. With a shape parameter κ and a mean parameter μ = kθ = α/β.
        * 
        * This class only uses the parameterization described in #1.
        */

        /// <summary>
        /// Constructs a Gamma distribution with scale θ = 10 and shape κ = 2.
        /// </summary>
        public GammaDistribution()
        {
            SetParameters(10d, 2d);
        }

        /// <summary>
        /// Constructs a Gamma distribution with given parameters θ and shape κ.
        /// </summary>
        /// <param name="scale">The scale parameter θ (theta).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GammaDistribution(double scale, double shape)
        {
            SetParameters(scale, shape);
        }

        private double _theta; // scale (θ)
        private double _kappa; // shape (κ)

        /// <summary>
        /// Gets and sets the scale parameter θ (theta).
        /// </summary>
        public double Theta
        {
            get { return _theta; }
            set
            {
                _parametersValid = ValidateParameters(value, Kappa, false) is null;
                _theta = value;
            }
        }

        /// <summary>
        /// Gets the inverse scale parameter of the distribution, β = 1/θ.
        /// </summary>
        public double Rate
        {
            get { return 1.0d / _theta; }
        }

        /// <summary>
        /// Gets and sets the shape parameter κ (kappa).
        /// </summary>
        public double Kappa
        {
            get { return _kappa; }
            set
            {
                _parametersValid = ValidateParameters(Theta, value, false) is null;
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
            get { return UnivariateDistributionType.GammaDistribution; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Gamma"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "G2"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Scale (θ)";
                parmString[1, 0] = "Shape (κ)";
                parmString[0, 1] = Theta.ToString();
                parmString[1, 1] = Kappa.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["θ", "κ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Theta), nameof(Kappa)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Theta, Kappa]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Kappa * Theta; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            // There is no closed form solution
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (Kappa > 1d)
                {
                    return (Kappa - 1d) * Theta;
                }
                else
                {
                    // For shape <= 1 the density is nonincreasing on (0, inf), so the mode is the lower bound.
                    return 0d;
                }
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Theta * Math.Sqrt(Kappa); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 2.0d / Math.Sqrt(Kappa); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + 6.0d / Kappa; }
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
            DistributionNumerics.ValidateSample(sample, 4, true);
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
            var newDistribution = new GammaDistribution(Theta, Kappa);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="scale">The scale parameter θ (theta).</param>
        /// <param name="shape">The shape parameter k.</param>
        public void SetParameters(double scale, double shape)
        {
            // Validate parameters
            _parametersValid = ValidateParameters(scale, shape, false) is null;
            _theta = scale;
            _kappa = shape;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(parameters), "Exactly two parameters are required.");
            SetParameters(parameters[0], parameters[1]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="scale">The scale parameter θ (theta).</param>
        /// <param name="shape">The shape parameter k.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double scale, double shape, bool throwException)
        {
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Theta), "The scale parameter θ (theta) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Theta), "The scale parameter θ (theta) must be positive.");
            }
            if (double.IsNaN(shape) || double.IsInfinity(shape) || shape <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Kappa), "The shape parameter κ (kappa) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "Exactly two parameters are required.");
                if (throwException) throw exception;
                return exception;
            }
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var parms = new double[NumberOfParameters];
            parms[0] = moments[1] * (moments[1] / moments[0]);
            double ratio = moments[0] / moments[1];
            parms[1] = ratio * ratio;
            return parms;
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new GammaDistribution();
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
            double A1 = -0.3080d;
            double A2 = -0.05812d;
            double A3 = 0.01765d;
            double B1 = 0.7213d;
            double B2 = -0.5947d;
            double B3 = -2.1817d;
            double B4 = 1.2113d;
            double L1 = moments[0];
            double L2 = moments[1];
            double CV = L2 / L1;
            double T, theta, kappa;
            if (CV < 0.5)
            {
                T = Math.PI * CV * CV;
                kappa = (1d + A1 * T) / (T * (1d + T * (A2 + T * A3)));
            }
            else
            {
                T = 1d - CV;
                kappa = T * (B1 + T * B2) / (1d + T * (B3 + T * B4));
            }
            theta = L1 / kappa;     
            return [theta, kappa];
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The L-scale λ2 = θ Γ(κ + 1/2) / (√π Γ(κ)) is evaluated from the ratio Γ(κ) / Γ(κ + 1/2) of
        /// <see cref="Gamma.RatioHalf(double)"/>, which is formed from ln B(κ, 1/2) for κ ≥ 1 and does not overflow,
        /// so λ2 stays finite for large shapes (small skews), where the separate gamma functions overflow beyond
        /// κ ≈ 171. Against mpmath the relative error of λ2 is below 1E-15 for shapes from 0.01 to 1E6. τ3 and τ4 are
        /// the rational-function approximations of Hosking (routine LMRGAM), accurate to about 1E-6.
        /// </remarks>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            double alpha = parameters[1];
            double beta = parameters[0];
            double L1 = alpha * beta;
            double L2 = Math.Abs(beta / (Math.Sqrt(Math.PI) * Gamma.RatioHalf(alpha)));
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
        /// <remarks>Requires at least four finite, nonconstant observations. Preserves the
        /// legacy moment initialization and family-specific bounds whenever they are finite, ordered, and
        /// contain the initial values, including usable samples containing nonpositive observations.
        /// Otherwise, the exceptional-input fallback requires positive observations and evaluates the moment
        /// initialization in unit coordinates to preserve small and large scales.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample or a finite feasible initialization is invalid.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            return DistributionNumerics.PreferLegacyConstraints(
                () => GetLegacyParameterConstraints(sample), () => GetRobustParameterConstraints(sample));
        }

        /// <summary>Preserves the established initialization and family-specific prior envelope.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>The legacy initial values and lower and upper bounds.</returns>
        private Tuple<double[], double[], double[]> GetLegacyParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            initialVals = ParametersFromMoments(Statistics.ProductMoments(sample));
            // Get bounds of scale
            lowerVals[0] = Tools.DoubleMachineEpsilon;
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d));
            // Get bounds of shape
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Retains the established constraint initializer arithmetic for ordinary samples.</summary>
        /// <param name="moments">Sample moments.</param>
        /// <returns>The legacy initial parameter values.</returns>
        private double[] LegacyConstraintParametersFromMoments(IList<double> moments)
        {
            var parms = new double[NumberOfParameters];
            parms[0] = 1d / (moments[0] / Math.Pow(moments[1], 2d));
            parms[1] = Math.Pow(moments[0], 2d) / Math.Pow(moments[1], 2d);
            return parms;
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4, true);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            double normalization = DistributionNumerics.InitializationScale(sample);
            var normalized = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) normalized[i] = sample[i] / normalization;
            var initialVals = ParametersFromMoments(Statistics.ProductMoments(normalized));
            initialVals[0] *= normalization;
            DistributionNumerics.PositiveParameterBounds(initialVals[0], out lowerVals[0], out upperVals[0]);
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
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
                var GM = new GammaDistribution();
                GM.SetParameters(x);
                return GM.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <summary>
        /// Estimates parameters using a Newton-Raphson method.
        /// </summary>
        /// <param name="sample">Array of sample data.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is insufficient, constant, or contains a nonpositive or nonfinite observation.</exception>
        /// <exception cref="ArgumentException">The transformed observation statistic is not-a-number.</exception>
        public void MLE_NR(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2, true);
            double lnsum = 0d;
            for (int i = 0; i < sample.Count; i++)
                lnsum += Math.Log(sample[i]);
            double mean = Statistics.Mean(sample);
            double s = Math.Log(mean) - lnsum / sample.Count;
            if (double.IsNaN(s))
            {
                throw new ArgumentException("Observation vector contains negative values.", "observations");
            }
            double K_hat = (3d - s + Math.Sqrt((s - 3d) * (s - 3d) + 24d * s)) / (12d * s);
            double newK = NewtonRaphson.Solve((k) => Math.Log(k) - Gamma.Digamma(k) - s, (k) => 1d / k - Gamma.Trigamma(k), K_hat, 0.000001d, 100, true);
            double newTheta = mean / newK;
            SetParameters(newTheta, newK);
        }

        /// <summary>
        /// Estimates parameters using an approximation proposed by Bobee.
        /// </summary>
        /// <param name="sample">Array of sample data.</param>
        public void MLE_Bobee(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 2, true);
            double A = Statistics.Mean(sample);
            double G = Statistics.GeometricMean(sample);
            double U = Math.Log(A) - Math.Log(G);
            double La = (1d + Math.Sqrt(1d + 4d * U / 3d)) / (4d * U);
            double dL = 0.04475d * Math.Pow(0.26d, La);
            double K = La - dL;
            double T = A / K;
            SetParameters(T, K);
        }

        /// <inheritdoc/>
        public override double PDF(double X)
        {
            return Math.Exp(LogPDF(X));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated in log space, so far-tail densities that underflow <see cref="PDF(double)"/>
        /// keep a finite log density.
        /// When <c>X = 0</c> and the shape <c>κ &lt; 1</c>, the Gamma density has a genuine
        /// integrable singularity and this method intentionally returns positive infinity.
        /// </remarks>
        public override double LogPDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Theta, Kappa, true);
            if (X < Minimum || double.IsPositiveInfinity(X)) return double.NegativeInfinity;
            if (X == 0) return Kappa == 1 ? -Math.Log(Theta) : Kappa < 1 ? double.PositiveInfinity : double.NegativeInfinity;
            double unit = X / Theta;
            double lf = unit == 0 && X > 0
                ? (Kappa - 1) * (Math.Log(X) - Math.Log(Theta)) - Gamma.LogGamma(Kappa) - Math.Log(Theta)
                : DistributionNumerics.GammaLogDensity(Kappa, unit) - Math.Log(Theta);
            return double.IsNaN(lf) ? double.NegativeInfinity : lf;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The regularized lower incomplete gamma function P(κ, x/θ) of <see cref="Gamma.LowerIncomplete(double, double)"/>;
        /// when x/θ underflows to zero for a positive x, the value is the exponential of <see cref="LogCDF(double)"/>.
        /// </remarks>
        public override double CDF(double X)
        {
            if (!_parametersValid) ValidateParameters(Theta, Kappa, true);
            if (X <= 0) return 0;
            double unit = X / Theta;
            return unit == 0 ? Math.Exp(LogCDF(X)) : Gamma.LowerIncomplete(Kappa, unit);
        }

        /// <inheritdoc/>
        public override double LogCDF(double X)
        {
            if (!_parametersValid) ValidateParameters(Theta, Kappa, true);
            if (X <= 0) return double.NegativeInfinity;
            double unit = X / Theta;
            return unit == 0 ? Kappa * (Math.Log(X) - Math.Log(Theta)) - Gamma.LogGammaOnePlus(Kappa)
                : DistributionNumerics.GammaLogCDF(Kappa, unit);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The regularized upper incomplete gamma function Q(κ, x/θ) of <see cref="Gamma.UpperIncomplete(double, double)"/>;
        /// when x/θ underflows to zero for a positive x, the value is the exponential of <see cref="LogCCDF(double)"/>.
        /// </remarks>
        public override double CCDF(double X)
        {
            if (!_parametersValid) ValidateParameters(Theta, Kappa, true);
            if (X <= 0) return 1;
            double unit = X / Theta;
            return unit == 0 ? Math.Exp(LogCCDF(X)) : Gamma.UpperIncomplete(Kappa, unit);
        }

        /// <inheritdoc/>
        public override double LogCCDF(double X)
        {
            if (!_parametersValid) ValidateParameters(Theta, Kappa, true);
            if (X <= 0) return 0;
            double unit = X / Theta;
            return unit == 0 ? Tools.Log1mExp(LogCDF(X)) : DistributionNumerics.GammaLogSurvival(Kappa, unit);
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Theta, Kappa, true);
            // Preserve the exact exponential identity, including subnormal probabilities whose
            // unnecessary log/exp round trip loses range on .NET Framework.
            if (Kappa == 1) return -Theta * Tools.Log1p(-probability);
            double unitQuantile = Gamma.InverseIncompleteRatio(Kappa, probability, false);
            // The lower-tail series has negligible higher-order terms at these unit quantiles.
            // Combine logs before multiplication when the scale can rescue an underflowed unit value.
            return unitQuantile < 1E-200
                ? Math.Exp(Math.Log(Theta) + (Math.Log(probability) + Gamma.LogGammaOnePlus(Kappa)) / Kappa)
                : unitQuantile * Theta;
        }

        /// <summary>
        /// Returns the inverse CDF using the modified Wilson-Hilferty transformation.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>The approximate gamma quantile in physical coordinates.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The probability is outside the closed unit interval or the distribution parameters are invalid.</exception>
        /// <remarks>
        /// Cornish-Fisher transformation (Fisher and Cornish, 1960) for abs(skew) less than or equal to 2. If abs(skew) > 2 then use Modified Wilson-Hilferty transformation (Kirby,1972).
        /// </remarks>
        public double WilsonHilfertyInverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Theta, Kappa, true);
            return Theta * (FrequencyFactorKp(Skewness, probability) * Math.Sqrt(Kappa) + Kappa);
        }

        /// <summary>
        /// Gets the K frequency factor given the skewness coefficient through Cornish-Fisher transformation (Fisher and Cornish, 1960) for abs(skew) less than or equal to 2.
        /// If abs(skew) > 2 then use Modified Wilson-Hilferty transformation (Kirby,1972).
        /// </summary>
        /// <param name="skewness">Coefficient of skewness.</param>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>The named approximate frequency factor, not the actual Gamma quantile.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="skewness"/> is nonfinite or <paramref name="probability"/> is outside the closed unit interval or not-a-number.</exception>
        /// <remarks>
        /// For abs(skew) at most 2 the frequency factor is the Cornish-Fisher polynomial in the skew through its sixth power,
        /// <c>U + Σₖ Cᵏ Pₖ(U)</c> for k = 1 to 6, with <c>U</c> the standard Normal quantile. It is analytic at zero skew,
        /// where it equals <c>U</c>, and its first four coefficients and the sixth are those of the exact skew series of the
        /// standardized gamma quantile, <c>(U² - 1)/6</c>, <c>U(U² - 7)/144</c>, <c>-(3U⁴ + 7U² - 16)/6480</c>,
        /// <c>U(9U⁴ + 256U² - 433)/622080</c> and <c>-U(3753U⁶ + 4353U⁴ - 289517U² - 289717)/9405849600</c>.
        /// The fifth-order coefficient is <c>(12U⁶ - 143U⁴ - 923U² + 1472)/6531840</c>. The exact series term is
        /// <c>(12U⁶ - 243U⁴ - 923U² + 1472)/6531840</c>, but the established coefficient is kept deliberately: against the
        /// actual gamma quantile at 11 probabilities from 1E-6 to 1 - 1E-6, the largest error of the polynomial is 1.4E-4
        /// (versus 1.6E-4 with the exact term) for abs(skew) up to 0.5, 1.2E-2 (versus 2.0E-2) up to 1 and 0.68 (versus 0.93)
        /// up to 2, because the established term partly offsets the truncation after the sixth power. Its published source
        /// has not been verified, so it is treated as an empirically tuned coefficient.
        /// For abs(skew) below 1E-4 the Normal variate is returned at probabilities zero and one as the endpoint limit.
        /// Large skew magnitude is capped at 9.75 before evaluating its powers. Negative
        /// skew uses the reflected probability and sign of the corresponding positive-skew approximation.
        /// The Cornish-Fisher coefficients and the Wilson-Hilferty constants for A, B and G are attributed to Hoshi and
        /// Burges (1981). They are documented as implemented and have not been checked against Kirby (1972), Kite (1977),
        /// Hoshi and Burges (1981), PeakFQ, HEC-SSP or Bulletin 17B/17C reference code.
        /// The frequency factor is an approximation by design and does not return the support limits of the gamma
        /// distribution at the probability endpoints. At probability one it is positive infinity for zero skew and for
        /// abs(skew) below 1E-4, not-a-number for abs(skew) from 1E-4 through 2, where the polynomial is evaluated at an
        /// infinite Normal variate, and for abs(skew) above 2 the Wilson-Hilferty value: positive infinity for positive skew
        /// and the clamped bound <c>-2/C</c> (with C capped at 9.75 in magnitude) for negative skew. At probability zero the
        /// Normal variate is finite (about -38.4674), so the polynomial and the Wilson-Hilferty form are evaluated there;
        /// for skew above 2 the Wilson-Hilferty clamp gives <c>-2/C</c>.
        /// </remarks>
        public static double FrequencyFactorKp(double skewness, double probability)
        {
            if (!Tools.IsFinite(skewness)) throw new ArgumentOutOfRangeException(nameof(skewness));
            if (!(probability >= 0 && probability <= 1)) throw new ArgumentOutOfRangeException(nameof(probability));
            double C = skewness;
            double absC = Math.Abs(C);
            // The series equals the standard Normal variate at zero skew; returning it directly keeps the endpoint limits.
            // Below absolute skew 1E-4 the variate is also the limit at probabilities zero and one.
            if (C == 0d || (absC < 1E-4 && (probability == 0d || probability == 1d))) return Normal.StandardZ(probability);

            // If abs(skew) is less than or equal to 2, use Cornish-Fisher transformation (Fisher and Cornish, 1960)
            if (absC <= 2d) return CornishFisherFrequencyFactor(C, Normal.StandardZ(probability), out _);
            else
            {
                // If abs(skew) is greater than 2, use Modified Wilson-Hilferty transformation (Kirby, 1972)
                // Only, valid if abs(skew) <= 9.75. Enforce limits.
                double sign = Math.Sign(C);
                C = Math.Min(absC, 9.75);
                absC = Math.Abs(C);

                // Hoshi and Burges (1981b) gave polynomial expressions for 1/A, B, G and H^3 as a function of Cs
                // Compute skew orders
                double C2 = Math.Pow(absC, 2d);
                double C3 = Math.Pow(absC, 3d);
                double C4 = Math.Pow(absC, 4d);
                double C5 = Math.Pow(absC, 5d);
                // Compute A
                double a0 = 0.00199447d;
                double a1 = 0.48489d;
                double a2 = 0.0230935d;
                double a3 = -0.0152435d;
                double a4 = 0.00160597d;
                double a5 = -0.000055869d;
                double A = 1d / (a0 + a1 * C + a2 * C2 + a3 * C3 + a4 * C4 + a5 * C5);
                // Compute B
                double b0 = 0.990562d;
                double b1 = 0.0319647d;
                double b2 = -0.0274231d;
                double b3 = 0.00777405d;
                double b4 = -0.000571184d;
                double b5 = 0.0000142077d;
                double B = b0 + b1 * C + b2 * C2 + b3 * C3 + b4 * C4 + b5 * C5;
                // Compute G
                double g0 = -0.00385205d;
                double g1 = 1.00426d;
                double g2 = 0.00651207d;
                double g3 = -0.0149166d;
                double g4 = 0.00163945d;
                double g5 = -0.0000583804d;
                double G = g0 + g1 * C + g2 * C2 + g3 * C3 + g4 * C4 + g5 * C5;
                // Compute H
                double H = Math.Pow(B - 2.0d / absC / A, 1d / 3d);
                // Reflect the normal score directly; forming 1-p would round a tiny p to one.
                return sign * A * (Math.Pow(Math.Max(H, 1.0d - Math.Pow(G / 6.0d, 2d) + G / 6.0d * sign * Normal.StandardZ(probability)), 3d) - B);
            }
        }

        /// <summary>
        /// Gets the partial derivative of the frequency factor Kp with respect to skew.
        /// </summary>
        /// <param name="skewness">Coefficient of skewness.</param>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <returns>The partial derivative of the frequency factor with respect to skewness.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="skewness"/> is nonfinite or <paramref name="probability"/> is not finite and strictly between zero and one.</exception>
        /// <remarks>
        /// For abs(skew) at most 2 this is the derivative of the Cornish-Fisher polynomial of
        /// <see cref="FrequencyFactorKp(double, double)"/>, analytic through zero skew. It uses the same deliberately tuned
        /// fifth-order coefficient <c>-143U⁴</c> in place of the exact series term <c>-243U⁴</c>, which measured more accurate
        /// against the actual gamma quantile (largest errors 1.4E-4, 1.2E-2 and 0.68 versus 1.6E-4, 2.0E-2 and 0.93 for abs(skew)
        /// up to 0.5, 1 and 2). Above 2 it is a numerical derivative of the Wilson-Hilferty form. The probability must be
        /// strictly between zero and one. The coefficients share the provenance described for
        /// <see cref="FrequencyFactorKp(double, double)"/>: they are documented as implemented, without a verified published
        /// source.
        /// </remarks>
        public static double PartialKp(double skewness, double probability)
        {
            if (!Tools.IsFinite(skewness)) throw new ArgumentOutOfRangeException(nameof(skewness));
            DistributionNumerics.ValidateProbability(probability);
            double C = skewness;
            double absC = Math.Abs(C);
            // If abs(skew) is less than or equal to 2, use Cornish-Fisher transformation (Fisher and Cornish, 1960)
            if (absC <= 2d)
            {
                CornishFisherFrequencyFactor(C, Normal.StandardZ(probability), out double derivative);
                return derivative;
            }
            else
            {
                // If abs(skew) is greater than 2, use Modified Wilson-Hilferty transformation (Kirby, 1972)
                return NumericalDerivative.Derivative((x) => FrequencyFactorKp(x, probability), skewness, 0.0001d);
            }
        }

        /// <summary>Evaluates the Cornish-Fisher frequency factor polynomial and its skew derivative.</summary>
        /// <param name="skew">The coefficient of skewness, at most 2 in absolute value.</param>
        /// <param name="u">The standard Normal quantile.</param>
        /// <param name="derivative">The derivative of the frequency factor with respect to the skew.</param>
        /// <returns>The frequency factor <c>u + Σₖ skewᵏ Pₖ(u)</c> for k = 1 to 6.</returns>
        /// <remarks>The polynomials in <c>u</c> are those of the Cornish-Fisher transformation (Fisher and Cornish,
        /// 1960) with their powers of two folded into the denominators, and both sums are evaluated in Horner form.</remarks>
        private static double CornishFisherFrequencyFactor(double skew, double u, out double derivative)
        {
            double u2 = u * u, u4 = u2 * u2, u6 = u4 * u2;
            double p1 = (u2 - 1d) / 6d;
            double p2 = u * (u2 - 7d) / 144d;
            double p3 = -(6d * u4 + 14d * u2 - 32d) / 12960d;
            double p4 = u * (9d * u4 + 256d * u2 - 433d) / 622080d;
            double p5 = (12d * u6 - 143d * u4 - 923d * u2 + 1472d) / 6531840d;
            double p6 = -u * (3753d * u6 + 4353d * u4 - 289517d * u2 - 289717d) / 9405849600d;
            derivative = p1 + skew * (2d * p2 + skew * (3d * p3 + skew * (4d * p4 + skew * (5d * p5 + skew * 6d * p6))));
            return u + skew * (p1 + skew * (p2 + skew * (p3 + skew * (p4 + skew * (p5 + skew * p6)))));
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override UnivariateDistributionBase Clone()
        {
            return new GammaDistribution(Theta, Kappa);
        }

        /// <inheritdoc/>
        /// <remarks>Uses scale-theta and shape-kappa coordinates. MoM is the sample-mean/sample-variance
        /// sandwich; MLE is the expected Fisher inverse. Both require positive sample size and valid parameters.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters or sample size are invalid.</exception>
        /// <exception cref="NotImplementedException">The estimation method is unsupported.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments &&
                estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_theta, _kappa, true);
            // Compute covariance in user-facing (θ=scale, κ=shape) parameterization.
            double t = Theta, k = Kappa;
            var covar = new double[2, 2];
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                // Algebraic simplification of the same sample-mean/sample-variance sandwich.
                double scaled = t / Math.Sqrt(sampleSize);
                covar[0, 0] = (scaled * scaled) * (2 + 3 / k);
                covar[1, 1] = 2 * (k / sampleSize) * (k + 1);
                covar[0, 1] = -2 * (t / sampleSize) * (k + 1);
                covar[1, 0] = covar[0, 1];
            }
            else
            {
                // MLE: Fisher information inverse in (θ, κ) space.
                // Transformed from (α=1/θ, κ) via delta method.
                double logResidual = Math.Log(DistributionNumerics.GammaScaledFisherResidual(k));
                double logScale = Math.Log(t), logShape = Math.Log(k), logCount = Math.Log(sampleSize);
                covar[0, 0] = Math.Exp(2 * logScale - logCount + DistributionNumerics.LogSum(-logResidual, -logShape)); // Var(θ̂)
                covar[1, 1] = Math.Exp(2 * logShape - logCount - logResidual); // Var(κ̂)
                covar[0, 1] = -Math.Exp(logScale + logShape - logCount - logResidual); // Cov(θ̂, κ̂), negative
                covar[1, 0] = covar[0, 1];
            }
            return covar;
        }

        /// <summary>
        /// The quantile variance given probability and sample size.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="estimationMethod">The distribution parameter estimation method.</param>
        /// <returns>The covariance quadratic form for the actual inverse CDF.</returns>
        /// <remarks>Both supported methods differentiate the actual Gamma quantile. The named
        /// Wilson-Hilferty and frequency-factor approximations are not used here. An algebraically
        /// equivalent sum of two nonnegative terms preserves the mean-direction variance at large
        /// shape; physical scale is restored in logarithms. Probability must be finite and strictly interior.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The probability, sample size, or distribution parameters are invalid.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not maximum likelihood or method of moments.</exception>
        /// <exception cref="InvalidOperationException">A gamma tail or quantile calculation does not converge.</exception>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Theta, Kappa, true);
            DistributionNumerics.ValidateSampleSize(sampleSize);
            double logCoefficient;
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
                logCoefficient = -Math.Log(DistributionNumerics.GammaScaledFisherResidual(Kappa));
            else if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
                logCoefficient = Math.Log(2) + DistributionNumerics.LogSum(0, -Math.Log(Kappa));
            else throw new NotImplementedException();

            double unitQuantile = Gamma.InverseIncompleteRatio(Kappa, probability, false);
            double logQuantile, logDifference;
            if (unitQuantile < 1E-200)
            {
                // The same lower-tail log quantile and derivative used by the public gradient,
                // with kappa cancelled before forming q-kappa*dq/dkappa.
                logQuantile = (Math.Log(probability) + Gamma.LogGammaOnePlus(Kappa)) / Kappa;
                if (double.IsNegativeInfinity(logQuantile)) return 0;
                logDifference = logQuantile + Math.Log(Math.Abs(1 - Gamma.Digamma(1 + Kappa) + logQuantile));
            }
            else
            {
                logQuantile = Math.Log(unitQuantile);
                double derivative = DistributionNumerics.GammaQuantileShapeDerivative(Kappa, unitQuantile);
                logDifference = Math.Log(Math.Abs(unitQuantile - Kappa * derivative));
            }
            // Exactly q^2/kappa + A*(q-kappa*dq/dkappa)^2, before theta^2/n.
            double logVariance = DistributionNumerics.LogSum(2 * logQuantile - Math.Log(Kappa), logCoefficient + 2 * logDifference);
            return Math.Exp(2 * Math.Log(Theta) - Math.Log(sampleSize) + logVariance);
        }

        /// <inheritdoc/>
        /// <remarks>Returns [unit-scale quantile, theta times its implicit shape derivative]. The
        /// convergent incomplete-Gamma equations are differentiated at the actual quantile.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(_theta, _kappa, true);
            double unitQuantile = Gamma.InverseIncompleteRatio(Kappa, probability, false);
            double shapeDerivative;
            if (unitQuantile < 1E-200)
            {
                double logQuantile = (Math.Log(probability) + Gamma.LogGammaOnePlus(Kappa)) / Kappa;
                double logarithmicDerivative = (Gamma.Digamma(1 + Kappa) - logQuantile) / Kappa;
                shapeDerivative = Math.Exp(Math.Log(Theta) + logQuantile + Math.Log(logarithmicDerivative));
            }
            else shapeDerivative = Theta * DistributionNumerics.GammaQuantileShapeDerivative(Kappa, unitQuantile);
            return [unitQuantile, shapeDerivative];
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

        /// <inheritdoc/>
        public override double[] ConditionalMoments(double a, double b)
        {
            if (a >= b)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            double alpha = Kappa; // shape > 0
            double beta = Theta;  // scale > 0
            if (!(alpha > 0.0) || !(beta > 0.0))
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Effective bounds within support [0, ∞)
            double A = Math.Max(0.0, a);
            double B = b;

            // Standardized bounds in z = x / beta
            double zA = A / beta;
            double zB = double.IsPositiveInfinity(B) ? double.PositiveInfinity : (B / beta);

            // Regularized lower incomplete gamma CDF P(s, z)
            static double Pg(double s, double z)
            {
                if (double.IsNegativeInfinity(z) || z <= 0.0) return 0.0;
                if (double.IsPositiveInfinity(z)) return 1.0;
                return Mathematics.SpecialFunctions.Gamma.LowerIncomplete(s, z); // P(s,z)
            }

            // Normalizing probability over (a,b): P0 = P(alpha, zB) - P(alpha, zA)
            double P0 = Pg(alpha, zB) - Pg(alpha, zA);
            if (P0 <= 1e-15)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // For stability: Γ(α+k)/Γ(α) via log-gamma differences
            double lgAlpha = Mathematics.SpecialFunctions.Gamma.LogGamma(alpha);
            double[] gr = new double[5]; // gr[k] = Γ(α+k)/Γ(α)
            for (int k = 0; k <= 4; k++)
                gr[k] = Math.Exp(Mathematics.SpecialFunctions.Gamma.LogGamma(alpha + k) - lgAlpha);

            // Truncated raw moments E[X^k | a < X < b], k = 1..4
            double[] EX = new double[5];
            EX[0] = 1.0; // by definition under conditioning
            for (int k = 1; k <= 4; k++)
            {
                double deltaPk = Pg(alpha + k, zB) - Pg(alpha + k, zA);
                EX[k] = Math.Pow(beta, k) * gr[k] * (deltaPk / P0);
            }

            // Central moments about the *unconditional* mean μ = αβ
            double mu = alpha * beta;
            double m1 = EX[1];
            double m2 = EX[2] - 2.0 * mu * EX[1] + mu * mu;
            double m3 = EX[3] - 3.0 * mu * EX[2] + 3.0 * mu * mu * EX[1] - mu * mu * mu;
            double m4 = EX[4]
                      - 4.0 * mu * EX[3]
                      + 6.0 * mu * mu * EX[2]
                      - 4.0 * mu * mu * mu * EX[1]
                      + mu * mu * mu * mu;

            return new[] { m1, m2, m3, m4 };
        }

    }
}
