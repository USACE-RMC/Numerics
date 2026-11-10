// This file uses functions from the Accord Math Library; the required license header follows:
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
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
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
    public class StudentT : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IBootstrappable
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
        private double _degreesOfFreedom;
        private double _densityConstant;
        private double _sqrtNuPi;
        private double _logDensityConstant;

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
        public double DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
            set
            {
                _parametersValid = ValidateParameters(Mu, Sigma, value, false) is null;
                _degreesOfFreedom = value;
                UpdateDensityConstants();
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
        /// <param name="v">The degrees of freedom ν (nu). Range: ν ≥ 1.</param>
        public void SetParameters(double mu, double sigma, double v)
        {
            _mu = mu;
            _sigma = sigma;
            _degreesOfFreedom = v;
            UpdateDensityConstants();
            _parametersValid = ValidateParameters(mu, sigma, v, false) is null;
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
        public ArgumentOutOfRangeException? ValidateParameters(double location, double scale, double v, bool throwException)
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
            if (double.IsNaN(v) || double.IsInfinity(v) || v < 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
                return new ArgumentOutOfRangeException(nameof(DegreesOfFreedom), "The degrees of freedom ν (nu) must greater than or equal to one.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Only <see cref="ParameterEstimationMethod.MaximumLikelihood"/> is supported, with the location, scale and
        /// degrees of freedom all estimated. The method of moments is not provided because its estimate of ν exists
        /// only for samples with positive excess kurtosis and assumes ν &gt; 4.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than four values, a nonfinite value, or no variation.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not maximum likelihood.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSample(sample, 4);
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
            var newDistribution = new StudentT(Mu, Sigma, DegreesOfFreedom);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The location starts at the sample median. The degrees of freedom start at the moment estimate
        /// ν = 4 + 6/κ from the sample excess kurtosis κ, limited to [2.5, 100], or at 30 when κ is not positive.
        /// The scale starts at s √((ν − 2)/ν), where s is the sample standard deviation.
        /// </para>
        /// <para>
        /// The location is bounded by ±10^(⌈log10(max(|median|, s))⌉ + 1), the scale by
        /// [min(ε, s/10), 10^(⌈log10(s)⌉ + 1)], and the degrees of freedom by [1, 1000]. Above ν = 1000 the
        /// likelihood is practically indistinguishable from that of the Normal distribution.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample has fewer than four values, a nonfinite value, or no variation.</exception>
        /// <exception cref="ArgumentException">The sample moments are not finite.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            var moments = Statistics.ProductMoments(sample);
            double location = Statistics.Percentile(sample, 0.5d);
            double standardDeviation = moments[1];
            if (!(standardDeviation > 0d) || double.IsInfinity(standardDeviation) || !Tools.IsFinite(location))
                throw new ArgumentException("Sample moments must be finite with positive dispersion.", nameof(sample));
            double v = moments[3] > 0d ? Tools.Clamp(4d + 6d / moments[3], 2.5d, 100d) : 30d;
            double scale = standardDeviation * Math.Sqrt((v - 2d) / v);
            double locationBound = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Max(Math.Abs(location), standardDeviation)) + 1d));
            double scaleBound = Math.Pow(10d, Math.Ceiling(Math.Log10(standardDeviation) + 1d));
            var initialVals = new double[] { location, scale, v };
            var lowerVals = new double[] { -locationBound, Math.Min(Tools.DoubleMachineEpsilon, scale / 10d), 1d };
            var upperVals = new double[] { locationBound, scaleBound, 1000d };
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The log-likelihood is maximized over the location, scale and degrees of freedom with the Nelder-Mead
        /// method, within the bounds from <see cref="GetParameterConstraints(IList{double})"/>.
        /// </remarks>
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
                var T = new StudentT();
                T.SetParameters(x);
                return T.LogLikelihood(sample);
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
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double inverseSigma = 1.0d / Sigma;
            double Z = (x - Mu) * inverseSigma;
            if (DegreesOfFreedom >= 100000000.0d)
            {
                return Normal.StandardPDF(Z) * inverseSigma;
            }
            return _densityConstant * Math.Pow(1.0d + Z * Z / DegreesOfFreedom, -0.5d * (DegreesOfFreedom + 1.0d)) * inverseSigma / _sqrtNuPi;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated in log space as ln Γ((ν + 1)/2) − ln Γ(ν/2) − ln √(νπ) − ln σ − ((ν + 1)/2) ln(1 + z²/ν),
        /// with z = (x − µ)/σ, so far-tail densities that underflow <see cref="PDF(double)"/> keep a finite log
        /// density. For ν ≥ 10^8 the Normal log density is used, as in <see cref="PDF(double)"/>.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double Z = (x - Mu) / Sigma;
            if (DegreesOfFreedom >= 100000000.0d)
            {
                return Normal.StandardLogPDF(Z) - Math.Log(Sigma);
            }
            // ln(1 + z²/ν), written as 2 ln|z| − ln ν when z²/ν overflows.
            double ratio = Z * Z / DegreesOfFreedom;
            double logKernel = double.IsInfinity(ratio) ? 2.0d * Math.Log(Math.Abs(Z)) - Math.Log(DegreesOfFreedom) : Tools.Log1p(ratio);
            return _logDensityConstant - 0.5d * (DegreesOfFreedom + 1.0d) * logKernel - Math.Log(Sigma);
        }

        /// <summary>
        /// Caches the density factors that depend only on ν: Γ((ν + 1)/2) / Γ(ν/2), evaluated as the
        /// exponential of ln Γ(1/2) - ln B(ν/2, 1/2) so that no large ln Γ values are subtracted, √(νπ), and the
        /// log of their ratio.
        /// </summary>
        private void UpdateDensityConstants()
        {
            double logGammaRatio = Gamma.LogGamma(0.5d) - Beta.LogFunction(_degreesOfFreedom / 2.0d, 0.5d);
            _densityConstant = Math.Exp(logGammaRatio);
            _sqrtNuPi = Math.Sqrt(_degreesOfFreedom * Math.PI);
            _logDensityConstant = logGammaRatio - Math.Log(_sqrtNuPi);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// With z = (x − µ)/σ &lt; 0, the CDF is the regularized incomplete beta function I_w(ν/2, ν/2) with
        /// w = (z + √(z² + ν)) / (2√(z² + ν)), evaluated in the algebraically equal form
        /// w = 1 / (2√(1 + s²)(s + √(1 + s²))), s = |z|/√ν, which avoids the cancellation in z + √(z² + ν), so the
        /// lower tail keeps full relative precision. For z² &lt; ½ the equal form ½ − ½ I_y(1/2, ν/2), y = z²/(ν + z²),
        /// is used, which is more accurate at large ν. For z &gt; 0 the CDF is 1 − F(−z) by symmetry, and F(0) = ½.
        /// Beyond |z| = 10^150 √ν the probability is the exponential of <see cref="LogCDF(double)"/>.
        /// </para>
        /// <para>
        /// For ν ≥ 30 and ν/(ν + z²) ≥ 0.7 the lower tail is instead half of I_x(ν/2, 1/2), x = ν/(ν + z²), from the
        /// asymptotic expansion in the large shape ν/2 (routine BGRAT of algorithm 708, DiDonato and Morris, 1992;
        /// <see cref="Beta.AsymptoticIncomplete"/>), whose leading term is erfc(√((ν/2 − 1/4) ln(1 + z²/ν))), with
        /// z²/ν carried in double-double arithmetic. This covers every ν, so no Normal approximation is used. On the
        /// 97-point mpmath grid of the unit tests (ν from 1 to 10^10, |z| from 10^-6 to 10^6), the relative error of the
        /// two-sided tail 2 CCDF(|z|) is at most 3.4E-15 where it is at least 0.01, 1.1E-14 between 10^-10 and 0.01, and
        /// 8E-13 below 10^-10, where the tail is ill-conditioned in z; the unit tests also check ν up to 10^30 to 1E-12.
        /// </para>
        /// <para>
        /// CDF(+∞) = 1, CDF(−∞) = 0 and CDF(NaN) = NaN.
        /// </para>
        /// </remarks>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            //
            double Z = (x - Mu) / Sigma;
            if (double.IsNaN(Z)) return double.NaN;
            return StandardCdf(Z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as the lower tail at −z, z = (x − µ)/σ, by symmetry, so the upper tail keeps full relative
        /// precision instead of being formed as 1 − CDF.
        /// </remarks>
        public override double CCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double Z = (x - Mu) / Sigma;
            if (double.IsNaN(Z)) return double.NaN;
            return StandardCdf(-Z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// When the lower-tail probability is below about 10^-290, the logarithm is evaluated directly as
        /// ln(I_y(ν/2, 1/2)/2) with y = ν/(ν + z²), from the continued fractions of the incomplete beta function,
        /// so it stays finite after the probability itself underflows.
        /// </remarks>
        public override double LogCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double Z = (x - Mu) / Sigma;
            if (double.IsNaN(Z)) return double.NaN;
            return StandardLogCdf(Z);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as the log lower tail at −z by symmetry; see <see cref="LogCDF(double)"/>.
        /// </remarks>
        public override double LogCCDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            double Z = (x - Mu) / Sigma;
            if (double.IsNaN(Z)) return double.NaN;
            return StandardLogCdf(-Z);
        }

        /// <summary>
        /// Returns the CDF of the standard Student's t distribution with ν degrees of freedom.
        /// </summary>
        /// <param name="z">The standardized variate, not NaN.</param>
        /// <returns>The lower-tail probability at <paramref name="z"/>.</returns>
        private double StandardCdf(double z)
        {
            if (z == 0d) return 0.5d;
            if (z > 0d) return 1.0d - StandardCdf(-z);
            if (double.IsNegativeInfinity(z)) return 0.0d;
            double v = DegreesOfFreedom;
            double tail = AsymptoticTwoSidedTail(-z, false);
            if (!double.IsNaN(tail)) return 0.5d * tail;
            double s = -z / Math.Sqrt(v);
            // Near the center, F = ½ − ½ I_y(1/2, ν/2) with y = z²/(ν + z²) stays accurate for large ν, and F ≥ 0.3,
            // so the subtraction does not cancel.
            if (z * z < 0.5d) return 0.5d - 0.5d * Beta.Incomplete(0.5d, v / 2.0d, s * s / (1.0d + s * s));
            // s² would leave the floating-point range.
            if (s > 1E150) return Math.Exp(StandardLogLowerTail(z));
            double q = Math.Sqrt(1.0d + s * s);
            return Beta.Incomplete(v / 2.0d, v / 2.0d, 1.0d / (2d * q * (s + q)));
        }

        /// <summary>
        /// Returns the two-sided tail P(|T| ≥ s) = I_x(ν/2, 1/2), x = ν/(ν + s²), or its logarithm, by the asymptotic
        /// expansion in the large shape ν/2, where it applies.
        /// </summary>
        /// <param name="s">The magnitude |z|, positive.</param>
        /// <param name="logarithm">True to return the natural logarithm, which stays finite after the tail underflows.</param>
        /// <returns>The tail or its logarithm when ν ≥ 30 and x ≥ 0.7; NaN otherwise.</returns>
        /// <remarks>
        /// u = s²/ν is formed in double-double arithmetic from the exact square s², and ln x = −ln(1 + u) carries its
        /// trailing part, so the exponent of the tail does not take the rounding of x.
        /// </remarks>
        private double AsymptoticTwoSidedTail(double s, bool logarithm)
        {
            double v = DegreesOfFreedom;
            double a = 0.5d * v;
            double s2 = s * s;
            if (a < Beta.AsymptoticMinimumShape || !(s2 < double.PositiveInfinity)) return double.NaN;
            double s2Low = Beta.ProductError(s, s, s2);
            double u = s2 / v;
            double uLow = (s2 - u * v - Beta.ProductError(u, v, u * v) + s2Low) / v;
            double logX = -Tools.Log1p(u);
            if (logX < Beta.AsymptoticMinimumLogVariable) return double.NaN;
            return Beta.AsymptoticIncomplete(a, 0.5d, logX, -uLow / (1.0d + u), logarithm);
        }

        /// <summary>
        /// Returns the log CDF of the standard Student's t distribution with ν degrees of freedom.
        /// </summary>
        /// <param name="z">The standardized variate, not NaN.</param>
        /// <returns>The natural logarithm of the lower-tail probability at <paramref name="z"/>.</returns>
        private double StandardLogCdf(double z)
        {
            if (z > 0d) return Tools.Log1p(-StandardCdf(-z));
            double probability = StandardCdf(z);
            if (probability >= 1E-290) return Math.Log(probability);
            if (double.IsNegativeInfinity(z)) return double.NegativeInfinity;
            return StandardLogLowerTail(z);
        }

        /// <summary>
        /// Returns ln F(z) = ln(I_y(ν/2, 1/2)/2), y = ν/(ν + z²), for a far lower-tail standardized variate.
        /// </summary>
        /// <param name="z">The standardized variate, at most −1.</param>
        /// <returns>The natural logarithm of the lower-tail probability at <paramref name="z"/>.</returns>
        /// <remarks>
        /// Uses the leading factor y^a (1 − y)^b / (a B(a, b)) of the incomplete beta function in log space, with
        /// a = ν/2 and b = 1/2, times the continued fraction <see cref="Beta.Incbcf(double, double, double)"/> or
        /// <see cref="Beta.Incbd(double, double, double)"/> selected as in <see cref="Beta.Incomplete(double, double, double)"/>.
        /// For z ≤ −1, y ≤ a/(a + b), so no tail reflection is needed. ln y is formed from ln(√ν/|z|), so the result
        /// stays finite when y underflows. Where the asymptotic expansion of <see cref="CDF(double)"/> applies, its
        /// logarithmic form is used instead.
        /// </remarks>
        private double StandardLogLowerTail(double z)
        {
            double logTail = AsymptoticTwoSidedTail(-z, true);
            if (!double.IsNaN(logTail)) return logTail - Tools.Log2;
            double a = DegreesOfFreedom / 2.0d;
            const double b = 0.5d;
            double logR = 0.5d * Math.Log(DegreesOfFreedom) - Math.Log(-z);
            double r2 = Math.Exp(2.0d * logR);
            double logOneMinusY = -Tools.Log1p(r2);
            double logY = 2.0d * logR + logOneMinusY;
            double y = Math.Exp(logY);
            double fraction = y * (a + b - 2.0d) - (a - 1.0d) < 0.0d ? Beta.Incbcf(a, b, y) : Beta.Incbd(a, b, y) / Math.Exp(logOneMinusY);
            return a * logY + b * logOneMinusY + Gamma.LogGamma(a + b) - Gamma.LogGamma(a) - Gamma.LogGamma(b) + Math.Log(fraction / a) - Math.Log(2.0d);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// This function is also know as the Quantile Function. A probability outside [0, 1], including NaN,
        /// throws <see cref="ArgumentOutOfRangeException"/>. For 1/4 &lt; p &lt; 3/4 the quantile follows from the inverse
        /// incomplete beta function of the two-sided centre probability |1 - 2p|, which keeps the relative accuracy of t
        /// near the median; for ν ≥ 1E16 the centre quantile equals the Normal quantile to double precision and is
        /// returned as such. In the tails with ν ≥ 30 the quantile is the root of ln CDF(t) = ln p (or of the upper tail for
        /// p &gt; 1/2), found by Brent's method between the Normal quantile and a widened bound, consistent with
        /// <see cref="CDF(double)"/>; below ν = 30 the tails use the inverse incomplete beta function of 2 min(p, 1 - p),
        /// except that the log CDF is also solved where 2 min(p, 1 - p) is subnormal, and for ν = 1 below 1E-150 the Cauchy
        /// quantile ∓1/(π min(p, 1 - p)) is returned (−∞ or +∞ beyond the double range).
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
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, DegreesOfFreedom, true);
            if (probability > 0.25d && probability < 0.75d)
            {
                if (probability == 0.5d)
                {
                    return Mu;
                }

                // Here |t| < 0.68 and t = z (1 + (z² + 1)/(4ν) + ...) with z the Normal quantile, so for ν ≥ 1E16 the
                // correction is below half an ulp, while t²/ν would lose precision as a subnormal number.
                if (DegreesOfFreedom >= 1E16d)
                {
                    return Mu + Sigma * Normal.StandardZ(probability);
                }

                double z = Beta.IncompleteInverse(0.5d, 0.5d * DegreesOfFreedom, Math.Abs(1.0d - 2.0d * probability));
                double t = Math.Sqrt(DegreesOfFreedom * z / (1.0d - z));
                if (probability < 0.5d)
                {
                    t = -t;
                }

                return Mu + Sigma * t;
            }
            else if (DegreesOfFreedom == 1d && Math.Min(probability, 1.0d - probability) < 1E-150)
            {
                // The Cauchy quantile -cot(πq) = -1/(πq) (1 - (πq)²/3 + ...), where the incomplete beta root (πq)² would
                // lose precision as a subnormal number; the correction term is below 1E-299 relative.
                double t = 1.0d / Math.PI / Math.Min(probability, 1.0d - probability);
                return Mu + Sigma * (probability < 0.5d ? -t : t);
            }
            else if (DegreesOfFreedom >= LargeDegreesOfFreedom || Math.Min(probability, 1.0d - probability) < SubnormalTail)
            {
                return Mu + Sigma * TailQuantile(probability);
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
                double t = Math.Sqrt(DegreesOfFreedom) * Math.Sqrt(1.0d - z) / Math.Sqrt(z);
                return Mu + Sigma * (rflg * t);
            }
        }

        /// <summary>
        /// The degrees of freedom from which <see cref="InverseCDF(double)"/> solves the log CDF in the tails: ν/2 is the
        /// smallest shape of the asymptotic expansion of the CDF.
        /// </summary>
        private const double LargeDegreesOfFreedom = 2.0d * Beta.AsymptoticMinimumShape;

        /// <summary>
        /// Half the smallest normal double: below it the two-sided tail 2 min(p, 1 - p) passed to the inverse incomplete
        /// beta function is subnormal and loses relative precision, so <see cref="InverseCDF(double)"/> solves the log CDF.
        /// </summary>
        private const double SubnormalTail = 1.1125369292536007E-308;

        /// <summary>
        /// Returns the standard quantile for ν ≥ 30 and p ≤ 1/4 or p ≥ 3/4 as the root of the log tail probability.
        /// </summary>
        /// <param name="probability">The probability, strictly between 0 and 1.</param>
        /// <returns>The standard quantile t with CDF(t) = <paramref name="probability"/>.</returns>
        /// <remarks>
        /// With q the smaller of p and 1 − p, the lower-tail quantile solves ln F(t) = ln q. The t tail is heavier than
        /// the Normal tail, so the root lies at or below the Normal quantile Φ⁻¹(q) &lt; 0; the lower bound is widened
        /// until ln F falls below ln q.
        /// </remarks>
        private double TailQuantile(double probability)
        {
            if (probability == 0.5d) return 0.0d;
            double q = probability < 0.5d ? probability : 1.0d - probability;
            double target = Math.Log(q);
            double upper = Normal.StandardZ(q);
            while (StandardLogCdf(upper) < target) upper *= 0.5d;
            double lower = 2.0d * upper - 1.0d;
            while (StandardLogCdf(lower) > target) lower *= 2.0d;
            double t = Brent.Solve(x => StandardLogCdf(x) - target, lower, upper, 1E-15 * Math.Abs(upper), 200, false);
            return probability < 0.5d ? t : -t;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new StudentT(Mu, Sigma, DegreesOfFreedom);
        }
              
    }
}