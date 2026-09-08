using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Generalized Extreme Value distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class GeneralizedExtremeValue : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
    
        /// <summary>
        /// Constructs a Generalized Extreme Value with a location of 100, scale of 10, and shape of 0.
        /// </summary>
        public GeneralizedExtremeValue()
        {
            SetParameters(100d, 10d, 0d);
        }

        /// <summary>
        /// Constructs a Generalized Extreme Value (GEV) distribution with the given parameters ξ, α, and κ.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GeneralizedExtremeValue(double location, double scale, double shape)
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
            get { return UnivariateDistributionType.GeneralizedExtremeValue; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Extreme Value"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GEV"; }
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
        /// <remarks>The mean exists for kappa &gt; -1. Log-Gamma arithmetic preserves finite scaled
        /// values even when a raw Gamma function overflows.</remarks>
        public override double Mean
        {
            get
            {
                if (Kappa <= -1) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return Xi + Alpha * SmallShapeStandardizedMoments()[0];
                double logarithm = Gamma.LogGamma(1 + Kappa);
                if (logarithm == 0) return Xi;
                double magnitude = logarithm > 0 ? logarithm + DistributionNumerics.Log1mExp(-logarithm)
                    : DistributionNumerics.Log1mExp(logarithm);
                return Xi - Math.Sign(Kappa) * Math.Sign(logarithm)
                    * Math.Exp(Math.Log(Alpha) + magnitude - Math.Log(Math.Abs(Kappa)));
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(.5); }
        }

        /// <inheritdoc/>
        /// <remarks>The mode is the finite upper endpoint when kappa is at least one; otherwise it is the interior stationary point.</remarks>
        public override double Mode
        {
            get
            {
                if (Kappa >= 1) return Maximum;
                double logarithm = Tools.Log1p(-Kappa);
                return Xi - Alpha * logarithm * DistributionNumerics.Exprel(Kappa * logarithm);
            }
        }

        /// <inheritdoc/>
        /// <remarks>The variance exists for kappa &gt; -1/2. Nonexistent moments return NaN.</remarks>
        public override double StandardDeviation
        {
            get
            {
                if (Kappa <= -.5) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return Alpha * SmallShapeStandardizedMoments()[1];
                if (Kappa == 1) return Alpha;
                double logVariance = LogPowerVariance(Kappa);
                return Math.Exp(Math.Log(Alpha) + .5 * logVariance - Math.Log(Math.Abs(Kappa)));
            }
        }

        /// <inheritdoc/>
        /// <remarks>The third moment exists for kappa &gt; -1/3; positive bounded shapes are not excluded.</remarks>
        public override double Skewness
        {
            get
            {
                if (Kappa <= -1d / 3d) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return SmallShapeStandardizedMoments()[2];
                if (Kappa == 1) return -2;
                double l1 = Gamma.LogGamma(1 + Kappa), l2 = Gamma.LogGamma(1 + 2 * Kappa), l3 = Gamma.LogGamma(1 + 3 * Kappa);
                if (double.IsPositiveInfinity(l3)) return double.NegativeInfinity;
                double largest = Math.Max(l3, Math.Max(l1 + l2, 3 * l1));
                double centered = Math.Exp(l3 - largest) - 3 * Math.Exp(l1 + l2 - largest) + 2 * Math.Exp(3 * l1 - largest);
                return centered == 0 ? 0 : -Math.Sign(Kappa) * Math.Sign(centered)
                    * Math.Exp(largest + Math.Log(Math.Abs(centered)) - 1.5 * LogPowerVariance(Kappa));
            }
        }

        /// <inheritdoc/>
        /// <remarks>Ordinary kurtosis requires kappa &gt; -1/4. Analytical normalized central moments are used without raw-moment overflow.</remarks>
        public override double Kurtosis
        {
            get
            {
                if (Kappa <= -.25) return double.NaN;
                if (Math.Abs(Kappa) <= .05) return SmallShapeStandardizedMoments()[3];
                if (Kappa == 1) return 9;
                double l1 = Gamma.LogGamma(1 + Kappa), l2 = Gamma.LogGamma(1 + 2 * Kappa);
                double l3 = Gamma.LogGamma(1 + 3 * Kappa), l4 = Gamma.LogGamma(1 + 4 * Kappa);
                if (double.IsPositiveInfinity(l4)) return double.PositiveInfinity;
                double largest = Math.Max(Math.Max(l4, l1 + l3), Math.Max(2 * l1 + l2, 4 * l1));
                double centered = Math.Exp(l4 - largest) - 4 * Math.Exp(l1 + l3 - largest)
                    + 6 * Math.Exp(2 * l1 + l2 - largest) - 3 * Math.Exp(4 * l1 - largest);
                return Math.Exp(largest + Math.Log(centered) - 2 * LogPowerVariance(Kappa));
            }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get
            {
                if (Kappa >= 0)
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
                if (Kappa <= 0)
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
            DistributionNumerics.ValidateSample(sample, 4);
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
            var newDistribution = new GeneralizedExtremeValue(Xi, Alpha, Kappa);
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
            if (parameters == null || parameters.Count != NumberOfParameters)
                throw new ArgumentOutOfRangeException(nameof(parameters), "Exactly three parameters are required.");
            SetParameters(parameters[0], parameters[1], parameters[2]);
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException? ValidateParameters(double location, double scale, double shape, bool throwException)
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
            return null!;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "Exactly three parameters are required.");
                if (throwException) throw exception;
                return exception;
            }
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
                a = Math.Sqrt(6d) / Math.PI * moments[1];
                x = moments[0] - a * Tools.Euler;
            }
            else
            {
                double U1 = Gamma.Function(1d + k);
                double U2 = Gamma.Function(1d + 2d * k);
                a = Math.Sqrt(moments[1] * moments[1] * k * k / (U2 - Math.Pow(U1, 2d)));
                x = moments[0] - a / k * (1d - U1);
            }

            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            // Solve for kappa
            double k = SolveForKappa(moments[2]);
            double a;
            double x;
            if (Math.Abs(k) <= NearZero)
            {
                a = Math.Sqrt(6d) / Math.PI * moments[1];
                x = moments[0] - a * Tools.Euler;
            }
            else
            {
                double U1 = Gamma.Function(1d + k);
                double U2 = Gamma.Function(1d + 2d * k);
                a = Math.Sqrt(moments[1] * moments[1] * k * k / (U2 - Math.Pow(U1, 2d)));
                x = moments[0] - a / k * (1d - U1);
            }

            // return parameters
            return [x, a, k];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new GeneralizedExtremeValue();
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
            if (skew > 1.14d && skew < 10d)
            {
                // Extreme Value Type II
                return 0.2858221d - 0.357983d * skew + 0.116659d * Math.Pow(skew, 2d) - 0.022725d * Math.Pow(skew, 3d) + 0.002604d * Math.Pow(skew, 4d) - 0.000161d * Math.Pow(skew, 5d) + 0.000004d * Math.Pow(skew, 6d);
            }
            else if (skew == 1.14d)
            {
                // Extreme Value Type I
                return 0d;
            }
            else if (skew >= 0d && skew < 1.14d)
            {
                // Extreme Value Type III
                // This regression equation works for -2 < Cs < 1.14
                return 0.277648d - 0.322016d * skew + 0.060278d * Math.Pow(skew, 2d) + 0.016759d * Math.Pow(skew, 3d) - 0.005873d * Math.Pow(skew, 4d) - 0.00244d * Math.Pow(skew, 5d) - 0.00005d * Math.Pow(skew, 6d);
            }
            else if (skew < 0d && skew >= -2)
            {
                // For negative values of skew, there exists two possible values of kappa.
                // Either extreme value type II or III
                // Therefor it must be solved. The Brent method is used here. 
                return Brent.Solve((x) =>
                {
                    double U1 = Gamma.Lanczos(1d + x);
                    double U2 = Gamma.Lanczos(1d + 2d * x);
                    double U3 = Gamma.Lanczos(1d + 3d * x);
                    double k = Math.Sign(x) * (-U3 + 3d * U1 * U2 - 2d * Math.Pow(U1, 3d)) / Math.Pow(U2 - Math.Pow(U1, 2d), 3d / 2d);
                    return k - skew;
                }, -(1d / 3d), 1d);
            }
            else if (skew < -2)
            {
                // EV2(3)
                // This regression equation works for -10 < Cs < 0
                return -0.50405d - 0.00861d * skew + 0.015497d * Math.Pow(skew, 2d) + 0.005613d * Math.Pow(skew, 3d) + 0.00087d * Math.Pow(skew, 4d) + 0.000065d * Math.Pow(skew, 5d);
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
            // The following approximation given by Hosking et al. (1985b) has accuracy better than 9x10-4 for abs(t3)<=0.5.
            if (Math.Abs(T3) <= 0.5d)
            {
                double c = 2d / (3d + T3) - Math.Log(2d) / Math.Log(3d);
                double kappa = 7.859d * c + 2.9554d * Math.Pow(c, 2d);
                double alpha = L2 * kappa / ((1d - Math.Pow(2d, -kappa)) * Gamma.Function(1d + kappa));
                double xi = L1 - alpha * (1d - Gamma.Function(1d + kappa)) / kappa;
                return [xi, alpha, kappa];
            }
            else
            {
                // Solve for kappa
                double kappa = Brent.Solve(x => T3 - (2.0d * (1.0d - Math.Pow(3.0d, -x)) / (1.0d - Math.Pow(2.0d, -x)) - 3.0d), -1, 10d);
                double alpha = L2 * kappa / ((1d - Math.Pow(2d, -kappa)) * Gamma.Function(1d + kappa));
                double xi = L1 - alpha * (1d - Gamma.Function(1d + kappa)) / kappa;
                return [xi, alpha, kappa];
            }
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            if (kappa <= -1.0d)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "L-moments can only be defined for kappa > -1.");
            double L1 = xi + alpha * (1.0d - Gamma.Function(1.0d + kappa)) / kappa;
            double L2 = alpha * (1.0d - Math.Pow(2.0d, -kappa)) * Gamma.Function(1.0d + kappa) / kappa;
            double T3 = 2.0d * (1.0d - Math.Pow(3.0d, -kappa)) / (1.0d - Math.Pow(2.0d, -kappa)) - 3.0d;
            double T4 = (5.0d * (1.0d - Math.Pow(4.0d, -kappa)) - 10.0d * (1.0d - Math.Pow(3.0d, -kappa)) + 6.0d * (1.0d - Math.Pow(2.0d, -kappa))) / (1.0d - Math.Pow(2.0d, -kappa));
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        /// <remarks>Requires at least four finite, nonconstant observations. Initialization uses
        /// the existing linear-moment estimator and preserves the fitting shape bounds.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample or a finite feasible initialization is invalid.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            double normalization = DistributionNumerics.InitializationScale(sample);
            var normalized = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) normalized[i] = sample[i] / normalization;
            var initialVals = ParametersFromLinearMoments(Statistics.LinearMoments(normalized));
            initialVals[0] *= normalization;
            initialVals[1] *= normalization;
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], initialVals[1], Statistics.Minimum(sample),
                Statistics.Maximum(sample), false, out lowerVals[0], out upperVals[0]);
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
            // Get bounds of shape
            lowerVals[2] = -10;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
            if (!DistributionNumerics.IsFinite(initialVals[2]) || initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
            // 
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
                var GEV = new GeneralizedExtremeValue();
                GEV.SetParameters(x);
                return GEV.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;

        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            return Math.Exp(LogPDF(x));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated in log space, so far-tail densities that underflow <see cref="PDF(double)"/>
        /// keep a finite log density.
        /// The exact shape controls support. At the finite upper endpoint, shape one has density
        /// 1/alpha and shapes above one have a genuine integrable infinite density.
        /// </remarks>
        public override double LogPDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, Alpha, Kappa, true);
            if (x < Minimum || x > Maximum || double.IsInfinity(x)) return double.NegativeInfinity;
            if (Kappa > 0 && x == Maximum)
                return Kappa < 1 ? double.NegativeInfinity : Kappa == 1 ? -Math.Log(Alpha) : double.PositiveInfinity;
            if (Kappa < 0 && x == Minimum) return double.NegativeInfinity;
            double y = TransformedValue(x);
            double lf = -(1d - Kappa) * y - Math.Exp(-y) - Math.Log(Alpha);
            return double.IsNaN(lf) ? double.NegativeInfinity : lf;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return Math.Exp(LogCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0;
            return -Math.Exp(-TransformedValue(x));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => -Tools.Expm1(LogCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return 0;
            if (x >= Maximum) return double.NegativeInfinity;
            double y = TransformedValue(x);
            double exponential = Math.Exp(-y);
            return exponential == 0 ? -y : DistributionNumerics.Log1mExp(-exponential);
        }

        /// <summary>Maps an interior observation to its Gumbel coordinate using the exact nonzero shape.</summary>
        private double TransformedValue(double x)
        {
            return DistributionNumerics.HoskingShapeTransform(x, Xi, Alpha, Kappa);
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
                ValidateParameters(Xi, Alpha, Kappa, true);
            double logarithm = Math.Log(-Math.Log(probability));
            double product = Kappa * logarithm;
            double unitQuantile = double.IsNegativeInfinity(product) ? 1 / Kappa : DistributionNumerics.ScaledExprelProduct(1, -logarithm, product);
            double displacement = double.IsNegativeInfinity(product) ? Alpha / Kappa : DistributionNumerics.ScaledExprelProduct(Alpha, -logarithm, product);
            return double.IsInfinity(displacement) && DistributionNumerics.IsFinite(unitQuantile)
                ? Alpha * (Xi / Alpha + unitQuantile) : Xi + displacement;
        }

        /// <summary>
        /// Gets the expected Fisher information matrix.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <returns>The full three-parameter expected information in location, scale and shape coordinates.</returns>
        /// <remarks>Regular information requires kappa &lt; 1/2. At zero shape the shape parameter
        /// remains estimated, so this is not the two-parameter Gumbel information.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Sample size, distribution parameters or regularity are invalid.</exception>
        /// <exception cref="InvalidOperationException">Information cannot be resolved or represented numerically.</exception>
        public Matrix ExpectedInformationMatrix(int sampleSize)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            var information = KappaExpectedInformation.ExpectedInformation(Kappa, 0, 3, out _, out _, out _);
            for (int i = 0; i < 3; i++)
            for (int j = i; j < 3; j++)
            {
                double value = information[i, j];
                if (value != 0)
                {
                    double logarithm = Math.Log(Math.Abs(value)) + Math.Log(sampleSize)
                        - (i < 2 ? Math.Log(Alpha) : 0) - (j < 2 ? Math.Log(Alpha) : 0);
                    value = Math.Sign(value) * Math.Exp(logarithm);
                    if (!DistributionNumerics.IsFinite(value))
                        throw new InvalidOperationException("Expected information is outside the finite floating-point range.");
                }
                information[i, j] = information[j, i] = value;
            }
            return new Matrix(information);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedExtremeValue(Xi, Alpha, Kappa);
        }

        /// <inheritdoc/>
        /// <remarks>Full three-parameter MLE covariance requires positive sample size and kappa &lt; 1/2.
        /// The zero-shape limit still estimates shape; uncertainty regularity does not restrict the distribution's valid shape domain.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters, sample size or information regularity are invalid.</exception>
        /// <exception cref="NotImplementedException">The method is not maximum likelihood.</exception>
        /// <exception cref="InvalidOperationException">The information or covariance is numerically unresolved.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            return KappaExpectedInformation.ParameterCovariance(Alpha, Kappa, 0, sampleSize, 3);
        }

        /// <inheritdoc/>
        /// <remarks>Analytical exponential divided differences retain the exact nonzero shape and
        /// the full [1,-L,-alpha*L squared/2] gradient at zero shape, where L=log(-log(p)).</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            double logarithm = Math.Log(-Math.Log(probability));
            double product = Kappa * logarithm;
            var gradient = new double[]
            {
                1.0d, // location
                double.IsNegativeInfinity(product) ? 1 / Kappa : DistributionNumerics.ScaledExprelProduct(1, -logarithm, product), // scale
                double.IsNegativeInfinity(product) ? -(Alpha / Kappa) / Kappa : -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, logarithm, product) // shape
            };
            return gradient;
        }

        /// <inheritdoc/>
        /// <remarks>Uses finite physical common-coordinate gradients and unit-scale covariance,
        /// permitting scalar variance evaluation independently of the physical covariance matrix's range.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            var unit = new GeneralizedExtremeValue(0, 1, Kappa);
            double logarithm = Math.Log(-Math.Log(probability)), product = Kappa * logarithm;
            double scaleGradient = DistributionNumerics.ScaledExprelProduct(Alpha, -logarithm, product);
            double shapeGradient = -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, logarithm, product);
            return DistributionNumerics.ScaledQuantileVariance(unit.ParameterCovariance(sampleSize, estimationMethod),
                [Alpha, scaleGradient, shapeGradient]);
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

        /// <summary>Returns log Var(T^power), T unit exponential, without forming raw gamma moments.</summary>
        /// <remarks>The caller establishes power &gt; -1/2 and uses a divided series near zero.</remarks>
        internal static double LogPowerVariance(double power)
        {
            double first = Gamma.LogGamma(1 + power), second = Gamma.LogGamma(1 + 2 * power);
            return double.IsPositiveInfinity(second) ? double.PositiveInfinity
                : second + DistributionNumerics.Log1mExp(2 * first - second);
        }

        /// <summary>Evaluates analytical standardized central moments with log-Gamma divided differences near zero shape.</summary>
        /// <remarks>Finite differences of log Gamma remove the cancelling powers before evaluation.
        /// Exponential polynomial identities then retain the second through fourth centered moments.
        /// The exact zero is the analytical Gumbel limit; every nonzero kappa remains in the series.</remarks>
        private double[] SmallShapeStandardizedMoments()
        {
            double k = Kappa;
            if (k == 0) return [Tools.Euler, Math.PI / Math.Sqrt(6), 1.1395470994046487, 5.4];
            double logarithm = DistributionNumerics.LogGammaOnePlus(k);
            double mean = -(logarithm / k) * DistributionNumerics.Exprel(logarithm);
            double a = NormalizedLogGammaDifference(k, 2);
            double b = NormalizedLogGammaDifference(k, 3);
            double c = NormalizedLogGammaDifference(k, 4);
            double a2 = k * k * a, b3 = k * k * k * b, c4 = k * k * k * k * c;
            double u = Tools.Expm1(a2), v = Tools.Expm1(b3);
            double variance = a * DistributionNumerics.Exprel(a2);
            double third = b * DistributionNumerics.Exprel(b3);
            double fourth = c * DistributionNumerics.Exprel(c4);
            double skew = -(k * variance * variance * (3 + u) + Math.Exp(3 * a2) * third)
                / (variance * Math.Sqrt(variance));
            double kurtosis = (variance * variance * (3 + u * (16 + u * (15 + u * (6 + u))))
                + 12 * k * third * Math.Exp(3 * a2) * a * DistributionNumerics.Exprel(3 * a2)
                + Math.Exp(6 * a2) * (k * k * third * third * (6 + v * (4 + v)) + Math.Exp(4 * b3) * fourth))
                / (variance * variance);
            return [mean, Math.Exp(logarithm) * Math.Sqrt(variance), skew, kurtosis];
        }

        /// <summary>Returns the order-r forward difference of log Gamma(1+j*kappa), divided by kappa to power r.</summary>
        private static double NormalizedLogGammaDifference(double kappa, int order)
        {
            double sum = 0, power = 1;
            for (int n = order; n <= 32; n++)
            {
                double factor = order == 2 ? Math.Pow(2, n) - 2
                    : order == 3 ? Math.Pow(3, n) - 3 * Math.Pow(2, n) + 3
                    : Math.Pow(4, n) - 4 * Math.Pow(3, n) + 6 * Math.Pow(2, n) - 4;
                sum += (n % 2 == 0 ? 1 : -1) * DistributionNumerics.ZetaInteger(n) * factor * power / n;
                power *= kappa;
            }
            return sum;
        }

    }
}
