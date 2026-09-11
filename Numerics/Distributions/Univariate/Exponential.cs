using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The exponential distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <see href = "https://en.wikipedia.org/wiki/Exponential_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Exponential : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs an Exponential distribution with a location of 100 and scale of 10.
        /// </summary>
        public Exponential()
        {
            SetParameters(100d, 10d);
        }

        /// <summary>
        /// Constructs an Exponential distribution with a given ξ and α.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public Exponential(double location, double scale)
        {
            SetParameters(location, scale);
        }

        /// <summary>
        /// Constructs an Exponential distribution with a location of 0 and a given scale.
        /// </summary>
        /// <param name="scale">The scale parameter α (alpha).</param>
        public Exponential(double scale) : this(0.0d, scale)
        {
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
            get { return UnivariateDistributionType.Exponential; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Exponential"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "EXP"; }
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
            get { return Xi + Alpha; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Xi - Math.Log(0.5d) * Alpha; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Xi; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Alpha; }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 2.0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 9.0d; }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return Xi; }
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
            DistributionNumerics.ValidateSample(sample, 4);
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
            var newDistribution = new Exponential(Xi, Alpha);
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
            _parametersValid = ValidateParameters([location, scale], false) is null;
            _xi = location;
            _alpha = scale;
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
        /// <param name="parameters">A list of parameters.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        /// <returns><see langword="null"/> when the parameter vector is valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="throwException"/> is <see langword="true"/> and the parameter count, location, or scale is invalid.</exception>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (parameters == null || parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "Exactly two parameters are required.");
                if (throwException) throw exception;
                return exception;
            }
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
                return new ArgumentOutOfRangeException(nameof(Xi), "The location parameter ξ (Xi) must be a number.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Alpha), "The scale parameter α (alpha) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            var parms = new double[NumberOfParameters];
            parms[0] = moments[0] - moments[1];
            parms[1] = moments[1];
            return parms;
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new Exponential();
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
            double alpha = 2.0d * L2;
            double xi = L1 - alpha;
            return [xi, alpha];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            double xi = parameters[0];
            double alpha = parameters[1];
            double L1 = xi + alpha;
            double L2 = 0.5d * alpha;
            double T3 = 1d / 3d;
            double T4 = 1d / 6d;
            return [L1, L2, T3, T4];
        }

        /// <inheritdoc/>
        /// <remarks>Requires at least four finite, nonconstant observations. Preserves the legacy initialization
        /// and family-specific bounds whenever they are finite, ordered, and contain the initial values.
        /// Otherwise, the exceptional-input fallback evaluates the estimator in unit coordinates before
        /// forming finite location and positive scale bounds.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample or a representable feasible initialization is invalid.</exception>
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
            var moments = Statistics.ProductMoments(sample);
            double minData = Statistics.Minimum(sample);
            initialVals[0] = (sample.Count * minData - moments[0]) / (sample.Count - 1);
            initialVals[1] = sample.Count * (moments[0] - minData) / (sample.Count - 1);

            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = initialVals[0] - Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0]))));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[0]) + 1d)); 

            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));

            // Correct initial values if necessary
            if (initialVals[0] <= lowerVals[0] || initialVals[0] >= upperVals[0])
            {
                initialVals[0] = Statistics.Mean([lowerVals[0], upperVals[0]]);
            }
            if (initialVals[1] <= lowerVals[1] || initialVals[1] >= upperVals[1])
            {
                initialVals[1] = Statistics.Mean([lowerVals[1], upperVals[1]]);
            }
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            double normalization = DistributionNumerics.InitializationScale(sample);
            var normalized = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) normalized[i] = sample[i] / normalization;
            // The existing bias-corrected start is evaluated in unit coordinates.
            var moments = Statistics.ProductMoments(normalized);
            double minData = Statistics.Minimum(sample);
            double unitMinimum = minData / normalization;
            double n = sample.Count;
            initialVals[0] = (unitMinimum - (moments[0] - unitMinimum) / (n - 1)) * normalization;
            initialVals[1] = (n / (n - 1)) * (moments[0] - unitMinimum) * normalization;
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], initialVals[1], minData,
                Statistics.Maximum(sample), true, out lowerVals[0], out upperVals[0]);
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
                var EXP = new Exponential();
                EXP.SetParameters(x);
                return EXP.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
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
        /// </remarks>
        public override double LogPDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, Alpha], true);
            if (X < Minimum || X > Maximum) return double.NegativeInfinity;
            double lf = -Math.Log(Alpha) - DistributionNumerics.Standardize(X, Xi, Alpha);
            return double.IsNaN(lf) ? double.NegativeInfinity : lf;
        }

        /// <inheritdoc/>
        public override double CDF(double X)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, Alpha], true);
            if (X <= Minimum) return 0d;
            if (X >= Maximum) return 1d;
            return -Tools.Expm1(-DistributionNumerics.Standardize(X, Xi, Alpha));
        }

        /// <inheritdoc/>
        /// <remarks>Computes the logarithm directly without subtracting a near-unit exponential.</remarks>
        public override double LogCDF(double X)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            if (X <= Xi) return double.NegativeInfinity;
            return DistributionNumerics.Log1mExp(-DistributionNumerics.Standardize(X, Xi, Alpha));
        }

        /// <inheritdoc/>
        public override double CCDF(double X) => Math.Exp(LogCCDF(X));

        /// <inheritdoc/>
        public override double LogCCDF(double X)
        {
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            return X <= Xi ? 0 : -DistributionNumerics.Standardize(X, Xi, Alpha);
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
                ValidateParameters([Xi, Alpha], true);
            double unitQuantile = -Tools.Log1p(-probability);
            double displacement = Alpha * unitQuantile;
            return double.IsInfinity(displacement) && Tools.IsFinite(unitQuantile)
                ? Alpha * (Xi / Alpha + unitQuantile) : Xi + displacement;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Exponential(Xi, Alpha);
        }

        /// <inheritdoc/>
        /// <remarks>The actual MLE uses the sample minimum and mean minus minimum. Its covariance is
        /// diagonal with alpha squared times [1/n squared, (n-1)/n squared]. The MoM covariance is unchanged.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or sample size is not positive.</exception>
        /// <exception cref="NotImplementedException">The requested estimation method is unsupported.</exception>
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
                ValidateParameters([Xi, _alpha], true);

            // Compute covariance
            double n = sampleSize;
            double a = Alpha / Math.Sqrt(n);
            var covar = new double[2, 2];
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                covar[0, 0] = a * a; // location
                covar[1, 1] = 2d * (a * a); // scale
                covar[0, 1] = -(a * a); // location & scale
                covar[1, 0] = covar[0, 1];
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                // Actual MLE: location=min(sample), scale=mean(sample)-min(sample).
                double locationScale = Alpha / n;
                covar[0, 0] = locationScale * locationScale;
                covar[1, 1] = (a * a) * ((n - 1) / n);
            }
            return covar;
        }

        /// <inheritdoc/>
        /// <remarks>Uses the same covariance quadratic form in normalized coordinates; probability must be finite and strictly between zero and one.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters([Xi, Alpha], true);
            var unit = new Exponential(0, 1);
            return DistributionNumerics.ScaledQuantileVariance(unit.ParameterCovariance(sampleSize, estimationMethod),
                unit.QuantileGradient(probability), Alpha);
        }

        /// <inheritdoc/>
        /// <remarks>Returns the derivative of the actual inverse CDF in location and scale coordinates.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Xi, _alpha], true);
            var gradient = new double[]
            {
                1.0d, // location
                -Tools.Log1p(-probability) // scale
            };
            return gradient;
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

            double xi = Xi;
            double alpha = Alpha;
            if (!(alpha > 0.0))
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Map to Y = X - xi, truncated on (A, B)
            // Note: support is y >= 0
            double A = Math.Max(0.0, a - xi);
            double B = b - xi;

            if (double.IsNaN(A) || double.IsNaN(B) || B <= 0.0) // interval entirely left of support or invalid
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Standardized limits t = y/alpha
            double tA = A / alpha;
            double tB = double.IsPositiveInfinity(B) ? double.PositiveInfinity : (B / alpha);

            // Normalizing probability Z = P(A < Y < B) = e^{-tA} - e^{-tB}
            double eA = Math.Exp(-tA);
            double eB = double.IsPositiveInfinity(tB) ? 0.0 : Math.Exp(-tB);
            double Z = eA - eB;
            if (Z <= 1e-15)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            // Helper: S_n(t) = e^{-t} * sum_{k=0}^n t^k/k!  (appears in γ(n+1, t) = n! * (1 - S_n(t)))
            static double Sn(double t, int n)
            {
                if (double.IsPositiveInfinity(t)) return 0.0;
                double term = 1.0; // t^0/0!
                double sum = term;
                for (int k = 1; k <= n; k++)
                {
                    term *= t / k;    // t^k/k!
                    sum += term;
                }
                return Math.Exp(-t) * sum;
            }

            // Precompute factorials for n = 0..4
            double[] fact = { 1.0, 1.0, 2.0, 6.0, 24.0 };

            // Raw truncated moments of Y: E[Y^n | A<Y<B] for n = 0..4
            // Using: ∫_A^B y^n (1/α) e^{-y/α} dy = α^n * [γ(n+1, B/α) - γ(n+1, A/α)]
            // and γ(n+1, t) = n! * (1 - S_n(t)).
            double[] EY = new double[5];
            EY[0] = 1.0; // by definition under conditioning

            for (int n = 1; n <= 4; n++)
            {
                double SnA = Sn(tA, n);
                double SnB = Sn(tB, n);
                double numer = Math.Pow(alpha, n) * fact[n] * (SnA - SnB); // α^n n! [S_n(tA) - S_n(tB)]
                EY[n] = numer / Z;
            }

            // Convert to raw moments of X via binomial expansion: E[X^k] = sum_{r=0}^k C(k,r) xi^(k-r) E[Y^r]
            double[] EX = new double[5];
            for (int k = 0; k <= 4; k++)
            {
                double sum = 0.0;
                for (int r = 0; r <= k; r++)
                {
                    double bc = Mathematics.SpecialFunctions.Factorial.BinomialCoefficient(k, r);
                    sum += bc * Math.Pow(xi, k - r) * EY[r];
                }
                EX[k] = sum;
            }

            // Central moments about the *unconditional* mean μ = xi + alpha
            double mu = xi + alpha;
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
