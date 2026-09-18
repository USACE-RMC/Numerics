using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.RootFinding;

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
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                if (Math.Abs(Kappa) >= 1) return double.NaN;
                if (Kappa == 0) return Xi;
                return Xi + (Math.Abs(Kappa) <= .05 ? -Alpha * Kappa * Polynomial(ReciprocalCoefficients, Kappa * Kappa, 1) : Alpha * StandardMean(Kappa));
            }
        }

        /// <inheritdoc/>
        public override double Median => InverseCDF(.5);

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                if (Kappa <= -1) return Minimum;
                if (Kappa >= 1) return Maximum;
                double z = Tools.Log1p(Kappa) - Tools.Log1p(-Kappa);
                return QuantileAtLatent(z);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                if (Math.Abs(Kappa) >= .5) return double.NaN;
                return Alpha * Math.Sqrt(StandardVariance(Kappa));
            }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get
            {
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                if (Math.Abs(Kappa) >= 1d / 3) return double.NaN;
                if (Kappa == 0) return 0;
                if (Math.Abs(Kappa) <= .05)
                    return Kappa * Polynomial(ThirdCoefficients, Kappa * Kappa, 2) / Math.Pow(StandardVariance(Kappa), 1.5);
                double b1 = ReciprocalSinc(Kappa), b2 = ReciprocalSinc(2 * Kappa), b3 = ReciprocalSinc(3 * Kappa);
                return Math.Sign(Kappa) * (-b3 + 3 * b1 * b2 - 2 * b1 * b1 * b1) / Math.Pow(b2 - b1 * b1, 1.5);
            }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get
            {
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                if (Math.Abs(Kappa) >= .25) return double.NaN;
                if (Kappa == 0) return 21d / 5;
                double variance = StandardVariance(Kappa);
                if (Math.Abs(Kappa) <= .05)
                    return Polynomial(FourthCoefficients, Kappa * Kappa, 2) / (variance * variance);
                double b1 = ReciprocalSinc(Kappa), b2 = ReciprocalSinc(2 * Kappa);
                double b3 = ReciprocalSinc(3 * Kappa), b4 = ReciprocalSinc(4 * Kappa);
                double v = b2 - b1 * b1;
                return (b4 - 4 * b1 * b3 + 6 * b1 * b1 * b2 - 3 * b1 * b1 * b1 * b1) / (v * v);
            }
        }

        private static readonly double[] ReciprocalCoefficients = BuildReciprocalCoefficients();
        private static readonly double[] VarianceCoefficients = BuildMomentCoefficients(2);
        private static readonly double[] ThirdCoefficients = BuildMomentCoefficients(3);
        private static readonly double[] FourthCoefficients = BuildMomentCoefficients(4);

        /// <summary>Coefficients of pi*k/sin(pi*k) as a power series in k squared.</summary>
        /// <returns>The reciprocal-sinc series coefficients indexed by powers of squared shape.</returns>
        private static double[] BuildReciprocalCoefficients()
        {
            var sinc = new double[15];
            var reciprocal = new double[15];
            sinc[0] = reciprocal[0] = 1;
            for (int n = 1; n < sinc.Length; n++)
            {
                sinc[n] = -sinc[n - 1] * Math.PI * Math.PI / (2 * n) / (2 * n + 1);
                for (int j = 1; j <= n; j++) reciprocal[n] -= sinc[j] * reciprocal[n - j];
            }
            return reciprocal;
        }

        /// <summary>Multiplies truncated power series used to remove exact central-moment zeros algebraically.</summary>
        /// <param name="left">The first coefficient vector.</param>
        /// <param name="right">The second coefficient vector of the same length.</param>
        /// <returns>The product truncated to the input vector length.</returns>
        private static double[] Multiply(double[] left, double[] right)
        {
            var result = new double[left.Length];
            for (int n = 0; n < result.Length; n++)
            for (int j = 0; j <= n; j++) result[n] += left[j] * right[n - j];
            return result;
        }

        /// <summary>Forms central-moment numerator coefficients before evaluation, avoiding cancellation near zero shape.</summary>
        /// <param name="order">The central-moment order, expected to be two, three, or four.</param>
        /// <returns>The numerator coefficients indexed by powers of squared shape.</returns>
        private static double[] BuildMomentCoefficients(int order)
        {
            var b1 = ReciprocalCoefficients;
            var b2 = new double[b1.Length];
            var b3 = new double[b1.Length];
            var b4 = new double[b1.Length];
            for (int n = 0; n < b1.Length; n++)
            {
                b2[n] = b1[n] * Math.Pow(4, n);
                b3[n] = b1[n] * Math.Pow(9, n);
                b4[n] = b1[n] * Math.Pow(16, n);
            }
            double[] b11 = Multiply(b1, b1), b12 = Multiply(b1, b2), b111 = Multiply(b11, b1);
            double[] b13 = Multiply(b1, b3), b112 = Multiply(b11, b2), b1111 = Multiply(b11, b11);
            var result = new double[b1.Length];
            for (int n = 0; n < result.Length; n++)
                result[n] = order == 2 ? b2[n] - b11[n]
                    : order == 3 ? -b3[n] + 3 * b12[n] - 2 * b111[n]
                    : b4[n] - 4 * b13[n] + 6 * b112[n] - 3 * b1111[n];
            return result;
        }

        /// <summary>Horner evaluation after dividing out the exact leading power of kappa squared.</summary>
        /// <param name="coefficients">The power-series coefficients.</param>
        /// <param name="squaredShape">The squared shape at which to evaluate the reduced series.</param>
        /// <param name="first">The first coefficient retained after removing the exact leading zero.</param>
        /// <returns>The reduced polynomial value.</returns>
        private static double Polynomial(double[] coefficients, double squaredShape, int first)
        {
            double value = 0;
            for (int n = coefficients.Length - 1; n >= first; n--) value = value * squaredShape + coefficients[n];
            return value;
        }

        /// <summary>Returns pi*k/sin(pi*k), including its removable singularity.</summary>
        /// <param name="k">The generalized-logistic shape.</param>
        /// <returns><c>pi*k/sin(pi*k)</c>, using its power series near zero.</returns>
        private static double ReciprocalSinc(double k) => Math.Abs(k) <= .05
            ? 1 + k * k * Polynomial(ReciprocalCoefficients, k * k, 1) : Math.PI * k / Math.Sin(Math.PI * k);

        /// <summary>Returns the standardized mean shift without subtracting nearly equal raw moments.</summary>
        /// <param name="k">The generalized-logistic shape.</param>
        /// <returns>The standardized mean displacement from location.</returns>
        private static double StandardMean(double k) => Math.Abs(k) <= .05
            ? -k * Polynomial(ReciprocalCoefficients, k * k, 1) : (1 - ReciprocalSinc(k)) / k;

        /// <summary>Returns standardized variance after dividing out its exact kappa-squared zero.</summary>
        /// <param name="k">The generalized-logistic shape.</param>
        /// <returns>The standardized variance.</returns>
        private static double StandardVariance(double k) => Math.Abs(k) <= .05
            ? Polynomial(VarianceCoefficients, k * k, 1)
            : (ReciprocalSinc(2 * k) - Math.Pow(ReciprocalSinc(k), 2)) / k / k;
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
                    return FiniteShapeEndpoint();
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
                    return FiniteShapeEndpoint();
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
            return ValidateParameters(parameters[0], parameters[1], parameters[2], throwException);
        }

        /// <summary>
        /// Gets the parameters using the direct method of moments. Moments are derived from the real-space data.
        /// </summary>
        /// <param name="moments">The array of sample moments.</param>
        /// <returns>The location, scale, and shape derived from the supplied product moments.</returns>
        /// <exception cref="ArgumentNullException"><paramref name="moments"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The moment vector is too short or contains an invalid mean, dispersion, or skewness.</exception>
        public double[] DirectMethodOfMoments(IList<double> moments)
        {
            if (moments == null) throw new ArgumentNullException(nameof(moments));
            if (moments.Count < 3 || !Tools.IsFinite(moments[0]) || !Tools.IsFinite(moments[1])
                || moments[1] <= 0 || !Tools.IsFinite(moments[2]))
                throw new ArgumentOutOfRangeException(nameof(moments));
            double k = SolveForKappa(moments[2]);
            double a = moments[1] / Math.Sqrt(StandardVariance(k));
            double x = moments[0] - a * StandardMean(k);
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
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="skew"/> is nonfinite.</exception>
        public double SolveForKappa(double skew)
        {
            if (!Tools.IsFinite(skew)) throw new ArgumentOutOfRangeException(nameof(skew));
            if (skew == 0) return 0;
            if (Math.Abs(skew) >= 10) return double.NaN;
            // Retain Brent and its convergence settings, evaluating finite moment values inside the open moment domain.
            double bound = 1d / 3 - Tools.DoubleMachineEpsilon;
            return Brent.Solve(k => new GeneralizedLogistic(0, 1, k).Skewness - skew, -bound, bound);
        }
        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            if (moments == null) throw new ArgumentNullException(nameof(moments));
            if (moments.Count < 3 || !Tools.IsFinite(moments[0]) || !Tools.IsFinite(moments[1])
                || moments[1] <= 0 || !Tools.IsFinite(moments[2]) || Math.Abs(moments[2]) >= 1)
                throw new ArgumentOutOfRangeException(nameof(moments));
            double kappa = -moments[2];
            double alpha = moments[1] / ReciprocalSinc(kappa);
            double xi = moments[0] - alpha * StandardMean(kappa);
            return [xi, alpha, kappa];
        }
        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0], alpha = parameters[1], kappa = parameters[2];
            ValidateParameters(xi, alpha, kappa, true);
            if (Math.Abs(kappa) >= 1) throw new ArgumentOutOfRangeException(nameof(parameters), "L-moments require absolute kappa below one.");
            return [xi + alpha * StandardMean(kappa), alpha * ReciprocalSinc(kappa), -kappa, (1 + 5 * kappa * kappa) / 6];
        }
        /// <inheritdoc/>
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
            // Estimate initial values using the method of moments (a.k.a product moments).
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            // initialVals = DirectMethodOfMoments(Statistics.ComputeProductMoments(sample))
            initialVals = LegacyConstraintParametersFromLinearMoments(Statistics.LinearMoments(sample));
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

        /// <summary>Retains the established constraint initializer arithmetic for ordinary samples.</summary>
        /// <param name="moments">Sample moments.</param>
        /// <returns>The legacy initial parameter values.</returns>
        private double[] LegacyConstraintParametersFromLinearMoments(IList<double> moments)
        {
            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];
            double T4 = moments[3];
            double kappa = -T3;
            double alpha;
            double xi;
            if (kappa == 0.0d)
            {
                alpha = L2;
                xi = L1;
            }
            else if (Math.Abs(kappa) <= NearZero)
            {
                double kappa2 = kappa * kappa;
                double pi2 = Math.PI * Math.PI;
                double sinc = 1.0d - pi2 * kappa2 / 6.0d + pi2 * pi2 * kappa2 * kappa2 / 120.0d;
                double reciprocalDifference = -pi2 * kappa / 6.0d - 7.0d * pi2 * pi2 * kappa * kappa2 / 360.0d;
                alpha = L2 * sinc;
                xi = L1 - alpha * reciprocalDifference;
            }
            else
            {
                alpha = L2 * Math.Sin(kappa * Math.PI) / (kappa * Math.PI);
                xi = L1 - alpha * (1.0d / kappa - Math.PI / Math.Sin(kappa * Math.PI));
            }
            return [xi, alpha, kappa];
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The sample is invalid, insufficient, or constant.</exception>
        /// <exception cref="InvalidOperationException">No finite supported initializer can be placed within finite bounds.</exception>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            // Estimate initial values using the method of moments (a.k.a product moments).
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            // initialVals = DirectMethodOfMoments(Statistics.ComputeProductMoments(sample))
            double magnitude = 0;
            foreach (double value in sample) magnitude = Math.Max(magnitude, Math.Abs(value));
            var scaled = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) scaled[i] = sample[i] / magnitude;
            double[] moments = Statistics.LinearMoments(scaled);
            initialVals = ParametersFromLinearMoments(moments);
            initialVals[0] *= magnitude;
            initialVals[1] *= magnitude;
            var candidate = new GeneralizedLogistic(initialVals[0], initialVals[1], initialVals[2]);
            if (!candidate.ParametersValid || !Tools.IsFinite(candidate.LogLikelihood(sample)))
                initialVals = [moments[0] * magnitude, moments[1] * magnitude, 0];
            // Get bounds of location
            double locationMagnitude = Math.Max(Math.Abs(initialVals[0]), initialVals[1]);
            lowerVals[0] = -FiniteDecimalBound(locationMagnitude);
            upperVals[0] = -lowerVals[0];
            // Get bounds of scale
            lowerVals[1] = Math.Min(Tools.DoubleMachineEpsilon, initialVals[1] / 10);
            upperVals[1] = FiniteDecimalBound(initialVals[1]);
            // Get bounds of shape
            lowerVals[2] = -10;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
            candidate.SetParameters(initialVals);
            if (!candidate.ParametersValid || !Tools.IsFinite(candidate.LogLikelihood(sample))
                || initialVals[0] <= lowerVals[0] || initialVals[0] >= upperVals[0]
                || initialVals[1] <= lowerVals[1] || initialVals[1] >= upperVals[1])
                throw new InvalidOperationException("The sample does not admit a finite supported generalized-logistic initializer within finite bounds.");
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
            if (solver.Status != OptimizationStatus.Success
                || ValidateParameters(solver.BestParameterSet.Values, false) != null
                || !Tools.IsFinite(new GeneralizedLogistic(solver.BestParameterSet.Values[0], solver.BestParameterSet.Values[1], solver.BestParameterSet.Values[2]).LogLikelihood(sample)))
                throw new InvalidOperationException($"Generalized logistic maximum likelihood estimation failed with optimizer status {solver.Status} or a nonfinite fit.");
            return solver.BestParameterSet.Values;
        }

        /// <summary>Retains decimal-order fitting bounds without overflow or a zero-centered collapse.</summary>
        /// <param name="value">The positive magnitude from which to select the next decimal order.</param>
        /// <returns>The next decimal-order bound, capped at the largest finite binary64 value.</returns>
        private static double FiniteDecimalBound(double value)
        {
            double bound = Math.Pow(10, Math.Ceiling(Math.Log10(value)) + 1);
            return double.IsPositiveInfinity(bound) ? double.MaxValue : bound;
        }

        /// <summary>Retains a finite support endpoint when an intermediate scale/shape quotient overflows.</summary>
        /// <returns>The finite-shape support endpoint in physical coordinates.</returns>
        private double FiniteShapeEndpoint()
        {
            double shift = Alpha / Kappa;
            return double.IsInfinity(shift) && Math.Sign(Xi) != Math.Sign(shift)
                ? (Xi * Kappa + Alpha) / Kappa : Xi + shift;
        }

        /// <summary>Inverts the exact Hosking transformation, including compensated endpoint residuals.</summary>
        /// <param name="x">An observation inside the distribution support.</param>
        /// <returns>The corresponding standard logistic variate.</returns>
        private double LatentLogistic(double x)
        {
            double y = DistributionNumerics.Standardize(x, Xi, Alpha);
            if (Kappa == 0) return y;
            double product = -Kappa * y;
            if (product < -.9) return -KappaFourBoundary.LogT(x, Xi, Alpha, Kappa);
            return DistributionNumerics.HoskingShapeTransform(x, Xi, Alpha, Kappa);
        }

        /// <inheritdoc/>
        public override double PDF(double x) => Math.Exp(LogPDF(x));

        /// <inheritdoc/>
        /// <remarks>Evaluates log density directly and preserves the one-sided infinite density at singular support endpoints.</remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsInfinity(x) || x < Minimum || x > Maximum) return double.NegativeInfinity;
            if (x == Minimum || x == Maximum)
                return Math.Abs(Kappa) < 1 ? double.NegativeInfinity : Math.Abs(Kappa) == 1 ? -Math.Log(Alpha) : double.PositiveInfinity;
            double z = LatentLogistic(x);
            return (z >= 0 ? (Kappa - 1) * z - 2 * Tools.Log1p(Math.Exp(-z))
                : (Kappa + 1) * z - 2 * Tools.Log1p(Math.Exp(z))) - Math.Log(Alpha);
        }

        /// <inheritdoc/>
        public override double CDF(double x) => Math.Exp(LogCDF(x));

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0;
            double z = LatentLogistic(x);
            return z <= 0 ? z - Tools.Log1p(Math.Exp(z)) : -Tools.Log1p(Math.Exp(-z));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return 0;
            if (x >= Maximum) return double.NegativeInfinity;
            double z = LatentLogistic(x);
            return z >= 0 ? -z - Tools.Log1p(Math.Exp(-z)) : -Tools.Log1p(Math.Exp(z));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            if (!(probability >= 0 && probability <= 1))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between zero and one.");
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (probability == 0) return Minimum;
            if (probability == 1) return Maximum;
            double z = Math.Log(probability) - Tools.Log1p(-probability);
            return QuantileAtLatent(z);
        }

        /// <summary>Combines shape exponentials and physical scale before exponentiation or affine addition.</summary>
        /// <param name="z">The standard logistic quantile.</param>
        /// <returns>The corresponding quantile in physical coordinates.</returns>
        private double QuantileAtLatent(double z)
        {
            double v = -Kappa * z;
            double standard = v > 50 ? -Math.Sign(Kappa) * Math.Exp(v - Math.Log(Math.Abs(Kappa)))
                : v < -50 ? 1 / Kappa : z * DistributionNumerics.Exprel(v);
            double offset = v > 50 ? -Math.Sign(Kappa) * Math.Exp(Math.Log(Alpha) + v - Math.Log(Math.Abs(Kappa)))
                : Alpha * standard;
            double value = Xi + offset;
            if (double.IsInfinity(value) && Tools.IsFinite(standard))
            {
                double combined = Xi / Alpha + standard;
                if (Tools.IsFinite(combined)) return Alpha * combined;
            }
            return value;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone() => new GeneralizedLogistic(Xi, Alpha, Kappa);

        /// <inheritdoc/>
        /// <remarks>
        /// Local asymptotic MLE covariance in xi, alpha, kappa order. Regular information
        /// requires absolute kappa below one half; this condition does not restrict distribution
        /// validity or assert existence of a finite global MLE.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Sample size, scale or information regularity is invalid.</exception>
        /// <exception cref="InvalidOperationException">Numerical information or its inversion cannot be resolved.</exception>
        /// <exception cref="NotImplementedException">The requested estimator is not maximum likelihood.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException("Generalized-logistic covariance is implemented only for local maximum-likelihood uncertainty.");
            return KappaExpectedInformation.ParameterCovariance(Alpha, Kappa, -1, sampleSize, 3);
        }

        /// <inheritdoc/>
        /// <remarks>Applies the local-MLE delta method in common physical quantile coordinates,
        /// avoiding underflow or overflow from forming physical covariance entries first.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            var unit = new GeneralizedLogistic(0, 1, Kappa);
            double[,] covariance = unit.ParameterCovariance(sampleSize, estimationMethod);
            double z = Math.Log(probability) - Tools.Log1p(-probability), argument = -Kappa * z;
            double scaleGradient = DistributionNumerics.ScaledExprelProduct(Alpha, z, argument);
            double shapeGradient = -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, z, argument);
            return DistributionNumerics.ScaledQuantileVariance(covariance, [Alpha, scaleGradient, shapeGradient]);
        }

        /// <inheritdoc/>
        /// <remarks>The analytic kappa derivative is continuous at zero and equals -alpha*logit(p)^2/2 there.</remarks>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            double z = Math.Log(probability) - Tools.Log1p(-probability), v = -Kappa * z;
            double shape = v < -50 ? -Math.Exp(Math.Log(Alpha) - 2 * Math.Log(Math.Abs(Kappa)))
                : v > 50 ? -Math.Exp(Math.Log(Alpha) + v + Math.Log(v - 1) - 2 * Math.Log(Math.Abs(Kappa)))
                : -Alpha * (z * z * DistributionNumerics.ExprelDerivative(v));
            double scale = v > 50 ? -Math.Sign(Kappa) * Math.Exp(v - Math.Log(Math.Abs(Kappa)))
                : v < -50 ? 1 / Kappa : z * DistributionNumerics.Exprel(v);
            return [1, scale, shape];
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
            => DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
    }
}
