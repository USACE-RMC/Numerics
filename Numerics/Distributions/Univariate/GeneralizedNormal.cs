using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;

namespace Numerics.Distributions
{

    /// <summary>
    /// The generalized normal distribution (LogNormal-3).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public class GeneralizedNormal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {

        /// <summary>
        /// Constructs a Generalized Normal distribution with a location of 100, scale of 10, and shape of 0.
        /// </summary>
        public GeneralizedNormal()
        {
            SetParameters(100d, 10d, 0d);
        }

        /// <summary>
        /// Constructs a Generalized Normal distribution with the given parameters ξ, α, and κ.
        /// </summary>
        /// <param name="location">The location parameter ξ (Xi).</param>
        /// <param name="scale">The scale parameter α (alpha).</param>
        /// <param name="shape">The shape parameter κ (kappa).</param>
        public GeneralizedNormal(double location, double scale, double shape)
        {
            SetParameters(location, scale, shape);
        }

        private double _xi; // location
        private double _alpha; // scale
        private double _kappa; // shape
        private bool _momentsComputed = false;
        private double[] u = [ double.NaN, double.NaN, double.NaN, double.NaN ];

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
                _momentsComputed = false;
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
                _momentsComputed = false;
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
                _momentsComputed = false;
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
            get { return UnivariateDistributionType.GeneralizedNormal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Generalized Normal"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "GNO"; }
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
            get { return [Xi, Alpha, Kappa ]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = AnalyticalMoments();
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

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
                return QuantileAtLatent(Kappa);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get
            {
                if (!_momentsComputed)
                {
                    u = AnalyticalMoments();
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
                    u = AnalyticalMoments();
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
                    u = AnalyticalMoments();
                    _momentsComputed = true;
                }
                return u[3];
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
            if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
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
            var newDistribution = new GeneralizedNormal(Xi, Alpha, Kappa);
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
            Alpha = scale;
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

        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            if (moments == null) throw new ArgumentNullException(nameof(moments));
            if (moments.Count < 3 || !DistributionNumerics.IsFinite(moments[0]) || !DistributionNumerics.IsFinite(moments[1])
                || moments[1] <= 0 || !DistributionNumerics.IsFinite(moments[2]) || Math.Abs(moments[2]) >= 1)
                throw new ArgumentOutOfRangeException(nameof(moments), "Finite L-moments require positive L-scale and absolute L-skewness below one.");
            double L1 = moments[0];
            double L2 = moments[1];
            double T3 = moments[2];

            double E0 = 2.0466534;
            double E1 = -3.6544371;
            double E2 = 1.8396733;
            double E3 = -0.20360244;
            double F1 = -2.0182173;
            double F2 = 1.2420401;
            double F3 = -0.21741801;

            double kappa = -T3 * (E0 + E1 * Math.Pow(T3, 2d) + E2 * Math.Pow(T3, 4d) + E3 * Math.Pow(T3, 6d)) / (1d + F1 * Math.Pow(T3, 2d) + F2 * Math.Pow(T3, 4d) + F3 * Math.Pow(T3, 6d));
            double alpha = L2 / NormalLScale(kappa);
            double xi = L1 - new GeneralizedNormal(0, alpha, kappa).Mean;
            return [xi, alpha, kappa];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double xi = parameters[0];
            double alpha = parameters[1];
            double kappa = parameters[2];
            ValidateParameters(xi, alpha, kappa, true);

            double A0 = 4.8860251 * Math.Pow(10, -1);
            double A1 = 4.4493076 * Math.Pow(10, -3); 
            double A2 = 8.8027039 * Math.Pow(10, -4); 
            double A3 = 1.1507084 * Math.Pow(10, -6); 
            double B1 = 6.4662924 * Math.Pow(10, -2); 
            double B2 = 3.3090406 * Math.Pow(10, -3); 
            double B3 = 7.4290680 * Math.Pow(10, -5); 
            double C0 = 1.8756590 * Math.Pow(10, -1); 
            double C1 = -2.5352147 * Math.Pow(10, -3); 
            double C2 = 2.6995102 * Math.Pow(10, -4); 
            double C3 = -1.8446680 * Math.Pow(10, -6); 
            double D1 = 8.2325617 * Math.Pow(10, -2); 
            double D2 = 4.2681448 * Math.Pow(10, -3); 
            double D3 = 1.1653690 * Math.Pow(10, -4); 
            double tau40 = 0.12260171954089095; // 30*asin(1/3)/pi - 9, the exact normal L-kurtosis limit.

            double L1 = new GeneralizedNormal(xi, alpha, kappa).Mean;
            double L2 = alpha * NormalLScale(kappa);
            double T3 = -kappa * (A0 + A1 * Math.Pow(kappa, 2d) + A2 * Math.Pow(kappa, 4d) + A3 * Math.Pow(kappa, 6d)) / (1d + B1 * Math.Pow(kappa, 2d) + B2 * Math.Pow(kappa, 4d) + B3 * Math.Pow(kappa, 6d));
            double T4 = tau40 + Math.Pow(kappa, 2d) * (C0 + C1 * Math.Pow(kappa, 2d) + C2 * Math.Pow(kappa, 4d) + C3 * Math.Pow(kappa, 6d)) / (1d + D1 * Math.Pow(kappa, 2d) + D2 * Math.Pow(kappa, 4d) + D3 * Math.Pow(kappa, 6d));
            return [L1, L2, T3, T4];
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
            initialVals = LegacyConstraintParametersFromLinearMoments(Statistics.LinearMoments(sample));
            // Get bounds of location
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of scale
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[1]))) + 1d);
            // Get bounds of shape
            lowerVals[2] = -10d;
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

            double E0 = 2.0466534;
            double E1 = -3.6544371;
            double E2 = 1.8396733;
            double E3 = -0.20360244;
            double F1 = -2.0182173;
            double F2 = 1.2420401;
            double F3 = -0.21741801;

            double kappa = -T3 * (E0 + E1 * Math.Pow(T3, 2d) + E2 * Math.Pow(T3, 4d) + E3 * Math.Pow(T3, 6d)) / (1d + F1 * Math.Pow(T3, 2d) + F2 * Math.Pow(T3, 4d) + F3 * Math.Pow(T3, 6d));
            double alpha = (L2 * kappa * Math.Exp(-(kappa * kappa) / 2d)) / (1d - 2 * Normal.StandardCDF(-kappa / Tools.Sqrt2));
            double xi = L1 - alpha * (1.0d - Math.Exp(kappa * kappa / 2d)) / kappa;
            return [xi, alpha, kappa];
        }

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
        {
            DistributionNumerics.ValidateSample(sample, 4);
            // Estimate initial values using the method of moments (a.k.a product moments).
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Get initial values
            double magnitude = 0;
            foreach (double value in sample) magnitude = Math.Max(magnitude, Math.Abs(value));
            var scaled = new double[sample.Count];
            for (int i = 0; i < sample.Count; i++) scaled[i] = sample[i] / magnitude;
            double[] moments = Statistics.LinearMoments(scaled);
            initialVals = ParametersFromLinearMoments(moments);
            initialVals[0] *= magnitude;
            initialVals[1] *= magnitude;
            var candidate = new GeneralizedNormal(initialVals[0], initialVals[1], initialVals[2]);
            if (!candidate.ParametersValid || !DistributionNumerics.IsFinite(candidate.LogLikelihood(sample)))
                initialVals = [moments[0] * magnitude, moments[1] * Math.Sqrt(Math.PI) * magnitude, 0];
            // Get bounds of location
            double locationMagnitude = Math.Max(Math.Abs(initialVals[0]), initialVals[1]);
            lowerVals[0] = -FiniteDecimalBound(locationMagnitude);
            upperVals[0] = -lowerVals[0];
            // Get bounds of scale
            lowerVals[1] = Math.Min(Tools.DoubleMachineEpsilon, initialVals[1] / 10);
            upperVals[1] = FiniteDecimalBound(initialVals[1]);
            // Get bounds of shape
            lowerVals[2] = -10d;
            upperVals[2] = 10d;
            // Correct initial value of kappa if necessary
            if (initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
            {
                initialVals[2] = 0d;
            }
            candidate.SetParameters(initialVals);
            if (!candidate.ParametersValid || !DistributionNumerics.IsFinite(candidate.LogLikelihood(sample))
                || initialVals[0] <= lowerVals[0] || initialVals[0] >= upperVals[0]
                || initialVals[1] <= lowerVals[1] || initialVals[1] >= upperVals[1])
                throw new InvalidOperationException("The sample does not admit a finite supported generalized-normal initializer within finite bounds.");
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
                var GLO = new GeneralizedNormal();
                GLO.SetParameters(x);
                return GLO.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            if (solver.Status != OptimizationStatus.Success
                || ValidateParameters(solver.BestParameterSet.Values, false) != null
                || !DistributionNumerics.IsFinite(new GeneralizedNormal(solver.BestParameterSet.Values[0], solver.BestParameterSet.Values[1], solver.BestParameterSet.Values[2]).LogLikelihood(sample)))
                throw new InvalidOperationException($"Generalized normal maximum likelihood estimation failed with optimizer status {solver.Status} or a nonfinite fit.");
            return solver.BestParameterSet.Values;
        }

        /// <summary>Analytical shifted-lognormal moments, preserving representable scale products.</summary>
        private double[] AnalyticalMoments()
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (Kappa == 0) return [Xi, Alpha, 0, 3];
            double v = Kappa * Kappa;
            double t = Tools.Expm1(v);
            double mean, sd, skew;
            if (v <= .5)
            {
                double relative = DistributionNumerics.Exprel(v);
                mean = Xi - Alpha * Kappa * .5 * DistributionNumerics.Exprel(.5 * v);
                sd = Alpha * Math.Exp(.5 * v) * Math.Sqrt(relative);
                skew = -Kappa * (t + 3) * Math.Sqrt(relative);
            }
            else
            {
                double logt = v > 36 ? v + Tools.Log1p(-Math.Exp(-v)) : Math.Log(t);
                double loghalf = v > 72 ? v / 2 + Tools.Log1p(-Math.Exp(-v / 2)) : Math.Log(Tools.Expm1(v / 2));
                mean = Xi - Math.Sign(Kappa) * Math.Exp(Math.Log(Alpha) + loghalf - Math.Log(Math.Abs(Kappa)));
                sd = Math.Exp(Math.Log(Alpha) + v / 2 + logt / 2 - Math.Log(Math.Abs(Kappa)));
                double logsum = v > 36 ? v + Tools.Log1p(2 * Math.Exp(-v)) : Math.Log(t + 3);
                skew = -Math.Sign(Kappa) * Math.Exp(logsum + logt / 2);
            }
            return [mean, sd, skew, 3 + t * (16 + t * (15 + t * (6 + t)))];
        }

        /// <summary>Returns unit-scale L-scale with the exact kappa=0 normal limit.</summary>
        private static double NormalLScale(double k)
        {
            double v = k * k;
            if (Math.Abs(k) < .5)
            {
                double sum = 1, power = 1;
                for (int n = 1; n < 24; n++)
                {
                    power *= -v / (4 * n);
                    double term = power / (2 * n + 1);
                    sum += term;
                    if (Math.Abs(term) < Math.Abs(sum) * 1E-17) break;
                }
                return Math.Exp(v / 2) * sum / Math.Sqrt(Math.PI);
            }
            return Math.Exp(v / 2) * (1 - 2 * Normal.StandardCDF(-Math.Abs(k) / Tools.Sqrt2)) / Math.Abs(k);
        }

        /// <summary>Retains the existing decimal-order bound construction while preventing overflow.</summary>
        private static double FiniteDecimalBound(double value)
        {
            double bound = Math.Pow(10, Math.Ceiling(Math.Log10(value)) + 1);
            return double.IsPositiveInfinity(bound) ? double.MaxValue : bound;
        }

        /// <summary>Retains a finite support endpoint when an intermediate scale/shape quotient overflows.</summary>
        private double FiniteShapeEndpoint()
        {
            double shift = Alpha / Kappa;
            return double.IsInfinity(shift) && Math.Sign(Xi) != Math.Sign(shift)
                ? (Xi * Kappa + Alpha) / Kappa : Xi + shift;
        }

        /// <summary>Inverts the Hosking shape transform, retaining its exact nonzero shape and support residual.</summary>
        private double LatentNormal(double x)
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
        /// <remarks>Evaluates the transformed-normal log density directly, including underflowing densities.</remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (double.IsNaN(x)) return double.NaN;
            if (double.IsInfinity(x) || x <= Minimum || x >= Maximum) return double.NegativeInfinity;
            double z = LatentNormal(x);
            return -Math.Log(Alpha) - Tools.LogSqrt2PI + z * (Kappa - z / 2);
        }

        /// <inheritdoc/>
        public override double CDF(double x) => Math.Exp(LogCDF(x));

        /// <inheritdoc/>
        public override double LogCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return double.NegativeInfinity;
            if (x >= Maximum) return 0;
            return DistributionNumerics.NormalLogCDF(LatentNormal(x));
        }

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return 0;
            if (x >= Maximum) return double.NegativeInfinity;
            return DistributionNumerics.NormalLogSurvival(LatentNormal(x));
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            if (!(probability >= 0 && probability <= 1))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between zero and one.");
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (probability == 0) return Minimum;
            if (probability == 1) return Maximum;
            return QuantileAtLatent(Normal.StandardZ(probability));
        }

        /// <summary>Combines shape exponentials and physical scale before exponentiation or affine addition.</summary>
        private double QuantileAtLatent(double z)
        {
            double v = -Kappa * z;
            double standard = v > 50 ? -Math.Sign(Kappa) * Math.Exp(v - Math.Log(Math.Abs(Kappa)))
                : v < -50 ? 1 / Kappa : z * DistributionNumerics.Exprel(v);
            double offset = v > 50 ? -Math.Sign(Kappa) * Math.Exp(Math.Log(Alpha) + v - Math.Log(Math.Abs(Kappa)))
                : Alpha * standard;
            double value = Xi + offset;
            if (double.IsInfinity(value) && DistributionNumerics.IsFinite(standard))
            {
                double combined = Xi / Alpha + standard;
                if (DistributionNumerics.IsFinite(combined)) return Alpha * combined;
            }
            return value;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone() => new GeneralizedNormal(Xi, Alpha, Kappa);
        /// <inheritdoc/>
        /// <remarks>
        /// Local asymptotic MLE covariance for all three estimated parameters in xi, alpha,
        /// kappa order. Information is finite for every finite kappa. This does not assert
        /// existence of a finite global MLE for the three-parameter lognormal likelihood.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">The sample size or distribution parameters are invalid.</exception>
        /// <exception cref="NotImplementedException">The requested estimator is not maximum likelihood.</exception>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException("Generalized-normal covariance is implemented only for local maximum-likelihood uncertainty.");
            NormalInformationFactors(Kappa, out double c, out double inverseR, out double logInverseR);
            double scale = Alpha / Math.Sqrt(sampleSize), shapedScale = Alpha * Kappa / Math.Sqrt(sampleSize);
            double logAlpha = Math.Log(Alpha), logN = Math.Log(sampleSize), logK = Math.Log(Math.Abs(Kappa));
            var covariance = new double[3, 3];
            covariance[0, 0] = scale * scale * (1 + c * c * inverseR);
            covariance[1, 1] = shapedScale * shapedScale + scale * scale / 2;
            covariance[2, 2] = (Kappa / 2 / sampleSize) * Kappa + inverseR / sampleSize;
            covariance[0, 1] = -scale * shapedScale;
            if (!DistributionNumerics.IsFinite(covariance[0, 1]) || covariance[0, 1] == 0)
                covariance[0, 1] = -Math.Sign(Kappa) * Math.Exp(2 * logAlpha + logK - logN);
            covariance[0, 2] = c * inverseR * Alpha / sampleSize;
            if (c > 0 && (covariance[0, 2] == 0 || !DistributionNumerics.IsFinite(covariance[0, 2])))
                covariance[0, 2] = Math.Exp(logAlpha + Math.Log(c) + logInverseR - logN);
            covariance[1, 2] = Alpha * Kappa / 2 / sampleSize;
            if (!DistributionNumerics.IsFinite(covariance[1, 2]) || covariance[1, 2] == 0)
                covariance[1, 2] = Math.Sign(Kappa) * Math.Exp(logAlpha + logK - Math.Log(2) - logN);
            covariance[1, 0] = covariance[0, 1];
            covariance[2, 0] = covariance[0, 2];
            covariance[2, 1] = covariance[1, 2];
            return covariance;
        }

        /// <summary>Evaluates the closed-form information factors, retaining log(1/R) after underflow.</summary>
        private static void NormalInformationFactors(double kappa, out double c, out double inverseR, out double logInverseR)
        {
            double v = kappa * kappa;
            c = .5 * DistributionNumerics.Exprel(-v / 2);
            if (v < .1)
            {
                double sum = 1.5, term = 1.5;
                for (int m = 2; m < 24; m++)
                {
                    term *= v * (m + 2d) / (m + 1) / (m + 1);
                    sum += term;
                    if (Math.Abs(term) < Math.Abs(sum) * 1E-17) break;
                }
                inverseR = 1 / sum;
                logInverseR = -Math.Log(sum);
            }
            else
            {
                logInverseR = double.IsPositiveInfinity(v) ? double.NegativeInfinity : 2 * Math.Log(v) - v
                    - (v > 350 ? Tools.Log1p(v) : Math.Log(1 + v - (1 + 2 * v) * Math.Exp(-v)));
                inverseR = Math.Exp(logInverseR);
            }
            if (double.IsPositiveInfinity(v)) c = 0;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Applies the local-MLE delta method in public parameter coordinates. An exact positive
        /// factorization of the covariance avoids cancellation near a finite shape endpoint.
        /// Scale and exponential factors are combined before the scalar variance is exponentiated.
        /// </remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            DistributionNumerics.ValidateSampleSize(sampleSize);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
                throw new NotImplementedException("Generalized-normal covariance is implemented only for local maximum-likelihood uncertainty.");
            NormalInformationFactors(Kappa, out double c, out _, out double logInverseR);
            double z = Normal.StandardZ(probability), argument = -Kappa * z;
            double logScale = Math.Log(Alpha) - .5 * Math.Log(sampleSize);
            // C = a*a' + b*b'/2 + d*d'/R, with a=[1,-k,0], b=[0,1,k], d=[c,0,1].
            // For g=[1,S,T], g'a=exp(-k*z) and g'b=z*exp(-k*z).
            double first = 2 * (logScale + argument) + Tools.Log1p(z * z / 2);
            if (double.IsPositiveInfinity(first)) return double.PositiveInfinity;
            double logResidual;
            if (Math.Abs(Kappa) < .1)
            {
                double residual = c - z * z * DistributionNumerics.ExprelDerivative(argument);
                logResidual = Math.Log(Math.Abs(residual));
            }
            else
            {
                // c+T = ((1-argument)*exp(argument)-exp(-k*k/2))/(k*k).
                // Retain this difference even when c and T individually round to +/-1/(k*k).
                double negative = -(Kappa / 2) * Kappa;
                double numerator;
                if (argument > 1)
                    numerator = DistributionNumerics.LogSum(argument + Math.Log(argument - 1), negative);
                else
                {
                    double positive = argument == 1 || double.IsNegativeInfinity(argument)
                        ? double.NegativeInfinity : argument + Tools.Log1p(-argument);
                    numerator = DistributionNumerics.LogDifference(Math.Max(positive, negative), Math.Min(positive, negative));
                }
                logResidual = numerator - 2 * Math.Log(Math.Abs(Kappa));
            }
            double second = 2 * (logScale + logResidual) + logInverseR;
            return Math.Exp(DistributionNumerics.LogSum(first, second));
        }

        /// <inheritdoc/>
        /// <remarks>The analytical shape derivative is continuous at kappa=0 and equals -alpha*z*z/2 there.</remarks>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            double z = Normal.StandardZ(probability), v = -Kappa * z;
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
