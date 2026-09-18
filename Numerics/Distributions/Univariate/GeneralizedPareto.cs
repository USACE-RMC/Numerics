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
        /// <remarks>The mean exists for kappa &gt; -1; bounded positive shapes are not excluded.</remarks>
        public override double Mean
        {
            get { return Kappa > -1 ? Xi + Alpha / (1 + Kappa) : double.NaN; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(.5); }
        }

        /// <inheritdoc/>
        /// <remarks>The lower endpoint is the mode for kappa &lt; 1, the uniform kappa=1 case has no
        /// unique mode and returns NaN, and kappa &gt; 1 has its mode at the upper endpoint.</remarks>
        public override double Mode
        {
            get { return Kappa < 1 ? Xi : Kappa == 1 ? double.NaN : Maximum; }
        }

        /// <inheritdoc/>
        /// <remarks>The variance exists for kappa &gt; -1/2; scale is applied after standardization.</remarks>
        public override double StandardDeviation
        {
            get { return Kappa > -.5 ? (Alpha / (1 + Kappa)) / Math.Sqrt(1 + 2 * Kappa) : double.NaN; }
        }

        /// <inheritdoc/>
        /// <remarks>The third moment exists for kappa &gt; -1/3.</remarks>
        public override double Skewness
        {
            get
            {
                if (Kappa > -1d / 3d)
                {
                    if (Kappa > 1)
                    {
                        double inverse = 1 / Kappa;
                        return (2 * (inverse - 1) / (inverse + 3)) * Math.Sqrt(Kappa) * Math.Sqrt(inverse + 2);
                    }
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
        /// <remarks>Ordinary kurtosis exists for kappa &gt; -1/4.</remarks>
        public override double Kurtosis
        {
            get
            {
                if (Kappa > -.25)
                {
                    if (Kappa > 1)
                    {
                        double inverse = 1 / Kappa;
                        double coefficient = 3 * (inverse + 2) * (3 * inverse * inverse - inverse + 2)
                            / ((inverse + 3) * (inverse + 4));
                        return Kappa * coefficient;
                    }
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
            _xi = location;
            _alpha = scale;
            _kappa = shape;
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
            ValidateParameters(parameters, true);
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
            ValidateParameters(parameters, true);
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
        /// <remarks>Requires at least four finite, nonconstant observations. Preserves the legacy linear-moment
        /// initialization and family-specific bounds whenever they are finite, ordered, and contain the initial
        /// values; the legacy upper location bound is the sample minimum plus <see cref="Tools.DoubleMachineEpsilon"/>.
        /// Otherwise, the exceptional-input fallback evaluates the linear-moment initialization in normalized
        /// coordinates and rescales it algebraically; its upper location bound is the sample minimum.</remarks>
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

        /// <summary>Handles samples whose legacy initialization or bounds are not finite or outside ordered bounds.</summary>
        /// <param name="sample">The validated observations.</param>
        /// <returns>Finite initial values and bounds from the hardened initialization path.</returns>
        private Tuple<double[], double[], double[]> GetRobustParameterConstraints(IList<double> sample)
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
            double minData = Statistics.Minimum(sample);
            DistributionNumerics.LocationParameterBounds(ref initialVals[0], initialVals[1], minData,
                Statistics.Maximum(sample), true, out lowerVals[0], out upperVals[0]);
            DistributionNumerics.PositiveParameterBounds(initialVals[1], out lowerVals[1], out upperVals[1]);
            // Get bounds of shape
            lowerVals[2] = -10d;
            upperVals[2] = 10d;
            // Correct initial values if necessary
            if (!Tools.IsFinite(initialVals[2]) || initialVals[2] <= lowerVals[2] || initialVals[2] >= upperVals[2])
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
            return Math.Exp(LogPDF(x));
        }

        /// <inheritdoc/>
        /// <remarks>Uses the exact nonzero shape and retains one-sided finite-endpoint density limits.</remarks>
        public override double LogPDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x < Minimum || x > Maximum || double.IsInfinity(x)) return double.NegativeInfinity;
            if (Kappa > 0 && x == Maximum)
                return Kappa < 1 ? double.NegativeInfinity : Kappa == 1 ? -Math.Log(Alpha) : double.PositiveInfinity;
            double value = -(1 - Kappa) * TransformedValue(x) - Math.Log(Alpha);
            return double.IsNaN(value) ? double.NegativeInfinity : value;
        }

        /// <summary>Maps an interior observation to its exponential coordinate without a shape-zero plateau.</summary>
        /// <param name="x">The observation in physical coordinates.</param>
        /// <returns>The corresponding unit-exponential coordinate.</returns>
        private double TransformedValue(double x)
        {
            return DistributionNumerics.HoskingShapeTransform(x, Xi, Alpha, Kappa);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            return -Tools.Expm1(LogCCDF(x));
        }

        /// <inheritdoc/>
        public override double LogCDF(double x) => DistributionNumerics.Log1mExp(LogCCDF(x));

        /// <inheritdoc/>
        public override double CCDF(double x) => Math.Exp(LogCCDF(x));

        /// <inheritdoc/>
        public override double LogCCDF(double x)
        {
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            if (x <= Minimum) return 0;
            if (x >= Maximum) return double.NegativeInfinity;
            return -TransformedValue(x);
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
            double logarithm = Tools.Log1p(-probability);
            double product = Kappa * logarithm;
            double unitQuantile = double.IsNegativeInfinity(product) ? 1 / Kappa : DistributionNumerics.ScaledExprelProduct(1, -logarithm, product);
            double displacement = double.IsNegativeInfinity(product) ? Alpha / Kappa : DistributionNumerics.ScaledExprelProduct(Alpha, -logarithm, product);
            return double.IsInfinity(displacement) && Tools.IsFinite(unitQuantile)
                ? Alpha * (Xi / Alpha + unitQuantile) : Xi + displacement;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new GeneralizedPareto(Xi, Alpha, Kappa) { Lambda = Lambda };
        }

        /// <inheritdoc/>
        /// <remarks>MLE uncertainty requires kappa &lt; 1/2; MoM requires kappa &gt; -1/4. Location
        /// covariance additionally requires sampleSize + 2*kappa &gt; 0. These restrictions do not narrow distribution validity or fitting bounds.</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters, sample size or the requested uncertainty domain are invalid.</exception>
        /// <exception cref="NotImplementedException">The requested estimator combination is unsupported.</exception>
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
                ValidateParameters(Xi, _alpha, Kappa, true);
            if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood && Kappa >= .5)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "Regular maximum-likelihood uncertainty requires kappa < 1/2.");
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments && Kappa <= -.25)
                throw new ArgumentOutOfRangeException(nameof(Kappa), "Method-of-moments uncertainty requires kappa > -1/4.");
            double a = Alpha;
            double k = Kappa;
            double N = sampleSize;
            if (!(N + 2 * k > 0))
                throw new ArgumentOutOfRangeException(nameof(sampleSize), "Location covariance requires sampleSize + 2*kappa > 0.");
            var covar = new double[3, 3];
            double locationScale = a / (N + k);
            covar[0, 0] = (N / (N + 2 * k)) * (locationScale * locationScale); // location
            double scaleVariance = a / Math.Sqrt(N);
            scaleVariance *= scaleVariance;
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                double num = Math.Pow(1d + k, 2d) * (1d + 6d * k + 12d * Math.Pow(k, 2d));
                double den = (1d + 2d * k) * (1d + 3d * k) * (1d + 4d * k);
                covar[1, 1] = 2d * scaleVariance * (num / den); // scale
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
                covar[1, 1] = (1d - k) * (2d * scaleVariance); // scale
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
        /// <remarks>Preserves the documented omission of order-n-to-the-minus-two location uncertainty.
        /// The scale/shape covariance is contracted in finite physical common coordinates. Probability must be finite and strictly interior.</remarks>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (!_parametersValid) ValidateParameters(Xi, Alpha, Kappa, true);
            // The variance in the location parameter is of order N-2 and thus
            // can be neglected relative to the variances of scale and shape (order N-1)
            // to obtain approximate confidence intervals. 
            var unit = new GeneralizedPareto(0, 1, Kappa);
            var covar = unit.ParameterCovariance(sampleSize, estimationMethod);
            double logarithm = Tools.Log1p(-probability), product = Kappa * logarithm;
            double scaleGradient = DistributionNumerics.ScaledExprelProduct(Alpha, -logarithm, product);
            double shapeGradient = -DistributionNumerics.ScaledExprelDerivativeProduct(Alpha, logarithm, product);
            var scaleShapeCovariance = new[,] { { covar[1, 1], covar[1, 2] }, { covar[2, 1], covar[2, 2] } };
            return DistributionNumerics.ScaledQuantileVariance(scaleShapeCovariance, [scaleGradient, shapeGradient]);
        }

        /// <inheritdoc/>
        /// <remarks>Uses analytical exponential divided differences, including the full nonzero
        /// shape sensitivity at kappa=0, with L=log(1-p).</remarks>
        /// <exception cref="ArgumentOutOfRangeException">Parameters are invalid or probability is not finite and strictly interior.</exception>
        public double[] QuantileGradient(double probability)
        {
            DistributionNumerics.ValidateProbability(probability);
            if (_parametersValid == false)
                ValidateParameters(Xi, _alpha, Kappa, true);
            double logarithm = Tools.Log1p(-probability);
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
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            return DistributionNumerics.QuantileJacobian(this, probabilities, out determinant);
        }

    }
}
