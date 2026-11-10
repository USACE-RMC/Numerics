using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.RootFinding;
using System.Globalization;
using System.Threading;
using System.Xml.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// A general truncated probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// The base-distribution probabilities at the truncation bounds are cached together with a bitwise snapshot of
    /// the base parameters. When the base distribution is changed through <see cref="BaseDistribution"/>, the next
    /// evaluation detects the change, recomputes the cached probabilities and moments, and refreshes
    /// <see cref="UnivariateDistributionBase.ParametersValid"/>.
    /// </para>
    /// <para>
    /// When the base CDF at the lower bound exceeds one half, the truncated mass, distribution, survival and quantile
    /// functions use differences of the base survival function (CCDF), which keep their relative accuracy in the
    /// upper tail where differences of CDF values near one cancel.
    /// </para>
    /// </remarks>
    [Serializable]
    public class TruncatedDistribution : UnivariateDistributionBase
    {

        /// <summary>
        /// Constructs a truncated distribution.
        /// </summary>
        /// <param name="basDistribution">The base distribution to truncate.</param>
        /// <param name="min">The minimum possible value of the distribution.</param>
        /// <param name="max">The maximum possible value of the distribution.</param>
        public TruncatedDistribution(UnivariateDistributionBase basDistribution, double min, double max)
        {
            _baseDist = basDistribution;
            _min = min;
            _max = max;
            var parameters = _baseDist.GetParameters.ToList();
            parameters.Add(min);
            parameters.Add(max);
            var state = CreateState(parameters);
            _parametersValid = state.Valid;
            _state = state;
        }

        private readonly UnivariateDistributionBase _baseDist;
        private double _min;
        private double _max;
        [NonSerialized] private TruncationState? _state;

        /// <summary>
        /// The cached truncation probabilities and moments for one bitwise state of the base parameters.
        /// </summary>
        /// <remarks>
        /// Instances are immutable apart from the lazily published moments, so concurrent readers share them safely.
        /// </remarks>
        private sealed class TruncationState
        {
            /// <summary>The allocation-free snapshot of the base distribution, or null when its type is not capturable.</summary>
            private readonly DistributionSnapshot? _snapshot;

            /// <summary>The base parameters, compared bitwise when no snapshot is available.</summary>
            private readonly double[] _baseParameters;

            /// <summary>The moments {mean, standard deviation, skewness, kurtosis}, or null until computed.</summary>
            private double[]? _moments;

            /// <summary>The base CDF at the lower bound.</summary>
            internal readonly double LowerProbability;

            /// <summary>The base survival probability at the upper bound, used when <see cref="UpperTail"/> is true.</summary>
            internal readonly double UpperSurvival;

            /// <summary>The base survival probability at the lower bound, used when <see cref="UpperTail"/> is true.</summary>
            internal readonly double LowerSurvival;

            /// <summary>True when the base CDF at the lower bound exceeds one half.</summary>
            internal readonly bool UpperTail;

            /// <summary>The truncated probability mass: F(max) − F(min), or S(min) − S(max) in the upper tail.</summary>
            internal readonly double Mass;

            /// <summary>True when the base parameters and bounds are valid and the interval has positive mass.</summary>
            internal readonly bool Valid;

            /// <summary>Initializes the cached state.</summary>
            /// <param name="baseDistribution">The live base distribution whose parameters are captured.</param>
            /// <param name="lowerProbability">The base CDF at the lower bound.</param>
            /// <param name="upperProbability">The base CDF at the upper bound.</param>
            /// <param name="lowerSurvival">The base survival probability at the lower bound.</param>
            /// <param name="upperSurvival">The base survival probability at the upper bound.</param>
            /// <param name="valid">Whether the parameters passed validation.</param>
            internal TruncationState(UnivariateDistributionBase baseDistribution, double lowerProbability, double upperProbability,
                double lowerSurvival, double upperSurvival, bool valid)
            {
                _snapshot = DistributionSnapshot.TryCapture(baseDistribution);
                _baseParameters = baseDistribution.GetParameters;
                LowerProbability = lowerProbability;
                LowerSurvival = lowerSurvival;
                UpperSurvival = upperSurvival;
                UpperTail = lowerProbability > 0.5d;
                Mass = UpperTail ? lowerSurvival - upperSurvival : upperProbability - lowerProbability;
                Valid = valid;
            }

            /// <summary>Tests whether the live base parameters are bitwise those captured by this state.</summary>
            /// <param name="baseDistribution">The live base distribution.</param>
            /// <returns>True when the base distribution is unchanged.</returns>
            internal bool Matches(UnivariateDistributionBase baseDistribution)
            {
                if (_snapshot is not null) return _snapshot.Matches(baseDistribution);
                var live = baseDistribution.GetParameters;
                if (live.Length != _baseParameters.Length) return false;
                for (int i = 0; i < live.Length; i++)
                    if (BitConverter.DoubleToInt64Bits(live[i]) != BitConverter.DoubleToInt64Bits(_baseParameters[i])) return false;
                return true;
            }

            /// <summary>Gets the published moments, or null when they have not been computed.</summary>
            internal double[]? Moments => Volatile.Read(ref _moments);

            /// <summary>Publishes the moments computed for this state.</summary>
            /// <param name="moments">The moments {mean, standard deviation, skewness, kurtosis}.</param>
            internal void PublishMoments(double[] moments) => Volatile.Write(ref _moments, moments);
        }

        /// <summary>
        /// Gets the base distribution. 
        /// </summary>
        public UnivariateDistributionBase BaseDistribution => _baseDist;

        /// <summary>
        /// Gets the min of the distribution.
        /// </summary>
        public double Min => _min;

        /// <summary>
        /// Gets the max of the distribution.
        /// </summary>
        public double Max => _max;

        /// <summary>
        /// The name of the XML attribute that holds the lower truncation bound.
        /// </summary>
        internal const string TruncationMinAttribute = "TruncationMin";

        /// <summary>
        /// The name of the XML attribute that holds the upper truncation bound.
        /// </summary>
        internal const string TruncationMaxAttribute = "TruncationMax";

        /// <inheritdoc/>
        /// <remarks>
        /// The element is the base distribution's element with the truncation bounds added as the <c>TruncationMin</c>
        /// and <c>TruncationMax</c> attributes, so the bounds cannot collide with base parameters named Min and Max.
        /// <see cref="UnivariateDistributionFactory.CreateDistribution(XElement)"/> restores the truncated distribution
        /// from it.
        /// </remarks>
        public override XElement ToXElement()
        {
            var result = _baseDist.ToXElement();
            result.SetAttributeValue(TruncationMinAttribute, _min.ToString("G17", CultureInfo.InvariantCulture));
            result.SetAttributeValue(TruncationMaxAttribute, _max.ToString("G17", CultureInfo.InvariantCulture));
            return result;
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return _baseDist.NumberOfParameters + 2; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type => _baseDist.Type;

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Truncated " + _baseDist.DisplayName; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Trunc. " +_baseDist.ShortDisplayName; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[_baseDist.NumberOfParameters + 2, 2];
                for (int i = 0; i < _baseDist.NumberOfParameters; i++)
                {
                    parmString[i, 0] = _baseDist.ParametersToString[i, 0];
                    parmString[i, 1] = _baseDist.ParametersToString[i, 1];
                }
                parmString[NumberOfParameters - 2, 0] = "Min";
                parmString[NumberOfParameters - 1, 0] = "Max";
                parmString[NumberOfParameters - 2, 1] = Min.ToString();
                parmString[NumberOfParameters - 1, 1] = Max.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get
            {
                var parms = _baseDist.ParameterNamesShortForm.ToList();
                parms.AddRange(new[] { "Min", "Max" });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get
            {
                var parms = _baseDist.GetParameterPropertyNames.ToList();
                parms.AddRange(new[] { nameof(Min), nameof(Max) });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get
            {
                var parms = _baseDist.GetParameters.ToList();
                parms.AddRange(new[] { Min, Max });                
                return parms.ToArray(); 
            }
        }


        /// <inheritdoc/>
        /// <remarks>Computed by adaptive integration over [<see cref="Minimum"/>, <see cref="Maximum"/>]; see <see cref="TruncatedMoments"/>.</remarks>
        public override double Mean
        {
            get { return TruncatedMoments()[0]; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5d); }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// When the base distribution is one whose density is unimodal with a closed-form mode (for example
        /// <see cref="Normal"/>, <see cref="GammaDistribution"/> with shape above one, or <see cref="Weibull"/>), and
        /// that mode lies strictly inside the base support, the density increases up to the base mode and decreases after
        /// it. The truncated mode is then the base mode when it lies in [<see cref="Minimum"/>, <see cref="Maximum"/>],
        /// and otherwise the nearer bound, and no search is needed.
        /// </para>
        /// <para>
        /// Otherwise the global mode is searched: the density is evaluated on 65 equally spaced points of
        /// [InverseCDF(0.001), InverseCDF(0.999)], the grid peaks within 0.9 of the highest are polished with a bounded
        /// Brent search over two grid cells. The base distribution's mode, when it lies inside [<see cref="Minimum"/>,
        /// <see cref="Maximum"/>], and the two bounds are then compared with it, and the one with the largest density is
        /// returned, so a mode at a truncation bound is found exactly. A peak narrower than a grid cell that is not the
        /// base mode can be missed.
        /// </para>
        /// </remarks>
        public override double Mode
        {
            get
            {
                var state = CurrentState();
                if (state.Valid == false)
                    ValidateParameters(GetParameters, true);
                if (UnimodalBases.Contains(_baseDist.GetType()))
                {
                    double peak = _baseDist.Mode;
                    if (peak > _baseDist.Minimum && peak < _baseDist.Maximum)
                    {
                        double minimum = Minimum, maximum = Maximum;
                        return peak < minimum ? minimum : peak > maximum ? maximum : peak;
                    }
                }
                double lower = InverseCDF(0.001), upper = InverseCDF(0.999);
                double mode = lower;
                if (upper > lower)
                {
                    // Inside the truncation interval the density is the base density over a constant mass, so the grid
                    // and the polish run on the base density, which has the same maximizer.
                    const int points = 65;
                    var grid = new double[points];
                    double step = (upper - lower) / (points - 1);
                    for (int i = 0; i < points; i++) grid[i] = _baseDist.PDF(lower + i * step);
                    mode = DistributionNumerics.GridMode(_baseDist.PDF, lower, upper, grid, exactGrid: true);
                }
                double density = PDF(mode);
                // A spike narrower than a grid cell can be missed; the base distribution's own mode (a global search for a
                // mixture) is a candidate whenever it lies inside the truncation interval.
                double baseMode = _baseDist.Mode;
                if (baseMode >= Minimum && baseMode <= Maximum)
                {
                    double baseModeDensity = PDF(baseMode);
                    if (baseModeDensity > density)
                    {
                        mode = baseMode;
                        density = baseModeDensity;
                    }
                }
                double lowerBound = Minimum, lowerBoundDensity = PDF(lowerBound);
                if (lowerBoundDensity > density)
                {
                    mode = lowerBound;
                    density = lowerBoundDensity;
                }
                double upperBound = Maximum;
                if (PDF(upperBound) > density) mode = upperBound;
                return mode;
            }
        }

        /// <summary>
        /// The base distribution types whose density is unimodal on its support and whose <see cref="Mode"/> is the
        /// closed-form location of that maximum.
        /// </summary>
        /// <remarks>
        /// A type is listed only when every valid parameter set gives a density that is nondecreasing below the reported
        /// mode and nonincreasing above it. The exact runtime type is matched, so a derived class that overrides the
        /// density or the mode takes the global search. Parameter sets whose mode is a support bound (for example a
        /// gamma shape at or below one) also take the global search.
        /// </remarks>
        private static readonly HashSet<Type> UnimodalBases = new HashSet<Type>
        {
            typeof(Normal), typeof(LogNormal), typeof(LnNormal), typeof(Gumbel), typeof(Logistic), typeof(Rayleigh),
            typeof(Weibull), typeof(GammaDistribution), typeof(ChiSquared), typeof(Triangular), typeof(Pert),
            typeof(TruncatedNormal), typeof(GeneralizedExtremeValue), typeof(GeneralizedLogistic),
            typeof(GeneralizedNormal), typeof(PearsonTypeIII), typeof(LogPearsonTypeIII), typeof(StudentT), typeof(Cauchy)
        };

        /// <inheritdoc/>
        /// <remarks>Computed by adaptive integration over [<see cref="Minimum"/>, <see cref="Maximum"/>]; see <see cref="TruncatedMoments"/>.</remarks>
        public override double StandardDeviation
        {
            get { return TruncatedMoments()[1]; }
        }

        /// <inheritdoc/>
        /// <remarks>Computed by adaptive integration over [<see cref="Minimum"/>, <see cref="Maximum"/>]; see <see cref="TruncatedMoments"/>.</remarks>
        public override double Skewness
        {
            get { return TruncatedMoments()[2]; }
        }

        /// <inheritdoc/>
        /// <remarks>Computed by adaptive integration over [<see cref="Minimum"/>, <see cref="Maximum"/>]; see <see cref="TruncatedMoments"/>.</remarks>
        public override double Kurtosis
        {
            get { return TruncatedMoments()[3]; }
        }

        /// <inheritdoc/>
        public override double Minimum => Math.Max(_baseDist.Minimum, Min);

        /// <inheritdoc/>
        public override double Maximum => Math.Min(_baseDist.Maximum, Max);

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get 
            {
                var parms = _baseDist.MinimumOfParameters.ToList();
                parms.AddRange(new[] { double.NegativeInfinity, Min });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get
            {
                var parms = _baseDist.MaximumOfParameters.ToList();
                parms.AddRange(new[] { Max, double.PositiveInfinity });
                return parms.ToArray();
            }
        }

        /// <inheritdoc/>
        /// <exception cref="ArgumentException">Thrown when the flattened parameter count does not match the base distribution and two truncation bounds.</exception>
        public override void SetParameters(IList<double> parameters)
        {
            if (parameters.Count != NumberOfParameters)
            {
                throw new ArgumentException("The length of the parameter array is invalid.", nameof(parameters));
            }
            _baseDist.SetParameters(parameters.ToArray().Subset(0, parameters.Count - 3));
            _min = parameters[parameters.Count - 2];
            _max = parameters[parameters.Count - 1];
            var state = CreateState(parameters);
            _parametersValid = state.Valid;
            Volatile.Write(ref _state, state);
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParametersCore(parameters, throwException, out _, out _, out _, out _);
        }

        /// <summary>
        /// Validates a flattened parameter vector and captures the resulting truncation state of the live base distribution.
        /// </summary>
        /// <param name="parameters">The base parameters followed by the lower and upper truncation bounds.</param>
        /// <returns>The cached state for the current base parameters.</returns>
        private TruncationState CreateState(IList<double> parameters)
        {
            var validation = ValidateParametersCore(parameters, false, out double lowerProbability, out double upperProbability,
                out double lowerSurvival, out double upperSurvival);
            return new TruncationState(_baseDist, lowerProbability, upperProbability, lowerSurvival, upperSurvival, validation is null);
        }

        /// <summary>
        /// Returns the cached state, rebuilding it and refreshing <see cref="UnivariateDistributionBase.ParametersValid"/>
        /// when the base distribution has changed since the state was captured.
        /// </summary>
        /// <returns>The state that matches the live base parameters.</returns>
        private TruncationState CurrentState()
        {
            var state = Volatile.Read(ref _state);
            if (state is not null && state.Matches(_baseDist)) return state;
            state = CreateState(GetParameters);
            _parametersValid = state.Valid;
            Volatile.Write(ref _state, state);
            return state;
        }

        /// <summary>
        /// Returns the moments {mean, standard deviation, skewness, kurtosis}, computed once per base-parameter state.
        /// </summary>
        /// <returns>The mean, standard deviation, skewness and kurtosis.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the parameters are invalid.</exception>
        /// <remarks>
        /// The density is integrated with <see cref="AdaptiveGaussKronrod"/> over the exact support
        /// [<see cref="Minimum"/>, <see cref="Maximum"/>], mapped to [-1, 1] about its midpoint, with relative and absolute
        /// tolerances of 1E-12 and 1E-15. The mean is integrated first, then the variance about it, then the third and
        /// fourth standardized central moments. A moment whose integration fails is returned as <see cref="double.NaN"/>.
        /// </remarks>
        private double[] TruncatedMoments()
        {
            var state = CurrentState();
            var cached = state.Moments;
            if (cached is not null) return cached;
            if (!state.Valid) ValidateParameters(GetParameters, true);
            double lower = Minimum, upper = Maximum;
            double center = lower / 2d + upper / 2d, half = upper / 2d - lower / 2d;
            double Integrate(Func<double, double> weight)
            {
                var agk = new AdaptiveGaussKronrod(t => weight(t) * PDF(center + half * t) * half, -1d, 1d)
                {
                    RelativeTolerance = 1E-12,
                    AbsoluteTolerance = 1E-15,
                    ReportFailure = false
                };
                agk.Integrate();
                return agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            }
            double meanT = Integrate(t => t);
            double sdT = Math.Sqrt(Integrate(t => (t - meanT) * (t - meanT)));
            double skew = Integrate(t => Math.Pow((t - meanT) / sdT, 3d));
            double kurt = Integrate(t => Math.Pow((t - meanT) / sdT, 4d));
            var moments = new[] { center + half * meanT, half * sdT, skew, kurt };
            state.PublishMoments(moments);
            return moments;
        }

        /// <summary>
        /// Validates a flattened base-distribution parameter vector and its truncation bounds.
        /// </summary>
        /// <param name="parameters">The base parameters followed by the lower and upper truncation bounds.</param>
        /// <param name="throwException">Whether to throw the first validation exception.</param>
        /// <param name="lowerProbability">The base-distribution CDF at the lower bound when valid; otherwise <see cref="double.NaN"/>.</param>
        /// <param name="upperProbability">The base-distribution CDF at the upper bound when valid; otherwise <see cref="double.NaN"/>.</param>
        /// <param name="lowerSurvival">The base-distribution CCDF at the lower bound when the lower-bound CDF exceeds one half; otherwise <see cref="double.NaN"/>.</param>
        /// <param name="upperSurvival">The base-distribution CCDF at the upper bound when the lower-bound CDF exceeds one half; otherwise <see cref="double.NaN"/>.</param>
        /// <returns><see langword="null"/> when all parameters are valid; otherwise, the validation exception.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="throwException"/> is true and validation fails.</exception>
        /// <remarks>
        /// The interval mass is the CDF difference when the lower-bound CDF is at most one half and the survival difference
        /// otherwise; either way the interval is rejected when its mass is below 1E-15.
        /// </remarks>
        private ArgumentOutOfRangeException? ValidateParametersCore(IList<double> parameters, bool throwException, out double lowerProbability, out double upperProbability,
            out double lowerSurvival, out double upperSurvival)
        {
            lowerProbability = double.NaN;
            upperProbability = double.NaN;
            lowerSurvival = double.NaN;
            upperSurvival = double.NaN;
            if (parameters.Count != NumberOfParameters)
            {
                var exception = new ArgumentOutOfRangeException(nameof(parameters), "The flattened parameter count is invalid.");
                if (throwException) throw exception;
                return exception;
            }
            var baseParameters = parameters.ToArray().Subset(0, parameters.Count - 3);
            var baseValidation = _baseDist.ValidateParameters(baseParameters, false);
            if (baseValidation is not null)
            {
                if (throwException) _baseDist.ValidateParameters(baseParameters, true);
                return baseValidation;
            }
            double min = parameters[parameters.Count - 2];
            double max = parameters[parameters.Count - 1];
            if (double.IsNaN(min) || double.IsNaN(max) || double.IsInfinity(min) || double.IsInfinity(max) || min >= max)
            {
                var exception = new ArgumentOutOfRangeException(nameof(Min), "The min must be finite and less than the max.");
                if (throwException) throw exception;
                return exception;
            }
            var candidate = _baseDist.Clone();
            candidate.SetParameters(baseParameters);
            lowerProbability = candidate.CDF(min);
            upperProbability = candidate.CDF(max);
            double mass = lowerProbability - upperProbability;
            if (lowerProbability > 0.5d)
            {
                lowerSurvival = candidate.CCDF(min);
                upperSurvival = candidate.CCDF(max);
                mass = lowerSurvival - upperSurvival;
            }
            if (Math.Abs(mass) < 1e-15)
            {
                var exception = new ArgumentOutOfRangeException(nameof(Min), "Truncation interval has zero probability mass.");
                if (throwException) throw exception;
                return exception;
            }
            return null;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            var state = CurrentState();
            if (state.Valid == false)
                ValidateParameters(GetParameters, true);
            if (x < Min || x > Max) return 0.0;
            return _baseDist.PDF(x) / state.Mass;
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            var state = CurrentState();
            if (state.Valid == false)
                ValidateParameters(GetParameters, true);
            if (x <= Min) return 0.0;
            if (x >= Max) return 1.0;
            if (state.UpperTail) return (state.LowerSurvival - _baseDist.CCDF(x)) / state.Mass;
            return (_baseDist.CDF(x) - state.LowerProbability) / state.Mass;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Equals 1 − <see cref="CDF(double)"/> when the base CDF at the lower bound is at most one half, and
        /// (S(x) − S(max)) / (S(min) − S(max)) with the base survival function S otherwise.
        /// </remarks>
        public override double CCDF(double x)
        {
            // Validate parameters
            var state = CurrentState();
            if (state.Valid == false)
                ValidateParameters(GetParameters, true);
            if (!state.UpperTail) return 1.0d - CDF(x);
            if (x <= Min) return 1.0;
            if (x >= Max) return 0.0;
            return (_baseDist.CCDF(x) - state.UpperSurvival) / state.Mass;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// When the base CDF at the lower bound exceeds one half, the quantile solves S(x) = S(max) + (1 − p)(S(min) − S(max))
        /// for x with Brent's method on [<see cref="Minimum"/>, <see cref="Maximum"/>], where S is the base survival function,
        /// because the base quantile of a probability near one cannot resolve the upper tail.
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate parameters
            var state = CurrentState();
            if (state.Valid == false)
                ValidateParameters(GetParameters, true);
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d) return Minimum;
            if (probability == 1.0d) return Maximum;
            if (state.UpperTail)
            {
                double target = state.UpperSurvival + (1d - probability) * state.Mass;
                double lower = Minimum, upper = Maximum;
                return Brent.Solve(x => _baseDist.CCDF(x) - target, lower, upper, 1E-15 * (upper - lower), 1000, false);
            }
            return _baseDist.InverseCDF(probability * state.Mass + state.LowerProbability);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The base distribution is cloned as well, so changing the parameters of the clone does not change this instance.
        /// </remarks>
        public override UnivariateDistributionBase Clone()
        {
            return new TruncatedDistribution(_baseDist.Clone(), Min, Max);
        }

    }
}
