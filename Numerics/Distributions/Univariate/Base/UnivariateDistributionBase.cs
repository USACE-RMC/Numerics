using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// Declares common functionality for all univariate distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class UnivariateDistributionBase : IUnivariateDistribution, IComparable<UnivariateDistributionBase>
    {
        /// <summary>
        /// Approximation for assessing if a parameter is near zero. 
        /// </summary>
        protected double NearZero = 1E-4;

        /// <summary>
        /// Protected parameter is valid property. 
        /// </summary>
        protected bool _parametersValid = true;

        /// <summary>
        /// Returns the continuous distribution type.
        /// </summary>
        public abstract UnivariateDistributionType Type { get; }

        /// <inheritdoc/>
        public abstract string DisplayName { get; }

        /// <inheritdoc/>
        public abstract string ShortDisplayName { get; }

        /// <inheritdoc/>
        public string DisplayLabel
        {
            get
            {
                StringBuilder label = new StringBuilder(ShortDisplayName + " (");
                var paramVals = GetParameters;
                for (int i = 0; i < paramVals.Length; i++)
                {
                    label.Append(string.Format("{0:0.#####}", paramVals[i]));
                    if (i != paramVals.Length - 1)
                        label.Append(", ");
                }
                label.Append(")");
                return label.ToString();
            }
        }

        /// <inheritdoc/>
        public abstract int NumberOfParameters { get; }

        /// <inheritdoc/>
        public abstract string[,] ParametersToString { get; }

        /// <summary>
        /// Returns the distribution parameter names as an array of string.
        /// </summary>
        public virtual string[] ParameterNames
        {
            get { return ParametersToString.GetColumn(0); }
        }

        /// <inheritdoc/>
        public abstract string[] ParameterNamesShortForm { get; }

        /// <inheritdoc/>
        public abstract double[] GetParameters { get; }

        /// <inheritdoc/>
        public abstract string[] GetParameterPropertyNames { get; }

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public abstract double Mean { get; }

        /// <inheritdoc/>
        public abstract double Median { get; }

        /// <inheritdoc/>
        public abstract double Mode { get; }

        /// <inheritdoc/>
        public double Variance
        {
            get { return Math.Pow(StandardDeviation, 2d); }
        }

        /// <inheritdoc/>
        public abstract double StandardDeviation { get; }

        /// <inheritdoc/>
        public double CoefficientOfVariation
        {
            get { return StandardDeviation / Mean; }
        }

        /// <inheritdoc/>
        public abstract double Skewness { get; }

        /// <inheritdoc/>
        public abstract double Kurtosis { get; }

        /// <inheritdoc/>
        public abstract double Minimum { get; }

        /// <inheritdoc/>
        public abstract double Maximum { get; }

        /// <inheritdoc/>
        public abstract double[] MinimumOfParameters { get; }

        /// <inheritdoc/>
        public abstract double[] MaximumOfParameters { get; }

        /// <inheritdoc/>
        public abstract void SetParameters(IList<double> parameters);

        /// <inheritdoc/>
        public abstract ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException);

        /// <summary>
        /// The log likelihood function.
        /// </summary>
        /// <param name="sample">Sample of observed data.</param>
        /// <returns>
        /// Log likelihood.
        /// </returns>
        public double LogLikelihood(IList<double> sample)
        {
            var logLH = default(double);
            for (int i = 0; i < sample.Count; i++)
            {        
                logLH += LogPDF(sample[i]);
            }
            if (double.IsNaN(logLH) || double.IsInfinity(logLH)) return double.NegativeInfinity;
            return logLH;
        }

        /// <summary>
        /// The Log-Likelihood function for a single data point.
        /// </summary>
        /// <param name="value">A single observed data point.</param>
        /// <returns>
        /// Log likelihood.
        /// </returns>
        public double LogLikelihood(double value)
        {
            return LogPDF(value);
        }

        /// <summary>
        /// The log likelihood function for left censored data.
        /// </summary>
        /// <param name="threshold">The threshold.</param>
        /// <param name="numberBelow">The number of data points below the threshold.</param>
        /// <returns>The censored log likelihood; zero when the category is empty.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The count is negative.</exception>
        public double LogLikelihood_LeftCensored(double threshold, long numberBelow)
        {
            if (numberBelow < 0) throw new ArgumentOutOfRangeException(nameof(numberBelow));
            if (numberBelow == 0) return 0;
            return numberBelow * LogCDF(threshold);
        }

        /// <summary>
        /// The log likelihood function for right censored data.
        /// </summary>
        /// <param name="threshold">The threshold.</param>
        /// <param name="numberAbove">The number of data points above the threshold.</param>
        /// <returns>The censored log likelihood; zero when the category is empty.</returns>
        /// <exception cref="ArgumentOutOfRangeException">The count is negative.</exception>
        public double LogLikelihood_RightCensored(double threshold, long numberAbove)
        {
            if (numberAbove < 0) throw new ArgumentOutOfRangeException(nameof(numberAbove));
            if (numberAbove == 0) return 0;
            return numberAbove * LogCCDF(threshold);
        }

        /// <summary>
        /// The log likelihood function for interval data.
        /// </summary>
        /// <param name="lowerLimit">The lower limit of the interval.</param>
        /// <param name="upperLimit">The upper limit of the interval.</param>
        /// <returns>The logarithm of the probability of (lowerLimit, upperLimit].</returns>
        /// <exception cref="ArgumentOutOfRangeException">The limits are reversed or NaN.</exception>
        /// <remarks>Uses lower or upper logarithmic tails to avoid subtracting rounded probabilities.
        /// The open lower and closed upper endpoint convention also applies to atoms.</remarks>
        public double LogLikelihood_Intervals(double lowerLimit, double upperLimit)
        {
            if (double.IsNaN(lowerLimit) || double.IsNaN(upperLimit) || lowerLimit > upperLimit)
                throw new ArgumentOutOfRangeException(nameof(upperLimit), "Interval limits must be ordered and not NaN.");
            if (lowerLimit == upperLimit) return double.NegativeInfinity;
            if (this is Mixture mixture) return mixture.LogIntervalProbability(lowerLimit, upperLimit);
            double lower = LogCDF(lowerLimit), upper = LogCDF(upperLimit);
            double result = lower < -0.6931471805599453
                ? DistributionNumerics.LogDifference(upper, lower)
                : DistributionNumerics.LogDifference(LogCCDF(lowerLimit), LogCCDF(upperLimit));
            return double.IsNegativeInfinity(result)
                ? DistributionNumerics.CollapsedContinuousLogInterval(this, lowerLimit, upperLimit) : result;
        }

        /// <inheritdoc/>
        public abstract double PDF(double x);

        /// <inheritdoc/>
        public virtual double LogPDF(double x)
        {
            double f = PDF(x);
            if (double.IsNaN(f) || double.IsInfinity(f) || f <= 0.0)
                return double.NegativeInfinity;
            return Math.Log(f);
        }

        /// <inheritdoc/>
        public abstract double CDF(double x);

        /// <inheritdoc/>
        public virtual double HF(double x)
        {
            return PDF(x) / CCDF(x);
        }

        /// <inheritdoc/>
        public virtual double LogCDF(double x)
        {
            double F = CDF(x);
            if (double.IsNaN(F) || double.IsInfinity(F) || F <= 0.0)
                return double.NegativeInfinity;
            return Math.Log(F);
        }

        /// <inheritdoc/>
        public virtual double CCDF(double x)
        {
            return 1.0d - CDF(x);
        }

        /// <inheritdoc/>
        public virtual double LogCCDF(double x)
        {
            double cF = CCDF(x);
            if (double.IsNaN(cF) || double.IsInfinity(cF) || cF <= 0.0)
                return double.NegativeInfinity;
            return Math.Log(cF);
        }

        /// <inheritdoc/>
        public abstract double InverseCDF(double probability);

        /// <summary>
        /// The inverse of the complementary cumulative distribution function (the survival or
        /// exceedance quantile): the value x at which CCDF(x) = <paramref name="probability"/>.
        /// </summary>
        /// <param name="probability">The survival (exceedance) probability q, between 0 and 1.</param>
        /// <returns>The value exceeded with probability q, the same quantile as InverseCDF(1 − q).</returns>
        /// <remarks>
        /// <para>
        /// The base implementation returns <see cref="InverseCDF(double)"/> at 1 − q, so its result is
        /// bitwise that of the existing quantile. Forming 1 − q rounds away the low-order digits of a
        /// small q, and every q below about 1.1E-16 rounds to the upper support limit. Distributions
        /// with a closed-form survival quantile override this member to evaluate it from q directly
        /// (for example <see cref="GeneralizedPareto.InverseCCDF(double)"/>).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the probability is not between 0 and 1.</exception>
        public virtual double InverseCCDF(double probability)
        {
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException(nameof(probability), "Probability must be between 0 and 1.");
            return InverseCDF(1.0d - probability);
        }

        /// <summary>
        /// The PDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        public double[] PDF(IList<double> xValues)
        {
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = PDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The CDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        public double[] CDF(IList<double> xValues)
        {
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = CDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The CCDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        /// <returns>One exceedance probability per x value, each equal to <see cref="CCDF(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xValues"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="CDF(IList{double})"/>; a caller crossing
        /// a language boundary (for example from Python) pays the crossing once instead of once per value.
        /// </remarks>
        public double[] CCDF(IList<double> xValues)
        {
            if (xValues == null) throw new ArgumentNullException(nameof(xValues));
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = CCDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The natural logarithm of the PDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        /// <returns>One log density per x value, each equal to <see cref="LogPDF(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xValues"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="PDF(IList{double})"/>.
        /// </remarks>
        public double[] LogPDF(IList<double> xValues)
        {
            if (xValues == null) throw new ArgumentNullException(nameof(xValues));
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = LogPDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The natural logarithm of the CDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        /// <returns>One log non-exceedance probability per x value, each equal to <see cref="LogCDF(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xValues"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="CDF(IList{double})"/>.
        /// </remarks>
        public double[] LogCDF(IList<double> xValues)
        {
            if (xValues == null) throw new ArgumentNullException(nameof(xValues));
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = LogCDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The natural logarithm of the CCDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        /// <returns>One log exceedance probability per x value, each equal to <see cref="LogCCDF(double)"/> at that value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="xValues"/> is null.</exception>
        /// <remarks>
        /// One call evaluates the whole list, in the same way as <see cref="CDF(IList{double})"/>.
        /// </remarks>
        public double[] LogCCDF(IList<double> xValues)
        {
            if (xValues == null) throw new ArgumentNullException(nameof(xValues));
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = LogCCDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The Inverse CDF evaluated over a list of probabilities.
        /// </summary>
        /// <param name="probabilities">List of probability values.</param>
        public double[] InverseCDF(IList<double> probabilities)
        {
            var result = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                result[i] = InverseCDF(probabilities[i]);
            return result;
        }

        /// <summary>
        /// The inverse CCDF (survival quantile) evaluated over a list of survival probabilities.
        /// </summary>
        /// <param name="probabilities">The survival (exceedance) probabilities, each between 0 and 1.</param>
        /// <returns>One quantile per probability, each equal to <see cref="InverseCCDF(double)"/> at that probability.</returns>
        /// <remarks>
        /// One call evaluates the whole list, so a caller crossing a language boundary (for example
        /// from Python) pays the crossing once instead of once per value.
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="probabilities"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a probability is not between 0 and 1.</exception>
        public double[] InverseCCDF(IList<double> probabilities)
        {
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            var result = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                result[i] = InverseCCDF(probabilities[i]);
            return result;
        }

        /// <summary>
        /// The Inverse CDF evaluated over a list of probabilities under each of many parameter sets of this distribution's family.
        /// </summary>
        /// <param name="parameterSets">The parameter sets, one row per set and one column per parameter, in the order of <see cref="GetParameters"/>.</param>
        /// <param name="probabilities">List of non-exceedance probabilities.</param>
        /// <param name="valid">On return, one flag per row: true when <see cref="ValidateParameters(IList{double}, bool)"/> accepts the row.</param>
        /// <returns>
        /// A 2-D array of quantiles with one row per parameter set and one column per probability; the rows of
        /// invalid parameter sets are <see cref="double.NaN"/>.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="parameterSets"/> or <paramref name="probabilities"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the number of columns differs from <see cref="NumberOfParameters"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown by <see cref="InverseCDF(double)"/> when a probability is outside [0, 1] and a row is valid.</exception>
        /// <exception cref="NotImplementedException">Thrown by distributions that do not support <see cref="SetParameters(IList{double})"/>,
        /// such as <see cref="EmpiricalDistribution"/> and <see cref="KernelDensity"/>, when a row is valid.</exception>
        /// <remarks>
        /// <para>
        /// Each row is evaluated on one private <see cref="Clone"/> of this distribution, so this instance is never
        /// modified and settings that are not parameters (for example, the logarithm base of a log-Pearson Type III)
        /// carry over to every row. A wrapper distribution that holds another distribution must therefore clone it
        /// deeply; a clone that shares the wrapped instance lets the rows modify the original. A row is valid when <c>ValidateParameters(row, false)</c> returns null; a valid
        /// row is applied with <c>SetParameters(row)</c> and evaluated with <see cref="InverseCDF(double)"/> at each
        /// probability. The arithmetic is therefore exactly that of the per-row loop
        /// <c>ValidateParameters</c> → <c>SetParameters</c> → <c>InverseCDF</c> on a clone, and the results are
        /// bitwise identical to it. Invalid rows are not applied and their quantiles are <see cref="double.NaN"/>,
        /// matching the convention of <see cref="BootstrapAnalysis.Quantiles(IList{double}, IUnivariateDistribution[])"/>.
        /// </para>
        /// <para>
        /// The rows are evaluated serially in order, so the result does not depend on the machine or the thread
        /// count. This member lets a caller evaluate a sampled parameter ensemble in one call rather than one call
        /// per row.
        /// </para>
        /// <para>
        /// A family may override this member with a faster evaluation (for example, computing each probability's
        /// transform once for all rows) only when its results stay bitwise identical to this per-row loop;
        /// <see cref="GeneralizedExtremeValue"/> does.
        /// </para>
        /// </remarks>
        public virtual double[,] ParameterSetQuantiles(double[,] parameterSets, IList<double> probabilities, out bool[] valid)
        {
            if (parameterSets == null) throw new ArgumentNullException(nameof(parameterSets));
            if (probabilities == null) throw new ArgumentNullException(nameof(probabilities));
            int count = parameterSets.GetLength(0);
            int arity = NumberOfParameters;
            if (parameterSets.GetLength(1) != arity)
                throw new ArgumentException("Each parameter set must have exactly " + arity + " parameters.", nameof(parameterSets));

            var result = new double[count, probabilities.Count];
            valid = new bool[count];
            var probe = Clone();
            var row = new double[arity];
            for (int i = 0; i < count; i++)
            {
                for (int j = 0; j < arity; j++)
                    row[j] = parameterSets[i, j];
                valid[i] = probe.ValidateParameters(row, false) is null;
                if (!valid[i])
                {
                    for (int k = 0; k < probabilities.Count; k++)
                        result[i, k] = double.NaN;
                    continue;
                }
                probe.SetParameters(row);
                for (int k = 0; k < probabilities.Count; k++)
                    result[i, k] = probe.InverseCDF(probabilities[k]);
            }
            return result;
        }

        /// <summary>
        /// Returns the central moments {Mean, Standard Deviation, Skew, and Kurtosis} of the distribution using adaptive numerical integration with the Gauss-Kronrod rule.
        /// </summary>
        /// <param name="tolerance">The desired relative tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.</param>
        /// <returns>Mean, Standard Deviation, Skew, and Kurtosis.</returns>
        /// <remarks>
        /// <para>
        /// <b>This overload takes a tolerance, not a step count.</b> There is a sibling overload,
        /// <see cref="CentralMoments(int)"/>, whose argument is a number of integration steps. The two are
        /// distinguished only by the type of the single argument, so the choice is invisible at the call
        /// site: <c>CentralMoments(1E-8)</c> requests a relative tolerance of 1E-8 from this adaptive
        /// routine, while <c>CentralMoments(1000)</c> selects the fixed-step bin-expectation overload with 1,000
        /// steps. Writing <c>CentralMoments(1)</c> when a tolerance of 1.0 was meant silently runs the other
        /// method with a single step.
        /// </para>
        /// <para>
        /// This overload integrates with <see cref="AdaptiveGaussKronrod"/> over the range
        /// [InverseCDF(1E-16), InverseCDF(1 − 1E-16)], subdividing until the requested relative tolerance is
        /// met. A moment whose integration fails is returned as <see cref="double.NaN"/> rather than
        /// throwing.
        /// </para>
        /// <para>This inherited approximation is independent of analytical moment properties overridden
        /// by individual distributions; use those properties when an analytical result is available.</para>
        /// </remarks>
        public virtual double[] CentralMoments(double tolerance = 1E-8)
        {
            double u1 = 0, u2 = 0, u3 = 0, u4 = 0;
            double a = InverseCDF(1E-16);
            double b = InverseCDF(1-1E-16);
            if (a >= b) return new double[] { a, double.NaN, double.NaN, double.NaN };

            // Mean
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u1 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            // Standard Deviation
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - u1, 2d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u2 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            u2 = Math.Sqrt(u2);
            // Skewness
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow((x - u1) / u2, 3d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u3 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            // Kurtosis
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow((x - u1) / u2, 4d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u4 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;

            return new double[] { u1, u2, u3, u4 };
        }

        /// <summary>
        /// Returns the central moments {Mean, Standard Deviation, Skew, and Kurtosis} of the distribution using a fixed-step discrete expectation over stratified bins.
        /// </summary>
        ///<param name="steps">Number of integration steps. Default = 300.</param>
        /// <returns>Mean, Standard Deviation, Skew, and Kurtosis.</returns>
        /// <remarks>
        /// <para>
        /// <b>This overload takes a step count, not a tolerance.</b> There is a sibling overload,
        /// <see cref="CentralMoments(double)"/>, whose argument is an integration tolerance. The two are
        /// distinguished only by the type of the single argument, so the choice is invisible at the call
        /// site: <c>CentralMoments(1000)</c> selects this fixed-step routine with 1,000 steps, while
        /// <c>CentralMoments(1E-8)</c> selects the adaptive overload with a relative tolerance of 1E-8.
        /// Writing <c>CentralMoments(1)</c> when a tolerance was meant silently runs this method with a
        /// single step.
        /// </para>
        /// <para>
        /// <b>What this overload computes:</b> it stratifies
        /// [InverseCDF(1E-8), InverseCDF(1 − 1E-8)] into <paramref name="steps"/> equal bins, takes each
        /// bin's probability mass ΔFᵢ as a difference of <see cref="CDF(double)"/> values, and accumulates
        /// the discrete expectation Σᵢ xᵢᵏ · ΔFᵢ. It is therefore a bin-probability (midpoint) expectation
        /// against the distribution, <b>not</b> a trapezoidal rule applied to the integrand x·f(x). The
        /// representative point xᵢ is the bin midpoint for interior bins, the upper bound for the first bin
        /// and the lower bound for the last. Moments are accumulated in local coordinates about the
        /// resulting mean to avoid cancellation between large raw first and second moments.
        /// </para>
        /// <para>
        /// The first and last bins carry the whole of their tails: ΔF₀ is CDF(upper bound of bin 0), taken
        /// from −∞, and the final ΔF is 1 − CDF(lower bound of the last bin). So despite the 1E-8 stratification
        /// endpoints the total probability sums to one and the effective range is <b>not</b> truncated. Cost
        /// and accuracy are both fixed by <paramref name="steps"/> and there is no convergence check.
        /// </para>
        /// <para>This legacy approximation does not dispatch to analytical moment properties. Its
        /// fixed tail representatives can be inaccurate for heavy tails. Prefer the distribution's
        /// moment properties when analytical values are available.</para>
        /// </remarks>
        public virtual double[] CentralMoments(int steps = 300)
        {
            double a = InverseCDF(1E-8);
            double b = InverseCDF(1 - 1E-8);
            if (a >= b) return [a, double.NaN, double.NaN, double.NaN];

            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));
            var dFx = new double[steps];
            var coordinates = new double[steps];
            double reference = a / 2 + b / 2;
            double scale = Math.Max(Math.Abs(a - reference), Math.Abs(b - reference));
            // Preserve the existing bins and probability masses, including endpoint tail representatives.
            // A bin whose lower bound is bitwise the previous upper bound reuses that CDF value.
            double previousUpper = bins[0].UpperBound;
            double previousCDF = CDF(previousUpper);
            dFx[0] = previousCDF;
            coordinates[0] = DistributionNumerics.Standardize(bins[0].UpperBound, reference, scale);
            for (int i = 1; i < steps - 1; i++)
            {
                double upperCDF = CDF(bins[i].UpperBound);
                dFx[i] = upperCDF - (SameBits(bins[i].LowerBound, previousUpper) ? previousCDF : CDF(bins[i].LowerBound));
                previousUpper = bins[i].UpperBound;
                previousCDF = upperCDF;
                coordinates[i] = DistributionNumerics.Standardize(bins[i].Midpoint, reference, scale);
            }
            double lastLower = bins.Last().LowerBound;
            dFx[steps - 1] = 1 - (steps > 1 && SameBits(lastLower, previousUpper) ? previousCDF : CDF(lastLower));
            coordinates[steps - 1] = DistributionNumerics.Standardize(bins.Last().LowerBound, reference, scale);
            double offset = 0;
            for (int i = 0; i < steps; i++) offset += coordinates[i] * dFx[i];
            double second = 0, third = 0, fourth = 0;
            for (int i = 0; i < steps; i++)
            {
                double centered = coordinates[i] - offset, square = centered * centered;
                second += square * dFx[i];
                third += square * centered * dFx[i];
                fourth += square * square * dFx[i];
            }
            return [reference + scale * offset, scale * Math.Sqrt(second), third / second / Math.Sqrt(second), fourth / second / second];
        }

        /// <summary>Tests whether two bin edges are the same floating-point value bit for bit.</summary>
        /// <param name="a">The first edge.</param>
        /// <param name="b">The second edge.</param>
        /// <returns><see langword="true"/> when the raw bits are equal, so a CDF value at one edge is exactly the value at the other.</returns>
        private static bool SameBits(double a, double b) => BitConverter.DoubleToInt64Bits(a) == BitConverter.DoubleToInt64Bits(b);

        /// <summary>
        /// Returns the conditional mean and the second to fourth conditional moments about the unconditional mean,
        /// given that X lies in [a, b].
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <returns>
        /// A 4-element array: m1 = E[X | a &lt;= X &lt;= b], then mk = E[(X - μ)^k | a &lt;= X &lt;= b] for k = 2, 3, 4,
        /// where μ is the unconditional <see cref="Mean"/>. All four are NaN when a &gt;= b or the interval has no
        /// probability mass.
        /// </returns>
        /// <remarks>
        /// m2 to m4 are centred about the unconditional mean μ, not about the conditional mean m1, so they are not
        /// the central moments of the truncated distribution. The truncated variance, for example, is
        /// m2 - (m1 - μ)^2. The expectation uses a 300-bin midpoint rule with the CDF mass of each bin.
        /// </remarks>
        public virtual double[] ConditionalMoments(double a, double b)
        {
            if (a >= b)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            int steps = 300;
            double l = CDF(a), u = CDF(b), totalProb = u - l;
            if (totalProb <= 0.0 || double.IsNaN(totalProb))
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            double mean = Mean;
            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));

            // Midpoint integration
            double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
            // A bin whose lower bound is bitwise the previous upper bound reuses that CDF value.
            double previousUpper = double.NaN, previousCDF = double.NaN;

            for (int i = 0; i < steps; i++)
            {
                double upperCDF = CDF(bins[i].UpperBound);
                double dF = upperCDF - (i > 0 && SameBits(bins[i].LowerBound, previousUpper) ? previousCDF : CDF(bins[i].LowerBound));
                previousUpper = bins[i].UpperBound;
                previousCDF = upperCDF;
                double xMid = bins[i].Midpoint;

                double delta = xMid - mean;
                double delta2 = delta * delta;
                double delta3 = delta2 * delta;
                double delta4 = delta2 * delta2;

                sum1 += xMid * dF;
                sum2 += delta2 * dF;
                sum3 += delta3 * dF;
                sum4 += delta4 * dF;
            }

            double m1 = sum1 / totalProb;
            double m2 = sum2 / totalProb;
            double m3 = sum3 / totalProb;
            double m4 = sum4 / totalProb;

            return new[] { m1, m2, m3, m4 };
        }

        /// <summary>
        /// Returns the conditional mean of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalMean(double a, double b, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional variance of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalVariance(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 2d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional skewness of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalSkewness(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 3d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional kurtosis of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalKurtosis(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 4d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional expected value at the given confidence level. 
        /// </summary>
        /// <param name="alpha">The threshold confidence level, α.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.</param>
        public virtual double ConditionalExpectedValue(double alpha, double tolerance = 1E-8)
        {
            double a = InverseCDF(alpha);
            double b = InverseCDF(1 - 1E-16);
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (1 - alpha);
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and probability density.</returns>
        public virtual double[,] CreatePDFGraph()
        {
            var XValues = Stratify.XValues(new StratificationOptions(InverseCDF(0.0005d), InverseCDF(0.9995d), 1000, false));
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = PDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Stratified X values.</param>
        /// <returns>A 2-column array of X and probability density.</returns>
        public double[,] CreatePDFGraph(List<StratificationBin> XValues)
        {
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = PDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Array of stratified X values.</param>
        /// <returns>A 2-column array of X and probability density, one row per bin, with the ranges concatenated in order.</returns>
        public double[,] CreatePDFGraph(List<StratificationBin>[] XValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < XValues.Length; i++)
                totalBinCount += XValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < XValues.Length; i++)
            {
                for (int j = 0; j < XValues[i].Count; j++)
                {
                    graph[binCounter, 0] = XValues[i][j].Midpoint;
                    graph[binCounter, 1] = PDF(XValues[i][j].Midpoint);
                    binCounter += 1;
                }
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateCDFGraph()
        {
            var XValues = Stratify.XValues(new StratificationOptions(InverseCDF(0.0005d), InverseCDF(0.9995d), 1000, false));
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = CDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Stratified X values.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateCDFGraph(List<StratificationBin> XValues)
        {
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = CDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Array of stratified X values.</param>
        /// <returns>A 2-column array of X and non-exceedance probability, one row per bin, with the ranges concatenated in order.</returns>
        public double[,] CreateCDFGraph(List<StratificationBin>[] XValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < XValues.Length; i++)
                totalBinCount += XValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < XValues.Length; i++)
            {
                for (int j = 0; j < XValues[i].Count; j++)
                {
                    graph[binCounter, 0] = XValues[i][j].Midpoint;
                    graph[binCounter, 1] = CDF(XValues[i][j].Midpoint);
                    binCounter += 1;
                }
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph()
        {
            var PValues = Stratify.Probabilities(new StratificationOptions(CDF(InverseCDF(0.0005d)), CDF(InverseCDF(0.9995d)), 1000, true));
            var graph = new double[PValues.Count, 2];
            for (int i = 0; i < PValues.Count; i++)
            {
                graph[i, 0] = InverseCDF(PValues[i].Midpoint);
                graph[i, 1] = PValues[i].Midpoint;
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// </summary>
        /// <param name="PValues">Stratified probabilities.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph(List<StratificationBin> PValues)
        {
            var graph = new double[PValues.Count, 2];
            for (int i = 0; i < PValues.Count; i++)
            {
                graph[i, 0] = InverseCDF(PValues[i].Midpoint);
                graph[i, 1] = PValues[i].Midpoint;
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// </summary>
        /// <param name="PValues">Array of stratified probabilities.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph(List<StratificationBin>[] PValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < PValues.Length; i++)
                totalBinCount += PValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < PValues.Length; i++)
            {
                for (int j = 0; j < PValues[i].Count; j++)
                {
                    graph[binCounter, 0] = InverseCDF(PValues[i][j].Midpoint);
                    graph[binCounter, 1] = PValues[i][j].Midpoint;
                    binCounter += 1;
                }
            }
            return graph;
        }

        /// <summary>
        /// Generate random values of a distribution given a sample size.
        /// </summary>
        /// <param name="sampleSize">Size of random sample to generate.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <returns>
        /// Array of random values.
        /// </returns>
        public virtual double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            // Create PRNG for generating random numbers
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
                sample[i] = InverseCDF(rnd.NextDouble());
            // Return array of random values
            return sample;
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public abstract UnivariateDistributionBase Clone();

        /// <summary>
        /// Returns the distribution as XElement (XML).
        /// </summary>
        /// <returns>An element holding the type, each flattened parameter, and any non-default settings outside the parameters.</returns>
        /// <remarks>
        /// Settings that change evaluated values but are not parameters are written as optional attributes only when
        /// they differ from the family default: the logarithm <c>Base</c> of <see cref="LogNormal"/> and
        /// <see cref="LogPearsonTypeIII"/>, and the natural-log coordinates <c>Mu</c> and <c>Sigma</c> of an
        /// <see cref="LnNormal"/> whose physical-moment surface is inactive. Default-configured distributions
        /// therefore serialize exactly as elements written without these attributes, and
        /// <see cref="UnivariateDistributionFactory.CreateDistribution(XElement)"/> restores both forms.
        /// </remarks>
        public virtual XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            var names = GetParameterPropertyNames;
            var parms = GetParameters;
            for (int i = 0; i < NumberOfParameters; i++)
            {
                result.SetAttributeValue(names[i], parms[i].ToString("G17", CultureInfo.InvariantCulture));
            }
            DistributionNumerics.WriteConfigurationAttributes(this, result);
            return result;
        }

        /// <summary>
        /// The equality operator to allow comparison of continuous distribution parameters. If distributions are equal the operator will return true, otherwise false.
        /// </summary>
        /// <param name="left">Continuous distribution on the left side of the operator.</param>
        /// <param name="right">Continuous distribution on the right side of the operator.</param>
        /// <returns>Boolean indicating if the equality operator results in a true statement.</returns>
        public static bool operator ==(UnivariateDistributionBase left, UnivariateDistributionBase right)
        {
            // Check for null arguments. Keep in mind null == null
            if (left is null && right is null)
            {
                return true;
            }
            else if (left is null)
            {
                return false;
            }
            else if (right is null)
            {
                return false;
            }
            // Check to see if they are the same distribution
            if (left.Type != right.Type)
                return false;
            // Test for equality of parameters
            // Compare parameters here. 
            // Need to have special cases for kernel density, univariate, and bivariate distributions. Everything else can use the GetParameters method.
            if (left.Type == UnivariateDistributionType.Empirical)
            {
                EmpiricalDistribution leftTyped = (EmpiricalDistribution)left;
                EmpiricalDistribution rightTyped = (EmpiricalDistribution)right;
                // Check X values for equality
                var leftValues = leftTyped.XValues;
                var rightValues = rightTyped.XValues;
                if (leftValues.Count != rightValues.Count)
                    return false;
                for (int i = 0; i < leftValues.Count; i++)
                {
                    if (Math.Abs(leftValues[i] - rightValues[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
                // Check P values for equality, no need to test count since it was tested with X values
                leftValues = leftTyped.ProbabilityValues;
                rightValues = rightTyped.ProbabilityValues;
                for (int i = 0; i < leftValues.Count; i++)
                {
                    if (Math.Abs(leftValues[i] - rightValues[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
            }
            else if (left.Type == UnivariateDistributionType.KernelDensity)
            {
                throw new NotImplementedException();
            }
            else
            {
                var leftParams = left.GetParameters;
                var rightParams = right.GetParameters;
                for (int i = 0; i < leftParams.Count(); i++)
                {
                    if (Math.Abs(leftParams[i] - rightParams[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
            }
            // All tests passed
            return true;
            // 
        }

        /// <summary>
        /// The non-equality operator to allow comparison of continuous distribution parameters. If distributions are not equal the operator will return true, otherwise false.
        /// </summary>
        /// <param name="left">Continuous distribution on the left side of the operator.</param>
        /// <param name="right">Continuous distribution on the right side of the operator.</param>
        /// <returns>Boolean indicating if the non-equality operator results in a true statement.</returns>
        public static bool operator !=(UnivariateDistributionBase left, UnivariateDistributionBase right)
        {
            return !(left == right);
        }

        /// <summary>
        /// Determines whether the specified object is equal to the current object.
        /// </summary>
        /// <param name="obj">The object to compare with the current object.</param>
        /// <returns>True if the specified object is equal to the current object; otherwise, False.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is UnivariateDistributionBase other)
            {
                return this == other;
            }
            return false;
        }

        /// <summary>
        /// Serves as the default hash function.
        /// </summary>
        /// <returns>A hash code for the current object.</returns>
        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 23 + Type.GetHashCode();

                if (Type == UnivariateDistributionType.Empirical)
                {
                    var empirical = this as EmpiricalDistribution;
                    if (empirical is not null)
                    {
                        foreach (var x in empirical.XValues)
                        {
                            hash = hash * 23 + x.GetHashCode();
                        }
                        foreach (var p in empirical.ProbabilityValues)
                        {
                            hash = hash * 23 + p.GetHashCode();
                        }
                    }
                }
                else if (Type != UnivariateDistributionType.KernelDensity)
                {
                    var parameters = GetParameters;
                    foreach (var param in parameters)
                    {
                        hash = hash * 23 + param.GetHashCode();
                    }
                }

                return hash;
            }
        }

        /// <summary>
        /// Implementation of IComparable to sort distributions by their Mean.
        /// </summary>
        /// <param name="other">Another UnivariateDistributionBase instance.</param>
        /// <returns>
        /// A negative number if this instance's mean is less than the other; zero if equal; a positive number if greater.
        /// </returns>
        public int CompareTo(UnivariateDistributionBase? other)
        {
            if (other is null)
                return 1; // non-null instance is considered greater than null

            return this.Mean.CompareTo(other.Mean);
        }

    }
}
