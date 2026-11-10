using System;
using System.Collections.Generic;
using Numerics.Data.Statistics;
using Numerics.Mathematics.RootFinding;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Poisson distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Poisson_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Poisson : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, IBootstrappable
    {
 
        /// <summary>
        /// Constructs a Poisson distribution with λ = 1.
        /// </summary>
        public Poisson()
        {
            SetParameters([1.0d]);
        }

        /// <summary>
        /// Constructs a Poisson distribution with a given rate λ (lambda).
        /// </summary>
        /// <param name="rate">The rate parameter λ (lambda). Range: λ > 0.</param>
        public Poisson(double rate)
        {
            SetParameters([rate]);
        }

        private double _lambda;

        /// <summary>
        /// Gets and sets the Poisson's rate parameter λ (lambda). Range: λ > 0.
        /// </summary>
        public double Lambda
        {
            get { return _lambda; }

            set
            {
                _parametersValid = ValidateParameters([value], false) is null;
                _lambda = value;
            }
        }

        /// <inheritdoc/>
        public override int NumberOfParameters
        {
            get { return 1; }
        }

        /// <inheritdoc/>
        public override UnivariateDistributionType Type
        {
            get { return UnivariateDistributionType.Poisson; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Poisson"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "P"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[1, 2];
                parmString[0, 0] = "Rate (λ)";
                parmString[0, 1] = Lambda.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["λ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Lambda)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Lambda]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return Lambda; }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return InverseCDF(0.5); }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Math.Floor(Lambda); }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(Lambda); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 1.0d / Math.Sqrt(Lambda); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3.0d + 1.0d / Lambda; }
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
            get { return [0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The method of moments and maximum likelihood both give λ = mean(x).
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a value that is not a
        /// nonnegative integer, or has no positive count.</exception>
        /// <exception cref="NotImplementedException"><paramref name="estimationMethod"/> is not supported.</exception>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            ValidateSample(sample);
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(ParametersFromMoments([Statistics.Mean(sample)]));
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
            var newDistribution = new Poisson(Lambda);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            return [moments[0]];
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            ValidateParameters(parameters, true);
            var dist = new Poisson(parameters[0]);
            return [dist.Mean, dist.StandardDeviation, dist.Skewness, dist.Kurtosis];
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The initial value is the closed-form maximum likelihood estimate, the sample mean.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a value that is not a
        /// nonnegative integer, or has no positive count.</exception>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = MLE(sample);
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            DistributionNumerics.PositiveParameterBounds(initialVals[0], out lowerVals[0], out upperVals[0]);
            return new Tuple<double[], double[], double[]>(initialVals, lowerVals, upperVals);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The maximum likelihood estimate has the closed form λ = mean(x).
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a value that is not a
        /// nonnegative integer, or has no positive count.</exception>
        public double[] MLE(IList<double> sample)
        {
            ValidateSample(sample);
            return [Statistics.Mean(sample)];
        }

        /// <summary>
        /// Checks that a sample can be used to estimate the rate.
        /// </summary>
        /// <param name="sample">The observed counts.</param>
        /// <exception cref="ArgumentNullException"><paramref name="sample"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException">The sample is empty, contains a value that is not a
        /// nonnegative integer, or has no positive count.</exception>
        /// <remarks>
        /// The observations are counts, so each must be a finite nonnegative integer. A constant positive sample
        /// still determines λ and is accepted; an all-zero sample would give λ = 0 and is rejected.
        /// </remarks>
        private static void ValidateSample(IList<double> sample)
        {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            if (sample.Count < 1) throw new ArgumentOutOfRangeException(nameof(sample), "Insufficient observations to initialize the distribution.");
            bool positive = false;
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]) || sample[i] < 0d || sample[i] != Math.Floor(sample[i]))
                    throw new ArgumentOutOfRangeException(nameof(sample), "Observations must be nonnegative integer counts.");
                positive |= sample[i] > 0d;
            }
            if (!positive) throw new ArgumentOutOfRangeException(nameof(sample), "At least one count must be positive.");
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            Lambda = parameters[0];
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Lambda), "The rate (λ) must be positive.");
                return new ArgumentOutOfRangeException(nameof(Lambda), "The rate (λ) must be positive.");
            }
            return null!;
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The mass is zero away from the nonnegative integers. See <see cref="LogPDF(double)"/> for the evaluation.
        /// </remarks>
        public override double PDF(double k)
        {
            return Math.Exp(LogPDF(k));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// For rates below 1E5 and integer k up to <see cref="int.MaxValue"/> the log mass is −λ + k·ln(λ) − ln(k!), with
        /// ln(k!) from <see cref="Factorial.LogFactorial(int)"/>. Rates of 1E5 and above, and larger k, use the deviance
        /// form of the gamma kernel, ln(λ^k·e^(−λ)/Γ(k + 1)), which avoids cancelling the large terms.
        /// </para>
        /// <para>
        /// The result is negative infinity away from the nonnegative integers and at +∞, and not-a-number for a
        /// not-a-number argument.
        /// </para>
        /// </remarks>
        public override double LogPDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            return LogMass(Lambda, k);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// The CDF at k is the regularized upper incomplete gamma integral Q(⌊k⌋ + 1, λ). For λ below
        /// 100,000 it is evaluated with <see cref="Gamma.UpperIncomplete(double, double)"/>. From λ = 100,000
        /// upward the relative error of that routine grows from about 1e-10 (it rises roughly in proportion to λ),
        /// so the CDF is evaluated from the log upper gamma tail, which uses Temme's uniform asymptotic expansion
        /// (DLMF 8.12) near the mean and keeps a relative error near 1e-14.
        /// </para>
        /// <para>
        /// Non-integer k is floored. The result is 0 below zero, 1 at +∞, and not-a-number for a not-a-number
        /// argument.
        /// </para>
        /// </remarks>
        public override double CDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            return Distribution(Lambda, k);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The survival probability P(X &gt; k) is the regularized lower incomplete gamma integral P(⌊k⌋ + 1, λ),
        /// evaluated directly from the log lower gamma tail so far upper-tail probabilities are not lost to the
        /// subtraction 1 − CDF. The result is 1 below zero, 0 at +∞, and not-a-number for a not-a-number argument.
        /// </remarks>
        public override double CCDF(double k)
        {
            return Math.Exp(LogCCDF(k));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// Evaluated as ln P(⌊k⌋ + 1, λ), which stays finite after the survival probability underflows.
        /// </remarks>
        public override double LogCCDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            if (double.IsNaN(k)) return double.NaN;
            k = Math.Floor(k);
            if (k < Minimum)
                return 0.0d;
            if (k >= Maximum)
                return double.NegativeInfinity;
            return DistributionNumerics.GammaLogCDF(k + 1d, Lambda);
        }

        /// <inheritdoc/>
        /// <remarks>
        /// <para>
        /// Returns the smallest integer k with CDF(k) ≥ p. For λ below 100,000 and p in [1E-300, 1 − 1E-10] the
        /// search starts at the rounded Cornish-Fisher estimate λ + √λ z + (z² − 1) / 6, with z the standard normal
        /// quantile of p, and walks one step at a time, so it typically takes two CDF evaluations. Other
        /// probabilities below that rate walk one step at a time from ⌊λ⌋; both starts give the same k. For larger
        /// rates the search brackets the quantile with steps that start at √λ and double, then bisects on the
        /// integers, so it stays short for any finite rate.
        /// </para>
        /// <para>
        /// A probability outside [0, 1], or not-a-number, throws.
        /// </para>
        /// </remarks>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (!(probability >= 0.0d && probability <= 1.0d))
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([Lambda], true);
            return StartedQuantile(Lambda, probability);
        }

        /// <summary>
        /// The log probability mass for many rates, one rate per count.
        /// </summary>
        /// <param name="rates">The rate λ of each element. Each must be positive and finite.</param>
        /// <param name="k">The count at which each element is evaluated.</param>
        /// <returns>One log mass per element; element i equals <see cref="LogPDF(double)"/> at k[i] of a
        /// Poisson distribution with rate rates[i].</returns>
        /// <remarks>
        /// Each element is evaluated by the same private kernel as the instance member, so each result is bitwise
        /// identical to the instance result, and no distribution object is created per element.
        /// Arrays of 4096 or more elements are evaluated in parallel ranges; each element is independent, so the
        /// results do not depend on the scheduling or the processor count.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="rates"/> or <paramref name="k"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A rate is not positive and finite.</exception>
        public static double[] LogPDF(double[] rates, double[] k)
        {
            ValidateElementRates(rates, k, nameof(k));
            var result = new double[rates.Length];
            DistributionNumerics.ForEachElementRange(result.Length, (start, end) =>
            {
                for (int i = start; i < end; i++)
                    result[i] = LogMass(rates[i], k[i]);
            });
            return result;
        }

        /// <summary>
        /// The cumulative distribution function for many rates, one rate per count.
        /// </summary>
        /// <param name="rates">The rate λ of each element. Each must be positive and finite.</param>
        /// <param name="k">The count at which each element is evaluated.</param>
        /// <returns>One probability P(X ≤ k[i]) per element; element i equals <see cref="CDF(double)"/> at k[i] of a
        /// Poisson distribution with rate rates[i].</returns>
        /// <remarks>
        /// Each element is evaluated by the same private kernel as the instance member, so each result is bitwise
        /// identical to the instance result.
        /// Arrays of 4096 or more elements are evaluated in parallel ranges; each element is independent, so the
        /// results do not depend on the scheduling or the processor count.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="rates"/> or <paramref name="k"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A rate is not positive and finite.</exception>
        public static double[] CDF(double[] rates, double[] k)
        {
            ValidateElementRates(rates, k, nameof(k));
            var result = new double[rates.Length];
            DistributionNumerics.ForEachElementRange(result.Length, (start, end) =>
            {
                for (int i = start; i < end; i++)
                    result[i] = Distribution(rates[i], k[i]);
            });
            return result;
        }

        /// <summary>
        /// The inverse cumulative distribution function for many rates, one rate per probability.
        /// </summary>
        /// <param name="rates">The rate λ of each element. Each must be positive and finite.</param>
        /// <param name="probabilities">The probability of each element, between 0 and 1.</param>
        /// <returns>One quantile per element: the smallest integer k with P(X ≤ k) ≥ probabilities[i] under rate
        /// rates[i], equal to <see cref="InverseCDF(double)"/> of the corresponding instance.</returns>
        /// <remarks>
        /// The elements are found by the same private search as the instance member, which starts from a
        /// Cornish-Fisher estimate of the quantile, so each result is bitwise identical to the instance result. A
        /// probability of zero gives 0 and one gives +∞.
        /// Arrays of 4096 or more elements are evaluated in parallel ranges; each element is independent, so the
        /// results do not depend on the scheduling or the processor count.
        /// </remarks>
        /// <exception cref="ArgumentNullException"><paramref name="rates"/> or <paramref name="probabilities"/> is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A rate is not positive and finite, or a probability is
        /// outside [0, 1] or not-a-number.</exception>
        public static double[] InverseCDF(double[] rates, double[] probabilities)
        {
            ValidateElementRates(rates, probabilities, nameof(probabilities));
            for (int i = 0; i < probabilities.Length; i++)
            {
                if (!(probabilities[i] >= 0.0d && probabilities[i] <= 1.0d))
                    throw new ArgumentOutOfRangeException(nameof(probabilities), "Probability must be between 0 and 1.");
            }
            var result = new double[rates.Length];
            DistributionNumerics.ForEachElementRange(result.Length, (start, end) =>
            {
                for (int i = start; i < end; i++)
                    result[i] = StartedQuantile(rates[i], probabilities[i]);
            });
            return result;
        }

        /// <summary>
        /// Checks the per-element arguments of the static array members.
        /// </summary>
        /// <param name="rates">The rates, each of which must be positive and finite.</param>
        /// <param name="values">The second argument array, which must match <paramref name="rates"/> in length.</param>
        /// <param name="valuesName">The parameter name of <paramref name="values"/>.</param>
        /// <exception cref="ArgumentNullException">An array is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A rate is not positive and finite.</exception>
        private static void ValidateElementRates(double[] rates, double[] values, string valuesName)
        {
            if (rates == null) throw new ArgumentNullException(nameof(rates));
            if (values == null) throw new ArgumentNullException(valuesName);
            if (rates.Length != values.Length)
                throw new ArgumentException("The rates and " + valuesName + " must have the same length.", valuesName);
            for (int i = 0; i < rates.Length; i++)
            {
                if (!Tools.IsFinite(rates[i]) || rates[i] <= 0.0d)
                    throw new ArgumentOutOfRangeException(nameof(rates), "The rate (λ) must be positive.");
            }
        }

        /// <summary>
        /// The log probability mass at k for a validated rate.
        /// </summary>
        /// <param name="rate">A positive, finite rate λ.</param>
        /// <param name="k">The count.</param>
        /// <returns>The log mass; see <see cref="LogPDF(double)"/> for the evaluation and the special values.</returns>
        private static double LogMass(double rate, double k)
        {
            if (double.IsNaN(k)) return double.NaN;
            if (k < 0.0d || double.IsPositiveInfinity(k)) return double.NegativeInfinity;
            // The probability mass is zero away from the integers.
            if (k != Math.Floor(k)) return double.NegativeInfinity;
            if (k > int.MaxValue || rate >= LargeRate)
                return DistributionNumerics.GammaLogDensity(k + 1d, rate);
            return -rate + k * Math.Log(rate) - Factorial.LogFactorial((int)k);
        }

        /// <summary>
        /// The cumulative probability at k for a validated rate.
        /// </summary>
        /// <param name="rate">A positive, finite rate λ.</param>
        /// <param name="k">The count, floored when not an integer.</param>
        /// <returns>P(X ≤ ⌊k⌋); see <see cref="CDF(double)"/> for the evaluation and the special values.</returns>
        private static double Distribution(double rate, double k)
        {
            if (double.IsNaN(k)) return double.NaN;
            k = Math.Floor(k);
            if (k < 0.0d)
                return 0.0d;
            if (double.IsPositiveInfinity(k))
                return 1.0d;
            if (rate >= LargeRate)
                return Math.Exp(DistributionNumerics.GammaLogSurvival(k + 1d, rate));
            return Gamma.UpperIncomplete(k + 1d, rate);
        }

        /// <summary>
        /// The quantile at a probability in [0, 1] for a validated rate.
        /// </summary>
        /// <param name="rate">A positive, finite rate λ.</param>
        /// <param name="probability">A probability between 0 and 1.</param>
        /// <returns>The smallest integer k with P(X ≤ k) ≥ <paramref name="probability"/>; see
        /// <see cref="InverseCDF(double)"/> for the search.</returns>
        private static double Quantile(double rate, double probability)
        {
            if (probability == 0.0d)
                return 0.0d;
            if (probability == 1.0d)
                return double.PositiveInfinity;
            if (rate >= LargeRate)
                return LargeRateQuantile(rate, probability);

            int k = (int)Math.Max(0, Math.Floor(rate));  // start near the mean

            // Move downward if needed
            while (k > 0 && Distribution(rate, k - 1) >= probability)
                k--;

            // Move upward to find first k such that CDF(k) >= probability
            while (Distribution(rate, k) < probability)
                k++;

            return k;
        }

        /// <summary>
        /// The quantile of <see cref="Quantile"/>, searched from a Cornish-Fisher start.
        /// </summary>
        /// <param name="rate">A positive, finite rate λ.</param>
        /// <param name="probability">A probability between 0 and 1.</param>
        /// <returns>The smallest integer k with P(X ≤ k) ≥ <paramref name="probability"/>.</returns>
        /// <remarks>
        /// <para>
        /// For rates below 1E5 and probabilities in [1E-300, 1 − 1E-10], the search starts at the rounded two-term
        /// Cornish-Fisher expansion λ + √λ z + (z² − 1) / 6, with z the standard normal quantile of p (the start of
        /// R's qpois), and then steps one integer at a time with <see cref="Distribution"/>, as <see cref="Quantile"/>
        /// does from ⌊λ⌋. The start is usually the quantile itself or next to it, so two CDF evaluations typically
        /// suffice.
        /// </para>
        /// <para>
        /// The smallest k with CDF(k) ≥ p is independent of the start wherever the computed CDF does not decrease
        /// between the two starts. In that probability range the mass near the quantile is at least about 1E-12 of
        /// the CDF, far above its rounding, so the result equals that of <see cref="Quantile"/>. Outside the range,
        /// where the mass can fall to the rounding of the CDF, and for rates of 1E5 or more, <see cref="Quantile"/> is
        /// used directly.
        /// </para>
        /// </remarks>
        private static double StartedQuantile(double rate, double probability)
        {
            if (rate >= LargeRate || !(probability >= MinimumStartedProbability && probability <= MaximumStartedProbability))
                return Quantile(rate, probability);
            double z = Normal.StandardZ(probability);
            double start = Math.Floor(rate + Math.Sqrt(rate) * z + (z * z - 1d) / 6d + 0.5d);
            int k = start > 0d ? (int)start : 0;
            if (Distribution(rate, k) >= probability)
            {
                while (k > 0 && Distribution(rate, k - 1) >= probability)
                    k--;
                return k;
            }
            do
                k++;
            while (Distribution(rate, k) < probability);
            return k;
        }

        /// <summary>The smallest probability for which <see cref="StartedQuantile"/> uses its Cornish-Fisher start.</summary>
        private const double MinimumStartedProbability = 1E-300;

        /// <summary>The largest probability for which <see cref="StartedQuantile"/> uses its Cornish-Fisher start.</summary>
        private const double MaximumStartedProbability = 1d - 1E-10;

        /// <summary>
        /// The rate from which the mass uses the gamma deviance form, the CDF uses the log gamma tail and the quantile
        /// search uses bisection.
        /// </summary>
        private const double LargeRate = 1E5;

        /// <summary>
        /// Finds the smallest integer k with CDF(k) ≥ p by bracketing and bisection.
        /// </summary>
        /// <param name="rate">A rate λ of at least 1E5.</param>
        /// <param name="probability">A probability strictly between 0 and 1.</param>
        /// <returns>The smallest integer k with CDF(k) ≥ <paramref name="probability"/>.</returns>
        /// <remarks>
        /// The search index is a double, so rates beyond the range of <see cref="int"/> are handled. The bracket
        /// keeps CDF(lower) &lt; p ≤ CDF(upper), starting from ⌊λ⌋ with a step of √λ that doubles until the
        /// quantile is enclosed; bisection then narrows it to adjacent integers, or to adjacent representable
        /// values where the spacing of doubles exceeds one.
        /// </remarks>
        private static double LargeRateQuantile(double rate, double probability)
        {
            double start = Math.Floor(rate);
            double step = Math.Ceiling(Math.Sqrt(rate));
            double lower, upper;
            if (Distribution(rate, start) >= probability)
            {
                upper = start;
                lower = upper - step;
                while (lower >= 0d && Distribution(rate, lower) >= probability)
                {
                    upper = lower;
                    step *= 2d;
                    lower = upper - step;
                }
                if (lower < 0d) lower = -1d;
            }
            else
            {
                lower = start;
                upper = lower + step;
                while (Distribution(rate, upper) < probability)
                {
                    lower = upper;
                    step *= 2d;
                    upper = lower + step;
                }
            }
            while (upper - lower > 1d)
            {
                double middle = Math.Floor(lower + (upper - lower) / 2d);
                if (middle <= lower || middle >= upper) break;
                if (Distribution(rate, middle) >= probability)
                    upper = middle;
                else
                    lower = middle;
            }
            return upper;
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Poisson(Lambda);
        }
    
    }
}