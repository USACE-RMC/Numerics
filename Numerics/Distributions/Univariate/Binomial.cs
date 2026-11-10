using System;
using System.Collections.Generic;
using Numerics.Mathematics.SpecialFunctions;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Binomial distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Binomial_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class Binomial : UnivariateDistributionBase
    {

        
        /// <summary>
        /// Constructs a Binomial distribution with p=0.5 and n=10.
        /// </summary>
        public Binomial()
        {
            SetParameters([0.5d, 10d]);
        }

        /// <summary>
        /// Constructs a Binomial distribution with a given probability and sample size.
        /// </summary>
        /// <param name="probability">The success probability (p) in each trial. Range: 0 ≤ p ≤ 1.</param>
        /// <param name="numberOfTrials">The number of trials (n). Range: n ≥ 0.</param>
        public Binomial(double probability, int numberOfTrials)
        {
            SetParameters([probability, numberOfTrials]);
        }
     
        private double _probabilityOfSuccess;
        private int _numberOfTrials;

        /// <summary>
        /// Gets and sets the success probability in each trial. Range: 0 ≤ p ≤ 1.
        /// </summary>
        public double ProbabilityOfSuccess
        {
            get { return _probabilityOfSuccess; }
            set
            {
                _parametersValid = ValidateParameters(new[] { value, NumberOfTrials }, false) is null;
                _probabilityOfSuccess = value;
            }
        }

        /// <summary>
        /// Gets the complement of the probability.
        /// </summary>
        public double Complement
        {
            get { return 1d - ProbabilityOfSuccess; }
        }

        /// <summary>
        /// Gets and sets the number of trials. Range: n ≥ 0.
        /// </summary>
        public int NumberOfTrials
        {
            get { return _numberOfTrials; }
            set
            {
                _parametersValid = ValidateParameters(new[] { ProbabilityOfSuccess, value }, false) is null;
                _numberOfTrials = value;
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
            get { return UnivariateDistributionType.Binomial; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Binomial"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Bin"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Probability of Success (p)";
                parmString[1, 0] = "Number of Trials (n)";
                parmString[0, 1] = ProbabilityOfSuccess.ToString();
                parmString[1, 1] = NumberOfTrials.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["p", "n"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(ProbabilityOfSuccess), nameof(NumberOfTrials)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [ProbabilityOfSuccess, NumberOfTrials]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get { return NumberOfTrials * ProbabilityOfSuccess; }
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The median is the smallest k with CDF(k) ≥ 1/2. Every median lies in [floor(np), ceil(np)], and it equals np
        /// when np is an integer (Kaas and Buhrman, 1980), so the result is floor(np) when CDF(floor(np)) ≥ 1/2 and
        /// ceil(np) otherwise. For p = 1/2 and odd n the distribution is symmetric about n/2, so CDF((n − 1)/2) is exactly
        /// 1/2 and the median is (n − 1)/2; this case is taken from the symmetry rather than from the rounded CDF.
        /// </remarks>
        public override double Median
        {
            get
            {
                // Validate parameters
                if (_parametersValid == false)
                    ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
                if (ProbabilityOfSuccess == 0.5d && NumberOfTrials % 2 == 1)
                    return (NumberOfTrials - 1) / 2;
                double np = NumberOfTrials * ProbabilityOfSuccess;
                double lower = Math.Floor(np);
                return CDF(lower) >= 0.5d ? lower : Math.Ceiling(np);
            }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get
            {
                if (ProbabilityOfSuccess == 1.0d)
                    return NumberOfTrials;
                if (ProbabilityOfSuccess == 0.0d)
                    return 0.0d;
                return (int)Math.Floor((NumberOfTrials + 1) * ProbabilityOfSuccess);
            }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Math.Sqrt(NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return (1.0d - 2.0d * ProbabilityOfSuccess) / Math.Sqrt(NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d + (1.0d - 6d * Complement * ProbabilityOfSuccess) / (NumberOfTrials * ProbabilityOfSuccess * Complement); }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return NumberOfTrials; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [0.0d, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [1.0d, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            _probabilityOfSuccess = parameters[0];
            _numberOfTrials = (int)parameters[1];
            _parametersValid = ValidateParameters(parameters, false) is null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException? ValidateParameters(IList<double> parameters, bool throwException)
        {
            // Validate probability
            if (double.IsNaN(parameters[0]) || double.IsInfinity(parameters[0]) || parameters[0] < 0.0d || parameters[0] > 1.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "Probability must be between 0 and 1.");
                return new ArgumentOutOfRangeException(nameof(ProbabilityOfSuccess), "Probability must be between 0 and 1.");
            }
            if (double.IsNaN(parameters[1]) || double.IsInfinity(parameters[1]) || parameters[1] <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(NumberOfTrials), "The number of trials (n) must be positive.");
                return new ArgumentOutOfRangeException(nameof(NumberOfTrials), "The number of trials (n) must be positive.");
            }
            return null!;            
        }

        /// <inheritdoc/>
        public override double PDF(double k)
        {
            // Exponentiate the log-mass so the result stays finite when the binomial coefficient overflows.
            return Math.Exp(LogPDF(k));
        }

        /// <inheritdoc/>
        /// <remarks>
        /// The mass is zero, and the log mass negative infinity, away from the integers.
        /// </remarks>
        public override double LogPDF(double k)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
            return LogMass(ProbabilityOfSuccess, NumberOfTrials, k);
        }

        /// <inheritdoc/>
        public override double CDF(double k)
        {
            if (_parametersValid == false)
                ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
            return Distribution(ProbabilityOfSuccess, NumberOfTrials, k);
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
                ValidateParameters([ProbabilityOfSuccess, NumberOfTrials], true);
            double k = 0d;
            for (int i = 0; i <= NumberOfTrials; i++)
            {
                if (Distribution(ProbabilityOfSuccess, NumberOfTrials, i) >= probability)
                {
                    k = i;
                    break;
                }
            }
            return k;
        }

        /// <summary>
        /// The log probability mass for many parameter sets, one set per count.
        /// </summary>
        /// <param name="successProbabilities">The success probability p of each element, between 0 and 1.</param>
        /// <param name="trials">The number of trials n of each element. Each must be positive.</param>
        /// <param name="k">The number of successes at which each element is evaluated.</param>
        /// <returns>One log mass per element; element i equals <see cref="LogPDF(double)"/> at k[i] of a binomial
        /// distribution with success probability successProbabilities[i] and trials[i] trials.</returns>
        /// <remarks>
        /// Each element is evaluated by the same private kernel as the instance member, so each result is bitwise
        /// identical to the instance result, and no distribution object is created per element.
        /// Arrays of 4096 or more elements are evaluated in parallel ranges; each element is independent, so the
        /// results do not depend on the scheduling or the processor count.
        /// </remarks>
        /// <exception cref="ArgumentNullException">An array is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A success probability is outside [0, 1] or not-a-number, or a
        /// number of trials is not positive.</exception>
        public static double[] LogPDF(double[] successProbabilities, int[] trials, double[] k)
        {
            ValidateElementParameters(successProbabilities, trials, k);
            var result = new double[k.Length];
            DistributionNumerics.ForEachElementRange(result.Length, (start, end) =>
            {
                for (int i = start; i < end; i++)
                    result[i] = LogMass(successProbabilities[i], trials[i], k[i]);
            });
            return result;
        }

        /// <summary>
        /// The cumulative distribution function for many parameter sets, one set per count.
        /// </summary>
        /// <param name="successProbabilities">The success probability p of each element, between 0 and 1.</param>
        /// <param name="trials">The number of trials n of each element. Each must be positive.</param>
        /// <param name="k">The number of successes at which each element is evaluated.</param>
        /// <returns>One probability P(X ≤ k[i]) per element; element i equals <see cref="CDF(double)"/> at k[i] of the
        /// corresponding binomial distribution.</returns>
        /// <remarks>
        /// Each element is evaluated by the same private kernel as the instance member, so each result is bitwise
        /// identical to the instance result.
        /// Arrays of 4096 or more elements are evaluated in parallel ranges; each element is independent, so the
        /// results do not depend on the scheduling or the processor count.
        /// </remarks>
        /// <exception cref="ArgumentNullException">An array is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A success probability is outside [0, 1] or not-a-number, or a
        /// number of trials is not positive.</exception>
        public static double[] CDF(double[] successProbabilities, int[] trials, double[] k)
        {
            ValidateElementParameters(successProbabilities, trials, k);
            var result = new double[k.Length];
            DistributionNumerics.ForEachElementRange(result.Length, (start, end) =>
            {
                for (int i = start; i < end; i++)
                    result[i] = Distribution(successProbabilities[i], trials[i], k[i]);
            });
            return result;
        }

        /// <summary>
        /// Checks the per-element arguments of the static array members.
        /// </summary>
        /// <param name="successProbabilities">The success probabilities, each between 0 and 1.</param>
        /// <param name="trials">The numbers of trials, each positive.</param>
        /// <param name="k">The evaluation points, which must match the parameter arrays in length.</param>
        /// <exception cref="ArgumentNullException">An array is <see langword="null"/>.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">A parameter is invalid.</exception>
        private static void ValidateElementParameters(double[] successProbabilities, int[] trials, double[] k)
        {
            if (successProbabilities == null) throw new ArgumentNullException(nameof(successProbabilities));
            if (trials == null) throw new ArgumentNullException(nameof(trials));
            if (k == null) throw new ArgumentNullException(nameof(k));
            if (trials.Length != successProbabilities.Length || k.Length != successProbabilities.Length)
                throw new ArgumentException("The success probabilities, trials and k must have the same length.", nameof(k));
            for (int i = 0; i < k.Length; i++)
            {
                double p = successProbabilities[i];
                if (!(p >= 0.0d && p <= 1.0d))
                    throw new ArgumentOutOfRangeException(nameof(successProbabilities), "Probability must be between 0 and 1.");
                if (trials[i] <= 0)
                    throw new ArgumentOutOfRangeException(nameof(trials), "The number of trials (n) must be positive.");
            }
        }

        /// <summary>
        /// The log probability mass at k for validated parameters.
        /// </summary>
        /// <param name="probability">The success probability p, between 0 and 1.</param>
        /// <param name="trials">The positive number of trials n.</param>
        /// <param name="k">The number of successes.</param>
        /// <returns>ln C(n, k) + k ln(p) + (n − k) ln(1 − p), with the exact limits at p = 0 and p = 1 and negative
        /// infinity away from the integers in [0, n].</returns>
        private static double LogMass(double probability, int trials, double k)
        {
            // The probability mass is zero away from the integers.
            if (k != Math.Floor(k) && !double.IsNaN(k)) return double.NegativeInfinity;
            if (k < 0.0d || k > trials) return double.NegativeInfinity;
            if (probability == 0.0d)
                return k == 0d ? 0.0d : double.NegativeInfinity;
            if (probability == 1.0d)
                return k == trials ? 0.0d : double.NegativeInfinity;
            // ln C(n, k) + k ln(p) + (n - k) ln(1 - p)
            int i = (int)k;
            double logChoose = Factorial.LogFactorial(trials) - Factorial.LogFactorial(i) - Factorial.LogFactorial(trials - i);
            return logChoose + k * Math.Log(probability) + (trials - k) * Math.Log(1d - probability);
        }

        /// <summary>
        /// The cumulative probability at k for validated parameters.
        /// </summary>
        /// <param name="probability">The success probability p, between 0 and 1.</param>
        /// <param name="trials">The positive number of trials n.</param>
        /// <param name="k">The number of successes, floored when not an integer.</param>
        /// <returns>P(X ≤ ⌊k⌋), the regularized incomplete beta function I(1 − p; n − ⌊k⌋, ⌊k⌋ + 1).</returns>
        private static double Distribution(double probability, int trials, double k)
        {
            k = Math.Floor(k);
            if (k < 0.0d)
                return 0.0d;
            if (k >= trials)
                return 1.0d;
            return Beta.Incomplete(trials - k, k + 1d, 1d - probability);
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Binomial(ProbabilityOfSuccess, NumberOfTrials);
        }
     
    }
}