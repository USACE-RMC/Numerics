using System;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Multinomial distribution, a generalization of the binomial distribution to K categories.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The multinomial distribution models the number of outcomes in K categories from N independent trials,
    /// where each trial results in exactly one of the K categories with fixed probabilities p₁, p₂, ..., pₖ.
    /// The probability mass function gives the probability of observing x₁ in category 1, x₂ in category 2,
    /// etc., where Σxᵢ = N.
    /// </para>
    /// <para>
    /// Key applications include:
    /// <list type="bullet">
    /// <item><description>MCMC sampling: selecting trajectory states in NUTS (No-U-Turn Sampler) via weighted multinomial draws.</description></item>
    /// <item><description>Categorical data modeling: modeling counts from multiple categories in LifeSim.</description></item>
    /// <item><description>Bayesian inference: likelihood function for categorical observations.</description></item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Johnson, N.L., Kotz, S. and Balakrishnan, N. (1997). "Discrete Multivariate Distributions." Wiley.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Multinomial_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Multinomial : MultivariateDistribution
    {

        /// <summary>
        /// Private parameterless constructor for use by Clone().
        /// </summary>
        private Multinomial() { }

        /// <summary>
        /// Constructs a multinomial distribution with the specified number of trials and category probabilities.
        /// </summary>
        /// <param name="numberOfTrials">The number of trials N. Must be positive.</param>
        /// <param name="probabilities">The probability vector p₁, ..., pₖ. Each must be in [0, 1] and they must sum to 1.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when N &lt; 1, probabilities have fewer than 2 elements,
        /// any probability is negative, or probabilities don't sum to 1.</exception>
        public Multinomial(int numberOfTrials, double[] probabilities)
        {
            if (numberOfTrials < 1)
                throw new ArgumentOutOfRangeException(nameof(numberOfTrials), "The number of trials must be positive.");
            if (probabilities == null || probabilities.Length < 2)
                throw new ArgumentOutOfRangeException(nameof(probabilities), "The probability vector must have at least 2 elements.");

            double sum = 0;
            for (int i = 0; i < probabilities.Length; i++)
            {
                if (double.IsNaN(probabilities[i]) || probabilities[i] < 0 || probabilities[i] > 1)
                    throw new ArgumentOutOfRangeException(nameof(probabilities), $"Probability p[{i}] must be in [0, 1].");
                sum += probabilities[i];
            }
            if (Math.Abs(sum - 1.0) > 1e-10)
                throw new ArgumentOutOfRangeException(nameof(probabilities), "Probabilities must sum to 1.");

            _n = numberOfTrials;
            _p = (double[])probabilities.Clone();
            _dimension = probabilities.Length;
        }

        private int _n;
        private double[] _p = null!;
        private int _dimension;

        /// <summary>
        /// Gets the number of trials N.
        /// </summary>
        public int NumberOfTrials
        {
            get { return _n; }
        }

        /// <summary>
        /// Gets the probability vector p.
        /// </summary>
        public double[] Probabilities
        {
            get { return (double[])_p.Clone(); }
        }

        /// <inheritdoc/>
        public override int Dimension
        {
            get { return _dimension; }
        }

        /// <inheritdoc/>
        public override MultivariateDistributionType Type
        {
            get { return MultivariateDistributionType.Multinomial; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Multinomial"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Mult"; }
        }

        /// <inheritdoc/>
        public override bool ParametersValid
        {
            get
            {
                if (_p == null || _p.Length < 2 || _n < 1) return false;
                double sum = 0;
                for (int i = 0; i < _p.Length; i++)
                {
                    if (double.IsNaN(_p[i]) || _p[i] < 0 || _p[i] > 1) return false;
                    sum += _p[i];
                }
                return Math.Abs(sum - 1.0) <= 1e-10;
            }
        }

        /// <summary>
        /// Gets the mean vector. Mean[i] = N * p[i].
        /// </summary>
        public double[] Mean
        {
            get
            {
                var mean = new double[_dimension];
                for (int i = 0; i < _dimension; i++)
                    mean[i] = _n * _p[i];
                return mean;
            }
        }

        /// <summary>
        /// Gets the marginal variance vector. Var[i] = N * p[i] * (1 - p[i]).
        /// </summary>
        public double[] Variance
        {
            get
            {
                var variance = new double[_dimension];
                for (int i = 0; i < _dimension; i++)
                    variance[i] = _n * _p[i] * (1.0 - _p[i]);
                return variance;
            }
        }

        /// <summary>
        /// Gets the covariance between components i and j: Cov(Xᵢ, Xⱼ) = -N * pᵢ * pⱼ.
        /// </summary>
        /// <param name="i">The first component index (0-based).</param>
        /// <param name="j">The second component index (0-based).</param>
        /// <returns>The covariance. Negative for i != j.</returns>
        public double Covariance(int i, int j)
        {
            if (i < 0 || i >= _dimension || j < 0 || j >= _dimension)
                throw new ArgumentOutOfRangeException("Index out of range.");
            if (i == j) return _n * _p[i] * (1.0 - _p[i]);
            return -_n * _p[i] * _p[j];
        }

        /// <summary>
        /// Computes the probability mass function (PMF) for a given count vector x.
        /// </summary>
        /// <param name="x">The count vector. Each xᵢ must be a non-negative integer and Σxᵢ = N.</param>
        /// <returns>The probability P(X = x). Returns 0 if x is not a valid count vector.</returns>
        /// <remarks>
        /// <para>
        /// The PMF is: P(X = x) = N! / (x₁! x₂! ... xₖ!) * p₁^x₁ * p₂^x₂ * ... * pₖ^xₖ
        /// Computed in log-space for numerical stability.
        /// </para>
        /// </remarks>
        public override double PDF(double[] x)
        {
            return Math.Exp(LogPMF(x));
        }

        /// <summary>
        /// Computes the log of the probability mass function.
        /// </summary>
        /// <param name="x">The count vector. Each xᵢ must be a non-negative integer and Σxᵢ = N.</param>
        /// <returns>The log probability. Returns <see cref="double.NegativeInfinity"/> if x is not valid.</returns>
        public double LogPMF(double[] x)
        {
            if (x == null || x.Length != _dimension) return double.NegativeInfinity;

            // Validate: all non-negative integers that sum to N
            int sum = 0;
            for (int i = 0; i < _dimension; i++)
            {
                int xi = (int)Math.Round(x[i]);
                if (xi < 0 || Math.Abs(x[i] - xi) > 1e-10) return double.NegativeInfinity;
                sum += xi;
            }
            if (sum != _n) return double.NegativeInfinity;

            // Compute in log-space: log(N!) - sum(log(xi!)) + sum(xi * log(pi))
            double logPmf = Factorial.LogFactorial(_n);
            for (int i = 0; i < _dimension; i++)
            {
                int xi = (int)Math.Round(x[i]);
                logPmf -= Factorial.LogFactorial(xi);
                if (xi > 0)
                {
                    if (_p[i] <= 0) return double.NegativeInfinity;
                    logPmf += xi * Math.Log(_p[i]);
                }
            }
            return logPmf;
        }

        /// <inheritdoc/>
        public override double LogPDF(double[] x)
        {
            return LogPMF(x);
        }

        /// <summary>
        /// The CDF of the multinomial distribution is not available in closed form.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        /// <returns>This method always throws <see cref="NotImplementedException"/>.</returns>
        /// <exception cref="NotImplementedException">Always thrown.</exception>
        public override double CDF(double[] x)
        {
            throw new NotImplementedException("The CDF of the multinomial distribution does not have a closed-form expression.");
        }

        /// <inheritdoc/>
        public override MultivariateDistribution Clone()
        {
            var clone = new Multinomial();
            clone._n = _n;
            clone._p = (double[])_p.Clone();
            clone._dimension = _dimension;
            return clone;
        }

        /// <summary>
        /// Generates random samples from the multinomial distribution using sequential binomial sampling.
        /// </summary>
        /// <param name="sampleSize">The number of samples to generate.</param>
        /// <param name="seed">Optional seed for reproducibility. Use -1 for a random seed.</param>
        /// <returns>A 2D array of shape [sampleSize, Dimension]. Each row sums to N.</returns>
        /// <remarks>
        /// <para>
        /// For each sample, draws are made sequentially: for category i, draw from
        /// Binomial(n_remaining, p_i / p_remaining). This is exact and efficient.
        /// </para>
        /// </remarks>
        public double[,] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            var rng = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize, _dimension];

            for (int s = 0; s < sampleSize; s++)
            {
                int nRemaining = _n;
                double pRemaining = 1.0;

                for (int i = 0; i < _dimension - 1; i++)
                {
                    if (nRemaining == 0 || pRemaining <= 0)
                    {
                        sample[s, i] = 0;
                        continue;
                    }

                    double conditionalP = _p[i] / pRemaining;
                    if (conditionalP >= 1.0)
                    {
                        sample[s, i] = nRemaining;
                        nRemaining = 0;
                    }
                    else
                    {
                        // Draw from Binomial(nRemaining, conditionalP)
                        int count = BinomialSample(nRemaining, conditionalP, rng);
                        sample[s, i] = count;
                        nRemaining -= count;
                    }
                    pRemaining -= _p[i];
                }

                // Last category gets the remainder
                sample[s, _dimension - 1] = nRemaining;
            }

            return sample;
        }

        /// <summary>
        /// Samples a single index from a discrete set of categories weighted by the given probabilities.
        /// </summary>
        /// <param name="weights">Unnormalized weights. Must be non-negative with at least one positive value.</param>
        /// <param name="rng">The random number generator.</param>
        /// <returns>The 0-based index of the selected category.</returns>
        /// <remarks>
        /// <para>
        /// This is a weighted categorical sampling method used by the NUTS algorithm
        /// to select a trajectory state weighted by exp(H).
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentException">Thrown if all weights are zero or the array is empty.</exception>
        public static int Sample(double[] weights, Random rng)
        {
            if (weights == null || weights.Length == 0)
                throw new ArgumentException("Weights array must be non-empty.", nameof(weights));

            double totalWeight = 0;
            for (int i = 0; i < weights.Length; i++)
            {
                if (weights[i] < 0) throw new ArgumentException("Weights must be non-negative.", nameof(weights));
                totalWeight += weights[i];
            }
            if (totalWeight <= 0)
                throw new ArgumentException("At least one weight must be positive.", nameof(weights));

            double u = rng.NextDouble() * totalWeight;
            double cumulative = 0;
            for (int i = 0; i < weights.Length; i++)
            {
                cumulative += weights[i];
                if (u <= cumulative)
                    return i;
            }

            // Should not reach here, but return last index as a safeguard
            return weights.Length - 1;
        }

        /// <summary>
        /// Maps an array of uniform variates to category indices through a cumulative weight table:
        /// each index is the first k with <c>cumulativeWeights[k] &gt; u</c>, clipped to the last category.
        /// </summary>
        /// <param name="cumulativeWeights">The cumulative category weights c[0] ≤ c[1] ≤ ... ≤ c[n-1].
        /// Must be nonempty, finite, and non-decreasing. Normalization is not required.</param>
        /// <param name="uniforms">The variates to map, on the scale of <paramref name="cumulativeWeights"/>
        /// (for a table normalized to end at one, uniforms on [0, 1)). Must not contain NaN.</param>
        /// <returns>One 0-based category index in [0, n-1] per variate, in input order.</returns>
        /// <remarks>
        /// <para>
        /// The semantics equal numpy <c>clip(searchsorted(cumulativeWeights, uniforms, side="right"), 0, n - 1)</c>.
        /// A variate exactly on a boundary c[k] selects category k + 1 (strictly greater), every variate
        /// in [c[0], c[n-1]) skips a zero-weight category (c[k] = c[k-1]), a variate below c[0] selects
        /// category 0, and a variate at or above c[n-1] (including +∞) selects the last category.
        /// </para>
        /// <para>
        /// Each variate is located by binary search over the first n - 1 entries, so the cost is
        /// O(m log n) for m variates, with one O(n) check of the table. The whole array is mapped in one
        /// call, so a caller crossing a language boundary pays the crossing once. Unlike
        /// <see cref="Sample(double[], Random)"/>, which accepts a variate on a boundary into the lower
        /// category, this method follows the right-continuous inverse of the cumulative table.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Devroye, L. (1986). Non-Uniform Random Variate Generation, Chapter III. Springer-Verlag, New York.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either array is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the table is empty or decreases.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a table entry is not finite or a variate is NaN.</exception>
        public static int[] CategoricalIndices(double[] cumulativeWeights, double[] uniforms)
        {
            return CategoricalIndices(cumulativeWeights, uniforms, false);
        }

        /// <summary>
        /// Maps an array of uniform variates to category indices through a cumulative weight table, selecting either
        /// the first k with <c>cumulativeWeights[k] &gt; u</c> or the first k with <c>cumulativeWeights[k] ≥ u</c>,
        /// clipped to the last category.
        /// </summary>
        /// <param name="cumulativeWeights">The cumulative category weights c[0] ≤ c[1] ≤ ... ≤ c[n-1].
        /// Must be nonempty, finite, and non-decreasing. Normalization is not required.</param>
        /// <param name="uniforms">The variates to map, on the scale of <paramref name="cumulativeWeights"/>.
        /// Must not contain NaN.</param>
        /// <param name="inclusive">If false, each index is min{k : c[k] &gt; u}, as in
        /// <see cref="CategoricalIndices(double[], double[])"/>. If true, each index is min{k : c[k] ≥ u},
        /// the generalized inverse F⁻¹(u) = min{k : F(k) ≥ u} of a discrete CDF table F.</param>
        /// <returns>One 0-based category index in [0, n-1] per variate, in input order.</returns>
        /// <remarks>
        /// <para>
        /// With <paramref name="inclusive"/> true the semantics equal numpy
        /// <c>clip(searchsorted(cumulativeWeights, uniforms, side="left"), 0, n - 1)</c>. A variate exactly on a
        /// boundary c[k] selects category k (the lowest such k when entries repeat), a variate at or below c[0]
        /// selects category 0 (even when c[0] is zero), and a variate above c[n-1] (including +∞) selects the last
        /// category. The two forms differ only for variates equal to a table entry.
        /// </para>
        /// <para>
        /// This is the inverse-transform rule for a discrete distribution: for a table of CDF values F(0), F(1), ...
        /// of a count distribution such as the Poisson, the index is the smallest count whose CDF reaches u.
        /// Each variate is located by binary search, so the cost is O(m log n) for m variates, with one O(n)
        /// check of the table.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Devroye, L. (1986). Non-Uniform Random Variate Generation, Chapter III. Springer-Verlag, New York.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either array is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the table is empty or decreases.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a table entry is not finite or a variate is NaN.</exception>
        public static int[] CategoricalIndices(double[] cumulativeWeights, double[] uniforms, bool inclusive)
        {
            if (cumulativeWeights == null) throw new ArgumentNullException(nameof(cumulativeWeights));
            if (uniforms == null) throw new ArgumentNullException(nameof(uniforms));
            int n = cumulativeWeights.Length;
            if (n == 0) throw new ArgumentException("The cumulative weights must not be empty.", nameof(cumulativeWeights));
            for (int k = 0; k < n; k++)
            {
                if (!Tools.IsFinite(cumulativeWeights[k]))
                    throw new ArgumentOutOfRangeException(nameof(cumulativeWeights), "The cumulative weights must be finite.");
                if (k > 0 && cumulativeWeights[k] < cumulativeWeights[k - 1])
                    throw new ArgumentException("The cumulative weights must be non-decreasing.", nameof(cumulativeWeights));
            }

            var indices = new int[uniforms.Length];
            int last = n - 1;
            for (int i = 0; i < uniforms.Length; i++)
            {
                double u = uniforms[i];
                if (double.IsNaN(u))
                    throw new ArgumentOutOfRangeException(nameof(uniforms), "The uniform variates must not be NaN.");
                // First k in [0, n - 1) with c[k] > u (or c[k] >= u when inclusive), else the last category.
                int low = 0, high = last;
                if (inclusive)
                {
                    while (low < high)
                    {
                        int middle = (low + high) >> 1;
                        if (cumulativeWeights[middle] >= u) high = middle;
                        else low = middle + 1;
                    }
                }
                else
                {
                    while (low < high)
                    {
                        int middle = (low + high) >> 1;
                        if (cumulativeWeights[middle] > u) high = middle;
                        else low = middle + 1;
                    }
                }
                indices[i] = low;
            }
            return indices;
        }

        /// <summary>
        /// Samples from a Binomial(n, p) distribution using the inverse CDF method for small n
        /// and the BTPE algorithm for large n.
        /// </summary>
        private static int BinomialSample(int n, double p, MersenneTwister rng)
        {
            if (n <= 0 || p <= 0) return 0;
            if (p >= 1.0) return n;

            // For small n, use direct simulation
            if (n < 25)
            {
                int count = 0;
                for (int i = 0; i < n; i++)
                {
                    if (rng.NextDouble() < p)
                        count++;
                }
                return count;
            }

            // For larger n, use the inverse CDF via the existing Binomial distribution
            var binom = new Binomial(p, n);
            return (int)Math.Round(binom.InverseCDF(rng.NextDouble()));
        }

    }
}
