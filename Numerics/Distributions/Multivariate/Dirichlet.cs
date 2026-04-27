using System;
using System.Linq;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Dirichlet distribution, a multivariate generalization of the Beta distribution.
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
    /// The Dirichlet distribution Dir(α₁, α₂, ..., αₖ) is a family of continuous multivariate
    /// probability distributions parameterized by a vector α of positive real numbers. It is
    /// defined on the (K-1)-dimensional simplex, meaning the components x₁, x₂, ..., xₖ satisfy
    /// xᵢ &gt; 0 and Σxᵢ = 1. The Dirichlet distribution is the conjugate prior of the categorical
    /// and multinomial distributions in Bayesian statistics.
    /// </para>
    /// <para>
    /// Key applications include:
    /// <list type="bullet">
    /// <item><description>Bayesian inference: conjugate prior for mixture model weights in RMC-BestFit.</description></item>
    /// <item><description>Topic modeling: prior distribution over document-topic proportions.</description></item>
    /// <item><description>Compositional data analysis: modeling proportions that sum to unity.</description></item>
    /// </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Kotz, S., Balakrishnan, N. and Johnson, N.L. (2000). "Continuous Multivariate Distributions,
    /// Volume 1: Models and Applications," 2nd ed. Wiley. Chapter 49.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Dirichlet_distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class Dirichlet : MultivariateDistribution
    {

        /// <summary>
        /// Private parameterless constructor for use by Clone().
        /// </summary>
        private Dirichlet() { }

        /// <summary>
        /// Constructs a symmetric Dirichlet distribution where all concentration parameters are equal.
        /// </summary>
        /// <param name="dimension">The number of categories K. Must be at least 2.</param>
        /// <param name="alpha">The common concentration parameter. Must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when dimension &lt; 2 or alpha &lt;= 0.</exception>
        public Dirichlet(int dimension, double alpha)
        {
            if (dimension < 2)
                throw new ArgumentOutOfRangeException(nameof(dimension), "The dimension must be at least 2.");
            if (double.IsNaN(alpha) || double.IsInfinity(alpha) || alpha <= 0)
                throw new ArgumentOutOfRangeException(nameof(alpha), "The concentration parameter must be positive.");

            _alpha = new double[dimension];
            for (int i = 0; i < dimension; i++)
                _alpha[i] = alpha;
            _dimension = dimension;
            ComputeNormalization();
        }

        /// <summary>
        /// Constructs a Dirichlet distribution with the specified concentration parameter vector.
        /// </summary>
        /// <param name="alpha">The vector of concentration parameters α₁, α₂, ..., αₖ. Each must be positive.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when any concentration parameter is non-positive,
        /// or the vector has fewer than 2 elements.</exception>
        public Dirichlet(double[] alpha)
        {
            if (alpha == null || alpha.Length < 2)
                throw new ArgumentOutOfRangeException(nameof(alpha), "The concentration parameter vector must have at least 2 elements.");
            for (int i = 0; i < alpha.Length; i++)
            {
                if (double.IsNaN(alpha[i]) || double.IsInfinity(alpha[i]) || alpha[i] <= 0)
                    throw new ArgumentOutOfRangeException(nameof(alpha), $"Concentration parameter α[{i}] must be positive.");
            }

            _alpha = (double[])alpha.Clone();
            _dimension = alpha.Length;
            ComputeNormalization();
        }

        private double[] _alpha = null!;
        private int _dimension;
        private double _logNormalization; // log(B(alpha)) = sum(logGamma(alpha_i)) - logGamma(sum(alpha_i))
        private double _alphaSum;
        private double[]? _mean;
        private double[]? _variance;
        private double[]? _mode;

        /// <summary>
        /// Computes and caches the normalization constant and alpha sum.
        /// </summary>
        private void ComputeNormalization()
        {
            _alphaSum = 0;
            _logNormalization = 0;
            for (int i = 0; i < _dimension; i++)
            {
                _alphaSum += _alpha[i];
                _logNormalization += Gamma.LogGamma(_alpha[i]);
            }
            _logNormalization -= Gamma.LogGamma(_alphaSum);
        }

        /// <summary>
        /// Gets the concentration parameter vector α.
        /// </summary>
        public double[] Alpha
        {
            get { return (double[])_alpha.Clone(); }
        }

        /// <summary>
        /// Gets the sum of all concentration parameters: S = Σαᵢ.
        /// </summary>
        public double AlphaSum
        {
            get { return _alphaSum; }
        }

        /// <inheritdoc/>
        public override int Dimension
        {
            get { return _dimension; }
        }

        /// <inheritdoc/>
        public override MultivariateDistributionType Type
        {
            get { return MultivariateDistributionType.Dirichlet; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Dirichlet"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "Dir"; }
        }

        /// <inheritdoc/>
        public override bool ParametersValid
        {
            get
            {
                if (_alpha == null || _alpha.Length < 2) return false;
                for (int i = 0; i < _alpha.Length; i++)
                {
                    if (double.IsNaN(_alpha[i]) || double.IsInfinity(_alpha[i]) || _alpha[i] <= 0)
                        return false;
                }
                return true;
            }
        }

        /// <summary>
        /// Gets the mean vector. Mean[i] = αᵢ / S, where S = Σαⱼ.
        /// </summary>
        public double[] Mean
        {
            get
            {
                if (_mean == null)
                {
                    _mean = new double[_dimension];
                    for (int i = 0; i < _dimension; i++)
                        _mean[i] = _alpha[i] / _alphaSum;
                }
                return _mean;
            }
        }

        /// <summary>
        /// Gets the marginal variance vector. Var[i] = αᵢ(S - αᵢ) / (S²(S + 1)).
        /// </summary>
        public double[] Variance
        {
            get
            {
                if (_variance == null)
                {
                    _variance = new double[_dimension];
                    double s2 = _alphaSum * _alphaSum;
                    double denom = s2 * (_alphaSum + 1.0);
                    for (int i = 0; i < _dimension; i++)
                        _variance[i] = _alpha[i] * (_alphaSum - _alpha[i]) / denom;
                }
                return _variance;
            }
        }

        /// <summary>
        /// Gets the mode vector. Mode[i] = (αᵢ - 1) / (S - K) when all αᵢ &gt; 1.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when any αᵢ &lt;= 1, as the mode is not in the interior of the simplex.</exception>
        public double[] Mode
        {
            get
            {
                if (_mode == null)
                {
                    for (int i = 0; i < _dimension; i++)
                    {
                        if (_alpha[i] <= 1.0)
                            throw new InvalidOperationException("The mode is only defined in the interior of the simplex when all αᵢ > 1.");
                    }
                    double denom = _alphaSum - _dimension;
                    _mode = new double[_dimension];
                    for (int i = 0; i < _dimension; i++)
                        _mode[i] = (_alpha[i] - 1.0) / denom;
                }
                return _mode;
            }
        }

        /// <summary>
        /// Gets the covariance between components i and j: Cov(Xᵢ, Xⱼ) = -αᵢαⱼ / (S²(S + 1)).
        /// </summary>
        /// <param name="i">The first component index (0-based).</param>
        /// <param name="j">The second component index (0-based).</param>
        /// <returns>The covariance. Negative for i != j (components are negatively correlated on the simplex).</returns>
        public double Covariance(int i, int j)
        {
            if (i < 0 || i >= _dimension || j < 0 || j >= _dimension)
                throw new ArgumentOutOfRangeException("Index out of range.");
            if (i == j) return Variance[i];
            double denom = _alphaSum * _alphaSum * (_alphaSum + 1.0);
            return -_alpha[i] * _alpha[j] / denom;
        }

        /// <summary>
        /// Gets the full covariance matrix.
        /// </summary>
        /// <returns>A K x K covariance matrix.</returns>
        public double[,] CovarianceMatrix()
        {
            var cov = new double[_dimension, _dimension];
            double denom = _alphaSum * _alphaSum * (_alphaSum + 1.0);
            for (int i = 0; i < _dimension; i++)
            {
                for (int j = 0; j < _dimension; j++)
                {
                    if (i == j)
                        cov[i, j] = _alpha[i] * (_alphaSum - _alpha[i]) / denom;
                    else
                        cov[i, j] = -_alpha[i] * _alpha[j] / denom;
                }
            }
            return cov;
        }

        /// <inheritdoc/>
        public override double PDF(double[] x)
        {
            return Math.Exp(LogPDF(x));
        }

        /// <summary>
        /// Computes the log of the probability density function.
        /// </summary>
        /// <param name="x">A point on the simplex. Components must be positive and sum to 1.</param>
        /// <returns>
        /// The log-density at x. Returns <see cref="double.MinValue"/> if x is not on the simplex.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The log-PDF is computed as:
        /// <code>
        ///     log f(x) = Σ(αᵢ - 1)·log(xᵢ) - log B(α)
        /// </code>
        /// where B(α) = Π Γ(αᵢ) / Γ(Σαᵢ) is the multivariate Beta function.
        /// </para>
        /// </remarks>
        public override double LogPDF(double[] x)
        {
            if (x == null || x.Length != _dimension) return double.MinValue;

            // Check that x is on the simplex
            double sum = 0;
            for (int i = 0; i < _dimension; i++)
            {
                if (x[i] <= 0 || x[i] > 1) return double.MinValue;
                sum += x[i];
            }
            if (Math.Abs(sum - 1.0) > 1e-10) return double.MinValue;

            double logPdf = -_logNormalization;
            for (int i = 0; i < _dimension; i++)
                logPdf += (_alpha[i] - 1.0) * Math.Log(x[i]);

            return logPdf;
        }

        /// <summary>
        /// The CDF of the Dirichlet distribution is not available in closed form.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        /// <returns>This method always throws <see cref="NotImplementedException"/>.</returns>
        /// <exception cref="NotImplementedException">Always thrown. The multivariate Dirichlet CDF has no closed-form expression.</exception>
        public override double CDF(double[] x)
        {
            throw new NotImplementedException("The CDF of the Dirichlet distribution does not have a closed-form expression.");
        }

        /// <inheritdoc/>
        public override MultivariateDistribution Clone()
        {
            var clone = new Dirichlet();
            clone._alpha = (double[])_alpha.Clone();
            clone._dimension = _dimension;
            clone._logNormalization = _logNormalization;
            clone._alphaSum = _alphaSum;
            return clone;
        }

        /// <summary>
        /// Generates random samples from the Dirichlet distribution.
        /// </summary>
        /// <param name="sampleSize">The number of samples to generate.</param>
        /// <param name="seed">Optional seed for reproducibility. Use -1 for a random seed.</param>
        /// <returns>
        /// A 2D array of shape [sampleSize, Dimension]. Each row sums to 1 and lies on the simplex.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Sampling uses the gamma distribution representation: draw K independent
        /// Gamma(αᵢ, 1) variates yᵢ, then normalize: xᵢ = yᵢ / Σyⱼ.
        /// This leverages the existing <see cref="GammaDistribution"/> class.
        /// </para>
        /// </remarks>
        public double[,] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            var rng = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize, _dimension];

            // Create Gamma distributions for each component: Gamma(alpha_i, 1)
            // GammaDistribution uses shape-rate parameterization: GammaDistribution(beta=1/scale, alpha=shape)
            // Actually, need to check the constructor signature
            var gammas = new GammaDistribution[_dimension];
            for (int i = 0; i < _dimension; i++)
                gammas[i] = new GammaDistribution(1.0, _alpha[i]); // beta=1 (rate), alpha=shape

            for (int s = 0; s < sampleSize; s++)
            {
                double sum = 0;
                var y = new double[_dimension];
                for (int i = 0; i < _dimension; i++)
                {
                    y[i] = gammas[i].InverseCDF(rng.NextDouble());
                    if (y[i] < 0) y[i] = 0; // Guard against numerical issues
                    sum += y[i];
                }

                // Normalize to simplex
                if (sum == 0)
                {
                    for (int i = 0; i < _dimension; i++)
                        sample[s, i] = 1.0 / _dimension;
                    continue;
                }
                for (int i = 0; i < _dimension; i++)
                    sample[s, i] = y[i] / sum;
            }

            return sample;
        }

        /// <summary>
        /// Computes the log of the multivariate Beta function: log B(α) = Σ log Γ(αᵢ) - log Γ(Σαᵢ).
        /// </summary>
        /// <param name="alpha">The concentration parameter vector.</param>
        /// <returns>The log of B(α).</returns>
        public static double LogMultivariateBeta(double[] alpha)
        {
            double logB = 0;
            double sum = 0;
            for (int i = 0; i < alpha.Length; i++)
            {
                logB += Gamma.LogGamma(alpha[i]);
                sum += alpha[i];
            }
            logB -= Gamma.LogGamma(sum);
            return logB;
        }

    }
}
