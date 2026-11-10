using System;
using System.Collections.Generic;
using System.Linq;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Multivariate Student's t-distribution.
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
    /// The multivariate Student's t-distribution is a generalization of the univariate Student's t-distribution
    /// to multiple dimensions. It arises naturally as the marginal distribution of a multivariate normal
    /// distribution with an unknown variance scaled by an inverse chi-squared random variable. It is
    /// parameterized by a degrees of freedom ν, a location vector μ, and a positive-definite scale matrix Σ.
    /// The scale matrix Σ is not the covariance matrix; the covariance is ν/(ν−2)·Σ for ν > 2.
    /// </para>
    /// <para>
    /// The distribution has heavier tails than the multivariate normal, controlled by the degrees of freedom ν.
    /// As ν → ∞, the distribution converges to the multivariate normal. For small ν, the tails are substantially
    /// heavier, making it useful for robust statistical inference and confidence interval construction
    /// where normal approximations underestimate tail probabilities.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <list type="bullet">
    /// <item><description>
    /// Kotz, S. and Nadarajah, S. (2004). "Multivariate t Distributions and Their Applications."
    /// Cambridge University Press.
    /// </description></item>
    /// <item><description>
    /// Wikipedia contributors, "Multivariate t-distribution." Wikipedia, The Free Encyclopedia.
    /// Available at: <see href="https://en.wikipedia.org/wiki/Multivariate_t-distribution"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public class MultivariateStudentT : MultivariateDistribution
    {

        /// <summary>
        /// Private parameterless constructor for use by Clone().
        /// </summary>
        private MultivariateStudentT() { }

        /// <summary>
        /// Constructs a standard multivariate Student's t-distribution with zero location vector,
        /// identity scale matrix, and the specified degrees of freedom.
        /// </summary>
        /// <param name="dimension">The number of dimensions in the distribution.</param>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than zero.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="degreesOfFreedom"/> is not positive,
        /// or <paramref name="dimension"/> is not positive.</exception>
        public MultivariateStudentT(int dimension, double degreesOfFreedom)
        {
            var location = new double[dimension];
            SetParameters(degreesOfFreedom, location, Matrix.Identity(dimension).ToArray());
        }

        /// <summary>
        /// Constructs a multivariate Student's t-distribution with the specified location vector,
        /// identity scale matrix, and degrees of freedom.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than zero.</param>
        /// <param name="location">The location vector μ (mu) for the distribution.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="degreesOfFreedom"/> is not positive,
        /// or <paramref name="location"/> is null.</exception>
        public MultivariateStudentT(double degreesOfFreedom, double[] location)
        {
            SetParameters(degreesOfFreedom, location, Matrix.Identity(location.Length).ToArray());
        }

        /// <summary>
        /// Constructs a multivariate Student's t-distribution with the specified degrees of freedom,
        /// location vector, and scale matrix.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than zero.</param>
        /// <param name="location">The location vector μ (mu) for the distribution.</param>
        /// <param name="scaleMatrix">The scale matrix Σ (sigma) for the distribution. Must be positive definite.
        /// Note: this is the scale matrix, not the covariance matrix. The covariance is ν/(ν−2)·Σ for ν > 2.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when parameters are invalid.</exception>
        public MultivariateStudentT(double degreesOfFreedom, double[] location, double[,] scaleMatrix)
        {
            SetParameters(degreesOfFreedom, location, scaleMatrix);
        }

        private bool _parametersValid = true;
        private int _dimension;
        private double _degreesOfFreedom;
        private double[] _location = null!;
        private Matrix _scaleMatrix = null!;

        private CholeskyDecomposition _cholesky = null!;
        private double _lnconstant;
        private double[]? _variance;
        private double[]? _standardDeviation;
        private Random _MVNUNI = new MersenneTwister(MultivariateNormal.DefaultMVNUNISeed);
        private int _maxEvaluations = 25000;
        private double _absoluteError = 1E-5;
        private double _relativeError = 0d;
        private double[]? _correl;
        private MultivariateNormal? _integrator;

        /// <summary>
        /// The maximum number of integrand evaluations of the lattice rule used by <see cref="CDF(double[])"/>
        /// above two dimensions. Default = 25,000.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the value is not positive.</exception>
        public int MaxEvaluations
        {
            get { return _maxEvaluations; }
            set
            {
                if (value <= 0) throw new ArgumentOutOfRangeException(nameof(MaxEvaluations), "The maximum number of evaluations must be positive.");
                _maxEvaluations = value;
            }
        }

        /// <summary>
        /// The absolute error tolerance of the lattice rule used by <see cref="CDF(double[])"/> above two
        /// dimensions. Default = 1E-5.
        /// </summary>
        /// <remarks>
        /// The rule stops when its 99% error estimate is at most the larger of this value and
        /// <see cref="RelativeError"/> times the probability, or when <see cref="MaxEvaluations"/> is reached,
        /// whichever comes first. With the default budget of 25,000 evaluations the error in three to five
        /// dimensions is typically of order 1E-5; a tolerance of 1E-7 typically needs about 10^6 evaluations.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the value is negative or NaN.</exception>
        public double AbsoluteError
        {
            get { return _absoluteError; }
            set
            {
                if (!(value >= 0d)) throw new ArgumentOutOfRangeException(nameof(AbsoluteError), "The absolute error must be non-negative.");
                _absoluteError = value;
            }
        }

        /// <summary>
        /// The relative error tolerance of the lattice rule used by <see cref="CDF(double[])"/> above two
        /// dimensions. Default = 0, so only <see cref="AbsoluteError"/> applies.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the value is negative or NaN.</exception>
        public double RelativeError
        {
            get { return _relativeError; }
            set
            {
                if (!(value >= 0d)) throw new ArgumentOutOfRangeException(nameof(RelativeError), "The relative error must be non-negative.");
                _relativeError = value;
            }
        }

        /// <summary>
        /// The uniform(0,1) random number generator used by the lattice rule of <see cref="CDF(double[])"/> for
        /// dimensions greater than two.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Above two dimensions <see cref="CDF(double[])"/> is the Genz-Bretz randomized Korobov lattice rule,
        /// which draws its lattice shifts from this generator. The CDF therefore carries a small stochastic
        /// error above two dimensions and only reproduces when the generator is seeded. At one and two
        /// dimensions the CDF is deterministic and never touches this generator.
        /// </para>
        /// <para>
        /// The default is a Mersenne twister seeded with <see cref="MultivariateNormal.DefaultMVNUNISeed"/>,
        /// which makes the default result reproducible, but every instance left at the default replays the
        /// identical lattice shifts, so the quadrature errors of separate instances are correlated rather
        /// than independent and do not average out when many evaluations are aggregated. Assign a seeded
        /// generator to tie results to a caller's own seed and to decorrelate the error across instances;
        /// this mirrors <see cref="MultivariateNormal.MVNUNI"/>.
        /// </para>
        /// <para>
        /// The lattice rule advances the generator, so successive CDF evaluations on the same instance consume
        /// successive shifts and are not bit-identical to one another above two dimensions; two freshly
        /// constructed instances with the same parameters and the same seed are.
        /// </para>
        /// <para>
        /// <b>Not thread-safe.</b> <see cref="MersenneTwister"/> has no internal synchronization, so above
        /// two dimensions a single instance must not have <see cref="CDF(double[])"/> called concurrently
        /// from several threads. Give each thread its own instance, or assign each thread's instance a
        /// generator of its own (<c>mvt.MVNUNI = new MersenneTwister(seedForThisThread)</c>).
        /// <see cref="Clone"/> is not a remedy: the property is typed <see cref="Random"/>, which exposes
        /// no general deep copy, so a clone shares one generator with the original and races exactly as
        /// the original would.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when the assigned generator is null.</exception>
        public Random MVNUNI
        {
            get { return _MVNUNI; }
            set { _MVNUNI = value ?? throw new ArgumentNullException(nameof(MVNUNI)); }
        }

        /// <summary>
        /// Gets the number of variables for the distribution.
        /// </summary>
        public override int Dimension
        {
            get { return _dimension; }
        }

        /// <summary>
        /// Returns the multivariate distribution type.
        /// </summary>
        public override MultivariateDistributionType Type
        {
            get { return MultivariateDistributionType.MultivariateStudentT; }
        }

        /// <summary>
        /// Returns the display name of the distribution type as a string.
        /// </summary>
        public override string DisplayName
        {
            get { return "Multivariate Student's t"; }
        }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public override string ShortDisplayName
        {
            get { return "Multi-T"; }
        }

        /// <summary>
        /// Determines whether the parameters are valid or not.
        /// </summary>
        public override bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <summary>
        /// Gets the degrees of freedom ν (nu) of the distribution.
        /// </summary>
        public double DegreesOfFreedom
        {
            get { return _degreesOfFreedom; }
        }

        /// <summary>
        /// Gets the location vector μ (mu) of the distribution.
        /// </summary>
        public double[] Location
        {
            get { return _location; }
        }

        /// <summary>
        /// Gets the mean vector of the distribution. Equal to the location vector μ when ν > 1.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when ν ≤ 1, as the mean is undefined.</exception>
        public double[] Mean
        {
            get
            {
                if (_degreesOfFreedom <= 1.0)
                    throw new InvalidOperationException("The mean is undefined for degrees of freedom ν ≤ 1.");
                return _location;
            }
        }

        /// <summary>
        /// Gets the median vector of the distribution. Equal to the location vector μ.
        /// </summary>
        public double[] Median
        {
            get { return _location; }
        }

        /// <summary>
        /// Gets the mode vector of the distribution. Equal to the location vector μ.
        /// </summary>
        public double[] Mode
        {
            get { return _location; }
        }

        /// <summary>
        /// Gets the marginal variance vector of the distribution: ν/(ν−2) · diag(Σ).
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when ν ≤ 2, as the variance is undefined.</exception>
        public double[] Variance
        {
            get
            {
                if (_degreesOfFreedom <= 2.0)
                    throw new InvalidOperationException("The variance is undefined for degrees of freedom ν ≤ 2.");
                if (_variance == null)
                {
                    double scale = _degreesOfFreedom / (_degreesOfFreedom - 2.0);
                    var diag = _scaleMatrix.Diagonal();
                    _variance = new double[Dimension];
                    for (int i = 0; i < Dimension; i++)
                        _variance[i] = scale * diag[i];
                }
                return _variance;
            }
        }

        /// <summary>
        /// Gets the marginal standard deviation vector of the distribution.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when ν ≤ 2, as the variance is undefined.</exception>
        public double[] StandardDeviation
        {
            get
            {
                if (_standardDeviation == null)
                {
                    _standardDeviation = new double[Dimension];
                    var v = Variance;
                    for (int i = 0; i < Dimension; i++)
                        _standardDeviation[i] = Math.Sqrt(v[i]);
                }
                return _standardDeviation;
            }
        }

        /// <summary>
        /// Gets the scale matrix Σ (sigma) of the distribution.
        /// </summary>
        /// <remarks>
        /// The scale matrix is not the covariance matrix. The covariance matrix is ν/(ν−2)·Σ for ν > 2.
        /// </remarks>
        public double[,] ScaleMatrix
        {
            get { return _scaleMatrix.ToArray(); }
        }

        /// <summary>
        /// Gets the covariance matrix of the distribution: ν/(ν−2) · Σ.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when ν ≤ 2, as the covariance is undefined.</exception>
        public double[,] Covariance
        {
            get
            {
                if (_degreesOfFreedom <= 2.0)
                    throw new InvalidOperationException("The covariance is undefined for degrees of freedom ν ≤ 2.");
                double scale = _degreesOfFreedom / (_degreesOfFreedom - 2.0);
                var cov = _scaleMatrix.Clone();
                for (int i = 0; i < Dimension; i++)
                    for (int j = 0; j < Dimension; j++)
                        cov[i, j] *= scale;
                return cov.ToArray();
            }
        }

        /// <summary>
        /// Determines if the scale matrix is positive definite.
        /// </summary>
        public bool IsPositiveDefinite => _cholesky.IsPositiveDefinite;

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu). Must be greater than zero.</param>
        /// <param name="location">The location vector μ (mu) for the distribution.</param>
        /// <param name="scaleMatrix">The scale matrix Σ (sigma) for the distribution. Must be positive definite.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when parameters are invalid.</exception>
        public void SetParameters(double degreesOfFreedom, double[] location, double[,] scaleMatrix)
        {
            // Validate parameters. The validation factors its own copy of the scale matrix, and that
            // factor is reused rather than recomputed from identical values.
            ValidateParameters(degreesOfFreedom, location, scaleMatrix, true, out var cholesky);

            _degreesOfFreedom = degreesOfFreedom;
            _dimension = location.Length;
            _location = location;
            _scaleMatrix = new Matrix(scaleMatrix);
            _cholesky = cholesky!;

            // Precompute the log of the normalizing constant for the PDF:
            // ln C = LogGamma((ν+p)/2) - LogGamma(ν/2) - (p/2)*ln(νπ) - (1/2)*ln|Σ|
            double lndet = _cholesky.LogDeterminant();
            _lnconstant = Gamma.LogGamma((_degreesOfFreedom + _dimension) / 2.0)
                        - Gamma.LogGamma(_degreesOfFreedom / 2.0)
                        - (_dimension / 2.0) * Math.Log(_degreesOfFreedom * Math.PI)
                        - 0.5 * lndet;

            // Reset cached derived properties
            _variance = null;
            _standardDeviation = null;
            _correl = null;
            _integrator = null;
        }

        /// <summary>
        /// Validate the distribution parameters.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu).</param>
        /// <param name="location">The location vector μ (mu).</param>
        /// <param name="scaleMatrix">The scale matrix Σ (sigma).</param>
        /// <param name="throwException">Determines whether to throw an exception or return it.</param>
        /// <returns>
        /// An <see cref="ArgumentOutOfRangeException"/> if validation fails; otherwise null.
        /// </returns>
        public ArgumentOutOfRangeException? ValidateParameters(double degreesOfFreedom, double[] location, double[,] scaleMatrix, bool throwException)
        {
            return ValidateParameters(degreesOfFreedom, location, scaleMatrix, throwException, out _);
        }

        /// <summary>
        /// Validate the distribution parameters, handing back the factorization of the scale matrix built
        /// along the way so that <see cref="SetParameters"/> does not have to recompute it.
        /// </summary>
        /// <param name="degreesOfFreedom">The degrees of freedom ν (nu).</param>
        /// <param name="location">The location vector μ (mu).</param>
        /// <param name="scaleMatrix">The scale matrix Σ (sigma).</param>
        /// <param name="throwException">Determines whether to throw an exception or return it.</param>
        /// <param name="cholesky">On return, the Cholesky decomposition of the scale matrix when the
        /// parameters are valid; null otherwise.</param>
        /// <returns>
        /// An <see cref="ArgumentOutOfRangeException"/> if validation fails; otherwise null.
        /// </returns>
        private ArgumentOutOfRangeException? ValidateParameters(double degreesOfFreedom, double[] location, double[,] scaleMatrix, bool throwException, out CholeskyDecomposition? cholesky)
        {
            cholesky = null;
            if (double.IsNaN(degreesOfFreedom) || double.IsInfinity(degreesOfFreedom) || degreesOfFreedom <= 0)
            {
                var ex = new ArgumentOutOfRangeException(nameof(degreesOfFreedom), "Degrees of freedom must be greater than zero.");
                if (throwException) throw ex; else return ex;
            }
            if (location == null)
            {
                var ex = new ArgumentOutOfRangeException(nameof(location), "Location vector must not be null.");
                if (throwException) throw ex; else return ex;
            }
            if (scaleMatrix == null)
            {
                var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale matrix must not be null.");
                if (throwException) throw ex; else return ex;
            }
            for (int i = 0; i < location.Length; i++)
            {
                if (double.IsNaN(location[i]) || double.IsInfinity(location[i]))
                {
                    var ex = new ArgumentOutOfRangeException(nameof(location), "Location values must be finite.");
                    if (throwException) throw ex; else return ex;
                }
            }
            var m = new Matrix(scaleMatrix);
            if (!m.IsSquare)
            {
                var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale matrix must be square.");
                if (throwException) throw ex; else return ex;
            }
            if (m.NumberOfRows != location.Length)
            {
                var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Location vector length must match scale matrix dimension.");
                if (throwException) throw ex; else return ex;
            }
            for (int i = 0; i < m.NumberOfRows; i++)
            {
                for (int j = 0; j < m.NumberOfColumns; j++)
                {
                    if (double.IsNaN(m[i, j]) || double.IsInfinity(m[i, j]))
                    {
                        var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale-matrix values must be finite.");
                        if (throwException) throw ex; else return ex;
                    }
                }
            }
            try
            {
                var chol = new CholeskyDecomposition(m);
                if (!chol.IsPositiveDefinite)
                {
                    var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale matrix is not positive-definite.");
                    if (throwException) throw ex; else return ex;
                }
                cholesky = chol;
            }
            catch (ArgumentOutOfRangeException) { throw; }
            catch (Exception)
            {
                var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale matrix is not positive-definite.");
                if (throwException) throw ex; else return ex;
            }
            return null;
        }

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        /// <returns>The probability density at point x.</returns>
        /// <remarks>
        /// <para>
        /// The PDF of the multivariate Student's t-distribution is:
        /// </para>
        /// <para>
        /// f(x) = C · [1 + (x−μ)'Σ⁻¹(x−μ)/ν]^(−(ν+p)/2)
        /// </para>
        /// <para>
        /// where C = Γ((ν+p)/2) / [Γ(ν/2) · (νπ)^(p/2) · |Σ|^(1/2)].
        /// </para>
        /// </remarks>
        public override double PDF(double[] x)
        {
            return Math.Exp(LogPDF(x));
        }

        /// <summary>
        /// Returns the natural log of the PDF evaluated at a point x.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        /// <returns>The natural logarithm of the probability density at point x.</returns>
        /// <remarks>
        /// <para>
        /// Uses the numerically stable formulation:
        /// </para>
        /// <para>
        /// log f(x) = ln C − ((ν+p)/2) · ln(1 + Q/ν)
        /// </para>
        /// <para>
        /// where Q = (x−μ)'Σ⁻¹(x−μ) is the squared Mahalanobis distance, and ln C is the precomputed
        /// log normalizing constant.
        /// </para>
        /// </remarks>
        public override double LogPDF(double[] x)
        {
            double Q = Mahalanobis(x);
            double f = _lnconstant - ((_degreesOfFreedom + _dimension) / 2.0) * Math.Log(1.0 + Q / _degreesOfFreedom);
            if (double.IsNaN(f) || double.IsInfinity(f)) return double.NegativeInfinity;
            return f;
        }

        /// <summary>
        /// Gets the squared Mahalanobis distance between a sample point and this distribution:
        /// (x−μ)'Σ⁻¹(x−μ).
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        /// <returns>The squared Mahalanobis distance.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the vector dimension does not match
        /// the distribution dimension.</exception>
        public double Mahalanobis(double[] x)
        {
            if (x.Length != Dimension)
                throw new ArgumentOutOfRangeException(nameof(x), "The vector must be the same dimension as the distribution.");
            var z = new double[_location.Length];
            for (int i = 0; i < x.Length; i++)
                z[i] = x[i] - _location[i];
            var a = _cholesky.Solve(new Vector(z));
            double b = 0.0;
            for (int i = 0; i < z.Length; i++)
                b += a[i] * z[i];
            return b;
        }

        /// <summary>
        /// Returns the radial factor r = sqrt(W / ν) of the multivariate t integrand, where W is the χ²(ν) quantile of u.
        /// </summary>
        /// <param name="halfNu">Half the degrees of freedom, a = ν / 2 &gt; 0.</param>
        /// <param name="logGammaHalfNu">ln Γ(a), computed once by the caller.</param>
        /// <param name="u">The uniform coordinate in [0, 1).</param>
        /// <returns>sqrt(x / a), where x solves P(a, x) = u; zero for u = 0.</returns>
        /// <remarks>
        /// W / 2 has the Gamma(a, 1) distribution, so r = sqrt(x / a) with x the inverse of the regularized lower
        /// incomplete gamma function, found by <see cref="Gamma.InverseIncompleteRatio(double, double, bool, double)"/>.
        /// </remarks>
        internal static double RadialFactor(double halfNu, double logGammaHalfNu, double u)
        {
            return Math.Sqrt(Gamma.InverseIncompleteRatio(halfNu, u, false, logGammaHalfNu) / halfNu);
        }

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        /// <returns>The cumulative probability P(X ≤ x).</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the vector dimension does not match the
        /// distribution dimension.</exception>
        /// <remarks>
        /// <para>
        /// The CDF is computed with the Genz-Bretz method, a translation of Alan Genz's Fortran <c>MVTDST</c>.
        /// A coordinate at +∞ is unconstrained and drops out, a coordinate at −∞ gives zero, and NaN gives NaN.
        /// The remaining limits are standardized by the scale-matrix diagonal and sorted, and the scale
        /// correlation is factored, by the same routine as <see cref="MultivariateNormal.CDF(double[])"/>.
        /// </para>
        /// <para>
        /// One remaining variable uses the univariate Student's t CDF. Two remaining variables are
        /// deterministic: integer ν up to 10,000 uses the Dunnett and Sobel (1954) closed form, which is exact
        /// to rounding, and any other ν integrates the exact bivariate normal probability over the χ(ν)/√ν
        /// mixing variable by adaptive Gauss-Kronrod quadrature, which agrees with exact references to about 1E-14.
        /// </para>
        /// <para>
        /// Three or more remaining variables use the randomized Korobov lattice rule of MVNDST with one extra
        /// coordinate for the mixing variable. The rule refines until its 99% error estimate meets
        /// <see cref="AbsoluteError"/> or <see cref="RelativeError"/> times the value, or until
        /// <see cref="MaxEvaluations"/> integrand evaluations are spent.
        /// </para>
        /// <para>
        /// <b>Above two dimensions this method is stochastic, stateful, and not thread-safe.</b> The lattice
        /// rule draws its shifts from <see cref="MVNUNI"/>, which is instance state; up to two dimensions the
        /// CDF is deterministic, touches no generator or integration workspace, and may be called concurrently
        /// on one instance. Each call advances the generator, so repeated calls
        /// at the same point on the same instance return slightly different values; it is two freshly
        /// constructed instances with the same parameters and seed that agree bit for bit. Because
        /// <see cref="MersenneTwister"/> has no internal synchronization, calling this method concurrently
        /// on a single shared instance is a data race. Use one instance per thread, or a distinct generator
        /// per thread assigned through <see cref="MVNUNI"/>; <see cref="Clone"/> shares the generator by
        /// reference and does not help.
        /// </para>
        /// <para>
        /// References: Genz, A. and Bretz, F. (2009). "Computation of Multivariate Normal and t Probabilities."
        /// Lecture Notes in Statistics, Vol. 195. Springer. Dunnett, C. W. and Sobel, M. (1954). "A bivariate
        /// generalization of Student's t-distribution with tables for certain special cases." Biometrika 41,
        /// 153-169.
        /// </para>
        /// </remarks>
        public override double CDF(double[] x)
        {
            if (x.Length != Dimension)
                throw new ArgumentOutOfRangeException(nameof(x), "The vector must be the same dimension as the distribution.");

            if (Dimension == 1)
            {
                // Delegate to univariate Student's t CDF
                double sigma = Math.Sqrt(_scaleMatrix[0, 0]);
                var univT = new StudentT(_location[0], sigma, _degreesOfFreedom);
                return univT.CDF(x[0]);
            }

            var lower = new double[Dimension];
            var upper = new double[Dimension];
            var infin = new int[Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                if (double.IsNaN(x[i])) return double.NaN;
                if (double.IsNegativeInfinity(x[i])) return 0d;
                if (double.IsPositiveInfinity(x[i]))
                {
                    infin[i] = -1;
                }
                else
                {
                    upper[i] = (x[i] - _location[i]) / Math.Sqrt(_scaleMatrix[i, i]);
                }
            }

            // The packed correlations are filled before they are published, so a concurrent caller never reads a
            // partly filled array.
            var correl = _correl;
            if (correl == null)
            {
                correl = new double[Dimension * (Dimension - 1) / 2];
                int k = 0;
                for (int i = 1; i < Dimension; i++)
                    for (int j = 0; j < i; j++)
                        correl[k++] = _scaleMatrix[i, j] / Math.Sqrt(_scaleMatrix[i, i] * _scaleMatrix[j, j]);
                _correl = correl;
            }
            var integrator = _integrator ??= new MultivariateNormal(new double[Dimension], _scaleMatrix.ToArray());
            integrator.MVNUNI = _MVNUNI;

            double error = 0d, value = 0d;
            int inform = 0;
            integrator.MVTDST(Dimension, _degreesOfFreedom, lower, upper, infin, correl, _maxEvaluations, _absoluteError, _relativeError, ref error, ref value, ref inform);
            return value;
        }

        /// <summary>
        /// Generate a matrix of random values from the multivariate Student's t-distribution.
        /// </summary>
        /// <param name="sampleSize">Size of random sample to generate.</param>
        /// <param name="seed">Optional seed for the random number generator. Use -1 for a random seed.</param>
        /// <returns>
        /// A 2D array of random values. The number of rows equals the sample size.
        /// The number of columns equals the dimension of the distribution.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Sampling uses the representation: X = μ + L·z / √(W/ν),
        /// where L is the Cholesky factor of Σ, z ~ N(0, I_p), and W ~ χ²(ν).
        /// The χ² variate is generated via the equivalent Gamma(ν/2, 2) distribution
        /// to support non-integer degrees of freedom.
        /// </para>
        /// </remarks>
        public double[,] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            return GenerateRandomValues(sampleSize, seed > 0 ? new MersenneTwister(seed) : new MersenneTwister());
        }

        /// <summary>
        /// Generate a matrix of random values from the multivariate Student's t-distribution, drawing from a
        /// supplied generator.
        /// </summary>
        /// <param name="sampleSize">Size of random sample to generate.</param>
        /// <param name="random">The generator. Each row draws <see cref="Dimension"/> values for the normal
        /// component and then one for the χ²(ν) mixing variate.</param>
        /// <returns>
        /// A 2D array of random values. The number of rows equals the sample size.
        /// The number of columns equals the dimension of the distribution.
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="random"/> is null.</exception>
        /// <remarks>
        /// The sample is a deterministic function of the generator state, with no clock fallback. Passing
        /// <c>new MersenneTwister(seed)</c> with a positive seed reproduces
        /// <see cref="GenerateRandomValues(int, int)"/> with that seed exactly.
        /// </remarks>
        public double[,] GenerateRandomValues(int sampleSize, Random random)
        {
            if (random == null) throw new ArgumentNullException(nameof(random));
            var rnd = random;
            var sample = new double[sampleSize, Dimension];

            // Use Gamma(ν/2, 2) to generate χ²(ν) variates, supporting non-integer ν
            var gamma = new GammaDistribution(2.0, _degreesOfFreedom / 2.0);

            for (int i = 0; i < sampleSize; i++)
            {
                // z ~ N(0, I_p)
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                    z[j] = Normal.StandardZ(rnd.NextDouble());

                // W ~ χ²(ν) via Gamma(ν/2, 2)
                double w = gamma.InverseCDF(rnd.NextDouble());

                // scale = √(ν / W)
                double scale = Math.Sqrt(_degreesOfFreedom / w);

                // x = μ + L·z · scale
                var Lz = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = _location[j] + Lz[j] * scale;
            }
            return sample;
        }

        /// <summary>
        /// Generate random values using Latin Hypercube Sampling (LHS) for improved space-filling properties.
        /// </summary>
        /// <param name="sampleSize">Size of random sample to generate.</param>
        /// <param name="seed">Seed for the random number generator.</param>
        /// <returns>
        /// A 2D array of random values with improved coverage of the distribution space.
        /// </returns>
        /// <remarks>
        /// <para>
        /// Latin Hypercube Sampling ensures that each marginal dimension is evenly stratified.
        /// The LHS uniforms are used for the normal component z ~ N(0, I_p), while a separate
        /// independent random stream generates the χ²(ν) mixing variates. This preserves the
        /// space-filling properties of LHS in the correlated normal dimensions while maintaining
        /// correct marginal t-distribution behavior through the independent χ² scaling.
        /// </para>
        /// </remarks>
        public double[,] LatinHypercubeRandomValues(int sampleSize, int seed)
        {
            var r = LatinHypercube.Random(sampleSize, Dimension, seed);
            var sample = new double[sampleSize, Dimension];

            // Separate PRNG stream for χ² variates (independent of the LHS grid)
            var rndChi = new MersenneTwister(seed + 1);
            var gamma = new GammaDistribution(2.0, _degreesOfFreedom / 2.0);

            for (int i = 0; i < sampleSize; i++)
            {
                // z ~ N(0, I_p) via LHS uniforms
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                    z[j] = Normal.StandardZ(r[i, j]);

                // W ~ χ²(ν) via independent stream
                double w = gamma.InverseCDF(rndChi.NextDouble());
                double scale = Math.Sqrt(_degreesOfFreedom / w);

                // x = μ + L·z · scale
                var Lz = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = _location[j] + Lz[j] * scale;
            }
            return sample;
        }

        /// <summary>
        /// Returns a 2D array of stratified random variates. The first dimension is stratified,
        /// and the remaining dimensions are sampled randomly.
        /// </summary>
        /// <param name="stratificationBins">A list of stratification bins defining the stratified dimension.</param>
        /// <param name="seed">Seed for the random number generator.</param>
        /// <returns>
        /// A 2D array of stratified random values.
        /// </returns>
        /// <remarks>
        /// <para>
        /// The first marginal dimension uses stratified uniform quantiles (bin midpoints),
        /// while the remaining dimensions use independent random uniforms. The χ²(ν) mixing
        /// variate is drawn from a separate random stream.
        /// </para>
        /// </remarks>
        public double[,] StratifiedRandomValues(List<StratificationBin> stratificationBins, int seed)
        {
            int samplesize = stratificationBins.Count;
            var rnd = new MersenneTwister(seed);
            var rndChi = new MersenneTwister(seed + 1);
            var gamma = new GammaDistribution(2.0, _degreesOfFreedom / 2.0);
            var sample = new double[samplesize, Dimension];

            for (int i = 0; i < samplesize; i++)
            {
                // z ~ N(0, I_p) with first dimension stratified
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                {
                    if (j == 0)
                        z[j] = Normal.StandardZ(stratificationBins[i].Midpoint);
                    else
                        z[j] = Normal.StandardZ(rnd.NextDouble());
                }

                // W ~ χ²(ν) via independent stream
                double w = gamma.InverseCDF(rndChi.NextDouble());
                double scale = Math.Sqrt(_degreesOfFreedom / w);

                // x = μ + L·z · scale
                var Lz = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = _location[j] + Lz[j] * scale;
            }
            return sample;
        }

        /// <summary>
        /// The inverse cumulative distribution function (InverseCDF).
        /// </summary>
        /// <param name="probabilities">
        /// Array of p+1 uniform probabilities in (0, 1). The first p entries correspond to the
        /// correlated normal dimensions, and the last entry is the probability for the χ²(ν)
        /// mixing variable that controls the tail thickness.
        /// </param>
        /// <returns>A p-dimensional sample point from the distribution.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the length of
        /// <paramref name="probabilities"/> is not equal to Dimension + 1.</exception>
        /// <remarks>
        /// <para>
        /// Uses the representation X = μ + L·z · √(ν/W), where z_j = Φ⁻¹(p_j) for j = 1..p
        /// and W = F_{χ²(ν)}⁻¹(p_{p+1}). L is the Cholesky factor of the scale matrix Σ.
        /// </para>
        /// <para>
        /// The extra probability (compared to <see cref="MultivariateNormal.InverseCDF"/>) controls
        /// the χ² mixing variate. Values near 0 produce extreme tail samples (large scaling),
        /// while values near 1 produce samples close to the multivariate normal.
        /// </para>
        /// </remarks>
        public double[] InverseCDF(double[] probabilities)
        {
            if (probabilities.Length != Dimension + 1)
                throw new ArgumentOutOfRangeException(nameof(probabilities),
                    $"The probabilities array must have length {Dimension + 1} (Dimension + 1 for the χ² mixing variable).");

            // Convert first p probabilities to standard normal variates
            var z = new double[Dimension];
            for (int j = 0; j < Dimension; j++)
                z[j] = Normal.StandardZ(probabilities[j]);

            // Convert last probability to χ²(ν) variate via Gamma(ν/2, 2)
            var gamma = new GammaDistribution(2.0, _degreesOfFreedom / 2.0);
            double w = gamma.InverseCDF(probabilities[Dimension]);
            double scale = Math.Sqrt(_degreesOfFreedom / w);

            // x = μ + L·z · scale
            var Lz = _cholesky.L * z;
            var sample = new double[Dimension];
            for (int j = 0; j < Dimension; j++)
                sample[j] = _location[j] + Lz[j] * scale;

            return sample;
        }

        /// <summary>
        /// Creates a deep copy of this distribution.
        /// </summary>
        /// <returns>A new <see cref="MultivariateStudentT"/> instance with identical parameters.</returns>
        /// <remarks>
        /// The parameters, and the factorization built from them, are copied deeply. <see cref="MVNUNI"/>
        /// is the one exception: the property is typed <see cref="Random"/>, which exposes no general deep
        /// copy, so it is copied by reference and the clone shares a single generator with the original.
        /// Cloning is therefore not a way to make concurrent <see cref="CDF(double[])"/> calls safe above
        /// two dimensions, and the clone's CDF results are not reproducible independently of the
        /// original's; assign the clone its own generator when either matters.
        /// </remarks>
        public override MultivariateDistribution Clone()
        {
            var clone = new MultivariateStudentT()
            {
                _parametersValid = this._parametersValid,
                _dimension = this._dimension,
                _degreesOfFreedom = this._degreesOfFreedom,
                _location = this._location.ToArray(),
                _scaleMatrix = this._scaleMatrix.Clone(),
                _cholesky = new CholeskyDecomposition(this._scaleMatrix.Clone()),
                _lnconstant = this._lnconstant,
                _MVNUNI = this._MVNUNI,
                _maxEvaluations = this._maxEvaluations,
                _absoluteError = this._absoluteError,
                _relativeError = this._relativeError,
            };
            return clone;
        }

    }
}
