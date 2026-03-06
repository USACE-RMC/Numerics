/*
* NOTICE:
* The U.S. Army Corps of Engineers, Risk Management Center (USACE-RMC) makes no guarantees about
* the results, or appropriateness of outputs, obtained from Numerics.
*
* LIST OF CONDITIONS:
* Redistribution and use in source and binary forms, with or without modification, are permitted
* provided that the following conditions are met:
* ● Redistributions of source code must retain the above notice, this list of conditions, and the
* following disclaimer.
* ● Redistributions in binary form must reproduce the above notice, this list of conditions, and
* the following disclaimer in the documentation and/or other materials provided with the distribution.
* ● The names of the U.S. Government, the U.S. Army Corps of Engineers, the Institute for Water
* Resources, or the Risk Management Center may not be used to endorse or promote products derived
* from this software without specific prior written permission. Nor may the names of its contributors
* be used to endorse or promote products derived from this software without specific prior
* written permission.
*
* DISCLAIMER:
* THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY CORPS OF ENGINEERS RISK MANAGEMENT CENTER
* (USACE-RMC) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
* THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL USACE-RMC BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
* THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
            // Validate parameters
            ValidateParameters(degreesOfFreedom, location, scaleMatrix, true);

            _degreesOfFreedom = degreesOfFreedom;
            _dimension = location.Length;
            _location = location;
            _scaleMatrix = new Matrix(scaleMatrix);
            _cholesky = new CholeskyDecomposition(_scaleMatrix);

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
            if (degreesOfFreedom <= 0)
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
            try
            {
                var chol = new CholeskyDecomposition(m);
                if (!chol.IsPositiveDefinite)
                {
                    var ex = new ArgumentOutOfRangeException(nameof(scaleMatrix), "Scale matrix is not positive-definite.");
                    if (throwException) throw ex; else return ex;
                }
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
            if (double.IsNaN(f) || double.IsInfinity(f)) return double.MinValue;
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
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point x.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        /// <returns>The cumulative probability P(X ≤ x).</returns>
        /// <remarks>
        /// <para>
        /// For dimension 1, this delegates to the univariate Student's t CDF.
        /// </para>
        /// <para>
        /// For dimensions ≥ 2, the CDF is computed via stratified numerical integration over the
        /// mixing representation: P(X ≤ x) = E_W[ Φ_p((x−μ)·√(W/ν); 0, Σ) ] where W ~ χ²(ν),
        /// and Φ_p is the multivariate normal CDF. The expectation is computed by evaluating the
        /// MVN CDF at K=200 stratified quantiles of the χ²(ν) distribution and averaging.
        /// </para>
        /// <para>
        /// Reference: Genz, A. and Bretz, F. (2009). "Computation of Multivariate Normal and t Probabilities."
        /// Lecture Notes in Statistics, Vol. 195. Springer.
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

            // For D >= 2, use stratified quantile integration over the χ²(ν) mixing variable.
            // The multivariate-t CDF can be written as:
            //   P(X ≤ x) = E_W[ Φ_MVN((x−μ)·√(W/ν); 0, Σ) ]   where W ~ χ²(ν)
            //
            // We evaluate this by computing the MVN CDF at K equally-spaced quantiles
            // of χ²(ν) and averaging. This is deterministic and works for any ν.

            const int K = 200;
            var gamma = new GammaDistribution(2.0, _degreesOfFreedom / 2.0);

            // Centered point for MVN CDF evaluation
            var zVec = new double[Dimension];
            for (int i = 0; i < Dimension; i++)
                zVec[i] = x[i] - _location[i];

            // Create MVN with zero mean and the scale matrix Σ for CDF evaluation
            var mvn = new MultivariateNormal(new double[Dimension], _scaleMatrix.ToArray());

            double sum = 0.0;
            for (int k = 0; k < K; k++)
            {
                // Use midpoint of each stratum for stratified integration
                double p = (k + 0.5) / K;
                double w = gamma.InverseCDF(p);
                double scaleFactor = Math.Sqrt(w / _degreesOfFreedom);

                // Scale the centered vector: (x-μ) · √(W/ν)
                var scaledZ = new double[Dimension];
                for (int i = 0; i < Dimension; i++)
                    scaledZ[i] = zVec[i] * scaleFactor;

                sum += mvn.CDF(scaledZ);
            }

            double result = sum / K;

            // Clamp to [0, 1]
            return Math.Max(0.0, Math.Min(1.0, result));
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
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
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
            };
            return clone;
        }

    }
}