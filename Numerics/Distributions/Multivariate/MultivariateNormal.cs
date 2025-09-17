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
    /// The Multivariate Normal distribution.
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
    /// <see href = "https://en.wikipedia.org/wiki/Multivariate_normal_distribution" />
    /// </para>
    /// </remarks>
    [Serializable]
    public class MultivariateNormal : MultivariateDistribution
    {
        private MultivariateNormal() { }

        /// <summary>
        /// Constructs a multivariate Gaussian distribution with zero mean vector and identity covariance matrix.
        /// </summary>
        /// <param name="dimension">The number of dimensions in the distribution.</param>
        public MultivariateNormal(int dimension)
        {
            var mean = new double[dimension];
            SetParameters(mean, Matrix.Identity(dimension).ToArray());
        }

        /// <summary>
        /// Constructs a new Multivariate Normal distribution with an identity covariance matrix.
        /// </summary>
        /// <param name="mean">The mean vector μ (mu) for the distribution.</param>
        public MultivariateNormal(double[] mean)
        {
            SetParameters(mean, Matrix.Identity(mean.Length).ToArray());
        }

        /// <summary>
        /// Constructs a new Multivariate Normal distribution.
        /// </summary>
        /// <param name="mean">The mean vector μ (mu) for the distribution.</param>
        /// <param name="covariance">The covariance matrix Σ (sigma) for the distribution.</param>
        public MultivariateNormal(double[] mean, double[,] covariance)
        {
            SetParameters(mean, covariance);
        }

        private bool _parametersValid = true;
        private int _dimension = 0;
        private double[] _mean;
        private Matrix _covariance;
        
        private CholeskyDecomposition _cholesky;
        private double _lnconstant;
        private double[] _variance;
        private double[] _standardDeviation;

        // variables required for the multivariate CDF
        private Matrix _correlation;
        private double[] _correl;
        private Random _MVNUNI = new MersenneTwister();
        private int _maxEvaluations = 100000;
        private double _absoluteError = 1E-4;
        private double _relativeError = 1E-4;
        private double[] _lower;
        private double[] _upper;
        private int[] _infin;
        private bool _correlationMatrixCreated = false;
        private bool _covSRTed = false;

        /// <summary>
        /// The uniform(0,1) random number generator required to compute the multivariate CDF for dimensions greater than 2.
        /// </summary>
        public Random MVNUNI
        {
           get { return _MVNUNI; }
           set { _MVNUNI = value; }
        }

        /// <summary>
        /// The maximum number of function evaluations allowed when computing the multivariate CDF. Default = 1,000 x D.
        /// </summary>
        public int MaxEvaluations
        {
            get { return _maxEvaluations; }
            set { _maxEvaluations = value; }
        }

        /// <summary>
        /// The absolute error tolerance when computing the multivariate CDF. Default = 1E-3.
        /// </summary>
        public double AbsoluteError
        {
            get { return _absoluteError; }
            set { _absoluteError = value; }
        }

        /// <summary>
        /// The relative error tolerance when computing the multivariate CDF. Default = 1E-8.
        /// </summary>
        public double RelativeError
        {
            get { return _relativeError; }
            set { _relativeError = value; }
        }

        /// <summary>
        /// Gets the number of variables for the distribution.
        /// </summary>
        public override int Dimension {
            get { return _dimension; }
        }

        /// <summary>
        /// Returns the multivariate distribution type.
        /// </summary>
        public override MultivariateDistributionType Type
        {
            get { return MultivariateDistributionType.MultivariateNormal; }
        }

        /// <summary>
        /// Returns the display name of the distribution type as a string.
        /// </summary>
        public override string DisplayName
        {
            get { return "Multivariate Normal"; }
        }

        /// <summary>
        /// Returns the short display name of the distribution as a string.
        /// </summary>
        public override string ShortDisplayName
        {
            get { return "Multi-N"; }
        }

        /// <summary>
        /// Determines whether the parameters are valid or not.
        /// </summary>
        public override bool ParametersValid
        {
            get { return _parametersValid; }
        }

        /// <summary>
        /// Gets the mean vector of the distribution.
        /// </summary>
        public double[] Mean
        {
            get { return _mean; }
        }

        /// <summary>
        /// Gets the median vector of the distribution.
        /// </summary>
        public double[] Median
        {
            get { return _mean; }
        }

        /// <summary>
        /// Gets the mode vector of the distribution.
        /// </summary>
        public double[] Mode
        {
            get { return _mean; }
        }

        /// <summary>
        /// Gets the variance vector of the distribution.
        /// </summary>
        public double[] Variance
        {
            get
            {
                if (_variance == null)
                    _variance = _covariance.Diagonal();
                return _variance;
            }
        }

        /// <summary>
        /// Gets the standard deviation vector of the distribution.
        /// </summary>
        public double[] StandardDeviation
        {
            get
            {
                if (_standardDeviation == null)
                {
                    _standardDeviation = new double[Dimension];
                    for (int i = 0; i < Dimension; i++)
                        _standardDeviation[i] = Math.Sqrt(Variance[i]);
                }
                return _standardDeviation;
            }
        }

        /// <summary>
        /// Gets the Variance-Covariance matrix for the distribution.
        /// </summary>
        public double[,] Covariance
        {
            get { return _covariance.ToArray(); }
        }

        /// <summary>
        /// Determines if the covariance matrix is positive definite.
        /// </summary>
        public bool IsPositiveDefinite => _cholesky.IsPositiveDefinite;

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="mean">The mean vector μ (mu) for the distribution.</param>
        /// <param name="covariance">The covariance matrix Σ (sigma) for the distribution.</param>
        public void SetParameters(double[] mean, double[,] covariance)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(mean, covariance, true);

            _dimension = mean.Length;      
            _mean = mean;
            _covariance = new Matrix(covariance);
            _cholesky = new CholeskyDecomposition(_covariance);
            double lndet = _cholesky.LogDeterminant();
            _lnconstant = -(Math.Log(2d * Math.PI) * _mean.Length + lndet) * 0.5d;

            // Set up parameters for MVN CDF
            _correlationMatrixCreated = false;
            _maxEvaluations = 1000 * _dimension;
            _lower = new double[Dimension];
            _upper = new double[Dimension];
            _infin = new int[Dimension];
            for (int i = 0; i < Dimension; i++)
            {
                _lower[i] = double.NegativeInfinity;
                _infin[i] = 0;
            }
            // Set up arrays
            NL = (int)(Dimension * (Dimension - 1) / 2d);
            // Create empty arrays
            _correl = new double[NL];
            _correlation = new Matrix(Dimension, Dimension);

        }

        /// <summary>
        /// Create the correlation arrays required for computing the CDF.
        /// </summary>
        private void CreateCorrelationMatrix()
        {
            // Save correlation matrix
            var D = Matrix.Diagonal(_covariance);
            D.Sqrt();
            var invSqrtD = !D;
            _correlation = (invSqrtD * _covariance) * invSqrtD;
            // Create collapsed matrix
            _correl = new double[NL];
            int t = 0;
            for (int i = 1; i < Dimension; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    _correl[t] = _correlation[i, j];
                    t++;
                }
            }
            _correlationMatrixCreated = true;
        }


        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="mean">The mean vector μ (mu) for the distribution.</param>
        /// <param name="covariance">The covariance matrix Σ (sigma) for the distribution.</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double[] mean, double[,] covariance, bool throwException)
        {
            // Check if matrix is square
            if (_covariance.IsSquare == false)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Covariance), "Covariance matrix must be square.");
                return new ArgumentOutOfRangeException(nameof(Covariance), "Covariance matrix must be square.");
            }
            // 
            // Check if cholesky decomp is positive definite
            if (_cholesky.IsPositiveDefinite == false)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Covariance), "Covariance matrix is not positive-definite.");
                return new ArgumentOutOfRangeException(nameof(Covariance), "Covariance matrix is not positive-definite.");
            }
            // 
            return null;
        }

        /// <summary>
        /// The Probability Density Function (PDF) of the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        public override double PDF(double[] x)
        {
            return Math.Exp(-0.5d * Mahalanobis(x) + _lnconstant);
        }

        /// <summary>
        /// Returns the natural log of the PDF.
        /// </summary>
        /// <param name="x">The vector of x values.</param>
        public override double LogPDF(double[] x)
        {
            double f = -0.5d * Mahalanobis(x) + _lnconstant;
            if (double.IsNaN(f) || double.IsInfinity(f)) return double.MinValue;
            return f;
        }

        /// <summary>
        /// Gets the Mahalanobis distance between a sample and this distribution.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        public double Mahalanobis(double[] x)
        {
            if (x.Length != Dimension)
                throw new ArgumentOutOfRangeException(nameof(x), "The vector must be the same dimension as the distribution.");
            // 
            var z = new double[_mean.Length];
            for (int i = 0; i < x.Length; i++)
                z[i] = x[i] - _mean[i];
            var a = _cholesky.Solve(new Vector(z));
            double b = 0d;
            for (int i = 0; i < z.Length; i++)
                b += a[i] * z[i];
            return b;
        }

        /// <summary>
        /// The Cumulative Distribution Function (CDF) for the distribution evaluated at a point X.
        /// </summary>
        /// <param name="x">A point in the distribution space.</param>
        public override double CDF(double[] x)
        {
            if (Dimension == 1)
            {
                double stdDev = Math.Sqrt(Covariance[0, 0]);
                if (stdDev == 0d)
                    return x[0] == _mean[0] ? 1 : 0;
                double z = (x[0] - _mean[0]) / stdDev;
                return Normal.StandardCDF(z);
            }
            else if (Dimension == 2)
            {
                double sigma1 = Math.Sqrt(Covariance[0, 0]);
                double sigma2 = Math.Sqrt(Covariance[1, 1]);
                double rho = Covariance[0, 1] / (sigma1 * sigma2);
                if (double.IsNaN(rho)) return x.Equals(_mean) ? 1 : 0;
                double z = (x[0] - _mean[0]) / sigma1;
                double w = (x[1] - _mean[1]) / sigma2;
                return BivariateCDF(-z, -w, rho);
            }
            else
            {
                // Construct inputs for GenzBretz method
                for (int i = 0; i < x.Length; i++)
                    _upper[i] = (x[i] - _mean[i]) / Math.Sqrt(Covariance[i, i]);

                double ERROR = 0;
                double VAL = 0;
                int INFORM = 0;

                if (_correlationMatrixCreated == false)
                    CreateCorrelationMatrix();
                MVNDST(Dimension, _lower, _upper, _infin, _correl, _maxEvaluations, _absoluteError, _relativeError, ref ERROR, ref VAL, ref INFORM);
                return VAL;
            }
        }

        public double Interval(double[] lower, double[] upper)
        {
            // Construct inputs for GenzBretz method
            var infin = new int[upper.Length];
            for (int i = 0; i < upper.Length; i++)
            {
                lower[i] = (lower[i] - _mean[i]) / Math.Sqrt(Covariance[i, i]);
                upper[i] = (upper[i] - _mean[i]) / Math.Sqrt(Covariance[i, i]);
                infin[i] = 2;
            }
                
            double ERROR = 0;
            double VAL = 0;
            int INFORM = 0;

            if (_correlationMatrixCreated == false)
                CreateCorrelationMatrix();
            MVNDST(Dimension, lower, upper, infin, _correl, _maxEvaluations, _absoluteError, _relativeError, ref ERROR, ref VAL, ref INFORM);
            return VAL;
        }

        /// <summary>
        /// The inverse cumulative distribution function (InverseCDF).
        /// </summary>
        /// <param name="probabilities">Array of probabilities.</param>
        public double[] InverseCDF(double[] probabilities)
        {
            var sample = new double[Dimension];
            // Create vector z of standard normal variates for each dimension
            var z = new double[Dimension];
            for (int j = 0; j < Dimension; j++)
                z[j] = Normal.StandardZ(probabilities[j]);
            // x = A*z + mu
            var Az = _cholesky.L * z;
            for (int j = 0; j < Dimension; j++)
                sample[j] = Az[j] + _mean[j];
            return sample;
        }

        /// <summary>
        /// Returns an independent univariate Normal distribution for the given index.
        /// </summary>
        /// <param name="index">The zero-based index of the distribution.</param>
        public Normal IndependentNormal(int index)
        {
            return new Normal(Mean[index], StandardDeviation[index]);
        }

        /// <summary>
        /// Returns a new univariate Normal distribution.
        /// </summary>
        /// <param name="mu">Mean of the distribution.</param>
        /// <param name="sigma">Standard deviation of the distribution.</param>
        public static MultivariateNormal Univariate(double mu, double sigma)
        {
            return new MultivariateNormal(new[] { mu }, new[,] { { sigma * sigma } });
        }

        /// <summary>
        /// Returns a new bivariate Normal distribution.
        /// </summary>
        /// <param name="mu1">Mean of the first variate of the distribution.</param>
        /// <param name="mu2">Mean of the second variate of the distribution.</param>
        /// <param name="sigma1">Standard deviation of the first variate of the distribution.</param>
        /// <param name="sigma2">Standard deviation of the second variate of the distribution.</param>
        /// <param name="rho">The correlation coefficient between the two distributions.</param>
        public static MultivariateNormal Bivariate(double mu1, double mu2, double sigma1, double sigma2, double rho)
        {
            var mu = new double[] { mu1, mu2 };
            var covariance = new double[,] { { sigma1 * sigma1, sigma1 * sigma2 * rho }, { sigma1 * sigma2 * rho, sigma2 * sigma2 } };
            return new MultivariateNormal(mu, covariance);
        }

        /// <summary>
        /// Generate random values of a distribution given a sample size.
        /// </summary>
        /// <param name="sampleSize"> Size of random sample to generate. </param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <returns>
        /// Array of random values. The number of rows are equal to the sample size.
        /// The number of columns are equal to the dimensions of this distribution.
        /// </returns>
        public double[,] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            // Create PRNG for generating random numbers
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize, Dimension];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
            {
                // Create vector z of standard normal variates for each dimension
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                    z[j] = Normal.StandardZ(rnd.NextDouble());
                // x = A*z + mu
                var Az = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = Az[j] + _mean[j];
            }
            // Return array of random values
            return sample;
        }

        /// <summary>
        /// Use Latin hypercube method to generate random values of a distribution given a sample size and a user-defined seed.
        /// </summary>
        /// <param name="sampleSize"> Size of random sample to generate. </param>
        /// <param name="seed"> Seed for random number generator. </param>
        /// <returns>
        /// Array of random values.
        /// </returns>
        public double[,] LatinHypercubeRandomValues(int sampleSize, int seed)
        {
            var r = LatinHypercube.Random(sampleSize, Dimension, seed);
            var sample = new double[sampleSize, Dimension];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
            {
                // Create vector z of standard normal variates for each dimension
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                    z[j] = Normal.StandardZ(r[i, j]);
                // x = A*z + mu
                var Az = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = Az[j] + _mean[j];
            }
            // Return array of random values
            return sample;
        }

        /// <summary>
        /// Returns a 2D array of stratified z-variates. The first dimension is stratified, and the remaining are sampled randomly. 
        /// </summary>
        /// <param name="stratificationBins">A list of stratification bins.</param>
        /// <param name="seed"> Seed for random number generator. </param>
        public double[,] StratifiedRandomValues(List<StratificationBin> stratificationBins, int seed)
        {
            int samplesize = stratificationBins.Count;
            var r = new MersenneTwister(seed);
            var sample = new double[samplesize, Dimension];
            // Generate values
            for (int i = 0; i < samplesize; i++)
            {
                // Create vector z of standard normal variates for each dimension
                var z = new double[Dimension];
                for (int j = 0; j < Dimension; j++)
                {
                    if (j == 0)
                    {
                        z[j] = Normal.StandardZ(stratificationBins[i].Midpoint);
                    }
                    else
                    {
                        z[j] = Normal.StandardZ(r.NextDouble());
                    }               
                }               
                // x = A*z + mu
                var Az = _cholesky.L * z;
                for (int j = 0; j < Dimension; j++)
                    sample[i, j] = Az[j] + _mean[j];
            }
            // Return array of random values
            return sample;
        }

        #region Cumulative Distribution Support

    
        //****************************************************************************80

        public static double bivnor(double ah, double ak, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVNOR computes the bivariate normal CDF.
        //
        //  Discussion:
        //
        //    BIVNOR computes the probability for two normal variates X and Y
        //    whose correlation is R, that AH <= X and AK <= Y.
        //
        //  Licensing:
        //
        //    This code is distributed under the MIT license.
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Thomas Donnelly.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Thomas Donnelly,
        //    Algorithm 462: Bivariate Normal Distribution,
        //    Communications of the ACM,
        //    October 1973, Volume 16, Number 10, page 638.
        //
        //  Parameters:
        //
        //    Input, double AH, AK, the lower limits of integration.
        //
        //    Input, double R, the correlation between X and Y.
        //
        //    Output, double BIVNOR, the bivariate normal CDF.
        //
        //  Local Parameters:
        //
        //    Local, int IDIG, the number of significant digits
        //    to the right of the decimal point desired in the answer.
        //
        {
            double a2;
            double ap;
            double b;
            double cn;
            double con;
            double conex;
            double ex;
            double g2;
            double gh;
            double gk;
            double gw = 0;
            double h2;
            double h4;
            int i;
            int idig = 15;
            int i5 = 0;
            double rr;
            double s1;
            double s2;
            double sgn;
            double sn;
            double sp;
            double sqr;
            double t;
            double twopi = 6.283185307179587;
            double w2=0;
            double wh = 0;
            double wk = 0;

            b = 0.0;

            gh = gauss(-ah) / 2.0;
            gk = gauss(-ak) / 2.0;

            if (r == 0.0)
            {
                b = 4.00 * gh * gk;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }

            rr = (1.0 + r) * (1.0 - r);

            if (rr < 0.0)
            {
                //cerr << "\n";
                //cerr << "BIVNOR - Fatal error!\n";
                //cerr << "  1 < |R|.\n";
                //exit(0);
            }

            if (rr == 0.0)
            {
                if (r < 0.0)
                {
                    if (ah + ak < 0.0)
                    {
                        b = 2.0 * (gh + gk) - 1.0;
                    }
                }
                else
                {
                    if (ah - ak < 0.0)
                    {
                        b = 2.0 * gk;
                    }
                    else
                    {
                        b = 2.0 * gh;
                    }
                }
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }

            sqr = Math.Sqrt(rr);

            if (idig == 15)
            {
                con = twopi * 1.0E-15 / 2.0;
            }
            else
            {
                con = twopi / 2.0;
                for (i = 1; i <= idig; i++)
                {
                    con = con / 10.0;
                }
            }
            //
            //  (0,0)
            //
            if (ah == 0.0 && ak == 0.0)
            {
                b = 0.25 + Math.Asin(r) / twopi;
                b = Math.Max(b, 0.0);
                b = Math.Min(b, 1.0);
                return b;
            }
            //
            //  (0,nonzero)
            //
            if (ah == 0.0 && ak != 0.0)
            {
                b = gk;
                wh = -ak;
                wk = (ah / ak - r) / sqr;
                gw = 2.0 * gk;
                i5 = 1;
            }
            //
            //  (nonzero,0)
            //
            else if (ah != 0.0 && ak == 0.0)
            {
                b = gh;
                wh = -ah;
                wk = (ak / ah - r) / sqr;
                gw = 2.0 * gh;
                i5 = -1;
            }
            //
            //  (nonzero,nonzero)
            //
            else if (ah != 0.0 && ak != 0.0)
            {
                b = gh + gk;
                if (ah * ak < 0.0)
                {
                    b = b - 0.5;
                }
                wh = -ah;
                wk = (ak / ah - r) / sqr;
                gw = 2.0 * gh;
                i5 = -1;
            }

            for (; ; )
            {
                sgn = -1.0;
                t = 0.0;

                if (wk != 0.0)
                {
                    if (Math.Abs(wk) == 1.0)
                    {
                        t = wk * gw * (1.0 - gw) / 2.0;
                        b = b + sgn * t;
                    }
                    else
                    {
                        if (1.0 < Math.Abs(wk))
                        {
                            sgn = -sgn;
                            wh = wh * wk;
                            g2 = gauss(wh);
                            wk = 1.0 / wk;

                            if (wk < 0.0)
                            {
                                b = b + 0.5;
                            }
                            b = b - (gw + g2) / 2.0 + gw * g2;
                        }
                        h2 = wh * wh;
                        a2 = wk * wk;
                        h4 = h2 / 2.0;
                        ex = Math.Exp(-h4);
                        w2 = h4 * ex;
                        ap = 1.0;
                        s2 = ap - ex;
                        sp = ap;
                        s1 = 0.0;
                        sn = s1;
                        conex = Math.Abs(con / wk);

                        for (; ; )
                        {
                            cn = ap * s2 / (sn + sp);
                            s1 = s1 + cn;

                            if (Math.Abs(cn) <= conex)
                            {
                                break;
                            }
                            sn = sp;
                            sp = sp + 1.0;
                            s2 = s2 - w2;
                            w2 = w2 * h4 / sp;
                            ap = -ap * a2;
                        }
                        t = (Math.Atan(wk) - wk * s1) / twopi;
                        b = b + sgn * t;
                    }
                }
                if (0 <= i5 )
                {
                    break;
                }
                if (ak == 0.0)
                {
                    break;
                }
                wh = -ak;
                wk = (ah / ak - r) / sqr;
                gw = 2.0 * gk;
                i5 = 1;
            }

            b = Math.Max(b, 0.0);
            b = Math.Min(b, 1.0);

            return b;
        }
        //****************************************************************************80

        private static double gauss(double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAUSS returns the area of the lower tail of the normal curve.
        //
        //  Licensing:
        //
        //    This code is distributed under the MIT license.
        //
        //  Modified:
        //
        //    13 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the evaluation point.
        //
        //    Output, double GAUSS, the lower normal tail area.
        //
        {
            double value;

            value = (1.0 + Erf.Function(t / Math.Sqrt(2.0))) / 2.0;

            return value;
        }


        /// <summary>
        /// A function for computing the bivariate normal CDF.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This code was copied from Accord Math Library, http://accord-framework.net
        /// </para>
        /// <para>
        /// This method is based on the work done by Alan Genz, Department of
        /// Mathematics, Washington State University. Pullman, WA 99164-3113
        /// Email: alangenz@wsu.edu. This work was shared under a 3-clause BSD
        /// license. Please see source file for more details and the actual
        /// license text.</para>
        /// <para>
        /// This function is based on the method described by Drezner, Z and G.O.
        /// Wesolowsky, (1989), On the computation of the bivariate normal integral,
        /// Journal of Statist. Comput. Simul. 35, pp. 101-107, with major modifications
        /// for double precision, and for |R| close to 1.</para>
        /// </remarks>
        /// <param name="z1">The z variate for the first Normal.</param>
        /// <param name="z2">The z variate for the second Normal.</param>
        /// <param name="r">The correlation coefficient.</param>
        public static double BivariateCDF(double z1, double z2, double r)
        {

            // Copyright (C) 2013, Alan Genz,  All rights reserved.               
            // 
            // Redistribution and use in source and binary forms, with or without
            // modification, are permitted provided the following conditions are met:
            // 1. Redistributions of source code must retain the above copyright
            // notice, this list of conditions and the following disclaimer.
            // 2. Redistributions in binary form must reproduce the above copyright
            // notice, this list of conditions and the following disclaimer in 
            // the documentation and/or other materials provided with the 
            // distribution.
            // 3. The contributor name(s) may not be used to endorse or promote 
            // products derived from this software without specific prior 
            // written permission.
            // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
            // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
            // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
            // FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
            // COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
            // INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
            // BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
            // OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
            // ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
            // TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
            // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

            const double TWOPI = 2.0d * Math.PI;
            double[] x;
            double[] w;
            if (Math.Abs(r) < 0.3d)
            {
                // Gauss Legendre Points and Weights N =  6
                x = BVND_XN6;
                w = BVND_WN6;
            }
            else if (Math.Abs(r) < 0.75d)
            {
                // Gauss Legendre Points and Weights N =  12
                x = BVND_XN12;
                w = BVND_WN12;
            }
            else
            {
                // Gauss Legendre Points and Weights N =  20
                x = BVND_XN20;
                w = BVND_WN20;
            }

            double h = z1;
            double k = z2;
            double hk = h * k;
            double bvn = 0d;
            if (Math.Abs(r) < 0.925d)
            {
                if (Math.Abs(r) > 0d)
                {
                    double sh = (h * h + k * k) / 2d;
                    double asr = Math.Asin(r);
                    for (int i = 0; i < x.Length;  i++)
                    {
                        for (int j = -1; j <= 1; j += 2)
                        {
                            double sn = Math.Sin(asr * (j * x[i] + 1d) / 2d);
                            bvn = bvn + w[i] * Math.Exp((sn * hk - sh) / (1d - sn * sn));
                        }
                    }

                    bvn = bvn * asr / (2d * TWOPI);
                }

                return bvn + MVNPHI(-h) * MVNPHI(-k);
            }

            if (r < 0d)
            {
                k = -k;
                hk = -hk;
            }

            if (Math.Abs(r) < 1d)
            {
                double sa = (1d - r) * (1d + r);
                double A = Math.Sqrt(sa);
                double sb = h - k;
                sb = sb * sb;
                double c = (4d - hk) / 8d;
                double d = (12d - hk) / 16d;
                double asr = -(sb / sa + hk) / 2d;
                if (asr > -100)
                    bvn = A * Math.Exp(asr) * (1d - c * (sb - sa) * (1d - d * sb / 5d) / 3d + c * d * sa * sa / 5d);
                if (-hk < 100d)
                {
                    double B = Math.Sqrt(sb);
                    bvn = bvn - Math.Exp(-hk / 2d) * Math.Sqrt(TWOPI) * MVNPHI(-B / A) * B * (1d - c * sb * (1d - d * sb / 5d) / 3d);
                }

                A = A / 2d;
                for (int i = 0; i < x.Length; i++)
                {
                    for (int j = -1; j <= 1; j += 2)
                    {
                        double xs = A * (j * x[i] + 1d);
                        xs = xs * xs;
                        double rs = Math.Sqrt(1d - xs);
                        asr = -(sb / xs + hk) / 2d;
                        if (asr > -100)
                        {
                            bvn = bvn + A * w[i] * Math.Exp(asr) * (Math.Exp(-hk * xs / (2d * (1d + rs) * (1d + rs))) / rs - (1d + c * xs * (1d + d * xs)));
                        }
                    }
                }

                bvn = -bvn / TWOPI;
            }

            if (r > 0d)
                return bvn + MVNPHI(-Math.Max(h, k));
            bvn = -bvn;
            if (k <= h)
                return bvn;
            if (h < 0d)
                return bvn + MVNPHI(k) - MVNPHI(h);
            return bvn + MVNPHI(-h) - MVNPHI(-k);
        }

        private static readonly double[] BVND_WN20 = new[] { 0.017614007139152121d, 0.040601429800386939d, 0.062672048334109054d, 0.083276741576704755d, 0.1019301198172404d, 0.1181945319615184d, 0.13168863844917661d, 0.1420961093183821d, 0.14917298647260371d, 0.15275338713072589d };
        private static readonly double[] BVND_XN20 = new[] { -0.99312859918509488d, -0.96397192727791381d, -0.912234428251326d, -0.83911697182221878d, -0.7463319064601508d, -0.636053680726515d, -0.51086700195082713d, -0.3737060887154196d, -0.2277858511416451d, -0.076526521133497324d };
        private static readonly double[] BVND_WN12 = new[] { 0.047175336386511772d, 0.1069393259953183d, 0.16007832854334639d, 0.2031674267230659d, 0.2334925365383547d, 0.24914704581340291d };
        private static readonly double[] BVND_XN12 = new[] { -0.98156063424671913d, -0.904117256370475d, -0.769902674194305d, -0.58731795428661715d, -0.36783149899818018d, -0.12523340851146919d };
        private static readonly double[] BVND_WN6 = new[] { 0.1713244923791705d, 0.36076157304813838d, 0.46791393457269043d };
        private static readonly double[] BVND_XN6 = new[] { -0.93246951420315216d, -0.6612093864662647d, -0.238619186083197d };


        private static double[] MVNPHI_A = new double[]
        {
            6.10143081923200417926465815756E-1,
            -4.34841272712577471828182820888E-1,
            1.76351193643605501125840298123E-1,
            -6.0710795609249414860051215825E-2,
            1.7712068995694114486147141191E-2,
            -4.321119385567293818599864968E-3,
            8.54216676887098678819832055E-4,
            -1.27155090609162742628893940E-4,
            1.1248167243671189468847072E-5,
            3.13063885421820972630152E-7,
            -2.70988068537762022009086E-7,
            3.0737622701407688440959E-8,
            2.515620384817622937314E-9,
            -1.028929921320319127590E-9,
            2.9944052119949939363E-11,
            2.6051789687266936290E-11,
            -2.634839924171969386E-12,
            -6.43404509890636443E-13,
            1.12457401801663447E-13,
            1.7281533389986098E-14,
            -4.264101694942375E-15,
            -5.45371977880191E-16,
            1.58697607761671E-16,
            2.0899837844334E-17,
            -5.900526869409E-18,
            -9.41893387554E-19,
            2.14977356470E-19,
            4.6660985008E-20,
            -7.243011862E-21,
            -2.387966824E-21,
            1.91177535E-22,
            1.20482568E-22,
            -6.72377E-25,
            -5.747997E-24,
            -4.28493E-25,
            2.44856E-25,
            4.3793E-26,
            -8.151E-27,
            -3.089E-27,
            9.3E-29,
            1.74E-28,
            1.6E-29,
            -8.0E-30,
            -2.0E-30
        };

        /// <summary>
        /// Gauss Legendre Points and Weights, N = 6, 12, and 20
        /// </summary>
        private static double[,] X = new double[10, 3]
        {
            {-0.932469514203152,-0.9815606342467192,-0.9931285991850949},
            {-0.6612093864662645,-0.9041172563704749,-0.9639719272779138},
            {-0.2386191860831969,-0.7699026741943047,-0.912234428251326},
            {0,-0.5873179542866175,-0.8391169718222188},
            {0,-0.3678314989981802,-0.7463319064601508},
            {0,-0.1252334085114689,-0.636053680726515},
            {0,0,-0.5108670019508271},
            {0,0,-0.37370608871541955},
            {0,0,-0.22778585114164507},
            {0,0,-0.07652652113349734}
        };

        private static double[,] W = new double[10, 3]
        {
            {0.17132449237917036,0.04717533638651183,0.0176140071391521},
            {0.3607615730481386,0.10693932599531843,0.0406014298003869},
            {0.46791393457269104,0.16007832854334622,0.062672048334109},
            {0,0.20316742672306592,0.0832767415767047},
            {0,0.2334925365383548,0.10193011981724},
            {0,0.24914704581340277,0.118194531961518},
            {0,0,0.131688638449176},
            {0,0,0.142096109318382},
            {0,0,0.149172986472603},
            {0,0,0.152753387130725}
        };

        //SAVE A, B, INFI, COV
        // NL = 500
        private int NL = 500;
        private double[] MVNDFN_A = new double[500];
        private double[] MVNDFN_B = new double[500];
        private int[] MVNDFN_INFI = new int[500];
        private double[] MVNDFN_COV = new double[125250]; //NL*(NL+1)/2

        // Optimal Parameters for Lattice Rules
        private int[] DKBVRC_P = new int[28] { 31, 47, 73, 113, 173, 263, 397, 593, 907, 1361, 2053, 3079, 4621, 6947, 10427, 15641, 23473, 35221, 52837, 79259, 118891, 178349, 267523, 401287, 601943, 902933, 1354471, 2031713 };
        private int[,] DKBVRC_C = new int[28, 99]
                {
                    {12, 9, 9, 13, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3, 12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3, 12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3, 12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 12, 3, 3, 3, 12, 7, 7, 12, 12, 12, 12, 12, 12, 12, 12, 7, 3, 3, 3, 7, 7, 7, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,},
                    {13, 11, 17, 10, 15, 15, 15, 15, 15, 15, 22, 15, 15, 6, 6, 6, 15, 15, 9, 13, 2, 2, 2, 13, 11, 11, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 6, 6, 6, 15, 15, 9, 13, 2, 2, 2, 13, 11, 11, 10, 15, 15, 15, 15, 15, 15, 15, 15, 15, 6, 6, 6, 15, 15, 9, 13, 2, 2, 2, 13, 11, 11, 10, 10, 15, 15, 15, 15, 15, 15, 15, 15, 6, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,},
                    {27, 28, 10, 11, 11, 20, 11, 11, 28, 13, 13, 28, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 31, 31, 5, 5, 5, 31, 13, 11, 11, 11, 11, 11, 11, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 31, 31, 5, 5, 5, 11, 13, 11, 11, 11, 11, 11, 11, 11, 13, 13, 11, 13, 5, 5, 5, 5, 14, 13, 5, 5, 5, 5, 5, 5, 5, 5,},
                    {35, 27, 27, 36, 22, 29, 29, 20, 45, 5, 5, 5, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 29, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 21, 27, 3, 3, 3, 24, 27, 27, 17, 29, 29, 29, 17, 5, 5, 5, 5, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 17, 17, 17, 6, 17, 17, 6, 3, 6, 6, 3, 3, 3, 3, 3,},
                    { 64, 66, 28, 28, 44, 44, 55, 67, 10, 10, 10, 10, 10, 10, 38, 38, 10, 10, 10, 10, 10, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 38, 38, 31, 4, 4, 31, 64, 4, 4, 4, 64, 45, 45, 45, 45, 45, 45, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 11, 66, 66, 66, 66, 66, 66, 66, 66, 66, 45, 11, 7, 3, 2, 2, 2, 27, 5, 3, 3, 5, 5, 2, 2, 2, 2, 2, 2, 2,},
                    { 111, 42, 54, 118, 20, 31, 31, 72, 17, 94, 14, 14, 11, 14, 14, 14, 94, 10, 10, 10, 10, 14, 14, 14, 14, 14, 14, 14, 11, 11, 11, 8, 8, 8, 8, 8, 8, 8, 18, 18, 18, 18, 18, 113, 62, 62, 45, 45, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 63, 63, 53, 63, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 51, 51, 51, 51, 51, 12, 51, 12, 51, 5, 3, 3, 2, 2, 5,},
                    { 163, 154, 83, 43, 82, 92, 150, 59, 76, 76, 47, 11, 11, 100, 131, 116, 116, 116, 116, 116, 116, 138, 138, 138, 138, 138, 138, 138, 138, 138, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 116, 116, 116, 116, 116, 116, 100, 100, 100, 100, 100, 138, 138, 138, 138, 138, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 38, 38, 38, 38, 38, 38, 38, 38, 3, 3, 3, 3, 3,},
                    { 246, 189, 242, 102, 250, 250, 102, 250, 280, 118, 196, 118, 191, 215, 121, 121, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 161, 161, 161, 161, 161, 161, 161, 161, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 10, 10, 10, 10, 10, 10, 103, 10, 10, 10, 10, 5,},
                    { 347, 402, 322, 418, 215, 220, 339, 339, 339, 337, 218, 315, 315, 315, 315, 167, 167, 167, 167, 361, 201, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 231, 231, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 48, 48, 48, 48, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 243, 243, 243, 243, 243, 243, 243, 243, 243, 243, 283, 283, 283, 283, 283, 283, 283, 283, 283, 16, 283, 16, 283, 283,},
                    { 505, 220, 601, 644, 612, 160, 206, 206, 206, 422, 134, 518, 134, 134, 518, 652, 382, 206, 158, 441, 179, 441, 56, 559, 559, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 101, 101, 56, 101, 101, 101, 101, 101, 101, 101, 101, 193, 193, 193, 193, 193, 193, 193, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 101, 101, 101, 101,},
                    {794, 325, 960, 528, 247, 247, 338, 366, 847, 753, 753, 236, 334, 334, 461, 711, 652, 381, 381, 381, 652, 381, 381, 381, 381, 381, 381, 381, 226, 326, 326, 326, 326, 326, 326, 326, 126, 326, 326, 326, 326, 326, 326, 326, 326, 326, 326, 195, 195, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 195, 195, 195, 195, 195, 195, 195, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 387, 387, 387, 387, 387, 387, 387, 387, 387, 387, 387, 387, 387,},
                    { 1189, 888, 259, 1082, 725, 811, 636, 965, 497, 497, 1490, 1490, 392, 1291, 508, 508, 1291, 1291, 508, 1291, 508, 508, 867, 867, 867, 867, 934, 867, 867, 867, 867, 867, 867, 867, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 1284, 563, 563, 563, 563, 1010, 1010, 1010, 208, 838, 563, 563, 563, 759, 759, 564, 759, 759, 801, 801, 801, 801, 759, 759, 759, 759, 759, 563, 563, 563, 563, 563, 563, 563, 563, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226, 226,},
                    { 1763, 1018, 1500, 432, 1332, 2203, 126, 2240, 1719, 1284, 878, 1983, 266, 266, 266, 266, 747, 747, 127, 127, 2074, 127, 2074, 1400, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 1400, 1383, 1383, 1383, 1383, 1383, 1383, 1383, 507, 1073, 1073, 1073, 1073, 1990, 1990, 1990, 1990, 1990, 507, 507, 507, 507, 507, 507, 507, 507, 507, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 1073, 22, 22, 22, 22, 22, 22, 1073, 452, 452, 452, 452, 452, 452, 318, 301, 301, 301, 301, 86, 86, 15,},
                    { 2872, 3233, 1534, 2941, 2910, 393, 1796, 919, 446, 919, 919, 1117, 103, 103, 103, 103, 103, 103, 103, 2311, 3117, 1101, 3117, 3117, 1101, 1101, 1101, 1101, 1101, 2503, 2503, 2503, 2503, 2503, 2503, 2503, 2503, 429, 429, 429, 429, 429, 429, 429, 1702, 1702, 1702, 184, 184, 184, 184, 184, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 784, 784, 784, 784, 784, 784, 784, 784, 784, 784, 784, 784, 784,},
                    { 4309, 3758, 4034, 1963, 730, 642, 1502, 2246, 3834, 1511, 1102, 1102, 1522, 1522, 3427, 3427, 3928, 915, 915, 3818, 3818, 3818, 3818, 4782, 4782, 4782, 3818, 4782, 3818, 3818, 1327, 1327, 1327, 1327, 1327, 1327, 1327, 1387, 1387, 1387, 1387, 1387, 1387, 1387, 1387, 1387, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 2339, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 3148, 1776, 1776, 1776, 3354, 3354, 3354, 925, 3354, 3354, 925, 925, 925, 925, 925, 2133, 2133, 2133, 2133, 2133, 2133, 2133, 2133,},
                    {6610, 6977, 1686, 3819, 2314, 5647, 3953, 3614, 5115, 423, 423, 5408, 7426, 423, 423, 487, 6227, 2660, 6227, 1221, 3811, 197, 4367, 351, 1281, 1221, 351, 351, 351, 7245, 1984, 2999, 2999, 2999, 2999, 2999, 2999, 3995, 2063, 2063, 2063, 2063, 1644, 2063, 2077, 2512, 2512, 2512, 2077, 2077, 2077, 2077, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 754, 1097, 1097, 754, 754, 754, 754, 248, 754, 1097, 1097, 1097, 1097, 222, 222, 222, 222, 754, 1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982, 1982,},
                    { 9861, 3647, 4073, 2535, 3430, 9865, 2830, 9328, 4320, 5913, 10365, 8272, 3706, 6186, 7806, 7806, 7806, 8610, 2563, 11558, 11558, 9421, 1181, 9421, 1181, 1181, 1181, 9421, 1181, 1181, 10574, 10574, 3534, 3534, 3534, 3534, 3534, 2898, 2898, 2898, 3450, 2141, 2141, 2141, 2141, 2141, 2141, 2141, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 7055, 2831, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 8204, 4688, 4688, 4688, 2831, 2831, 2831, 2831, 2831, 2831, 2831, 2831,},
                    {10327, 7582, 7124, 8214, 9600, 10271, 10193, 10800, 9086, 2365, 4409, 13812, 5661, 9344, 9344, 10362, 9344, 9344, 8585, 11114, 13080, 13080, 13080, 6949, 3436, 3436, 3436, 13213, 6130, 6130, 8159, 8159, 11595, 8159, 3436, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 7096, 4377, 7096, 4377, 4377, 4377, 4377, 4377, 5410, 5410, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 4377, 440, 440, 1199, 1199, 1199,},
                    {19540, 19926, 11582, 11113, 24585, 8726, 17218, 419, 4918, 4918, 4918, 15701, 17710, 4037, 4037, 15808, 11401, 19398, 25950, 25950, 4454, 24987, 11719, 8697, 1452, 1452, 1452, 1452, 1452, 8697, 8697, 6436, 21475, 6436, 22913, 6434, 18497, 11089, 11089, 11089, 11089, 3036, 3036, 14208, 14208, 14208, 14208, 12906, 12906, 12906, 12906, 12906, 12906, 12906, 12906, 7614, 7614, 7614, 7614, 5021, 5021, 5021, 5021, 5021, 5021, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 10145, 4544, 4544, 4544, 4544, 4544, 4544, 8394, 8394, 8394, 8394,},
                    { 34566, 9579, 12654, 26856, 37873, 38806, 29501, 17271, 3663, 10763, 18955, 1298, 26560, 17132, 17132, 4753, 4753, 8713, 18624, 13082, 6791, 1122, 19363, 34695, 18770, 18770, 18770, 18770, 15628, 18770, 18770, 18770, 18770, 33766, 20837, 20837, 20837, 20837, 20837, 20837, 6545, 6545, 6545, 6545, 6545, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 30483, 30483, 30483, 30483, 30483, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 12138, 9305, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 11107, 9305, 9305,},
                    { 31929, 49367, 10982, 3527, 27066, 13226, 56010, 18911, 40574, 20767, 20767, 9686, 47603, 47603, 11736, 11736, 41601, 12888, 32948, 30801, 44243, 53351, 53351, 16016, 35086, 35086, 32581, 2464, 2464, 49554, 2464, 2464, 49554, 49554, 2464, 81, 27260, 10681, 2185, 2185, 2185, 2185, 2185, 2185, 2185, 18086, 18086, 18086, 18086, 18086, 17631, 17631, 18086, 18086, 18086, 37335, 37774, 37774, 37774, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 26401, 12982, 40398, 40398, 40398, 40398, 40398, 40398, 3518, 3518, 3518, 37799, 37799, 37799, 37799, 37799, 37799, 37799, 37799, 37799, 4721, 4721, 4721, 4721, 7067, 7067, 7067, 7067,},
                    { 40701, 69087, 77576, 64590, 39397, 33179, 10858, 38935, 43129, 35468, 35468, 5279, 61518, 61518, 27945, 70975, 70975, 86478, 86478, 20514, 20514, 73178, 73178, 43098, 43098, 4701, 59979, 59979, 58556, 69916, 15170, 15170, 4832, 4832, 43064, 71685, 4832, 15170, 15170, 15170, 27679, 27679, 27679, 60826, 60826, 6187, 6187, 4264, 4264, 4264, 4264, 4264, 45567, 32269, 32269, 32269, 32269, 62060, 62060, 62060, 62060, 62060, 62060, 62060, 62060, 62060, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 1803, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 51108, 55315, 55315, 54140, 54140, 54140, 54140, 54140, 13134,},
                    { 103650, 125480, 59978, 46875, 77172, 83021, 126904, 14541, 56299, 43636, 11655, 52680, 88549, 29804, 101894, 113675, 48040, 113675, 34987, 48308, 97926, 5475, 49449, 6850, 62545, 62545, 9440, 33242, 9440, 33242, 9440, 33242, 9440, 62850, 9440, 9440, 9440, 90308, 90308, 90308, 47904, 47904, 47904, 47904, 47904, 47904, 47904, 47904, 47904, 41143, 41143, 41143, 41143, 41143, 41143, 41143, 36114, 36114, 36114, 36114, 36114, 24997, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 65162, 47650, 47650, 47650, 47650, 47650, 47650, 47650, 40586, 40586, 40586, 40586, 40586, 40586, 40586, 38725, 38725, 38725, 38725, 88329, 88329, 88329, 88329, 88329,},
                    { 165843, 90647, 59925, 189541, 67647, 74795, 68365, 167485, 143918, 74912, 167289, 75517, 8148, 172106, 126159, 35867, 35867, 35867, 121694, 52171, 95354, 113969, 113969, 76304, 123709, 123709, 144615, 123709, 64958, 64958, 32377, 193002, 193002, 25023, 40017, 141605, 189165, 189165, 141605, 189165, 189165, 141605, 141605, 141605, 189165, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, 80822, 80822, 80822, 80822, 80822, 80822, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, 7114, 131661,},
                    {130365, 236711, 110235, 125699, 56483, 93735, 234469, 60549, 1291, 93937, 245291, 196061, 258647, 162489, 176631, 204895, 73353, 172319, 28881, 136787, 122081, 122081, 275993, 64673, 211587, 211587, 211587, 282859, 282859, 211587, 242821, 256865, 256865, 256865, 122203, 291915, 122203, 291915, 291915, 122203, 25639, 25639, 291803, 245397, 284047, 245397, 245397, 245397, 245397, 245397, 245397, 245397, 94241, 66575, 66575, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453, 94453,},
                    { 333459, 375354, 102417, 383544, 292630, 41147, 374614, 48032, 435453, 281493, 358168, 114121, 346892, 238990, 317313, 164158, 35497, 70530, 70530, 434839, 24754, 24754, 24754, 393656, 118711, 118711, 148227, 271087, 355831, 91034, 417029, 417029, 91034, 91034, 417029, 91034, 299843, 299843, 413548, 413548, 308300, 413548, 413548, 413548, 308300, 308300, 308300, 413548, 308300, 308300, 308300, 308300, 308300, 15311, 15311, 15311, 15311, 176255, 176255, 23613, 23613, 23613, 23613, 23613, 23613, 172210, 204328, 204328, 204328, 204328, 121626, 121626, 121626, 121626, 121626, 200187, 200187, 200187, 200187, 200187, 121551, 121551, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, 13942, 13942, 13942, 13942, 13942,},
                    {500884, 566009, 399251, 652979, 355008, 430235, 328722, 670680, 405585, 405585, 424646, 670180, 670180, 641587, 215580, 59048, 633320, 81010, 20789, 389250, 389250, 638764, 638764, 389250, 389250, 398094, 80846, 147776, 147776, 296177, 398094, 398094, 147776, 147776, 396313, 578233, 578233, 578233, 19482, 620706, 187095, 620706, 187095, 126467, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 241663, 321632, 23210, 23210, 394484, 394484, 394484, 78101, 78101, 78101, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, 277743, 277743, 277743, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259,},
                    { 858339, 918142, 501970, 234813, 460565, 31996, 753018, 256150, 199809, 993599, 245149, 794183, 121349, 150619, 376952, 809123, 809123, 804319, 67352, 969594, 434796, 969594, 804319, 391368, 761041, 754049, 466264, 754049, 754049, 466264, 754049, 754049, 282852, 429907, 390017, 276645, 994856, 250142, 144595, 907454, 689648, 687580, 687580, 687580, 687580, 978368, 687580, 552742, 105195, 942843, 768249, 307142, 307142, 307142, 307142, 880619, 880619, 880619, 880619, 880619, 880619, 880619, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 117185, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 60731, 178309, 178309, 178309, 178309, 74373, 74373, 74373, 74373, 74373, 74373, 74373, 74373, 214965, 214965, 214965,},
                };

        private int DKBVRC_SAMPLS;
        private int DKBVRC_NP;
        private double DKBVRC_VAREST;

        /// <summary>
        /// A subroutine for computing multivariate normal probabilities.
        ///     This subroutine uses an algorithm given in the paper
        ///     "Numerical Computation of Multivariate Normal Probabilities", in
        ///     J.of Computational and Graphical Stat., 1(1992), pp. 141-149, by
        /// Alan Genz
        /// Department of Mathematics
        ///          Washington State University
        /// Pullman, WA 99164-3113
        ///          Email : AlanGenz @wsu.edu
        /// </summary>
        /// <param name="N">the number of variables.</param>
        /// <param name="LOWER">array of lower integration limits.</param>
        /// <param name="UPPER">array of upper integration limits.</param>
        /// <param name="INFIN">array of integration limits flags:
        ///            if INFIN(I) &gt; 0, Ith limits are (-infinity, infinity);
        ///            if INFIN(I) = 0, Ith limits are(-infinity, UPPER(I)];
        ///            if INFIN(I) = 1, Ith limits are[LOWER(I), infinity);
        ///            if INFIN(I) = 2, Ith limits are[LOWER(I), UPPER(I)].</param>
        /// <param name="CORREL">array of correlation coefficients; the correlation coefficient in row I column J of the correlation matrix
        ///  should be stored in CORREL(J + ((I-2)*(I-1))/2 ), for J less than I. The correlation matrix must be positive semidefinite.</param>
        /// <param name="MAXPTS">maximum number of function values allowed. This parameter can be used to limit the time.A sensible strategy is to start with MAXPTS = 1000 * N, and then increase MAXPTS if ERROR is too large.</param>
        /// <param name="ABSEPS">absolute error tolerance.</param>
        /// <param name="RELEPS">relative error tolerance.</param>
        /// <param name="ERROR">estimated absolute error, with 99% confidence level.</param>
        /// <param name="VALUE">estimated value for the integral.</param>
        /// <param name="INFORM">termination status parameter: if INFORM = 0, normal completion with ERROR &lt; EPS; if INFORM = 1, completion with ERROR > EPS and MAXPTS function values used; increase MAXPTS to decrease ERROR; if INFORM = 2, N > 500 or N less than 1.</param>
        public void MVNDST(int N, double[] LOWER, double[] UPPER, int[] INFIN, double[] CORREL, int MAXPTS, double ABSEPS, double RELEPS, ref double ERROR, ref double VALUE, ref int INFORM)
        {
            int INFIS = 0, IVLS;
            double D = 0, E = 0;
            double[] Y = new double[500];

            if (N > 500 || N < 1)
            {
                INFORM = 2;
                VALUE = 0;
                ERROR = 1;
            }
            else
            {
                INFORM = (int)MVNDNT(N, CORREL, LOWER, UPPER, INFIN, ref INFIS, ref D, ref E, Y);
                if (N - INFIS == 0)
                {
                    VALUE = 1;
                    ERROR = 0;
                }
                else if (N - INFIS == 1)
                {
                    VALUE = E - D;
                    // This can sometimes produce very small, but negative values (approx. -5E-73)
                    if (VALUE < 0) VALUE = 0;
                    ERROR = 2E-16;
                }
                else
                {
                    // Call the lattice rule integration subroutine.
                    IVLS = 0;
                    DKBVRC(N - INFIS - 1, IVLS, MAXPTS, MVNDFN, ABSEPS, RELEPS, ref ERROR, ref VALUE, ref INFORM);
                }
            }
        }

        /// <summary>
        /// Integrand subroutine.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="W"></param>
        /// <returns></returns>
        private double MVNDFN(int N, double[] W)
        {

            double SUM, AI = 0, BI = 0, DI = 0, EI = 0;

            var Y = new double[NL];
            //double[] Y = new double[500];

            double result = 1;
            int INFA = 0;
            int INFB = 0;
            int IK = 0;
            int IJ = 0;
            for (int i = 0; i < N + 1; i++)
            {
                SUM = 0;
                for (int j = 0; j <= i - 1; j++)
                {
                    if (j < IK) { SUM += MVNDFN_COV[IJ] * Y[j]; }
                    IJ++;
                }

                if (MVNDFN_INFI[i] != 0)
                {
                    if (INFA == 1)
                    { AI = Math.Max(AI, MVNDFN_A[i] - SUM); }
                    else
                    {
                        AI = MVNDFN_A[i] - SUM;
                        INFA = 1;
                    }
                }
                if (MVNDFN_INFI[i] != 1)
                {
                    if (INFB == 1)
                    {
                        BI = Math.Min(BI, MVNDFN_B[i] - SUM);
                    }
                    else
                    {
                        BI = MVNDFN_B[i] - SUM;
                        INFB = 1;
                    }
                }
                IJ++;

                if (i == (N) || MVNDFN_COV[IJ + IK + 1] > 0)
                {
                    MVNLMS(AI, BI, 2 * INFA + INFB - 1, ref DI, ref EI);

                    if (DI >= EI)
                    {
                        result = 0;
                        break;
                    }
                    else
                    {
                        result = result * (EI - DI);
                        if (i <= N) { Y[IK] = PHINVS(DI + W[IK] * (EI - DI)); }
                        IK = IK + 1;
                        INFA = 0;
                        INFB = 0;
                    }
                }

            }
            return result;
        }

        private double MVNDNT(int N, double[] CORREL, double[] LOWER, double[] UPPER, int[] INFIN, ref int INFIS, ref double D, ref double E, double[] Y)
        {
            double result = 0;

            COVSRT(N, LOWER, UPPER, CORREL, INFIN, Y, ref INFIS, MVNDFN_A, MVNDFN_B, MVNDFN_COV, MVNDFN_INFI);
         
            if (N - INFIS == 1)
            {
                MVNLMS(MVNDFN_A[0], MVNDFN_B[0], MVNDFN_INFI[0], ref D, ref E);
            }
            else if (N - INFIS == 2)
            {
                if (Math.Abs(MVNDFN_COV[2]) > 0)
                {
                    D = Math.Sqrt(1 + Math.Pow(MVNDFN_COV[1], 2));
                    if (MVNDFN_INFI[1] != 0) { MVNDFN_A[1] = MVNDFN_A[1] / D; }
                    if (MVNDFN_INFI[1] != 1) { MVNDFN_B[1] = MVNDFN_B[1] / D; }
                    E = BVNMVN(MVNDFN_A, MVNDFN_B, MVNDFN_INFI, MVNDFN_COV[1] / D);
                    D = 0;
                }
                else
                {
                    if (MVNDFN_INFI[0] != 0)
                    { if (MVNDFN_INFI[1] != 0) { MVNDFN_A[0] = Math.Max(MVNDFN_A[0], MVNDFN_A[1]); } }
                    else
                    { if (MVNDFN_INFI[1] != 0) { MVNDFN_A[0] = MVNDFN_A[1]; } }

                    if (MVNDFN_INFI[0] != 1)
                    { if (MVNDFN_INFI[1] != 1) { MVNDFN_B[0] = Math.Min(MVNDFN_B[0], MVNDFN_B[1]); } }
                    else
                    { if (MVNDFN_INFI[1] != 1) { MVNDFN_B[0] = MVNDFN_B[1]; } }

                    if (MVNDFN_INFI[0] != MVNDFN_INFI[1]) { MVNDFN_INFI[0] = 2; }
                    MVNLMS(MVNDFN_A[0], MVNDFN_B[0], MVNDFN_INFI[0], ref D, ref E);
                }
                INFIS++;
            }

            return result;
        }

        private void MVNLMS(double A, double B, int INFIN, ref double LOWER, ref double UPPER)
        {
            LOWER = 0;
            UPPER = 1;
            if (INFIN >= 0)
            {
                if (INFIN != 0) { LOWER = MVNPHI(A); }
                if (INFIN != 1) { UPPER = MVNPHI(B); }
            }
            UPPER = Math.Max(UPPER, LOWER);
        }

        /// <summary>
        /// Subroutine to sort integration limits and determine Cholesky factor.
        /// </summary>
        /// <param name="N"></param>
        /// <param name="LOWER"></param>
        /// <param name="UPPER"></param>
        /// <param name="CORREL"></param>
        /// <param name="INFIN"></param>
        /// <param name="Y"></param>
        /// <param name="INFIS"></param>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="COV"></param>
        /// <param name="INFI"></param>
        private void COVSRT(int N, double[] LOWER, double[] UPPER, double[] CORREL, int[] INFIN, double[] Y, ref int INFIS, double[] A, double[] B, double[] COV, int[] INFI)
        {
            double SQTWPI = 2.506628274631001;
            double EPS = 1E-10;

            int IJ = 0;
            int II = 0;
            INFIS = 0;

            for (int i = 0; i < N; i++)
            {
                A[i] = 0;
                B[i] = 0;
                INFI[i] = INFIN[i];
                if (INFI[i] < 0)
                {
                    INFIS++;
                }
                else
                {
                    if (INFI[i] != 0) { A[i] = LOWER[i]; }
                    if (INFI[i] != 1) { B[i] = UPPER[i]; }
                }

                for (int j = 0; j < i; j++)
                {
                    COV[IJ] = CORREL[II];
                    II++;
                    IJ++;
                }
                COV[IJ] = 1;
                IJ++;
            }

            // First move any doubly infinite limits to innermost positions.
            if (INFIS < N)
            {
                for (int i = N - 1; i >= N - INFIS; i--)
                {
                    if (INFI[i] >= 0)
                    {
                        for (int j = 0; j < i - 1; j++)
                        {
                            if (INFI[j] < 0)
                            {
                                RCSWP(j, i, A, B, INFI, N, COV);
                                goto ten;
                            }
                        }
                    }
                }

            ten:

                // Sort remaining limits and determine Cholesky factor.
                II = 0;
                double AMIN = 0, BMIN = 0, DMIN, EMIN, CVDIAG, SUMSQ = 0, SUM, AJ = 0, BJ = 0, D = 0, E = 0, YL, YU;
                int JMIN, IL, M;
                for (int i = 0; i < N - INFIS; i++)
                {
                    // Determine the integration limits for variable with minimum expected probability and interchange that variable with Ith.

                    DMIN = 0;
                    EMIN = 1;
                    JMIN = i;
                    CVDIAG = 0;
                    IJ = II;

                    for (int j = i; j < N - INFIS; j++)
                    {


                        if (COV[IJ + j] > EPS)
                        {
                            SUMSQ = Math.Sqrt(COV[IJ + j]);
                            SUM = 0;
                            for (int k = 0; k < i; k++)
                            {
                                SUM += COV[IJ + k] * Y[k];
                            }

                            AJ = (A[j] - SUM) / SUMSQ;
                            BJ = (B[j] - SUM) / SUMSQ;

                            MVNLMS(AJ, BJ, INFI[j], ref D, ref E);

                            if (EMIN + D >= E + DMIN)
                            {
                                JMIN = j;
                                AMIN = AJ;
                                BMIN = BJ;
                                DMIN = D;
                                EMIN = E;
                                CVDIAG = SUMSQ;
                            }

                        }
                        IJ += (j + 1);
                    }

                    if (JMIN > i) { RCSWP(i, JMIN, A, B, INFI, N, COV); }
                    COV[II + i] = CVDIAG;


                    // Compute Ith column of Cholesky factor.
                    // Compute expected value for Ith integration variable and
                    // scale Ith covariance matrix row and limits.

                    if (CVDIAG > 0)
                    {
                        IL = II + (i + 1);
                        for (int l = i + 2; l <= N - INFIS; l++)
                        {
                            COV[IL + i] = COV[IL + i] / CVDIAG;
                            IJ = II + i;

                            for (int j = i + 1; j < l; j++)
                            {
                                COV[IL + j] = COV[IL + j] - COV[IL + i] * COV[IJ + i + 1];
                                IJ += j + 1;
                            }
                            IL += l;
                        }

                        if (EMIN > DMIN + EPS)
                        {
                            YL = 0;
                            YU = 0;
                            if (INFI[i] != 0) { YL = -Math.Exp(-Math.Pow(AMIN, 2) / 2) / SQTWPI; }
                            if (INFI[i] != 1) { YU = -Math.Exp(-Math.Pow(BMIN, 2) / 2) / SQTWPI; }
                            Y[i] = (YU - YL) / (EMIN - DMIN);
                        }
                        else
                        {
                            if (INFI[i] == 0) Y[i] = BMIN;
                            if (INFI[i] == 1) Y[i] = AMIN;
                            if (INFI[i] == 2) Y[i] = (AMIN + BMIN) / 2;
                        }
                        for (int j = 0; j <= i; j++)
                        {
                            COV[II] = COV[II] / CVDIAG;
                            II += 1;
                        }
                        A[i] = A[i] / CVDIAG;
                        B[i] = B[i] / CVDIAG;
                    }
                    else
                    {
                        IL = II + i;
                        for (int l = i + 1; l < N - INFIS; l++)
                        {
                            COV[IL + i] = 0;
                            IL += l;
                        }

                        // If the covariance matrix diagonal entry is zero, permute limits and/ or rows, if necessary.

                        for (int j = i - 1; j < 0; j--)
                        {
                            if (Math.Abs(COV[II + j]) > EPS)
                            {
                                A[i] = A[i] / COV[II + j];
                                B[i] = B[i] / COV[II + j];
                                if (COV[II + j] < 0)
                                {
                                    Tools.Swap(ref A[i], ref B[i]);
                                    if (INFI[i] != 2) INFI[i] = 1 - INFI[i];
                                }
                                for (int l = 0; l < j; l++)
                                {
                                    COV[II + l] = COV[II + l] / COV[II + j];
                                }
                                for (int l = j + 1; l < i - 1; l++)
                                {
                                    if (COV[(l - 1) * l / 2 + j + 1] > 0)
                                    {
                                        IJ = II;
                                        for (int k = i - 1; k < l; k--)
                                        {
                                            for (int m = 0; m < k; m++)
                                            {
                                                Tools.Swap(ref COV[IJ - k + m], ref COV[IJ + m]);
                                            }
                                            Tools.Swap(ref A[k], ref A[k + 1]);
                                            Tools.Swap(ref B[k], ref B[k + 1]);
                                            M = INFI[k];
                                            INFI[k] = INFI[k + 1];
                                            INFI[k + 1] = M;
                                            IJ = IJ - k;
                                        }
                                        goto twenty;
                                    }
                                }
                                goto twenty;
                            }
                            COV[II + j] = 0;
                        }
                    twenty:
                        II = II + i;
                        Y[i] = 0;
                    }

                }

            }

        }

        /// <summary>
        /// Swaps rows and columns P and Q in situ, with P &lt;= Q.
        /// </summary>
        /// <param name="P">Rows</param>
        /// <param name="Q">Columns</param>
        /// <param name="A">Array 1</param>
        /// <param name="B">Array 2</param>
        /// <param name="INFIN"></param>
        /// <param name="N"></param>
        /// <param name="C"></param>
        private void RCSWP(int P, int Q, double[] A, double[] B, int[] INFIN, int N, double[] C)
        {
            int II, JJ;
            Tools.Swap(ref A[P], ref A[Q]);
            Tools.Swap(ref B[P], ref B[Q]);
            int J = INFIN[P];
            INFIN[P] = INFIN[Q];
            INFIN[Q] = J;
            JJ = (P * (P + 1)) / 2;
            II = (Q * (Q + 1)) / 2;

            Tools.Swap(ref C[JJ + P], ref C[II + Q]);

            for (int i = 0; i < P; i++)
            {
                Tools.Swap(ref C[JJ + i], ref C[II + i]);
            }

            JJ += (P + 1);
            for (int i = P + 1; i < Q; i++)
            {
                Tools.Swap(ref C[JJ + P], ref C[II + i]);
                JJ += (i + 1);
            }

            II += (Q + 1);
            for (int i = Q + 1; i < N; i++)
            {
                Tools.Swap(ref C[II + P], ref C[II + Q]);
                II += (i + 1);
            }
        }

        /// <summary>
        /// Automatic Multidimensional Integration Subroutine
        ///          AUTHOR: Alan Genz
        ///                 Department of Mathematics
        /// Washington State University
        ///                 Pulman, WA 99164-3113
        ///                 Email: AlanGenz @wsu.edu
        ///
        ///         Last Change: 7/3/7
        ///
        ///  DKBVRC computes an approximation to the integral
        ///
        ///      1  1     1
        ///     I I ... I F(X) dx(NDIM)...dx(2)dx(1)
        ///      0  0     0
        ///
        ///  DKBVRC uses randomized Korobov rules for the first 100 variables.
        /// The primary references are
        ///   "Randomization of Number Theoretic Methods for Multiple Integration"
        ///    R.Cranley and T.N.L.Patterson, SIAM J Numer Anal, 13, pp. 904-14,
        ///  and
        ///   "Optimal Parameters for Multidimensional Integration", 
        ///    P.Keast, SIAM J Numer Anal, 10, pp.831-838.
        ///  If there are more than 100 variables, the remaining variables are
        ///  integrated using the rules described in the reference
        ///   "On a Number-Theoretical Integration Method"
        /// H.Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11.
        /// </summary>
        /// <param name="NDIM">Number of variables, must exceed 1, but not exceed 40.</param>
        /// <param name="MINVLS">minimum number of function evaluations allowed. MINVLS must not exceed MAXVLS.If MINVLS &lt; 0 then the routine assumes a previous call has been made with the same integrand and continues that calculation.</param>
        /// <param name="MAXVLS">maximum number of function evaluations allowed.</param>
        /// <param name="FUNCTN">EXTERNALly declared user defined function to be integrated. It must have parameters (NDIM, Z), where Z is a real array of dimension NDIM.</param>
        /// <param name="ABSEPS">Required absolute accuracy.</param>
        /// <param name="RELEPS">Required relative accuracy.</param>
        /// <param name="ABSERR">Estimated absolute accuracy of FINEST.</param>
        /// <param name="FINEST">Estimated value of integral.</param>
        /// <param name="INFORM">INFORM = 0 for normal exit, when ABSERR &lt;= MAX(ABSEPS, RELEPS* ABS(FINEST)) and INTVLS &lt;= MAXCLS. 
        ///                      INFORM = 1 If MAXVLS was too small to obtain the required accuracy.In this case a value FINEST is returned with estimated absolute accuracy ABSERR.</param>
        private void DKBVRC(int NDIM, int MINVLS, int MAXVLS, Func<int, double[], double> FUNCTN, double ABSEPS, double RELEPS, ref double ABSERR, ref double FINEST, ref int INFORM)
        {

            int PLIM = 28, NLIM = 1000, KLIM = 100, KLIMI, K, INTVLS, MINSMP = 8;

            //int[] P = new int[PLIM];
            //int[,] C = new int[PLIM, KLIM - 1];
            double DIFINT, FINVAL, VARSQR, VARPRD, VALUE = 0;
            double[] X = new double[2 * NLIM];
            double[] VK = new double[NLIM];

            // DKBVRC SAVE P, C, SAMPLS, NP, VAREST
            INFORM = 1;
            INTVLS = 0;
            KLIMI = KLIM;
            if (MINVLS >= 0)
            {
                FINEST = 0;
                DKBVRC_VAREST = 0;
                DKBVRC_SAMPLS = MINSMP;
                for (int i = Math.Min(NDIM, 10) - 1; i < PLIM; i++)
                {
                    DKBVRC_NP = i;
                    if (MINVLS < 2 * DKBVRC_SAMPLS * DKBVRC_P[i]) { goto ten; }
                }
                DKBVRC_SAMPLS = Math.Max(MINSMP, MINVLS / (2 * DKBVRC_P[DKBVRC_NP]));
            }

        ten:

            VK[0] = 1.0 / DKBVRC_P[DKBVRC_NP];
            K = 1;
            for (int i = 1; i < NDIM; i++)
            {
                if (i < KLIM)
                {
                    K = (int)((DKBVRC_C[DKBVRC_NP, Math.Min(NDIM - 1, KLIM - 1) - 1] * (double)K) % (double)DKBVRC_P[DKBVRC_NP]);
                    VK[i] = K * VK[0];
                }
                else
                {
                    VK[i] = (int)(Math.Pow(DKBVRC_P[DKBVRC_NP] * 2, ((double)(i - KLIM) / (NDIM - KLIM + 1))));
                    VK[i] = VK[i] / DKBVRC_P[DKBVRC_NP] % 1;
                }
            }

            FINVAL = 0;
            VARSQR = 0;

            for (int i = 1; i <= DKBVRC_SAMPLS; i++)
            {
                DKSMRC(NDIM, KLIMI, ref VALUE, DKBVRC_P[DKBVRC_NP], ref VK, FUNCTN, ref X);
                DIFINT = (VALUE - FINVAL) / i;
                FINVAL += DIFINT;
                VARSQR = (i - 2) * VARSQR / i + Math.Pow(DIFINT, 2);
            }

            INTVLS = INTVLS + 2 * DKBVRC_SAMPLS * DKBVRC_P[DKBVRC_NP];
            VARPRD = DKBVRC_VAREST * VARSQR;
            FINEST = FINEST + (FINVAL - FINEST) / (1 + VARPRD);
            if (VARSQR > 0) DKBVRC_VAREST = (1 + VARPRD) / VARSQR;
            ABSERR = 7 * Math.Sqrt(VARSQR / (1 + VARPRD)) / 2;

            if (ABSERR > Math.Max(ABSEPS, Math.Abs(FINEST) * RELEPS))
            {
                if (DKBVRC_NP < PLIM - 1)
                {
                    DKBVRC_NP = DKBVRC_NP + 1;
                }
                else
                {
                    DKBVRC_SAMPLS = Math.Min(3 * DKBVRC_SAMPLS / 2, (MAXVLS - INTVLS) / (2 * DKBVRC_P[DKBVRC_NP]));
                    DKBVRC_SAMPLS = Math.Max(MINSMP, DKBVRC_SAMPLS);
                }
                if (INTVLS + 2 * DKBVRC_SAMPLS * DKBVRC_P[DKBVRC_NP] <= MAXVLS) { goto ten; }
            }
            else
            {
                INFORM = 0;
            }

            MINVLS = INTVLS;
        }

        private void DKSMRC(int NDIM, int KLIM, ref double SUMKRO, int PRIME, ref double[] VK, Func<int, double[], double> FUNCTN, ref double[] X)
        {
            double sampleValue;
            double XT;
            int JP;
            SUMKRO = 0;
            int NK = Math.Min(NDIM, KLIM);
            for (int i = 0; i < NK - 1; i++)
            {
                sampleValue = MVNUNI.NextDouble();
                JP = (int)(i + sampleValue * (NK - i));
                XT = VK[i];
                VK[i] = VK[JP];
                VK[JP] = XT;
            }

            for (int i = 0; i < NDIM; i++)
            {
                X[NDIM + i] = MVNUNI.NextDouble();
            }


            for (int i = 1; i <= PRIME; i++)
            {

                for (int j = 0; j < NDIM; j++)
                {
                    X[j] = Math.Abs(2d * ((i * VK[j] + X[NDIM + j]) % 1d) - 1);
                }

                SUMKRO += (FUNCTN(NDIM, X) - SUMKRO) / (2 * i - 1);

                for (int j = 0; j < NDIM; j++)
                {
                    X[j] = 1 - X[j];
                }
                SUMKRO += (FUNCTN(NDIM, X) - SUMKRO) / (2 * i);
            }

        }

        /// <summary>
        /// A function for computing bivariate normal probabilities.
        /// 
        /// </summary>
        /// <param name="LOWER">array of lower integration limits.</param>
        /// <param name="UPPER">array of upper integration limits.</param>
        /// <param name="INFIN">array of integration limits flags: if INFIN(I) = 0, Ith limits are(-infinity, UPPER(I)]; if INFIN(I) = 1, Ith limits are[LOWER(I), infinity); if INFIN(I) = 2, Ith limits are[LOWER(I), UPPER(I)].</param>
        /// <param name="CORREL">correlation coefficient.</param>
        /// <returns></returns>
        private double BVNMVN(double[] LOWER, double[] UPPER, int[] INFIN, double CORREL)
        {

            double result = 0;
            if (INFIN[0] == 2 && INFIN[1] == 2)
            {
                result = BVU(LOWER[0], LOWER[1], CORREL) - BVU(UPPER[0], LOWER[1], CORREL) - BVU(LOWER[0], UPPER[1], CORREL) + BVU(UPPER[0], UPPER[1], CORREL);
            }
            else if (INFIN[0] == 2 && INFIN[1] == 1)
            {
                result = BVU(LOWER[0], LOWER[1], CORREL) - BVU(UPPER[0], LOWER[1], CORREL);
            }
            else if (INFIN[0] == 1 && INFIN[1] == 2)
            {
                result = BVU(LOWER[0], LOWER[1], CORREL) - BVU(LOWER[0], UPPER[1], CORREL);
            }
            else if (INFIN[0] == 2 && INFIN[1] == 0)
            {
                result = BVU(-UPPER[0], -UPPER[1], CORREL) - BVU(-LOWER[0], -UPPER[1], CORREL);
            }
            else if (INFIN[0] == 0 && INFIN[1] == 2)
            {
                result = BVU(-UPPER[0], -UPPER[1], CORREL) - BVU(-UPPER[0], -LOWER[1], CORREL);
            }
            else if (INFIN[0] == 1 && INFIN[1] == 0)
            {
                result = BVU(LOWER[0], -UPPER[1], -CORREL);
            }
            else if (INFIN[0] == 0 && INFIN[1] == 1)
            {
                result = BVU(-UPPER[0], LOWER[1], -CORREL);
            }
            else if (INFIN[0] == 1 && INFIN[1] == 1)
            {
                result = BVU(LOWER[0], LOWER[1], CORREL);
            }
            else if (INFIN[0] == 0 && INFIN[1] == 0)
            {
                result = BVU(-UPPER[0], -UPPER[1], CORREL);
            }

            // This can sometimes produce very small, but negative values (approx. -5E-73)
            if (result < 0d) result = 0d;
            return result;
        }

        public static double MVNPHI(double Z)
        {
            //     
            //     Normal distribution probabilities accurate to 1d-15.
            //     Reference: J.L.Schonfelder, Math Comp 32(1978), pp 1232-1240. 
            //
            double RTWO = 1.414213562373095048801688724209;
            int IM = 24;
            double BM, B, BP = 0, P, T, XA;

            //     
            XA = Math.Abs(Z) / RTWO;
            if (XA > 100)
            {
                P = 0;
            }
            else
            {
                T = (8 * XA - 30) / (4 * XA + 15);
                BM = 0;
                B = 0;
                for (int i = IM; i >= 0; i--)
                {
                    BP = B;
                    B = BM;
                    BM = T * B - BP + MVNPHI_A[i];
                }

                P = Math.Exp(-XA * XA) * (BM - BP) / 4;
            }

            if (Z > 0) P = 1 - P;

            return P;
        }

        //    Coefficients for P close to 0.5
        const double A0 = 3.3871328727963666080;
        const double A1 = 1.3314166789178437745E+2;
        const double A2 = 1.9715909503065514427E+3;
        const double A3 = 1.3731693765509461125E+4;
        const double A4 = 4.5921953931549871457E+4;
        const double A5 = 6.7265770927008700853E+4;
        const double A6 = 3.3430575583588128105E+4;
        const double A7 = 2.5090809287301226727E+3;
        const double B1 = 4.2313330701600911252E+1;
        const double B2 = 6.8718700749205790830E+2;
        const double B3 = 5.3941960214247511077E+3;
        const double B4 = 2.1213794301586595867E+4;
        const double B5 = 3.9307895800092710610E+4;
        const double B6 = 2.8729085735721942674E+4;
        const double B7 = 5.2264952788528545610E+3;
        // HASH SUM AB    55.88319 28806 14901 4439

        //     Coefficients for P not close to 0, 0.5 or 1.    
        const double C0 = 1.42343711074968357734;
        const double C1 = 4.63033784615654529590;
        const double C2 = 5.76949722146069140550;
        const double C3 = 3.64784832476320460504;
        const double C4 = 1.27045825245236838258;
        const double C5 = 2.41780725177450611770E-1;
        const double C6 = 2.27238449892691845833E-2;
        const double C7 = 7.74545014278341407640E-4;
        const double D1 = 2.05319162663775882187;
        const double D2 = 1.67638483018380384940;
        const double D3 = 6.89767334985100004550E-1;
        const double D4 = 1.48103976427480074590E-1;
        const double D5 = 1.51986665636164571966E-2;
        const double D6 = 5.47593808499534494600E-4;
        const double D7 = 1.05075007164441684324E-9;
        //     HASH SUM CD    49.33206 50330 16102 89036

        //	Coefficients for P near 0 or 1.
        const double E0 = 6.65790464350110377720;
        const double E1 = 5.46378491116411436990;
        const double E2 = 1.78482653991729133580;
        const double E3 = 2.96560571828504891230E-1;
        const double E4 = 2.65321895265761230930E-2;
        const double E5 = 1.24266094738807843860E-3;
        const double E6 = 2.71155556874348757815E-5;
        const double E7 = 2.01033439929228813265E-7;
        const double F1 = 5.99832206555887937690E-1;
        const double F2 = 1.36929880922735805310E-1;
        const double F3 = 1.48753612908506148525E-2;
        const double F4 = 7.86869131145613259100E-4;
        const double F5 = 1.84631831751005468180E-5;
        const double F6 = 1.42151175831644588870E-7;
        const double F7 = 2.04426310338993978564E-15;
        //     HASH SUM EF    47.52583 31754 92896 71629

        /// <summary>
        /// ALGORITHM AS241  APPL.STATIST. (1988) VOL. 37, NO. 3
        /// 
        /// Produces the normal deviate Z corresponding to a given lower tail area of p.
        /// </summary>
        /// <param name="P">Lower tail area.</param>
        /// <returns></returns>
        public static double PHINVS(double P)
        {
            double SPLIT1 = 0.425, SPLIT2 = 5, CONST1 = 0.180625, CONST2 = 1.6;

            double Q = (2 * P - 1) / 2;
            double R;
            if (Math.Abs(Q) < SPLIT1)
            {
                R = CONST1 - Q * Q;
                return Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) * R + A2) * R + A1) * R + A0) / (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3) * R + B2) * R + B1) * R + 1);
            }
            else
            {
                double result;
                R = Math.Min(P, 1 - P);
                if (R > 0)
                {
                    R = Math.Sqrt(-Math.Log(R));
                    if (R < SPLIT2)
                    {
                        R -= CONST2;
                        result = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3) * R + C2) * R + C1) * R + C0) / (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3) * R + D2) * R + D1) * R + 1);
                    }
                    else
                    {
                        R -= SPLIT2;
                        result = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3) * R + E2) * R + E1) * R + E0) / (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3) * R + F2) * R + F1) * R + 1);
                    }
                }
                else
                {
                    result = 9;
                }
                if (Q < 0) { result = -result; }
                return result;
            }
        }

        /// <summary>
        /// Adapted From:
        /// A function for computing bivariate normal probabilities.
        ///
        ///       Yihong Ge
        ///       Department of Computer Science and Electrical Engineering
        /// Washington State University
        ///       Pullman, WA 99164-2752
        ///     and
        /// Alan Genz
        /// Department of Mathematics
        ///       Washington State University
        /// Pullman, WA 99164-3113
        ///       Email : alangenz @wsu.edu
        ///
        /// BVN - calculate the probability that X is larger than SH and Y is
        ///       larger than SK.
        /// </summary>
        /// <param name="SH">integration limit</param>
        /// <param name="SK">integration limit</param>
        /// <param name="R">correlation coefficient</param>
        /// <returns></returns>
        public static double BVU(double SH, double SK, double R)
        {
            //
            double BVN;
            int LG, NG;
            double TWOPI = 6.283185307179586;
            double AS, A, B, C, D, RS, XS;
            double SN, ASR, H, K, BS, HS, HK;

            if (Math.Abs(R) < .3)
            {
                NG = 0;
                LG = 3;
            }
            else if (Math.Abs(R) < .75)
            {
                NG = 1;
                LG = 6;
            }
            else
            {
                NG = 2;
                LG = 10;
            }

            H = SH;
            K = SK;
            HK = H * K;
            BVN = 0;
            if (Math.Abs(R) < .925)
            {
                HS = (H * H + K * K) / 2;
                ASR = Math.Asin(R);
                for (int i = 0; i < LG; i++)
                {
                    SN = Math.Sin(ASR * (X[i, NG] + 1) / 2);
                    BVN += W[i, NG] * Math.Exp((SN * HK - HS) / (1 - SN * SN));
                    SN = Math.Sin(ASR * (-X[i, NG] + 1) / 2);
                    BVN += W[i, NG] * Math.Exp((SN * HK - HS) / (1 - SN * SN));
                }

                BVN = BVN * ASR / (2 * TWOPI) + MVNPHI(-H) * MVNPHI(-K);
            }
            else
            {
                if (R < 0)
                {
                    K = -K;
                    HK = -HK;
                }
                if (Math.Abs(R) < 1)
                {
                    AS = (1 - R) * (1 + R);
                    A = Math.Sqrt(AS);
                    BS = Math.Pow(H - K, 2);
                    C = (4 - HK) / 8;
                    D = (12 - HK) / 16;
                    BVN = A * Math.Exp(-(BS / AS + HK) / 2) * (1 - C * (BS - AS) * (1 - D * BS / 5) / 3 + C * D * AS * AS / 5);
                    if (HK > -160)
                    {
                        B = Math.Sqrt(BS);
                        BVN = BVN - Math.Exp(-HK / 2) * Math.Sqrt(TWOPI) * MVNPHI(-B / A) * B * (1 - C * BS * (1 - D * BS / 5) / 3);
                    }
                    A = A / 2;
                    for (int i = 0; i < LG; i++)
                    {
                        XS = Math.Pow(A * (X[i, NG] + 1), 2);
                        RS = Math.Sqrt(1 - XS);
                        BVN = BVN + A * W[i, NG] * (Math.Exp(-BS / (2 * XS) - HK / (1 + RS)) / RS - Math.Exp(-(BS / XS + HK) / 2) * (1 + C * XS * (1 + D * XS)));
                        XS = AS * Math.Pow(-X[i, NG] + 1, 2) / 4;
                        RS = Math.Sqrt(1 - XS);
                        BVN = BVN + A * W[i, NG] * Math.Exp(-(BS / XS + HK) / 2) * (Math.Exp(-HK * XS / (2 * Math.Pow(1 + RS, 2))) / RS - (1 + C * XS * (1 + D * XS)));
                    }

                    BVN = -BVN / TWOPI;
                }

                if (R > 0)
                {
                    BVN = BVN + MVNPHI(-Math.Max(H, K));
                }
                else
                {
                    BVN = -BVN;
                    if (K > H)
                    {
                        if (H < 0)
                        {
                            BVN = BVN + MVNPHI(K) - MVNPHI(H);
                        }
                        else
                        {
                            BVN = BVN + MVNPHI(-H) - MVNPHI(-K);
                        }
                    }
                }
            }

            // This can sometimes produce very small, but negative values (approx. -5E-73)
            if (BVN < 0d) BVN = 0d;
            return BVN;
        }


        #endregion


        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public override MultivariateDistribution Clone()
        {
            var mvn = new MultivariateNormal()
            {
                _parametersValid = this.ParametersValid,
                _dimension = this.Dimension,
                _mean = this.Mean.ToArray(),
                _covariance = this._covariance.Clone(),
                _cholesky = this._cholesky,
                _lnconstant = this._lnconstant,
                _variance = this.Variance.ToArray(),
                _standardDeviation = this.StandardDeviation.ToArray(),
                _correlation = this._correlation.Clone(),
                _correl = this._correl.ToArray(),
                _MVNUNI = this._MVNUNI,
                _maxEvaluations = this._maxEvaluations,
                _absoluteError = this._absoluteError,
                _relativeError = this._relativeError,
                _lower = this._lower.ToArray(),
                _upper = this._upper.ToArray(),
                _infin = this._infin.ToArray(),
                _covSRTed = this._covSRTed,
                _correlationMatrixCreated = this._correlationMatrixCreated,
            };

            return mvn;

        }

    }
}