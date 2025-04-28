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
using System.Threading.Tasks;
using Numerics.Data.Statistics;
using Numerics.Mathematics.Optimization;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;

namespace Numerics.Distributions
{

    /// <summary>
    /// The Normal (Gaussian) probability distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Wikipedia contributors, "Normal distribution,". Wikipedia, The Free
    /// Encyclopedia. Available at: <see href="https://en.wikipedia.org/wiki/Normal_distribution"/>
    /// </description></item>
    /// <item><description>
    /// Accord Math Library, <see href="http://accord-framework.net"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    [Serializable]
    public sealed class Normal : UnivariateDistributionBase, IEstimation, IMaximumLikelihoodEstimation, IMomentEstimation, ILinearMomentEstimation, IStandardError, IBootstrappable
    {
      
        /// <summary>
        /// Constructs a Normal (Gaussian) distribution with a mean of 0 and standard deviation of 1.
        /// </summary>
        public Normal()
        {
            SetParameters(0d, 1d);
        }

        /// <summary>
        /// Constructs a Normal (Gaussian) distribution with a given mean and a standard deviation of 1.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        public Normal(double mean)
        {
            SetParameters(mean, 1d);
        }

        /// <summary>
        /// Constructs a Normal (Gaussian) distribution with given mean and standard deviation.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="standardDeviation">Standard deviation of the distribution.</param>
        public Normal(double mean, double standardDeviation)
        {
            SetParameters(mean, standardDeviation);
        }

        private double _mu;
        private double _sigma;

        /// <summary>
        /// Gets and sets the location parameter µ (Mu).
        /// </summary>
        public double Mu
        {
            get { return _mu; }
            set
            {
                _parametersValid = ValidateParameters(value, Sigma, false) is null;
                _mu = value;
            }
        }

        /// <summary>
        /// Gets and sets the scale parameter σ (sigma).
        /// </summary>
        public double Sigma
        {
            get { return _sigma; }
            set
            {
                if (value < 1E-16 && Math.Sign(value) != -1) value = 1E-16;
                _parametersValid = ValidateParameters(Mu, value, false) is null;
                _sigma = value;
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
            get { return UnivariateDistributionType.Normal; }
        }

        /// <inheritdoc/>
        public override string DisplayName
        {
            get { return "Normal"; }
        }

        /// <inheritdoc/>
        public override string ShortDisplayName
        {
            get { return "N"; }
        }

        /// <inheritdoc/>
        public override string[,] ParametersToString
        {
            get
            {
                var parmString = new string[2, 2];
                parmString[0, 0] = "Mean (µ)";
                parmString[1, 0] = "Std Dev (σ)";
                parmString[0, 1] = Mu.ToString();
                parmString[1, 1] = Sigma.ToString();
                return parmString;
            }
        }

        /// <inheritdoc/>
        public override string[] ParameterNamesShortForm
        {
            get { return ["µ", "σ"]; }
        }

        /// <inheritdoc/>
        public override string[] GetParameterPropertyNames
        {
            get { return [nameof(Mu), nameof(Sigma)]; }
        }

        /// <inheritdoc/>
        public override double[] GetParameters
        {
            get { return [Mu, Sigma]; }
        }

        /// <inheritdoc/>
        public override double Mean
        {
            get
            {
                return Mu;
            }
        }

        /// <inheritdoc/>
        public override double Median
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
        public override double Mode
        {
            get { return Mu; }
        }

        /// <inheritdoc/>
        public override double StandardDeviation
        {
            get { return Sigma; }
        }

        /// <inheritdoc/>
        public override double Skewness
        {
            get { return 0.0d; }
        }

        /// <inheritdoc/>
        public override double Kurtosis
        {
            get { return 3d; }
        }

        /// <inheritdoc/>
        public override double Minimum
        {
            get { return double.NegativeInfinity; }
        }

        /// <inheritdoc/>
        public override double Maximum
        {
            get { return double.PositiveInfinity; }
        }

        /// <inheritdoc/>
        public override double[] MinimumOfParameters
        {
            get { return [double.NegativeInfinity, 0.0d]; }
        }

        /// <inheritdoc/>
        public override double[] MaximumOfParameters
        {
            get { return [double.PositiveInfinity, double.PositiveInfinity]; }
        }

        /// <inheritdoc/>
        public void Estimate(IList<double> sample, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod == ParameterEstimationMethod.MethodOfMoments)
            {
                SetParameters(Statistics.ProductMoments(sample));
            }
            else if (estimationMethod == ParameterEstimationMethod.MethodOfLinearMoments)
            {
                SetParameters(ParametersFromLinearMoments(Statistics.LinearMoments(sample)));
            }
            else if (estimationMethod == ParameterEstimationMethod.MaximumLikelihood)
            {
                // The only difference between this method and MOM is that the population standard deviation is used
                // rather than the sample; Note: the maximum likelihood estimate of the standard deviation
                // is a biased estimator.
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
            var newDistribution = new Normal(Mu, Sigma);
            var sample = newDistribution.GenerateRandomValues(sampleSize, seed);
            newDistribution.Estimate(sample, estimationMethod);
            if (newDistribution.ParametersValid == false)
                throw new Exception("Bootstrapped distribution parameters are invalid.");
            return newDistribution;
        }

        /// <summary>
        /// Set the distribution parameters.
        /// </summary>
        /// <param name="location">The location parameter µ (Mu).</param>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        public void SetParameters(double location, double scale)
        {
            Mu = location;
            Sigma = scale;
        }

        /// <inheritdoc/>
        public override void SetParameters(IList<double> parameters)
        {
            SetParameters(parameters[0], parameters[1]);
        }

        /// <inheritdoc/>
        public double[] ParametersFromMoments(IList<double> moments)
        {
            return moments.ToArray().Subset(0, 1);
        }

        /// <inheritdoc/>
        public double[] MomentsFromParameters(IList<double> parameters)
        {
            var dist = new Normal();
            dist.SetParameters(parameters);
            var m1 = dist.Mean;
            var m2 = dist.StandardDeviation;
            var m3 = dist.Skewness;
            var m4 = dist.Kurtosis;
            return [m1, m2, m3, m4];
        }

        /// <inheritdoc/>
        public double[] ParametersFromLinearMoments(IList<double> moments)
        {
            double mu = moments[0];
            double sigma = moments[1] * Math.Sqrt(Math.PI);
            return [mu, sigma];
        }

        /// <inheritdoc/>
        public double[] LinearMomentsFromParameters(IList<double> parameters)
        {
            double L1 = parameters[0];
            double L2 = parameters[1] * Math.Pow(Math.PI, -0.5);
            double T3 = 0d;
            double T4 = 30d * Math.Pow(Math.PI, -1d) * Math.Atan(Tools.Sqrt2) - 9d;
            return [L1, L2, T3, T4];
        }

        /// <summary>
        /// Validate the parameters.
        /// </summary>
        /// <param name="location">The location parameter µ (Mu).</param>
        /// <param name="scale">The scale parameter σ (sigma).</param>
        /// <param name="throwException">Determines whether to throw an exception or not.</param>
        public ArgumentOutOfRangeException ValidateParameters(double location, double scale, bool throwException)
        {
            if (double.IsNaN(location) || double.IsInfinity(location))
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Mu), "The mean must be a number.");
                return new ArgumentOutOfRangeException(nameof(Mu), "The mean must be a number.");
            }
            if (double.IsNaN(scale) || double.IsInfinity(scale) || scale <= 0.0d)
            {
                if (throwException)
                    throw new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
                return new ArgumentOutOfRangeException(nameof(Sigma), "Standard deviation must be positive.");
            }
            return null;
        }

        /// <inheritdoc/>
        public override ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException)
        {
            return ValidateParameters(parameters[0], parameters[1], throwException);
        }

        /// <inheritdoc/>
        public Tuple<double[], double[], double[]> GetParameterConstraints(IList<double> sample)
        {
            var initialVals = new double[NumberOfParameters];
            var lowerVals = new double[NumberOfParameters];
            var upperVals = new double[NumberOfParameters];
            // Estimate initial values using the method of moments (a.k.a product moments).
            var moments = Statistics.ProductMoments(sample);
            initialVals[0] = moments[0];
            initialVals[1] = moments[1];
            // Get bounds of mean
            if (initialVals[0] == 0d) initialVals[0] = Tools.DoubleMachineEpsilon;
            lowerVals[0] = -Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            upperVals[0] = Math.Pow(10d, Math.Ceiling(Math.Log10(Math.Abs(initialVals[0])) + 1d));
            // Get bounds of standard deviation
            lowerVals[1] = Tools.DoubleMachineEpsilon;
            upperVals[1] = Math.Pow(10d, Math.Ceiling(Math.Log10(initialVals[1]) + 1d));
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

            // Solve using Nelder-Mead (Downhill Simplex)
            double logLH(double[] x)
            {
                var N = new Normal();
                N.SetParameters(x);
                return N.LogLikelihood(sample);
            }
            var solver = new NelderMead(logLH, NumberOfParameters, Initials, Lowers, Uppers);
            solver.ReportFailure = true;
            solver.Maximize();
            return solver.BestParameterSet.Values;
        }

        /// <inheritdoc/>
        public override double PDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            double z = (x - Mu) / Sigma;
            return Math.Exp(-0.5d * z * z) / (Tools.Sqrt2PI * Sigma);
        }

        /// <inheritdoc/>
        public override double CDF(double x)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            return 0.5d * (1.0d + Erf.Function((x - Mu) / (Sigma * Tools.Sqrt2)));
        }

        private static readonly double[] a = [3.3871328727963666080, 1.3314166789178437745e+2, 1.9715909503065514427e+3, 1.3731693765509461125e+4, 4.5921953931549871457e+4, 6.7265770927008700853e+4, 3.3430575583588128105e+4, 2.5090809287301226727e+3];
        private static readonly double[] b = [1.0, 4.2313330701600911252e+1, 6.8718700749205790830e+2, 5.3941960214247511077e+3, 2.1213794301586595867e+4, 3.9307895800092710610e+4, 2.8729085735721942674e+4, 5.2264952788528545610e+3];
        private static readonly double[] c = [1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550, 3.64784832476320460504, 1.27045825245236838258, 2.41780725177450611770e-1, 2.27238449892691845833e-2, 7.74545014278341407640e-4];
        private static readonly double[] d = [1.0, 2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1, 1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4, 1.05075007164441684324e-9];
        private static readonly double[] e = [6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580, 2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3, 2.71155556874348757815e-5, 2.01033439929228813265e-7];
        private static readonly double[] f = [1.0, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2, 7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7, 2.04426310338993978564e-15];

        /// <summary>
        /// R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
        /// </summary>
        /// <param name="p">The probability value.</param>
        private static double r8_normal_01_cdf_inverse(double p)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_NORMAL_01_CDF_INVERSE inverts the standard normal CDF.
            //
            //  Discussion:
            //
            //    The result is accurate to about 1 part in 10**16.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 March 2010
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Michael Wichura.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Michael Wichura,
            //    The Percentage Points of the Normal Distribution,
            //    Algorithm AS 241,
            //    Applied Statistics,
            //    Volume 37, Number 3, pages 477-484, 1988.
            //
            //  Parameters:
            //
            //    Input, double P, the value of the cumulative probability 
            //    density function.  0 < P < 1.  If P is outside this range, an "infinite"
            //    value is returned.
            //
            //    Output, double R8_NORMAL_01_CDF_INVERSE, the normal deviate value 
            //    with the property that the probability of a standard normal deviate being 
            //    less than or equal to this value is P.
            //

            double q;
            double r;
            double value;

            if (p <= 0.0)
            {
                return double.NegativeInfinity;
            }
            if (1.0 <= p)
            {
                return double.PositiveInfinity;
            }

            q = p - 0.5;

            // 0.075 <= p <= 0.925
            if (Math.Abs(q) <= 0.425)
            {
                r = 0.180625 - q * q;
                value = q * r8poly_value(8, a, r) / r8poly_value(8, b, r);
            }
            else
            {
                if (q < 0.0)
                {
                    r = p;
                }
                else
                {
                    r = 1.0 - p;
                }

                r = Math.Sqrt(-Math.Log(r));

                if (r <= 5.0)
                {
                    r = r - 1.6;
                    value = r8poly_value(8, c, r) / r8poly_value(8, d, r);
                }
                else
                {
                    r = r - 5.0;
                    value = r8poly_value(8, e, r) / r8poly_value(8, f, r);
                }

                if (q < 0.0)
                {
                    value = -value;
                }

            }

            return value;
        }

        /// <summary>
        /// R8POLY_VALUE evaluates a double precision polynomial.
        /// </summary>
        /// <param name="n">The number of coefficients.</param>
        /// <param name="a">The coefficients.</param>
        /// <param name="x">The point to evaluate.</param>
        private static double r8poly_value(int n, double[] a, double x)
        {
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VALUE evaluates a double precision polynomial.
            //
            //  Discussion:
            //
            //    For sanity's sake, the value of N indicates the NUMBER of 
            //    coefficients, or more precisely, the ORDER of the polynomial,
            //    rather than the DEGREE of the polynomial.  The two quantities
            //    differ by 1, but cause a great deal of confusion.
            //
            //    Given N and A, the form of the polynomial is:
            //
            //      p(x) = a[0] + a[1] * x + ... + a[n-2] * x^(n-2) + a[n-1] * x^(n-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double A[N], the coefficients of the polynomial.
            //    A[0] is the constant term.
            //
            //    Input, double X, the point at which the polynomial is to be evaluated.
            //
            //    Output, double R8POLY_VALUE, the value of the polynomial at X.
            //

            int i;
            double value = 0.0;
            for (i = n - 1; 0 <= i; i--)
            {
                value = value * x + a[i];
            }

            return value;
        }

        /// <inheritdoc/>
        public override double InverseCDF(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            if (probability == 0.0d)
                return Minimum;
            if (probability == 1.0d)
                return Maximum;
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, Sigma, true);
            return Mu + Sigma * r8_normal_01_cdf_inverse(probability);
        }


        /// <summary>
        /// Returns the PDF for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="Z">The Z variate of a standard Normal.</param>
        public static double StandardPDF(double Z)
        {
            return Math.Exp(-0.5d * Z * Z) / Tools.Sqrt2PI;
        }

        /// <summary>
        /// Returns the log PDF for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="Z">The Z variate of a standard Normal.</param>
        public static double StandardLogPDF(double Z)
        {
            return -0.5d * Z * Z - Tools.LogSqrt2PI;
        }

        /// <summary>
        /// Returns the PDF for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="zValues">A list of Z variates.</param>
        public static double[] StandardPDF(IList<double> zValues)
        {
            var result = new double[zValues.Count];
            for (int i = 0; i < zValues.Count; i++)
                result[i] = StandardPDF(zValues[i]);
            return result;
        }

        /// <summary>
        /// Returns the CDF for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="Z">The Z variate of a standard Normal.</param>
        public static double StandardCDF(double Z)
        {
            return MultivariateNormal.MVNPHI(Z);
        }

        /// <summary>
        /// Returns the CDF for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="zValues">A list of Z variates.</param>
        /// <returns>An array of probabilities.</returns>
        public static double[] StandardCDF(IList<double> zValues)
        {
            var result = new double[zValues.Count];
            for (int i = 0; i < zValues.Count; i++)
                result[i] = StandardCDF(zValues[i]);
            return result;
        }

        /// <summary>
        /// Returns the Z variate for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="probability">Probability between 0 and 1.</param>
        public static double StandardZ(double probability)
        {
            // Validate probability
            if (probability < 0.0d || probability > 1.0d)
                throw new ArgumentOutOfRangeException("probability", "Probability must be between 0 and 1.");
            double limit = 1E-16;
            //if (probability <= limit)
            //    return -8.2220822161304348;
            if (probability >= 1d - limit)
                return 8.2095361516013856;
            return r8_normal_01_cdf_inverse(probability);
        }

        /// <summary>
        /// Returns the Z variate for a standard Normal distribution, where mean of 0 and standard deviation of 1.
        /// </summary>
        /// <param name="probabilities">A list of probabilities.</param>
        /// <returns>An array of Z variates.</returns>
        public static double[] StandardZ(IList<double> probabilities)
        {
            var result = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                result[i] = StandardZ(probabilities[i]);
            return result;
        }

        /// <summary>
        /// Get confidence intervals based on the Normal approximation.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="quantiles">List of nonexceedance probabilities for output frequency curves.</param>
        /// <param name="percentiles">List of confidence percentiles for confidence interval output.</param>
        /// <remarks>
        /// References: Stedinger, J. Confidence Intervals for Design Events. Journal of Hydraulic Engineering. 1983.
        /// </remarks>
        public double[,] NormalConfidenceIntervals(int sampleSize, IList<double> quantiles, IList<double> percentiles)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Dimension output array
            int q = quantiles.Count;
            int p = percentiles.Count;
            var Output = new double[q, p];
            for (int i = 0; i < q; i++)
            {
                double Xp = InverseCDF(quantiles[i]);
                double Zp = StandardZ(quantiles[i]);
                // Record percentiles for user-defined probabilities
                for (int j = 0; j < p; j++)
                {
                    double Z = StandardZ(percentiles[j]);
                    double VARp = Sigma * Sigma / sampleSize * (1.0d + 0.5d * Zp * Zp);
                    Output[i, j] = Xp + Z * Math.Sqrt(VARp);
                }
            }
            // Return confidence percentile output
            return Output;
        }

        /// <summary>
        /// Get exact confidence intervals based on the Noncentral-T distribution.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="quantiles">List of nonexceedance probabilities for output frequency curves.</param>
        /// <param name="percentiles">List of confidence percentiles for confidence interval output.</param>
        /// <remarks>
        /// References: Stedinger, J. Confidence Intervals for Design Events. Journal of Hydraulic Engineering. 1983.
        /// </remarks>
        public double[,] NoncentralTConfidenceIntervals(int sampleSize, IList<double> quantiles, IList<double> percentiles)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Dimension output array
            int q = quantiles.Count;
            int p = percentiles.Count;
            var Output = new double[q, p];
            int df = sampleSize - 1;
            double SE = Sigma / Math.Sqrt(sampleSize);
            for (int i = 0; i < q; i++)
            {
                double Z = StandardZ(quantiles[i]);
                var NCT = new NoncentralT(df, Z * Math.Sqrt(sampleSize));
                // Record percentiles for user-defined probabilities
                for (int j = 0; j < p; j++)
                    Output[i, j] = Mu + NCT.InverseCDF(percentiles[j]) * SE;
            }
            // Return confidence percentile output
            return Output;
        }

        /// <summary>
        /// Get confidence intervals using Monte Carlo simulation.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="realizations">The number of Monte Carlo realizations.</param>
        /// <param name="quantiles">List of nonexceedance probabilities for output frequency curves.</param>
        /// <param name="percentiles">List of confidence percentiles for confidence interval output.</param>
        /// <remarks>
        /// This is the same sampling approach as used in HEC-FDA.
        /// </remarks>
        public double[,] MonteCarloConfidenceIntervals(int sampleSize, int realizations, IList<double> quantiles, IList<double> percentiles)
        {
            // validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            // Dimension output array
            int q = quantiles.Count;
            int p = percentiles.Count;
            var Output = new double[q, p];

            // Variables
            double OriginalMean = Mu;
            double OriginalStdDev = Sigma;

            // Create random numbers for mean and standard deviation
            var r = new MersenneTwister(12345);
            var rndMean = r.NextDoubles(realizations);
            r = new MersenneTwister(45678);
            var rndStdDev = r.NextDoubles(realizations);

            // Create list of Monte Carlo distributions
            var MonteCarloDistributions = new UnivariateDistributionBase[realizations];
            Parallel.For(0, realizations, idx =>
            {

                // Generate new mean
                var Normal = new Normal(OriginalMean, OriginalStdDev / Math.Sqrt(sampleSize));
                double NewMu = Normal.InverseCDF(rndMean[idx]);
                // Generate new standard deviation
                var Chi = new ChiSquared(sampleSize - 1);
                double NewSigma = Math.Sqrt((sampleSize - 1) * Math.Pow(OriginalStdDev, 2d) / Chi.InverseCDF(rndStdDev[idx]));
                // Create a new distribution with the new parameters
                MonteCarloDistributions[idx] = new Normal(NewMu, NewSigma);
            });
            // Create confidence intervals
            for (int i = 0; i < q; i++)
            {
                // Create array of X values across user-defined probabilities
                var XValues = new double[realizations];
                // Record X values
                Parallel.For(0, realizations, idx => XValues[idx] = MonteCarloDistributions[idx].InverseCDF(quantiles[i]));
                // Record percentiles for user-defined probabilities
                for (int j = 0; j < p; j++)
                    Output[i, j] = Statistics.Percentile(XValues, percentiles[j]);
            }
            // Return confidence percentile output
            return Output;
        }

        /// <summary>
        /// Gets the expected probability given an exceedance probability.
        /// </summary>
        /// <param name="sampleSize">The data sample size N used for computing the standard error.</param>
        /// <param name="probability">Exceedance probability.</param>
        public double ExpectedProbability(int sampleSize, double probability)
        {
            int N = sampleSize;
            var T = new StudentT(N - 1);
            return T.CDF(StandardZ(probability) * Math.Sqrt(N / (double)(N + 1)));
        }

        /// <inheritdoc/>
        public override UnivariateDistributionBase Clone()
        {
            return new Normal(Mu, Sigma);
        }

        /// <inheritdoc/>
        public double[,] ParameterCovariance(int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            if (estimationMethod != ParameterEstimationMethod.MethodOfMoments &&
                estimationMethod != ParameterEstimationMethod.MaximumLikelihood)
            {
                throw new NotImplementedException();
            }
            // Compute covariance
            double u2 = Sigma;
            var covar = new double[2, 2];
            covar[0, 0] = Math.Pow(u2, 2d) / sampleSize; // location
            covar[1, 1] = 2d * Math.Pow(u2, 4d) / sampleSize; // scale
            covar[0, 1] = 0.0;
            covar[1, 0] = covar[0, 1];
            return covar;
        }

        /// <inheritdoc/>
        public double QuantileVariance(double probability, int sampleSize, ParameterEstimationMethod estimationMethod)
        {
            var covar = ParameterCovariance(sampleSize, estimationMethod);
            var grad = QuantileGradient(probability);
            double varA = covar[0, 0];
            double varB = covar[1, 1];
            double covAB = covar[1, 0];
            double dQx1 = grad[0];
            double dQx2 = grad[1];
            return Math.Pow(dQx1, 2d) * varA + Math.Pow(dQx2, 2d) * varB + 2d * dQx1 * dQx2 * covAB;
        }

        /// <inheritdoc/>
        public double[] QuantileGradient(double probability)
        {
            // Validate parameters
            if (_parametersValid == false)
                ValidateParameters(Mu, _sigma, true);
            double u2 = Sigma;
            double z = StandardZ(probability);
            var gradient = new double[]
            {
                1.0d, // location
                z / (2d * u2) // scale
            };
            return gradient;
        }

        /// <inheritdoc/>
        public double[,] QuantileJacobian(IList<double> probabilities, out double determinant)
        {
            if (probabilities.Count != NumberOfParameters)
            {
                throw new ArgumentOutOfRangeException(nameof(probabilities), "The number of probabilities must be the same length as the number of distribution parameters.");
            }
            // Get gradients
            var dQp1 = QuantileGradient(probabilities[0]);
            var dQp2 = QuantileGradient(probabilities[1]);
            // Compute determinant
            // |a b|
            // |c d|
            // |A| = ad − bc
            double a = dQp1[0];
            double b = dQp1[1];
            double c = dQp2[0];
            double d = dQp2[1];
            determinant = a * d - b * c;
            // Return Jacobian
            var jacobian = new double[,] { { a, b }, { c, d } };
            return jacobian;
        }

    }
}