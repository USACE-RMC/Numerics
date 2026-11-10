using System;
using System.Collections.Generic;
using Numerics.Distributions;

namespace Numerics.Data.Statistics
{

    /// <summary>
    /// Computes the autocovariance, autocorrelation, or partial autocorrelation function.
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    ///<b> Description: </b>
    /// Autocovariance is a function that gives the covariance (measure of the joint variability of two random variables)
    /// of the process with itself at pairs of time points.
    /// Autocorrelation is the similarity of between observations of a random variable as a function of the time lag between them.
    /// The partial autocorrelation function gives the partial correlation of a stationary time series with it own lagged values of the 
    /// time series at all short lags. This is in contrast with the autocorrelation function that does not control for other lags.
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description> 
    /// <see href="https://en.wikipedia.org/wiki/Autocovariance"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Autocorrelation"/>
    /// </description></item>
    /// <item><description>
    /// 
    /// </description></item>
    /// </list>
    /// </remarks>
    public class Autocorrelation
    {
        /// <summary>
        /// Enumeration of the type of autocorrelation.
        /// </summary>
        public enum Type
        {
            /// <summary>
            /// Autocorrelation
            /// </summary>
            Correlation,
            /// <summary>
            /// Autocovariance
            /// </summary>
            Covariance,
            /// <summary>
            /// Partial autocorrelation
            /// </summary>
            Partial
        }

        /// <summary>
        /// Computes the autocovariance, autocorrelation, or partial autocorrelation function.
        /// </summary>
        /// <param name="data">The list of data to evaluate.</param>
        /// <param name="lagMax">The maximum lag at which to estimate the function. Default is 10*log10(N/m) where N is the number of observations.
        /// Lags at or beyond the number of observations give zero.</param>
        /// <param name="type">The type of function to be computed. Default is correlation. </param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the function evaluated at the given values.
        /// </returns>
        public static double[,]? Function(IList<double> data, int lagMax = -1, Type type = Type.Correlation)
        {
            return Function(data, lagMax, type, true);
        }

        /// <summary>
        /// Computes the autocovariance, autocorrelation, or partial autocorrelation function, about the sample mean or
        /// about zero.
        /// </summary>
        /// <param name="data">The list of data to evaluate.</param>
        /// <param name="lagMax">The maximum lag at which to estimate the function. A negative value selects 10*log10(N/m) where N is the number of observations.
        /// Lags at or beyond the number of observations give zero.</param>
        /// <param name="type">The type of function to be computed.</param>
        /// <param name="centered">True to remove the sample mean, as <see cref="Function(IList{double}, int, Type)"/> does;
        /// false for a process whose mean is known to be zero, so the autocovariance at lag k is
        /// Σ x_t x_(t-k) / n and the autocorrelation is Σ x_t x_(t-k) / Σ x_t².</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the function evaluated at the given values.
        /// </returns>
        public static double[,]? Function(IList<double> data, int lagMax, Type type, bool centered)
        {
            if (type == Type.Correlation)
            {
                return Correlation(data, lagMax, centered);
            }
            else if (type == Type.Covariance)
            {
                return Covariance(data, lagMax, centered);
            }
            else if (type == Type.Partial)
            {
                return Partial(data, lagMax, centered);
            }

            return null;
        }

        /// <summary>
        /// Computes the autocovariance, autocorrelation, or partial autocorrelation function.
        /// </summary>
        /// <param name="timeSeries">The time-series to evaluate.</param>
        /// <param name="lagMax">The maximum lag at which to estimate the function. Default is 10*log10(N/m) where N is the number of observations.
        /// Lags at or beyond the number of observations give zero.</param>
        /// <param name="type">The type of function to be computed. Default is correlation.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the function evaluated at the given values.
        /// </returns>
        public static double[,]? Function(TimeSeries timeSeries, int lagMax = -1, Type type = Type.Correlation)
        {
            return Function(timeSeries, lagMax, type, true);
        }

        /// <summary>
        /// Computes the autocovariance, autocorrelation, or partial autocorrelation function, about the sample mean or
        /// about zero.
        /// </summary>
        /// <param name="timeSeries">The time-series to evaluate.</param>
        /// <param name="lagMax">The maximum lag at which to estimate the function. A negative value selects 10*log10(N/m) where N is the number of observations.
        /// Lags at or beyond the number of observations give zero.</param>
        /// <param name="type">The type of function to be computed.</param>
        /// <param name="centered">True to remove the sample mean, as <see cref="Function(TimeSeries, int, Type)"/> does;
        /// false for a process whose mean is known to be zero (see <see cref="Function(IList{double}, int, Type, bool)"/>).</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the function evaluated at the given values.
        /// </returns>
        public static double[,]? Function(TimeSeries timeSeries, int lagMax, Type type, bool centered)
        {
            if (type == Type.Correlation)
            {
                return Correlation(timeSeries, lagMax, centered);
            }
            else if (type == Type.Covariance)
            {
                return Covariance(timeSeries, lagMax, centered);
            }
            else if (type == Type.Partial)
            {
                return Partial(timeSeries, lagMax, centered);
            }

            return null;
        }

        /// <summary>
        /// Compute the autocovariance function.
        /// </summary>
        /// <param name="data">The list of data to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the covariance of the given values.
        /// </returns>
        private static double[,]? Covariance(IList<double> data, int lagMax, bool centered)
        {
            int n = data.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            if (!centered) return UncenteredCovariance(LaggedProductSums(data as double[] ?? [.. data], lagMax), n);
            var acvf = new double[lagMax + 1, 2];
            double mean = centered ? Statistics.Mean(data) : 0d;
            for (int lag = 0; lag <= lagMax; lag++)
            {
                double cov = 0d;
                for (int t = lag; t < n; t++)
                    cov += (data[t] - mean) * (data[t - lag] - mean);
                acvf[lag, 0] = lag;
                acvf[lag, 1] = cov / n;
            }

            return acvf;
        }

        /// <summary>
        /// Compute the autocovariance function.
        /// </summary>
        /// <param name="timeSeries">The time-series to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the covariance of the given values.
        /// </returns>
        private static double[,]? Covariance(TimeSeries timeSeries, int lagMax, bool centered)
        {
            int n = timeSeries.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            if (!centered) return UncenteredCovariance(LaggedProductSums(timeSeries.ValuesToArray(), lagMax), n);
            var acvf = new double[lagMax + 1, 2];
            double mean = centered ? timeSeries.MeanValue() : 0d;
            for (int lag = 0; lag <= lagMax; lag++)
            {
                double cov = 0d;
                for (int t = lag; t < n; t++)
                    cov += (timeSeries[t].Value - mean) * (timeSeries[t - lag].Value - mean);
                acvf[lag, 0] = lag;
                acvf[lag, 1] = cov / n;
            }

            return acvf;
        }

        /// <summary>
        /// Compute the autocorrelation function.
        /// </summary>
        /// <param name="data">The list of data to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the autocorrelation of the given values.
        /// </returns>
        private static double[,]? Correlation(IList<double> data, int lagMax, bool centered)
        {
            int n = data.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            if (!centered) return UncenteredCorrelation(LaggedProductSums(data as double[] ?? [.. data], lagMax));
            var acf = Covariance(data, lagMax, centered);
            if (acf == null) return null;
            double den = acf[0, 1];
            if (den == 0) return null;
            for (int i = 0; i < acf.GetLength(0); i++)
                acf[i, 1] /= den;
            return acf;
        }

        /// <summary>
        /// Compute the autocorrelation function.
        /// </summary>
        /// <param name="timeSeries">The time-series to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the autocorrelation of the given values.
        /// </returns>
        private static double[,]? Correlation(TimeSeries timeSeries, int lagMax, bool centered)
        {
            int n = timeSeries.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            if (!centered) return UncenteredCorrelation(LaggedProductSums(timeSeries.ValuesToArray(), lagMax));
            var acf = Covariance(timeSeries, lagMax, centered);
            if (acf == null) return null;
            double den = acf[0, 1];
            if (den == 0) return null;
            for (int i = 0; i < acf.GetLength(0); i++)
                acf[i, 1] /= den;
            return acf;
        }

        /// <summary>
        /// Compute the partial autocorrelation function.
        /// </summary>
        /// <param name="data">The list of data to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the partial autocorrelation of the given values.
        /// </returns>
        private static double[,]? Partial(IList<double> data, int lagMax, bool centered)
        {
            int n = data.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            // First compute the ACVF
            var acvf = Covariance(data, lagMax, centered);
            if (acvf == null) return null;
            // Then compute PACF using the Durbin-Levinson algorithm
            int i, j;
            var phis = new double[lagMax + 1];
            var phis2 = new double[lagMax + 1];
            var pacf = new double[lagMax, 2];
            double vi, phinn;
            phis[0] = acvf[1, 1] / acvf[0, 1];
            pacf[0, 0] = 1d;
            pacf[0, 1] = phis[0];
            vi = acvf[0, 1];
            vi *= 1d - phis[0] * phis[0];
            for (i = 2; i <= lagMax; i++)
            {
                for (j = 0; j < i - 1; j++)
                    phis2[j] = phis[i - j - 2];
                phinn = acvf[i, 1];
                for (j = 1; j < i; j++)
                    phinn -= phis[j - 1] * acvf[i - j, 1];
                phinn /= vi;
                for (j = 0; j < i - 1; j++)
                    phis[j] -= phinn * phis2[j];
                vi *= 1d - phinn * phinn;
                phis[i - 1] = phinn;
                pacf[i - 1, 0] = i;
                pacf[i - 1, 1] = phis[i - 1];
            }

            return pacf;
        }

        /// <summary>
        /// Compute the partial autocorrelation function.
        /// </summary>
        /// <param name="timeSeries">The time-series to assess.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <param name="centered">True to remove the sample mean; false to evaluate about zero.</param>
        /// <returns>
        /// An n x 2 matrix, with n being the number of given input data points. The first column contains the lag and the 
        /// second column contains the partial autocorrelation of the given values.
        /// </returns>
        private static double[,]? Partial(TimeSeries timeSeries, int lagMax, bool centered)
        {
            int n = timeSeries.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2) return null;
            // First compute the ACVF
            var acvf = Covariance(timeSeries, lagMax, centered);
            if (acvf == null) return null;
            // Then compute PACF using the Durbin-Levinson algorithm
            int i, j;
            var phis = new double[lagMax + 1];
            var phis2 = new double[lagMax + 1];
            var pacf = new double[lagMax, 2];
            double vi, phinn;
            phis[0] = acvf[1, 1] / acvf[0, 1];
            pacf[0, 0] = 1d;
            pacf[0, 1] = phis[0];
            vi = acvf[0, 1];
            vi *= 1d - phis[0] * phis[0];
            for (i = 2; i <= lagMax; i++)
            {
                for (j = 0; j < i - 1; j++)
                    phis2[j] = phis[i - j - 2];
                phinn = acvf[i, 1];
                for (j = 1; j < i; j++)
                    phinn -= phis[j - 1] * acvf[i - j, 1];
                phinn /= vi;
                for (j = 0; j < i - 1; j++)
                    phis[j] -= phinn * phis2[j];
                vi *= 1d - phinn * phinn;
                phis[i - 1] = phinn;
                pacf[i - 1, 0] = i;
                pacf[i - 1, 1] = phis[i - 1];
            }

            return pacf;
        }

        /// <summary>
        /// Returns the lagged product sums S_k = Σ_(t=k)^(n-1) x_t x_(t-k) for k = 0 to <paramref name="lagMax"/>.
        /// </summary>
        /// <param name="x">The series values.</param>
        /// <param name="lagMax">The maximum lag.</param>
        /// <returns>The sums S_0 to S_lagMax; a lag at or beyond the length of <paramref name="x"/> has no products and a sum of zero.</returns>
        /// <remarks>
        /// Each sum adds the products in increasing t, the order of the centered autocovariance loop, so S_k / n has
        /// the same bits as that loop evaluated about zero. Four lags are accumulated together over one sweep of a
        /// contiguous array, so the four sums are independent dependency chains and each x_t is read once per block;
        /// the order of every individual sum is unchanged.
        /// </remarks>
        private static double[] LaggedProductSums(double[] x, int lagMax)
        {
            int n = x.Length;
            var sums = new double[lagMax + 1];
            int last = Math.Min(lagMax, n - 1);
            int lag = 0;
            for (; lag + 3 <= last; lag += 4)
            {
                double s0 = 0d, s1 = 0d, s2 = 0d, s3 = 0d;
                // Lead-in: the first terms of the lower lags, before every lag of the block has a term.
                s0 += x[lag] * x[0];
                s0 += x[lag + 1] * x[1];
                s1 += x[lag + 1] * x[0];
                s0 += x[lag + 2] * x[2];
                s1 += x[lag + 2] * x[1];
                s2 += x[lag + 2] * x[0];
                for (int t = lag + 3; t < n; t++)
                {
                    double xt = x[t];
                    s0 += xt * x[t - lag];
                    s1 += xt * x[t - lag - 1];
                    s2 += xt * x[t - lag - 2];
                    s3 += xt * x[t - lag - 3];
                }
                sums[lag] = s0;
                sums[lag + 1] = s1;
                sums[lag + 2] = s2;
                sums[lag + 3] = s3;
            }
            for (; lag <= last; lag++)
            {
                double sum = 0d;
                for (int t = lag; t < n; t++)
                    sum += x[t] * x[t - lag];
                sums[lag] = sum;
            }
            return sums;
        }

        /// <summary>
        /// Returns the autocovariance about zero, S_k / n, as a lag table.
        /// </summary>
        /// <param name="sums">The lagged product sums S_k.</param>
        /// <param name="n">The number of observations.</param>
        /// <returns>The lags in the first column and S_k / n in the second.</returns>
        private static double[,] UncenteredCovariance(double[] sums, int n)
        {
            var acvf = new double[sums.Length, 2];
            for (int lag = 0; lag < sums.Length; lag++)
            {
                acvf[lag, 0] = lag;
                acvf[lag, 1] = sums[lag] / n;
            }
            return acvf;
        }

        /// <summary>
        /// Returns the autocorrelation about zero, S_k / S_0, as a lag table, or null when S_0 is zero.
        /// </summary>
        /// <param name="sums">The lagged product sums S_k.</param>
        /// <returns>The lags in the first column and S_k / S_0 in the second, or null for an all-zero series.</returns>
        /// <remarks>
        /// Each lag sum is divided by S_0 once. Dividing the autocovariances S_k / n by S_0 / n would add two
        /// roundings per lag.
        /// </remarks>
        private static double[,]? UncenteredCorrelation(double[] sums)
        {
            double den = sums[0];
            if (den == 0) return null;
            var acf = new double[sums.Length, 2];
            for (int lag = 0; lag < sums.Length; lag++)
            {
                acf[lag, 0] = lag;
                acf[lag, 1] = sums[lag] / den;
            }
            return acf;
        }

        /// <summary>
        /// Get confidence interval for ACF and PACF rho values.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="interval">The confidence interval width. Default = 0.95, or 95%.</param>
        /// <returns>
        /// An array that contains the lower and upper bound of the confidence interval.
        /// </returns>
        public static double[] CorrelationConfidenceInterval(int sampleSize, double interval = 0.95d)
        {
            double alpha = 0.5d * (1.0d - interval);
            double lo = Normal.StandardZ(alpha) / Math.Sqrt(sampleSize);
            double hi = Normal.StandardZ(1d - alpha) / Math.Sqrt(sampleSize);
            return new[] { lo, hi };
        }
    }
}