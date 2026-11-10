using System;
using System.Collections.Generic;
using System.Linq;
using Numerics.Data.Statistics;

namespace Numerics.Mathematics
{

    /// <summary>
    /// Contains Fourier transform methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors:  </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains shared functions for Fast Fourier Transform (FFT), an algorithm that computes the Discrete Fourier Transform (DFT) which is obtained by 
    /// decomposing a sequence of values into components of different frequencies. The DFT converts a finite sequence of equally-spaced samples of the discrete-time Fourier 
    /// transform (DTFT), a complex valued function of frequency. A FFT rapidly computes these transformations by factorizing the DFT matrix into a product of mostly zero factors.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Fast_Fourier_transform"/>
    /// </description></item>
    /// <item><description>
    /// <see href="https://en.wikipedia.org/wiki/Discrete_Fourier_transform"/>
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class Fourier
    {
        /// <summary>
        /// Performs the fast Fourier transform (FFT) on a complex data array.
        /// </summary>
        /// <param name="data">Data is a complex array of length n stored as a real array of length 2*n. n must be an integer power of 2.</param>
        /// <param name="inverse">Determines whether to perform an inverse transform.</param>
        /// <remarks>
        /// Replaces the data array by its discrete Fourier transform. If inverse is true, then
        /// this routine replaces the data array by n times its inverse discrete Fourier Transform.
        /// The forward transform uses the Numerical Recipes sign convention, H(k) = sum h(j) exp(+2 pi i j k / n),
        /// which is the complex conjugate of the numpy.fft.fft convention for real input.
        /// </remarks>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the length of data is odd, or data.Length / 2 is
        /// not a power of 2 of at least 2.</exception>
        public static void FFT(double[] data, bool inverse = false)
        {
            // create variables
            int n, nn, mmax, m, j, istep, i, isign;
            double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
            if (data.Length % 2 != 0)
                throw new ArgumentOutOfRangeException(nameof(data), "The data array must hold interleaved real and imaginary parts, so its length must be even.");
            n = (int)(data.Count() / 2d);
            // Check if n is a power of 2.
            if (n < 2 || Tools.IsPowerOfTwo(n) == false)
            {
                throw new ArgumentOutOfRangeException(nameof(data), "The data array length must be a power of 2.");
            }
            nn = n << 1;
            isign = inverse == false ? 1 : -1;
            j = 1;
            // This is the bit-reversal section of the routine.
            for (i = 1; i < nn; i += 2)
            {
                if (j > i)
                {
                    // swap values
                    tempr = data[j - 1];
                    tempi = data[j];
                    data[j - 1] = data[i - 1];
                    data[j] = data[i];
                    data[i - 1] = tempr;
                    data[i] = tempi;
                }
                m = n;
                while (m >= 2 && j > m)
                {
                    j -= m;
                    m >>= 1;
                }
                j += m;
            }
            // Here begins the Danielson-Lanczos section of the routine.
            mmax = 2;
            // Outer loop executed log2 n times
            while (nn > mmax)
            {
                istep = mmax << 1;
                // Initialize the trigonometric recurrence.
                theta = isign * (2d * Math.PI / mmax);
                wpr = -2.0d * Math.Pow(Math.Sin(0.5d * theta), 2d);
                wpi = Math.Sin(theta);
                wr = 1.0d;
                wi = 0.0d;
                // Here are the two nested inner loops.
                for (m = 1; m < mmax; m += 2)
                {
                    for (i = m; istep >= 0 ? i <= nn : i >= nn; i += istep)
                    {
                        j = i + mmax;
                        // This is the Danielson-Lanczos formula:
                        tempr = wr * data[j - 1] - wi * data[j];
                        tempi = wr * data[j] + wi * data[j - 1];
                        data[j - 1] = data[i - 1] - tempr;
                        data[j] = data[i] - tempi;
                        data[i - 1] += tempr;
                        data[i] += tempi;
                    }
                    // Trigonometric recurrence
                    wtemp = wr;
                    wr = wr * wpr - wi * wpi + wr;
                    wi = wi * wpr + wtemp * wpi + wi;
                }

                mmax = istep;
            }
        }

        /// <summary>
        /// Calculates the Fourier transform of a set of n real-valued data points.
        /// </summary>
        /// <param name="data">Data is a real array of length n. n must be an integer power of 2.</param>
        /// <param name="inverse">Determines whether to perform an inverse transform.</param>
        /// <remarks>
        /// Replaces these data by the positive frequency half of their complex Fourier transform.
        /// The real-valued first and last components of the complex transform are returned as
        /// elements data(0) and data(1), respectively. This routine also calculates the inverse transform of a complex data array
        /// if it is the transform of real data. (Result in this case must be multiplied by 2/n.)
        /// </remarks>
        public static void RealFFT(double[] data, bool inverse = false)
        {
            int i, i1, i2, i3, i4, n, isign;
            double c1, c2, h1r, h1i, h2r, h2i, wr, wi, wpr, wpi, wtemp, theta;
            n = data.Count();
            isign = inverse == false ? 1 : -1;
            // Initial the recurrence
            theta = Math.PI / (n >> 1);
            c1 = 0.5d;
            if (isign == 1)
            {
                c2 = -0.5d;
                // The forward transform is here.
                FFT(data);
            }
            else
            {
                // Otherwise set up for an inverse transform.
                c2 = 0.5d;
                theta = -theta;
            }
            wpr = -2.0d * Math.Pow(Math.Sin(0.5d * theta), 2d);
            wpi = Math.Sin(theta);
            wr = 1.0d + wpr;
            wi = wpi;
            // Case i = 0 done separately below.
            for (i = 1; i < (n >> 2); i++)
            {
                i1 = i + i;
                i2 = 1 + i1;
                i3 = n - i1;
                i4 = 1 + i3;
                // The two separate transforms are separated out of data.
                h1r = c1 * (data[i1] + data[i3]);
                h1i = c1 * (data[i2] - data[i4]);
                h2r = -c2 * (data[i2] + data[i4]);
                h2i = c2 * (data[i1] - data[i3]);
                // Here they are recombined to form the true transform of 
                // the original real data.
                data[i1] = h1r + wr * h2r - wi * h2i;
                data[i2] = h1i + wr * h2i + wi * h2r;
                data[i3] = h1r - wr * h2r + wi * h2i;
                data[i4] = -h1i + wr * h2i + wi * h2r;
                // The recurrence.
                wtemp = wr;
                wr = wr * wpr - wi * wpi + wr;
                wi = wi * wpr + wtemp * wpi + wi;
            }

            if (isign == 1)
            {
                // Squeeze the first and last data together to get them 
                // all within the original array.
                h1r = data[0];
                data[0] = h1r + data[1];
                data[1] = h1r - data[1];
            }
            else
            {
                // This is the inverse transform for the case isign = -1.
                h1r = data[0];
                data[0] = c1 * (h1r + data[1]);
                data[1] = c1 * (h1r - data[1]);
                FFT(data, true);
            }
        }

        /// <summary>
        /// Performs the discrete Fourier transform of a complex array of any length, using Bluestein's algorithm
        /// when the length is not a power of 2.
        /// </summary>
        /// <param name="data">A complex array of length n, stored as a real array of length 2*n with interleaved real
        /// and imaginary parts. n can be any positive integer. The array is replaced by its transform.</param>
        /// <param name="inverse">Determines whether to perform an inverse transform.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="data"/> is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the array length is zero or odd.</exception>
        /// <remarks>
        /// <para>
        /// Uses the same conventions as <see cref="FFT(double[], bool)"/>: the forward transform is
        /// H(k) = sum h(j) exp(+2 pi i j k / n), and the inverse transform returns n times the inverse DFT. When n is a
        /// power of 2 (n &gt;= 2) the call is passed to <see cref="FFT(double[], bool)"/> unchanged.
        /// </para>
        /// <para>
        /// Otherwise the identity j k = (j^2 + k^2 - (k - j)^2) / 2 turns the DFT into a convolution with a chirp
        /// exp(-+i pi m^2 / n), which is evaluated with power-of-2 FFTs of length at least 2n - 1 (Bluestein's
        /// algorithm). The cost is O(n log n), about four to six times that of a power-of-2 FFT of similar length, and
        /// the round-off error grows like that of an FFT of the padded length. The chirp phase is reduced modulo 2n in
        /// integer arithmetic before scaling, so it stays accurate for large n.
        /// </para>
        /// <para>
        /// <b> References: </b> Bluestein, L. I. (1970). A linear filtering approach to the computation of discrete Fourier
        /// transform. IEEE Transactions on Audio and Electroacoustics, 18(4), 451-455. https://doi.org/10.1109/TAU.1970.1162132
        /// </para>
        /// </remarks>
        public static void BluesteinFFT(double[] data, bool inverse = false)
        {
            if (data == null) throw new ArgumentNullException(nameof(data));
            if (data.Length == 0 || data.Length % 2 != 0)
                throw new ArgumentOutOfRangeException(nameof(data), "The data array must hold n >= 1 interleaved complex values, so its length must be even and positive.");
            int n = data.Length / 2;
            if (n == 1) return;
            if (Tools.IsPowerOfTwo(n))
            {
                FFT(data, inverse);
                return;
            }

            // Chirp c(j) = exp(sign * i * pi * j^2 / n), with the sign of the transform kernel.
            double sign = inverse ? -1d : 1d;
            var chirpRe = new double[n];
            var chirpIm = new double[n];
            long twoN = 2L * n;
            for (int j = 0; j < n; j++)
            {
                long phase = (long)j * j % twoN;
                double angle = Math.PI * phase / n;
                chirpRe[j] = Math.Cos(angle);
                chirpIm[j] = sign * Math.Sin(angle);
            }

            int m = 1;
            while (m < 2 * n - 1)
                m <<= 1;
            // a(j) = h(j) * c(j), zero-padded to m.
            var a = new double[2 * m];
            for (int j = 0; j < n; j++)
            {
                double hr = data[2 * j], hi = data[2 * j + 1];
                a[2 * j] = hr * chirpRe[j] - hi * chirpIm[j];
                a[2 * j + 1] = hr * chirpIm[j] + hi * chirpRe[j];
            }
            // b(k) = conj(c(k)) for k = -(n - 1)..(n - 1), stored circularly.
            var b = new double[2 * m];
            b[0] = chirpRe[0];
            b[1] = -chirpIm[0];
            for (int k = 1; k < n; k++)
            {
                b[2 * k] = b[2 * (m - k)] = chirpRe[k];
                b[2 * k + 1] = b[2 * (m - k) + 1] = -chirpIm[k];
            }

            // Circular convolution of a and b by power-of-2 FFTs; the unnormalized inverse is divided by m.
            FFT(a);
            FFT(b);
            for (int k = 0; k < m; k++)
            {
                double ar = a[2 * k], ai = a[2 * k + 1], br = b[2 * k], bi = b[2 * k + 1];
                a[2 * k] = ar * br - ai * bi;
                a[2 * k + 1] = ar * bi + ai * br;
            }
            FFT(a, true);

            // H(k) = c(k) * (a * b)(k).
            for (int k = 0; k < n; k++)
            {
                double cr = a[2 * k] / m, ci = a[2 * k + 1] / m;
                data[2 * k] = cr * chirpRe[k] - ci * chirpIm[k];
                data[2 * k + 1] = cr * chirpIm[k] + ci * chirpRe[k];
            }
        }

        /// <summary>
        /// Computes the full linear convolution of two real sequences by FFT.
        /// </summary>
        /// <param name="data1">The first real sequence, of any length n1 &gt;= 1.</param>
        /// <param name="data2">The second real sequence, of any length n2 &gt;= 1.</param>
        /// <returns>
        /// The convolution c(k) = sum over j of data1(j) * data2(k - j), for k = 0..n1 + n2 - 2 (length n1 + n2 - 1),
        /// the same as numpy.convolve(data1, data2, mode="full").
        /// </returns>
        /// <exception cref="ArgumentNullException">Thrown when either sequence is null.</exception>
        /// <exception cref="ArgumentException">Thrown when either sequence is empty.</exception>
        /// <remarks>
        /// <para>
        /// Both sequences are zero-padded to a power of 2 of at least n1 + n2 - 1, so there is no wraparound, and are
        /// transformed with <see cref="RealFFT(double[], bool)"/>, which costs about half of a complex FFT of the same
        /// length. The cost is O(N log N) for N = n1 + n2 - 1, against O(n1 n2) for the direct sum.
        /// </para>
        /// <para>
        /// The error is absolute: about machine epsilon times log2(N) times the largest |c(k)|, rather than relative
        /// to each term. Entries far smaller than the largest output (for example deep tails of a convolved
        /// probability mass function) can therefore carry round-off of order 1E-16 times the peak, and exact zeros can
        /// come back as tiny positive or negative values. Clamp or use the direct sum when that matters.
        /// </para>
        /// </remarks>
        public static double[] Convolve(IList<double> data1, IList<double> data2)
        {
            if (data1 == null) throw new ArgumentNullException(nameof(data1));
            if (data2 == null) throw new ArgumentNullException(nameof(data2));
            if (data1.Count == 0 || data2.Count == 0)
                throw new ArgumentException("Both sequences must contain at least one value.");
            int length = data1.Count + data2.Count - 1;
            // RealFFT needs at least two complex points, that is four real values.
            int n = 4;
            while (n < length)
                n <<= 1;
            var x = new double[n];
            var y = new double[n];
            for (int i = 0; i < data1.Count; i++) x[i] = data1[i];
            for (int i = 0; i < data2.Count; i++) y[i] = data2[i];
            RealFFT(x);
            RealFFT(y);
            // Multiply the packed half spectra; data[0] and data[1] hold the real DC and Nyquist terms.
            // The inverse RealFFT returns n/2 times the sequence, so scale by 2/n.
            double scale = 2d / n;
            x[0] *= y[0] * scale;
            x[1] *= y[1] * scale;
            for (int i = 2; i < n; i += 2)
            {
                double xr = x[i], xi = x[i + 1], yr = y[i], yi = y[i + 1];
                x[i] = (xr * yr - xi * yi) * scale;
                x[i + 1] = (xr * yi + xi * yr) * scale;
            }
            RealFFT(x, true);
            var result = new double[length];
            Array.Copy(x, result, length);
            return result;
        }

        /// <summary>
        /// Computes the correlation of two real datasets. This is NOT a normalized correlation coefficient.
        /// </summary>
        /// <param name="data1">Data 1 is a real array of length n. n must be an integer power of 2.</param>
        /// <param name="data2">Data 2 is a real array of length n. n must be an integer power of 2.</param>
        /// <remarks>
        /// The answer is returned in wraparound order, i.e., correlations at increasing negative lags
        /// are in [n-1] on down to  [n/2], while correlations at increasingly positive lags are in [0]
        /// up to [n/2-1].
        /// </remarks>
        /// <returns>
        /// An array of the correlation between the two datasets
        /// </returns>
        public static double[] Correlation(double[] data1, double[] data2)
        {
            int no2, i;
            int n = data1.Count();
            double tmp;
            var temp = new double[n];
            var ans = new double[n];
            for (i = 0; i < n; i++)
            {
                ans[i] = data1[i];
                temp[i] = data2[i];
            }
            // Transform both data vectors
            RealFFT(ans);
            RealFFT(temp);
            // Normalization for inverse FFT.
            // Multiply to find FFT of their correlation.
            no2 = n >> 1;
            for (i = 2; i < n; i += 2)
            {
                tmp = ans[i];
                ans[i] = (ans[i] * temp[i] + ans[i + 1] * temp[i + 1]) / no2;
                ans[i + 1] = (ans[i + 1] * temp[i] - tmp * temp[i + 1]) / no2;
            }
            ans[0] = ans[0] * temp[0] / no2;
            ans[1] = ans[1] * temp[1] / no2;
            // Inverse transform gives correlation.
            RealFFT(ans, true);
            return ans;
        }

        /// <summary>
        /// Computes the autocorrelation function (ACF) given a series of data.
        /// </summary>
        /// <param name="series">The series of data to assess.</param>
        /// <param name="lagMax">The maximum lag. The first lag begins at zero.</param>
        /// <returns>
        /// A 2-column array with the first column contains the lags, and the second the autocorrelation.
        /// </returns>
        /// <remarks>
        /// Returns the standard (linear) sample autocorrelation, r(k) = sum over t of (x_t - mean)(x_(t+k) - mean)
        /// divided by the lag-0 sum, computed by FFT on a zero-padded copy of length at least 2n, so the circular correlation equals the
        /// linear one at every lag.
        /// </remarks>
        public static double[,]? Autocorrelation(IList<double> series, int lagMax = -1)
        {
            int n = series.Count;
            if (lagMax < 0)
                lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2)
                return null!;
            // Zero-pad to a power of 2 of at least 2n, so the FFT's circular correlation equals the linear
            // correlation at every lag and no wraparound reaches the reported lags.
            int newLength = 1;
            while (newLength < 2 * n)
                newLength <<= 1;
            // Normalize the data series
            var normalizedSeries = new double[newLength];
            double mean = Statistics.Mean(series);
            for (int i = 0; i < newLength; i++)
            {
                if (i < series.Count)
                {
                    normalizedSeries[i] = series[i] - mean;
                }
            }
            // Get correlation function
            var correlation = Correlation(normalizedSeries, normalizedSeries);
            // Convert to normalized correlation coefficient
            double maxValue = correlation.Max();
            var acf = new double[lagMax + 1, 2];
            for (int i = 0; i <= lagMax; i++)
            {
                acf[i, 0] = i;
                acf[i, 1] = correlation[i] / maxValue;
            }

            return acf;
        }

    }
}