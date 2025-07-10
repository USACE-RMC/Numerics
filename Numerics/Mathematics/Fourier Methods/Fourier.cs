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
        /// </remarks>
        public static void FFT(double[] data, bool inverse = false)
        {
            // create variables
            int n, nn, mmax, m, j, istep, i, isign;
            double wtemp, wr, wpr, wpi, wi, theta, tempr, tempi;
            n = (int)(data.Count() / 2d);
            // Check if n is a power of 2.
            if (n < 2 || IsPowerOfTwo(n) == false)
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
        public static double[,] Autocorrelation(IList<double> series, int lagMax = -1)
        {
            int n = series.Count;
            if (lagMax < 0)
                lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            if (lagMax < 1 || n < 2)
                return null;
            // Pad the length to be the power of 2 to facilitate FFT speed.
            int newLength = Convert.ToInt32(Math.Pow(2d, Math.Ceiling(Math.Log(series.Count, 2d))));
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

        /// <summary>
        /// Helper function to determine if integer is power of 2.
        /// </summary>
        /// <param name="n">Integer to test.</param>
        /// <returns> 
        /// True if the integer is a power of 2
        /// </returns>
        private static bool IsPowerOfTwo(int n)
        {
            return (n & n - 1) == 0;
        }
    }
}