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

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains methods for determining the correlation of two data sets.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// This class contains a shared functions for Pearson and Spearman Correlation Coefficients. The Pearson correlation coefficient measures linear correlation
    /// between two sets of data. The Spearman correlation coefficient is a nonparametric measure of rank correlation (i.e. measures the strength and direction of
    /// association between two ranked variables).
    /// </para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item>
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// </item>
    /// <item>
    /// <see href="https://en.wikipedia.org/wiki/Pearson_correlation_coefficient"/>
    /// </item>
    /// <item>
    /// <see href="https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient"/>
    /// </item>
    /// </list>
    /// </remarks>
    public sealed class Correlation
    {

        /// <summary>
        /// Computes the Pearson correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Pearson correlation coefficient.</returns>
        public static double Pearson(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The samples arrays must be the same length.");
            }
            // 
            int i, n = sample1.Count;
            // Find means.
            double ax = 0d;
            double ay = 0d;
            for (i = 0; i <= n - 1; i++)
            {
                ax += sample1[i];
                ay += sample2[i];
            }
            ax /= n;
            ay /= n;
            // Compute the correlation coefficient
            double sxx = 0d, syy = 0d, sxy = 0d;
            double xt, yt;
            for (i = 0; i <= n - 1; i++)
            {
                xt = sample1[i] - ax;
                yt = sample2[i] - ay;
                sxx += xt * xt;
                syy += yt * yt;
                sxy += xt * yt;
            }
            return sxy / Math.Sqrt(sxx * syy);
        }

        /// <summary>
        /// Computes the Spearman ranked correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Spearman ranked correlation coefficient.</returns>
        public static double Spearman(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The sample arrays must be the same length.");
            }

            var rank1 = Statistics.RanksInPlace(sample1.ToArray());
            var rank2 = Statistics.RanksInPlace(sample2.ToArray());
            return Pearson(rank1, rank2);
        }

        /// <summary>
        /// Computes Kendall's Tau ranked correlation coefficient.
        /// </summary>
        /// <param name="sample1">Sample data 1.</param>
        /// <param name="sample2">Sample data 2.</param>
        /// <returns>The Kendall's Tau ranked correlation coefficient.</returns>
        public static double KendallsTau(IList<double> sample1, IList<double> sample2)
        {
            // Check if the arrays are the same length
            if (sample2.Count != sample1.Count)
            {
                throw new ArgumentOutOfRangeException(nameof(sample1), "The sample arrays must be the same length.");
            }

            int n = sample1.Count;
            int j, k, i = 0, n1 = 0, n2 = 0;
            double aa, a1, a2;
            // Loop over first member of pair and second member.
            for (j = 0; j <= n - 2; j++)
            {
                for (k = j + 1; k <= n - 1; k++)
                {
                    a1 = sample1[j] - sample1[k];
                    a2 = sample2[j] - sample2[k];
                    aa = a1 * a2;
                    if (aa != 0.0d)
                    {
                        // Neither has a tie
                        n1 += 1;
                        n2 += 1;
                        i = aa > 0.0 ? ++i : --i;
                    } else {
                        // One or both arrays have ties.
                        if (a1 != 0.0d) n1 += 1;
                        if (a2 != 0.0d) n2 += 1;
                    }
                }
            }
            return i / (Math.Sqrt(n1) * Math.Sqrt(n2));
        }
      
    }
}