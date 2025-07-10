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

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// Contains plotting position formulas.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// The choice of quantiles from a theoretical distribution can depend upon context and purpose.
    /// Several different formulas have been used or proposed as affine symmetrical plotting positions.
    /// Such formulas have the form (k − a) / (n + 1 − 2a) for some value of a in the range from 0 to 1,
    /// which gives a range between k / (n + 1) and (k − 1) / (n - 1).
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot"/>
    /// </remarks>
    public class PlottingPositions
    {
        /// <summary>
        /// The general plotting position formula. This formula assumes the data is uncensored and complete.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <param name="alpha">The alpha coefficient. Range (0,1).</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Function(int N, double alpha)
        {
            if (N <= 2)
            {
                throw new ArgumentOutOfRangeException(nameof(N), "The sample size N must be greater than 2.");
            }

            if (alpha < 0d || alpha > 1d)
            {
                throw new ArgumentOutOfRangeException(nameof(alpha), "The alpha coefficient must be between 0 and 1.");
            }

            var PP = new double[N];
            for (int i = 1; i <= N; i++)
                PP[i - 1] = (i - alpha) / (N + 1 - 2d * alpha);
            return PP;
        }

        /// <summary>
        /// The general plotting position formula. This formula assumes the data is uncensored and complete.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <param name="plottingPostionType">The plotting position formula type.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Function(int N, PlottingPositions.PlottingPostionType plottingPostionType)
        {
            if (plottingPostionType == PlottingPostionType.Blom)
            {
                return Blom(N);
            }
            else if (plottingPostionType == PlottingPostionType.Cunnane)
            {
                return Cunnane(N);
            }
            else if (plottingPostionType == PlottingPostionType.Gringorten)
            {
                return Gringorten(N);
            }
            else if (plottingPostionType == PlottingPostionType.Hazen)
            {
                return Hazen(N);
            }
            else if (plottingPostionType == PlottingPostionType.KaplanMeier)
            {
                return KaplanMeier(N);
            }
            else if (plottingPostionType == PlottingPostionType.Median)
            {
                return Median(N);
            }
            else if (plottingPostionType == PlottingPostionType.Weibull)
            {
                return Weibull(N);
            }
            return null;
        }

        /// <summary>
        /// The Weibull plotting position formula (alpha = 0.0). Recommended for uniform distribution.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Weibull(int N)
        {
            return Function(N, 0d);
        }

        /// <summary>
        /// The Median plotting position formula (alpha = 0.3175). Provides median exceedance probabilities for all distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Median(int N)
        {
            return Function(N, 0.3175d);
        }

        /// <summary>
        /// The Blom (1958) plotting position formula (alpha = 0.375). Recommended for Normal, Gamma, 2-parameter Log Normal,
        /// 3-parameter Log Normal, and Log Pearson Type III distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Blom(int N)
        {
            return Function(N, 0.375d);
        }

        /// <summary>
        /// The Cunnane (1978) plotting position formula (alpha = 0.40). Recommended for GEV and Log-Gumbel distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Cunnane(int N)
        {
            return Function(N, 0.4d);
        }

        /// <summary>
        /// The Gringorten (1963) plotting position formula (alpha = 0.44). Recommended for Exponential, Gumbel and Weibull distributions.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Gringorten(int N)
        {
            return Function(N, 0.44d);
        }

        /// <summary>
        /// The Hazen plotting position formula (alpha = 0.50). Recommended when the parameters of the parent distribution
        /// are unknown.
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] Hazen(int N)
        {
            return Function(N, 0.5d);
        }

        /// <summary>
        /// The Kaplan-Meier plotting position formula for uncensored data (k/n).
        /// </summary>
        /// <param name="N">The sample size.</param>
        /// <returns>An array of plotting positions of size N.</returns>
        public static double[] KaplanMeier(int N)
        {
            if (N <= 2)
            {
                throw new ArgumentOutOfRangeException(nameof(N), "The sample size N must be greater than 2.");
            }
            var PP = new double[N];
            for (int i = 1; i <= N; i++)
                PP[i - 1] = i / (double)N;
            return PP;
        }

        /// <summary>
        /// Enumeration of plotting position types.
        /// </summary>
        public enum PlottingPostionType
        {
            /// <summary>
            /// The Weibull plotting position formula (alpha = 0.0).
            /// </summary>
            Weibull,

            /// <summary>
            /// The Median plotting position formula (alpha = 0.3175).
            /// </summary>
            Median,

            /// <summary>
            /// The Blom (1958) plotting position formula (alpha = 0.375).
            /// </summary>
            Blom,

            /// <summary>
            /// The Cunnane (1978) plotting position formula (alpha = 0.40).
            /// </summary>
            Cunnane,

            /// <summary>
            /// The Gringorten (1963) plotting position formula (alpha = 0.44).
            /// </summary>
            Gringorten,

            /// <summary>
            /// The Hazen plotting position formula (alpha = 0.50).
            /// </summary>
            Hazen,

            /// <summary>
            /// The Kaplan-Meier plotting position formula for uncensored data (k/n).
            /// </summary>
            KaplanMeier
        }
    }
}