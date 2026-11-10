using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using Numerics.Data.Statistics;

namespace Numerics.Distributions.Regional
{
    /// <summary>
    /// The sites of a region for regional frequency analysis by L-moments: each site's name, record length and sample
    /// L-moments.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Each site holds its record length n, its sample mean l1 (the index flood), its L-CV t = l2 / l1 and its
    /// L-moment ratios t3, t4 and t5, the layout of the regional data of Hosking's LMOMENTS Fortran package and of the
    /// R package lmomRFA. The data are built either from each site's raw sample, whose L-moments are the unbiased
    /// sample L-moments of <see cref="Statistics.LinearMoments(IList{double}, int)"/>, or from precomputed L-moments.
    /// </para>
    /// <para>
    /// <see cref="RegionalAverage()"/> is the record-length weighted regional average of the L-moment ratios
    /// (Hosking and Wallis, 1997, section 4.1), and <see cref="Discordancy()"/> the discordancy measure D_i of
    /// section 3.2.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item><description>
    /// Hosking, J. R. M. and Wallis, J. R. (1993). Some statistics useful in regional frequency analysis. Water
    /// Resources Research, 29(2), 271-281. https://doi.org/10.1029/92WR01980
    /// </description></item>
    /// <item><description>
    /// Hosking, J. R. M. and Wallis, J. R. (1997). Regional Frequency Analysis: An Approach Based on L-Moments.
    /// Cambridge University Press. https://doi.org/10.1017/CBO9780511529443
    /// </description></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class RegionalData
    {
        /// <summary>
        /// The upper 10% critical values of the discordancy measure for regions of 1 to 18 sites, as tabulated by
        /// Hosking and Wallis (1997, Table 3.1) and in lmomRFA; 3 from 15 sites on.
        /// </summary>
        private static readonly double[] Critical10 =
            [3, 3, 3, 3, 1.3330, 1.6481, 1.9166, 2.1401, 2.3287, 2.4906, 2.6321, 2.7573, 2.8694, 2.9709, 3, 3, 3, 3];

        /// <summary>
        /// The upper 5% critical values of the discordancy measure for regions of 1 to 18 sites, as in lmomRFA; 4 from
        /// 19 sites on.
        /// </summary>
        private static readonly double[] Critical5 =
            [4, 4, 4, 4, 1.3333, 1.6648, 1.9821, 2.2728, 2.5337, 2.7666, 2.9748, 3.1620, 3.3310, 3.4844, 3.6246, 3.7532, 3.8718, 3.9816];

        private readonly string[] _names;
        private readonly int[] _recordLengths;
        private readonly double[] _mean, _lcv, _t3, _t4, _t5;

        /// <summary>
        /// Constructs regional data from each site's precomputed L-moments.
        /// </summary>
        /// <param name="names">The site names.</param>
        /// <param name="recordLengths">The record length of each site, at least 1.</param>
        /// <param name="mean">The sample mean l1 of each site (the index flood).</param>
        /// <param name="lCV">The sample L-CV t = l2 / l1 of each site.</param>
        /// <param name="tau3">The sample L-skewness t3 of each site.</param>
        /// <param name="tau4">The sample L-kurtosis t4 of each site.</param>
        /// <param name="tau5">Optional. The sample t5 of each site; null leaves every t5 as NaN.</param>
        /// <exception cref="ArgumentNullException">Thrown when a required list is null.</exception>
        /// <exception cref="ArgumentException">Thrown when there are no sites or the lists differ in length.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a record length is below 1, a mean, L-CV, t3 or
        /// t4 is not finite, an L-CV is not positive, or a t5 is infinite.</exception>
        /// <remarks>
        /// The values are copied. A t5 may be NaN, which leaves that site out of the regional average t5. A null
        /// name is replaced by the site's one-based position.
        /// </remarks>
        public RegionalData(IList<string> names, IList<int> recordLengths, IList<double> mean, IList<double> lCV,
            IList<double> tau3, IList<double> tau4, IList<double>? tau5 = null)
        {
            if (names == null) throw new ArgumentNullException(nameof(names));
            if (recordLengths == null) throw new ArgumentNullException(nameof(recordLengths));
            if (mean == null) throw new ArgumentNullException(nameof(mean));
            if (lCV == null) throw new ArgumentNullException(nameof(lCV));
            if (tau3 == null) throw new ArgumentNullException(nameof(tau3));
            if (tau4 == null) throw new ArgumentNullException(nameof(tau4));
            int count = names.Count;
            if (count == 0) throw new ArgumentException("A region needs at least one site.", nameof(names));
            if (recordLengths.Count != count || mean.Count != count || lCV.Count != count || tau3.Count != count ||
                tau4.Count != count || (tau5 != null && tau5.Count != count))
                throw new ArgumentException("Every site list must have one entry per site.", nameof(names));

            _names = new string[count];
            _recordLengths = new int[count];
            _mean = new double[count];
            _lcv = new double[count];
            _t3 = new double[count];
            _t4 = new double[count];
            _t5 = new double[count];
            for (int i = 0; i < count; i++)
            {
                if (recordLengths[i] < 1)
                    throw new ArgumentOutOfRangeException(nameof(recordLengths), "Every record length must be at least 1.");
                if (!Tools.IsFinite(mean[i]) || !Tools.IsFinite(lCV[i]) || !Tools.IsFinite(tau3[i]) || !Tools.IsFinite(tau4[i]))
                    throw new ArgumentOutOfRangeException(nameof(mean), $"The mean, L-CV, t3 and t4 of site {i} must be finite.");
                if (lCV[i] <= 0d)
                    throw new ArgumentOutOfRangeException(nameof(lCV), $"The L-CV of site {i} must be positive.");
                if (tau5 != null && double.IsInfinity(tau5[i]))
                    throw new ArgumentOutOfRangeException(nameof(tau5), $"The t5 of site {i} must be finite or NaN.");
                _names[i] = names[i] ?? (i + 1).ToString(CultureInfo.InvariantCulture);
                _recordLengths[i] = recordLengths[i];
                _mean[i] = mean[i];
                _lcv[i] = lCV[i];
                _t3[i] = tau3[i];
                _t4[i] = tau4[i];
                _t5[i] = tau5 == null ? double.NaN : tau5[i];
            }
        }

        /// <summary>
        /// Constructs regional data from each site's raw sample.
        /// </summary>
        /// <param name="samples">The sample of each site. NaN values mark missing years and are left out.</param>
        /// <param name="names">Optional. The site names; null names the sites "1", "2", and so on.</param>
        /// <returns>The regional data, with each site's L-moments from <see cref="Statistics.LinearMoments(IList{double}, int)"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="samples"/> or a sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when there are no samples, the names differ in number from the
        /// samples, or a sample holds an infinity, has fewer than four values, has a nonpositive mean or is constant
        /// (an L-CV of zero).</exception>
        /// <remarks>
        /// The record length of a site is its number of non-NaN values. The L-moments are the unbiased sample L-moments
        /// of Hosking (1990), as in the SAMLMU routine of LMOMENTS and in lmomRFA's regsamlmu; a site of exactly four
        /// values has t5 = NaN.
        /// </remarks>
        public static RegionalData FromSamples(IList<IList<double>> samples, IList<string>? names = null)
        {
            if (samples == null) throw new ArgumentNullException(nameof(samples));
            int count = samples.Count;
            if (count == 0) throw new ArgumentException("A region needs at least one site.", nameof(samples));
            if (names != null && names.Count != count)
                throw new ArgumentException("There must be one name per sample.", nameof(names));
            var siteNames = new string[count];
            var n = new int[count];
            var l1 = new double[count];
            var lcv = new double[count];
            var t3 = new double[count];
            var t4 = new double[count];
            var t5 = new double[count];
            for (int i = 0; i < count; i++)
            {
                var sample = samples[i] ?? throw new ArgumentNullException(nameof(samples), $"The sample of site {i} is null.");
                var values = new List<double>(sample.Count);
                for (int j = 0; j < sample.Count; j++)
                {
                    double x = sample[j];
                    if (double.IsNaN(x)) continue;
                    if (double.IsInfinity(x)) throw new ArgumentException($"The sample of site {i} holds an infinity.", nameof(samples));
                    values.Add(x);
                }
                if (values.Count < 4)
                    throw new ArgumentException($"The sample of site {i} has fewer than four values.", nameof(samples));
                var moments = Statistics.LinearMoments(values, values.Count >= 5 ? 5 : 4);
                if (!(moments[0] > 0d))
                    throw new ArgumentException($"The sample of site {i} has a nonpositive mean, so its L-CV is undefined.", nameof(samples));
                siteNames[i] = names?[i] ?? (i + 1).ToString(CultureInfo.InvariantCulture);
                n[i] = values.Count;
                l1[i] = moments[0];
                lcv[i] = moments[1] / moments[0];
                t3[i] = moments[2];
                t4[i] = moments[3];
                t5[i] = moments.Length > 4 ? moments[4] : double.NaN;
            }
            return new RegionalData(siteNames, n, l1, lcv, t3, t4, t5);
        }

        /// <summary>
        /// Gets the number of sites.
        /// </summary>
        public int SiteCount => _names.Length;

        /// <summary>
        /// Gets the site names.
        /// </summary>
        public ReadOnlyCollection<string> Names => Array.AsReadOnly(_names);

        /// <summary>
        /// Gets the record length of each site.
        /// </summary>
        public ReadOnlyCollection<int> RecordLengths => Array.AsReadOnly(_recordLengths);

        /// <summary>
        /// Gets the sample mean l1 of each site, its index flood.
        /// </summary>
        public ReadOnlyCollection<double> Mean => Array.AsReadOnly(_mean);

        /// <summary>
        /// Gets the sample L-CV t = l2 / l1 of each site.
        /// </summary>
        public ReadOnlyCollection<double> LCV => Array.AsReadOnly(_lcv);

        /// <summary>
        /// Gets the sample L-skewness t3 of each site.
        /// </summary>
        public ReadOnlyCollection<double> Tau3 => Array.AsReadOnly(_t3);

        /// <summary>
        /// Gets the sample L-kurtosis t4 of each site.
        /// </summary>
        public ReadOnlyCollection<double> Tau4 => Array.AsReadOnly(_t4);

        /// <summary>
        /// Gets the sample t5 of each site; NaN where it is not available.
        /// </summary>
        public ReadOnlyCollection<double> Tau5 => Array.AsReadOnly(_t5);

        /// <summary>
        /// Returns the sample L-moments of one site.
        /// </summary>
        /// <param name="site">The zero-based site index.</param>
        /// <returns>{l1, l2, t3, t4, t5}, with l2 = l1 t.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="site"/> is not a site index.</exception>
        public double[] SiteLinearMoments(int site)
        {
            if (site < 0 || site >= SiteCount) throw new ArgumentOutOfRangeException(nameof(site));
            return [_mean[site], _mean[site] * _lcv[site], _t3[site], _t4[site], _t5[site]];
        }

        /// <summary>
        /// Returns the record-length weighted regional average L-moments.
        /// </summary>
        /// <returns>{1, t^R, t3^R, t4^R, t5^R}: the regional L-moments scaled to l1 = 1, so the second entry is both
        /// the regional L-CV and the regional L-scale.</returns>
        /// <remarks>
        /// Hosking and Wallis (1997), equation 4.3: t^R = Σ n_i t^(i) / Σ n_i, and likewise for t3, t4 and t5, as in
        /// lmomRFA's regavlmom and the REGTST routine. Sites whose t5 is NaN are left out of t5^R, which is NaN when
        /// no site has a t5.
        /// </remarks>
        public double[] RegionalAverage()
        {
            var weights = new double[SiteCount];
            for (int i = 0; i < SiteCount; i++) weights[i] = _recordLengths[i];
            return WeightedAverage(weights);
        }

        /// <summary>
        /// Returns the regional average L-moments with the given site weights.
        /// </summary>
        /// <param name="weights">One nonnegative weight per site, not all zero.</param>
        /// <returns>{1, t^R, t3^R, t4^R, t5^R}, as for <see cref="RegionalAverage()"/>.</returns>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="weights"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the weights differ in number from the sites, a weight is
        /// negative or not finite, or every weight is zero.</exception>
        public double[] RegionalAverage(IList<double> weights)
        {
            if (weights == null) throw new ArgumentNullException(nameof(weights));
            if (weights.Count != SiteCount) throw new ArgumentException("There must be one weight per site.", nameof(weights));
            var w = new double[SiteCount];
            double total = 0d;
            for (int i = 0; i < SiteCount; i++)
            {
                if (!Tools.IsFinite(weights[i]) || weights[i] < 0d)
                    throw new ArgumentException("Every weight must be finite and nonnegative.", nameof(weights));
                w[i] = weights[i];
                total += w[i];
            }
            if (!(total > 0d)) throw new ArgumentException("At least one weight must be positive.", nameof(weights));
            return WeightedAverage(w);
        }

        /// <summary>
        /// Forms the weighted means of the ratios; t5 is averaged over the sites that have one.
        /// </summary>
        /// <param name="w">The site weights.</param>
        /// <returns>{1, t^R, t3^R, t4^R, t5^R}.</returns>
        private double[] WeightedAverage(double[] w)
        {
            double sum = 0d, cv = 0d, s3 = 0d, s4 = 0d, s5 = 0d, sum5 = 0d;
            for (int i = 0; i < SiteCount; i++)
            {
                sum += w[i];
                cv += w[i] * _lcv[i];
                s3 += w[i] * _t3[i];
                s4 += w[i] * _t4[i];
                if (!double.IsNaN(_t5[i]))
                {
                    s5 += w[i] * _t5[i];
                    sum5 += w[i];
                }
            }
            return [1d, cv / sum, s3 / sum, s4 / sum, sum5 > 0d ? s5 / sum5 : double.NaN];
        }

        /// <summary>
        /// Returns the discordancy measure D_i of each site.
        /// </summary>
        /// <returns>D_i by site; every D_i is 1 for three or fewer sites, and every D_i is NaN when the sum of
        /// squares and cross-products matrix A is singular.</returns>
        /// <remarks>
        /// <para>
        /// Hosking and Wallis (1997), section 3.2.3: with u_i = (t^(i), t3^(i), t4^(i)), their unweighted mean ū and
        /// A = Σ (u_i − ū)(u_i − ū)ᵀ, D_i = N/3 (u_i − ū)ᵀ A⁻¹ (u_i − ū). The mean of the D_i is 1, and no D_i
        /// exceeds (N − 1)/3.
        /// </para>
        /// <para>
        /// A is inverted by the three symmetric sweeps of Hosking's REGTST routine, which fail when a pivot is not
        /// positive (the sites' ratios lie in a plane). A singular A then gives NaN for every D_i. This differs
        /// deliberately from REGTST, which leaves D_i = 0 for every site, so that a failed inversion cannot be read as
        /// a region without discordant sites; lmomRFA's regtst also reports NA, with a warning.
        /// </para>
        /// </remarks>
        public double[] Discordancy()
        {
            int sites = SiteCount;
            var d = new double[sites];
            if (sites <= 3)
            {
                for (int i = 0; i < sites; i++) d[i] = 1d;
                return d;
            }
            double m2 = 0d, m3 = 0d, m4 = 0d;
            for (int i = 0; i < sites; i++)
            {
                m2 += _lcv[i];
                m3 += _t3[i];
                m4 += _t4[i];
            }
            m2 /= sites;
            m3 /= sites;
            m4 /= sites;
            var work = new double[sites, 3];
            for (int i = 0; i < sites; i++)
            {
                work[i, 0] = _lcv[i] - m2;
                work[i, 1] = _t3[i] - m3;
                work[i, 2] = _t4[i] - m4;
            }
            var s = new double[3, 3];
            for (int j = 0; j < 3; j++)
                for (int k = j; k < 3; k++)
                    for (int i = 0; i < sites; i++)
                        s[j, k] += work[i, j] * work[i, k];

            // Invert the symmetric 3 x 3 matrix in place by three sweeps of its upper triangle.
            for (int k = 0; k < 3; k++)
            {
                if (!(s[0, 0] > 0d))
                {
                    for (int i = 0; i < sites; i++) d[i] = double.NaN;
                    return d;
                }
                double temp0 = 1d / s[0, 0];
                double temp1 = -s[0, 1] * temp0;
                double temp2 = -s[0, 2] * temp0;
                if (k > 1) temp1 = -temp1;
                if (k > 0) temp2 = -temp2;
                s[0, 0] = s[1, 1] + temp1 * s[0, 1];
                s[0, 1] = s[1, 2] + temp1 * s[0, 2];
                s[1, 1] = s[2, 2] + temp2 * s[0, 2];
                s[0, 2] = temp1;
                s[1, 2] = temp2;
                s[2, 2] = temp0;
            }
            s[1, 0] = s[0, 1];
            s[2, 0] = s[0, 2];
            s[2, 1] = s[1, 2];

            double factor = sites / 3d;
            for (int i = 0; i < sites; i++)
            {
                double q = 0d;
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        q += work[i, j] * work[i, k] * s[j, k];
                d[i] = q * factor;
            }
            return d;
        }

        /// <summary>
        /// Returns the critical value of the discordancy measure for a region.
        /// </summary>
        /// <param name="siteCount">The number of sites, at least 1.</param>
        /// <param name="significance">Optional. The significance level, 0.10 (the default) or 0.05.</param>
        /// <returns>The value at or above which a site is declared discordant.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="siteCount"/> is below 1, or
        /// <paramref name="significance"/> is neither 0.10 nor 0.05.</exception>
        /// <remarks>
        /// The 10% values are those of Hosking and Wallis (1997), Table 3.1, rising from 1.333 at 5 sites to 3 at 15
        /// sites and more; the 5% values, as in lmomRFA, rise to 4 from 19 sites. For four or fewer sites, which
        /// cannot produce a discordant site (D_i is at most (N − 1)/3), the values are 3 and 4, as in lmomRFA.
        /// </remarks>
        public static double DiscordancyCriticalValue(int siteCount, double significance = 0.10)
        {
            if (siteCount < 1) throw new ArgumentOutOfRangeException(nameof(siteCount), "The number of sites must be at least 1.");
            double[] table;
            double cap;
            if (significance == 0.10)
            {
                table = Critical10;
                cap = 3d;
            }
            else if (significance == 0.05)
            {
                table = Critical5;
                cap = 4d;
            }
            else
            {
                throw new ArgumentOutOfRangeException(nameof(significance), "The significance level must be 0.10 or 0.05.");
            }
            return siteCount > table.Length ? cap : table[siteCount - 1];
        }
    }
}
