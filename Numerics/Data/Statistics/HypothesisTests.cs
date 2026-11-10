using Numerics.Distributions;
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.SpecialFunctions;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.Data.Statistics
{
    /// <summary>
    /// A class for performing statistical hypothesis tests. 
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// <para>
    /// <b> Description: </b>
    /// This is a method of statistical inference used to decide whether the data sufficiently support a particular hypothesis, typically
    /// through the calculation of some test statistic. The null hypothesis (H0) is assumed to be true until proven otherwise through the tests. In this class
    /// the p value is calculated to determine whether or not to reject the null hypothesis. If the p value is less than the specified significance level,
    /// the null hypothesis is rejected in favor of the alternative hypothesis (Ha).
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Statistical_hypothesis_test"/>
    /// </remarks>
    public class HypothesisTests
    {
        /// <summary>
        /// The one sample t-Test compares the mean of the sample to a hypothesized population mean.
        /// One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <param name="populationMean">Optional. The hypothesized mean. Default = 0.</param>
        /// <returns>Returns the 2-sided p-value of the test statistic.</returns>
        public static double OneSampleTtest(IList<double> sample, double populationMean = 0d)
        {
            if (sample.Count < 2) throw new ArgumentException("Sample must have at least 2 observations.", nameof(sample));
            var meanVar = Statistics.MeanVariance(sample);
            int N = sample.Count;
            double se = Math.Sqrt(meanVar.Item2) / Math.Sqrt(N);
            double t = Math.Abs((meanVar.Item1 - populationMean) / se);
            return 2d * new StudentT(N - 1).CCDF(t);
        }

        /// <summary>
        /// The t-test determines if there is a significant difference between the means of two samples drawn from populations with equal variances.
        /// One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample1">Data sample 1.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the 2-sided p-value of the test statistic.</returns>
        public static double EqualVarianceTtest(IList<double> sample1, IList<double> sample2)
        {
            if (sample1.Count + sample2.Count < 3) throw new ArgumentException("Combined sample size must be at least 3.");
            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            int N1 = sample1.Count;
            int N2 = sample2.Count;
            double df = N1 + N2 - 2;
            // A single-value sample contributes no squared deviations to the pooled sum; its sample variance is NaN.
            double ss1 = N1 == 1 ? 0d : (N1 - 1) * meanVar1.Item2;
            double ss2 = N2 == 1 ? 0d : (N2 - 1) * meanVar2.Item2;
            double svar = (ss1 + ss2) / df;
            double t = Math.Abs((meanVar1.Item1 - meanVar2.Item1) / Math.Sqrt(svar * (1.0d / N1 + 1.0d / N2)));
            return 2d * new StudentT(df).CCDF(t);
        }

        /// <summary>
        /// The t-test determines if there is a significant difference between the means of two samples drawn from populations with unequal variances.
        /// One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample1">Data sample 1.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the 2-sided p-value of the test statistic.</returns>
        public static double UnequalVarianceTtest(IList<double> sample1, IList<double> sample2)
        {
            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            double ave1 = meanVar1.Item1, var1 = meanVar1.Item2, ave2 = meanVar2.Item1, var2 = meanVar2.Item2;
            int n1 = sample1.Count;
            int n2 = sample2.Count;
            double t = Math.Abs((ave1 - ave2) / Math.Sqrt(var1 / n1 + var2 / n2));
            double df = Tools.Sqr(var1 / n1 + var2 / n2) / (Tools.Sqr(var1 / n1) / (n1 - 1) + Tools.Sqr(var2 / n2) / (n2 - 1));
            return 2d * new StudentT(df).CCDF(t);
        }

        /// <summary>
        /// The paired t-test determines whether the mean difference between two sets of observations is zero.
        /// One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample1">Data sample 1.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the 2-sided p-value of the test statistic.</returns>
        public static double PairedTtest(IList<double> sample1, IList<double> sample2)
        {
            if (sample1.Count != sample2.Count) throw new ArgumentException("The two data samples must be the sample length.");

            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            int n = sample1.Count;
            double ave1 = meanVar1.Item1, var1 = meanVar1.Item2, ave2 = meanVar2.Item1, var2 = meanVar2.Item2, sd, df, cov = 0.0;
            for (int j = 0; j < n; j++) cov += (sample1[j] - ave1) * (sample2[j] - ave2);
            cov /= (df = n - 1);
            sd = Math.Sqrt((var1 + var2 - 2.0 * cov) / n);
            double t = Math.Abs((ave1 - ave2) / sd);
            return 2d * new StudentT(df).CCDF(t);
        }

        /// <summary>
        /// The F-test for significantly different variances. One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample1">Data sample 1.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double Ftest(IList<double> sample1, IList<double> sample2)
        {
            if (sample1.Count < 2 || sample2.Count < 2) throw new ArgumentException("Each sample must have at least 2 observations.");
            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            int n1 = sample1.Count, n2 = sample2.Count;
            double ave1 = meanVar1.Item1, var1 = meanVar1.Item2, ave2 = meanVar2.Item1, var2 = meanVar2.Item2, df1, df2, f, pVal;
            if (var1 == 0 && var2 == 0) return 1.0;
            if (var1 > var2)
            {
                f = var1 / var2;
                df1 = n1 - 1;
                df2 = n2 - 1;
            }
            else
            {
                f = var2 / var1;
                df1 = n2 - 1;
                df2 = n1 - 1;
            }
            pVal = 2.0 * Beta.Incomplete(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f));
            if (pVal > 1.0) pVal = 2d - pVal;
            return pVal;
        }

        /// <summary>
        /// The F-test comparing two models. The null, H0, states that the restricted and full models are equal. 
        /// One rejects the null if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sseRestricted">The sum of squared errors (SSE) of the restricted model.</param>
        /// <param name="sseFull">The sum of squared errors (SSE) of the full model.</param>
        /// <param name="dfRestricted">The residual degrees of freedom of the restricted model.</param>
        /// <param name="dfFull">The residual degrees of freedom of the full model, fewer than the restricted model's.</param>
        /// <param name="fStat">The F-test statistic.</param>
        /// <param name="pValue">The p-value of the test statistic.</param>
        /// <exception cref="ArgumentException">Thrown when <paramref name="dfFull"/> is not positive or is not less
        /// than <paramref name="dfRestricted"/>.</exception>
        /// <remarks>
        /// F = ((SSE_R − SSE_F) / (df_R − df_F)) / (SSE_F / df_F), and the p-value is its one-sided upper tail
        /// P(F(df_R − df_F, df_F) &gt; F) = I_x(df_F / 2, (df_R − df_F) / 2) with x = df_F / (df_F + (df_R − df_F) F),
        /// as in R's anova for nested linear models. A statistic of zero or less (no improvement) gives 1.
        /// </remarks>
        public static void FtestModels(double sseRestricted, double sseFull, int dfRestricted, int dfFull, out double fStat, out double pValue)
        {
            if (dfRestricted == dfFull) throw new ArgumentException("Restricted and full model degrees of freedom cannot be equal.");
            if (dfFull <= 0) throw new ArgumentException("Full model degrees of freedom must be positive.", nameof(dfFull));
            if (dfRestricted < dfFull) throw new ArgumentException("The restricted model must have more residual degrees of freedom than the full model.", nameof(dfRestricted));
            int dfNumerator = dfRestricted - dfFull;
            fStat = ((sseRestricted - sseFull) / dfNumerator) / (sseFull / dfFull);
            if (double.IsNaN(fStat))
                pValue = double.NaN;
            else if (fStat <= 0d)
                pValue = 1d;
            else
                pValue = Beta.Incomplete(0.5 * dfFull, 0.5 * dfNumerator, dfFull / (dfFull + dfNumerator * fStat));
        }

        /// <summary>
        /// The Jarque-Bera test for normality. One rejects H0 of normality if the p-value of the JB statistic is less than the significance level.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        /// <remarks>
        /// JB = (n / 6) (S² + K² / 4), with S = m3 / m2^(3/2) and K = m4 / m2² − 3 the moment (not bias-corrected)
        /// skewness and excess kurtosis, m_k = (1/n) Σ (x − x̄)^k. Under normality JB is asymptotically chi-squared with
        /// 2 degrees of freedom, as in Jarque and Bera (1987), R's moments::jarque.test and scipy.stats.jarque_bera.
        /// </remarks>
        public static double JarqueBeraTest(IList<double> sample)
        {
            int N = sample.Count;
            double mean = Statistics.Mean(sample);
            double m2 = 0, m3 = 0, m4 = 0;
            for (int i = 0; i < N; i++)
            {
                double d = sample[i] - mean;
                double d2 = d * d;
                m2 += d2;
                m3 += d2 * d;
                m4 += d2 * d2;
            }
            m2 /= N;
            m3 /= N;
            m4 /= N;
            double S = m3 / Math.Pow(m2, 1.5);
            double K = m4 / (m2 * m2) - 3d;
            double JB = N / 6.0d * (S * S + K * K / 4.0d);
            var chi = new ChiSquared(2);
            return 1d - chi.CDF(JB);
        }

        /// <summary>
        /// The Wald and Wolfowitz test for independence and stationarity (trend). 
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <returns>Returns the 2-sided p-value of the test statistic.</returns>
        public static double WaldWolfowitzTest(IList<double> sample)
        {
            double xn = sample.Last(), x1 = sample.First();
            double R = 0, eR, varR, U, S1 = 0, S2 = 0, S3 = 0, S4 = 0;
            int n = sample.Count;
            for (int i = 0; i < n; i++)
            {

                S1 += sample[i];
                S2 += Tools.Sqr(sample[i]);
                S3 += Tools.Pow(sample[i], 3);
                S4 += Tools.Pow(sample[i], 4);

                if (i < n - 1)
                {
                    R += sample[i] * sample[i + 1];
                }
            }
            R += xn * x1;
            double s12 = Tools.Sqr(S1);
            double s22 = Tools.Sqr(S2);
            double s14 = Tools.Pow(S1, 4);

            eR = (s12 - S2) / (n - 1);
            varR = (s22 - S4) / (n - 1) - Tools.Sqr(eR);
            varR += (s14 - 4 * s12 * S2 + 4 * S1 * S3 + s22 - 2 * S4) / ((n - 1) * (n - 2));
            U = Math.Abs((R - eR) / Math.Sqrt(varR));
            return (1 - Normal.StandardCDF(U)) * 2d;
        }

        /// <summary>
        /// The Ljung-Box test whether the autocorrelations of the data are different from zero.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <param name="lagMax">The max lag to evaluate.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double LjungBoxTest(IList<double> sample, int lagMax = -1)
        {
            int n = sample.Count;
            if (lagMax < 0) lagMax = (int)Math.Floor(Math.Min(10d * Math.Log10(n), n - 1));
            var acf = Autocorrelation.Function(sample, lagMax, Autocorrelation.Type.Correlation);
            if (acf == null) return double.NaN;
            double Q = 0;
            for (int k = 1; k <= lagMax; k++)
                Q += Tools.Sqr(acf[k, 1]) / (n - k);
            Q *= n * (n + 2);
            var chi2 = new ChiSquared(lagMax);
            return 1d - chi2.CDF(Q);
        }

        /// <summary>
        /// Returns the size of each group of tied values in the sample.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <returns>One entry per group of two or more equal values, holding the group size.</returns>
        private static List<double> TieGroupSizes(IList<double> sample)
        {
            var sorted = sample.ToArray();
            Array.Sort(sorted);
            var groups = new List<double>();
            int run = 1;
            for (int i = 1; i < sorted.Length; i++)
            {
                if (sorted[i] == sorted[i - 1])
                {
                    run++;
                }
                else
                {
                    if (run > 1) groups.Add(run);
                    run = 1;
                }
            }
            if (run > 1) groups.Add(run);
            return groups;
        }

        /// <summary>
        /// The Mann-Whitney test for homogeneity and stationarity (jump). 
        /// </summary>
        /// <param name="sample1">Data sample 1. Must be less than or equal in length to sample 2.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double MannWhitneyTest(IList<double> sample1, IList<double> sample2)
        {
            int n1 = sample1.Count, n2 = sample2.Count, n = n1 + n2;
            if (n1 > n2) throw new ArgumentException("The first sample must have a length less than or equal to sample 2.");
            if (n1 <= 3 || n2 <= 3) throw new ArgumentException("Each sample must have a length greater than 3.");
            if (n <= 20) throw new ArgumentException("The combined sample 1 & 2 must have a length greater than 20.");

            var sample = new List<double>();
            sample.AddRange(sample1.ToList());
            sample.AddRange(sample2.ToList());
            double R = 0, T = 0;

            var ranks = Statistics.RanksInPlace(sample.ToArray(), out _);

            for (int i = 0; i < sample1.Count; i++) R += ranks[i];
            // The variance correction sums over tie groups by their full size
            foreach (double g in TieGroupSizes(sample))
                T += (Tools.Pow(g, 3) - g) / (n * (n - 1));

            double V = R - n1 * (n1 + 1d) / 2d;
            double W = n1 * n2 - V;
            double U = Math.Min(V, W);
            double eU = n1 * n2 / 2d;
            double varU = n1 * n2 / 12d * ((n + 1) - T);
            double z = Math.Abs((U - eU) / Math.Sqrt(varU));
            return (1 - Normal.StandardCDF(z)) * 2d;
        }

        /// <summary>
        /// The Mann-Kendall test for homogeneity and stationarity (trend).
        /// </summary>
        /// <param name="sample">Data sample.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double MannKendallTest(IList<double> sample)
        {
            int i, j, n = sample.Count;
            if (n < 10) throw new ArgumentException("The sample size must be greater than or equal to 10.");

            double S = 0, T = 0, varS, z;

            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    S += Math.Sign(sample[j] - sample[i]);
                }
            }

            // The variance correction sums over tie groups by their full size
            foreach (double g in TieGroupSizes(sample))
                T += g * (g - 1) * (2 * g + 5);
            varS = (n * (n - 1) * (2 * n + 5) - T) / 18;
            z = Math.Abs((S - Math.Sign(S)) / Math.Sqrt(varS));

            return (1 - Normal.StandardCDF(z)) * 2d;
        }

        /// <summary>
        /// The linear trend test for stationarity (trend). 
        /// </summary>
        /// <param name="indices">Time series indices.</param>
        /// <param name="sample">Time series sample data.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double LinearTrendTest(IList<double> indices, IList<double> sample)
        {
            if (indices.Count != sample.Count)
                throw new ArgumentException("The indices must be the same length as the sample data");
            var xVals = new Matrix(indices.ToArray());
            var yVals = new Vector(sample.ToArray());
            var lm = new LinearRegression(xVals, yVals, true);
            return 2d * new StudentT(lm.DegreesOfFreedom).CCDF(Math.Abs(lm.Parameters[1] / lm.ParameterStandardErrors[1]));
        }

        /// <summary>
        /// Computes the exact two-sample Kolmogorov-Smirnov statistic between two samples sorted in
        /// ascending order: the supremum absolute difference of their empirical distribution functions.
        /// </summary>
        /// <param name="first">The first sample, sorted ascending. Must be nonempty and finite.</param>
        /// <param name="second">The second sample, sorted ascending. Must be nonempty and finite.</param>
        /// <returns>The statistic sup |F1(t) - F2(t)|, in [0, 1].</returns>
        /// <remarks>
        /// <para>
        /// The samples are swept once in merged order. At each distinct value both empirical
        /// distributions advance through every tied copy of that value before the gap is compared,
        /// so ties within or across the samples are counted as single jumps and the statistic is
        /// exact for discrete as well as continuous data. The cost is O(n1 + n2) after the
        /// O(n1 + n2) argument check; an array argument is not copied.
        /// </para>
        /// <para>
        /// The gap is evaluated in integer arithmetic: with i and j the counts of first and second
        /// values at or below the current value, the statistic is D / (n1 n2), where
        /// D = max |i n2 - j n1| is accumulated in 64-bit integers and divided once in double
        /// precision. The products are at most n1 n2 &lt; 2^62 for any array lengths, so they cannot
        /// overflow. When n1 n2 &lt;= 2^53 both D and n1 n2 are exact doubles and the result is the
        /// correctly rounded value of the rational statistic; beyond that D and n1 n2 are each rounded
        /// once and the result is within two ulps of it.
        /// </para>
        /// <para>
        /// The inputs must already be sorted ascending; the method checks the order and throws
        /// rather than sorting, so a caller that already holds sorted samples pays no sort. This is
        /// the statistic <see cref="GlobalSensitivity.Pawn"/> computes for each bin, evaluated by the same sweep.
        /// No p-value is returned. The result equals scipy.stats.ks_2samp(first, second).statistic where
        /// scipy reports the statistic rounded from its integer count (its exact mode in scipy 1.17); a
        /// statistic taken as the floating-point difference i/n1 - j/n2 can differ by an ulp or so.
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// Conover, W. J. (1999). Practical Nonparametric Statistics, 3rd ed., Chapter 6. Wiley.
        /// </para>
        /// </remarks>
        /// <exception cref="ArgumentNullException">Thrown when either sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when a sample is empty or is not sorted ascending.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a sample value is not finite.</exception>
        public static double TwoSampleKolmogorovSmirnov(IList<double> first, IList<double> second)
        {
            ValidateSortedSample(first, nameof(first));
            ValidateSortedSample(second, nameof(second));
            return TwoSampleKolmogorovSmirnovSorted(first as double[] ?? first.ToArray(), second as double[] ?? second.ToArray());
        }

        /// <summary>
        /// Checks that a sample for <see cref="TwoSampleKolmogorovSmirnov"/> is non-null, nonempty,
        /// finite, and sorted ascending.
        /// </summary>
        /// <param name="sample">The sample to check.</param>
        /// <param name="name">The argument name reported by the exceptions.</param>
        private static void ValidateSortedSample(IList<double> sample, string name)
        {
            if (sample == null) throw new ArgumentNullException(name);
            if (sample.Count == 0) throw new ArgumentException("The sample must not be empty.", name);
            for (int i = 0; i < sample.Count; i++)
            {
                if (!Tools.IsFinite(sample[i]))
                    throw new ArgumentOutOfRangeException(name, "The sample values must be finite.");
                if (i > 0 && sample[i] < sample[i - 1])
                    throw new ArgumentException("The sample must be sorted in ascending order.", name);
            }
        }

        /// <summary>
        /// Computes the exact two-sample Kolmogorov-Smirnov statistic between two sorted samples,
        /// evaluating the gap only after advancing both empirical distributions through tied values.
        /// </summary>
        /// <param name="first">The first sorted sample.</param>
        /// <param name="second">The second sorted sample.</param>
        /// <returns>The supremum absolute difference of the empirical distributions.</returns>
        /// <remarks>
        /// The gap |i/n1 - j/n2| is tracked as the integer |i n2 - j n1|, which is exact in 64-bit
        /// arithmetic for any array lengths, and the maximum is divided by n1 n2 once, so the result
        /// is the correctly rounded rational statistic whenever n1 n2 &lt;= 2^53.
        /// </remarks>
        internal static double TwoSampleKolmogorovSmirnovSorted(double[] first, double[] second)
        {
            int n1 = first.Length, n2 = second.Length;
            int i = 0, j = 0;
            long supremum = 0;
            while (i < n1 && j < n2)
            {
                double value = Math.Min(first[i], second[j]);
                while (i < n1 && first[i] == value) i++;
                while (j < n2 && second[j] == value) j++;
                long gap = Math.Abs((long)i * n2 - (long)j * n1);
                if (gap > supremum) supremum = gap;
            }
            // The remaining tail of either sample only shrinks toward equality at one, and the gap
            // at the crossover was already recorded, so the sweep is complete.
            return supremum / (double)((long)n1 * n2);
        }

        /// <summary>
        /// Computes Hartigan's dip statistic, the largest distance between the empirical distribution
        /// function of the sample and the closest unimodal distribution function.
        /// </summary>
        /// <param name="sample">The data sample. The order of the values does not matter.</param>
        /// <returns>The dip statistic, between 1/(2n) and 1/4 for n &gt;= 1, and 0 for an empty sample.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the sample contains NaN or infinite values.</exception>
        /// <remarks>
        /// <para>
        /// The dip is computed with the algorithm of Hartigan (1985), Applied Statistics algorithm AS 217,
        /// which alternates the greatest convex minorant and least concave majorant of the empirical
        /// distribution function over a shrinking modal interval. Values agree with R's
        /// <c>diptest::dip</c>. Tied values are allowed; a constant sample has the minimum dip 1/(2n).
        /// </para>
        /// <para>
        /// <b> References: </b>
        /// </para>
        /// <para>
        /// Hartigan, J. A. and Hartigan, P. M. (1985). The dip test of unimodality. The Annals of
        /// Statistics, 13(1), 70-84.
        /// </para>
        /// <para>
        /// Hartigan, P. M. (1985). Algorithm AS 217: Computation of the dip statistic to test for
        /// unimodality. Journal of the Royal Statistical Society, Series C (Applied Statistics), 34, 320-325.
        /// </para>
        /// </remarks>
        public static double HartiganDip(IList<double> sample)
        {
            if (sample == null) throw new ArgumentNullException(nameof(sample));
            var x = new double[sample.Count];
            for (int i = 0; i < x.Length; i++)
            {
                x[i] = sample[i];
                if (!Tools.IsFinite(x[i]))
                    throw new ArgumentException("The sample must not contain NaN or infinite values.", nameof(sample));
            }
            Array.Sort(x);
            return SortedDip(x);
        }

        /// <summary>
        /// Hartigan's dip test of unimodality. One rejects the null hypothesis of a unimodal distribution
        /// if the p-value is less than the significance level.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <param name="realizations">Optional. The number of simulated uniform samples. Default = 2,000.</param>
        /// <param name="seed">Optional. The seed of the simulation, which makes the p-value reproducible. Default = 12345.</param>
        /// <returns>Returns the p-value of the dip statistic. Samples of three or fewer values return 1.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the sample is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the sample contains NaN or infinite values.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the number of realizations is less than one.</exception>
        /// <remarks>
        /// The uniform distribution is the least favorable unimodal null, so the p-value is the fraction
        /// of <paramref name="realizations"/> uniform samples of the same size whose dip is at least the
        /// observed dip, as in R's <c>diptest::dip.test(x, simulate.p.value = TRUE)</c>. Its Monte Carlo
        /// standard error is about sqrt(p(1 - p) / realizations). The test is sensitive to any departure
        /// from a single mode, including multimodality in skewed or heavy-tailed data, and makes no
        /// assumption about the shape of the unimodal distribution. See <see cref="HartiganDip"/>.
        /// </remarks>
        public static double HartiganDipTest(IList<double> sample, int realizations = 2000, int seed = 12345)
        {
            if (realizations < 1) throw new ArgumentOutOfRangeException(nameof(realizations), "The number of realizations must be at least one.");
            double dip = HartiganDip(sample);
            int n = sample.Count;
            if (n <= 3) return 1d;

            var seeds = new MersenneTwister(seed).NextIntegers(1, int.MaxValue, realizations);
            var simulated = new double[realizations];
            Parallel.For(0, realizations, b =>
            {
                var u = new MersenneTwister(seeds[b]).NextDoubles(n);
                Array.Sort(u);
                simulated[b] = SortedDip(u);
            });

            int count = 0;
            for (int b = 0; b < realizations; b++)
                if (dip <= simulated[b]) count++;
            return (double)count / realizations;
        }

        /// <summary>
        /// Computes the dip statistic of sorted, finite data with algorithm AS 217.
        /// </summary>
        /// <param name="sorted">The data in ascending order.</param>
        /// <returns>The dip statistic.</returns>
        /// <remarks>
        /// Indices are 1-based to follow the published algorithm. The modal interval [low, high] starts
        /// at the full sample; each pass builds the greatest convex minorant (from high down to low) and
        /// the least concave majorant (from low up to high), finds the largest vertical distance between
        /// them, takes the dips of the minorant and majorant outside the new interval, and stops when the
        /// interval no longer shrinks or its distance falls below the current dip. Distances are counted
        /// in observations and scaled by 1/(2n) at the end.
        /// </remarks>
        private static double SortedDip(double[] sorted)
        {
            int n = sorted.Length;
            if (n == 0) return 0d;
            var x = new double[n + 1];
            for (int i = 0; i < n; i++)
                x[i + 1] = sorted[i];

            double dip = 1d;
            if (n < 2 || x[n] == x[1])
                return dip / (2d * n);

            // Convex minorant pointers: mn[j] is the predecessor of j on the minorant of x[1..j]
            var mn = new int[n + 1];
            mn[1] = 1;
            for (int j = 2; j <= n; j++)
            {
                mn[j] = j - 1;
                while (true)
                {
                    int mnj = mn[j];
                    int mnmnj = mn[mnj];
                    if (mnj == 1 || (x[j] - x[mnj]) * (mnj - mnmnj) < (x[mnj] - x[mnmnj]) * (j - mnj))
                        break;
                    mn[j] = mnmnj;
                }
            }

            // Concave majorant pointers: mj[k] is the successor of k on the majorant of x[k..n]
            var mj = new int[n + 1];
            mj[n] = n;
            for (int k = n - 1; k >= 1; k--)
            {
                mj[k] = k + 1;
                while (true)
                {
                    int mjk = mj[k];
                    int mjmjk = mj[mjk];
                    if (mjk == n || (x[k] - x[mjk]) * (mjk - mjmjk) < (x[mjk] - x[mjmjk]) * (k - mjk))
                        break;
                    mj[k] = mjmjk;
                }
            }

            var gcm = new int[n + 1];
            var lcm = new int[n + 1];
            int low = 1, high = n;
            while (true)
            {
                // Change points of the convex minorant from high down to low
                int ic = 1;
                gcm[1] = high;
                while (gcm[ic] > low)
                {
                    gcm[ic + 1] = mn[gcm[ic]];
                    ic++;
                }
                int lGcm = ic, ig = lGcm, ix = lGcm - 1;

                // Change points of the concave majorant from low up to high
                ic = 1;
                lcm[1] = low;
                while (lcm[ic] < high)
                {
                    lcm[ic + 1] = mj[lcm[ic]];
                    ic++;
                }
                int lLcm = ic, ih = lLcm, iv = 2;

                // Largest distance between the minorant and the majorant over [low, high]
                double d = 0d;
                if (lGcm != 2 || lLcm != 2)
                {
                    do
                    {
                        int gcmix = gcm[ix];
                        int lcmiv = lcm[iv];
                        if (gcmix > lcmiv)
                        {
                            int gcmi1 = gcm[ix + 1];
                            double dx = (lcmiv - gcmi1 + 1) - (x[lcmiv] - x[gcmi1]) * (gcmix - gcmi1) / (x[gcmix] - x[gcmi1]);
                            iv++;
                            if (dx >= d)
                            {
                                d = dx;
                                ig = ix + 1;
                                ih = iv - 1;
                            }
                        }
                        else
                        {
                            int lcmiv1 = lcm[iv - 1];
                            double dx = (x[gcmix] - x[lcmiv1]) * (lcmiv - lcmiv1) / (x[lcmiv] - x[lcmiv1]) - (gcmix - lcmiv1 - 1);
                            ix--;
                            if (dx >= d)
                            {
                                d = dx;
                                ig = ix + 1;
                                ih = iv;
                            }
                        }
                        if (ix < 1) ix = 1;
                        if (iv > lLcm) iv = lLcm;
                    } while (gcm[ix] != lcm[iv]);
                }
                else
                {
                    d = 1d;
                }

                if (d < dip)
                    break;

                // Dip of the convex minorant between its change points ig..lGcm
                double dipL = 0d;
                for (int j = ig; j < lGcm; j++)
                {
                    double maxT = 1d;
                    int jb = gcm[j + 1], je = gcm[j];
                    if (je - jb > 1 && x[je] != x[jb])
                    {
                        double c = (je - jb) / (x[je] - x[jb]);
                        for (int jj = jb; jj <= je; jj++)
                        {
                            double t = (jj - jb + 1) - (x[jj] - x[jb]) * c;
                            if (maxT < t) maxT = t;
                        }
                    }
                    if (dipL < maxT) dipL = maxT;
                }

                // Dip of the concave majorant between its change points ih..lLcm
                double dipU = 0d;
                for (int j = ih; j < lLcm; j++)
                {
                    double maxT = 1d;
                    int jb = lcm[j], je = lcm[j + 1];
                    if (je - jb > 1 && x[je] != x[jb])
                    {
                        double c = (je - jb) / (x[je] - x[jb]);
                        for (int jj = jb; jj <= je; jj++)
                        {
                            double t = (x[jj] - x[jb]) * c - (jj - jb - 1);
                            if (maxT < t) maxT = t;
                        }
                    }
                    if (dipU < maxT) dipU = maxT;
                }

                double dipNew = Math.Max(dipL, dipU);
                if (dip < dipNew) dip = dipNew;

                // Stop when the modal interval no longer shrinks
                if (low == gcm[ig] && high == lcm[ih])
                    break;
                low = gcm[ig];
                high = lcm[ih];
            }

            return dip / (2d * n);
        }

        /// <summary>
        /// The Gaussian mixture likelihood-ratio test of a single Gaussian against a two-component
        /// Gaussian mixture, with an approximate chi-square p-value.
        /// </summary>
        /// <param name="sample">The data sample, with at least 10 values.</param>
        /// <returns>Returns the p-value of the test statistic, or NaN when a mixture fit fails.</returns>
        /// <exception cref="ArgumentException">Thrown when the sample has fewer than 10 values.</exception>
        /// <remarks>
        /// <para>
        /// The statistic is 2(log L2 - log L1), where L1 and L2 are the maximized likelihoods of one- and
        /// two-component Gaussian mixtures fitted by expectation-maximization with k-Means++ starts (seed
        /// 12345). The p-value refers it to a chi-square distribution with 3 degrees of freedom, the
        /// number of extra parameters (a mean, a variance and a weight).
        /// </para>
        /// <para>
        /// The test asks whether a sample is better described by one normal component or by a mixture of
        /// two, which makes it a useful check for mixed populations, such as floods from two generating
        /// mechanisms. The chi-square reference is an approximation: under a single Gaussian the second
        /// component is not identifiable, so the usual likelihood-ratio theory does not hold exactly. In
        /// simulation (200 samples per case, nominal 5% level) normal samples were rejected about 17% of the
        /// time at n = 30 and 11% at n = 100. Skewed or heavy-tailed unimodal samples
        /// (lognormal, gamma, Student t) are often rejected because a two-component mixture fits their
        /// shape better than one Gaussian, so a small p-value indicates "not a single Gaussian" rather
        /// than "more than one mode". Use <see cref="UnimodalityTest(IList{double}, int, int)"/> for a
        /// calibrated p-value, and <see cref="HartiganDipTest"/> to test for a single mode without
        /// assuming normality.
        /// </para>
        /// </remarks>
        public static double UnimodalityTest(IList<double> sample)
        {
            int n = sample.Count;
            if (n < 10) throw new ArgumentException("The sample size must be greater than or equal to 10.");

            try
            {
                double testStat = MixtureLikelihoodRatio(sample.ToArray(), out _, out _);
                int df = 3;

                // Compute p-value from chi-square distribution
                var chiSquared = new ChiSquared(df);
                double pval = 1.0 - chiSquared.CDF(testStat);

                return pval;
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"Error in hypothesis testing: {ex.Message}");
                return double.NaN;
            }

        }

        /// <summary>
        /// The Gaussian mixture likelihood-ratio test of a single Gaussian against a two-component
        /// Gaussian mixture, with a parametric bootstrap p-value.
        /// </summary>
        /// <param name="sample">The data sample, with at least 10 values.</param>
        /// <param name="realizations">The number of bootstrap samples, for example 1,000.</param>
        /// <param name="seed">Optional. The seed of the bootstrap, which makes the p-value reproducible. Default = 12345.</param>
        /// <returns>Returns the bootstrap p-value of the test statistic, or NaN when the mixture fit of the sample fails.</returns>
        /// <exception cref="ArgumentException">Thrown when the sample has fewer than 10 values.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the number of realizations is less than one.</exception>
        /// <remarks>
        /// The statistic is the same as in <see cref="UnimodalityTest(IList{double})"/>. Its null
        /// distribution is simulated by drawing <paramref name="realizations"/> samples of the same size from
        /// the fitted single Gaussian and refitting both models to each, so the p-value
        /// (1 + #{T_b &gt;= T}) / (1 + B) has the nominal size for normal data at any sample size
        /// (McLachlan, 1987). Bootstrap samples whose mixture fit fails are left out of the count. The cost
        /// is about 2B mixture fits.
        /// <para>
        /// McLachlan, G. J. (1987). On bootstrapping the likelihood ratio test statistic for the number of
        /// components in a normal mixture. Journal of the Royal Statistical Society, Series C (Applied
        /// Statistics), 36(3), 318-324.
        /// </para>
        /// </remarks>
        public static double UnimodalityTest(IList<double> sample, int realizations, int seed = 12345)
        {
            int n = sample.Count;
            if (n < 10) throw new ArgumentException("The sample size must be greater than or equal to 10.");
            if (realizations < 1) throw new ArgumentOutOfRangeException(nameof(realizations), "The number of realizations must be at least one.");

            double observed, mean, standardDeviation;
            try
            {
                observed = MixtureLikelihoodRatio(sample.ToArray(), out mean, out standardDeviation);
            }
            catch (Exception ex)
            {
                Debug.WriteLine($"Error in hypothesis testing: {ex.Message}");
                return double.NaN;
            }

            var seeds = new MersenneTwister(seed).NextIntegers(1, int.MaxValue, realizations);
            var simulated = new double[realizations];
            Parallel.For(0, realizations, b =>
            {
                var replicate = new Normal(mean, standardDeviation).GenerateRandomValues(n, seeds[b]);
                try
                {
                    simulated[b] = MixtureLikelihoodRatio(replicate, out _, out _);
                }
                catch (Exception ex)
                {
                    // A replicate whose mixture fit fails is excluded from the bootstrap count.
                    Debug.WriteLine($"Bootstrap replicate {b} failed: {ex.Message}");
                    simulated[b] = double.NaN;
                }
            });

            int valid = 0, exceed = 0;
            for (int b = 0; b < realizations; b++)
            {
                if (double.IsNaN(simulated[b])) continue;
                valid++;
                if (simulated[b] >= observed) exceed++;
            }
            return (exceed + 1d) / (valid + 1d);
        }

        /// <summary>
        /// Returns the likelihood-ratio statistic 2(log L2 - log L1) of two- against one-component
        /// Gaussian mixtures, fitted with k-Means++ starts and seed 12345.
        /// </summary>
        /// <param name="data">The data sample.</param>
        /// <param name="mean">The mean of the fitted single Gaussian.</param>
        /// <param name="standardDeviation">The standard deviation of the fitted single Gaussian.</param>
        /// <returns>The likelihood-ratio statistic.</returns>
        private static double MixtureLikelihoodRatio(double[] data, out double mean, out double standardDeviation)
        {
            // Fit 1-component GMM (unimodal)
            var gmm1 = new GaussianMixtureModel(data, 1);
            gmm1.Train(12345, true);
            var logLH1 = gmm1.LogLikelihood;
            mean = gmm1.Means[0, 0];
            standardDeviation = Math.Sqrt(gmm1.Sigmas[0][0, 0]);

            // Fit 2-component GMM (potentially bimodal)
            var gmm2 = new GaussianMixtureModel(data, 2);
            gmm2.Train(12345, true);
            var logLH2 = gmm2.LogLikelihood;

            // Compute test statistic: Likelihood Ratio Test
            return 2 * (logLH2 - logLH1);
        }

    }


}