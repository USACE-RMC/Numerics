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
using Numerics.Distributions;
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.SpecialFunctions;

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
            var meanVar = Statistics.MeanVariance(sample);
            int N = sample.Count;
            double se = Math.Sqrt(meanVar.Item2) / Math.Sqrt(N);
            double t = Math.Abs((meanVar.Item1 - populationMean) / se);
            var tdist = new StudentT(N - 1);
            return (1d - tdist.CDF(t)) * 2d;
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
            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            int N1 = sample1.Count;
            int N2 = sample2.Count;
            double df = N1 + N2 - 2;
            double svar = ((N1 - 1) * meanVar1.Item2 + (N2 - 1) * meanVar2.Item2) / df;
            double t = Math.Abs((meanVar1.Item1 - meanVar2.Item1) / Math.Sqrt(svar * (1.0d / N1 + 1.0d / N2)));
            var tdist = new StudentT(df);
            return (1d - tdist.CDF(t)) * 2d;
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
            var tdist = new StudentT((int)df);
            return (1d - tdist.CDF(t)) * 2d;
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
            var tdist = new StudentT((int)df);
            return (1d - tdist.CDF(t)) * 2d;
        }

        /// <summary>
        /// The F-test for significantly different variances. One rejects H0 if the p-value of the statistic is less than the significance level.
        /// </summary>
        /// <param name="sample1">Data sample 1.</param>
        /// <param name="sample2">Data sample 2.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double Ftest(IList<double> sample1, IList<double> sample2)
        {
            var meanVar1 = Statistics.MeanVariance(sample1);
            var meanVar2 = Statistics.MeanVariance(sample2);
            int n1 = sample1.Count, n2 = sample2.Count;
            double ave1 = meanVar1.Item1, var1 = meanVar1.Item2, ave2 = meanVar2.Item1, var2 = meanVar2.Item2, df1, df2, f, pVal;
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
        /// <param name="dfRestricted">The degrees of freedom of the restricted model.</param>
        /// <param name="dfFull">The degrees of freedom of the full model. </param>
        /// <param name="fStat">The F-test statistic.</param>
        /// <param name="pValue">The p-value of the test statistic.</param>
        public static void FtestModels(double sseRestricted, double sseFull, int dfRestricted, int dfFull, out double fStat, out double pValue)
        {
            fStat = ((sseRestricted - sseFull) / (dfRestricted - dfFull)) / (sseFull / dfFull);
            pValue = 2.0 * Beta.Incomplete(0.5 * dfFull, 0.5 * dfRestricted, dfFull / (dfFull + dfRestricted * fStat));
            if (pValue > 1.0) pValue = 2d - pValue;
        }

        /// <summary>
        /// The Jarque-Bera test for normality. One rejects H0 of normality if the p-value of the JB statistic is less than the significance level.
        /// </summary>
        /// <param name="sample">The data sample.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        /// <remarks>The Jarque-Bera normality test statistics is a chi-squared random variable with 2 degrees of freedom.</remarks>
        public static double JarqueBeraTest(IList<double> sample)
        {
            var moments = Statistics.ProductMoments(sample);
            double S2 = moments[2] * moments[2];
            double K2 = moments[3] * moments[3];
            int N = sample.Count;
            double JB = N / 6.0d * (S2 + K2 / 4.0d);
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
            double Q = 0;
            for (int k = 1; k <= lagMax; k++)
                Q += Tools.Sqr(acf[k, 1]) / (n - k);
            Q *= n * (n + 2);
            var chi2 = new ChiSquared(lagMax);
            return 1d - chi2.CDF(Q);
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
            var ties = new double[n1]; double R = 0, T = 0;

            var ranks = Statistics.RanksInPlace(sample.ToArray(), out ties);

            for (int i = 0; i < sample1.Count; i++) R += ranks[i];
            for (int i = 0; i < ties.Length; i++) T += (Tools.Pow(ties[i], 3) - ties[i]) / (n * (n - 1));

            double V = R - n1 * (n1 + 1d) / 2d;
            double W = n1 * n2 - V;
            double U = Math.Min(V, W);
            double eU = n1 * n2 / 2;
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

            var ties = new double[n];
            var R = Statistics.RanksInPlace(sample.ToArray(), out ties);
            for (i = 0; i < ties.Length; i++) T += ties[i] * (ties[i] - 1) * (2 * ties[i] + 5);
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
            var tdist = new StudentT(lm.DegreesOfFreedom);
            double d = Math.Abs(lm.Parameters[1] / lm.ParameterStandardErrors[1]);
            return (1 - tdist.CDF(Math.Abs(lm.Parameters[1] / lm.ParameterStandardErrors[1]))) * 2;
        }

        /// <summary>
        /// The Gaussian Mixture Model Unimodality Test.
        /// </summary>
        /// <param name="sample">Time series sample data.</param>
        /// <returns>Returns the p-value of the test statistic.</returns>
        public static double UnimodalityTest(IList<double> sample)
        {
            int n = sample.Count;
            if (n < 10) throw new ArgumentException("The sample size must be greater than or equal to 10.");

            // Fit 1-component GMM (unimodal)
            var gmm1 = new GaussianMixtureModel(sample.ToArray(), 1);
            gmm1.Train();
            var logLH1 = gmm1.LogLikelihood;

            // Fit 2-component GMM (potentially bimodal)
            var gmm2 = new GaussianMixtureModel(sample.ToArray(), 2);
            gmm2.Train();
            var logLH2 = gmm2.LogLikelihood;

            // Compute test statistic: Likelihood Ratio Test
            double testStat = 2 * (logLH2 - logLH1);
            int df = 3;

            // Compute p-value from chi-square distribution
            var chiSquared = new ChiSquared(df);
            double pval = 1.0 - chiSquared.CDF(testStat);

            return pval;

        }


    }


}