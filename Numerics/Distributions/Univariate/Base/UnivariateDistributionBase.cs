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

using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Xml.Linq;

namespace Numerics.Distributions
{

    /// <summary>
    /// Declares common functionality for all univariate distributions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    [Serializable]
    public abstract class UnivariateDistributionBase : IUnivariateDistribution, IComparable<UnivariateDistributionBase>
    {
        /// <summary>
        /// Approximation for assessing if a parameter is near zero. 
        /// </summary>
        protected double NearZero = 1E-4;

        /// <summary>
        /// Protected parameter is valid property. 
        /// </summary>
        protected bool _parametersValid = true;

        /// <summary>
        /// Returns the continuous distribution type.
        /// </summary>
        public abstract UnivariateDistributionType Type { get; }

        /// <inheritdoc/>
        public abstract string DisplayName { get; }

        /// <inheritdoc/>
        public abstract string ShortDisplayName { get; }

        /// <inheritdoc/>
        public string DisplayLabel
        {
            get
            {
                StringBuilder label = new StringBuilder(ShortDisplayName + " (");
                var paramVals = GetParameters;
                for (int i = 0; i < paramVals.Length; i++)
                {
                    label.Append(string.Format("{0:0.#####}", paramVals[i]));
                    if (i != paramVals.Length - 1)
                        label.Append(", ");
                }
                label.Append(")");
                return label.ToString();
            }
        }

        /// <inheritdoc/>
        public abstract int NumberOfParameters { get; }

        /// <inheritdoc/>
        public abstract string[,] ParametersToString { get; }

        /// <summary>
        /// Returns the distribution parameter names as an array of string.
        /// </summary>
        public virtual string[] ParameterNames
        {
            get { return ParametersToString.GetColumn(0); }
        }

        /// <inheritdoc/>
        public abstract string[] ParameterNamesShortForm { get; }

        /// <inheritdoc/>
        public abstract double[] GetParameters { get; }

        /// <inheritdoc/>
        public abstract string[] GetParameterPropertyNames { get; }

        /// <inheritdoc/>
        public bool ParametersValid => _parametersValid;

        /// <inheritdoc/>
        public abstract double Mean { get; }

        /// <inheritdoc/>
        public abstract double Median { get; }

        /// <inheritdoc/>
        public abstract double Mode { get; }

        /// <inheritdoc/>
        public double Variance
        {
            get { return Math.Pow(StandardDeviation, 2d); }
        }

        /// <inheritdoc/>
        public abstract double StandardDeviation { get; }

        /// <inheritdoc/>
        public double CoefficientOfVariation
        {
            get { return StandardDeviation / Mean; }
        }

        /// <inheritdoc/>
        public abstract double Skewness { get; }

        /// <inheritdoc/>
        public abstract double Kurtosis { get; }

        /// <inheritdoc/>
        public abstract double Minimum { get; }

        /// <inheritdoc/>
        public abstract double Maximum { get; }

        /// <inheritdoc/>
        public abstract double[] MinimumOfParameters { get; }

        /// <inheritdoc/>
        public abstract double[] MaximumOfParameters { get; }

        /// <inheritdoc/>
        public abstract void SetParameters(IList<double> parameters);

        /// <inheritdoc/>
        public abstract ArgumentOutOfRangeException ValidateParameters(IList<double> parameters, bool throwException);

        /// <summary>
        /// The log likelihood function.
        /// </summary>
        /// <param name="sample">Sample of observed data.</param>
        /// <returns>
        /// Log likelihood.
        /// </returns>
        public double LogLikelihood(IList<double> sample)
        {
            var logLH = default(double);
            for (int i = 0; i < sample.Count; i++)
            {        
                logLH += LogPDF(sample[i]);
            }
            if (double.IsNaN(logLH) || double.IsInfinity(logLH)) return double.MinValue;
            return logLH;
        }

        /// <summary>
        /// The Log-Likelihood function for a single data point.
        /// </summary>
        /// <param name="value">A single observed data point.</param>
        /// <returns>
        /// Log likelihood.
        /// </returns>
        public double LogLikelihood(double value)
        {
            return LogPDF(value);
        }

        /// <summary>
        /// The log likelihood function for left censored data.
        /// </summary>
        /// <param name="threshold">The threshold.</param>
        /// <param name="numberBelow">The number of data points below the threshold.</param>
        public double LogLikelihood_LeftCensored(double threshold, long numberBelow)
        {
            return numberBelow * LogCDF(threshold);
        }

        /// <summary>
        /// The log likelihood function for right censored data.
        /// </summary>
        /// <param name="threshold">The threshold.</param>
        /// <param name="numberAbove">The number of data points above the threshold.</param>
        public double LogLikelihood_RightCensored(double threshold, long numberAbove)
        {
            return numberAbove * LogCCDF(threshold);
        }

        /// <summary>
        /// The log likelihood function for interval data.
        /// </summary>
        /// <param name="lowerLimit">The lower limit of the interval.</param>
        /// <param name="upperLimit">The upper limit of the interval.</param>
        public double LogLikelihood_Intervals(double lowerLimit, double upperLimit)
        {
            double interval = CDF(upperLimit) - CDF(lowerLimit);
            return Math.Log(interval);
        }

        /// <inheritdoc/>
        public abstract double PDF(double x);

        /// <inheritdoc/>
        public virtual double LogPDF(double x)
        {
            double f = PDF(x);
            return Math.Log(f);
        }

        /// <inheritdoc/>
        public abstract double CDF(double x);

        /// <inheritdoc/>
        public virtual double HF(double x)
        {
            return PDF(x) / CCDF(x);
        }

        /// <inheritdoc/>
        public virtual double LogCDF(double x)
        {
            double F = CDF(x);
            return Math.Log(F);
        }

        /// <inheritdoc/>
        public virtual double CCDF(double x)
        {
            return 1.0d - CDF(x);
        }

        /// <inheritdoc/>
        public virtual double LogCCDF(double x)
        {
            double cF = CCDF(x);
            return Math.Log(cF);
        }

        /// <inheritdoc/>
        public abstract double InverseCDF(double probability);

        /// <summary>
        /// The PDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        public double[] PDF(IList<double> xValues)
        {
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = PDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The CDF evaluated over a list of x values.
        /// </summary>
        /// <param name="xValues">List of x values.</param>
        public double[] CDF(IList<double> xValues)
        {
            var result = new double[xValues.Count];
            for (int i = 0; i < xValues.Count; i++)
                result[i] = CDF(xValues[i]);
            return result;
        }

        /// <summary>
        /// The Inverse CDF evaluated over a list of probabilities.
        /// </summary>
        /// <param name="probabilities">List of probability values.</param>
        public double[] InverseCDF(IList<double> probabilities)
        {
            var result = new double[probabilities.Count];
            for (int i = 0; i < probabilities.Count; i++)
                result[i] = InverseCDF(probabilities[i]);
            return result;
        }

        /// <summary>
        /// Returns the central moments {Mean, Standard Deviation, Skew, and Kurtosis} of the distribution using numerical integration with Adaptive Simpson's rule. 
        /// </summary>
        /// <param name="tolerance">The desired tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.</param>
        /// <returns>Mean, Standard Deviation, Skew, and Kurtosis.</returns>
        public virtual double[] CentralMoments(double tolerance = 1E-8)
        {
            double u1 = 0, u2 = 0, u3 = 0, u4 = 0;
            double a = InverseCDF(1E-16);
            double b = InverseCDF(1-1E-16);
            if (a >= b) return new double[] { a, double.NaN, double.NaN, double.NaN };

            // Mean
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u1 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            // Standard Deviation
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - u1, 2d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u2 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            u2 = Math.Sqrt(u2);
            // Skewness
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow((x - u1) / u2, 3d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u3 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            // Kurtosis
            agk = new AdaptiveGaussKronrod(x => { return Math.Pow((x - u1) / u2, 4d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            u4 = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;

            return new double[] { u1, u2, u3, u4 };
        }

        /// <summary>
        /// Returns the central moments {Mean, Standard Deviation, Skew, and Kurtosis} of the distribution using numerical integration with Trapezoidal rule. 
        /// </summary>
        ///<param name="steps">Number of integration steps. Default = 300.</param>
        public virtual double[] CentralMoments(int steps = 300)
        {
            double a = InverseCDF(1E-8);
            double b = InverseCDF(1 - 1E-8);
            if (a >= b) return [a, double.NaN, double.NaN, double.NaN];

            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));
            var dFx = new double[steps];
            double u1, u2, u3, u4;
            double sumU1 = 0;
            double sumU2 = 0;
            double sumU3 = 0;
            double sumU4 = 0;

            // First compute the mean and standard deviation
            dFx[0] = CDF(bins[0].UpperBound);
            sumU1 += bins[0].UpperBound * dFx[0];
            sumU2 += Math.Pow(bins[0].UpperBound, 2d) * dFx[0];
            for (int i = 1; i < steps - 1; i++)
            {
                dFx[i] = CDF(bins[i].UpperBound) - CDF(bins[i].LowerBound);
                sumU1 += bins[i].Midpoint * dFx[i];
                sumU2 += Math.Pow(bins[i].Midpoint, 2d) * dFx[i];
            }
            dFx[steps - 1] = 1 - CDF(bins.Last().LowerBound);
            sumU1 += bins.Last().LowerBound * dFx[steps - 1];
            sumU2 += Math.Pow(bins.Last().LowerBound, 2d) * dFx[steps - 1];
            u1 = sumU1;
            u2 = Math.Sqrt(sumU2 - Math.Pow(u1, 2d));

            // Then compute skewness and kurtosis
            sumU3 += Math.Pow((bins[0].UpperBound - u1) / u2, 3d) * dFx[0];
            sumU4 += Math.Pow((bins[0].UpperBound - u1) / u2, 4d) * dFx[0];
            for (int i = 1; i < steps - 1; i++)
            {
                sumU3 += Math.Pow((bins[i].Midpoint - u1) / u2, 3d) * dFx[i];
                sumU4 += Math.Pow((bins[i].Midpoint - u1) / u2, 4d) * dFx[i];
            }
            sumU3 += Math.Pow((bins.Last().LowerBound - u1) / u2, 3d) * dFx[steps - 1];
            sumU4 += Math.Pow((bins.Last().LowerBound - u1) / u2, 4d) * dFx[steps - 1];
            u3 = sumU3;
            u4 = sumU4;
            return [u1, u2, u3, u4];
        }

        /// <summary>
        /// Returns the conditional moments of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <summary>
        /// Returns conditional central moments (up to 4th order) between [a, b].
        /// </summary>
        public virtual double[] ConditionalMoments(double a, double b)
        {
            if (a >= b)
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            int steps = 300;
            double l = CDF(a), u = CDF(b), totalProb = u - l;
            if (totalProb <= 0.0 || double.IsNaN(totalProb))
                return new[] { double.NaN, double.NaN, double.NaN, double.NaN };

            double mean = Mean;
            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));

            // Midpoint integration
            double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;

            for (int i = 0; i < steps; i++)
            {
                double dF = CDF(bins[i].UpperBound) - CDF(bins[i].LowerBound);
                double xMid = bins[i].Midpoint;

                double delta = xMid - mean;
                double delta2 = delta * delta;
                double delta3 = delta2 * delta;
                double delta4 = delta2 * delta2;

                sum1 += xMid * dF;
                sum2 += delta2 * dF;
                sum3 += delta3 * dF;
                sum4 += delta4 * dF;
            }

            double m1 = sum1 / totalProb;
            double m2 = sum2 / totalProb;
            double m3 = sum3 / totalProb;
            double m4 = sum4 / totalProb;

            return new[] { m1, m2, m3, m4 };
        }

        /// <summary>
        /// Returns the conditional mean of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalMean(double a, double b, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional variance of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalVariance(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 2d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional skewness of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalSkewness(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 3d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional kurtosis of the distribution.
        /// </summary>
        /// <param name="a">The lower integration limit, a.</param>
        /// <param name="b">The upper integration limit, b.</param>
        /// <param name="mean">The mean for computing the central moment. Default = 0.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = 1E-8.</param>
        public virtual double ConditionalKurtosis(double a, double b, double mean = 0, double tolerance = 1E-8)
        {
            double l = CDF(a);
            double u = CDF(b);
            var agk = new AdaptiveGaussKronrod(x => { return Math.Pow(x - mean, 4d) * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (u - l);
        }

        /// <summary>
        /// Returns the conditional expected value at the given confidence level. 
        /// </summary>
        /// <param name="alpha">The threshold confidence level, α.</param>
        /// <param name="tolerance">The desired tolerance for the solution. Default = ~Sqrt(Machine Epsilon), or 1E-8.</param>
        public virtual double ConditionalExpectedValue(double alpha, double tolerance = 1E-8)
        {
            double a = InverseCDF(alpha);
            double b = InverseCDF(1 - 1E-16);
            var agk = new AdaptiveGaussKronrod(x => { return x * PDF(x); }, a, b) { RelativeTolerance = tolerance, ReportFailure = false };
            agk.Integrate();
            double e = agk.Status != IntegrationStatus.Failure ? agk.Result : double.NaN;
            return e / (1 - alpha);
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and probability density.</returns>
        public virtual double[,] CreatePDFGraph()
        {
            var XValues = Stratify.XValues(new StratificationOptions(InverseCDF(0.0005d), InverseCDF(0.9995d), 1000, false));
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = PDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Stratified X values.</param>
        /// <returns>A 2-column array of X and probability density.</returns>
        public double[,] CreatePDFGraph(List<StratificationBin> XValues)
        {
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = PDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a PDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Array of stratified X values.</param>
        /// <returns>A 2-column array of X and probability density.</returns>
        public double[,] CreatePDFGraph(List<StratificationBin>[] XValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < XValues.Length; i++)
                totalBinCount += XValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < XValues.Length; i++)
            {
                for (int j = 0; j < XValues[i].Count; j++)
                {
                    graph[binCounter, 0] = XValues[i][j].Midpoint;
                    graph[binCounter, 1] = PDF(XValues[i][j].Midpoint);
                    binCounter += 1;
                }
                binCounter += 1;
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateCDFGraph()
        {
            var XValues = Stratify.XValues(new StratificationOptions(InverseCDF(0.0005d), InverseCDF(0.9995d), 1000, false));
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = CDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Stratified X values.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateCDFGraph(List<StratificationBin> XValues)
        {
            var graph = new double[XValues.Count, 2];
            for (int i = 0; i < XValues.Count; i++)
            {
                graph[i, 0] = XValues[i].Midpoint;
                graph[i, 1] = CDF(XValues[i].Midpoint);
            }
            return graph;
        }

        /// <summary>
        /// Create a CDF table for graphing purposes.
        /// </summary>
        /// <param name="XValues">Array of stratified X values.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateCDFGraph(List<StratificationBin>[] XValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < XValues.Length; i++)
                totalBinCount += XValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < XValues.Length; i++)
            {
                for (int j = 0; j < XValues[i].Count; j++)
                {
                    graph[binCounter, 0] = XValues[i][j].Midpoint;
                    graph[binCounter, 1] = CDF(XValues[i][j].Midpoint);
                    binCounter += 1;
                }
                binCounter += 1;
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// The bounds of the table are automatically determined.
        /// </summary>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph()
        {
            var PValues = Stratify.Probabilities(new StratificationOptions(CDF(InverseCDF(0.0005d)), CDF(InverseCDF(0.9995d)), 1000, true));
            var graph = new double[PValues.Count, 2];
            for (int i = 0; i < PValues.Count; i++)
            {
                graph[i, 0] = InverseCDF(PValues[i].Midpoint);
                graph[i, 1] = PValues[i].Midpoint;
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// </summary>
        /// <param name="PValues">Stratified probabilities.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph(List<StratificationBin> PValues)
        {
            var graph = new double[PValues.Count, 2];
            for (int i = 0; i < PValues.Count; i++)
            {
                graph[i, 0] = InverseCDF(PValues[i].Midpoint);
                graph[i, 1] = PValues[i].Midpoint;
            }
            return graph;
        }

        /// <summary>
        /// Create an Inverse CDF table for graphing purposes.
        /// </summary>
        /// <param name="PValues">Array of stratified probabilities.</param>
        /// <returns>A 2-column array of X and non-exceedance probability.</returns>
        public double[,] CreateInverseCDFGraph(List<StratificationBin>[] PValues)
        {
            var totalBinCount = default(int);
            for (int i = 0; i < PValues.Length; i++)
                totalBinCount += PValues[i].Count;
            var graph = new double[totalBinCount, 2];
            int binCounter = 0;
            for (int i = 0; i < PValues.Length; i++)
            {
                for (int j = 0; j < PValues[i].Count; j++)
                {
                    graph[binCounter, 0] = InverseCDF(PValues[i][j].Midpoint);
                    graph[binCounter, 1] = PValues[i][j].Midpoint;
                    binCounter += 1;
                }
            }
            return graph;
        }

        /// <summary>
        /// Generate random values of a distribution given a sample size.
        /// </summary>
        /// <param name="sampleSize">Size of random sample to generate.</param>
        /// <param name="seed">Optional. The prng seed. If negative or zero, then the computer clock is used as a seed.</param>
        /// <returns>
        /// Array of random values.
        /// </returns>
        public virtual double[] GenerateRandomValues(int sampleSize, int seed = -1)
        {
            // Create PRNG for generating random numbers
            var rnd = seed > 0 ? new MersenneTwister(seed) : new MersenneTwister();
            var sample = new double[sampleSize];
            // Generate values
            for (int i = 0; i < sampleSize; i++)
                sample[i] = InverseCDF(rnd.NextDouble());
            // Return array of random values
            return sample;
        }

        /// <summary>
        /// Creates a copy of the distribution.
        /// </summary>
        public abstract UnivariateDistributionBase Clone();

        /// <summary>
        /// Returns the distribution as XElement (XML). 
        /// </summary>
        public virtual XElement ToXElement()
        {
            var result = new XElement("Distribution");
            result.SetAttributeValue(nameof(Type), Type.ToString());
            var names = GetParameterPropertyNames;
            var parms = GetParameters;
            for (int i = 0; i < NumberOfParameters; i++)
            {
                result.SetAttributeValue(names[i], parms[i].ToString("G17", CultureInfo.InvariantCulture));
            }
            return result;
        }

        /// <summary>
        /// The equality operator to allow comparison of continuous distribution parameters. If distributions are equal the operator will return true, otherwise false.
        /// </summary>
        /// <param name="left">Continuous distribution on the left side of the operator.</param>
        /// <param name="right">Continuous distribution on the right side of the operator.</param>
        /// <returns>Boolean indicating if the equality operator results in a true statement.</returns>
        public static bool operator ==(UnivariateDistributionBase left, UnivariateDistributionBase right)
        {
            // Check for null arguments. Keep in mind null == null
            if (left is null && right is null)
            {
                return true;
            }
            else if (left is null)
            {
                return false;
            }
            else if (right is null)
            {
                return false;
            }
            // Check to see if they are the same distribution
            if (left.Type != right.Type)
                return false;
            // Test for equality of parameters
            // Compare parameters here. 
            // Need to have special cases for kernel density, univariate, and bivariate distributions. Everything else can use the GetParameters method.
            if (left.Type == UnivariateDistributionType.Empirical)
            {
                EmpiricalDistribution leftTyped = (EmpiricalDistribution)left;
                EmpiricalDistribution rightTyped = (EmpiricalDistribution)right;
                // Check X values for equality
                var leftValues = leftTyped.XValues;
                var rightValues = rightTyped.XValues;
                if (leftValues.Count != rightValues.Count)
                    return false;
                for (int i = 0; i < leftValues.Count; i++)
                {
                    if (Math.Abs(leftValues[i] - rightValues[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
                // Check P values for equality, no need to test count since it was tested with X values
                leftValues = leftTyped.ProbabilityValues;
                rightValues = rightTyped.ProbabilityValues;
                for (int i = 0; i < leftValues.Count; i++)
                {
                    if (Math.Abs(leftValues[i] - rightValues[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
            }
            else if (left.Type == UnivariateDistributionType.KernelDensity)
            {
                throw new NotImplementedException();
            }
            else
            {
                var leftParams = left.GetParameters;
                var rightParams = right.GetParameters;
                for (int i = 0; i < leftParams.Count(); i++)
                {
                    if (Math.Abs(leftParams[i] - rightParams[i]) > Tools.DoubleMachineEpsilon)
                        return false;
                }
            }
            // All tests passed
            return true;
            // 
        }

        /// <summary>
        /// The non-equality operator to allow comparison of continuous distribution parameters. If distributions are not equal the operator will return true, otherwise false.
        /// </summary>
        /// <param name="left">Continuous distribution on the left side of the operator.</param>
        /// <param name="right">Continuous distribution on the right side of the operator.</param>
        /// <returns>Boolean indicating if the non-equality operator results in a true statement.</returns>
        public static bool operator !=(UnivariateDistributionBase left, UnivariateDistributionBase right)
        {
            return !(left == right);
        }

        /// <summary>
        /// Determines whether the specified object is equal to the current object.
        /// </summary>
        /// <param name="obj">The object to compare with the current object.</param>
        /// <returns>True if the specified object is equal to the current object; otherwise, False.</returns>
        public override bool Equals(object? obj)
        {
            if (obj is UnivariateDistributionBase other)
            {
                return this == other;
            }
            return false;
        }

        /// <summary>
        /// Serves as the default hash function.
        /// </summary>
        /// <returns>A hash code for the current object.</returns>
        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                hash = hash * 23 + Type.GetHashCode();

                if (Type == UnivariateDistributionType.Empirical)
                {
                    var empirical = this as EmpiricalDistribution;
                    if (empirical is not null)
                    {
                        foreach (var x in empirical.XValues)
                        {
                            hash = hash * 23 + x.GetHashCode();
                        }
                        foreach (var p in empirical.ProbabilityValues)
                        {
                            hash = hash * 23 + p.GetHashCode();
                        }
                    }
                }
                else if (Type != UnivariateDistributionType.KernelDensity)
                {
                    var parameters = GetParameters;
                    foreach (var param in parameters)
                    {
                        hash = hash * 23 + param.GetHashCode();
                    }
                }

                return hash;
            }
        }

        /// <summary>
        /// Implementation of IComparable to sort distributions by their Mean.
        /// </summary>
        /// <param name="other">Another UnivariateDistributionBase instance.</param>
        /// <returns>
        /// A negative number if this instance's mean is less than the other; zero if equal; a positive number if greater.
        /// </returns>
        public int CompareTo(UnivariateDistributionBase? other)
        {
            if (other is null)
                return 1; // non-null instance is considered greater than null

            return this.Mean.CompareTo(other.Mean);
        }

    }
}