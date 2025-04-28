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

using Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Numerics.Data.Statistics
{

    /// <summary>
    /// A class containing goodness-of-fit measures.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class GoodnessOfFit
    {

        /// <summary>
        /// Gets the Akaike Information Criterion (AIC) used for model selection among a finite set of models; the model with the lowest AIC is preferred.
        /// When comparing multiple model fits, additional model parameters often yield larger, optimized log-likelihood values.
        /// Unlike the optimized log-likelihood value, AIC penalizes for more complex models, i.e., models with additional parameters.
        /// </summary>
        /// <param name="numberOfParameters">The number of model parameters</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The AIC value for the model</returns>
        public static double AIC(int numberOfParameters, double logLikelihood)
        {
            return (-2d * logLikelihood) + (2d * numberOfParameters);
        }

        /// <summary>
        /// Gets the Akaike Information Criterion (AIC), corrected for small sample sizes, used for model selection among a finite set of models; the model with the lowest AIC is preferred.
        /// When comparing multiple model fits, additional model parameters often yield larger, optimized log-likelihood values.
        /// Unlike the optimized log-likelihood value, AIC penalizes for more complex models, i.e., models with additional parameters.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="numberOfParameters">The number of model parameters</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The AIC value for the model</returns>
        public static double AICc(int sampleSize, int numberOfParameters, double logLikelihood)
        {
            double _aic = AIC(numberOfParameters, logLikelihood);
            // Make adjustment for small samples
            return _aic + (2d * Tools.Sqr(numberOfParameters) + 2d * numberOfParameters) / (sampleSize - numberOfParameters - 1);
        }

        /// <summary>
        /// Gets the Bayesian information criterion (BIC) used for model selection among a finite set of models; the model with the lowest BIC is preferred.
        /// Like AIC, BIC uses the optimal log-likelihood function value and penalizes for more complex models, i.e., models with additional parameters.
        /// The penalty of BIC is a function of the sample size, and so is typically more severe than that of AIC.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="numberOfParameters">The number of model parameters</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The BIC value for the model</returns>
        public static double BIC(int sampleSize, int numberOfParameters, double logLikelihood)
        {
            return (-2d * logLikelihood) + (numberOfParameters * Math.Log(sampleSize));
        }

        /// <summary>
        /// Returns an array of weights based on a list of model AIC values.
        /// </summary>
        /// <param name="aicValues">The list of model AIC values.</param>
        public static double[] AICWeights(IList<double> aicValues)
        {
            var min = Tools.Min(aicValues);
            var weights = new double[aicValues.Count];
            var num = new double[aicValues.Count];
            double sum = 0;
            for (int i = 0; i < aicValues.Count; i++)
            {
                num[i] = Math.Exp(-0.5 * (aicValues[i] - min));
                sum += num[i];
            }
            for (int i = 0; i < aicValues.Count; i++)
            {
                weights[i] = num[i] / sum;
            }
            return weights;
        }

        /// <summary>
        /// Gets the Root Mean Square Error (RMSE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <param name="k">Number of model parameters. Default = 0.</param>
        /// <returns>The RMSE of the model</returns>
        public static double RMSE(IList<double> observedValues, IList<double> modeledValues, int k = 0)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");
            
            int n = observedValues.Count - k;
            double sse = 0d;
            for (int i = 0; i < n; i++)
                sse += Tools.Sqr(modeledValues[i] - observedValues[i]);
            return Math.Sqrt(sse / n);
        }

        /// <summary>
        /// Gets the Root Mean Square Error (RMSE) of the model compared to the observed data. Weibull plotting positions are assumed.
        /// </summary>inverse
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The RMSE of the model</returns>
        public static double RMSE(IList<double> observedValues, UnivariateDistributionBase model)
        {
            var observed = observedValues.ToArray();
            Array.Sort(observed);
            var pp = PlottingPositions.Weibull(observed.Length);
            var modeled = model.InverseCDF(pp);
            return RMSE(observed, modeled, model.NumberOfParameters);
        }

        /// <summary>
        /// Gets the Root Mean Square Error (RMSE) of the model compared to the observed data. 
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="plottingPositions">The plotting positions of the observed values.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The RMSE of the model</returns>
        public static double RMSE(IList<double> observedValues, IList<double> plottingPositions, UnivariateDistributionBase model)
        {
            var modeled = model.InverseCDF(plottingPositions);
            return RMSE(observedValues, modeled, model.NumberOfParameters);
        }

        /// <summary>
        /// Returns an array of weights based on a list of model RMSE values. Weights are derived using inverse-MSE weighting.
        /// </summary>
        /// <param name="rmseValues">The list of model RMSE values.</param>
        public static double[] RMSEWeights(IList<double> rmseValues)
        {
            var weights = new double[rmseValues.Count];
            var invMSE = new double[rmseValues.Count];
            double sum = 0;
            for (int i = 0; i < rmseValues.Count; i++)
            {
                invMSE[i] = 1 / (rmseValues[i] * rmseValues[i]);
                sum += invMSE[i];
            }
            for (int i = 0; i < rmseValues.Count; i++)
            {
                weights[i] = invMSE[i] / sum;
            }
            return weights;
        }

        /// <summary>
        /// Gets the R^2, the square of the correlation of the model data compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The R^2 value pf the model</returns>
        public static double RSquared(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            var corr = Correlation.Pearson(observedValues, modeledValues);
            return corr * corr;

        }

        /// <summary>
        /// Gets the accuracy of the modeled data compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>Th accuracy as percent.</returns>
        public static double Accuracy(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            double accuracy = 0;
            for (int i = 0; i < observedValues.Count; ++i)
            {
                if (observedValues[i]  == modeledValues[i])
                {
                    accuracy++;
                }
            }
            return 100 * accuracy / observedValues.Count;
        }

        /// <summary>
        /// Computes the Kolmogorov-Smirnov test statistic, which is used to decide if a sample comes from a population with a specified distribution.
        /// </summary>
        /// <remarks>
        /// <b> References: </b> 
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm"/>
        /// </remarks>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The Kolmogorov Smirnov test statistic</returns>
        public static double KolmogorovSmirnov(IList<double> observedValues, UnivariateDistributionBase model)
        {
            // Check if the lists are the same size
            if (observedValues.Count < 1)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "There must be more than one observed value.");

            int n = observedValues.Count;
            double D = double.MinValue;
            for (int i = 1; i <= n; i++)
            {
                double left = model.CDF(observedValues[i - 1]) - ((double)i - 1) / n;
                double right = (double)i / n - model.CDF(observedValues[i - 1]);
                D = Math.Max(D, Math.Max(left, right));
            }
            return D;
        }

        /// <summary>
        /// Computes the Chi-Squared test statistic, which is used to decide if a sample of data comes from a population with a specified distribution.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm"/>
        /// </remarks>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The Chi-Squared test statistic.</returns>
        public static double ChiSquared(IList<double> observedValues, UnivariateDistributionBase model)
        {
            // Check if the lists are the same size
            if (observedValues.Count < 1)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "There must be more than one observed value.");

            int n = observedValues.Count;
            var hist = new Histogram(observedValues);
            double x2 = 0;
            for (int i = 0; i < hist.NumberOfBins; i++)
            {
                double e = n * (model.CDF(hist[i].UpperBound) - model.CDF(hist[i].LowerBound));
                x2 += Math.Pow(hist[i].Frequency - e, 2) / e;
            }
            return x2;
        }

        /// <summary>
        /// Computes the Anderson-Darling test statistic, which is used to decide if a sample of data comes from a population with a specified distribution.
        /// </summary>
        /// <remarks>
        /// <b> References: </b>
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm"/>
        /// </remarks>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The Anderson-Darling test statistic.</returns>
        public static double AndersonDarling(IList<double> observedValues, UnivariateDistributionBase model)
        {
            // Check if the lists are the same size
            if (observedValues.Count < 1)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "There must be more than one observed value.");

            int n = observedValues.Count;
            double S = 0;
            for (int i = 1; i <= n; i++)
                S += (2 * (double)i - 1) * (model.LogCDF(observedValues[i - 1]) + model.LogCCDF(observedValues[n - i]));
            
            return -n - S / n;
        }

    }
}