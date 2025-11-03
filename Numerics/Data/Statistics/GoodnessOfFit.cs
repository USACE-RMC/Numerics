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
    /// A class containing goodness-of-fit measures for evaluating model performance.
    /// Includes information-theoretic criteria (AIC, BIC), error metrics (RMSE, MAE, MSE), 
    /// efficiency coefficients (Nash-Sutcliffe, Kling-Gupta), bias metrics (PBIAS), 
    /// and statistical tests (Kolmogorov-Smirnov, Anderson-Darling, Chi-Squared).
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    ///     <b> References: </b>
    ///     <list type="bullet">
    ///     <item><description>
    ///     Moriasi, D.N., Arnold, J.G., Van Liew, M.W., Bingner, R.L., Harmel, R.D., Veith, T.L. (2007). 
    ///     Model evaluation guidelines for systematic quantification of accuracy in watershed simulations. 
    ///     Transactions of the ASABE, 50(3), 885-900.
    ///     </description></item>
    ///     <item><description>
    ///     Moriasi, D.N., Gitau, M.W., Pai, N., Daggupati, P. (2015). 
    ///     Hydrologic and water quality models: Performance measures and evaluation criteria. 
    ///     Transactions of the ASABE, 58(6), 1763-1785.
    ///     </description></item>
    ///     <item><description>
    ///     Nash, J.E., Sutcliffe, J.V. (1970). 
    ///     River flow forecasting through conceptual models part I — A discussion of principles. 
    ///     Journal of Hydrology, 10(3), 282-290.
    ///     </description></item>
    ///     <item><description>
    ///     Gupta, H.V., Kling, H., Yilmaz, K.K., Martinez, G.F. (2009). 
    ///     Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling. 
    ///     Journal of Hydrology, 377(1-2), 80-91.
    ///     </description></item>
    ///     <item><description>
    ///     Kling, H., Fuchs, M., Paulin, M. (2012). 
    ///     Runoff conditions in the upper Danube basin under an ensemble of climate change scenarios. 
    ///     Journal of Hydrology, 424-425, 264-277.
    ///     </description></item>
    ///     <item><description>
    ///     Legates, D.R., McCabe, G.J. (1999). 
    ///     Evaluating the use of "goodness-of-fit" measures in hydrologic and hydroclimatic model validation. 
    ///     Water Resources Research, 35(1), 233-241.
    ///     </description></item>
    ///     </list>
    /// </para>
    /// </remarks>
    public class GoodnessOfFit
    {

        #region "Information Criteria"

        /// <summary>
        /// Computes the Akaike Information Criterion (AIC) used for model selection among a finite set of models.
        /// The model with the lowest AIC is preferred.
        /// </summary>
        /// <param name="numberOfParameters">The number of model parameters.</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The AIC value for the model.</returns>
        /// <remarks>
        /// AIC penalizes model complexity to prevent overfitting. When comparing multiple model fits, 
        /// additional model parameters often yield larger, optimized log-likelihood values. 
        /// Unlike the optimized log-likelihood value alone, AIC penalizes for more complex models, 
        /// i.e., models with additional parameters.
        /// </remarks>
        public static double AIC(int numberOfParameters, double logLikelihood)
        {
            return (-2d * logLikelihood) + (2d * numberOfParameters);
        }

        /// <summary>
        /// Computes the Akaike Information Criterion (AIC), corrected for small sample sizes (AICc), 
        /// used for model selection among a finite set of models. The model with the lowest AICc is preferred.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="numberOfParameters">The number of model parameters.</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The AICc value for the model.</returns>
        /// <remarks>
        /// AICc includes a correction term that becomes negligible as sample size increases. 
        /// It is recommended to use AICc instead of AIC when the sample size is small relative to 
        /// the number of parameters (n/k less than 40).
        /// </remarks>
        public static double AICc(int sampleSize, int numberOfParameters, double logLikelihood)
        {
            double _aic = AIC(numberOfParameters, logLikelihood);
            // Make adjustment for small samples
            return _aic + (2d * Tools.Sqr(numberOfParameters) + 2d * numberOfParameters) / (sampleSize - numberOfParameters - 1);
        }

        /// <summary>
        /// Computes the Bayesian Information Criterion (BIC) used for model selection among a finite set of models.
        /// The model with the lowest BIC is preferred.
        /// </summary>
        /// <param name="sampleSize">The sample size.</param>
        /// <param name="numberOfParameters">The number of model parameters.</param>
        /// <param name="logLikelihood">The maximum log-likelihood.</param>
        /// <returns>The BIC value for the model.</returns>
        /// <remarks>
        /// Like AIC, BIC uses the optimal log-likelihood function value and penalizes for more complex models. 
        /// The penalty of BIC is a function of the sample size, and is typically more severe than that of AIC, 
        /// especially for large sample sizes.
        /// </remarks>
        public static double BIC(int sampleSize, int numberOfParameters, double logLikelihood)
        {
            return (-2d * logLikelihood) + (numberOfParameters * Math.Log(sampleSize));
        }

        /// <summary>
        /// Computes an array of weights based on a list of model AIC values using Akaike weights.
        /// </summary>
        /// <param name="aicValues">The list of model AIC values.</param>
        /// <returns>An array of normalized weights that sum to 1.</returns>
        /// <remarks>
        /// Akaike weights provide a measure of the relative likelihood of each model being the best model 
        /// among the set of candidate models. The weight for model i is interpreted as the probability 
        /// that model i is the best model given the data and the set of candidate models.
        /// </remarks>
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

        #endregion

        #region "Error Metrics"

        /// <summary>
        /// Computes the Root Mean Square Error (RMSE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <param name="k">Number of model parameters. Default = 0.</param>
        /// <returns>The RMSE of the model.</returns>
        /// <remarks>
        /// RMSE is the square root of the average of squared differences between observed and modeled values. 
        /// It is sensitive to large errors due to the squaring operation and is expressed in the same units 
        /// as the observed data. Lower values indicate better model performance.
        /// </remarks>
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
        /// Computes the Root Mean Square Error (RMSE) of the model compared to the observed data. 
        /// Weibull plotting positions are assumed.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The RMSE of the model.</returns>
        /// <remarks>
        /// This method sorts the observed values and compares them to the model's inverse CDF 
        /// evaluated at Weibull plotting positions.
        /// </remarks>
        public static double RMSE(IList<double> observedValues, UnivariateDistributionBase model)
        {
            var observed = observedValues.ToArray();
            Array.Sort(observed);
            var pp = PlottingPositions.Weibull(observed.Length);
            var modeled = model.InverseCDF(pp);
            return RMSE(observed, modeled, model.NumberOfParameters);
        }

        /// <summary>
        /// Computes the Root Mean Square Error (RMSE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="plottingPositions">The plotting positions of the observed values.</param>
        /// <param name="model">The univariate continuous distribution.</param>
        /// <returns>The RMSE of the model.</returns>
        /// <remarks>
        /// This method allows for custom plotting positions to be specified for comparing 
        /// observed values to the model's inverse CDF.
        /// </remarks>
        public static double RMSE(IList<double> observedValues, IList<double> plottingPositions, UnivariateDistributionBase model)
        {
            var modeled = model.InverseCDF(plottingPositions);
            return RMSE(observedValues, modeled, model.NumberOfParameters);
        }

        /// <summary>
        /// Computes an array of weights based on a list of model RMSE values using inverse-MSE weighting.
        /// </summary>
        /// <param name="rmseValues">The list of model RMSE values.</param>
        /// <returns>An array of normalized weights that sum to 1.</returns>
        /// <remarks>
        /// Weights are derived using inverse mean squared error weighting. Models with lower RMSE 
        /// receive higher weights, providing a probabilistic interpretation of model performance.
        /// </remarks>
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
        /// Computes the Mean Squared Error (MSE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The MSE of the model.</returns>
        /// <remarks>
        /// MSE is the average of the squared differences between observed and modeled values. 
        /// It is sensitive to outliers due to the squaring operation and penalizes larger errors more heavily. 
        /// MSE is commonly used as a loss function in optimization algorithms because it is differentiable. 
        /// Units are the square of the original data units.
        /// </remarks>
        public static double MSE(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sse = 0d;
            for (int i = 0; i < n; i++)
                sse += Tools.Sqr(modeledValues[i] - observedValues[i]);
            return sse / n;
        }

        /// <summary>
        /// Computes the Mean Absolute Error (MAE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The MAE of the model.</returns>
        /// <remarks>
        /// MAE is the average of the absolute differences between observed and modeled values. 
        /// Unlike MSE or RMSE, MAE is less sensitive to outliers because it does not square the errors. 
        /// It treats all errors equally on a linear scale and is expressed in the same units as the observed data. 
        /// MAE is more robust to outliers and provides a direct, interpretable measure of average error magnitude.
        /// </remarks>
        public static double MAE(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sumAbsError = 0d;
            for (int i = 0; i < n; i++)
                sumAbsError += Math.Abs(modeledValues[i] - observedValues[i]);
            return sumAbsError / n;
        }

        /// <summary>
        /// Computes the Mean Absolute Percentage Error (MAPE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The MAPE as a percentage.</returns>
        /// <remarks>
        /// <para>
        /// MAPE expresses the average absolute error as a percentage of the observed values. 
        /// It is scale-independent, making it useful for comparing model performance across datasets 
        /// with different scales or units.
        /// </para>
        /// <para>
        /// Formula: MAPE = (100/n) × Σ|O - M| / |O|
        /// </para>
        /// <para>
        /// Interpretation:
        /// <list type="bullet">
        /// <item><description>MAPE = 0: Perfect predictions</description></item>
        /// <item><description>MAPE &lt; 10%: Highly accurate forecasting</description></item>
        /// <item><description>MAPE 10-20%: Good forecasting</description></item>
        /// <item><description>MAPE 20-50%: Reasonable forecasting</description></item>
        /// <item><description>MAPE &gt; 50%: Inaccurate forecasting</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// <b>Important Limitations:</b>
        /// <list type="bullet">
        /// <item><description>
        /// MAPE is undefined when observed values are zero and produces infinite or very large errors 
        /// when observed values are close to zero.
        /// </description></item>
        /// <item><description>
        /// MAPE has an asymmetric penalty: it penalizes over-predictions more heavily than under-predictions.
        /// </description></item>
        /// <item><description>
        /// Not recommended for data with zero or near-zero values. Consider using symmetric MAPE (sMAPE) 
        /// or weighted MAPE for such cases.
        /// </description></item>
        /// </list>
        /// </para>
        /// <para>
        /// This method will throw an exception if any observed value is zero.
        /// </para>
        /// </remarks>
        public static double MAPE(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sumPercentError = 0d;

            for (int i = 0; i < n; i++)
            {
                // Check for zero values which would cause division by zero
                if (Math.Abs(observedValues[i]) < double.Epsilon)
                    throw new ArgumentException("MAPE cannot be calculated when observed values contain zero or near-zero values.", nameof(observedValues));

                sumPercentError += Math.Abs((observedValues[i] - modeledValues[i]) / observedValues[i]);
            }

            return 100.0 * sumPercentError / n;
        }

        /// <summary>
        /// Computes the Symmetric Mean Absolute Percentage Error (sMAPE) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The sMAPE as a percentage.</returns>
        /// <remarks>
        /// <para>
        /// sMAPE is a variation of MAPE that addresses some of its limitations by using the average 
        /// of observed and modeled values in the denominator, making it more symmetric and bounded.
        /// </para>
        /// <para>
        /// Formula: sMAPE = (100/n) × Σ(2|O - M|) / (|O| + |M|)
        /// </para>
        /// <para>
        /// Advantages over MAPE:
        /// <list type="bullet">
        /// <item><description>Bounded between 0% and 200% (typically reported as 0% to 100%)</description></item>
        /// <item><description>More symmetric treatment of over- and under-predictions</description></item>
        /// <item><description>Less sensitive to zero values (though still problematic when both O and M are zero)</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Interpretation (when expressed as 0-100%):
        /// <list type="bullet">
        /// <item><description>sMAPE = 0%: Perfect predictions</description></item>
        /// <item><description>sMAPE &lt; 10%: Excellent forecasting</description></item>
        /// <item><description>sMAPE 10-20%: Good forecasting</description></item>
        /// <item><description>sMAPE 20-50%: Reasonable forecasting</description></item>
        /// <item><description>sMAPE &gt; 50%: Poor forecasting</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// This method will throw an exception if both observed and modeled values are zero for any data point.
        /// </para>
        /// </remarks>
        public static double sMAPE(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sumSymmetricPercentError = 0d;

            for (int i = 0; i < n; i++)
            {
                double denominator = Math.Abs(observedValues[i]) + Math.Abs(modeledValues[i]);

                // Check for case where both values are zero
                if (denominator < double.Epsilon)
                    throw new ArgumentException("sMAPE cannot be calculated when both observed and modeled values are zero.", nameof(observedValues));

                sumSymmetricPercentError += Math.Abs(observedValues[i] - modeledValues[i]) / denominator;
            }

            return 200.0 * sumSymmetricPercentError / n;
        }

        #endregion

        #region "Efficiency Coefficients"

        /// <summary>
        /// Computes the Nash-Sutcliffe Efficiency (NSE) coefficient.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The NSE coefficient.</returns>
        /// <remarks>
        /// <para>
        /// NSE ranges from -∞ to 1, where:
        /// <list type="bullet">
        /// <item><description>NSE = 1: Perfect model fit</description></item>
        /// <item><description>NSE = 0: Model performs as well as using the mean of observed data</description></item>
        /// <item><description>NSE &lt; 0: Mean of observed data is a better predictor than the model</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Performance ratings (Moriasi et al., 2007):
        /// <list type="bullet">
        /// <item><description>Very good: 0.75 &lt; NSE ≤ 1.00</description></item>
        /// <item><description>Good: 0.65 &lt; NSE ≤ 0.75</description></item>
        /// <item><description>Satisfactory: 0.50 &lt; NSE ≤ 0.65</description></item>
        /// <item><description>Unsatisfactory: NSE ≤ 0.50</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// NSE is sensitive to extreme values and systematic over- or under-prediction. 
        /// It is widely used in hydrological modeling for evaluating streamflow predictions.
        /// </para>
        /// </remarks>
        public static double NashSutcliffeEfficiency(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double observedMean = Statistics.Mean(observedValues);

            double numerator = 0d;
            double denominator = 0d;
            for (int i = 0; i < n; i++)
            {
                numerator += Tools.Sqr(observedValues[i] - modeledValues[i]);
                denominator += Tools.Sqr(observedValues[i] - observedMean);
            }

            return 1.0 - (numerator / denominator);
        }

        /// <summary>
        /// Computes the Log Nash-Sutcliffe Efficiency, which emphasizes low flow conditions.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <param name="epsilon">Small constant to add before taking logarithm (default = 0.01).</param>
        /// <returns>The Log-NSE coefficient.</returns>
        /// <remarks>
        /// <para>
        /// Log-NSE applies logarithmic transformation before calculating NSE, giving more weight 
        /// to low flow conditions and less weight to peak flows.
        /// </para>
        /// <para>
        /// Formula: Log-NSE = NSE(log(O + ε), log(M + ε))
        /// </para>
        /// <para>
        /// Use case: When accurate simulation of low flows is critical (e.g., drought analysis, 
        /// environmental flows, water quality)
        /// </para>
        /// <para>
        /// Note: Epsilon value should be small relative to typical flow values. Default is mean(observed)/100.
        /// </para>
        /// </remarks>
        public static double LogNashSutcliffeEfficiency(IList<double> observedValues,
            IList<double> modeledValues, double epsilon = -1)
        {
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues),
                    "The number of observed values must equal the number of modeled values.");

            // If epsilon not provided, use Pushpalatha et al. (2012) recommendation
            if (epsilon < 0)
                epsilon = Statistics.Mean(observedValues) / 100.0;

            var logObserved = new double[observedValues.Count];
            var logModeled = new double[modeledValues.Count];

            for (int i = 0; i < observedValues.Count; i++)
            {
                logObserved[i] = Math.Log(observedValues[i] + epsilon);
                logModeled[i] = Math.Log(modeledValues[i] + epsilon);
            }

            return NashSutcliffeEfficiency(logObserved, logModeled);
        }

        /// <summary>
        /// Computes the Kling-Gupta Efficiency (KGE) using the 2009 formulation.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The KGE coefficient, or a large negative value if calculation is undefined.</returns>
        /// <remarks>
        /// <para>
        /// KGE decomposes model performance into three components:
        /// <list type="bullet">
        /// <item><description>r: Linear correlation coefficient</description></item>
        /// <item><description>α (alpha): Ratio of standard deviations (modeled/observed)</description></item>
        /// <item><description>β (beta): Ratio of means (modeled/observed)</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// KGE ranges from -∞ to 1, where:
        /// <list type="bullet">
        /// <item><description>KGE = 1: Perfect model fit</description></item>
        /// <item><description>KGE &gt; -0.41: Model outperforms the mean flow benchmark</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// <b>Special Cases:</b>
        /// <list type="bullet">
        /// <item><description>If modeled values have zero variance (all constant), returns -10.0</description></item>
        /// <item><description>If observed values have zero variance, returns -10.0</description></item>
        /// <item><description>These cases indicate the model cannot capture any variability</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// KGE addresses known limitations of NSE by providing a more balanced evaluation of model performance. 
        /// It is less sensitive to systematic bias and better captures the variability of observations.
        /// Reference: Gupta et al. (2009)
        /// </para>
        /// </remarks>
        public static double KlingGuptaEfficiency(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;

            // First pass: compute means
            double sumObs = 0d;
            double sumMod = 0d;
            for (int i = 0; i < n; i++)
            {
                sumObs += observedValues[i];
                sumMod += modeledValues[i];
            }
            double obsMean = sumObs / n;
            double modMean = sumMod / n;

            // Second pass: compute standard deviations and correlation
            double varObs = 0d;
            double varMod = 0d;
            double covariance = 0d;

            for (int i = 0; i < n; i++)
            {
                double obsDeviation = observedValues[i] - obsMean;
                double modDeviation = modeledValues[i] - modMean;

                varObs += obsDeviation * obsDeviation;
                varMod += modDeviation * modDeviation;
                covariance += obsDeviation * modDeviation;
            }

            // Use sample standard deviation (N-1)
            double obsStdDev = Math.Sqrt(varObs / (n - 1));
            double modStdDev = Math.Sqrt(varMod / (n - 1));

            // Check for degenerate cases (zero variance)
            // When either standard deviation is zero, correlation is undefined
            // and the model fundamentally fails to capture variability
            if (obsStdDev < double.Epsilon || modStdDev < double.Epsilon)
            {
                // Return a very poor KGE value instead of NaN
                // This indicates the model is far worse than the mean baseline
                return -10.0;
            }

            // Correlation coefficient
            double r = covariance / (n - 1) / (obsStdDev * modStdDev);

            // Compute alpha (variability ratio)
            double alpha = modStdDev / obsStdDev;

            // Compute beta (bias ratio)
            double beta = modMean / obsMean;

            // Compute KGE
            double ed = Math.Sqrt(Tools.Sqr(r - 1.0) + Tools.Sqr(alpha - 1.0) + Tools.Sqr(beta - 1.0));
            return 1.0 - ed;
        }

        /// <summary>
        /// Computes the modified Kling-Gupta Efficiency (KGE') using the 2012 formulation.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The modified KGE' coefficient, or a large negative value if calculation is undefined.</returns>
        /// <remarks>
        /// <para>
        /// The 2012 modification replaces the variability ratio (α) with the variability ratio 
        /// based on coefficients of variation (γ, gamma):
        /// γ = (CV_modeled / CV_observed) = (σ_modeled/μ_modeled) / (σ_observed/μ_observed)
        /// </para>
        /// <para>
        /// <b>Special Cases:</b>
        /// <list type="bullet">
        /// <item><description>If modeled values have zero variance (all constant), returns -10.0</description></item>
        /// <item><description>If observed values have zero variance, returns -10.0</description></item>
        /// <item><description>These cases indicate the model cannot capture any variability</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// This modification provides better performance when dealing with biased flows and is particularly 
        /// useful for low-flow conditions. The interpretation of KGE' is the same as KGE.
        /// Reference: Kling et al. (2012)
        /// </para>
        /// </remarks>
        public static double KlingGuptaEfficiencyMod(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;

            // First pass: compute means
            double sumObs = 0d;
            double sumMod = 0d;
            for (int i = 0; i < n; i++)
            {
                sumObs += observedValues[i];
                sumMod += modeledValues[i];
            }
            double obsMean = sumObs / n;
            double modMean = sumMod / n;

            // Second pass: compute standard deviations and correlation
            double varObs = 0d;
            double varMod = 0d;
            double covariance = 0d;

            for (int i = 0; i < n; i++)
            {
                double obsDeviation = observedValues[i] - obsMean;
                double modDeviation = modeledValues[i] - modMean;

                varObs += obsDeviation * obsDeviation;
                varMod += modDeviation * modDeviation;
                covariance += obsDeviation * modDeviation;
            }

            // Use sample standard deviation (N-1)
            double obsStdDev = Math.Sqrt(varObs / (n - 1));
            double modStdDev = Math.Sqrt(varMod / (n - 1));

            // Check for degenerate cases (zero variance)
            if (obsStdDev < double.Epsilon || modStdDev < double.Epsilon)
            {
                // Return a very poor KGE' value instead of NaN
                return -10.0;
            }

            // Correlation coefficient
            double r = covariance / (n - 1) / (obsStdDev * modStdDev);

            // Compute gamma (coefficient of variation ratio)
            double cvObs = obsStdDev / obsMean;
            double cvMod = modStdDev / modMean;
            double gamma = cvMod / cvObs;

            // Compute beta (bias ratio)
            double beta = modMean / obsMean;

            // Compute KGE'
            double ed = Math.Sqrt(Tools.Sqr(r - 1.0) + Tools.Sqr(gamma - 1.0) + Tools.Sqr(beta - 1.0));
            return 1.0 - ed;
        }

        #endregion

        #region "Bias Metrics"

        /// <summary>
        /// Computes the Percent Bias (PBIAS) of the model compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The PBIAS as a percentage.</returns>
        /// <remarks>
        /// <para>
        /// PBIAS measures the average tendency of modeled values to be larger or smaller than observed values. 
        /// It is expressed as a percentage:
        /// <list type="bullet">
        /// <item><description>PBIAS = 0: Perfect model with no bias</description></item>
        /// <item><description>PBIAS &gt; 0: Model underestimates (negative bias)</description></item>
        /// <item><description>PBIAS &lt; 0: Model overestimates (positive bias)</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Performance ratings for streamflow (Moriasi et al., 2007):
        /// <list type="bullet">
        /// <item><description>Very good: PBIAS &lt; ±10%</description></item>
        /// <item><description>Good: ±10% ≤ PBIAS &lt; ±15%</description></item>
        /// <item><description>Satisfactory: ±15% ≤ PBIAS &lt; ±25%</description></item>
        /// <item><description>Unsatisfactory: PBIAS ≥ ±25%</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Note: Performance thresholds vary by constituent (sediment, nutrients have more lenient criteria).
        /// </para>
        /// </remarks>
        public static double PBIAS(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sumDiff = 0d;
            double sumObs = 0d;

            for (int i = 0; i < n; i++)
            {
                sumDiff += modeledValues[i] - observedValues[i];
                sumObs += observedValues[i];
            }

            return 100.0 * sumDiff / sumObs;
        }

        /// <summary>
        /// Computes the Ratio of RMSE to Standard Deviation of observations (RSR).
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The RSR value.</returns>
        /// <remarks>
        /// <para>
        /// RSR standardizes RMSE using the standard deviation of observations, making it a dimensionless metric. 
        /// It combines the benefits of error index statistics and incorporates a normalization factor:
        /// <list type="bullet">
        /// <item><description>RSR = 0: Perfect model (zero RMSE)</description></item>
        /// <item><description>RSR = 1: Model performance equals the standard deviation of observations</description></item>
        /// <item><description>RSR &gt; 1: Poor model performance</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Performance ratings (Moriasi et al., 2007):
        /// <list type="bullet">
        /// <item><description>Very good: 0.00 ≤ RSR ≤ 0.50</description></item>
        /// <item><description>Good: 0.50 &lt; RSR ≤ 0.60</description></item>
        /// <item><description>Satisfactory: 0.60 &lt; RSR ≤ 0.70</description></item>
        /// <item><description>Unsatisfactory: RSR &gt; 0.70</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Lower RSR values indicate better model performance.
        /// </para>
        /// <para>
        /// <b>Note:</b> This method uses sample standard deviation (N-1 denominator) for the observed values,
        /// which is consistent with the R hydroGOF package implementation.
        /// </para>
        /// </remarks>
        public static double RSR(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;

            // Single pass to compute RMSE and mean of observed
            double sse = 0d;
            double sumObs = 0d;
            for (int i = 0; i < n; i++)
            {
                sse += Tools.Sqr(modeledValues[i] - observedValues[i]);
                sumObs += observedValues[i];
            }

            double rmse = Math.Sqrt(sse / n);
            double obsMean = sumObs / n;

            // Second pass to compute sample standard deviation (using N-1)
            double variance = 0d;
            for (int i = 0; i < n; i++)
            {
                variance += Tools.Sqr(observedValues[i] - obsMean);
            }
            double stdDev = Math.Sqrt(variance / n);

            return rmse / stdDev;
        }

        #endregion

        #region "Correlation and Determination"

        /// <summary>
        /// Computes the coefficient of determination (R²), which is the square of the Pearson correlation coefficient.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The R² value.</returns>
        /// <remarks>
        /// <para>
        /// R² represents the proportion of variance in the observed data that is predictable from the modeled data.
        /// <list type="bullet">
        /// <item><description>R² = 1: Perfect correlation</description></item>
        /// <item><description>R² = 0: No linear relationship</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Note: R² does not indicate whether the model is biased or whether predictions are systematically 
        /// over or under-predicted. It should be used in combination with other metrics like PBIAS or visual 
        /// inspection of residuals. R² is equivalent to NSE when used in the context of comparing predictions 
        /// to observations in a regression setting.
        /// </para>
        /// </remarks>
        public static double RSquared(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            var corr = Correlation.Pearson(observedValues, modeledValues);
            return corr * corr;
        }

        #endregion

        #region "Index of Agreement Metrics"

        /// <summary>
        /// Computes the Index of Agreement (d), also known as Willmott's Index of Agreement.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The Index of Agreement (d).</returns>
        /// <remarks>
        /// <para>
        /// The Index of Agreement (d) was developed by Willmott (1981) as a standardized measure of 
        /// the degree of model prediction error.
        /// </para>
        /// <para>
        /// Range: 0 to 1
        /// <list type="bullet">
        /// <item><description>d = 1: Perfect agreement</description></item>
        /// <item><description>d = 0: No agreement</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Formula: d = 1 - [Σ(M - O)²] / [Σ(|M - Ō| + |O - Ō|)²]
        /// </para>
        /// <para>
        /// Advantages:
        /// <list type="bullet">
        /// <item><description>Can detect additive and proportional differences in means and variances</description></item>
        /// <item><description>Bounded and dimensionless</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Limitations:
        /// <list type="bullet">
        /// <item><description>Overly sensitive to extreme values due to squared differences</description></item>
        /// <item><description>Can be insensitive to proportional differences</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// References: 
        /// Willmott, C.J. (1981). On the validation of models. Physical Geography, 2, 184-194.
        /// </para>
        /// </remarks>
        public static double IndexOfAgreement(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double observedMean = Statistics.Mean(observedValues);

            double numerator = 0d;
            double denominator = 0d;

            for (int i = 0; i < n; i++)
            {
                numerator += Tools.Sqr(modeledValues[i] - observedValues[i]);
                denominator += Tools.Sqr(Math.Abs(modeledValues[i] - observedMean) + Math.Abs(observedValues[i] - observedMean));
            }

            return 1.0 - (numerator / denominator);
        }

        /// <summary>
        /// Computes the Modified Index of Agreement (d1), which uses absolute values instead of squares.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The Modified Index of Agreement (d1).</returns>
        /// <remarks>
        /// <para>
        /// The modified index (d1) was proposed by Willmott et al. (1985) to address the 
        /// over-sensitivity to extreme values in the original index.
        /// </para>
        /// <para>
        /// Formula: d1 = 1 - [Σ|M - O|] / [Σ(|M - Ō| + |O - Ō|)]
        /// </para>
        /// <para>
        /// Range: 0 to 1 (same as original d)
        /// </para>
        /// <para>
        /// Advantages over original d:
        /// <list type="bullet">
        /// <item><description>Less sensitive to outliers (uses absolute differences instead of squares)</description></item>
        /// <item><description>More robust evaluation of model performance</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// References:
        /// Willmott, C.J., Ackleson, S.G., Davis, R.E., et al. (1985). Statistics for the evaluation 
        /// and comparison of models. Journal of Geophysical Research, 90(C5), 8995-9005.
        /// </para>
        /// </remarks>
        public static double ModifiedIndexOfAgreement(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double observedMean = Statistics.Mean(observedValues);

            double numerator = 0d;
            double denominator = 0d;

            for (int i = 0; i < n; i++)
            {
                numerator += Math.Abs(modeledValues[i] - observedValues[i]);
                denominator += Math.Abs(modeledValues[i] - observedMean) + Math.Abs(observedValues[i] - observedMean);
            }

            return 1.0 - (numerator / denominator);
        }

        /// <summary>
        /// Computes the Refined Index of Agreement (dr), which has an extended range from -1 to 1.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The Refined Index of Agreement (dr).</returns>
        /// <remarks>
        /// <para>
        /// The refined index (dr) was proposed by Willmott et al. (2012) to provide a broader range 
        /// of values that better represents the variety of forms that predictions can differ from observations.
        /// </para>
        /// <para>
        /// Range: -1 to 1
        /// <list type="bullet">
        /// <item><description>dr = 1: Perfect agreement</description></item>
        /// <item><description>dr = 0: Model performs as well as the mean</description></item>
        /// <item><description>dr &lt; 0: Model performs worse than the mean</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Formula depends on the ratio of Σ|M - O| to 2×Σ|O - Ō|:
        /// <list type="bullet">
        /// <item><description>When Σ|M - O| ≤ 2×Σ|O - Ō|: dr = 1 - [Σ|M - O|] / [2×Σ|O - Ō|]</description></item>
        /// <item><description>When Σ|M - O| > 2×Σ|O - Ō|: dr = [2×Σ|O - Ō|] / [Σ|M - O|] - 1</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Advantages:
        /// <list type="bullet">
        /// <item><description>More rationally related to model accuracy than d or d1</description></item>
        /// <item><description>Can indicate when model is worse than using the mean (negative values)</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Limitation: Does not indicate direction of bias (over- or under-estimation)
        /// </para>
        /// <para>
        /// References:
        /// Willmott, C.J., Robeson, S.M., Matsuura, K. (2012). A refined index of model performance. 
        /// International Journal of Climatology, 32(13), 2088-2094.
        /// </para>
        /// </remarks>
        public static double RefinedIndexOfAgreement(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double observedMean = Statistics.Mean(observedValues);

            double sumAbsError = 0d;
            double sumAbsDeviation = 0d;

            for (int i = 0; i < n; i++)
            {
                sumAbsError += Math.Abs(modeledValues[i] - observedValues[i]);
                sumAbsDeviation += Math.Abs(observedValues[i] - observedMean);
            }

            double c = 2.0 * sumAbsDeviation;

            if (sumAbsError <= c)
            {
                return 1.0 - (sumAbsError / c);
            }
            else
            {
                return (c / sumAbsError) - 1.0;
            }
        }

        /// <summary>
        /// Computes the Volumetric Efficiency (VE), which measures the agreement between cumulative 
        /// volumes of observed and modeled values.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The Volumetric Efficiency (VE).</returns>
        /// <remarks>
        /// <para>
        /// VE was proposed by Criss and Winston (2008) as an alternative to NSE. It measures the 
        /// fractional error in the model's volumetric (or total mass) estimates.
        /// </para>
        /// <para>
        /// Range: -∞ to 1
        /// <list type="bullet">
        /// <item><description>VE = 1: Perfect volumetric match</description></item>
        /// <item><description>VE = 0: Model total equals mean of observations</description></item>
        /// <item><description>VE &lt; 0: Model performs worse than using the mean</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Formula: VE = 1 - [Σ|M - O|] / [Σ|O|]
        /// </para>
        /// <para>
        /// Advantages:
        /// <list type="bullet">
        /// <item><description>Directly interpretable in terms of volume or mass errors</description></item>
        /// <item><description>Less sensitive to outliers than NSE (uses absolute values)</description></item>
        /// <item><description>Focuses on overall water balance, important for water resources management</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Particularly useful for hydrological models where accurate water balance is critical.
        /// </para>
        /// <para>
        /// References:
        /// Criss, R.E., Winston, W.E. (2008). Do Nash values have value? Discussion and alternate proposals. 
        /// Hydrological Processes, 22(14), 2723-2725.
        /// </para>
        /// </remarks>
        public static double VolumetricEfficiency(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            int n = observedValues.Count;
            double sumAbsError = 0d;
            double sumAbsObserved = 0d;

            for (int i = 0; i < n; i++)
            {
                sumAbsError += Math.Abs(modeledValues[i] - observedValues[i]);
                sumAbsObserved += Math.Abs(observedValues[i]);
            }

            // Check for zero denominator
            if (Math.Abs(sumAbsObserved) < double.Epsilon)
                throw new ArgumentException("VE cannot be calculated when sum of absolute observed values is zero.", nameof(observedValues));

            return 1.0 - (sumAbsError / sumAbsObserved);
        }

        #endregion

        #region "Classification Metrics"

        /// <summary>
        /// Computes the accuracy of the modeled data compared to the observed data.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against.</param>
        /// <param name="modeledValues">The list of modeled values to compare against the observed values.</param>
        /// <returns>The accuracy as a percentage.</returns>
        /// <remarks>
        /// Accuracy is the proportion of exact matches between observed and modeled values. 
        /// This metric is primarily useful for classification problems or categorical data. 
        /// For continuous data, exact matches are rare, making this metric less meaningful.
        /// </remarks>
        public static double Accuracy(IList<double> observedValues, IList<double> modeledValues)
        {
            // Check if the lists are the same size
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues), "The number of observed values must equal the number of modeled values.");

            double accuracy = 0;
            for (int i = 0; i < observedValues.Count; ++i)
            {
                if (observedValues[i] == modeledValues[i])
                {
                    accuracy++;
                }
            }
            return 100 * accuracy / observedValues.Count;
        }


        /// <summary>
        /// Computes the confusion matrix components for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>A tuple containing (TP, TN, FP, FN).</returns>
        private static (int TP, int TN, int FP, int FN) GetConfusionMatrixComponents(
            IList<double> observedValues, IList<double> modeledValues)
        {
            if (observedValues.Count != modeledValues.Count)
                throw new ArgumentOutOfRangeException(nameof(observedValues),
                    "The number of observed values must equal the number of modeled values.");

            int TP = 0, TN = 0, FP = 0, FN = 0;

            for (int i = 0; i < observedValues.Count; i++)
            {
                bool observed = Math.Abs(observedValues[i] - 1.0) < 0.01;
                bool modeled = Math.Abs(modeledValues[i] - 1.0) < 0.01;

                if (observed && modeled) TP++;
                else if (!observed && !modeled) TN++;
                else if (!observed && modeled) FP++;
                else if (observed && !modeled) FN++;
            }

            return (TP, TN, FP, FN);
        }

        /// <summary>
        /// Computes the Precision (Positive Predictive Value) for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>The Precision value (0 to 1).</returns>
        /// <remarks>
        /// <para>
        /// Precision measures the proportion of positive predictions that are actually correct.
        /// Formula: Precision = TP / (TP + FP)
        /// </para>
        /// <para>
        /// Range: 0 to 1 (1 = perfect precision)
        /// </para>
        /// <para>
        /// Interpretation:
        /// <list type="bullet">
        /// <item><description>Precision = 1: All positive predictions are correct (no false positives)</description></item>
        /// <item><description>High precision: Model is conservative in making positive predictions</description></item>
        /// <item><description>Low precision: Many false alarms</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Use case: Prioritize when false positives are costly (e.g., spam detection, medical testing for benign conditions)
        /// </para>
        /// </remarks>
        public static double Precision(IList<double> observedValues, IList<double> modeledValues)
        {
            var (TP, TN, FP, FN) = GetConfusionMatrixComponents(observedValues, modeledValues);

            if (TP + FP == 0)
                throw new ArgumentException("Cannot calculate Precision when there are no positive predictions.");

            return (double)TP / (TP + FP);
        }

        /// <summary>
        /// Computes the Recall (Sensitivity, True Positive Rate) for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>The Recall value (0 to 1).</returns>
        /// <remarks>
        /// <para>
        /// Recall measures the proportion of actual positives that are correctly identified.
        /// Formula: Recall = TP / (TP + FN)
        /// </para>
        /// <para>
        /// Range: 0 to 1 (1 = perfect recall)
        /// </para>
        /// <para>
        /// Interpretation:
        /// <list type="bullet">
        /// <item><description>Recall = 1: All actual positives are identified (no false negatives)</description></item>
        /// <item><description>High recall: Model captures most positive cases</description></item>
        /// <item><description>Low recall: Many missed detections</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Use case: Prioritize when false negatives are costly (e.g., disease detection, fraud detection)
        /// </para>
        /// </remarks>
        public static double Recall(IList<double> observedValues, IList<double> modeledValues)
        {
            var (TP, TN, FP, FN) = GetConfusionMatrixComponents(observedValues, modeledValues);

            if (TP + FN == 0)
                throw new ArgumentException("Cannot calculate Recall when there are no actual positive cases.");

            return (double)TP / (TP + FN);
        }

        /// <summary>
        /// Computes the F1 Score (harmonic mean of Precision and Recall) for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>The F1 Score (0 to 1).</returns>
        /// <remarks>
        /// <para>
        /// F1 Score is the harmonic mean of Precision and Recall, providing a balanced measure 
        /// that considers both false positives and false negatives.
        /// </para>
        /// <para>
        /// Formula: F1 = 2 × (Precision × Recall) / (Precision + Recall) = 2TP / (2TP + FP + FN)
        /// </para>
        /// <para>
        /// Range: 0 to 1 (1 = perfect F1)
        /// </para>
        /// <para>
        /// Characteristics:
        /// <list type="bullet">
        /// <item><description>F1 = 1 only when both Precision and Recall are 1</description></item>
        /// <item><description>F1 approaches the lower of Precision and Recall when they differ greatly</description></item>
        /// <item><description>Particularly useful for imbalanced datasets</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Interpretation:
        /// <list type="bullet">
        /// <item><description>F1 &gt; 0.9: Excellent model</description></item>
        /// <item><description>F1 0.7-0.9: Good model</description></item>
        /// <item><description>F1 0.5-0.7: Moderate model</description></item>
        /// <item><description>F1 &lt; 0.5: Poor model</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Use case: When you need balanced performance with both types of errors equally important
        /// </para>
        /// </remarks>
        public static double F1Score(IList<double> observedValues, IList<double> modeledValues)
        {
            var (TP, TN, FP, FN) = GetConfusionMatrixComponents(observedValues, modeledValues);

            if (2 * TP + FP + FN == 0)
                throw new ArgumentException("Cannot calculate F1 Score when there are no positive predictions or actual positives.");

            return (2.0 * TP) / (2 * TP + FP + FN);
        }

        /// <summary>
        /// Computes the Specificity (True Negative Rate) for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>The Specificity value (0 to 1).</returns>
        /// <remarks>
        /// <para>
        /// Specificity measures the proportion of actual negatives that are correctly identified.
        /// Formula: Specificity = TN / (TN + FP)
        /// </para>
        /// <para>
        /// Range: 0 to 1 (1 = perfect specificity)
        /// </para>
        /// <para>
        /// Interpretation:
        /// <list type="bullet">
        /// <item><description>Specificity = 1: All actual negatives are correctly identified</description></item>
        /// <item><description>High specificity: Few false positives</description></item>
        /// <item><description>Low specificity: Many false alarms</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Complements Recall (Sensitivity). Together they provide complete picture of model performance.
        /// </para>
        /// </remarks>
        public static double Specificity(IList<double> observedValues, IList<double> modeledValues)
        {
            var (TP, TN, FP, FN) = GetConfusionMatrixComponents(observedValues, modeledValues);

            if (TN + FP == 0)
                throw new ArgumentException("Cannot calculate Specificity when there are no actual negative cases.");

            return (double)TN / (TN + FP);
        }

        /// <summary>
        /// Computes the Balanced Accuracy for binary classification.
        /// </summary>
        /// <param name="observedValues">The list of observed binary values (0 or 1).</param>
        /// <param name="modeledValues">The list of modeled binary values (0 or 1).</param>
        /// <returns>The Balanced Accuracy value (0 to 1).</returns>
        /// <remarks>
        /// <para>
        /// Balanced Accuracy is the average of Recall (Sensitivity) and Specificity.
        /// Formula: Balanced Accuracy = (Recall + Specificity) / 2
        /// </para>
        /// <para>
        /// Range: 0 to 1 (1 = perfect balanced accuracy)
        /// </para>
        /// <para>
        /// Advantages:
        /// <list type="bullet">
        /// <item><description>Accounts for class imbalance by giving equal weight to both classes</description></item>
        /// <item><description>More appropriate than regular accuracy for imbalanced datasets</description></item>
        /// <item><description>Random guessing gives 0.5 balanced accuracy regardless of class distribution</description></item>
        /// </list>
        /// </para>
        /// <para>
        /// Use case: Essential for evaluating models on imbalanced datasets
        /// </para>
        /// </remarks>
        public static double BalancedAccuracy(IList<double> observedValues, IList<double> modeledValues)
        {
            double recall = Recall(observedValues, modeledValues);
            double specificity = Specificity(observedValues, modeledValues);

            return (recall + specificity) / 2.0;
        }

        #endregion

        #region "Statistical Tests"

        /// <summary>
        /// Computes the Kolmogorov-Smirnov (K-S) test statistic, which measures the maximum vertical distance 
        /// between the empirical cumulative distribution function (ECDF) of the sample and the theoretical CDF.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution to test against.</param>
        /// <returns>The Kolmogorov-Smirnov test statistic (D).</returns>
        /// <remarks>
        /// <para>
        /// The K-S test statistic D measures the maximum absolute difference between the ECDF and the theoretical CDF. 
        /// Smaller values indicate better fit.
        /// </para>
        /// <para>
        /// The K-S test is distribution-free and is sensitive to both location and shape of the distribution. 
        /// However, it is generally less powerful than the Anderson-Darling test, especially in the tails.
        /// </para>
        /// <para>
        /// <b>Reference:</b> 
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm"/>
        /// </para>
        /// </remarks>
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
        /// Computes the Chi-Squared (χ²) test statistic for goodness-of-fit, which measures how well 
        /// the observed frequencies match the expected frequencies from the theoretical distribution.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution to test against.</param>
        /// <returns>The Chi-Squared test statistic (χ²).</returns>
        /// <remarks>
        /// <para>
        /// The Chi-Squared test divides the data into bins and compares observed frequencies with expected 
        /// frequencies based on the model. Smaller values indicate better fit.
        /// </para>
        /// <para>
        /// This test is best suited for large sample sizes (typically n > 50) and requires that expected 
        /// frequencies in each bin are at least 5 for the chi-square approximation to be valid.
        /// </para>
        /// <para>
        /// <b>Reference:</b>
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm"/>
        /// </para>
        /// </remarks>
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
        /// Computes the Anderson-Darling (A-D) test statistic, which is a weighted goodness-of-fit measure 
        /// that gives more weight to discrepancies in the tails of the distribution.
        /// </summary>
        /// <param name="observedValues">The list of observed values to measure against. Must be sorted in ascending order.</param>
        /// <param name="model">The univariate continuous distribution to test against.</param>
        /// <returns>The Anderson-Darling test statistic (A²).</returns>
        /// <remarks>
        /// <para>
        /// The Anderson-Darling test is more sensitive to discrepancies in the tails of the distribution 
        /// compared to the Kolmogorov-Smirnov test. Smaller values indicate better fit.
        /// </para>
        /// <para>
        /// This test is particularly useful for hydrological applications where extreme values (floods, droughts) 
        /// are of primary interest. It is generally more powerful than the K-S test for detecting departures 
        /// from the specified distribution.
        /// </para>
        /// <para>
        /// <b>Reference:</b>
        /// <see href="https://www.itl.nist.gov/div898/handbook/eda/section3/eda35e.htm"/>
        /// </para>
        /// </remarks>
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

        #endregion

    }
}