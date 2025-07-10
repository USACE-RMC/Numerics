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

using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;

namespace Numerics.Data
{

    /// <summary>
    /// A class for performing linear regression. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// Linear regression works to estimate the linear relationship between a single, scalar response (dependent) variable and one or more explanatory 
    /// (independent) variables. This implementation estimates Y = α + βX + e, where e ~ N(0,σ), and fits the linear regression model using singular value decomposition (SVD).
    /// </para>
    /// <b> References: </b>
    /// <see href="https://en.wikipedia.org/wiki/Linear_regression"/>
    /// </remarks>
    public class LinearRegression
    {

        /// <summary>
        /// Estimates the model Y = α + βX + e, where e ~ N(0,σ).
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        /// <param name="y">The response vector.</param>
        /// <param name="hasIntercept">Determines if an intercept should be estimate. Default = true.</param>
        public LinearRegression(Matrix x, Vector y, bool hasIntercept = true)
        {

            if (y.Length != x.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (y.Length <= 2) throw new ArithmeticException("There must be at least three data points.");
            if (x.NumberOfColumns > y.Length) throw new ArithmeticException($"A regression of the requested order requires at least {x.NumberOfColumns} data points. Only {y.Length} data points have been provided.");

            // Set inputs
            Y = y;
            X = x;
            this.hasIntercept = hasIntercept;
            ParameterNames = new List<string>();

            // Set model name
            if (Y.Header == null || Y.Header.Length == 0)
            {
                Y.Header = "Y Data";
            }

            // Set parameter names for summary report.
            if (hasIntercept) ParameterNames.Add("Intercept");
            if (X.Header != null && X.Header.Length == x.NumberOfColumns)
            {
                ParameterNames.AddRange(X.Header);
            }
            else
            {
                for (int i = 1; i <= x.NumberOfColumns; i++)
                {
                    ParameterNames.Add("β" + i);
                }
            }

            // Estimate the linear model.
            FitSVD();
        }

        private bool hasIntercept;

        /// <summary>
        /// The vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// The list of estimated parameter values.
        /// </summary>
        public List<double> Parameters { get; private set; }

        /// <summary>
        /// The list of the estimated parameter names. 
        /// </summary>
        public List<string> ParameterNames { get; private set; }

        /// <summary>
        /// The list of the estimated parameter standard errors. 
        /// </summary>
        public List<double> ParameterStandardErrors { get; private set; }

        /// <summary>
        /// The list of the estimated parameter t-statistics.
        /// </summary>
        public List<double> ParameterTStats { get; private set; }

        /// <summary>
        /// The estimate parameter covariance matrix. 
        /// </summary>
        public Matrix Covariance { get; private set; }

        /// <summary>
        /// The model standard error.
        /// </summary>
        public double StandardError { get; private set; }

        /// <summary>
        /// The data sample size. 
        /// </summary>
        public int SampleSize
        {
            get { return Y.Length; }
        }

        /// <summary>
        /// The model degrees of freedom.
        /// </summary>
        public int DegreesOfFreedom { get; private set; }

        /// <summary>
        /// The Coefficient of Determination (or R-squared). 
        /// </summary>
        public double RSquared { get; private set; }

        /// <summary>
        /// Adjusted R-squared.
        /// </summary>
        public double AdjRSquared { get; private set; }

        /// <summary>
        /// The residuals of the fitted linear model. 
        /// </summary>
        public double[] Residuals { get; private set; }

        /// <summary>
        /// Provides a standard summary output table in a list of strings. 
        /// </summary>
        public List<string> Summary()
        {

            // This summary output is designed to look just like the summary output from R. 
            var text = new List<string>();

            // Model parameters
            text.Add("");
            text.Add("Model for predicting " + Y.Header + ":");
            text.Add("Parameters:");
            text.Add($"{"", -15}" + $"{"Estimate",12}" + $"{"Std. Error", 12}" + $"{"t value", 12}" + $"{"Pr(>|t|)", 12}");
            for (int i = 0; i < Parameters.Count; i++)
            {
                var name = ParameterNames[i].Length < 15 ? ParameterNames[i] : ParameterNames[i].Substring(0, 14);
                // Get T-stats
                var tdist = new StudentT(DegreesOfFreedom) ;
                double tval = Parameters[i] / ParameterStandardErrors[i];
                double pval = (1 - tdist.CDF(Math.Abs(tval))) * 2;
                string pvalStrg = pval > 1E-4 ? pval.ToString("N4") : pval < 1E-15 ? "< 1E-15" : pval.ToString("E2");
                string sig =  pval < 1E-3 ? " ***" : pval < 1E-2 ? " **" : pval < 0.05 ? " *" : pval < 0.1 ? " ." : "  ";

                text.Add($"{name, -15}" + $"{Parameters[i].ToString("N5"), 12}" + $"{ParameterStandardErrors[i].ToString("N5"), 12}" +
                         $"{tval.ToString("N3"), 12}" + $"{pvalStrg, 12}" + $"{sig, -2}");
            }
            text.Add("---");
            text.Add("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1");
            text.Add("");

            // Goodness of fit
            text.Add("Residual standard error: " + StandardError.ToString("N4") + " on " + DegreesOfFreedom.ToString("N0") + " degrees of freedom");
            text.Add("Multiple R-squared: " + RSquared.ToString("N4") + ",  Adjusted R-squared: " + AdjRSquared.ToString("N4"));

            // F-Test
            // Estimate a model with just an intercept
            var lm = new LinearRegression(new Matrix(new Vector(Y.Length, 1)), Y, false);
            double F = 0, fP = 0, sseR = lm.StandardError * lm.StandardError * lm.DegreesOfFreedom, sseF = StandardError * StandardError * DegreesOfFreedom;
            int dfR = lm.DegreesOfFreedom, dfF = DegreesOfFreedom;
            HypothesisTests.FtestModels(sseR, sseF, dfR, dfF, out F, out fP);
            string fPStrg = fP > 1E-4 ? fP.ToString("N4") : fP < 1E-15 ? "< 1E-15" : fP.ToString("E2");
            text.Add("F-statistic: " + F.ToString("N1") + ", on " + (hasIntercept ? Parameters.Count - 1 : Parameters.Count).ToString("N0") + " and " + DegreesOfFreedom.ToString("N0") + " DF, p-value: " + fPStrg);
            text.Add("");

            // Residuals
            text.Add("Residuals:");
            text.Add($"{"Min",10}" + $"{"1Q",10}" + $"{"Median",10}" + $"{"3Q",10}" + $"{"Max",10}");
            var res = Statistics.Statistics.FiveNumberSummary(Residuals);
            text.Add($"{res[0].ToString("N4"),10}" + $"{res[1].ToString("N4"),10}" + $"{res[2].ToString("N4"),10}" + $"{res[3].ToString("N4"),10}" + $"{res[4].ToString("N4"),10}");
            text.Add("");

            return text;
        }

        /// <summary>
        /// Estimate the model using Singular Value Decomposition. 
        /// </summary>
        private void FitSVD()
        {
            double meanY = Statistics.Statistics.Mean(Y.ToArray());
            var xx = X;
            // see if we need to add an intercept column
            if (hasIntercept)
            {
                xx = new Matrix(X.NumberOfRows, X.NumberOfColumns + 1);
                for (int r = 0; r < X.NumberOfRows; r++)
                {
                    xx[r, 0] = 1;
                    for (int c = 0; c < X.NumberOfColumns; c++)
                    {
                        xx[r, c + 1] = X[r, c];
                    }
                }
            }

            int i, j, k, n = xx.NumberOfRows, m = xx.NumberOfColumns;
            double sse = 0.0, sst = 0.0, sum = 0.0;
            DegreesOfFreedom = n - m;
            Residuals = new double[n];
            Covariance = new Matrix(m);

            // Estimate coefficients
            var svd = new SingularValueDecomposition(xx);
            double tol = 1E-12;
            double thresh = (tol > 0d ? tol * svd.W[0] : -1d);
            var betas = svd.Solve(Y, thresh); // vector of fitted coefficients
       
            // Estimate uncertainty in the coefficients
            for (i = 0; i < n; i++)
            {
                sum = 0.0;
                for (j = 0; j < m; j++) sum += xx[i, j] * betas[j];
                Residuals[i] = Y[i] - sum;
                sse += Tools.Sqr(Residuals[i]);
                sst += Tools.Sqr(Y[i] - meanY);
            }
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < i + 1; j++)
                {
                    sum = 0.0;
                    for (k = 0; k < m; k++) if (svd.W[k] > svd.Threshold)
                            sum += svd.V[i, k] * svd.V[j, k] / Tools.Sqr(svd.W[k]);
                    Covariance[j, i] = Covariance[i, j] = sum;
                }
            }
            double SE = Math.Sqrt(sse / DegreesOfFreedom);

            // Set the output
            Parameters = betas.ToList();
            ParameterStandardErrors = new List<double>();
            for (j = 0; j < m; j++) ParameterStandardErrors.Add(Math.Sqrt(Covariance[j, j]) * SE);
            StandardError = SE;
            RSquared = 1 - sse / sst;
            AdjRSquared = 1 - (1 - RSquared) * (n - 1) / (DegreesOfFreedom);

        }


        /// <summary>
        /// Returns the predicted Y values given the X-values. 
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        public double[] Predict(Matrix x)
        {
            var result = new double[x.NumberOfRows];
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                if (hasIntercept == true)
                {
                    var values = new List<double>() { 1 };
                    values.AddRange(x.Row(i));
                    result[i] = Tools.SumProduct(Parameters, values);
                }
                else
                {
                    result[i] = Tools.SumProduct(Parameters, x.Row(i));
                }
            }
                
            return result;
        }

        /// <summary>
        /// Returns the prediction intervals for Y in a 2D array with columns: lower, upper, mean.
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        /// <param name="alpha">The confidence level; Default = 0.1, which will result in the 90% prediction intervals.</param>
        public double[,] PredictionIntervals(Matrix x, double alpha = 0.1)
        {
            var percentiles = new double[] { alpha / 2d, 1d - alpha / 2d };
            var result = new double[x.NumberOfRows, 3]; // lower, upper, mean
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                double mu = 0;
                if (hasIntercept == true)
                {
                    var values = new List<double>() { 1 };
                    values.AddRange(x.Row(i));
                    mu = Tools.SumProduct(Parameters, values);

                }
                else
                {
                    mu = Tools.SumProduct(Parameters, x.Row(i));
                }
                var s = StandardError;
                var s2 = s * s;
                var n = DegreesOfFreedom;
                var t = new StudentT(mu, Math.Sqrt(s2 / n + s2), n);
                result[i, 0] = t.InverseCDF(percentiles[0]);
                result[i, 1] = t.InverseCDF(percentiles[1]);
                result[i, 2] = mu;
            }

            return result;
        }

    }
}