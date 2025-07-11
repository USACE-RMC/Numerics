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

using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// A class for performing generalized linear regression. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// </remarks>
    public class GeneralizedLinearModel
    {

        public GeneralizedLinearModel(Matrix x, Vector y, bool hasIntercept = true) 
        {
            if (y.Length != x.NumberOfRows) throw new ArgumentException("X and Y must have the same number of rows.");
            if (y.Length <= 2) throw new ArithmeticException("There must be at least three data points.");
            if (x.NumberOfColumns > y.Length) throw new ArithmeticException($"A regression of the requested order requires at least {x.NumberOfColumns} data points. Only {y.Length} data points have been provided.");

            // Set inputs
            Y = y;
            X = hasIntercept ? AddInterceptColumn(x) : x;
            this.HasIntercept = hasIntercept;
            DegreesOfFreedom = y.Length - x.NumberOfRows;

            // Set model name
            if (Y.Header == null || Y.Header.Length == 0)
            {
                Y.Header = "Y Data";
            }

            // Set parameter names for summary report.
            ParameterNames = new List<string>();
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
        }

        /// <summary>
        /// Determines if the linear model has an intercept. 
        /// </summary>
        public bool HasIntercept { get; private set; }

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
        /// The residuals of the fitted linear model. 
        /// </summary>
        public double[] Residuals { get; private set; }

        /// <summary>
        /// The model standard error.
        /// </summary>
        public double StandardError { get; private set; }

        /// <summary>
        /// The data sample size. 
        /// </summary>
        public int SampleSize => Y.Length;

        /// <summary>
        /// The model degrees of freedom.
        /// </summary>
        public int DegreesOfFreedom { get; private set; }

        /// <summary>
        /// The Coefficient of Determination (or R-squared). 
        /// </summary>
        public double RSquared { get; private set; }

        /// <summary>
        /// the Akaike Information Criterion (AIC). 
        /// </summary>
        public double AIC { get; private set; }

        /// <summary>
        /// the Akaike Information Criterion corrected for small sample sizes (AICc)
        /// </summary>
        public double AICc { get; private set; }

        /// <summary>
        /// The Bayesian information criterion (BIC). 
        /// </summary>
        public double BIC { get; private set; }

        /// <summary>
        /// Estimate robust standard errors.
        /// </summary>
        public bool UseRobustSE { get; set; } = false;

        private Func<double, double> _inverseLink;
        private Func<double, double> _inverseLinkDerivative;
        private Func<(double mu, double y), double> _logLikelihoodTerm;
        private Func<(double eta, double y), double> _logLikelihoodGradient;

        /// <summary>
        /// Helper method to add an intercept column to the covariate matrix.
        /// </summary>
        /// <param name="x">The matrix of predictor values.</param>
        private static Matrix AddInterceptColumn(Matrix x)
        {
            Matrix result = new Matrix(x.NumberOfRows, x.NumberOfColumns + 1);
            for (int i = 0; i < x.NumberOfRows; i++)
            {
                result[i, 0] = 1.0;
                for (int j = 0; j < x.NumberOfColumns; j++)
                    result[i, j + 1] = x[i, j];
            }
            return result;
        }

    }
}
