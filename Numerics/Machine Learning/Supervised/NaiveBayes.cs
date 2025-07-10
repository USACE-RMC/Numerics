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
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Linq;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// A class for Gaussian Naive Bayes Classification.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The Naive Bayes classifier is a "probabilistic classifier" based on 
    /// applying Bayes' theorem with strong (naive) independence assumption between features.
    /// </para>
    /// <para>
    /// The class implements Gaussian Naive Bayes, which assumes the features associated with each class
    /// are continuous and distributed according to a Normal (Gaussian) distribution. 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// <list type="bullet">
    /// <item> <see href="https://en.wikipedia.org/wiki/Naive_Bayes_classifier"/></item>
    /// </list>
    /// </para>
    /// </remarks>
    public class NaiveBayes
    {

        #region Construction

        /// <summary>
        /// Create new Naive Bayes classifier.
        /// </summary>
        /// <param name="x">The training 1D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        public NaiveBayes(double[] x, double[] y)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Classes = y.Distinct().ToArray();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
            if (Classes.Length < 1) throw new ArgumentException("There must be at least 1 class to predict.");

        }

        /// <summary>
        /// Create new Naive Bayes classifier.
        /// </summary>
        /// <param name="x">The training 2D array of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        public NaiveBayes(double[,] x, double[] y)
        {
            // Set inputs
            Y = new Vector(y);
            X = new Matrix(x);
            Classes = y.Distinct().ToArray();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
            if (Classes.Length < 1) throw new ArgumentException("There must be at least 1 class to predict.");
        }

        /// <summary>
        /// Create new Naive Bayes classifier.
        /// </summary>
        /// <param name="x">The training matrix of predictor values.</param>
        /// <param name="y">The training response vector.</param>
        public NaiveBayes(Matrix x, Vector y)
        {
            // Set inputs
            Y = y;
            X = x;
            Classes = y.Array.Distinct().ToArray();

            if (Y.Length != X.NumberOfRows) throw new ArgumentException("The y vector must be the same length as the x matrix.");
            if (Y.Length < 10) throw new ArgumentException("There must be at least ten training data points.");
            if (Classes.Length < 1) throw new ArgumentException("There must be at least 1 class to predict.");
        }

        #endregion

        #region Members

        /// <summary>
        /// The training vector of response values.
        /// </summary>
        public Vector Y { get; private set; }

        /// <summary>
        /// The training matrix of predictor values. 
        /// </summary>
        public Matrix X { get; private set; }

        /// <summary>
        /// Returns the list of distinct classes from the training set. 
        /// </summary>
        public double[] Classes { get; private set; }

        /// <summary>
        /// The means of each feature given each class. 
        /// </summary>
        public double[,] Means { get; private set; }

        /// <summary>
        /// The standard deviations of each feature given each class.
        /// </summary>
        public double[,] StandardDeviations { get; private set; }

        /// <summary>
        /// The prior probability for each class.
        /// </summary>
        public double[] Priors { get; private set; }

        /// <summary>
        /// Determines if the classifier has been trained.
        /// </summary>
        public bool IsTrained { get; private set; } = false;

        #endregion

        #region Methods

        /// <summary>
        /// Train the Naive Bayes classifier. 
        /// </summary>
        public void Train()
        {
            // Set up training outputs
            IsTrained = false;
            int nSamples = X.NumberOfRows;
            int nFeatures = X.NumberOfColumns;
            int nClasses = Classes.Length;
            Means = new double[nClasses, nFeatures];
            StandardDeviations = new double[nClasses, nFeatures];
            Priors = new double[nClasses];

            for (int i = 0; i < nClasses; i++)
            {
                // Compute priors as relative frequency of classes
                for (int k = 0; k < nSamples; k++)
                {
                    if (Y[k] == Classes[i])
                    {
                        Priors[i]++;
                    }
                }
                Priors[i] /= nSamples;

                // Compute the mean and standard deviation of each feature j given the class i
                for (int j = 0; j < nFeatures; j++)
                {
                    double x = 0;     // sum
                    double x2 = 0;    // sum of X^2
                    double u1, u2;
                    double n = 0;
                    // Compute sums
                    for (int k = 0; k < nSamples; k++)
                    {
                        if (Y[k] == Classes[i])
                        {
                            x += X[k, j];
                            x2 += Math.Pow(X[k, j], 2);
                            n++;
                        }
                    }
                    // Compute averages
                    u1 = x / n;
                    u2 = x2 / n;
                    // Set means
                    Means[i, j] = u1;
                    // Set standard deviations
                    StandardDeviations[i, j] = Math.Sqrt((u2 - Math.Pow(u1, 2d)) * (n / (n - 1)));
                }

            }

            IsTrained = true;
        }

        /// <summary>
        /// Returns the prediction of the Naive Bayes classifier
        /// </summary>
        /// <param name="X">The 1D array of predictors.</param>
        public double[] Predict(double[] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction of the Naive Bayes classifier
        /// </summary>
        /// <param name="X">The 2D array of predictors.</param>
        public double[] Predict(double[,] X)
        {
            return Predict(new Matrix(X));
        }

        /// <summary>
        /// Returns the prediction of the Naive Bayes classifier.
        /// </summary>
        /// <param name="X">The matrix of predictors.</param>
        public double[] Predict(Matrix X)
        {
            if (!IsTrained || X.NumberOfColumns != this.X.NumberOfColumns) return null;
            var result = new double[X.NumberOfRows];
            for (int i = 0; i < X.NumberOfRows; i++)
            {
                result[i] = MAP(X.Row(i));
            }
            return result;
        }

        /// <summary>
        /// Returns the class with the maximum posterior probability.
        /// </summary>
        /// <param name="x">A single vector of features for a prediction point.</param>
        private double MAP(double[] x)
        {
            int nFeatures = X.NumberOfColumns;
            int nClasses = Classes.Length;
            double max = double.MinValue;
            int maxIdx = 0;

            // P(y)  = The prior probability is the relative frequency of each class
            // P(x|y) = The conditional probability is from the PDF of Normal distribution.

            for (int i = 0; i < nClasses; i++)
            {
                // Compute the log-likelihood of each class
                double logLH = Math.Log(Priors[i]);
                for (int j = 0; j < nFeatures; j++)
                {
                    var norm = new Normal(Means[i, j], StandardDeviations[i, j]);
                    logLH += norm.LogPDF(x[j]);
                }
                // keep track of the maximum
                if (logLH > max)
                {
                    max = logLH;
                    maxIdx = i;
                }
            }
            return Classes[maxIdx];
        }

        #endregion
    }
}
 
