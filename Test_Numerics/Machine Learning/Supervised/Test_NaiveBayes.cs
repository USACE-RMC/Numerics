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

using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using System.Collections.Generic;

namespace MachineLearning
{
    /// <summary>
    /// Unit tests for Naive Bayes classification.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b> 
    /// </para>
    /// <para> <b> References: </b> </para>
    /// <para>
    /// <see href = "https://www.geeksforgeeks.org/naive-bayes-classifier-in-r-programming/" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_NaiveBayes
    {
        /// <summary>
        /// Naive Bayes is tested against the Iris dataset in R.
        /// </summary>
        [TestMethod]
        public void Test_NaiveBayes_Iris()
        {
            // Create Training data based on 70% split
            var sepalLengthTrain = new double[] { 5.1, 4.7, 5, 5.4, 5, 4.9, 5.4, 4.8, 5.8, 5.7, 5.1, 5.1, 5.4, 4.6, 4.8, 5, 5.2, 4.7, 4.8, 5.2, 4.9, 5, 4.9, 5.1, 5, 4.4, 5.1, 4.8, 4.6, 5, 7, 6.9, 6.5, 5.7, 4.9, 5.2, 5, 6, 5.6, 6.7, 5.8, 5.6, 5.9, 6.3, 6.4, 6.6, 6.7, 5.7, 5.5, 5.8, 5.4, 6, 6.3, 5.5, 5.5, 5.8, 5.6, 5.7, 6.2, 5.7, 6.3, 7.1, 6.5, 7.6, 7.3, 7.2, 6.5, 6.8, 5.8, 6.4, 7.7, 6, 6.9, 7.7, 6.7, 7.2, 6.1, 7.2, 7.4, 6.4, 6.1, 7.7, 6.4, 6.9, 6.7, 5.8, 6.7, 6.7, 6.5, 5.9 };
            var sepalWidthTrain = new double[] { 3.5, 3.2, 3.6, 3.9, 3.4, 3.1, 3.7, 3, 4, 4.4, 3.5, 3.8, 3.4, 3.6, 3.4, 3, 3.5, 3.2, 3.1, 4.1, 3.1, 3.2, 3.6, 3.4, 3.5, 3.2, 3.8, 3, 3.2, 3.3, 3.2, 3.1, 2.8, 2.8, 2.4, 2.7, 2, 2.2, 2.9, 3.1, 2.7, 2.5, 3.2, 2.5, 2.9, 3, 3, 2.6, 2.4, 2.7, 3, 3.4, 2.3, 2.5, 2.6, 2.6, 2.7, 3, 2.9, 2.8, 3.3, 3, 3, 3, 2.9, 3.6, 3.2, 3, 2.8, 3.2, 3.8, 2.2, 3.2, 2.8, 3.3, 3.2, 3, 3, 2.8, 2.8, 2.6, 3, 3.1, 3.1, 3.1, 2.7, 3.3, 3, 3, 3 };
            var petalLengthTrain = new double[] { 1.4, 1.3, 1.4, 1.7, 1.5, 1.5, 1.5, 1.4, 1.2, 1.5, 1.4, 1.5, 1.7, 1, 1.9, 1.6, 1.5, 1.6, 1.6, 1.5, 1.5, 1.2, 1.4, 1.5, 1.3, 1.3, 1.9, 1.4, 1.4, 1.4, 4.7, 4.9, 4.6, 4.5, 3.3, 3.9, 3.5, 4, 3.6, 4.4, 4.1, 3.9, 4.8, 4.9, 4.3, 4.4, 5, 3.5, 3.8, 3.9, 4.5, 4.5, 4.4, 4, 4.4, 4, 4.2, 4.2, 4.3, 4.1, 6, 5.9, 5.8, 6.6, 6.3, 6.1, 5.1, 5.5, 5.1, 5.3, 6.7, 5, 5.7, 6.7, 5.7, 6, 4.9, 5.8, 6.1, 5.6, 5.6, 6.1, 5.5, 5.4, 5.6, 5.1, 5.7, 5.2, 5.2, 5.1 };
            var petalWidthTrain = new double[] { 0.2, 0.2, 0.2, 0.4, 0.2, 0.1, 0.2, 0.1, 0.2, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.1, 0.2, 0.2, 0.1, 0.2, 0.3, 0.2, 0.4, 0.3, 0.2, 0.2, 1.4, 1.5, 1.5, 1.3, 1, 1.4, 1, 1, 1.3, 1.4, 1, 1.1, 1.8, 1.5, 1.3, 1.4, 1.7, 1, 1.1, 1.2, 1.5, 1.6, 1.3, 1.3, 1.2, 1.2, 1.3, 1.2, 1.3, 1.3, 2.5, 2.1, 2.2, 2.1, 1.8, 2.5, 2, 2.1, 2.4, 2.3, 2.2, 1.5, 2.3, 2, 2.1, 1.8, 1.8, 1.6, 1.9, 2.2, 1.4, 2.3, 1.8, 2.1, 2.4, 1.9, 2.5, 2.3, 2, 1.8 };
            var speciesTrain = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

            var trainList = new List<double[]> { sepalLengthTrain, sepalWidthTrain, petalLengthTrain, petalWidthTrain };
            var Y_training = new Vector(speciesTrain) { Header = "Species" };
            var X_training = new Matrix(trainList) { Header = new string[] { "Sepal Length", "Sepal Width", "Petal Length", "Petal Width" } };

            var naiveBayes = new NaiveBayes(X_training, Y_training);
            naiveBayes.Train();

            // Test Means
            // Setosa - class 1
            Assert.AreEqual(5.016667, naiveBayes.Means[0, 0], 1E-6);
            Assert.AreEqual(3.456667, naiveBayes.Means[0, 1], 1E-6);
            Assert.AreEqual(1.466667, naiveBayes.Means[0, 2], 1E-6);
            Assert.AreEqual(0.220000, naiveBayes.Means[0, 3], 1E-6);
            // Versicolor - class 2
            Assert.AreEqual(5.916667, naiveBayes.Means[1, 0], 1E-6);
            Assert.AreEqual(2.750000, naiveBayes.Means[1, 1], 1E-6);
            Assert.AreEqual(4.220000, naiveBayes.Means[1, 2], 1E-6);
            Assert.AreEqual(1.303333, naiveBayes.Means[1, 3], 1E-6);
            // Virginica - class 3
            Assert.AreEqual(6.740000, naiveBayes.Means[2, 0], 1E-6);
            Assert.AreEqual(3.033333, naiveBayes.Means[2, 1], 1E-6);
            Assert.AreEqual(5.680000, naiveBayes.Means[2, 2], 1E-6);
            Assert.AreEqual(2.063333, naiveBayes.Means[2, 3], 1E-6);

            // Test Standard Deviations
            // Setosa - class 1
            Assert.AreEqual(0.3097088, naiveBayes.StandardDeviations[0, 0], 1E-6);
            Assert.AreEqual(0.3490710, naiveBayes.StandardDeviations[0, 1], 1E-6);
            Assert.AreEqual(0.1881550, naiveBayes.StandardDeviations[0, 2], 1E-6);
            Assert.AreEqual(0.08051558, naiveBayes.StandardDeviations[0, 3], 1E-6);
            // Versicolor - class 2
            Assert.AreEqual(0.5414434, naiveBayes.StandardDeviations[1, 0], 1E-6);
            Assert.AreEqual(0.3202908, naiveBayes.StandardDeviations[1, 1], 1E-6);
            Assert.AreEqual(0.4373904, naiveBayes.StandardDeviations[1, 2], 1E-6);
            Assert.AreEqual(0.20924055, naiveBayes.StandardDeviations[1, 3], 1E-6);
            // Virginica - class 3
            Assert.AreEqual(0.5763141, naiveBayes.StandardDeviations[2, 0], 1E-6);
            Assert.AreEqual(0.2928261, naiveBayes.StandardDeviations[2, 1], 1E-6);
            Assert.AreEqual(0.5019960, naiveBayes.StandardDeviations[2, 2], 1E-6);
            Assert.AreEqual(0.29182344, naiveBayes.StandardDeviations[2, 3], 1E-6);

            // Test Priors
            Assert.AreEqual(30d / X_training.NumberOfRows, naiveBayes.Priors[0], 1E-6);
            Assert.AreEqual(30d / X_training.NumberOfRows, naiveBayes.Priors[1], 1E-6);
            Assert.AreEqual(30d / X_training.NumberOfRows, naiveBayes.Priors[2], 1E-6);

            // Create Test data based on 30% split
            var sepalLengthTest = new double[] { 4.9, 4.6, 4.6, 4.4, 4.8, 4.3, 5.4, 5.7, 5.1, 5.1, 5, 5.2, 5.4, 5.5, 5.5, 4.4, 4.5, 5, 5.1, 5.3, 6.4, 5.5, 6.3, 6.6, 5.9, 6.1, 5.6, 6.2, 6.1, 6.1, 6.8, 6, 5.5, 6, 6.7, 5.6, 6.1, 5, 5.7, 5.1, 5.8, 6.3, 4.9, 6.7, 6.4, 5.7, 6.5, 7.7, 5.6, 6.3, 6.2, 6.4, 7.9, 6.3, 6.3, 6, 6.9, 6.8, 6.3, 6.2 };
            var sepalWidthTest = new double[] { 3, 3.1, 3.4, 2.9, 3.4, 3, 3.9, 3.8, 3.7, 3.3, 3.4, 3.4, 3.4, 4.2, 3.5, 3, 2.3, 3.5, 3.8, 3.7, 3.2, 2.3, 3.3, 2.9, 3, 2.9, 3, 2.2, 2.8, 2.8, 2.8, 2.9, 2.4, 2.7, 3.1, 3, 3, 2.3, 2.9, 2.5, 2.7, 2.9, 2.5, 2.5, 2.7, 2.5, 3, 2.6, 2.8, 2.7, 2.8, 2.8, 3.8, 2.8, 3.4, 3, 3.1, 3.2, 2.5, 3.4 };
            var petalLengthTest = new double[] { 1.4, 1.5, 1.4, 1.4, 1.6, 1.1, 1.3, 1.7, 1.5, 1.7, 1.6, 1.4, 1.5, 1.4, 1.3, 1.3, 1.3, 1.6, 1.6, 1.5, 4.5, 4, 4.7, 4.6, 4.2, 4.7, 4.5, 4.5, 4, 4.7, 4.8, 4.5, 3.7, 5.1, 4.7, 4.1, 4.6, 3.3, 4.2, 3, 5.1, 5.6, 4.5, 5.8, 5.3, 5, 5.5, 6.9, 4.9, 4.9, 4.8, 5.6, 6.4, 5.1, 5.6, 4.8, 5.1, 5.9, 5, 5.4 };
            var petalWidthTest = new double[] { 0.2, 0.2, 0.3, 0.2, 0.2, 0.1, 0.4, 0.3, 0.4, 0.5, 0.4, 0.2, 0.4, 0.2, 0.2, 0.2, 0.3, 0.6, 0.2, 0.2, 1.5, 1.3, 1.6, 1.3, 1.5, 1.4, 1.5, 1.5, 1.3, 1.2, 1.4, 1.5, 1, 1.6, 1.5, 1.3, 1.4, 1, 1.3, 1.1, 1.9, 1.8, 1.7, 1.8, 1.9, 2, 1.8, 2.3, 2, 1.8, 1.8, 2.1, 2, 1.5, 2.4, 1.8, 2.3, 2.3, 1.9, 2.3 };
            var speciesTest = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };
            
            var testList = new List<double[]> { sepalLengthTest, sepalWidthTest, petalLengthTest, petalWidthTest };
            var Y_test = new Vector(speciesTest) { Header = "Species" };
            var X_test = new Matrix(testList) { Header = new string[] { "Sepal Length", "Sepal Width", "Petal Length", "Petal Width" } };

            // Make predictions
            var prediction = naiveBayes.Predict(X_test);
            var truePred = new double[] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3 };
            for (int i = 0; i < prediction.Length; i++)
            {
                Assert.AreEqual(truePred[i], prediction[i]);
            }

        }
    }
}
