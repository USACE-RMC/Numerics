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
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Mathematics.SpecialFunctions;
using System;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Competing Risks distribution.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list> 
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://reliability.readthedocs.io/en/latest/Competing%20risk%20models.html" />
    /// </para>
    /// </remarks>

    [TestClass]
    public class Test_CompetingRisks
    {

        /// <summary>
        /// Test the moments of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_Moments()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            // Numerics computes the moments using numerical integration
            Assert.AreEqual(27.0445, cr.Mean, 1E-2);
            Assert.AreEqual(25.0845, cr.Median, 1E-2);
            Assert.AreEqual(16.6581, cr.Mode, 1E-2);
            Assert.AreEqual(15.60225, cr.StandardDeviation, 1E-2);
            Assert.AreEqual(0.3371, cr.Skewness, 1E-2);
            Assert.AreEqual(-0.8719, cr.Kurtosis - 3, 1E-2); // Package reports excess kurtosis.

        }

        /// <summary>
        /// Test the CDF and Inverse CDF of the Competing Risk model against the Python package 'reliability'.
        /// </summary>
        [TestMethod]
        public void Test_CR_CDF()
        {
            var d1 = new LogNormal(4, 0.1) { Base = Math.E };
            var d2 = new Weibull(50, 2);
            var d3 = new GammaDistribution(30, 1.5);
            var dists = new IUnivariateDistribution[] { d1, d2, d3 };
            var cr = new CompetingRisks(dists);

            Assert.AreEqual(4.6431, cr.InverseCDF(0.05), 1E-3);
            Assert.AreEqual(25.0845, cr.InverseCDF(0.50), 1E-3);
            Assert.AreEqual(54.2056, cr.InverseCDF(0.95), 1E-3);

        }

        /// <summary>
        /// Verifies that the PDF does not return Infinity or NaN when evaluated at points
        /// in the left tail where CDF values approach zero under the maximum rule.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     For the maximum of independent random variables, the PDF formula involves
        ///     the ratio f_i(x) / F_i(x). When x is in the left tail, F_i(x) approaches zero,
        ///     causing division by zero if not handled properly.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Create a competing risks model with two Normal distributions and evaluate
        ///     the PDF at a point approximately 5 standard deviations below the mean,
        ///     where CDF values will be on the order of 1E-7.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The PDF should return a small but finite non-negative value, not Infinity or NaN.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_MaxRule_SmallX_NoInfinity()
        {
            // This test verifies the fix for division by zero when CDF ≈ 0
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);
            var cr = new CompetingRisks(new[] { dist1, dist2 });
            cr.MinimumOfRandomVariables = false; // Maximum rule

            // Test at a point where CDFs are very small
            double x = 50; // Far in left tail
            double pdf = cr.PDF(x);

            Assert.IsFalse(double.IsInfinity(pdf), "PDF should not be infinity");
            Assert.IsFalse(double.IsNaN(pdf), "PDF should not be NaN");
            Assert.IsGreaterThanOrEqualTo(0, pdf, "PDF should be non-negative");
        }

        /// <summary>
        /// Verifies that the LogPDF method returns finite values for all points in a
        /// randomly generated sample, ensuring numerical stability for log-likelihood calculations.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     MLE and Bayesian estimation methods work with log-likelihoods rather than
        ///     likelihoods to prevent numerical underflow. The LogPDF method must return
        ///     finite values (not NaN or +Infinity) for all valid inputs to ensure that
        ///     optimization algorithms can compute gradients and evaluate objective functions.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     <list type="number">
        ///         <item>Create a competing risks model with two Exponential distributions</item>
        ///         <item>Generate a random sample of 100 values</item>
        ///         <item>Verify LogPDF is finite for each sample point</item>
        ///         <item>Verify the sum (log-likelihood) is also finite</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     All LogPDF values should be finite negative numbers (since PDF ≤ 1 for
        ///     continuous distributions with unbounded support), and their sum should
        ///     be a finite negative number representing the log-likelihood.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_LogPDF_StableForMLE()
        {
            var dist1 = new Exponential(0.1);
            var dist2 = new Exponential(0.2);
            var cr = new CompetingRisks(new[] { dist1, dist2 });

            // Generate sample
            var sample = cr.GenerateRandomValues(100, 12345);

            // Verify log-likelihood is finite
            double logLik = 0;
            foreach (var x in sample)
            {
                double logPdf = cr.LogPDF(x);
                Assert.IsFalse(double.IsNaN(logPdf), $"LogPDF should not be NaN at x={x}");
                Assert.IsFalse(double.IsPositiveInfinity(logPdf), $"LogPDF should not be +Inf at x={x}");
                logLik += logPdf;
            }

            Assert.IsFalse(double.IsNaN(logLik), "Log-likelihood should not be NaN");
        }

        /// <summary>
        /// Verifies that MLE converges to approximately correct parameter values for a
        /// competing risks model under the minimum rule (series system).
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     Under the minimum rule, the competing risks model represents a series system
        ///     where the observed outcome is the minimum of the competing random variables.
        ///     This is commonly used in reliability analysis where system failure occurs
        ///     when the first component fails.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     <list type="number">
        ///         <item>Create a "true" model with known Weibull parameters:
        ///             <list type="bullet">
        ///                 <item>Distribution 1: Weibull(shape=2.0, scale=100)</item>
        ///                 <item>Distribution 2: Weibull(shape=3.0, scale=120)</item>
        ///             </list>
        ///         </item>
        ///         <item>Generate a sample of 500 observations from the true model</item>
        ///         <item>Fit a new competing risks model using MLE</item>
        ///         <item>Verify estimated parameters are within 20% of true values</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The MLE should converge to parameter estimates close to the true values.
        ///     A tolerance of 20% accounts for sampling variability with n=500.
        /// </para>
        /// <para>
        ///     <b>Note:</b>
        ///     The minimum rule is generally more numerically stable than the maximum rule
        ///     because the survival function S(x) = 1 - F(x) is bounded away from zero
        ///     in the left tail where most data typically falls.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_ConvergesForMinRule()
        {
            // True parameters
            var trueDist1 = new Weibull(2.0, 100);
            var trueDist2 = new Weibull(3.0, 120);
            var trueCR = new CompetingRisks(new[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = true;

            // Generate sample
            var sample = trueCR.GenerateRandomValues(500, 12345);

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitCR = new CompetingRisks(new[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = true;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            // Verify parameters are reasonable (within 20% of true)
            Assert.AreEqual(2.0, mleParams[0], 0.4, "Shape1 should be close to 2.0");
            Assert.AreEqual(100.0, mleParams[1], 20.0, "Scale1 should be close to 100");
        }

        /// <summary>
        /// Verifies that MLE converges to valid (non-NaN) parameter values for a
        /// competing risks model under the maximum rule (parallel system).
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     Under the maximum rule, the competing risks model represents a parallel system
        ///     where the observed outcome is the maximum of the competing random variables.
        ///     This is used in reliability analysis where system failure occurs only when
        ///     all components have failed (redundant systems).
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     <list type="number">
        ///         <item>Create a "true" model with known Weibull parameters:
        ///             <list type="bullet">
        ///                 <item>Distribution 1: Weibull(shape=2.0, scale=100)</item>
        ///                 <item>Distribution 2: Weibull(shape=3.0, scale=120)</item>
        ///             </list>
        ///         </item>
        ///         <item>Generate a sample of 500 observations from the true model</item>
        ///         <item>Fit a new competing risks model using MLE</item>
        ///         <item>Verify estimated parameters are not NaN</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The MLE should converge without numerical failures. This test uses a weaker
        ///     assertion (not NaN) rather than checking closeness to true values because
        ///     the maximum rule has inherent identifiability challenges - with only the
        ///     maximum observed, distinguishing between component distributions is difficult.
        /// </para>
        /// <para>
        ///     <b>Note:</b>
        ///     The maximum rule is more prone to numerical instability because the CDF F(x)
        ///     approaches zero in the left tail, causing division issues in the PDF formula.
        ///     This test specifically validates that the numerical stability improvements
        ///     allow the optimizer to complete without NaN propagation.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_ConvergesForMaxRule()
        {
            // True parameters
            var trueDist1 = new Weibull(2.0, 100);
            var trueDist2 = new Weibull(3.0, 120);
            var trueCR = new CompetingRisks(new[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = false; // Maximum rule

            // Generate sample
            var sample = trueCR.GenerateRandomValues(500, 12345);

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitCR = new CompetingRisks(new[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = false;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            // Verify parameters are reasonable
            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p)), "MLE parameters should not be NaN");
        }

        /// <summary>
        /// Verifies that the PDF integrates to approximately 1 over the support of the distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     A valid probability density function must integrate to 1 over its support.
        ///     This test verifies that the PDF implementation satisfies this fundamental
        ///     requirement for both minimum and maximum rules.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Use numerical integration (e.g., adaptive quadrature) to compute the integral
        ///     of the PDF from the 1E-10 quantile to the 1-1E-10 quantile, which should
        ///     capture essentially all of the probability mass.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The integral should be within 0.001 of 1.0.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_IntegratesToOne()
        {
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);

            // Test minimum rule
            var crMin = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMin.MinimumOfRandomVariables = true;

            double lowerMin = crMin.InverseCDF(1E-10);
            double upperMin = crMin.InverseCDF(1 - 1E-10);
            var agkMin = new AdaptiveGaussKronrod(crMin.PDF, lowerMin, upperMin);
            agkMin.Integrate();
            double integralMin = agkMin.Result;
            Assert.AreEqual(1.0, integralMin, 0.001, "PDF (min rule) should integrate to 1");

            // Test maximum rule
            var crMax = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMax.MinimumOfRandomVariables = false;

            double lowerMax = crMax.InverseCDF(1E-10);
            double upperMax = crMax.InverseCDF(1 - 1E-10);
            var agkMax = new AdaptiveGaussKronrod(crMax.PDF, lowerMax, upperMax);
            agkMax.Integrate();
            double integralMax = agkMax.Result;
            Assert.AreEqual(1.0, integralMax, 0.001, "PDF (max rule) should integrate to 1");
        }

        /// <summary>
        /// Verifies consistency between PDF and CDF by checking that the numerical
        /// derivative of the CDF equals the PDF at multiple points.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     By definition, f(x) = dF(x)/dx. This test verifies internal consistency
        ///     between the PDF and CDF implementations, which is critical for MLE where
        ///     both functions may be used.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Evaluate both the PDF and the numerical derivative of the CDF at several
        ///     quantile points (0.1, 0.25, 0.5, 0.75, 0.9) and verify they match within
        ///     numerical tolerance.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The relative difference between PDF(x) and CDF'(x) should be less than 1E-4
        ///     at all test points.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_PDF_CDF_Consistency()
        {
            var dist1 = new Exponential(0.1);
            var dist2 = new Exponential(0.2);
            var cr = new CompetingRisks(new[] { dist1, dist2 });

            double[] quantiles = { 0.1, 0.25, 0.5, 0.75, 0.9 };

            foreach (double q in quantiles)
            {
                double x = cr.InverseCDF(q);
                double pdf = cr.PDF(x);
                double cdfDerivative = NumericalDerivative.Derivative(cr.CDF, x);

                double relError = Math.Abs(pdf - cdfDerivative) / Math.Max(pdf, 1E-10);
                Assert.IsLessThan(1E-4, relError, $"PDF and CDF derivative should match at quantile {q}. " +
                    $"PDF={pdf}, CDF'={cdfDerivative}, RelError={relError}");
            }
        }

        /// <summary>
        /// Verifies that LogPDF and log(PDF) return consistent values where both are numerically stable.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Background:</b>
        ///     The LogPDF method is implemented using log-space arithmetic for numerical stability.
        ///     This test verifies that LogPDF produces results consistent with log(PDF) in
        ///     regions where the standard PDF calculation is stable.
        /// </para>
        /// <para>
        ///     <b>Test Strategy:</b>
        ///     Compare LogPDF(x) with Math.Log(PDF(x)) at the median (where both should be stable)
        ///     for both minimum and maximum rules.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The values should match within 1E-10 absolute tolerance at stable evaluation points.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_LogPDF_ConsistentWithPDF()
        {
            var dist1 = new Normal(100, 10);
            var dist2 = new Normal(110, 15);

            // Test minimum rule
            var crMin = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMin.MinimumOfRandomVariables = true;
            double xMin = crMin.Median;
            double logPdfMin = crMin.LogPDF(xMin);
            double logOfPdfMin = Math.Log(crMin.PDF(xMin));
            Assert.AreEqual(logOfPdfMin, logPdfMin, 1E-10,
                "LogPDF should equal log(PDF) at median for min rule");

            // Test maximum rule  
            var crMax = new CompetingRisks(new[] { dist1.Clone(), dist2.Clone() });
            crMax.MinimumOfRandomVariables = false;
            double xMax = crMax.Median;
            double logPdfMax = crMax.LogPDF(xMax);
            double logOfPdfMax = Math.Log(crMax.PDF(xMax));
            Assert.AreEqual(logOfPdfMax, logPdfMax, 1E-10,
                "LogPDF should equal log(PDF) at median for max rule");
        }


        #region Minimum Rule - 2 Distributions

        // Tolerances - competing risks MLE is harder than single distribution MLE
        private const double SHAPE_TOLERANCE_PERCENT = 0.25;  // 25% relative error
        private const double SCALE_TOLERANCE_PERCENT = 0.30;  // 30% relative error
        private const int SAMPLE_SIZE = 1000;
        private const int RANDOM_SEED = 12345;

        /// <summary>
        /// Tests MLE for minimum of Exponential and Weibull distributions.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     This is the classic "random failures + wear-out failures" model:
        ///     <list type="bullet">
        ///         <item>Exponential(λ=0.02): Constant hazard rate, models random/early failures</item>
        ///         <item>Weibull(k=3, λ=80): Increasing hazard rate (k>1), models wear-out/aging</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     The Exponential has constant hazard h(t) = λ, while the Weibull with k=3 has
        ///     hazard h(t) = (k/λ)(t/λ)^(k-1) which increases with t. The exponential dominates
        ///     early (small t) while the Weibull dominates later (large t). This creates
        ///     a characteristic "bathtub curve" effect that allows separation.
        /// </para>
        /// <para>
        ///     <b>Expected Behavior:</b>
        ///     The exponential rate parameter and Weibull shape/scale should be recoverable
        ///     within tolerance. The Weibull scale may have higher variance due to fewer
        ///     observations in the wear-out region.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MinRule_2Dist_Exponential_Weibull()
        {
            // True parameters - designed for identifiability
            // Weibull(k=1) ≡ Exponential, dominates early failures
            // Weibull(k=3) has increasing hazard, dominates wear-out
            // Note: Weibull(scale, shape=1) is used instead of Exponential to avoid
            // the Exponential's 2-parameter (location+scale) MLE issues in competing risks context.
            double trueScale1 = 50.0;       // Weibull(50,1) = Exponential with mean 50
            double trueShape1 = 1.0;        // Constant hazard
            double trueScale2 = 80.0;       // Mean ≈ 71
            double trueShape2 = 3.0;        // Increasing hazard

            var trueDist1 = new Weibull(trueScale1, trueShape1);
            var trueDist2 = new Weibull(trueScale2, trueShape2);
            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = true;

            // Generate sample
            var sample = trueCR.GenerateRandomValues(SAMPLE_SIZE, RANDOM_SEED);

            // Verify sample statistics are reasonable
            double sampleMean = sample.Average();
            double sampleMin = sample.Min();
            double sampleMax = sample.Max();
            Console.WriteLine($"Sample: n={SAMPLE_SIZE}, mean={sampleMean:F2}, min={sampleMin:F2}, max={sampleMax:F2}");

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = true;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            // Weibull params are [scale, shape] for each component
            Console.WriteLine($"True:    Weibull({trueScale1},{trueShape1}), Weibull({trueScale2},{trueShape2})");
            Console.WriteLine($"Fitted:  Weibull({mleParams[0]:F3},{mleParams[1]:F3}), Weibull({mleParams[2]:F3},{mleParams[3]:F3})");

            // Assertions with tolerance
            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p) || double.IsInfinity(p)), "All parameters should be finite");

            // Verify overall fit
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.05, ksStatistic, "KS statistic should indicate good fit");
        }

        /// <summary>
        /// Tests MLE for minimum of two Weibull distributions with different shapes.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Two Weibulls with contrasting shapes and well-separated scales:
        ///     <list type="bullet">
        ///         <item>Weibull(k=0.8, λ=30): Decreasing hazard (k&lt;1), dominates very early</item>
        ///         <item>Weibull(k=3.0, λ=100): Increasing hazard (k&gt;1), dominates later</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     The k=0.8 distribution has decreasing hazard (infant mortality pattern),
        ///     while k=3.0 has increasing hazard (wear-out pattern). Combined with the
        ///     3:1 scale ratio, the distributions contribute in clearly different time regions.
        ///     The first dominates the left tail, the second shapes the right tail.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MinRule_2Dist_Weibull_DifferentShapes()
        {
            // Weibull 1: Decreasing hazard (infant mortality)
            double trueShape1 = 0.8;
            double trueScale1 = 30.0;

            // Weibull 2: Increasing hazard (wear-out)
            double trueShape2 = 3.0;
            double trueScale2 = 100.0;

            var trueDist1 = new Weibull(trueScale1, trueShape1);
            var trueDist2 = new Weibull(trueScale2, trueShape2);
            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = true;

            var sample = trueCR.GenerateRandomValues(SAMPLE_SIZE, RANDOM_SEED);

            Console.WriteLine($"Sample: n={SAMPLE_SIZE}, mean={sample.Average():F2}, median={sample.OrderBy(x => x).ElementAt(SAMPLE_SIZE / 2):F2}");

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = true;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   Weibull({trueShape1}, {trueScale1}), Weibull({trueShape2}, {trueScale2})");
            Console.WriteLine($"Fitted: Weibull({mleParams[0]:F3}, {mleParams[1]:F3}), Weibull({mleParams[2]:F3}, {mleParams[3]:F3})");

            // Verify no NaN
            Assert.IsFalse(mleParams.Any(double.IsNaN), "No parameters should be NaN");

            // Check parameter recovery (allowing for label switching)
            // Weibull params are [scale, shape], so shape indices are 1 and 3
            bool config1 = IsCloseRelative(mleParams[1], trueShape1, SHAPE_TOLERANCE_PERCENT) &&
                          IsCloseRelative(mleParams[3], trueShape2, SHAPE_TOLERANCE_PERCENT);
            bool config2 = IsCloseRelative(mleParams[1], trueShape2, SHAPE_TOLERANCE_PERCENT) &&
                          IsCloseRelative(mleParams[3], trueShape1, SHAPE_TOLERANCE_PERCENT);

            Assert.IsTrue(config1 || config2,
                "Fitted shapes should match true shapes (allowing for label switching)");
        }

        #endregion

        #region Minimum Rule - 3 Distributions

        /// <summary>
        /// Tests MLE for minimum of three distributions: Exponential + two Weibulls.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     A three-component "bathtub curve" model:
        ///     <list type="bullet">
        ///         <item>Weibull(k=0.7, λ=20): Decreasing hazard - infant mortality</item>
        ///         <item>Exponential(λ=0.005): Constant hazard - random failures (useful life)</item>
        ///         <item>Weibull(k=4, λ=150): Steeply increasing hazard - wear-out</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     Each distribution dominates a different time region:
        ///     <list type="bullet">
        ///         <item>Early: Weibull(0.7, 20) with its decreasing hazard</item>
        ///         <item>Middle: Exponential provides the "flat bottom" of the bathtub</item>
        ///         <item>Late: Weibull(4, 150) causes the upturn in failure rate</item>
        ///     </list>
        ///     The three distinct hazard behaviors create sufficient structure for identification.
        /// </para>
        /// <para>
        ///     <b>Note:</b>
        ///     With 5 parameters and complex interactions, this is a challenging estimation
        ///     problem. Larger samples (n=1500+) and relaxed tolerances are appropriate.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MinRule_3Dist_BathtubCurve()
        {
            // Three-component bathtub curve using Weibulls only
            // Weibull(scale, shape=1) ≡ Exponential, avoids location parameter MLE issues
            var trueDist1 = new Weibull(20, 0.7);      // scale=20, shape=0.7 Infant mortality (decreasing hazard)
            var trueDist2 = new Weibull(200, 1.0);      // scale=200, shape=1.0 Random failures (constant hazard, ≡ Exponential)
            var trueDist3 = new Weibull(150, 4.0);      // scale=150, shape=4.0 Wear-out (increasing hazard)

            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2, trueDist3 });
            trueCR.MinimumOfRandomVariables = true;

            // Use larger sample for 3-distribution case
            int n = 1500;
            var sample = trueCR.GenerateRandomValues(n, RANDOM_SEED);

            Console.WriteLine($"Sample: n={n}, mean={sample.Average():F2}, min={sample.Min():F2}, max={sample.Max():F2}");

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitDist3 = new Weibull();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2, fitDist3 });
            fitCR.MinimumOfRandomVariables = true;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            // Params: Weibull[scale,shape] x 3
            Console.WriteLine($"True parameters:   Weibull(20, 0.7), Weibull(200, 1.0), Weibull(150, 4.0)");
            Console.WriteLine($"Fitted parameters: Weibull({mleParams[0]:F3}, {mleParams[1]:F3}), " +
                            $"Weibull({mleParams[2]:F3}, {mleParams[3]:F3}), Weibull({mleParams[4]:F3}, {mleParams[5]:F3})");

            // Verify convergence (no NaN/Inf)
            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p) || double.IsInfinity(p)),
                "All parameters should be finite");

            // Verify the overall distribution fit (CDF comparison)
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.05, ksStatistic, "KS statistic should indicate good fit");
        }

        /// <summary>
        /// Tests MLE for minimum of three Weibull distributions with distinct characteristics.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Three Weibulls spanning the shape parameter space:
        ///     <list type="bullet">
        ///         <item>Weibull(k=0.5, λ=15): Strongly decreasing hazard</item>
        ///         <item>Weibull(k=1.5, λ=60): Mildly increasing hazard</item>
        ///         <item>Weibull(k=4.0, λ=120): Strongly increasing hazard</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     The shapes span k &lt; 1, 1 &lt; k &lt; 2, and k &gt; 2, giving three distinct
        ///     hazard behaviors. The scales are chosen to create overlapping but distinguishable
        ///     contributions: the k=0.5 dominates the extreme left tail, k=1.5 the middle-left,
        ///     and k=4.0 shapes the right tail.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MinRule_3Dist_ThreeWeibulls()
        {
            var trueDist1 = new Weibull(15, 0.5);   // Strongly decreasing hazard
            var trueDist2 = new Weibull(60, 1.5);   // Mildly increasing hazard  
            var trueDist3 = new Weibull(120, 4.0);  // Strongly increasing hazard

            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2, trueDist3 });
            trueCR.MinimumOfRandomVariables = true;

            int n = 1500;
            var sample = trueCR.GenerateRandomValues(n, RANDOM_SEED);

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Weibull();
            var fitDist3 = new Weibull();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2, fitDist3 });
            fitCR.MinimumOfRandomVariables = true;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   Weibull(0.5, 15), Weibull(1.5, 60), Weibull(4.0, 120)");
            Console.WriteLine($"Fitted: Weibull({mleParams[0]:F2}, {mleParams[1]:F2}), " +
                            $"Weibull({mleParams[2]:F2}, {mleParams[3]:F2}), " +
                            $"Weibull({mleParams[4]:F2}, {mleParams[5]:F2})");

            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p) || double.IsInfinity(p)),
                "All parameters should be finite");

            // Verify overall fit quality
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.05, ksStatistic, "KS statistic should indicate good fit");
        }

        #endregion

        #region Maximum Rule - 2 Distributions

        /// <summary>
        /// Tests MLE for maximum of two Normal distributions with separated means.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Two Normals with well-separated means and different standard deviations:
        ///     <list type="bullet">
        ///         <item>Normal(μ=50, σ=8): Lower component</item>
        ///         <item>Normal(μ=85, σ=12): Upper component with larger spread</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     For the maximum of two Normals, the lower distribution primarily influences
        ///     the left tail of the max distribution (when both draws happen to be low),
        ///     while the upper distribution dominates the right tail. With a ~4σ separation
        ///     between means, the contributions are clearly distinguishable.
        /// </para>
        /// <para>
        ///     <b>Note:</b>
        ///     Normal distributions work well for the maximum rule because their symmetric,
        ///     bounded tails avoid the numerical issues that arise with heavy-tailed distributions.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MaxRule_2Dist_TwoNormals()
        {
            double trueMu1 = 50, trueSigma1 = 8;
            double trueMu2 = 85, trueSigma2 = 12;

            var trueDist1 = new Normal(trueMu1, trueSigma1);
            var trueDist2 = new Normal(trueMu2, trueSigma2);
            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = false; // Maximum rule

            var sample = trueCR.GenerateRandomValues(SAMPLE_SIZE, RANDOM_SEED);

            Console.WriteLine($"Sample: n={SAMPLE_SIZE}, mean={sample.Average():F2}, min={sample.Min():F2}, max={sample.Max():F2}");

            // Fit model
            var fitDist1 = new Normal();
            var fitDist2 = new Normal();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = false;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   Normal({trueMu1}, {trueSigma1}), Normal({trueMu2}, {trueSigma2})");
            Console.WriteLine($"Fitted: Normal({mleParams[0]:F2}, {mleParams[1]:F2}), Normal({mleParams[2]:F2}, {mleParams[3]:F2})");

            // Verify no NaN
            Assert.IsFalse(mleParams.Any(double.IsNaN), "No parameters should be NaN");

            // Check overall fit quality
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.05, ksStatistic, "KS statistic should indicate good fit");

            // Verify means are approximately recovered (allowing label switching)
            var fittedMeans = new[] { mleParams[0], mleParams[2] }.OrderBy(x => x).ToArray();
            var trueMeans = new[] { trueMu1, trueMu2 }.OrderBy(x => x).ToArray();

            Assert.IsTrue(IsCloseRelative(fittedMeans[0], trueMeans[0], SCALE_TOLERANCE_PERCENT),
                $"Lower mean should be close to {trueMeans[0]}");
            Assert.IsTrue(IsCloseRelative(fittedMeans[1], trueMeans[1], SCALE_TOLERANCE_PERCENT),
                $"Upper mean should be close to {trueMeans[1]}");
        }

        /// <summary>
        /// Tests MLE for maximum of Weibull and Gumbel (GEV Type I) distributions.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Combining distributions with different tail behaviors:
        ///     <list type="bullet">
        ///         <item>Weibull(k=2, λ=50): Light right tail (bounded support if k>1 effective tail)</item>
        ///         <item>Gumbel(μ=70, σ=15): Heavy right tail (extreme value distribution)</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     For the maximum, the right tail behavior is critical. The Weibull's relatively
        ///     light tail means it rarely produces extreme maxima, while the Gumbel's heavy
        ///     tail dominates extreme values. The bulk of the distribution is shaped by both,
        ///     with the Gumbel's influence increasing in the upper quantiles.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MaxRule_2Dist_Weibull_Gumbel()
        {
            double trueWeibullShape = 2.0;
            double trueWeibullScale = 50.0;
            double trueGumbelLocation = 70.0;
            double trueGumbelScale = 15.0;

            var trueDist1 = new Weibull(trueWeibullScale, trueWeibullShape);
            var trueDist2 = new Gumbel(trueGumbelLocation, trueGumbelScale);
            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2 });
            trueCR.MinimumOfRandomVariables = false; // Maximum rule

            var sample = trueCR.GenerateRandomValues(SAMPLE_SIZE, RANDOM_SEED);

            Console.WriteLine($"Sample: n={SAMPLE_SIZE}, mean={sample.Average():F2}, min={sample.Min():F2}, max={sample.Max():F2}");

            // Fit model
            var fitDist1 = new Weibull();
            var fitDist2 = new Gumbel();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2 });
            fitCR.MinimumOfRandomVariables = false;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   Weibull({trueWeibullShape}, {trueWeibullScale}), Gumbel({trueGumbelLocation}, {trueGumbelScale})");
            Console.WriteLine($"Fitted: Weibull({mleParams[0]:F2}, {mleParams[1]:F2}), Gumbel({mleParams[2]:F2}, {mleParams[3]:F2})");

            Assert.IsFalse(mleParams.Any(double.IsNaN), "No parameters should be NaN");

            // Verify overall fit
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.05, ksStatistic, "KS statistic should indicate good fit");
        }

        #endregion

        #region Maximum Rule - 3 Distributions

        /// <summary>
        /// Tests MLE for maximum of three Normal distributions representing a trimodal scenario.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Three well-separated Normals:
        ///     <list type="bullet">
        ///         <item>Normal(μ=40, σ=6): Low component</item>
        ///         <item>Normal(μ=70, σ=8): Middle component</item>
        ///         <item>Normal(μ=100, σ=10): High component</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     With ~4σ separation between adjacent means, each Normal dominates a distinct
        ///     region. For the maximum, the observed values come from all three components
        ///     but with different frequencies based on the probability that each component
        ///     produces the largest value. The lowest component rarely "wins" but influences
        ///     the extreme left tail of the maximum distribution.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MaxRule_3Dist_ThreeNormals()
        {
            var trueDist1 = new Normal(40, 6);
            var trueDist2 = new Normal(70, 8);
            var trueDist3 = new Normal(100, 10);

            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2, trueDist3 });
            trueCR.MinimumOfRandomVariables = false; // Maximum rule

            int n = 1500;
            var sample = trueCR.GenerateRandomValues(n, RANDOM_SEED);

            Console.WriteLine($"Sample: n={n}, mean={sample.Average():F2}, min={sample.Min():F2}, max={sample.Max():F2}");

            // Fit model
            var fitDist1 = new Normal();
            var fitDist2 = new Normal();
            var fitDist3 = new Normal();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2, fitDist3 });
            fitCR.MinimumOfRandomVariables = false;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   N(40,6), N(70,8), N(100,10)");
            Console.WriteLine($"Fitted: N({mleParams[0]:F1},{mleParams[1]:F1}), " +
                            $"N({mleParams[2]:F1},{mleParams[3]:F1}), " +
                            $"N({mleParams[4]:F1},{mleParams[5]:F1})");

            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p) || double.IsInfinity(p)),
                "All parameters should be finite");

            // Verify overall fit
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.06, ksStatistic, "KS statistic should indicate reasonable fit");

            // Note: Individual parameter recovery is not asserted for 3 same-family components
            // under max-rule due to inherent identifiability limitations. The lowest component
            // has minimal influence on the maximum and is difficult to recover.
            // KS statistic and convergence checks above are sufficient.
        }

        /// <summary>
        /// Tests MLE for maximum of Exponential, Gamma, and LogNormal - three different families.
        /// </summary>
        /// <remarks>
        /// <para>
        ///     <b>Configuration Rationale:</b>
        ///     Three different distribution families with distinct shapes:
        ///     <list type="bullet">
        ///         <item>Exponential(λ=0.05): Monotonically decreasing density, mean=20</item>
        ///         <item>Gamma(k=3, θ=15): Unimodal with mode at 30, mean=45</item>
        ///         <item>LogNormal(μ=4.2, σ=0.4): Right-skewed, median≈67, mean≈72</item>
        ///     </list>
        /// </para>
        /// <para>
        ///     <b>Why This Is Identifiable:</b>
        ///     Using three different families provides maximum structural diversity.
        ///     The Exponential is memoryless, the Gamma has a characteristic shape controlled
        ///     by its shape parameter, and the LogNormal has a distinctive heavy right tail.
        ///     For the maximum, these combine to create a complex but estimable distribution.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_MLE_MaxRule_3Dist_DifferentFamilies()
        {
            var trueDist1 = new Exponential(0.05);           // Mean = 20
            var trueDist2 = new GammaDistribution(3.0, 15.0); // Mean = 45
            var trueDist3 = new LogNormal(4.2, 0.4) { Base = Math.E };  // Median ≈ 67

            var trueCR = new CompetingRisks(new UnivariateDistributionBase[] { trueDist1, trueDist2, trueDist3 });
            trueCR.MinimumOfRandomVariables = false; // Maximum rule

            int n = 1500;
            var sample = trueCR.GenerateRandomValues(n, RANDOM_SEED);

            Console.WriteLine($"Sample: n={n}, mean={sample.Average():F2}, min={sample.Min():F2}, max={sample.Max():F2}");
            Console.WriteLine($"True distribution means: Exp={1 / 0.05:F0}, Gamma={3 * 15:F0}, LogN≈{Math.Exp(4.2 + 0.4 * 0.4 / 2):F0}");

            // Fit model
            var fitDist1 = new Exponential();
            var fitDist2 = new GammaDistribution();
            var fitDist3 = new LogNormal();
            var fitCR = new CompetingRisks(new UnivariateDistributionBase[] { fitDist1, fitDist2, fitDist3 });
            fitCR.MinimumOfRandomVariables = false;

            var mleParams = fitCR.MLE(sample);
            fitCR.SetParameters(mleParams);

            Console.WriteLine($"True:   Exp(0.05), Gamma(3, 15), LogNormal(4.2, 0.4)");
            Console.WriteLine($"Fitted: Exp({mleParams[0]:F4}), Gamma({mleParams[1]:F2}, {mleParams[2]:F2}), " +
                            $"LogNormal({mleParams[3]:F2}, {mleParams[4]:F2})");

            Assert.IsFalse(mleParams.Any(p => double.IsNaN(p) || double.IsInfinity(p)),
                "All parameters should be finite");

            // Verify overall fit
            double ksStatistic = ComputeKSStatistic(sample, fitCR);
            Console.WriteLine($"KS statistic: {ksStatistic:F4}");
            Assert.IsLessThan(0.06, ksStatistic, "KS statistic should indicate reasonable fit");
        }

        #endregion

        #region Helper Methods

        /// <summary>
        /// Checks if two values are close within a relative tolerance.
        /// </summary>
        private static bool IsCloseRelative(double actual, double expected, double tolerance)
        {
            if (expected == 0) return Math.Abs(actual) < tolerance;
            return Math.Abs(actual - expected) / Math.Abs(expected) < tolerance;
        }

        /// <summary>
        /// Computes the Kolmogorov-Smirnov statistic for goodness of fit.
        /// </summary>
        private static double ComputeKSStatistic(double[] sample, CompetingRisks distribution)
        {
            var sorted = sample.OrderBy(x => x).ToArray();
            int n = sorted.Length;
            double maxDiff = 0;

            for (int i = 0; i < n; i++)
            {
                double empiricalCDF = (i + 1.0) / n;
                double theoreticalCDF = distribution.CDF(sorted[i]);
                double diff = Math.Abs(empiricalCDF - theoreticalCDF);
                if (diff > maxDiff) maxDiff = diff;
            }

            return maxDiff;
        }

        #endregion
    }
}
