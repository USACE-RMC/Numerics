using System;
using System.Diagnostics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Generalized Pareto distribution algorithm.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     <item> Tiki Gonzalez, USACE Risk Management Center, julian.t.gonzalez@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// <see href = "https://github.com/mathnet/mathnet-numerics/blob/master/src/Numerics.Tests/DistributionTests" />
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_GeneralizedPareto
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 8.3.1 White River at Mt. Carmel (Threshold = 50,000 cfs)
        private double[] sample = new double[] { 126000d, 148000d, 66000d, 156000d, 136000d, 122000d, 183000d, 162000d, 85200d, 56800d, 56600d, 138000d, 81000d, 51800d, 90700d, 139000d, 160000d, 118000d, 50600d, 137000d, 151000d, 172000d, 52600d, 248000d, 152000d, 64500d, 61500d, 143000d, 108000d, 53000d, 134000d, 115000d, 84100d, 105000d, 85400d, 76900d, 99100d, 73700d, 122000d, 62500d, 54300d, 58000d, 144000d, 55800d, 127000d, 55800d, 107000d, 56400d, 128000d, 106000d, 110000d, 232000d, 60400d, 60400d, 50800d, 100000d, 55900d, 167000d, 53700d, 56700d, 126000d, 100000d, 59500d, 164000d, 81800d, 56400d, 124000d, 64600d, 77300d, 65900d, 72500d, 65100d, 80800d, 69800d, 53000d, 195000d, 128000d, 114000d, 110000d, 149000d, 74100d, 75900d, 99300d, 168000d, 70800d, 104000d, 125000d, 77300d, 97300d, 140000d, 54900d, 66000d, 199000d, 99800d, 105000d, 93700d, 277000d, 85700d, 77300d, 122000d, 106000d, 93300d, 79000d, 130000d, 126000d, 57800d, 64700d, 162000d, 71900d, 63500d, 81500d, 51000d, 84400d, 108000d, 185000d, 55800d, 94600d, 82800d, 146000d, 66500d, 57700d, 78700d, 85100d, 129000d, 75700d, 104000d, 139000d, 50600d, 53500d, 178000d, 110000d, 50800d, 76000d, 130000d, 67300d, 149000d, 78400d, 96600d, 83300d, 68400d, 84300d, 56400d, 112000d, 76400d, 116000d, 51400d, 59800d, 63900d, 81900d, 88200d, 62300d, 162000d, 67200d, 85500d, 51000d, 286000d, 73800d, 61300d, 60800d, 91300d, 134000d, 106000d, 70800d, 106000d, 122000d, 149000d, 53700d, 85300d, 144000d, 54800d, 116000d, 67500d, 56500d, 86700d, 91500d, 105000d, 134000d, 97300d, 84000d, 141000d, 52600d, 124000d, 196000d, 84200d, 54500d, 74500d, 104000d, 57200d, 61000d, 155000d, 96500d, 89100d, 77900d, 70500d, 73400d, 180000d, 83700d, 302000d, 133000d, 92100d, 105000d, 235000d, 213000d, 96100d, 77100d, 73900d, 55400d, 55200d, 87800d, 52600d, 106000d, 93000d, 147000d, 61800d, 101000d, 154000d, 52000d, 121000d, 86700d, 57300d, 97500d, 112000d, 88500d, 76200d, 140000d, 87400d, 154000d, 95100d, 131000d, 131000d, 54900d, 78800d, 101000d, 224000d, 54800d, 50900d, 63500d, 63500d, 152000d, 51000d, 285000d, 114000d, 197000d, 106000d, 132000d, 83700d, 67200d, 110000d, 202000d, 127000d, 90600d, 126000d, 73900d, 86500d, 181000d, 141000d, 79700d, 97800d, 57300d, 77200d, 133000d, 82900d, 55000d, 62000d, 51700d, 54500d, 51600d, 103000d, 134000d, 71700d, 57000d, 63900d, 60700d, 81900d, 171000d, 111000d, 50400d, 50500d, 69700d, 88900d, 76600d };

        /// <summary>
        /// Verification of GPA Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 8.3.1 page 279.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_MOM_Fit()
        {
            var GPA = new GeneralizedPareto();
            GPA.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = GPA.Xi;
            double a = GPA.Alpha;
            double k = GPA.Kappa;
            double true_x = 50169.23d;
            double true_a = 55443d;
            double true_k = 0.0956d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x );
            Assert.IsLessThan(0.01d, (a - true_a) / true_a );
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Verification of GPA Distribution fit with method of linear moments.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var GPA = new GeneralizedPareto();
            GPA.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = GPA.Xi;
            double a = GPA.Alpha;
            double k = GPA.Kappa;
            double true_x = 1285.909d;
            double true_a = 589.7772d;
            double true_k = 0.6251903d;
            Assert.AreEqual(true_x, x, 0.001d);
            Assert.AreEqual(true_a, a, 0.001d);
            Assert.AreEqual(true_k, k, 0.001d);
            var lmom = GPA.LinearMomentsFromParameters(GPA.GetParameters);
            Assert.AreEqual(1648.806d, lmom[0],  0.001d);
            Assert.AreEqual(138.2366d, lmom[1],  0.001d);
            Assert.AreEqual(0.1033903d, lmom[2], 0.001d);
            Assert.AreEqual(0.03073215d, lmom[3], 0.001d);
        }

        /// <summary>
        /// Verification of GPA Distribution fit with modified method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 8.3.1 page 279.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_ModMOM_Fit()
        {
            var GPA = new GeneralizedPareto();
            GPA.SetParameters(GPA.ModifiedMethodOfMoments(sample));
            double x = GPA.Xi;
            double a = GPA.Alpha;
            double k = GPA.Kappa;
            double true_x = 50203.04d;
            double true_a = 55365.72d;
            double true_k = 0.0948d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x);
            Assert.IsLessThan(0.01d, (a - true_a) / true_a);
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Verification of Generalized Pareto Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 8.3.1 page 279.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_MLE_Fit()
        {
            var GPA = new GeneralizedPareto();
            GPA.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = GPA.Xi;
            double a = GPA.Alpha;
            double k = GPA.Kappa;
            double true_x = 50400d;
            double true_a = 55142.29d;
            double true_k = 0.0945d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x);
            Assert.IsLessThan(0.01d, (a - true_a) / true_a);
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Test the maximum likelihood fit with a known location against an independent reference.
        /// </summary>
        /// <remarks>
        /// The reference maximizes the same log-likelihood with SciPy (scipy.stats.genpareto.logpdf, shape
        /// c = −κ, location fixed at 0) by Nelder-Mead at xatol = 1E-12: α = 2.388710825197821,
        /// κ = −0.07721168825622311, log-likelihood −38.959310020261206. SciPy's own
        /// genpareto.fit(x, floc=0) gives α = 2.38873, κ = −0.07721. The log-likelihood is flat near its
        /// maximum, so the Nelder-Mead default convergence (relative change 1E-8 in the log-likelihood, as for
        /// every Numerics MLE) places the parameters within about 2E-4 relative of the maximum.
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_MLE_FixedLocation_Reference()
        {
            var excesses = new double[] { 1.2, 0.4, 3.9, 2.2, 0.15, 5.6, 1.05, 0.8, 7.3, 2.9, 0.33, 1.7, 4.4, 0.62, 2.5, 11.2, 0.09, 1.33, 3.1, 0.95 };
            var GPA = new GeneralizedPareto();
            var parameters = GPA.MLE(excesses, 0d);
            Assert.AreEqual(0d, parameters[0]);
            Assert.AreEqual(2.388710825197821, parameters[1], 1E-3 * 2.388710825197821);
            Assert.AreEqual(-0.07721168825622311, parameters[2], 1E-3);
            var fitted = new GeneralizedPareto(parameters[0], parameters[1], parameters[2]);
            Assert.AreEqual(-38.959310020261206, fitted.LogLikelihood(excesses), 1E-6);

            // Shifting the sample and the location together leaves the fit unchanged.
            var shifted = new double[excesses.Length];
            for (int i = 0; i < excesses.Length; i++) shifted[i] = excesses[i] + 10d;
            var moved = GPA.MLE(shifted, 10d);
            Assert.AreEqual(10d, moved[0]);
            Assert.AreEqual(parameters[1], moved[1], 1E-3 * parameters[1]);
            Assert.AreEqual(parameters[2], moved[2], 1E-3);

            // The fixed-location overload leaves the existing fit untouched.
            var minimumFit = new GeneralizedPareto().MLE(sample);
            Assert.AreEqual(Statistics.Minimum(sample), minimumFit[0]);
        }

        /// <summary>
        /// Test the maximum likelihood fit with a known location on the Rao and Hamed sample at its threshold.
        /// </summary>
        /// <remarks>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000, Example 8.3.1
        /// page 279, whose MLE has ξ = 50400 (the sample minimum), α = 55142.29 and κ = 0.0945.
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_MLE_FixedLocation_Threshold()
        {
            var parameters = new GeneralizedPareto().MLE(sample, 50400d);
            Assert.AreEqual(50400d, parameters[0]);
            Assert.AreEqual(55142.29d, parameters[1], 0.01 * 55142.29d);
            Assert.AreEqual(0.0945d, parameters[2], 0.01 * 0.0945d);
        }

        /// <summary>
        /// Test that the known-location fit recovers the parameters of a large simulated sample, including a
        /// bounded tail whose linear-moment start can fall outside the support.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_MLE_FixedLocation_Recovery()
        {
            foreach (var (alpha, kappa) in new[] { (10d, -0.2), (10d, 0d), (10d, 0.3) })
            {
                var truth = new GeneralizedPareto(5d, alpha, kappa);
                var draws = truth.GenerateRandomValues(20000, 12345);
                var parameters = new GeneralizedPareto().MLE(draws, 5d);
                Assert.AreEqual(5d, parameters[0]);
                Assert.AreEqual(alpha, parameters[1], 0.05 * alpha, $"kappa = {kappa}");
                Assert.AreEqual(kappa, parameters[2], 0.03, $"kappa = {kappa}");
                var fitted = new GeneralizedPareto(parameters[0], parameters[1], parameters[2]);
                double logLikelihood = fitted.LogLikelihood(draws);
                Assert.IsFalse(double.IsNaN(logLikelihood) || double.IsInfinity(logLikelihood));
            }
        }

        /// <summary>
        /// Test the argument checks of the known-location fit.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_MLE_FixedLocation_InvalidArguments()
        {
            var GPA = new GeneralizedPareto();
            Assert.ThrowsExactly<ArgumentNullException>(() => GPA.MLE(null!, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d, 2d }, double.NaN));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d }, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d, -0.5 }, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d, double.PositiveInfinity }, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 3d, 3d }, 3d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 5d, 5d, 5d }, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d, 2d }, double.NegativeInfinity));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 1d, double.NaN }, 0d));
        }

        /// <summary>
        /// Test the known-location fit at its boundaries: two observations, observations exactly at the location,
        /// and a sample whose linear-moment start leaves an observation outside the support (the exponential
        /// start is then used). Each fit is valid, keeps the location, and is not worse than its start.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_MLE_FixedLocation_Boundaries()
        {
            var cases = new (double[] Sample, double Location)[]
            {
                (new[] { 1d, 3d }, 0d),
                (new[] { 0d, 0d, 0.5, 2d, 4.5, 0d, 1.25 }, 0d),
                // Two tight clusters: the linear-moment start is κ = 4.07, α = 36.2, whose upper endpoint
                // α/κ = 8.91 lies below the largest observations, so the exponential start is used.
                (new[] { 5d, 5.1, 5.2, 5.3, 9.0, 9.1, 9.2, 9.3 }, 0d),
            };
            foreach (var (sample, location) in cases)
            {
                var parameters = new GeneralizedPareto().MLE(sample, location);
                Assert.AreEqual(location, parameters[0]);
                var fitted = new GeneralizedPareto(parameters[0], parameters[1], parameters[2]);
                Assert.IsTrue(fitted.ParametersValid);
                double logLikelihood = fitted.LogLikelihood(sample);
                Assert.IsFalse(double.IsNaN(logLikelihood) || double.IsInfinity(logLikelihood));
                double mean = 0;
                foreach (double x in sample) mean += x - location;
                mean /= sample.Length;
                var exponential = new GeneralizedPareto(location, mean, 0d);
                Assert.IsGreaterThanOrEqualTo(exponential.LogLikelihood(sample) - 1E-9, logLikelihood);
            }
        }

        /// <summary>
        /// Test the quantile function for the GPA Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 8.3.2 page 283.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_Quantile()
        {
            var GPA = new GeneralizedPareto(50203.04d, 55365.72d, 0.0948d);
            double q100 = GPA.InverseCDF(0.99d);
            double true_q100 = 256803d;
            Assert.IsLessThan(0.01d, (q100 - true_q100) / true_q100);
            double p = GPA.CDF(q100);
            double true_p = 0.99d;
            Assert.IsLessThan(0.01d, (p - true_p) / true_p);
        }

        /// <summary>
        /// The list quantile overload gives a growth curve bitwise equal to the scalar quantile.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_GrowthCurveMatchesScalarInverseCDF()
        {
            var probabilities = new double[] { 0d, 1E-12, 0.01, 0.1, 0.5, 0.9, 0.98, 0.99, 0.999, 1 - 1E-12, 1d };
            foreach (double k in new[] { -0.3, -1E-9, 0d, 0.2, 1.5 })
            {
                var gpa = new GeneralizedPareto(0.2d, 0.8d, k);
                double[] curve = gpa.InverseCDF(probabilities);
                for (int i = 0; i < probabilities.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(gpa.InverseCDF(probabilities[i])), BitConverter.DoubleToInt64Bits(curve[i]));
            }
        }

        /// <summary>
        /// Test the partial derivatives for the Generalized Pareto Distribution.
        /// </summary>
        [TestMethod()]
        public void Test_GPA_Partials()
        {
            var GPA = new GeneralizedPareto(50203.04d, 55365.72d, 0.0948d);
            double dQdLocation = GPA.QuantileGradient(0.99d)[0];
            double dQdScale = GPA.QuantileGradient(0.99d)[1];
            double dQdShape = GPA.QuantileGradient(0.99d)[2];
            double true_dLocation = 1.0d;
            double true_dScale = 3.7315488d;
            double true_dShape = -441209.53d;
            Assert.IsLessThan(0.01d, (dQdLocation - true_dLocation) / true_dLocation);
            Assert.IsLessThan(0.01d, (dQdScale - true_dScale) / true_dScale);
            Assert.IsLessThan(0.01d, (dQdShape - true_dShape) / true_dShape);
        }

        /// <summary>
        /// Test the standard error for the Generalized Pareto Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 8.3.3 page 286.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GPA_StandardError()
        {

            // Method of Moments
            var GPA = new GeneralizedPareto(50203.04d, 55365.72d, 0.0948d);
            double qVar99 = Math.Sqrt(GPA.QuantileVariance(0.99d, sample.Length, ParameterEstimationMethod.MethodOfMoments));
            double true_qVar99 = 16657d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99);

            // Maximum Likelihood
            GPA = new GeneralizedPareto(50400d, 55142.29d, 0.0945d);
            qVar99 = Math.Sqrt(GPA.QuantileVariance(0.99d, sample.Length, ParameterEstimationMethod.MaximumLikelihood));
            true_qVar99 = 15938d;
            Assert.IsLessThan(0.01d, (qVar99 - true_qVar99) / true_qVar99);
        }

        /// <summary>
        /// Checking if parameters can construct Generalized Pareto.
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(100,GPA.Xi);
            Assert.AreEqual(10,GPA.Alpha);
            Assert.AreEqual(0, GPA.Kappa);

            var GPA2 = new GeneralizedExtremeValue(-100, 1, 1);
            Assert.AreEqual(-100,GPA2.Xi);
            Assert.AreEqual(1, GPA2.Alpha);
            Assert.AreEqual(1, GPA2.Kappa);
        }

        /// <summary>
        /// Testing distribution with bad parameters.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var GPA = new GeneralizedPareto(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(GPA.ParametersValid);

            var GPA2 = new GeneralizedPareto(double.PositiveInfinity, double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(GPA2.ParametersValid);

            var GPA3 = new GeneralizedPareto(100, 0, 1);
            Assert.IsFalse(GPA3.ParametersValid);
        }

        /// <summary>
        /// Checking ParametersToString().
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual("Location (ξ)",GPA.ParametersToString[0, 0] );
            Assert.AreEqual("Scale (α)", GPA.ParametersToString[1, 0]);
            Assert.AreEqual("Shape (κ)", GPA.ParametersToString[2, 0]);
            Assert.AreEqual("100", GPA.ParametersToString[0, 1]);
            Assert.AreEqual("10", GPA.ParametersToString[1, 1]);
            Assert.AreEqual("0", GPA.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GeneralizedPareto(100, 10, -0.1);
            var mom = dist.CentralMoments(1E-8);
            Assert.AreEqual(dist.Mean, mom[0], 1E-2);
            Assert.AreEqual(dist.StandardDeviation, mom[1], 1E-2);
            Assert.AreEqual(dist.Skewness, mom[2], 1E-2);
            Assert.AreEqual(dist.Kurtosis, mom[3], 1E-2);
        }

        /// <summary>
        /// Testing mean function.
        /// </summary>
        [TestMethod()]
        public void Test_Mean()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(110, GPA.Mean);

            var GPA2 = new GeneralizedPareto(100, 10, 0.9);
            Assert.AreEqual(105.26315, GPA2.Mean,  1e-04);

            var GPA3 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(105, GPA3.Mean); // Uniform(100,110).
        }

        /// <summary>
        /// Testing median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(106.93147, GPA.Median,  1e-04);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(105, GPA2.Median);
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(100, GPA.Mode);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.IsTrue(double.IsNaN(GPA2.Mode)); // A uniform density has no unique mode.
        }

        /// <summary>
        /// Testing standard deviation function.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(10, GPA.StandardDeviation);

            var GPA2 = new GeneralizedPareto(100, 10, 0.25);
            Assert.AreEqual(6.531972, GPA2.StandardDeviation,  1e-04);

            var GPA3 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(2.8867513459481291, GPA3.StandardDeviation, 1E-14);
        }

        /// <summary>
        /// Testing the skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(2, GPA.Skewness);

            var GPA2 = new GeneralizedPareto(100, 10, 0.3);
            Assert.AreEqual(0.932039, GPA2.Skewness,  1e-04);

            var GPA3 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(0, GPA3.Skewness);
        }

        /// <summary>
        /// Testing Kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(9, GPA.Kurtosis);

            var GPA2 = new GeneralizedPareto(100, 10, 0.24);
            Assert.AreEqual(3.786748, GPA2.Kurtosis,  1e-04);

            var GPA3 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(1.8, GPA3.Kurtosis);
        }

        /// <summary>
        /// Checking minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(100, GPA.Minimum);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(double.PositiveInfinity, GPA.Maximum);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(110, GPA2.Maximum);
        }

        /// <summary>
        /// Testing PDF method at different locations and parameters.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(0.1, GPA.PDF(100), 3E-17); // exp(-log(10)) differs by two ulps from the decimal literal.
            Assert.AreEqual(4.53999e-06, GPA.PDF(200),  1e-10);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(0,GPA2.PDF(200));
            Assert.AreEqual(0, GPA2.PDF(50));
        }

        /// <summary>
        /// Testing CDF method at different locations and parameters.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(0, GPA.CDF(100));
            Assert.AreEqual(0, GPA.CDF(0), 1e-04);
            Assert.AreEqual(0.999954, GPA.CDF(200),  1e-06);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(0, GPA2.CDF(50));
            Assert.AreEqual(0, GPA2.CDF(20));
            Assert.AreEqual(1, GPA2.CDF(200));
        }

        /// <summary>
        /// Testing InverseCDF method at different probabilities.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(100, GPA.InverseCDF(0));
            Assert.AreEqual(double.PositiveInfinity,GPA.InverseCDF(1));
            Assert.AreEqual(106.93147, GPA.InverseCDF(0.5),  1e-04);

            var GPA2 = new GeneralizedPareto(100, 10, 1);
            Assert.AreEqual(103,GPA2.InverseCDF(0.3));
        }

        /// <summary>
        /// Method-of-moments parameters from the mean, standard deviation and skewness, and their round trip.
        /// Reference shape from scipy.stats.genpareto (c = -kappa) skewness inverted with scipy.optimize.brentq;
        /// the scale and location follow from the closed-form mean and variance.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersFromMoments()
        {
            var GPA = new GeneralizedPareto();
            var parameters = GPA.ParametersFromMoments(new[] { 5d, 2d, 1.5, 9d });
            Assert.AreEqual(2.8012552638060155, parameters[0], 1E-7);
            Assert.AreEqual(2.428095363983795, parameters[1], 1E-7);
            Assert.AreEqual(0.10430980186759424, parameters[2], 1E-7);
            foreach (var input in new[] { new[] { 5d, 2d, 1.5, 9d }, new[] { 100d, 30d, 0.5, 3d }, new[] { 10d, 1d, -1d, 3d }, new[] { 0d, 4d, 4d, 30d } })
            {
                var moments = GPA.MomentsFromParameters(GPA.ParametersFromMoments(input));
                Assert.AreEqual(input[0], moments[0], 1E-6 * Math.Max(1d, Math.Abs(input[0])));
                Assert.AreEqual(input[1], moments[1], 1E-6 * input[1]);
                Assert.AreEqual(input[2], moments[2], 1E-6);
            }
            Assert.IsInstanceOfType(GPA, typeof(IMomentEstimation));
        }

        /// <summary>
        /// Shapes for skews below the skewness at kappa = 1/3 (about 0.861). Reference values from scipy.stats.genpareto
        /// (c = -kappa) skewness inverted with scipy.optimize.brentq; skew 0 is the uniform case kappa = 1.
        /// </summary>
        [TestMethod()]
        public void Test_SolveForKappa_LowSkew()
        {
            var GPA = new GeneralizedPareto();
            Assert.AreEqual(0.544487317258348, GPA.SolveForKappa(0.5), 1E-7);
            Assert.AreEqual(1d, GPA.SolveForKappa(0d), 1E-7);
            Assert.AreEqual(2.8424436304492624, GPA.SolveForKappa(-1d), 1E-7);
            Assert.AreEqual(30.247749902892412, GPA.SolveForKappa(-5d), 1E-6);
            Assert.AreEqual(-9.9, new GeneralizedPareto(0, 1, GPA.SolveForKappa(-9.9)).Skewness, 1E-6);
            Assert.IsTrue(double.IsNaN(GPA.SolveForKappa(-10d)));
        }

        /// <summary>
        /// For skews at or above the skewness at kappa = 1/3 the shape is bitwise the Brent solution on [-1/3, 1/3]
        /// with the solver tolerances of SolveForKappa (absolute 1E-15, the minimum relative tolerance).
        /// </summary>
        [TestMethod()]
        public void Test_SolveForKappa_UpperRangeUnchanged()
        {
            var GPA = new GeneralizedPareto();
            double third = 1d / 3d;
            double boundary = 2d * (1d - third) * Math.Sqrt(1d + 2d * third) / (1d + 3d * third);
            foreach (double skew in new[] { boundary, 0.87, 1d, 1.5, 2d, 3d, 5d, 9.99 })
            {
                double expected = Numerics.Mathematics.RootFinding.Brent.Solve((x) =>
                {
                    double k = 2d * (1d - x) * Math.Sqrt(1d + 2d * x) / (1d + 3d * x);
                    return k - skew;
                }, -(1d / 3d), 1d / 3d, 1E-15, Numerics.Mathematics.RootFinding.Brent.MinimumRelativeTolerance);
                Assert.AreEqual(expected, GPA.SolveForKappa(skew), 0.0, $"skew {skew}");
            }
        }

        /// <summary>
        /// Known-location maximum likelihood rejects excesses over the location that are too small for the scale bounds.
        /// </summary>
        [TestMethod()]
        public void Test_KnownLocationMLE_TinyExcesses()
        {
            var GPA = new GeneralizedPareto();
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.MLE(new[] { 0d, 1E-300, 2E-300 }, 0d));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => GPA.Estimate(new[] { 0d, 1E-300, 2E-300 }, ParameterEstimationMethod.MaximumLikelihood, new[] { 0d, double.NaN, double.NaN }));
            var fit = GPA.MLE(new[] { 10.2, 11d, 10.5, 13d, 12.1, 10.9 }, 10d);
            Assert.AreEqual(10d, fit[0]);
            Assert.IsTrue(fit[1] > 0d && fit[2] > -10d && fit[2] < 10d);
        }

        /// <summary>
        /// Reads one bit-pattern cell of the survival-quantile fixture.
        /// </summary>
        /// <param name="row">The fixture row.</param>
        /// <param name="column">The column name.</param>
        /// <returns>The binary64 value of the big-endian hexadecimal bit pattern.</returns>
        private static double FixtureBits(System.Collections.Generic.Dictionary<string, string> row, string column)
        {
            return BitConverter.Int64BitsToDouble(long.Parse(row[column], System.Globalization.NumberStyles.HexNumber, System.Globalization.CultureInfo.InvariantCulture));
        }

        /// <summary>
        /// The survival quantile reproduces all 300 reference cases of docs/distributions/oracles/gpd-survival-quantile.py
        /// bit for bit and lies within 2 + |kappa| ln(1/q) ulp of the 200-bit mpmath value (at most 8 ulp here), while
        /// InverseCDF(1 - q) differs from the reference in 268 cases and from the exact value by more than 1E6 ulp.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_ReferenceVectors()
        {
            int rows = 0, differFromInverseCdf = 0;
            double worstSurvivalUlp = 0, worstComplementUlp = 0;
            foreach (var row in DistributionOracle.Read("gpd-survival-quantile.csv"))
            {
                var gpa = new GeneralizedPareto(FixtureBits(row, "location"), FixtureBits(row, "scale"), FixtureBits(row, "kappa"));
                double q = FixtureBits(row, "q"), reference = FixtureBits(row, "reference"), exact = FixtureBits(row, "exact");
                double survival = gpa.InverseCCDF(q);
                double complement = gpa.InverseCDF(1d - q);
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference), BitConverter.DoubleToInt64Bits(survival), $"row {rows}");
                double ulp = Math.Abs(BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(exact) + 1) - exact);
                // The rounding of ln q is amplified by the condition number |kappa| ln(1/q) of the heavy tail.
                double survivalUlp = Math.Abs(survival - exact) / ulp;
                Assert.IsLessThanOrEqualTo(2d + Math.Abs(gpa.Kappa) * -Math.Log(q), survivalUlp, $"row {rows}");
                worstSurvivalUlp = Math.Max(worstSurvivalUlp, survivalUlp);
                worstComplementUlp = Math.Max(worstComplementUlp, Math.Abs(complement - exact) / ulp);
                if (survival != complement) differFromInverseCdf++;
                rows++;
            }
            Assert.AreEqual(300, rows);
            Assert.AreEqual(268, differFromInverseCdf);
            Assert.IsLessThanOrEqualTo(8d, worstSurvivalUlp);
            Assert.IsGreaterThan(1E6, worstComplementUlp);
        }

        /// <summary>
        /// Known values and limits: the exponential case ξ + α ln(1/q), the bounded and heavy tails against the
        /// closed form, the support limits at q = 0 and q = 1, a finite quantile for q below the double epsilon where
        /// InverseCDF(1 - q) returns the support limit, the extended-range fallback, and the list overload.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_KnownValuesAndLimits()
        {
            var exponential = new GeneralizedPareto(10d, 2d, 0d);
            Assert.AreEqual(10d + 2d * Math.Log(100d), exponential.InverseCCDF(0.01), 1E-14);
            // scipy.stats.genpareto.isf(1e-20, c=0.2, loc=10, scale=2) = 10 + 2 (1e-20^(-0.2) - 1) / 0.2 = 100000.
            var heavy = new GeneralizedPareto(10d, 2d, -0.2);
            Assert.AreEqual(100000d, heavy.InverseCCDF(1E-20), 1E-9);
            Assert.AreEqual(double.PositiveInfinity, heavy.InverseCDF(1d - 1E-20));
            // A bounded tail approaches its upper limit ξ + α/κ = 20.
            var bounded = new GeneralizedPareto(10d, 2d, 0.2);
            Assert.AreEqual(10d + 2d * (1d - Math.Pow(0.5, 0.2)) / 0.2, bounded.InverseCCDF(0.5), 1E-14);
            Assert.AreEqual(bounded.Maximum, bounded.InverseCCDF(0d), 0d);
            Assert.AreEqual(bounded.Minimum, bounded.InverseCCDF(1d), 0d);
            Assert.AreEqual(double.PositiveInfinity, heavy.InverseCCDF(0d));
            // The closed form overflows here; the extended-range kernel still gives the finite quantile.
            // alpha (q^(-2) - 1) / 2 = 1e-300 (1e400 - 1) / 2 = 5e99, while expm1(2 ln 1e200) overflows.
            var extreme = new GeneralizedPareto(0d, 1E-300, -2d);
            Assert.AreEqual(5E99, extreme.InverseCCDF(1E-200), 1E87);
            Assert.AreEqual(double.PositiveInfinity, extreme.InverseCDF(1d - 1E-200));
            // CCDF(InverseCCDF(q)) = q.
            foreach (double q in new[] { 1E-300, 1E-12, 0.3, 0.999 })
                Assert.AreEqual(q, heavy.CCDF(heavy.InverseCCDF(q)), 1E-13 * q);
            var list = heavy.InverseCCDF(new[] { 0.5, 0.01 });
            Assert.AreEqual(heavy.InverseCCDF(0.5), list[0], 0d);
            Assert.AreEqual(heavy.InverseCCDF(0.01), list[1], 0d);
        }

        /// <summary>
        /// The base survival quantile of a family without an override is bitwise InverseCDF(1 - q), and invalid
        /// probabilities or parameters throw.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_BaseAndInvalidInputs()
        {
            var normal = new Normal(3d, 2d);
            foreach (double q in new[] { 0d, 1E-17, 1E-9, 0.25, 0.5, 0.9, 1d })
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(normal.InverseCDF(1d - q)), BitConverter.DoubleToInt64Bits(normal.InverseCCDF(q)));
            var gpa = new GeneralizedPareto(10d, 2d, 0.1);
            foreach (double bad in new[] { -0.1, 1.1, double.NaN })
            {
                Assert.AreEqual("probability", Assert.Throws<ArgumentOutOfRangeException>(() => gpa.InverseCCDF(bad)).ParamName);
                Assert.AreEqual("probability", Assert.Throws<ArgumentOutOfRangeException>(() => normal.InverseCCDF(bad)).ParamName);
            }
            Assert.Throws<ArgumentNullException>(() => gpa.InverseCCDF((double[])null!));
            var invalid = new GeneralizedPareto(10d, -2d, 0.1);
            Assert.Throws<ArgumentOutOfRangeException>(() => invalid.InverseCCDF(0.5));
        }

        /// <summary>
        /// mpmath (400-bit) values of the GEV to Poisson-GPD mapping ξ − α expm1(κ ln λ)/κ and α λ^κ, as
        /// { rate, location, scale, shape, GPD location, GPD scale }, over heavy, bounded, exponential and
        /// ±1E-12 shapes and rates from 1E-6 to 1E6.
        /// </summary>
        private static readonly double[][] PoissonGevCases =
        {
            new[] { 3.7, 12.8911, 2.4337, -0.15, 9.999975062407715, 2.0000312593611573 },
            new[] { 3.7, 12.8911, 2.4337, 0.15, 9.373089605069575, 2.9614015592395635 },
            new[] { 0.3, 2.0, 1.5, 0.2, 3.604976858025329, 1.1790046283949343 },
            new[] { 250.0, 19.5, 9.6, -0.45, -0.05511241277857893, 0.8001994142496389 },
            new[] { 1E-06, -6.5, 0.25, -0.1, 0.9526792638374316, 0.9952679263837432 },
            new[] { 10.0, 1.0756, 0.25, 0.0, 0.49995372675148847, 0.25 },
            new[] { 10.0, 1.0756, 0.25, 1E-12, 0.4999537267508257, 0.25000000000057565 },
            new[] { 10.0, 1.0756, 0.25, -1E-12, 0.4999537267521512, 0.24999999999942435 },
            new[] { 1E6, 3.5, 0.004, 0.4, 0.9981135684904192, 1.0047545726038323 },
        };

        /// <summary>
        /// The GEV to Poisson-GPD mapping agrees with the mpmath values to within 4 ulp of the terms it combines,
        /// keeps the shape exactly (±1E-12 is not collapsed to zero), gives ξ − α ln λ with scale α at κ = 0, and
        /// returns the GEV location and scale unchanged at a unit rate.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedExtremeValue_MpmathCases()
        {
            const double u = 1.1102230246251565E-16;
            foreach (var c in PoissonGevCases)
            {
                double rate = c[0], location = c[1], scale = c[2], shape = c[3];
                var gpd = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(rate, location, scale, shape);
                double shift = Math.Abs(location - c[4]);
                Assert.AreEqual(c[4], gpd.Xi, 4 * u * (Math.Abs(location) + shift), $"location, rate {rate}, shape {shape}");
                Assert.AreEqual(c[5], gpd.Alpha, 4 * u * c[5], $"scale, rate {rate}, shape {shape}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(shape), BitConverter.DoubleToInt64Bits(gpd.Kappa));
                Assert.IsTrue(gpd.ParametersValid);
            }
            var exponential = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(10d, 1.0756, 0.25, 0d);
            Assert.AreEqual(1.0756 - 0.25 * Math.Log(10d), exponential.Xi, 0d);
            Assert.AreEqual(0.25, exponential.Alpha, 0d);
            var unit = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(1d, 2d, 1.5d, -0.2d);
            Assert.AreEqual(2d, unit.Xi, 0d);
            Assert.AreEqual(1.5d, unit.Alpha, 0d);
        }

        /// <summary>
        /// Mapping a GPD to its annual-maximum GEV and back recovers the GPD, and mapping a GEV to its GPD and forward
        /// recovers the GEV, over shapes from −0.45 to 0.4 (including 0 and ±1E-12) and rates from 1E-6 to 1E6, to within
        /// 1E-14 of the magnitudes combined.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedExtremeValue_RoundTrip()
        {
            foreach (double rate in new[] { 1E-6, 0.3, 1d, 3.7, 250d, 1E6 })
            {
                foreach (double kappa in new[] { -0.45, -0.15, -1E-12, 0d, 1E-12, 0.15, 0.4 })
                {
                    var gev = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(rate, 10d, 2d, kappa);
                    var gpd = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(rate, gev.Xi, gev.Alpha, gev.Kappa);
                    Assert.AreEqual(10d, gpd.Xi, 1E-14 * (10d + 2d * Math.Abs(gev.Xi - 10d)), $"GPD location, rate {rate}, kappa {kappa}");
                    Assert.AreEqual(2d, gpd.Alpha, 1E-14 * 2d, $"GPD scale, rate {rate}, kappa {kappa}");
                    Assert.AreEqual(kappa, gpd.Kappa, 0d);

                    var inverse = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(rate, 5d, 3d, kappa);
                    var forward = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(rate, inverse.Xi, inverse.Alpha, inverse.Kappa);
                    Assert.AreEqual(5d, forward.Xi, 1E-14 * (5d + 2d * Math.Abs(inverse.Xi - 5d)), $"GEV location, rate {rate}, kappa {kappa}");
                    Assert.AreEqual(3d, forward.Alpha, 1E-14 * 3d, $"GEV scale, rate {rate}, kappa {kappa}");
                }
            }
        }

        /// <summary>
        /// Textbook identity (Madsen et al., 1997; Coles, 2001, Section 4.2): above the recovered threshold the GEV CDF
        /// equals exp(−λ S(x)) of the recovered GPD survival S, and the GEV non-exceedance probability of the threshold is
        /// exp(−λ). Invalid arguments throw, naming the argument.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedExtremeValue_AnnualMaximumIdentityAndInvalidInputs()
        {
            foreach (double rate in new[] { 0.2, 1d, 3.7, 50d })
            {
                foreach (double kappa in new[] { -0.3, -1E-9, 0d, 1E-9, 0.25 })
                {
                    var gev = new GeneralizedExtremeValue(100d, 20d, kappa);
                    var gpd = GeneralizedPareto.FromPoissonGeneralizedExtremeValue(rate, 100d, 20d, kappa);
                    Assert.AreEqual(Math.Exp(-rate), gev.CDF(gpd.Xi), 1E-13, $"rate {rate}, kappa {kappa}");
                    foreach (double p in new[] { 0.001, 0.1, 0.5, 0.9, 0.999 })
                    {
                        double x = gpd.InverseCDF(p);
                        Assert.AreEqual(Math.Exp(-rate * gpd.CCDF(x)), gev.CDF(x), 1E-13, $"rate {rate}, kappa {kappa}, p {p}");
                    }
                }
            }
            Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(0d, 1d, 1d, 0.1)).ParamName);
            Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(double.PositiveInfinity, 1d, 1d, 0.1)).ParamName);
            Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(double.NaN, 1d, 1d, 0.1)).ParamName);
            Assert.AreEqual("location", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(2d, double.NaN, 1d, 0.1)).ParamName);
            Assert.AreEqual("scale", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(2d, 1d, 0d, 0.1)).ParamName);
            Assert.AreEqual("scale", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(2d, 1d, double.PositiveInfinity, 0.1)).ParamName);
            Assert.AreEqual("shape", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(2d, 1d, 1d, double.NaN)).ParamName);
            Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedPareto.FromPoissonGeneralizedExtremeValue(1E300, 1d, 1d, 5d)).ParamName);
        }

        /// <summary>
        /// The direct method of moments is continuous through κ = 0 and solves the skewness equation to near machine
        /// precision: the scale α = σ(1 + κ)√(1 + 2κ) and the location ξ = μ - α/(1 + κ) are used for every shape.
        /// </summary>
        /// <remarks>
        /// The mean, standard deviation and skewness of GPA(ξ = 10, α = 2, κ) were computed with Python mpmath 1.3.0 at
        /// 60 significant digits and rounded to double. The expected κ is the exact root (mpmath findroot) of the skewness
        /// equation at the rounded skew, and the expected α and ξ follow from it and the rounded mean and standard
        /// deviation. The tolerances are 2E-15 absolute for κ and 2E-15 relative for α and ξ.
        /// </remarks>
        [TestMethod]
        [DataRow(0d, 12.0, 2.0, 2.0, 0d, 2.0, 10.0)]
        [DataRow(1E-12, 11.999999999998, 1.999999999996, 1.999999999994, 1.0000148857165327e-12, 2.0, 10.0)]
        [DataRow(-1E-12, 12.000000000002, 2.000000000004, 2.000000000006, -1.0000148857115326e-12, 1.9999999999999998, 10.0)]
        [DataRow(1E-9, 11.999999998, 1.9999999960000001, 1.999999994, 1.0000000112255027e-09, 2.0, 10.0)]
        [DataRow(-1E-9, 12.000000002, 2.000000004, 2.000000006, -1.0000000062255026e-09, 1.9999999999999998, 10.0)]
        [DataRow(1E-6, 11.999998000002, 1.999996000007, 1.999994000015, 9.999999999843747e-07, 2.0, 10.0)]
        [DataRow(-1E-6, 12.000002000002, 2.000004000007, 2.000006000015, -9.999999999847885e-07, 2.0, 10.0)]
        [DataRow(1E-4, 11.999800019998, 1.999600069988002, 1.9994001499570127, 0.00010000000000000334, 2.0, 10.0)]
        [DataRow(-1E-4, 12.000200020002, 2.000400070012002, 2.0006001500430126, -9.999999999998938e-05, 2.0, 10.0)]
        [DataRow(1.0001E-4, 11.999800000002, 1.9996000300019992, 1.9994000899870012, 0.00010001000000000417, 2.0, 10.0)]
        [DataRow(-1.0001E-4, 12.000200040006002, 2.0004001100260065, 2.000600210073027, -0.00010000999999997046, 2.0000000000000004, 10.0)]
        [DataRow(1E-3, 11.998001998001998, 1.9960069880207136, 1.9940149571263739, 0.0009999999999999976, 2.0, 10.0)]
        [DataRow(-1E-3, 12.002002002002001, 2.0040070120207867, 2.006015043127128, -0.000999999999999969, 2.0000000000000004, 9.999999999999998)]
        [DataRow(0.01, 11.98019801980198, 1.9606882039142066, 1.9414582308696033, 0.009999999999999998, 2.0, 10.0)]
        [DataRow(-0.01, 12.02020202020202, 2.0407122112165874, 2.0615443063871717, -0.010000000000000021, 2.0, 10.0)]
        [DataRow(0.1, 11.818181818181818, 1.6597653257732305, 1.5167701592450753, 0.10000000000000002, 2.0, 10.0)]
        [DataRow(-0.1, 12.222222222222221, 2.4845199749997664, 2.8110568859997356, -0.09999999999999999, 2.0, 10.0)]
        [DataRow(0.2, 11.666666666666666, 1.4085904245475276, 1.1832159566199232, 0.2, 2.0, 10.0)]
        [DataRow(-0.2, 12.5, 3.2274861218395143, 4.6475800154489, -0.2, 2.0, 10.0)]
        [DataRow(0.25, 11.6, 1.3063945294843617, 1.0497813183356477, 0.25, 2.0, 10.0)]
        [DataRow(-0.25, 12.666666666666666, 3.7712361663282534, 7.0710678118654755, -0.25, 2.0, 10.0)]
        [DataRow(0.3, 11.538461538461538, 1.2162606385262997, 0.9320397314180486, 0.3, 2.0, 10.0)]
        [DataRow(0.5, 11.333333333333334, 0.9428090415820634, 0.565685424949238, 0.5, 2.0, 10.0)]
        [DataRow(1.0, 11.0, 0.5773502691896257, 0.0, 1.0, 1.9999999999999998, 10.0)]
        [DataRow(2.0, 10.666666666666666, 0.29814239699997197, -0.6388765649999399, 2.0, 2.0, 10.0)]
        public void Test_DirectMethodOfMoments_ContinuousThroughZeroShape(double shape, double mean, double sd, double skew,
            double expectedKappa, double expectedAlpha, double expectedXi)
        {
            double[] p = new GeneralizedPareto().DirectMethodOfMoments(new[] { mean, sd, skew });
            Assert.AreEqual(expectedKappa, p[2], 2E-15, $"kappa at shape {shape:R}");
            Assert.AreEqual(expectedAlpha, p[1], 2E-15 * expectedAlpha, $"alpha at shape {shape:R}");
            Assert.AreEqual(expectedXi, p[0], 2E-15 * expectedXi, $"xi at shape {shape:R}");
        }
    }
}
