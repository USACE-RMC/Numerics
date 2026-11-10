using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Generalized Extreme Value distribution algorithm.
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
    public class Test_GeneralizedExtremeValue
    {

        // Reference: "Flood Frequency Analysis", A.R. Rao & K.H. Hamed, CRC Press, 2000.
        // Table 7.1.2 White River near Nora, IN
        private double[] sample = new double[] { 23200d, 2950d, 10300d, 23200d, 4540d, 9960d, 10800d, 26900d, 23300d, 20400d, 8480d, 3150d, 9380d, 32400d, 20800d, 11100d, 7270d, 9600d, 14600d, 14300d, 22500d, 14700d, 12700d, 9740d, 3050d, 8830d, 12000d, 30400d, 27000d, 15200d, 8040d, 11700d, 20300d, 22700d, 30400d, 9180d, 4870d, 14700d, 12800d, 13700d, 7960d, 9830d, 12500d, 10700d, 13200d, 14700d, 14300d, 4050d, 14600d, 14400d, 19200d, 7160d, 12100d, 8650d, 10600d, 24500d, 14400d, 6300d, 9560d, 15800d, 14300d, 28700d };

        /// <summary>
        /// Verification of GEV Distribution fit with method of moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 218.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_MOM_Fit()
        {
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MethodOfMoments);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            double true_x = 11012d;
            double true_a = 6209.4d;
            double true_k = 0.0736d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x );
            Assert.IsLessThan(0.01d, (a - true_a) / true_a );
            Assert.IsLessThan(0.01d, (k - true_k) / true_k);
        }

        /// <summary>
        /// Verification of GEV Distribution fit with method of linear moments.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 218.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_Fit()
        {
            var sample = new double[] { 1953d, 1939d, 1677d, 1692d, 2051d, 2371d, 2022d, 1521d, 1448d, 1825d, 1363d, 1760d, 1672d, 1603d, 1244d, 1521d, 1783d, 1560d, 1357d, 1673d, 1625d, 1425d, 1688d, 1577d, 1736d, 1640d, 1584d, 1293d, 1277d, 1742d, 1491d };
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MethodOfLinearMoments);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            // Exact root of the tau3 equation for the sample L-moments (exact rational PWMs), mpmath 1.3.0, 50 digits.
            double true_x = 1543.8805732587272924d;
            double true_a = 218.02814344122047698d;
            double true_k = 0.1063218096908424871d;
            Assert.AreEqual(true_x, x, 1E-11d * true_x);
            Assert.AreEqual(true_a, a, 1E-11d * true_a);
            Assert.AreEqual(true_k, k, 1E-12d);
            var lmom = GEV.LinearMomentsFromParameters(GEV.GetParameters);
            Assert.AreEqual(1648.8064516129032258d, lmom[0], 1E-9d);
            Assert.AreEqual(138.23655913978494624d, lmom[1], 1E-9d);
            Assert.AreEqual(0.10339033967770315216d, lmom[2], 1E-13d);
            Assert.AreEqual(0.12781548866149448569d, lmom[3], 1E-13d);
        }

        /// <summary>
        /// The L-moment shape against the exact root of the GEV tau3 equation.
        /// </summary>
        /// <remarks>
        /// The reference shapes solve 2(1 - 3^-k)/(1 - 2^-k) - 3 = tau3 to 50 digits (mpmath 1.3.0). After the
        /// Newton refinement of PELGEV the shape is within 3E-14 absolute, and within 4E-14 relative below tau3 = -0.8,
        /// including 0.5 &lt; |tau3| &lt; 0.8, where PELGEV alone is off by up to 3.2E-7.
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_ShapeMatchesExactRoot()
        {
            var cases = new (double t3, double k)[]
            {
                (-0.999, 10.949387972183680871), (-0.99, 7.5831530320736617664), (-0.9, 4.1055902133998764511),
                (-0.85, 3.4680027529558994851), (-0.7, 2.3489840323468421931), (-0.6, 1.8713722836967761566),
                (-0.8, 3.0084904219595737341), (-0.5, 1.4946403940688676395), (-0.2, 0.68084783782142647358),
                (0.0, 0.28377552616996779592), (0.1, 0.11189497450865134315), (0.1699, 3.8904229635559562514E-5),
                (0.3, -0.19281434131450358283), (0.5, -0.45744624399490275071), (0.6, -0.57791882116241275303),
                (0.7, -0.69173389648592054131), (0.8, -0.79962394689806202906), (0.95, -0.95167113080354617872),
                (0.999, -0.99904421373046579835)
            };
            var gev = new GeneralizedExtremeValue();
            foreach (var (t3, k) in cases)
            {
                double[] parameters = gev.ParametersFromLinearMoments([0d, 1d, t3]);
                Assert.AreEqual(k, parameters[2], t3 < -0.8 ? 4E-14 * k : 3E-14, "tau3 = " + t3);
            }
        }

        /// <summary>
        /// The L-moment estimator against Hosking's PELGEV as ported in the R package lmomco.
        /// </summary>
        /// <remarks>
        /// Reference values from lmomco 2.5.7 (R 4.4.3), <c>pargev(vec2lmom(c(0, 1, t3)))</c>, which is PELGEV without
        /// the Newton refinement, so the shapes agree to the PELGEV error of 3.2E-7 and the scale and location to
        /// the corresponding 1E-5 relative (largest near tau3 = 1, where Gamma(1 + k) grows without bound).
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_MatchesLmomcoPargev()
        {
            var cases = new (double t3, double xi, double alpha, double k)[]
            {
                (-0.99, 1.0051819830031103, 0.00045669004816578825, 7.583153032073672),
                (-0.9, 1.0239854128079684, 0.15472353378220055, 4.1055902133998767),
                (-0.85, 1.0005535987780667, 0.34262205355987913, 3.4680027529559014),
                (-0.8, 0.95360409302564286, 0.56648627614288716, 3.0084901085603337),
                (-0.5, 0.37963185055183174, 1.7493892642370792, 1.4946402656328543),
                (0.0, -0.62175370892758763, 1.7656695434633298, 0.28377530000000001),
                (0.1699, -0.83272069418978945, 1.4427466831279077, 3.8749417599763212E-05),
                (0.5, -1.0432720270934954, 0.7487981913091305, -0.4574461192885802),
                (0.95, -1.0174342062106454, 0.050535458173498138, -0.95167095945192781)
            };
            var gev = new GeneralizedExtremeValue();
            foreach (var (t3, xi, alpha, k) in cases)
            {
                double[] parameters = gev.ParametersFromLinearMoments([0d, 1d, t3]);
                Assert.AreEqual(k, parameters[2], 3.2E-7, "kappa at tau3 = " + t3);
                Assert.AreEqual(alpha, parameters[1], 1E-5 * alpha, "alpha at tau3 = " + t3);
                Assert.AreEqual(xi, parameters[0], 1E-5, "xi at tau3 = " + t3);
            }
        }

        /// <summary>
        /// The estimator is continuous through the Gumbel limit.
        /// </summary>
        /// <remarks>
        /// At tau3 = ln(9/8)/ln 2 the exact shape is zero (the double nearest it has the root -1.6E-17, mpmath 1.3.0),
        /// so alpha = lambda2/ln 2 and xi = lambda1 - gamma alpha. On either side, the tau3 of shapes -1E-9 and 1E-9
        /// (as doubles, with exact roots -1.0000000151927717E-9 and 9.9999998317556224E-10) give the mpmath 1.3.0 scale
        /// and location to rounding, with no cancellation near the limit.
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_GumbelLimit()
        {
            var gev = new GeneralizedExtremeValue();
            double t3 = 0.16992500144231236291;
            double[] parameters = gev.ParametersFromLinearMoments([10d, 2d, t3]);
            Assert.AreEqual(0d, parameters[2], 1E-14);
            Assert.AreEqual(2d / Math.Log(2d), parameters[1], 5E-14);
            Assert.AreEqual(10d - 0.57721566490153286061 * 2d / Math.Log(2d), parameters[0], 5E-14);
            var near = new (double t3, double k, double alpha, double xi)[]
            {
                (0.16992500208495936, -1.0000000151927717E-9, 2.8853900791124345, 8.334507644131017),
                (0.16992500079966538, 9.9999998317556224E-10, 2.885390084443419, 8.334507646761514)
            };
            foreach (var c in near)
            {
                double[] back = gev.ParametersFromLinearMoments([10d, 2d, c.t3]);
                Assert.AreEqual(c.k, back[2], 1E-14);
                Assert.AreEqual(c.alpha, back[1], 1E-14 * c.alpha);
                Assert.AreEqual(c.xi, back[0], 1E-14 * c.xi);
            }
        }

        /// <summary>
        /// Parameters to L-moments and back.
        /// </summary>
        /// <remarks>The L-moments of each parameter set are mpmath 1.3.0 values (50 digits).</remarks>
        [TestMethod()]
        public void Test_GEV_LMOM_RoundTrip()
        {
            var cases = new (double xi, double alpha, double k, double l1, double l2, double t3)[]
            {
                (1000, 200, -0.1, 1137.2574042386387111, 153.39836423353711137, 0.23582477228781528172),
                (10, 2, 0.3, 10.683528691291485424, 1.1233198171466950048, -0.0089961072644468849217),
                (0, 1, 0.8, 0.085770286274696606868, 0.49555533522784088599, -0.25241290490184556669),
                (5, 1.5, -0.5, 7.3173615527165480819, 2.2025232711764534254, 0.53465397595792068585)
            };
            var gev = new GeneralizedExtremeValue();
            foreach (var c in cases)
            {
                double[] lmom = gev.LinearMomentsFromParameters([c.xi, c.alpha, c.k]);
                Assert.AreEqual(c.l1, lmom[0], 1E-13 * Math.Abs(c.l1));
                Assert.AreEqual(c.l2, lmom[1], 1E-13 * c.l2);
                Assert.AreEqual(c.t3, lmom[2], 1E-14);
                double[] back = gev.ParametersFromLinearMoments(lmom);
                Assert.AreEqual(c.k, back[2], 1E-13);
                Assert.AreEqual(c.alpha, back[1], 1E-12 * c.alpha);
                Assert.AreEqual(c.xi, back[0], 1E-12 * c.alpha);
            }
        }

        /// <summary>
        /// Invalid L-moments are rejected.
        /// </summary>
        [TestMethod()]
        public void Test_GEV_LMOM_InvalidMoments()
        {
            var gev = new GeneralizedExtremeValue();
            Assert.ThrowsExactly<ArgumentNullException>(() => gev.ParametersFromLinearMoments(null!));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, 1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, 0d, 0.1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, -1d, 0.1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, 1d, 1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, 1d, -1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([double.NaN, 1d, 0.1d]));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gev.ParametersFromLinearMoments([0d, 1d, double.NaN]));
        }

        /// <summary>
        /// The list quantile overload gives a growth curve bitwise equal to the scalar quantile.
        /// </summary>
        [TestMethod()]
        public void Test_GEV_GrowthCurveMatchesScalarInverseCDF()
        {
            var probabilities = new double[] { 0d, 1E-12, 0.01, 0.1, 0.5, 0.9, 0.98, 0.99, 0.999, 1 - 1E-12, 1d };
            foreach (double k in new[] { -0.3, -1E-9, 0d, 0.2, 1.5 })
            {
                var gev = new GeneralizedExtremeValue(1d, 0.25d, k);
                double[] curve = gev.InverseCDF(probabilities);
                for (int i = 0; i < probabilities.Length; i++)
                    Assert.AreEqual(BitConverter.DoubleToInt64Bits(gev.InverseCDF(probabilities[i])), BitConverter.DoubleToInt64Bits(curve[i]));
            }
        }

        /// <summary>
        /// Verification of GEV Distribution fit with method of maximum likelihood.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.1 page 219.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_MLE_Fit()
        {
            var GEV = new GeneralizedExtremeValue();
            GEV.Estimate(sample, ParameterEstimationMethod.MaximumLikelihood);
            double x = GEV.Xi;
            double a = GEV.Alpha;
            double k = GEV.Kappa;
            double true_x = 10849d;
            double true_a = 5745.6d;
            Assert.IsLessThan(0.01d, (x - true_x) / true_x);
            Assert.IsLessThan(0.01d, (a - true_a) / true_a);
            // The likelihood is flat in the shape near zero, so the shape is checked against the maximum
            // likelihood estimate kappa = 0.0050428 (scipy 1.18.1 genextreme, Nelder-Mead with xatol 1E-10,
            // log-likelihood -634.66351948) to 2E-4, a range over which the log-likelihood changes by under 1E-5.
            Assert.AreEqual(0.0050428d, k, 2E-4d);
            Assert.AreEqual(-634.66351948d, GEV.LogLikelihood(sample), 1E-5d);
        }

        /// <summary>
        /// Test the quantile function for the GEV Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.2 page 221.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_Quantile()
        {
            var GEV = new GeneralizedExtremeValue(10849d, 5745.6d, 0.005d);
            double q100 = GEV.InverseCDF(0.99d);
            double true_q100 = 36977d;
            Assert.IsLessThan(0.01d, (q100 - true_q100) / true_q100);
            double p = GEV.CDF(q100);
            double true_p = 0.99d;
            Assert.IsLessThan(0.01d, (p - true_p) / true_p);
        }

        /// <summary>
        /// Test the standard error for the GEV Distribution.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Reference: "Flood Frequency Analysis", A.R. Rao &amp; K.H. Hamed, CRC Press, 2000.
        /// </para>
        /// <para>
        /// Example 7.1.3 page 226.
        /// </para>
        /// </remarks>
        [TestMethod()]
        public void Test_GEV_StandardError()
        {
            double u = 10849d;
            double a = 5745.6d;
            double k = 0.005d;
            double true_dXdU = 1.0d;
            double true_dxdA = 4.5472d;
            double true_dxdK = -59861;
            double true_VarU = 664669d;
            double true_VarA = 346400d;
            double true_VarK = 0.007655d;
            double true_CovarUA = 176180d;
            double true_CovarUK = 23.977d;
            double true_CovarAK = 13.8574d;
            double true_QVar = 26445364d;
            double true_QSigma = 5142d;
            var GEV = new GeneralizedExtremeValue(u, a, k);
            var partials = GEV.QuantileGradient(0.99d).ToArray();
            var covar = GEV.ParameterCovariance(sample.Length, ParameterEstimationMethod.MaximumLikelihood);
            double qVar = GEV.QuantileVariance(0.99d, sample.Length, ParameterEstimationMethod.MaximumLikelihood);
            double qSigma = Math.Sqrt(qVar);
            Assert.IsLessThan(0.01d,(partials[0] - true_dXdU) / true_dXdU);
            Assert.IsLessThan(0.01d, (partials[1] - true_dxdA) / true_dxdA);
            Assert.IsLessThan(0.01d, (partials[2] - true_dxdK) / true_dxdK);
            Assert.IsLessThan(0.01d, (covar[0, 0] - true_VarU) / true_VarU);
            Assert.IsLessThan(0.01d, (covar[1, 1] - true_VarA) / true_VarA);
            Assert.IsLessThan(0.01d, (covar[2, 2] - true_VarK) / true_VarK);
            Assert.IsLessThan(0.01d, (covar[0, 1] - true_CovarUA) / true_CovarUA);
            Assert.IsLessThan(0.01d, (covar[0, 2] - true_CovarUK) / true_CovarUK);
            Assert.IsLessThan(0.01d, (covar[1, 2] - true_CovarAK) / true_CovarAK);
            Assert.IsLessThan(0.01d, (qVar - true_QVar) / true_QVar);
            Assert.IsLessThan(0.01d, (qSigma - true_QSigma) / true_QSigma);
        }

        /// <summary>
        /// Testing GEV constructor
        /// </summary>
        [TestMethod()]
        public void Test_Construction()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(100,GEV.Xi);
            Assert.AreEqual(10,GEV.Alpha);
            Assert.AreEqual(0,GEV.Kappa);

            var GEV2 = new GeneralizedExtremeValue(-100, 1, 1);
            Assert.AreEqual(-100,GEV2.Xi);
            Assert.AreEqual(1, GEV2.Alpha);
            Assert.AreEqual(1, GEV2.Kappa);
        }

        /// <summary>
        /// Testing bad parameters on GEV distribution.
        /// </summary>
        [TestMethod()]
        public void Test_InvalidParameters()
        {
            var GEV = new GeneralizedExtremeValue(double.NaN, double.NaN, double.NaN);
            Assert.IsFalse(GEV.ParametersValid);

            var GEV2 = new GeneralizedExtremeValue(double.PositiveInfinity,double.PositiveInfinity, double.PositiveInfinity);
            Assert.IsFalse(GEV2.ParametersValid);

            var GEV3 = new GeneralizedExtremeValue(100, 0, 1);
            Assert.IsFalse(GEV3.ParametersValid);
        }

        /// <summary>
        /// Testing Parameters to string.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersToString()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual("Location (ξ)", GEV.ParametersToString[0, 0]);
            Assert.AreEqual("Scale (α)", GEV.ParametersToString[1, 0]);
            Assert.AreEqual("Shape (κ)", GEV.ParametersToString[2, 0]);
            Assert.AreEqual("100", GEV.ParametersToString[0, 1]);
            Assert.AreEqual("10", GEV.ParametersToString[1, 1]);
            Assert.AreEqual("0", GEV.ParametersToString[2, 1]);
        }

        /// <summary>
        /// Compare analytical moments against numerical integration.
        /// </summary>
        [TestMethod()]
        public void Test_Moments()
        {
            var dist = new GeneralizedExtremeValue(100, 10, -0.1);
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
            var GEV = new GeneralizedExtremeValue();
            var true_val = 100 + 10 * Tools.Euler;
            Assert.AreEqual(true_val, GEV.Mean);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.9);
            Assert.AreEqual(100.42482, GEV2.Mean, 1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 10);
            Assert.AreEqual(-3628699d, GEV3.Mean, 1E-7); // 100 + (1 - Gamma(11)); this bounded shape has a mean.
        }

        /// <summary>
        /// Testing Median function.
        /// </summary>
        [TestMethod()]
        public void Test_Median()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(103.66512, GEV.Median,  1e-04);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.9);
            Assert.AreEqual(103.12196423496414, GEV2.Median, 1E-11); // R: 100 - 10*expm1(.9*log(log(2)))/.9.
        }

        /// <summary>
        /// Testing mode function.
        /// </summary>
        [TestMethod()]
        public void Test_Mode()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(100, GEV.Mode);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(110, GEV2.Mode); // Reverse-exponential endpoint mode.
        }

        /// <summary>
        /// Testing standard deviation method.
        /// </summary>
        [TestMethod()]
        public void Test_StandardDeviation()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(12.825498, GEV.StandardDeviation,  1e-05);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.49);
            Assert.AreEqual(9.280898, GEV2.StandardDeviation, 1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(10, GEV3.StandardDeviation);
        }

        /// <summary>
        /// Testing skew function.
        /// </summary>
        [TestMethod()]
        public void Test_Skewness()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(1.1395470994046487, GEV.Skewness);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.3);
            Assert.AreEqual(-0.0690175, GEV2.Skewness,  1e-03);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(-2, GEV3.Skewness);
        }

        /// <summary>
        /// Testing kurtosis function.
        /// </summary>
        [TestMethod()]
        public void Test_Kurtosis()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(3 + 12d / 5d, GEV.Kurtosis);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 0.24);
            Assert.AreEqual(2.7659607, GEV2.Kurtosis,  1e-04);

            var GEV3 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(9, GEV3.Kurtosis);
        }

        /// <summary>
        /// Testing minimum function.
        /// </summary>
        [TestMethod()]
        public void Test_Minimum()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.NegativeInfinity, GEV.Minimum);

            var GEV2 = new GeneralizedExtremeValue(100, 10, -5);
            Assert.AreEqual(98, GEV2.Minimum);
        }

        /// <summary>
        /// Testing maximum function.
        /// </summary>
        [TestMethod()]
        public void Test_Maximum()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.PositiveInfinity, GEV.Maximum);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(110, GEV2.Maximum);
        }

        /// <summary>
        /// Testing PDF method.
        /// </summary>
        [TestMethod()]
        public void Test_PDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(0,GEV.PDF(0));
            Assert.AreEqual(0,GEV.PDF(1));

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(1.67017007902456E-06, GEV2.PDF(0), 1e-10);
        }

        /// <summary>
        /// Testing CDF method.
        /// </summary>
        [TestMethod()]
        public void Test_CDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(0.367879, GEV.CDF(100),  1e-04);
            Assert.AreEqual(0.9999546, GEV.CDF(200),  1e-07);

            var GEV2 = new GeneralizedExtremeValue(100, 10, 1);
            Assert.AreEqual(0.367879, GEV2.CDF(100),  1e-05);
            Assert.AreEqual(1, GEV2.CDF(200));
        }

        /// <summary>
        /// Testing InverseCDF method.
        /// </summary>
        [TestMethod()]
        public void Test_InverseCDF()
        {
            var GEV = new GeneralizedExtremeValue();
            Assert.AreEqual(double.NegativeInfinity, GEV.InverseCDF(0));
            Assert.AreEqual(103.66512, GEV.InverseCDF(0.5),  1e-05);
            Assert.AreEqual(double.PositiveInfinity,GEV.InverseCDF(1));
        }

        /// <summary>
        /// The shape solved from the skewness on the Brent branches: skew below -2 (kappa above one) and skew above 3.7.
        /// Reference values from scipy.stats.genextreme (c = kappa) skewness inverted with scipy.optimize.brentq.
        /// </summary>
        [TestMethod()]
        public void Test_SolveForKappa_BrentBranches()
        {
            var GEV = new GeneralizedExtremeValue();
            var skews = new[] { -2.5, -5d, -10d, 4d, 5d, 7d, 9d };
            var kappas = new[] { 1.1585142783230908, 1.7430590358305955, 2.382554224931509, -0.21549801240424302,
                -0.23952844567736997, -0.26724241387319125, -0.2825237373052092 };
            for (int i = 0; i < skews.Length; i++)
            {
                double kappa = GEV.SolveForKappa(skews[i]);
                Assert.AreEqual(kappas[i], kappa, 1E-7, $"skew {skews[i]}");
                Assert.AreEqual(skews[i], new GeneralizedExtremeValue(0, 1, kappa).Skewness, 1E-6 * Math.Abs(skews[i]), $"skew {skews[i]}");
            }
            Assert.IsTrue(double.IsNaN(GEV.SolveForKappa(10d)));
        }

        /// <summary>
        /// The regression polynomials are kept for 0 &lt;= skew &lt;= 3.7, where the Type II polynomial is within 2E-4 of the
        /// exact shape; the jump to the Brent branch at 3.7 is below that bound.
        /// </summary>
        [TestMethod()]
        public void Test_SolveForKappa_PolynomialRange()
        {
            var GEV = new GeneralizedExtremeValue();
            foreach (double s in new[] { 1.15, 2d, 3d, 3.5, 3.7 })
            {
                double polynomial = 0.2858221d - 0.357983d * s + 0.116659d * Math.Pow(s, 2d) - 0.022725d * Math.Pow(s, 3d) + 0.002604d * Math.Pow(s, 4d) - 0.000161d * Math.Pow(s, 5d) + 0.000004d * Math.Pow(s, 6d);
                Assert.AreEqual(polynomial, GEV.SolveForKappa(s), 0.0, $"skew {s}");
            }
            foreach (double s in new[] { 0d, 0.5, 1.1 })
            {
                double polynomial = 0.277648d - 0.322016d * s + 0.060278d * Math.Pow(s, 2d) + 0.016759d * Math.Pow(s, 3d) - 0.005873d * Math.Pow(s, 4d) - 0.00244d * Math.Pow(s, 5d) - 0.00005d * Math.Pow(s, 6d);
                Assert.AreEqual(polynomial, GEV.SolveForKappa(s), 0.0, $"skew {s}");
            }
            Assert.AreEqual(GEV.SolveForKappa(3.7), GEV.SolveForKappa(3.7 + 1E-9), 2E-4);
        }

        /// <summary>
        /// Method-of-moments parameters reproduce the mean, standard deviation and skewness on the Brent branches.
        /// </summary>
        [TestMethod()]
        public void Test_ParametersFromMoments_RoundTrip()
        {
            var GEV = new GeneralizedExtremeValue();
            foreach (double skew in new[] { -3d, 6d })
            {
                var moments = GEV.MomentsFromParameters(GEV.ParametersFromMoments(new[] { 100d, 20d, skew, 0d }));
                Assert.AreEqual(100d, moments[0], 1E-6);
                Assert.AreEqual(20d, moments[1], 1E-6);
                Assert.AreEqual(skew, moments[2], 1E-6);
            }
            Assert.IsInstanceOfType(GEV, typeof(IMomentEstimation));
        }

        /// <summary>
        /// The scalar parameter setters set the same validity state as SetParameters.
        /// </summary>
        [TestMethod()]
        public void Test_ParameterSetters_Validity()
        {
            var GEV = new GeneralizedExtremeValue(100, 10, 0.1);
            GEV.Xi = double.NaN;
            Assert.IsFalse(GEV.ParametersValid);
            GEV.Xi = 100;
            Assert.IsTrue(GEV.ParametersValid);
            GEV.Kappa = double.PositiveInfinity;
            Assert.IsFalse(GEV.ParametersValid);
            GEV.Kappa = -0.2;
            Assert.IsTrue(GEV.ParametersValid);
        }

        /// <summary>
        /// Poisson-GPD to GEV reference cases: rate, GPD location, scale, shape (Hosking sign), then the
        /// expected GEV location and scale bits. The inputs are the twelve peaks-over-threshold reference
        /// cases; every expected value lies within 1.3 ulp of the 200-bit mpmath value of
        /// ξ + α expm1(−κ ln λ)/(−κ) and α λ^(−κ).
        /// </summary>
        private static readonly long[][] PoissonGpdCases =
        {
            new[] { 0x4024000000000000L, 0x3FE0000000000000L, 0x3FD0000000000000L, 0x0000000000000000L, 0x3FF135D8DDDAAA8BL, 0x3FD0000000000000L },
            new[] { 0x4024000000000000L, 0x3FE0000000000000L, 0x3FD0000000000000L, unchecked((long)0xBD619799812DEA11UL), 0x3FF135D8DDDAB060L, 0x3FD0000000001441L },
            new[] { 0x4024000000000000L, 0x3FE0000000000000L, 0x3FD0000000000000L, 0x3D619799812DEA11L, 0x3FF135D8DDDAA4B6L, 0x3FCFFFFFFFFFD77EL },
            new[] { 0x3FF0000000000000L, 0x4000000000000000L, 0x3FF8000000000000L, unchecked((long)0xBFC999999999999AUL), 0x4000000000000000L, 0x3FF8000000000000L },
            new[] { 0x3FD3333333333333L, 0x4000000000000000L, 0x3FF8000000000000L, unchecked((long)0xBFC999999999999AUL), 0x3FD9480F24FC70A8L, 0x3FF2DD33F50C9F3CL },
            new[] { 0x400D99999999999AL, 0x4024000000000000L, 0x4000000000000000L, unchecked((long)0xBFC3333333333333UL), 0x4029C83B9AE8D2E5L, 0x40037823C3587E89L },
            new[] { 0x400D99999999999AL, 0x4024000000000000L, 0x4000000000000000L, 0x3FC3333333333333L, 0x4028C07726EF563AL, 0x3FFA4C3DD1473222L },
            new[] { 0x406F400000000000L, 0x0000000000000000L, 0x3FE999999999999AL, unchecked((long)0xBFDCCCCCCCCCCCCDUL), 0x40338CDC795D1DFAL, 0x402331F9A06D67C8L },
            new[] { 0x3EB0C6F7A0B5ED8DL, 0x3FF0000000000000L, 0x3FF0000000000000L, unchecked((long)0xBFB999999999999AUL), unchecked((long)0xC019F3D40B1599D0UL), 0x3FD0137987DD704CL },
            new[] { 0x412E848000000000L, 0x3FF0000000000000L, 0x3FF0000000000000L, 0x3FD999999999999AL, 0x400BEB9DEE006F16L, 0x3F704E74CC73EE86L },
            new[] { 0x4029000000000000L, 0x4008000000000000L, 0x3FE6666666666666L, unchecked((long)0xBD719799812DEA11UL), 0x401312713A107C92L, 0x3FE666666666A49BL },
            new[] { 0x4029000000000000L, 0x4008000000000000L, 0x3FE6666666666666L, unchecked((long)0xBD819799812DEA11UL), 0x401312713A108664L, 0x3FE666666666E2D0L },
        };

        /// <summary>
        /// The Poisson-GPD mapping reproduces the reference cases bit for bit, keeps the GPD shape exactly
        /// (a shape of 5e-13 is not collapsed to zero), and gives the Gumbel case ξ + α ln λ at κ = 0.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedPareto_ReferenceCases()
        {
            foreach (var c in PoissonGpdCases)
            {
                double rate = BitConverter.Int64BitsToDouble(c[0]), location = BitConverter.Int64BitsToDouble(c[1]);
                double scale = BitConverter.Int64BitsToDouble(c[2]), shape = BitConverter.Int64BitsToDouble(c[3]);
                var gev = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(rate, location, scale, shape);
                Assert.AreEqual(c[4], BitConverter.DoubleToInt64Bits(gev.Xi), $"location, rate {rate}, shape {shape}");
                Assert.AreEqual(c[5], BitConverter.DoubleToInt64Bits(gev.Alpha), $"scale, rate {rate}, shape {shape}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(shape), BitConverter.DoubleToInt64Bits(gev.Kappa));
                Assert.IsTrue(gev.ParametersValid);
            }
            // A unit rate returns the GPD location and scale unchanged.
            var unit = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(1d, 2d, 1.5d, -0.2d);
            Assert.AreEqual(2d, unit.Xi, 0d);
            Assert.AreEqual(1.5d, unit.Alpha, 0d);
            // The exponential case is the Gumbel annual maximum.
            var gumbel = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(10d, 0.5d, 0.25d, 0d);
            Assert.AreEqual(0.5d + 0.25d * Math.Log(10d), gumbel.Xi, 0d);
        }

        /// <summary>
        /// Textbook identity (Madsen et al., 1997; Coles, 2001, Section 4.2): above the threshold the annual
        /// maximum CDF exp(−λ S(x)) of a Poisson(λ) count of GPD exceedances with survival S equals the mapped
        /// GEV CDF, across heavy, bounded and exponential tails and rates above and below one.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedPareto_AnnualMaximumIdentity()
        {
            foreach (double rate in new[] { 0.2, 1d, 3.7, 50d })
            {
                foreach (double kappa in new[] { -0.3, -1E-9, 0d, 1E-9, 0.25 })
                {
                    var gpd = new GeneralizedPareto(10d, 2d, kappa);
                    var gev = GeneralizedExtremeValue.FromPoissonGeneralizedPareto(rate, 10d, 2d, kappa);
                    foreach (double p in new[] { 0.001, 0.1, 0.5, 0.9, 0.999 })
                    {
                        double x = gpd.InverseCDF(p);
                        double annual = Math.Exp(-rate * gpd.CCDF(x));
                        Assert.AreEqual(annual, gev.CDF(x), 1E-13, $"rate {rate}, kappa {kappa}, p {p}");
                    }
                    // The location is the level exceeded once per year on average: rate * S(location) = 1.
                    if (rate >= 1d) Assert.AreEqual(1d, rate * gpd.CCDF(gev.Xi), 1E-13);
                }
            }
        }

        /// <summary>
        /// The survival quantile is within 3 + |kappa g| ulp (g = log(-log1p(-q))) of the 400-bit mpmath values of
        /// docs/distributions/oracles/gev-survival-quantile.py over nine shapes and survival probabilities from 0.89
        /// down to 1E-300, while InverseCDF(1 - q) is off by more than 1E10 ulp on the same grid.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_MpmathGrid()
        {
            int rows = 0;
            double worstSurvival = 0, worstComplement = 0;
            foreach (var row in DistributionOracle.Read("gev-survival-quantile.csv"))
            {
                double kappa = BitConverter.Int64BitsToDouble(long.Parse(row["kappa"], System.Globalization.NumberStyles.HexNumber, System.Globalization.CultureInfo.InvariantCulture));
                double q = BitConverter.Int64BitsToDouble(long.Parse(row["q"], System.Globalization.NumberStyles.HexNumber, System.Globalization.CultureInfo.InvariantCulture));
                double exact = BitConverter.Int64BitsToDouble(long.Parse(row["exact"], System.Globalization.NumberStyles.HexNumber, System.Globalization.CultureInfo.InvariantCulture));
                var gev = new GeneralizedExtremeValue(10d, 2d, kappa);
                double ulp = Math.Abs(BitConverter.Int64BitsToDouble(BitConverter.DoubleToInt64Bits(exact) + 1) - exact);
                double survival = gev.InverseCCDF(q);
                double g = Math.Log(-Tools.Log1p(-q));
                double error = Math.Abs(survival - exact) / ulp;
                Assert.IsLessThanOrEqualTo(3d + Math.Abs(kappa * g), error, $"kappa {kappa}, q {q}");
                worstSurvival = Math.Max(worstSurvival, error);
                double complement = gev.InverseCDF(1d - q);
                if (Tools.IsFinite(complement)) worstComplement = Math.Max(worstComplement, Math.Abs(complement - exact) / ulp);
                rows++;
            }
            Assert.AreEqual(2196, rows);
            Assert.IsGreaterThan(1E10, worstComplement);
        }

        /// <summary>
        /// Survival quantile limits and identities: q = 0 and q = 1 give the support limits, CCDF(InverseCCDF(q)) = q,
        /// a q below the double epsilon gives a finite heavy-tail quantile where InverseCDF(1 - q) returns the upper
        /// limit, the list overload agrees with the scalar member, and invalid input throws.
        /// </summary>
        [TestMethod]
        public void Test_InverseCCDF_LimitsAndInvalidInputs()
        {
            var heavy = new GeneralizedExtremeValue(10d, 2d, -0.2);
            var bounded = new GeneralizedExtremeValue(10d, 2d, 0.2);
            Assert.AreEqual(heavy.Maximum, heavy.InverseCCDF(0d));
            Assert.AreEqual(heavy.Minimum, heavy.InverseCCDF(1d), 0d);
            Assert.AreEqual(bounded.Maximum, bounded.InverseCCDF(0d), 0d);
            Assert.AreEqual(bounded.Minimum, bounded.InverseCCDF(1d));
            foreach (double q in new[] { 1E-200, 1E-20, 1E-9, 0.3, 0.9 })
                Assert.AreEqual(q, heavy.CCDF(heavy.InverseCCDF(q)), 1E-13 * q);
            // Near the finite upper bound x resolves the survival probability only to about ulp(x) / alpha.
            foreach (double q in new[] { 1E-6, 0.3, 0.9 })
                Assert.AreEqual(q, bounded.CCDF(bounded.InverseCCDF(q)), 1E-9 * q);
            // 10 + 2 (1 - (1e-20)^(-0.2)) / (-0.2) = 10 + 10 (1e4 - 1) = 100000 to first order in q.
            Assert.AreEqual(100000d, heavy.InverseCCDF(1E-20), 1E-9);
            Assert.AreEqual(double.PositiveInfinity, heavy.InverseCDF(1d - 1E-20));
            var list = heavy.InverseCCDF(new[] { 0.5, 1E-6 });
            Assert.AreEqual(heavy.InverseCCDF(0.5), list[0], 0d);
            Assert.AreEqual(heavy.InverseCCDF(1E-6), list[1], 0d);
            foreach (double bad in new[] { -0.1, 1.1, double.NaN })
                Assert.AreEqual("probability", Assert.Throws<ArgumentOutOfRangeException>(() => heavy.InverseCCDF(bad)).ParamName);
            Assert.Throws<ArgumentOutOfRangeException>(() => new GeneralizedExtremeValue(10d, -1d, 0d).InverseCCDF(0.5));
        }

        /// <summary>
        /// Non-positive or non-finite rate and scale, non-finite location and shape, and a mapping that
        /// overflows the double range throw ArgumentOutOfRangeException naming the argument.
        /// </summary>
        [TestMethod]
        public void Test_FromPoissonGeneralizedPareto_InvalidInputs()
        {
            foreach (double bad in new[] { 0d, -1d, double.NaN, double.PositiveInfinity })
                Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedExtremeValue.FromPoissonGeneralizedPareto(bad, 0d, 1d, 0d)).ParamName);
            foreach (double bad in new[] { 0d, -1d, double.NaN, double.PositiveInfinity })
                Assert.AreEqual("scale", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedExtremeValue.FromPoissonGeneralizedPareto(2d, 0d, bad, 0d)).ParamName);
            Assert.AreEqual("location", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedExtremeValue.FromPoissonGeneralizedPareto(2d, double.NaN, 1d, 0d)).ParamName);
            Assert.AreEqual("shape", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedExtremeValue.FromPoissonGeneralizedPareto(2d, 0d, 1d, double.NegativeInfinity)).ParamName);
            Assert.AreEqual("rate", Assert.Throws<ArgumentOutOfRangeException>(() => GeneralizedExtremeValue.FromPoissonGeneralizedPareto(1E300, 0d, 1d, -5d)).ParamName);
        }

        /// <summary>
        /// The L-moments are continuous through κ = 0 and equal the Gumbel values there, and
        /// <see cref="GeneralizedExtremeValue.ParametersFromLinearMoments"/> inverts them.
        /// </summary>
        /// <remarks>
        /// Reference values for ξ = 10 and α = 2 computed with Python mpmath 1.3.0 at 50 significant digits from
        /// λ1 = ξ + α(1 - Γ(1 + κ))/κ, λ2 = αΓ(1 + κ)(1 - 2^-κ)/κ and the τ3 and τ4 ratios, with the exact Gumbel limits
        /// ξ + γα, α ln 2, 2 ln 3/ln 2 - 3 and 16 - 10 ln 3/ln 2 at κ = 0. The tolerances are 2E-15 relative for λ1 and
        /// λ2 and 1E-14 absolute for τ3 and τ4.
        /// </remarks>
        [TestMethod]
        [DataRow(0d, 11.154431329803066, 1.3862943611198906, 0.16992500144231237, 0.1503749927884382)]
        [DataRow(1E-12, 11.154431329801088, 1.38629436111861, 0.16992500144166972, 0.1503749927881857)]
        [DataRow(-1E-12, 11.154431329805044, 1.3862943611211713, 0.169925001442955, 0.15037499278869068)]
        [DataRow(1E-9, 11.154431327824954, 1.3862943598392468, 0.16992500079966538, 0.15037499253593725)]
        [DataRow(-1E-9, 11.154431331781177, 1.3862943624005344, 0.16992500208495936, 0.15037499304093913)]
        [DataRow(1E-6, 11.15442935169289, 1.3862930804778149, 0.16992435879548176, 0.15037474028788952)]
        [DataRow(-1E-6, 11.154433307916872, 1.3862956417654855, 0.16992564408946517, 0.1503752452897781)]
        [DataRow(1E-3, 11.152455030809067, 1.3850154749259203, 0.16928251552522589, 0.15012288732613316)]
        [DataRow(-1E-3, 11.156411258717297, 1.3875767662296468, 0.17056780955634154, 0.15062788953119546)]
        [DataRow(0.1, 10.972984602662537, 1.2741823011427935, 0.10724765961205163, 0.12892665962377894)]
        [DataRow(-0.1, 11.372574042386388, 1.5339836423353712, 0.2358247722878153, 0.17974345124238944)]
        [DataRow(0.5, 10.455092298188967, 1.0382794271800315, -0.11396711386776502, 0.10536947527156275)]
        [DataRow(-0.5, 13.089815403622064, 2.936697694901938, 0.5346539759579206, 0.39779793207587183)]
        public void Test_LinearMomentsFromParameters_ContinuousThroughZeroShape(double kappa, double l1, double l2, double t3, double t4)
        {
            var gev = new GeneralizedExtremeValue(10d, 2d, kappa);
            double[] lmom = gev.LinearMomentsFromParameters(new[] { 10d, 2d, kappa });
            Assert.AreEqual(l1, lmom[0], 2E-15 * l1);
            Assert.AreEqual(l2, lmom[1], 2E-15 * l2);
            Assert.AreEqual(t3, lmom[2], 1E-14);
            Assert.AreEqual(t4, lmom[3], 1E-14);

            double[] p = gev.ParametersFromLinearMoments(lmom);
            Assert.AreEqual(kappa, p[2], 2E-14);
            Assert.AreEqual(2d, p[1], 2E-14 * 2d);
            Assert.AreEqual(10d, p[0], 2E-15 * 10d);
        }
    }
}
