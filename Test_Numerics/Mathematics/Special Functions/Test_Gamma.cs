using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.SpecialFunctions;
using System;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Unit tests for the Gamma class
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> Description: </b>
    /// Some of these special functions are tested directly with MS Excel and with values provided by:
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// Others must be tested indirectly through the continuous distribution classes.
    /// </remarks>
    [TestClass]
    public class Test_Gamma
    {
        /// <summary>
        /// Test the Gamma Function (which uses the Stirling approximation)
        /// </summary>
        [TestMethod]
        public void Test_Function()
        {
            var testX = new double[] { 1d, 1.2d, 1.4d, 1.6d, 1.8d, 2d, 0.2d, 0.4d, 0.6d, 0.8d, -0.2d, -0.4d, -0.6d, -0.8d, 10d, 20d, 30d, 40d, 50d, 60d };
            var testValid = new double[] { 1d, 0.918168742399761d, 0.887263817503075d, 0.89351534928769d, 0.931383770980243d, 1d, 4.5908437119988d, 2.21815954375769d, 1.48919224881282d, 1.1642297137253d, -5.82114856862652d, -3.72298062203204d, -3.69693257292948d, -5.73855463999851d, 362880d, 121645100408832000d, 8.8417619937397E+30d, 2.03978820811974E+46d, 6.08281864034268E+62d, 1.3868311854569E+80d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Function(testX[i]);
                Assert.IsLessThan(0.01d, Math.Abs(testValid[i] - testResults[i]) / testValid[i]);
            }
        }

        /// <summary>
        /// Arguments one ulp below 1 and 2 are reduced upward through the recurrence, where rounding can carry the
        /// reduced argument to exactly 3; the result must still be close to Gamma(1) = Gamma(2) = 1.
        /// Reference values from scipy.special.gamma.
        /// </summary>
        [TestMethod]
        public void Test_Function_JustBelowIntegers()
        {
            var testX = new double[] { 1d - Math.Pow(2, -52), 2d - Math.Pow(2, -52), 0.5d - Math.Pow(2, -54), 1.5d - Math.Pow(2, -52), 3d - Math.Pow(2, -51) };
            var testValid = new double[] { 1.0000000000000004d, 1.0000000000000002d, 1.7724538509055165d, 0.8862269254527582d, 1.9999999999999993d };
            for (int i = 0; i < testX.Length; i++)
            {
                Assert.AreEqual(testValid[i], Gamma.Function(testX[i]), 1E-14 * testValid[i]);
            }
        }

        /// <summary>
        /// Test the Lanczos approximation to the Gamma Function
        /// </summary>
        [TestMethod]
        public void Test_Lanczos()
        {
            var testX = new double[] { 1d, 1.2d, 1.4d, 1.6d, 1.8d, 2d, 0.2d, 0.4d, 0.6d, 0.8d, -0.2d, -0.4d, -0.6d, -0.8d, 10d, 20d, 30d, 40d, 50d, 60d };
            var testValid = new double[] { 1d, 0.918168742399761d, 0.887263817503075d, 0.89351534928769d, 0.931383770980243d, 1d, 4.5908437119988d, 2.21815954375769d, 1.48919224881282d, 1.1642297137253d, -5.82114856862652d, -3.72298062203204d, -3.69693257292948d, -5.73855463999851d, 362880d, 121645100408832000d, 8.8417619937397E+30d, 2.03978820811974E+46d, 6.08281864034268E+62d, 1.3868311854569E+80d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Lanczos(testX[i]);
                Assert.IsLessThan(0.01d, Math.Abs(testValid[i] - testResults[i]) / testValid[i]);
            }
        }

        /// <summary>
        /// Test the digamma function
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Digamma()
        {
            double[] x = { 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };

            double[] expected =
            {
                -0.5772156649015329,
                -0.4237549404110768,
                -0.2890398965921883,
                -0.1691908888667997,
                -0.06138454458511615,
                 0.03648997397857652,
                 0.1260474527734763,
                 0.208547874873494,
                 0.2849914332938615,
                 0.3561841611640597,
                 0.4227843350984671,
            };

            for (int i = 0; i < x.Length; i++)
            {
                double xi = x[i];
                double expectedi = expected[i];
                double actual = Gamma.Digamma(xi);

                Assert.AreEqual(expectedi, actual, 1e-10);
            }
        }

        /// <summary>
        /// Test the trigamma function
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath at 30 digits, polygamma(1, x), and
        /// cross-checked with scipy 1.17 polygamma(1, x). The grid covers tiny arguments, both sides
        /// of the shift threshold (20) and large arguments. The relative tolerance of 1E-15 is about
        /// five ulps.
        /// </remarks>
        [TestMethod]
        public void Test_Trigamma()
        {
            var testX = new double[] { 1.2d, 1.4d, 1.6d, 1.8d, 0.2d, 0.4d, 0.6d, 0.8d,
                1E-8, 1E-4, 0.5, 1.0, 2.5, 5.0, 10.0, 19.9, 20.5, 100.0, 1E6 };
            var testValid = new double[] { 1.2673772054237793, 1.0253565905295976, 0.8584318931245798, 0.7369741375017003,
                26.267377205423777, 7.275356590529597, 3.6362096709023577, 2.2994741375017,
                1.0000000000000002E+16, 100000001.64469367, 4.934802200544679, 1.6449340668482264,
                0.49035775610023485, 0.22132295573711533, 0.10516633568168575, 0.05153498898307066,
                0.04998959242996099, 0.010050166663333571, 1.0000005000001667E-06 };

            for (int i = 0; i < testValid.Length; i++)
            {
                double actual = Gamma.Trigamma(testX[i]);
                Assert.AreEqual(testValid[i], actual, 1E-15 * testValid[i], $"x = {testX[i]}");
            }
        }

        /// <summary>
        /// Test the function that returns the natural logarithm of the gamma function
        /// </summary>
        /// <remarks>
        /// This function was tested against values calculated from MS Excel's "GAMMALN( )" function
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// <see href = "https://learn.microsoft.com/en-us/office/vba/api/excel.worksheetfunction.gammaln" />
        ///  </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_LogGamma()
        {
            var testX = new double[] {1.2d, 1.4d, 1.6d, 1.8d, 0.2d, 0.4d, 0.6d, 0.8d, 10d, 20d, 30d, 40d, 50d, 60d};
            var testValid = new double[] { -0.08537409d, -0.119612914d, -0.112591766d, -0.071083873d,
                1.524063822d, 0.796677818d, 0.398233858d, 0.152059678d, 12.80182748d, 39.33988419d, 71.25703897d,
                106.6317603d, 144.5657439d, 184.5338289d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.LogGamma(testX[i]);
                Assert.IsLessThan(0.01d,Math.Abs(testValid[i] - testResults[i]) / testValid[i]);
            }
        }

        /// <summary>
        /// LogGamma returns ln|Γ(z)| for negative non-integer z instead of NaN.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy 1.17 gammaln(-0.5) = 1.2655121234846454, gammaln(-2.5) and gammaln(-0.1).
        /// The absolute tolerance of 1E-14 is a few ulps of the reflection terms.
        /// </remarks>
        [TestMethod]
        public void Test_LogGamma_NegativeArguments()
        {
            Assert.AreEqual(1.2655121234846454, Gamma.LogGamma(-0.5), 1E-14);
            Assert.AreEqual(-0.05624371649767401, Gamma.LogGamma(-2.5), 1E-14);
            Assert.AreEqual(2.368961332728789, Gamma.LogGamma(-0.1), 1E-14);
        }

        /// <summary>
        /// LogGamma is accurate to a few ulps relative to 40-digit reference values on every interval of the
        /// approximation, at integers and half-integers, at both sides of the zeros at 1 and 2, and at tiny,
        /// large and negative non-integer arguments.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.3.0 at 40 significant digits, log|gamma(x)| at the
        /// double nearest each argument, and rounded to double. The relative tolerance of 1E-15 is about
        /// five ulps, so the test checks full relative accuracy on (0.5, 4) and next to the zeros at 1 and 2,
        /// where ln Γ is small and an absolute error bound would hide a loss of significant digits.
        /// </remarks>
        [TestMethod]
        public void Test_LogGamma_ReferenceValues()
        {
            var testX = new double[]
            {
                // Integers and half-integers
                3d, 4d, 5d, 10d, 20d, 30d, 60d, 0.5d, 1.5d, 2.5d, 10.5d, 99.5d, 199.5d,
                // Next to the zeros at 1 and 2
                1d + 1E-3, 1d - 1E-3, 1d + 1E-8, 1d - 1E-8, 1d + 1E-15, 1d - 1E-15,
                2d + 1E-3, 2d - 1E-3, 2d + 1E-8, 2d - 1E-8, 2d + 1E-14, 2d - 1E-14,
                // Tiny, intermediate and large arguments
                1E-300, 1E-10, 0.01d, 0.2d, 0.6d, 1.2d, 3.7d, 7.3d, 12.5d, 37.1d, 1E5, 1E15, 1E300,
                // Negative non-integers
                -0.5d, -2.5d, -0.1d, -1E-10, -10.3d, -100.7d
            };
            var testValid = new double[]
            {
                0.6931471805599453, 1.791759469228055, 3.1780538303479458, 12.801827480081469, 39.339884187199495,
                71.25703896716801, 184.53382886144948, 0.5723649429247001, -0.12078223763524522, 0.2846828704729192,
                13.940625219403763, 356.8353828236131, 855.2863892734525,
                -0.0005763935982833062, 0.0005780385328913802, -5.772156531688512E-09, 5.772156760265731E-09,
                -6.408381213479996E-16, 5.767543092132013E-16,
                0.000423106734800117, -0.0004224618006921073, 4.227843357536779E-09, -4.227843293043373E-09,
                4.318341110149868E-15, -4.224464129494371E-15,
                690.7755278982137, 23.025850929882736, 4.599479878042022, 1.5240638224307845, 0.39823385806923495, -0.08537409000331583,
                1.428072326665388, 7.147892523022248, 18.734347511936445, 96.07956574816038, 1051287.7089736569,
                3.3538776394910668E+16, 6.897755278982137E+302,
                1.2655121234846454, -0.056243716497674054, 2.3689613327287886, 23.025850929998178,
                -14.457515440024208, -365.61225576180533
            };

            for (int i = 0; i < testX.Length; i++)
            {
                // The reflection formula for negative arguments subtracts terms of order one, so the
                // tolerance there is relative to the larger of 1 and the result.
                double scale = testX[i] < 0d ? Math.Max(1d, Math.Abs(testValid[i])) : Math.Abs(testValid[i]);
                double actual = Gamma.LogGamma(testX[i]);
                Assert.AreEqual(testValid[i], actual, 1E-15 * scale, $"x = {testX[i]:R}");
            }

            // The zeros at 1 and 2 are exact
            Assert.AreEqual(0d, Gamma.LogGamma(1d));
            Assert.AreEqual(0d, Gamma.LogGamma(2d));
        }

        /// <summary>
        /// The shape-derivative form of LogIncompleteRatio returns the same logarithm as the value form, bit for bit, and a
        /// derivative that matches high-precision references in every region: the small-x series of Q, the power series of
        /// P, the continued fraction of Q and Temme's expansion, for both tails.
        /// </summary>
        /// <remarks>
        /// References computed with Python mpmath 1.4.1 at 90 significant digits: the series and continued fraction of the
        /// smaller ratio, the larger as log1p(-exp(·)), and mpmath.diff in the shape at fixed x.
        /// </remarks>
        [TestMethod]
        public void Test_LogIncompleteRatio_ShapeDerivative()
        {
            // { a, x, upper, ln P or ln Q, d/da }
            var table = new double[,]
            {
                { 0.3, 0.2, 0, -0.41930028795632035, -1.556456032228477 },
                { 0.3, 0.2, 1, -1.071503254376882, 2.988030363627726 },
                { 5, 2, 0, -2.944031732572081, -1.0809114750831792 },
                { 5, 9, 1, -2.901083376074, 0.839148885416154 },
                { 100, 105, 1, -1.2047909468838265, 0.1151664255845481 },
                { 100, 95, 0, -1.1477285508168387, -0.11340905103049531 },
                { 20000, 21000, 0, -1.714866608576869E-12, -8.536083365636196E-14 },
                { 0.5, 30, 1, -32.28898703270474, 5.3965135605268 },
                { 50, 10, 0, -43.13137931408243, -1.6241118212333785 },
                { 0.001, 1E-05, 1, -4.521108920080178, 994.6164157622044 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double a = table[i, 0], x = table[i, 1], log = table[i, 3], derivative = table[i, 4];
                bool upper = table[i, 2] == 1;
                double actual = Gamma.LogIncompleteRatio(a, x, upper, out double actualDerivative);
                Assert.AreEqual(Gamma.LogIncompleteRatio(a, x, upper), actual, 0d, $"value at a = {a:R}, x = {x:R}");
                Assert.AreEqual(log, actual, 1E-15 * Math.Max(1d, Math.Abs(log)), $"log at a = {a:R}, x = {x:R}");
                Assert.AreEqual(derivative, actualDerivative, 1E-14 * Math.Abs(derivative), $"derivative at a = {a:R}, x = {x:R}");
            }
        }

        /// <summary>
        /// RatioHalf returns Γ(x) / Γ(x + 1/2) to a few ulps below and above 1, including large x where Γ overflows.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.4.1 at 50 significant digits. The tolerance is 2E-15 relative,
        /// above the 8-ulp worst case measured up to 1E8.
        /// </remarks>
        [TestMethod]
        public void Test_RatioHalf_ReferenceValues()
        {
            // { x, Γ(x) / Γ(x + 1/2) }
            var table = new double[,]
            {
                { 0.001, 564.9713307845761 },
                { 0.25, 2.958675119188639 },
                { 0.75, 1.3519564801345694 },
                { 1.0, 1.1283791670955126 },
                { 2.5, 0.6646701940895685 },
                { 150.0, 0.08171772770686318 },
                { 10000.0, 0.010000125000781202 },
                { 100000000.0, 0.000100000000125 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double x = table[i, 0], expected = table[i, 1];
                Assert.AreEqual(expected, Gamma.RatioHalf(x), 2E-15 * expected, $"x = {x:R}");
            }
        }

        /// <summary>
        /// LogGammaOnePlus keeps full relative accuracy for small increments of either sign, where LogGamma(1 + a) loses
        /// the digits of a in the rounding of 1 + a, and matches LogGamma away from one.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.4.1 at 50 significant digits as loggamma(1 + a); the value at
        /// 1E-300 is -γ a. The relative tolerance of 1E-15 is about five ulps.
        /// </remarks>
        [TestMethod]
        public void Test_LogGammaOnePlus_ReferenceValues()
        {
            // { a, ln Γ(1 + a) }
            var table = new double[,]
            {
                { 1E-300, -5.772156649015329E-301 },
                { 1E-10, -5.7721566481928616E-11 },
                { -1E-08, 5.772156731262032E-09 },
                { 0.001, -0.0005763935982833696 },
                { 0.25, -0.09827183642181316 },
                { -0.3, 0.2608672465316665 },
                { -0.3203125, 0.2862428123190013 },
                { 0.5, -0.12078223763524522 },
                { -0.5, 0.5723649429247001 },
                { 0.75, -0.08440112102048555 },
                { -0.9, 2.252712651734206 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double a = table[i, 0], expected = table[i, 1];
                Assert.AreEqual(expected, Gamma.LogGammaOnePlus(a), 1E-15 * Math.Abs(expected), $"a = {a:R}");
            }
            Assert.AreEqual(0d, Gamma.LogGammaOnePlus(0d));
            Assert.AreEqual(0d, Gamma.LogGammaOnePlus(1d));
            Assert.AreEqual(Gamma.LogGamma(4.5d), Gamma.LogGammaOnePlus(3.5d), 0d);
        }

        /// <summary>
        /// The inverse incomplete gamma ratio returns a finite quantile wherever one exists, including tiny shapes,
        /// probabilities of 1E-300 in either tail and very large shapes, returns 0 where the quantile is below the smallest
        /// positive double, and returns NaN for invalid arguments.
        /// </summary>
        /// <remarks>
        /// Reference roots solved with Python mpmath 1.4.1 at 60 digits. The tolerance (5E-14 + 4E-16 / a) relative allows
        /// for the condition number 1/a of the quantile at small shapes.
        /// </remarks>
        [TestMethod]
        public void Test_InverseIncompleteRatio_TailsAndEdges()
        {
            // { a, p, upper, x }
            var table = new double[,]
            {
                { 0.001, 0.5, 0, 5.2442064082779784E-302 },
                { 0.001, 0.6, 0, 7.960233816826859E-223 },
                { 0.001, 0.9, 0, 9.821659644066881E-47 },
                { 0.001, 0.1, 1, 9.82165964406658E-47 },
                { 0.001, 0.25, 1, 6.470514765088826E-126 },
                { 2.0, 1E-300, 1, 697.3242113793526 },
                { 0.5, 1E-300, 1, 686.9363156111971 },
                { 1.7782794100389228, 1E-300, 0, 2.626385352621008E-169 },
                { 1000000, 1E-08, 0, 994398.1597861212 },
                { 1000000, 1E-300, 1, 1037505.6569794014 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double a = table[i, 0], p = table[i, 1], expected = table[i, 3];
                bool upper = table[i, 2] == 1;
                double actual = upper ? Gamma.InverseUpperIncomplete(a, p) : Gamma.InverseLowerIncomplete(a, p);
                Assert.AreEqual(expected, actual, (5E-14 + 4E-16 / a) * expected, $"a = {a:R}, p = {p:R}, upper = {upper}");
            }

            // Quantiles below the double range: ln x = (ln p + ln Γ(1 + a)) / a is below -745
            Assert.AreEqual(0d, Gamma.InverseLowerIncomplete(0.01, 1E-6));
            Assert.AreEqual(0d, Gamma.InverseUpperIncomplete(0.01, 1 - 1E-6));
            Assert.AreEqual(0d, Gamma.InverseLowerIncomplete(0.5, 1E-320));
            Assert.AreEqual(0d, Gamma.InverseLowerIncomplete(2.0, 0d));
            Assert.AreEqual(double.PositiveInfinity, Gamma.InverseLowerIncomplete(2.0, 1d));
            Assert.AreEqual(double.PositiveInfinity, Gamma.InverseUpperIncomplete(2.0, 0d));
            Assert.AreEqual(0d, Gamma.InverseUpperIncomplete(2.0, 1d));
            Assert.IsTrue(double.IsNaN(Gamma.InverseLowerIncomplete(-1d, 0.5)));
            Assert.IsTrue(double.IsNaN(Gamma.InverseLowerIncomplete(2d, 1.5)));
            Assert.IsTrue(double.IsNaN(Gamma.InverseUpperIncomplete(double.NaN, 0.5)));

            // A subnormal tail probability is inverted from its logarithm
            double x = Gamma.InverseUpperIncomplete(2.0, 1E-310);
            Assert.AreEqual(Math.Log(1E-310), Gamma.LogIncompleteRatio(2.0, x, true), 1E-12 * 713.8);
        }

        /// <summary>
        /// For tiny shapes a small upper tail Q(a, x) ≈ a E1(x) puts the quantile far below one, where the bracket of the
        /// root has 0 as its lower end. The quantile is found (the bound is squared while the lower end is 0), and 0 is
        /// returned only where the quantile is below the smallest positive double.
        /// </summary>
        /// <remarks>
        /// Reference roots from Python mpmath 1.4.1 at 32 digits, by bracketed root finding on ln T(a, x) - ln p in ln x.
        /// The condition number of the quantile is Q / (x f(x)) ≈ Q / a, about 100 at these points, so the tolerance is
        /// 1E-12 relative for a = 1E-8 and 1E-11 for a = 1E-300 (measured 3.3E-12 at a = 1E-300, Q = 1E-298).
        /// </remarks>
        [TestMethod]
        public void Test_InverseIncompleteRatio_TinyShapeSmallUpperTail()
        {
            Assert.AreEqual(2.0885675164425588E-44, Gamma.InverseLowerIncomplete(1E-8, 0.999999), 1E-12 * 2.0885675164425588E-44);
            Assert.AreEqual(2.0885675224483890E-44, Gamma.InverseUpperIncomplete(1E-8, 1E-6), 1E-12 * 2.0885675224483890E-44);
            Assert.AreEqual(3.9951639287605745E-218, Gamma.InverseUpperIncomplete(1E-8, 5E-6), 1E-12 * 3.9951639287605745E-218);
            Assert.AreEqual(5.5222627516570352E-305, Gamma.InverseUpperIncomplete(1E-8, 7E-6), 1E-12 * 5.5222627516570352E-305);
            Assert.AreEqual(2.0886719363262584E-44, Gamma.InverseUpperIncomplete(1E-300, 1E-298), 1E-11 * 2.0886719363262584E-44);
            Assert.AreEqual(4.0001609899617474E-218, Gamma.InverseUpperIncomplete(1E-300, 5E-298), 1E-11 * 4.0001609899617474E-218);

            // Quantiles below the smallest positive double (about e^-1000 and smaller) return exactly 0
            var underflow = new (double a, double p, bool upper)[]
            {
                (1E-8, 1E-5, true), (1E-20, 1E-10, true), (1E-20, 1E-5, true),
                (1E-300, 1E-100, true), (1E-300, 1E-20, true), (1E-300, 1E-10, true), (1E-300, 1E-5, true),
                (1E-20, 0.999999, false), (1E-20, 0.9999999999, false), (1E-20, 0.9999999999999999, false),
                (1E-300, 0.999999, false), (1E-300, 0.9999999999, false), (1E-300, 0.9999999999999999, false),
            };
            foreach (var (a, p, upper) in underflow)
            {
                double actual = upper ? Gamma.InverseUpperIncomplete(a, p) : Gamma.InverseLowerIncomplete(a, p);
                Assert.AreEqual(0d, actual, 0d, $"a = {a:R}, p = {p:R}, upper = {upper}");
            }
        }

        /// <summary>
        /// LogGamma returns positive infinity at the poles and at positive infinity, and NaN for NaN and
        /// negative infinity.
        /// </summary>
        /// <remarks>
        /// scipy.special.gammaln also returns positive infinity at zero, at the negative integers and at
        /// positive infinity.
        /// </remarks>
        [TestMethod]
        public void Test_LogGamma_SpecialValues()
        {
            var poles = new double[] { 0d, -0d, -1d, -2d, -50d, -1E20 };
            foreach (double x in poles)
            {
                Assert.AreEqual(double.PositiveInfinity, Gamma.LogGamma(x), $"x = {x:R}");
            }

            Assert.AreEqual(double.PositiveInfinity, Gamma.LogGamma(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(Gamma.LogGamma(double.NaN)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGamma(double.NegativeInfinity)));
        }

        /// <summary>
        /// The regularized incomplete gamma integral matches 50-digit reference values on every
        /// algorithm branch, including the continued fraction taken when the argument reaches the
        /// shape parameter.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.4.1 at 50 significant digits,
        /// gammainc(alpha, 0, X, regularized=True), and cross-checked with scipy 1.17.1
        /// scipy.special.gammainc, which agrees with mpmath to the last printed digit on every
        /// case below.
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete_ReferenceValues()
        {
            // Continued-fraction branch: X > 1 and X >= alpha
            Assert.AreEqual(0.9886559833621947, Gamma.Incomplete(1.5, 0.1), 1E-11);
            Assert.AreEqual(0.9932620530009145, Gamma.Incomplete(5.0, 1.0), 1E-11);
            Assert.AreEqual(0.9796589705830716, Gamma.Incomplete(12.0, 6.0), 1E-11);
            Assert.AreEqual(0.9991559189282155, Gamma.Incomplete(45.0, 26.0), 1E-11);
            Assert.AreEqual(0.5059471461707603, Gamma.Incomplete(500.0, 500.0), 1E-11);
            Assert.AreEqual(0.9729484747641512, Gamma.Incomplete(700.0, 650.0), 1E-11);

            // Series branch: X <= 1 or X < alpha
            Assert.AreEqual(0.30146464169666126, Gamma.Incomplete(0.075, 0.5), 1E-11);
            Assert.AreEqual(0.9083579897300342, Gamma.Incomplete(0.3, 0.1), 1E-11);

            // Hill approximation branch: alpha > 10000. With the series stopping on |TERM| the
            // branch agrees with the reference to about 1E-13 here.
            Assert.AreEqual(0.8366349011570255, Gamma.Incomplete(10500.0, 10400.0), 1E-12);
        }

        /// <summary>
        /// The Hill branch (alpha &gt; 10000) sums its alternating power series to convergence when
        /// X &gt; alpha instead of stopping after the first negative term.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy 1.17 gammainc(10001, 10500) and gammainc(1e5, 101000). The pre-fix values
        /// were off by about 7E-9 absolute; the fixed branch matches to the last digit, so 1E-12
        /// leaves margin for the approximation while catching the old error.
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete_HillBranchAboveShape()
        {
            Assert.AreEqual(0.9999995501972493, Gamma.Incomplete(10500.0, 10001.0), 1E-12);
            Assert.AreEqual(0.9991915784870744, Gamma.Incomplete(101000.0, 1E5), 1E-12);
            // Below the shape parameter the terms were already positive and stay accurate
            Assert.AreEqual(0.0007574199211747684, Gamma.Incomplete(99000.0, 1E5), 1E-15);
        }

        /// <summary>
        /// NaN arguments return NaN immediately and an infinite upper limit returns one.
        /// </summary>
        [TestMethod]
        public void Test_Incomplete_NaNAndInfinity()
        {
            Assert.IsTrue(double.IsNaN(Gamma.Incomplete(double.NaN, 2.0)));
            Assert.IsTrue(double.IsNaN(Gamma.Incomplete(1.0, double.NaN)));
            Assert.AreEqual(1.0, Gamma.Incomplete(double.PositiveInfinity, 2.0));
        }

        /// <summary>
        /// Edge values of the regularized incomplete gamma functions follow scipy.special.gammainc
        /// and gammaincc, so P + Q = 1 wherever both are defined.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy 1.17 gammainc(0, 1) = 1, gammaincc(0, 1) = 0, gammainc(0, inf) = 1,
        /// gammaincc(0, inf) = 0, gammainc(1, 0) = 0, gammaincc(1, 0) = 1, and NaN for
        /// gammainc(0, 0) and a &lt; 0. These are exact values. For x &lt; 0 with a &gt; 0 the library keeps
        /// its CDF convention P = 0, Q = 1 (scipy returns NaN), which distribution code relies on below the support.
        /// </remarks>
        [TestMethod]
        public void Test_IncompleteRegularized_EdgeValues()
        {
            Assert.AreEqual(1.0, Gamma.LowerIncomplete(0.0, 1.0));
            Assert.AreEqual(0.0, Gamma.UpperIncomplete(0.0, 1.0));
            Assert.AreEqual(1.0, Gamma.LowerIncomplete(0.0, double.PositiveInfinity));
            Assert.AreEqual(0.0, Gamma.UpperIncomplete(0.0, double.PositiveInfinity));
            Assert.AreEqual(0.0, Gamma.LowerIncomplete(1.0, 0.0));
            Assert.AreEqual(1.0, Gamma.UpperIncomplete(1.0, 0.0));

            Assert.IsTrue(double.IsNaN(Gamma.LowerIncomplete(0.0, 0.0)));
            Assert.IsTrue(double.IsNaN(Gamma.UpperIncomplete(0.0, 0.0)));
            Assert.IsTrue(double.IsNaN(Gamma.LowerIncomplete(-1.0, 1.0)));
            Assert.IsTrue(double.IsNaN(Gamma.UpperIncomplete(-1.0, 1.0)));
            Assert.AreEqual(0.0, Gamma.LowerIncomplete(1.0, -1.0));
            Assert.AreEqual(1.0, Gamma.UpperIncomplete(1.0, -1.0));
        }

        /// <summary>
        /// The inverse lower incomplete gamma keeps full relative precision for small probabilities
        /// and returns the exact limits at 0 and 1.
        /// </summary>
        /// <remarks>
        /// Oracle: scipy 1.17 gammaincinv(a, y). The relative tolerance of 1E-14 is a few dozen ulps;
        /// the pre-fix values were off by 1E-4 relative at y = 1E-12 and by 22% at y = 1E-15.
        /// </remarks>
        [TestMethod]
        public void Test_InverseLowerIncomplete_SmallProbabilities()
        {
            var cases = new (double a, double y, double x)[]
            {
                (2.0, 1E-12, 1.414214229040193E-06),
                (2.0, 1E-15, 4.4721360216662474E-08),
                (0.5, 1E-10, 7.85398163397448E-21),
                (10.0, 1E-8, 0.7696086530117634),
                (2.0, 0.3, 1.0973492107034917),
                (5.0, 0.1, 2.432591025962664),
                (2.0, 0.5, 1.6783469900166612),
                (2.0, 0.9, 3.889720169867429),
            };

            foreach (var (a, y, x) in cases)
            {
                double actual = Gamma.InverseLowerIncomplete(a, y);
                Assert.AreEqual(x, actual, 1E-14 * x, $"a = {a}, y = {y}");
            }

            Assert.AreEqual(0.0, Gamma.InverseLowerIncomplete(2.0, 0.0));
            Assert.AreEqual(double.PositiveInfinity, Gamma.InverseLowerIncomplete(2.0, 1.0));
        }

        /// <summary>
        /// Test the incomplete gamma integral
        /// </summary>
        /// <remarks>
        /// Test values of x and a, and exact values of the gamma function copied from the testing of the
        /// ACM TOMS algorithm 435 of the modified incomplete Gamma function originally written in FORTRAN77
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>        
        /// Wayne Fullerton,
        /// Algorithm 435: Modified Incomplete Gamma Function,
        /// Communications of the ACM,
        /// Volume 15, Number 11, November 1972, pages 993-995
        /// </description></item>
        /// <item><description>   
        /// <see href = "https://people.math.sc.edu/Burkardt/f77_src/toms435/toms435_prb_output.txt"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete()
        {
            var testX = new double[] { 0.03d, 0.3d ,1.5d, 0.075d, 0.75d, 3.5d, 0.1d, 1d, 5d, 0.1d, 1d, 5d,
                0.15d, 1.5d, 7d, 2.5d, 12d, 16d, 25d, 45d};
            var testA = new double[] { 0.1d, 0.1d, 0.1d, 0.5d, 0.5d, 0.5d, 1d, 1d, 1d, 1.1d, 1.1d, 1.1d,
                2d, 2d, 2d, 6d, 6d, 11d, 26d};
            var testValid = new double[] { 0.73823506d, 0.90835798d, 0.98865598d, 0.30146465d, 0.77932864d,
                0.99184901d, 0.95162585E-01d, 0.63212055d, 0.99326205d, 0.72059743E-01d, 0.58918095d,0.99153680d,
                0.10185828E-01d, 0.44217461d, 0.99270493d, 0.42021040E-01d, 0.97965896d, 0.92260396d, 0.44707859d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.Incomplete(testX[i], testA[i]);
                Assert.IsLessThan(0.01d, Math.Abs(testValid[i] - testResults[i]) / testValid[i]);
            }

        }

        /// <summary>
        /// Test the upper incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_UpperIncomplete()
        {
            double expected, actual;

            // Q(0, x) = 0 for x > 0 (scipy gammaincc(0, 2) = 0), so P + Q = 1
            actual = Gamma.UpperIncomplete(0.000000, 2);
            expected = 0.000000;
            Assert.AreEqual(expected, actual);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.250000, 2);
            expected = 0.017286;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.500000, 2);
            expected = 0.045500;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(0.750000, 2);
            expected = 0.085056;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.000000, 2);
            expected = 0.135335;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.250000, 2);
            expected = 0.194847;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.500000, 2);
            expected = 0.261464;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(1.750000, 2);
            expected = 0.332706;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.000000, 2);
            expected = 0.406006;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.250000, 2);
            expected = 0.478944;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.500000, 2);
            expected = 0.549416;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.UpperIncomplete(2.750000, 2);
            expected = 0.615734;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));
        }

        /// <summary>
        /// Test the lower incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// Test values of x and a, and exact values of the lower incomplete gamma function copied from the testing of the
        /// "gammaLower_reg(  )" function from the Cephes Math Library
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Cephes Math Library Release 2.8: June, 2000
        /// </description></item>
        /// <item><description>
        /// <see href = "https://www.codecogs.com/library/maths/special/gamma/gamma_lower_reg.php"/>
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_LowerIncomplete()
        {
            double expected, actual;

            actual = Gamma.LowerIncomplete(0.000000, 2);
            expected = 1.000000;
            Assert.AreEqual(expected, actual);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.250000, 2);
            expected = 0.982714;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.500000, 2);
            expected = 0.954500;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(0.750000, 2);
            expected = 0.914944;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.000000, 2);
            expected = 0.864665;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.250000, 2);
            expected = 0.805153;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.500000, 2);
            expected = 0.738536;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(1.750000, 2);
            expected = 0.667294;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.000000, 2);
            expected = 0.593994;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.250000, 2);
            expected = 0.521056;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.500000, 2);
            expected = 0.450584;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

            actual = Gamma.LowerIncomplete(2.750000, 2);
            expected = 0.384266;
            Assert.AreEqual(expected, actual, 1e-6);
            Assert.IsFalse(double.IsNaN(actual));

        }

        /// <summary>
        /// Test the inverse of the lower regularized Gamma function
        /// </summary>
        /// <remarks>
        /// The inverse lower incomplete gamma function was test against values calculated from MATLAB's "gammaincinv( )" function
        /// <para>
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// The MathWorks, Inc. (2022). MATLAB version: 9.13.0 (R2022b).
        /// </description></item>
        /// <item><description>
        /// <see href = "https://www.mathworks.com/help/matlab/ref/gammaincinv.html" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_InverseLowerIncomplete()
        {
            var testY = new double[] { 0.1d, 0.2d, 0.3d, 0.4d, 0.5d, 0.6d, 0.7d, 0.8d, 0.9d };
            var testValid = new double[] { 0.1054d, 0.8244d, 1.9138d, 3.2113d, 4.6709d, 6.2919d, 8.1110d, 10.2325d, 12.9947d};
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Gamma.InverseLowerIncomplete(i+1, testY[i]); ;
                Assert.AreEqual(testValid[i], testResults[i], 1e-4);
            }

        }

        /// <summary>
        /// Test the inverse of the upper incomplete regularized Gamma function
        /// </summary>
        /// <remarks>
        /// Part of this code was copied from the Accord Math Library
        /// <para> 
        /// </para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_InverseUpperIncomplete()
        {
            var vec = new double[100];
            vec[0] = 0.1;
            for (int i = 1; i < 100; i++)
            {
                vec[i] = vec[i - 1] + 0.1;
            }

            foreach (double lambda in vec)
            {
                for (int i = 1; i < 100; i++)
                {
                    double x = Gamma.UpperIncomplete(lambda, i);
                    double j = Gamma.InverseUpperIncomplete(lambda, x);

                    Assert.IsLessThan(1e-2 * Math.Abs(j),Math.Abs(i - j) );
                }
            }
        }

        /// <summary>
        /// Locks the values of the upper incomplete gamma function over a grid of shapes and arguments, so callers
        /// that depend on its exact results are protected from unintended change.
        /// </summary>
        [TestMethod]
        public void Test_UpperIncomplete_RegressionGrid()
        {
            (double a, double x, double q)[] grid =
            {
                (0.5, 0.1, 0.6547208460185769),
                (0.5, 0.25, 0.47950012218695337),
                (0.5, 0.5, 0.31731050786291404),
                (0.5, 0.5 + Math.Sqrt(0.5), 0.12023834065689726),
                (0.5, 6, 0.0005320055051392495),
                (1, 0.1, 0.9048374180359595),
                (1, 0.5, 0.6065306597126334),
                (1, 1, 0.36787944117144245),
                (1, 2, 0.1353352832366127),
                (1, 7, 0.0009118819655545161),
                (2.5, 0.1, 0.9991138612111875),
                (2.5, 1.25, 0.7764950711233227),
                (2.5, 2.5, 0.41588018699550783),
                (2.5, 2.5 + Math.Sqrt(2.5), 0.14751604793196355),
                (2.5, 10, 0.0012497305630313758),
                (10, 0.1, 1),
                (10, 5, 0.9681719426937951),
                (10, 10, 0.4579297144718521),
                (10, 10 + Math.Sqrt(10), 0.15535840360363729),
                (10, 25, 0.00022147663824878351),
                (100, 50, 0.9999999996799934),
                (100, 100, 0.48670120172085135),
                (100, 110, 0.1582786700600871),
                (100, 205, 1.3996754787754716E-16),
                (1000, 1000, 0.4957947558197845),
                (1000, 1000 + Math.Sqrt(1000), 0.15861581157127916),
                (1000, 2005, 5.575344815222742E-137),
                (10000, 10000, 0.4986701916600448),
                (10000, 10100, 0.15865124955282034),
                (100000, 100000, 0.4995794778896348),
                (100000, 100000 + Math.Sqrt(100000), 0.1586548515516751),
                (1000000, 1000000, 0.4998670192391274),
                (1000000, 1001000, 0.15865521363165966),
            };
            foreach (var (a, x, q) in grid)
                Assert.AreEqual(q, Gamma.UpperIncomplete(a, x), 0d);
        }

        /// <summary>
        /// LogGammaDifference returns ln Γ(x + c) - ln Γ(x) to full relative accuracy for small increments of either
        /// sign and large arguments, where the difference of the two log-gamma values cancels.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath 1.3.0 at 60 significant digits as loggamma(x + c) - loggamma(x).
        /// The tolerance is 2E-15 relative.
        /// </remarks>
        [TestMethod]
        public void Test_LogGammaDifference_ReferenceValues()
        {
            double[,] table =
            {
                { 1000.0, 0.001, 0.0069072556958987286 },
                { 1000.0, -0.001, -0.006907254695398562 },
                { 2.5, 0.001, 0.0007034017801652782 },
                { 0.5, 1e-06, -1.963507558623128e-06 },
                { 0.5, -0.4999, 8.637917715709373 },
                { 4001.0, 0.3, 2.4882636387808064 },
                { 1000000.0, -0.5, -6.907754903982012 },
                { 3.0, 10.0, 19.29406731510194 },
                { 20.0, -19.5, -38.767519244274794 },
                { 8.0, 1e-12, 2.0156414779556764e-12 },
                { 1.0005, -0.001, 0.0005772157650729545 },
                { 7.5, 0.75, 1.49882268284639 }
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double x = table[i, 0], c = table[i, 1], expected = table[i, 2];
                Assert.AreEqual(expected, Gamma.LogGammaDifference(x, c), 2E-15 * Math.Abs(expected), $"x = {x:R}, c = {c:R}");
            }
            Assert.AreEqual(0d, Gamma.LogGammaDifference(5d, 0d), 0d);
        }

        /// <summary>
        /// LogGammaDifference returns NaN outside its domain, where x or x + c is not a positive finite value,
        /// instead of looping or taking the logarithm of a negative number.
        /// </summary>
        [TestMethod]
        public void Test_LogGammaDifference_OutsideDomainIsNaN()
        {
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(1d, -1d)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(1d, -2d)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(1d, -1E20)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(1d, double.NegativeInfinity)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(0d, 1d)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(-3d, 0d)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(double.NaN, 1d)));
            Assert.IsTrue(double.IsNaN(Gamma.LogGammaDifference(double.PositiveInfinity, 1d)));
        }
    }
}
