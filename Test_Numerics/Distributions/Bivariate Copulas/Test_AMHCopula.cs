using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{

    /// <summary>
    /// Unit tests for the AMH Copula. Reference values come from the R 'copula' package for the original methods and from pyvinecopulib 0.7.6 and mpmath 1.4.1 for the conditional-distribution methods.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_AMHCopula
    {

        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new AMHCopula(0.8);
            Assert.AreEqual(0.6491167, copula.PDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(0.9630927, copula.PDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(1.158995, copula.PDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(1.259423, copula.PDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new AMHCopula(0.8);
            Assert.AreEqual(0.6491167, copula.PDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(0.9630927, copula.PDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(1.158995, copula.PDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(1.259423, copula.PDF(0.8, 0.2), 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new AMHCopula(0.8);
            Assert.AreEqual(0.1834862, copula.CDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(0.1626016, copula.CDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(0.1481481, copula.CDF(0.2, 0.8), 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(0.1398601, copula.CDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new AMHCopula(0.8);
            Assert.AreEqual(0.1834862, copula.CDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(0.1);
            Assert.AreEqual(0.1626016, copula.CDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(-0.5);
            Assert.AreEqual(0.1481481, copula.CDF(0.8, 0.2), 1E-6);

            copula = new AMHCopula(-0.9);
            Assert.AreEqual(0.1398601, copula.CDF(0.8, 0.2), 1E-6);

        }

        private double[] data1 = new double[] { 82.6207861517446, 122.99162719591, 85.5464857510852, 88.271076232573, 75.5957634925482, 115.666296571075, 88.1535968625871, 89.0486583921505, 95.6520223040253, 93.8704340350363, 94.8983409930905, 98.190034452137, 119.062097667798, 101.396326234446, 102.044102658515, 122.09278272123, 87.5930294017727, 89.467574539934, 68.759249687252, 97.3170891237622, 115.23726463545, 88.7249972292411, 104.808913056224, 81.4376985296622, 88.5496448465636, 85.9566777828142, 124.159708173032, 82.3989677011281, 81.0754908215474, 96.9227210098483, 81.9598285215311, 106.023194937944, 98.9429238885091, 99.6105440441767, 76.8967482790596, 93.8733402407522, 97.1593258552729, 83.8476766091024, 83.3372066079669, 84.2378120351813, 130.294157330837, 92.0183206123361, 91.0344208437507, 62.6399930777225, 100.724150023017, 105.381851336153, 106.052861517738, 105.695937045651, 98.7180560328115, 103.659278801543, 135.355695309546, 112.126374548073, 75.8894228126445, 113.303033282488, 108.572361018471, 84.4697433912218, 102.158309989318, 104.169634603778, 126.258771553477, 92.0530458784438, 121.851579456857, 112.009972090868, 82.566074288035, 100.197743002857, 109.653015549203, 103.193458190759, 120.818259910061, 109.584519216444, 102.469182012814, 107.385973965207, 97.2671206491186, 112.417149113877, 97.8701780074462, 86.9632866191297, 129.722705545036, 118.480752744013, 109.568435955988, 122.488504266404, 112.534825234329, 110.178551830369, 85.8393034892458, 128.148787670707, 88.8620889850668, 100.212554271105, 111.449083807286, 97.9922242053636, 96.0208140774823, 118.595431696388, 96.5101956950509, 89.0070522979225, 111.690497497928, 98.2256095188502, 109.751617939794, 113.262405920361, 89.1969936327915, 126.836051984505, 94.7671769565446, 101.621377036869, 85.9837125336272, 80.6830712615443 };
        private double[] data2 = new double[] { 23.3689941895676, 139.075605092248, 125.391972129225, 68.0466559628945, 21.3849501824245, 72.623739980515, 118.431110203811, 55.6433779771001, 61.7394973031962, 117.820425024964, 51.4915091268943, 90.6318678923286, 90.2225515717848, 82.9583115041602, 58.8208994489536, 89.3934020879786, 84.2687389199291, 52.8123343271364, 133.922245003475, 65.7460273014489, 94.0306797439919, 92.7622201239545, 67.4956202291755, 52.1748570092584, 76.289463467535, 30.6442121096422, 80.367744128015, 84.589322268581, 89.4931756387178, 99.8444599387594, 147.360635147879, 26.7990606303491, 75.9359921708942, 132.224006705134, 43.5713891672347, 69.4331898107351, 113.826594948864, 44.856181127981, 59.6818830690231, 34.8777185962398, 124.522970067359, 55.0531003117469, 112.641824278124, 53.3587463398865, 82.3818569319419, 131.306495042214, 92.5500215682676, 96.1266469742884, 126.22260773411, 105.118018987749, 98.5541508572697, 81.6528854407349, 77.1433525995117, 76.8338526759968, 89.723133164602, 55.8177007161922, 94.0035174340863, 164.923630526842, 133.423371650793, 89.0043059347403, 89.7587258229626, 88.4640664624267, 104.247452450998, 67.9322685363551, 83.4200659472935, 73.8431913350094, 17.7019731704449, 114.710483451006, 101.065426716036, 55.4556570293316, 67.8514034006699, 127.73989553804, 125.34525003517, 57.0971025319617, 85.1768621826704, 101.705832867004, 85.3021580666208, 89.4450158150175, 75.5941517929671, 124.424804607655, 83.3449371155293, 96.7642165625312, 57.9643876742987, 126.024606562933, 66.5425894700683, 125.964000028259, 45.1213432721552, 86.2256615507381, 101.231675185045, 27.2201147115773, 95.015298927854, 63.4674026223287, 86.224525433146, 116.818569536018, 63.082998559889, 73.7218123340297, 92.935473959462, 84.2370209154411, 49.9905902403315, 61.0021219814154 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new AMHCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(0.7965838, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Test fitting with the method of maximum pseudo likelihood.
        /// </summary>
        [TestMethod]
        public void Test_MPL_Fit()
        {
            // get ranks of data
            var rank1 = Statistics.RanksInPlace(data1);
            var rank2 = Statistics.RanksInPlace(data2);
            // get plotting positions
            for (int i = 0; i < data1.Length; i++)
            {
                rank1[i] = rank1[i] / (rank1.Length + 1d);
                rank2[i] = rank2[i] / (rank2.Length + 1d);
            }
            // Fit copula
            BivariateCopula copula = new AMHCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(0.8321504, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new AMHCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(0.8392506, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new AMHCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(0.8367644, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(100.2479, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(14.74416, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(83.76367, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(29.76202, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }

        /// <summary>
        /// Test the tail dependence coefficients.
        /// AMH copula: λ_U = λ_L = 0 (no tail dependence).
        /// </summary>
        [TestMethod]
        public void Test_TailDependence()
        {
            var copula = new AMHCopula(0.5);
            Assert.AreEqual(0.0, copula.UpperTailDependence, 1E-10, "AMH copula should have no upper tail dependence.");
            Assert.AreEqual(0.0, copula.LowerTailDependence, 1E-10, "AMH copula should have no lower tail dependence.");

            var copulaNeg = new AMHCopula(-0.5);
            Assert.AreEqual(0.0, copulaNeg.UpperTailDependence, 1E-10, "AMH copula should have no upper tail dependence for negative θ.");
            Assert.AreEqual(0.0, copulaNeg.LowerTailDependence, 1E-10, "AMH copula should have no lower tail dependence for negative θ.");
        }

        /// <summary>
        /// Test that ParametersValid tracks the dependency parameter's valid range.
        /// </summary>
        [TestMethod]
        public void Test_ParametersValid()
        {
            var copula = new AMHCopula(0.4);
            Assert.IsTrue(copula.ParametersValid);

            copula.Theta = 1.5; // above the maximum of 1
            Assert.IsFalse(copula.ParametersValid);

            copula.Theta = -0.5;
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test Clone produces an independent copy with deep-copied marginals.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new AMHCopula(0.4, new Normal(100, 10), new Gumbel(50, 5));
            var clone = copula.Clone() as AMHCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(copula.Theta, clone.Theta);
            Assert.IsTrue(clone.ParametersValid);

            // Marginals are deep-copied (distributions memoize lazily, so clones must
            // not share marginal instances), with identical quantiles
            Assert.AreNotSame(copula.MarginalDistributionX, clone.MarginalDistributionX);
            Assert.AreNotSame(copula.MarginalDistributionY, clone.MarginalDistributionY);
            Assert.AreEqual(copula.MarginalDistributionY.InverseCDF(0.9), clone.MarginalDistributionY.InverseCDF(0.9));

            // Mutating the clone should not affect the original
            clone.Theta = 0.8;
            Assert.AreEqual(0.4, copula.Theta);
        }

        /// <summary>
        /// Test the forward conditional CDF (h-function), h(v|u) = ∂C(u,v)/∂u, covering both
        /// positive and negative dependence. The Ali-Mikhail-Haq family is not available in
        /// pyvinecopulib, so reference values were computed with the Python package mpmath
        /// 1.4.1 as a 50-digit numerical ∂C/∂u of the AMH CDF C(u,v) = uv/(1 − θ(1−u)(1−v)),
        /// cross-checked at generation against the closed form v(1 − θ(1−v))/D²,
        /// D = 1 − θ(1−u)(1−v), to 1E-30. The analytic conditional is also cross-checked
        /// against a central finite difference of the closed-form CDF (ε = 1E-6; noise floor
        /// ~1E-9, so 1E-7 is asserted) and pinned at the exact upper edge h(1|u) = 1.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDF()
        {
            var copula = new AMHCopula(0.5);
            Assert.AreEqual(0.742798289691333, copula.ConditionalCDF(0.3, 0.7), 1E-12);
            Assert.AreEqual(0.24343809494085702, copula.ConditionalCDF(0.7, 0.3), 1E-12);
            Assert.AreEqual(0.9424018848037696, copula.ConditionalCDF(0.05, 0.9), 1E-12);
            Assert.AreEqual(0.028933391200115736, copula.ConditionalCDF(0.9, 0.05), 1E-12);

            copula = new AMHCopula(-0.7);
            Assert.AreEqual(0.6438083047470791, copula.ConditionalCDF(0.3, 0.7), 1E-12);
            Assert.AreEqual(0.3397666023871834, copula.ConditionalCDF(0.7, 0.3), 1E-12);
            Assert.AreEqual(0.8466512766037415, copula.ConditionalCDF(0.05, 0.9), 1E-12);
            Assert.AreEqual(0.07319181596808046, copula.ConditionalCDF(0.9, 0.05), 1E-12);

            foreach (double theta in new[] { 0.5, -0.7 })
            {
                var c = new AMHCopula(theta);
                foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    foreach (double v in new[] { 0.1, 0.5, 0.9 })
                    {
                        double fd = (c.CDF(u + 1E-6, v) - c.CDF(u - 1E-6, v)) / 2E-6;
                        Assert.AreEqual(fd, c.ConditionalCDF(u, v), 1E-7);
                    }
                    Assert.AreEqual(1d, c.ConditionalCDF(u, 1d), 1E-12);
                }
            }
        }

        /// <summary>
        /// Test the scalar inverse conditional CDF, covering both positive and negative
        /// dependence. Reference values were computed with the Python package mpmath 1.4.1
        /// as a 50-digit root solve of the conditional CDF (the AMH family is not available
        /// in pyvinecopulib). The closed-form round trip
        /// h(InverseConditionalCDF(u, t) | u) = t is asserted at 1E-10, and
        /// InverseCDF(u, t)[1] must recompose the scalar bit-for-bit.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDF()
        {
            var copula = new AMHCopula(0.5);
            Assert.AreEqual(0.6538170057580951, copula.InverseConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.08399461695722253, copula.InverseConditionalCDF(0.9, 0.05), 1E-10);

            copula = new AMHCopula(-0.7);
            Assert.AreEqual(0.750531409535383, copula.InverseConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.033999893497902646, copula.InverseConditionalCDF(0.9, 0.05), 1E-10);

            foreach (double theta in new[] { 0.5, -0.7 })
            {
                var c = new AMHCopula(theta);
                foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    foreach (double t in new[] { 0.1, 0.5, 0.9 })
                    {
                        double v = c.InverseConditionalCDF(u, t);
                        Assert.AreEqual(t, c.ConditionalCDF(u, v), 1E-10);
                        Assert.AreEqual(v, c.InverseCDF(u, t)[1], 0d);
                    }
                }
            }
        }

        /// <summary>
        /// Every generator member against mpmath 1.3.0 at 150 digits, for φ(t) = ln((1 − θ(1 − t))/t), the independence generator −ln t at θ = 0: Generator, GeneratorInverse (at
        /// the rounded generator value), GeneratorPrime, GeneratorPrime2 and GeneratorPrimeInverse (at the rounded slope),
        /// on t from 1E-12 to 1 - 1E-12 and θ including the independence region, to 1E-13 relative.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorMembersMatchMpmath()
        {
            // theta, t, phi(t), phi^-1(phi), phi'(t), phi''(t), (phi')^-1(phi').
            double[,] cases =
            {
            { -0.99, 1e-12, 28.31915575466445, 1.0000000000000014e-12, -1000000000000.4976, 1e+24, 1e-12 },
            { -0.99, 0.01, 5.288317534469703, 0.009999999999999998, -100.49997474875006, 9999.750025250612, 0.01 },
            { -0.99, 0.5, 1.095273387402595, 0.5, -2.6622073578595318, 3.561481415196698, 0.5 },
            { -0.99, 0.99, 0.01990165190387564, 0.99, -1.9903960888216754, 0.05932560924312536, 0.9900000000000018 },
            { -0.99, 0.999999999999, 1.9899559777769683e-12, 0.999999999999, -1.99000000000002, 0.01990000000394053, 0.9999999999989958 },
            { -1e-12, 1e-12, 27.63102111592955, 9.999999999999992e-13, -1000000000000.0, 1e+24, 1e-12 },
            { -1e-12, 0.01, 4.605170185989081, 0.010000000000000002, -100.000000000001, 10000.0, 0.01 },
            { -1e-12, 0.5, 0.6931471805604453, 0.5, -2.000000000001, 4.0, 0.5 },
            { -1e-12, 0.99, 0.01005033585351145, 0.99, -1.01010101010201, 1.0203040506070808, 0.9900000000000001 },
            { -1e-12, 0.999999999999, 9.999778782813784e-13, 0.999999999999, -1.000000000002, 1.000000000002, 0.999999999999 },
            { 0.0, 1e-12, 27.631021115928547, 1.000000000000001e-12, -1000000000000.0, 1e+24, 1e-12 },
            { 0.0, 0.01, 4.605170185988091, 0.010000000000000004, -100.0, 10000.0, 0.01 },
            { 0.0, 0.5, 0.6931471805599453, 0.5, -2.0, 4.0, 0.5 },
            { 0.0, 0.99, 0.01005033585350145, 0.99, -1.0101010101010102, 1.0203040506070808, 0.99 },
            { 0.0, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.000000000001, 1.000000000002, 0.9999999999989999 },
            { 0.5, 1e-12, 26.937873935369602, 1.000000000000001e-12, -999999999999.0, 1e+24, 1e-12 },
            { 0.5, 0.01, 3.9219733362813143, 0.009999999999999998, -99.00990099009901, 9999.019703950593, 0.01 },
            { 0.5, 0.5, 0.4054651081081644, 0.5, -1.3333333333333333, 3.5555555555555554, 0.5 },
            { 0.5, 0.99, 0.0050377940299571634, 0.99, -0.5075884472869397, 0.767785174821116, 0.99 },
            { 0.5, 0.999999999999, 4.999889391403143e-13, 0.999999999999, -0.50000000000075, 0.7500000000017499, 0.999999999999 },
            { 0.99, 1e-12, 23.025850930039457, 1.0000000000000006e-12, -999999999901.0, 1e+24, 1e-12 },
            { 0.99, 0.01, 0.6881346387364015, 0.01, -50.25125628140706, 7525.062498421759, 0.01 },
            { 0.99, 0.5, 0.009950330853168092, 0.5, -0.03960396039603964, 0.15684736790510748, 0.5 },
            { 0.99, 0.99, 0.00010100499983335884, 0.99, -0.010202010000010212, 0.020506040204101534, 0.99 },
            { 0.99, 0.999999999999, 9.999778782808743e-15, 0.999999999999, -0.010000000000019908, 0.01990000000005942, 0.999999999999 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new AMHCopula(cases[i, 0]);
                double t = cases[i, 1];
                AssertClose(cases[i, 2], copula.Generator(t), $"case {i}: Generator");
                if (cases[i, 2] > 0d && !double.IsInfinity(cases[i, 2]))
                    AssertClose(cases[i, 3], copula.GeneratorInverse(cases[i, 2]), $"case {i}: GeneratorInverse");
                AssertClose(cases[i, 4], copula.GeneratorPrime(t), $"case {i}: GeneratorPrime");
                AssertClose(cases[i, 5], copula.GeneratorPrime2(t), $"case {i}: GeneratorPrime2");
                if (cases[i, 4] < 0d && !double.IsInfinity(cases[i, 4]))
                    AssertClose(cases[i, 6], copula.GeneratorPrimeInverse(cases[i, 4]), $"case {i}: GeneratorPrimeInverse");
            }
        }

        private static void AssertClose(double expected, double actual, string message)
        {
            if (double.IsInfinity(expected) || expected == 0d)
                Assert.AreEqual(expected, actual, 1E-300, message);
            else
                Assert.AreEqual(expected, actual, 1E-13 * Math.Abs(expected), message);
        }

        /// <summary>
        /// InverseConditionalCDF against the root of h(v|u) = t found by mpmath 1.3.0 at 150 digits, including θ = 0
        /// (independence, where the quadratic degenerates), u = 0 and 1, and t = 1E-12, to 1E-13 relative.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDFMatchesMpmath()
        {
            double[,] cases =
            {
            { -0.99, 0.0, 1e-12, 1.98999999999803e-12 },
            { -0.99, 0.0, 0.5, 0.6655518394648829 },
            { -0.99, 0.0, 0.9, 0.9471179270227393 },
            { -0.99, 0.3, 1e-12, 1.4403261306526e-12 },
            { -0.99, 0.3, 0.5, 0.5870931406611762 },
            { -0.99, 0.3, 0.9, 0.9265696347374848 },
            { -0.99, 1.0, 1e-12, 5.02512562814196e-13 },
            { -0.99, 1.0, 0.5, 0.29436350457014576 },
            { -0.99, 1.0, 0.9, 0.6871895156048573 },
            { 0.0, 0.0, 1e-12, 1e-12 },
            { 0.0, 0.0, 0.5, 0.5 },
            { 0.0, 0.0, 0.9, 0.9 },
            { 0.0, 0.3, 1e-12, 1e-12 },
            { 0.0, 0.3, 0.5, 0.5 },
            { 0.0, 0.3, 0.9, 0.9 },
            { 0.0, 1.0, 1e-12, 1e-12 },
            { 0.0, 1.0, 0.5, 0.5 },
            { 0.0, 1.0, 0.9, 0.9 },
            { 0.5, 0.0, 1e-12, 5.0000000000025e-13 },
            { 0.5, 0.0, 0.5, 0.3333333333333333 },
            { 0.5, 0.0, 0.9, 0.8181818181818182 },
            { 0.5, 0.3, 1e-12, 8.450000000000549e-13 },
            { 0.5, 0.3, 0.5, 0.4496674680413865 },
            { 0.5, 0.3, 0.9, 0.8784373287140119 },
            { 0.5, 1.0, 1e-12, 1.9999999999959998e-12 },
            { 0.5, 1.0, 0.5, 0.6180339887498949 },
            { 0.5, 1.0, 0.9, 0.9317821063276354 },
            { 0.99, 0.0, 1e-12, 1.0000000000009908e-14 },
            { 0.99, 0.0, 0.5, 0.00990099009900991 },
            { 0.99, 0.0, 0.9, 0.08256880733944963 },
            { 0.99, 0.3, 1e-12, 9.424899991606978e-12 },
            { 0.99, 0.3, 0.5, 0.42000440550164764 },
            { 0.99, 0.3, 0.9, 0.8479796616237846 },
            { 0.99, 1.0, 1e-12, 9.999999900999993e-11 },
            { 0.99, 1.0, 0.5, 0.7056364954298543 },
            { 0.99, 1.0, 0.9, 0.9484254604002106 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new AMHCopula(cases[i, 0]);
                Assert.AreEqual(cases[i, 3], copula.InverseConditionalCDF(cases[i, 1], cases[i, 2]), 1E-13 * cases[i, 3], $"case {i}");
            }
        }

        /// <summary>
        /// t = 0 returns 0 at every θ, including θ = 1, where the quadratic degenerates; t = 1 returns 1.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDFAtZeroAndOne()
        {
            foreach (double theta in new[] { -1d, -0.5, 0d, 0.5, 1d })
            {
                var copula = new AMHCopula(theta);
                foreach (double u in new[] { 0d, 0.3, 1d })
                {
                    Assert.AreEqual(0d, copula.InverseConditionalCDF(u, 0d), 0d, $"theta={theta}, u={u}");
                    Assert.AreEqual(1d, copula.InverseConditionalCDF(u, 1d), 1E-15, $"theta={theta}, u={u}");
                }
            }
        }

        /// <summary>
        /// At θ = 1, where the generator is identically zero, the CDF and the conditional CDF use the closed forms
        /// C = uv / (u + v - uv) and h = v² / (u + v - uv)², and agree with θ = 1 - 1E-9.
        /// </summary>
        [TestMethod]
        public void Test_ThetaOne_ClosedForm()
        {
            var copula = new AMHCopula(1d);
            Assert.AreEqual(0.25d, copula.CDF(0.3d, 0.6d), 1E-15);
            Assert.AreEqual(0.36d / 0.5184d, copula.ConditionalCDF(0.3d, 0.6d), 1E-15);
            Assert.AreEqual(0d, copula.CDF(0d, 0.5d));
            Assert.AreEqual(0d, copula.ConditionalCDF(0.5d, 0d));
            Assert.AreEqual(1d, copula.ConditionalCDF(0.5d, 1d));
            var near = new AMHCopula(1d - 1E-9);
            foreach (double u in new[] { 0.01d, 0.2d, 0.5d, 0.8d, 0.99d })
            {
                foreach (double v in new[] { 0.01d, 0.2d, 0.5d, 0.8d, 0.99d })
                {
                    Assert.AreEqual(near.CDF(u, v), copula.CDF(u, v), 1E-7, $"u={u}, v={v}");
                    Assert.AreEqual(near.ConditionalCDF(u, v), copula.ConditionalCDF(u, v), 1E-7, $"u={u}, v={v}");
                }
            }
        }
    }
}
