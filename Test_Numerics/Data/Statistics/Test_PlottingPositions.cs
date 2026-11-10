using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;

namespace Data.Statistics
{
    /// <summary>
    /// Unit tests for the PlottingPositions class. These methods were tested against values calculated directly from their formulas.
    /// </summary>
    /// <remarks>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     <item>Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil</item>
    ///     </list>
    /// </remarks>
    [TestClass]
    public class Test_PlottingPositions
    {
        /// <summary>
        /// Test the basic plotting position function.
        /// </summary>
        [TestMethod]
        public void Test_Function()
        {
            int n = 30;
            double alpha = 0.1;
            var test = PlottingPositions.Function(n, alpha);
            var valid = new double[] { 0.0292207792207792, 0.0616883116883117, 0.0941558441558441, 0.126623376623377, 0.159090909090909, 0.191558441558442, 0.224025974025974, 0.256493506493506, 0.288961038961039, 0.321428571428571, 0.353896103896104, 0.386363636363636, 0.418831168831169, 0.451298701298701, 0.483766233766234, 0.516233766233766, 0.548701298701299, 0.581168831168831, 0.613636363636364, 0.646103896103896, 0.678571428571428, 0.711038961038961, 0.743506493506493, 0.775974025974026, 0.808441558441558, 0.840909090909091, 0.873376623376623, 0.905844155844156, 0.938311688311688, 0.970779220779221 };

            for(int i = 0; i < test.Length; i++)
            {
                Assert.AreEqual(valid[i], test[i], 1E-6);
            }
        }

        /// <summary>
        /// One and two values take the formula (i − α) / (N + 1 − 2α) directly: Weibull 1/2 for one value and
        /// 1/3, 2/3 for two; Hazen 1/2 and 1/4, 3/4; Kaplan-Meier 1 and 1/2, 1. Invalid sizes, a NaN alpha and the
        /// 0/0 case N = 1 with α = 1 are refused.
        /// </summary>
        [TestMethod]
        public void Test_SmallSamples()
        {
            CollectionAssert.AreEqual(new[] { 0.5 }, PlottingPositions.Weibull(1));
            CollectionAssert.AreEqual(new[] { 1d / 3d, 2d / 3d }, PlottingPositions.Weibull(2));
            CollectionAssert.AreEqual(new[] { 0.5 }, PlottingPositions.Hazen(1));
            CollectionAssert.AreEqual(new[] { 0.25, 0.75 }, PlottingPositions.Hazen(2));
            CollectionAssert.AreEqual(new[] { 1d }, PlottingPositions.KaplanMeier(1));
            CollectionAssert.AreEqual(new[] { 0.5, 1d }, PlottingPositions.KaplanMeier(2));
            CollectionAssert.AreEqual(new[] { 0.5 }, PlottingPositions.Function(1, PlottingPositions.PlottingPostionType.Gringorten));
            CollectionAssert.AreEqual(new[] { 0d, 1d }, PlottingPositions.Function(2, 1d));

            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(0, 0.4));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.KaplanMeier(0));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(5, double.NaN));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(1, 1d));
        }

        /// <summary>
        /// Test the Blom plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Blom()
        {
            int n = 30;
            var test1 = PlottingPositions.Blom(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Blom);
            var valid = new double[] { 0.0206611570247934, 0.0537190082644628, 0.0867768595041322, 0.119834710743802, 0.152892561983471, 0.18595041322314, 0.21900826446281, 0.252066115702479, 0.285123966942149, 0.318181818181818, 0.351239669421488, 0.384297520661157, 0.417355371900826, 0.450413223140496, 0.483471074380165, 0.516528925619835, 0.549586776859504, 0.582644628099174, 0.615702479338843, 0.648760330578512, 0.681818181818182, 0.714876033057851, 0.747933884297521, 0.78099173553719, 0.814049586776859, 0.847107438016529, 0.880165289256198, 0.913223140495868, 0.946280991735537, 0.979338842975207 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the Cunnane plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Cunnane()
        {
            int n = 30;
            var test1 = PlottingPositions.Cunnane(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Cunnane);
            var valid = new double[] { 0.0198675496688742, 0.0529801324503311, 0.0860927152317881, 0.119205298013245, 0.152317880794702, 0.185430463576159, 0.218543046357616, 0.251655629139073, 0.28476821192053, 0.317880794701987, 0.350993377483444, 0.384105960264901, 0.417218543046358, 0.450331125827815, 0.483443708609272, 0.516556291390728, 0.549668874172185, 0.582781456953642, 0.615894039735099, 0.649006622516556, 0.682119205298013, 0.71523178807947, 0.748344370860927, 0.781456953642384, 0.814569536423841, 0.847682119205298, 0.880794701986755, 0.913907284768212, 0.947019867549669, 0.980132450331126 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the Gringorten plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Gringorten()
        {
            int n = 30;
            var test1 = PlottingPositions.Gringorten(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Gringorten);
            var valid = new double[] { 0.0185922974767596, 0.051792828685259, 0.0849933598937583, 0.118193891102258, 0.151394422310757, 0.184594953519256, 0.217795484727756, 0.250996015936255, 0.284196547144754, 0.317397078353254, 0.350597609561753, 0.383798140770252, 0.416998671978752, 0.450199203187251, 0.48339973439575, 0.51660026560425, 0.549800796812749, 0.583001328021248, 0.616201859229748, 0.649402390438247, 0.682602921646746, 0.715803452855246, 0.749003984063745, 0.782204515272244, 0.815405046480744, 0.848605577689243, 0.881806108897742, 0.915006640106242, 0.948207171314741, 0.98140770252324 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the Hazen plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Hazen()
        {
            int n = 30;
            var test1 = PlottingPositions.Hazen(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Hazen);
            var valid = new double[] { 0.0166666666666667, 0.05, 0.0833333333333333, 0.116666666666667, 0.15, 0.183333333333333, 0.216666666666667, 0.25, 0.283333333333333, 0.316666666666667, 0.35, 0.383333333333333, 0.416666666666667, 0.45, 0.483333333333333, 0.516666666666667, 0.55, 0.583333333333333, 0.616666666666667, 0.65, 0.683333333333333, 0.716666666666667, 0.75, 0.783333333333333, 0.816666666666667, 0.85, 0.883333333333333, 0.916666666666667, 0.95, 0.983333333333333 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the Kaplan-Meier plotting position.
        /// </summary>
        [TestMethod]
        public void Test_KaplanMeier()
        {
            int n = 30;
            var test1 = PlottingPositions.KaplanMeier(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.KaplanMeier);
            var valid = new double[] { 0.0333333333333333, 0.0666666666666667, 0.1, 0.133333333333333, 0.166666666666667, 0.2, 0.233333333333333, 0.266666666666667, 0.3, 0.333333333333333, 0.366666666666667, 0.4, 0.433333333333333, 0.466666666666667, 0.5, 0.533333333333333, 0.566666666666667, 0.6, 0.633333333333333, 0.666666666666667, 0.7, 0.733333333333333, 0.766666666666667, 0.8, 0.833333333333333, 0.866666666666667, 0.9, 0.933333333333333, 0.966666666666667, 1 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the Median plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Median()
        {
            int n = 30;
            var test1 = PlottingPositions.Median(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Median);
            var valid = new double[] { 0.0224765354849333, 0.0554091882101103, 0.0883418409352874, 0.121274493660464, 0.154207146385641, 0.187139799110818, 0.220072451835995, 0.253005104561172, 0.285937757286349, 0.318870410011526, 0.351803062736703, 0.38473571546188, 0.417668368187057, 0.450601020912234, 0.483533673637412, 0.516466326362588, 0.549398979087766, 0.582331631812943, 0.61526428453812, 0.648196937263297, 0.681129589988474, 0.714062242713651, 0.746994895438828, 0.779927548164005, 0.812860200889182, 0.845792853614359, 0.878725506339536, 0.911658159064713, 0.94459081178989, 0.977523464515067 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }       
        }

        /// <summary>
        /// Test the Weibull plotting position.
        /// </summary>
        [TestMethod]
        public void Test_Weibull()
        {
            int n = 30;
            var test1 = PlottingPositions.Weibull(n);
            var test2 = PlottingPositions.Function(n, PlottingPositions.PlottingPostionType.Weibull);
            var valid = new double[] { 0.032258064516129, 0.0645161290322581, 0.0967741935483871, 0.129032258064516, 0.161290322580645, 0.193548387096774, 0.225806451612903, 0.258064516129032, 0.290322580645161, 0.32258064516129, 0.354838709677419, 0.387096774193548, 0.419354838709677, 0.451612903225806, 0.483870967741935, 0.516129032258065, 0.548387096774194, 0.580645161290323, 0.612903225806452, 0.645161290322581, 0.67741935483871, 0.709677419354839, 0.741935483870968, 0.774193548387097, 0.806451612903226, 0.838709677419355, 0.870967741935484, 0.903225806451613, 0.935483870967742, 0.967741935483871 };

            for (int i = 0; i < test1.Length; i++)
            {
                Assert.AreEqual(valid[i], test1[i], 1E-6);
                Assert.AreEqual(valid[i], test2[i], 1E-6);
            }
        }

        /// <summary>
        /// Test the weighted plotting positions: unit weights reproduce the unweighted formula bit for bit,
        /// integer weights and scaled weights follow the formula, and invalid weights are refused.
        /// </summary>
        [TestMethod]
        public void Test_WeightedFunction()
        {
            foreach (double alpha in new[] { 0.0, 0.3175, 0.375, 0.4, 0.44, 0.5 })
            {
                foreach (int n in new[] { 3, 10, 30, 1001 })
                {
                    var ones = new double[n];
                    for (int i = 0; i < n; i++) ones[i] = 1.0;
                    var weighted = PlottingPositions.Function(ones, alpha);
                    var plain = PlottingPositions.Function(n, alpha);
                    for (int i = 0; i < n; i++)
                        Assert.AreEqual(System.BitConverter.DoubleToInt64Bits(plain[i]), System.BitConverter.DoubleToInt64Bits(weighted[i]), $"alpha {alpha}, n {n}, i {i}");
                }
            }

            // Weights {1, 2, 1}, alpha = 0.5: W = 4, mean 4/3, denominator 4 + 4/3 - 4/3 = 4;
            // F = (1 - 0.5) / 4, (3 - 1) / 4, (4 - 0.5) / 4.
            var hazen = PlottingPositions.Function(new[] { 1.0, 2.0, 1.0 }, 0.5);
            CollectionAssert.AreEqual(new[] { 0.125, 0.5, 0.875 }, hazen);

            // Weibull with weights {1, 0, 3}: W = 4, mean 4/3, denominator 16/3; F = 3/16, 3/16, 12/16.
            var weibull = PlottingPositions.Function(new[] { 1.0, 0.0, 3.0 }, 0.0);
            Assert.AreEqual(3.0 / 16.0, weibull[0], 1E-15);
            Assert.AreEqual(weibull[0], weibull[1]);
            Assert.AreEqual(0.75, weibull[2], 1E-15);

            // Positions depend only on the relative weights.
            var scaled = PlottingPositions.Function(new[] { 0.25, 0.5, 0.25 }, 0.5);
            CollectionAssert.AreEqual(hazen, scaled);

            Assert.ThrowsExactly<System.ArgumentNullException>(() => PlottingPositions.Function(null!, 0.5));
            Assert.ThrowsExactly<System.ArgumentException>(() => PlottingPositions.Function(new double[0], 0.5));
            Assert.ThrowsExactly<System.ArgumentException>(() => PlottingPositions.Function(new[] { 0.0, 0.0 }, 0.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, -1.0 }, 0.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, double.NaN }, 0.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, 1.0 }, 1.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, 1.0 }, -0.1));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, 1.0 }, double.NaN));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, double.PositiveInfinity }, 0.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { double.MaxValue, double.MaxValue }, 0.5));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 1.0, 1.0 }, 0.5000001));
            Assert.ThrowsExactly<System.ArgumentOutOfRangeException>(() => PlottingPositions.Function(new[] { 2.0 }, 1.0));
        }

        /// <summary>
        /// Test the weighted plotting positions at their boundaries: a single weight, alpha = 0.5 with a zero last
        /// weight (the largest position the formula allows), and monotonicity within [0, 1) for random weights.
        /// </summary>
        [TestMethod]
        public void Test_WeightedFunction_Boundaries()
        {
            // One weight w: (w − α w) / (w + w − 2α w) = 1/2 for every α.
            foreach (double alpha in new[] { 0.0, 0.44, 0.5 })
                Assert.AreEqual(0.5, PlottingPositions.Function(new[] { 3.0 }, alpha)[0], 1E-15);

            // α = 0.5, weights {2, 2, 0}: W = 4, denominator 4; F = (2 − 1)/4, (4 − 1)/4, (4 − 0)/4. The last
            // position reaches one only in this limit (a zero last weight at α = 1/2).
            CollectionAssert.AreEqual(new[] { 0.25, 0.75, 1.0 }, PlottingPositions.Function(new[] { 2.0, 2.0, 0.0 }, 0.5));

            var random = new System.Random(20261023);
            var weights = new double[500];
            for (int i = 0; i < weights.Length; i++) weights[i] = random.NextDouble() < 0.1 ? 0.0 : random.NextDouble() * 10;
            foreach (double alpha in new[] { 0.0, 0.3175, 0.44, 0.5 })
            {
                var positions = PlottingPositions.Function(weights, alpha);
                Assert.IsGreaterThanOrEqualTo(0.0, positions[0]);
                Assert.IsLessThanOrEqualTo(1.0, positions[positions.Length - 1]);
                for (int i = 1; i < positions.Length; i++)
                    Assert.IsGreaterThanOrEqualTo(positions[i - 1], positions[i], $"alpha {alpha}, i {i}");
            }
        }
    }
}
