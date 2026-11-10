using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Frank Copula. Reference values come from the R 'copula' package for the original methods and from pyvinecopulib 0.7.6 and mpmath 1.4.1 for the conditional-distribution methods.
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
    public class Test_FrankCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new FrankCopula(10);
            Assert.AreEqual(0.02469702, copula.PDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(0.0001228828, copula.PDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(1.089996, copula.PDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(0.2408784, copula.PDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new FrankCopula(10);
            Assert.AreEqual(0.02469702, copula.PDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(0.0001228828, copula.PDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(1.089996, copula.PDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(0.2408784, copula.PDF(0.8, 0.2), 1E-6);
        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new FrankCopula(10);
            Assert.AreEqual(0.1998148, copula.CDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(0.1999997, copula.CDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(0.1534309, copula.CDF(0.2, 0.8), 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(0.1960338, copula.CDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new FrankCopula(10);
            Assert.AreEqual(0.1998148, copula.CDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(20);
            Assert.AreEqual(0.1999997, copula.CDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(-0.5);
            Assert.AreEqual(0.1534309, copula.CDF(0.8, 0.2), 1E-6);

            copula = new FrankCopula(5);
            Assert.AreEqual(0.1960338, copula.CDF(0.8, 0.2), 1E-6);
        }


        private double[] data1 = new double[] { 89.4465116086853, 91.7257825797634, 102.113218235982, 80.9258507082031, 83.4871920273991, 115.695292225839, 82.9627093424983, 88.6038265767081, 128.772275696669, 102.984538236351, 121.578135793465, 108.431406395934, 80.0293061353146, 111.431178502881, 79.2589288833566, 118.158376357189, 80.7620921129063, 108.376156927961, 95.8476352445177, 114.311862750292, 94.7747731262201, 118.875199050171, 93.2349488601981, 98.3959955105056, 108.649741186885, 122.182781815476, 77.0562107320218, 107.540344227161, 117.370406779471, 112.182445319761, 93.0289332852279, 90.529599567024, 120.201762251141, 114.278074578765, 104.568760050718, 99.7636591326382, 110.046851177132, 94.1493960120251, 104.412586917358, 98.7607116866476, 98.7463881469725, 105.1085034505, 121.269695545648, 113.731813064378, 107.800720302073, 106.36957073507, 91.6744161473404, 116.374681447127, 81.77133664558, 85.6156796128899, 107.491138655012, 102.732572344842, 121.586784481517, 116.993480194217, 124.032257210449, 99.1570143649983, 71.7128277313518, 54.6418351602958, 120.014404396034, 96.6470795565503, 107.428703111402, 99.3256773796286, 96.3098320017686, 115.238337503403, 72.0162782248935, 80.5473603316434, 85.8099029975343, 134.325558667819, 102.547186178411, 109.9543113498, 112.218655560109, 103.637934442487, 105.599530792916, 104.374711756055, 117.597527250648, 108.696922523792, 95.7635882353413, 113.556079933849, 91.7625855813032, 101.738578831165, 63.8120024793359, 63.7998925020047, 100.584112662855, 103.281311868372, 87.666605980906, 116.010424048388, 120.789094918399, 115.413149279569, 93.3814849565503, 91.1936357063422, 95.8835908172229, 104.01051720411, 118.117029822129, 101.092380522772, 105.948133484193, 101.511409332803, 91.6808235778427, 112.274907117709, 113.335593237504, 109.199879710188 };
        private double[] data2 = new double[] { 90.2725830142295, 50.9675107814588, 72.8091017436628, 37.834394607052, 66.0298908678443, 131.543666094728, 23.4646646758453, 38.387445258362, 90.0045371339763, 64.0836922901644, 113.713441516048, 86.2362874241136, 68.9809017958243, 83.2228976474751, 40.2100658509047, 126.376008778581, 79.489380523177, 105.171916388245, 92.2811081494569, 111.661707493053, 76.4955819232933, 64.4896701009032, 59.6898595730439, 95.893807000091, 84.416144203673, 119.239043275835, 41.0210587766723, 80.3403252114675, 98.987503709067, 117.127667708957, 75.0679115951468, 47.6762994744396, 98.631369967578, 145.07181315581, 73.404950592094, 77.2596264433046, 85.5789153296777, 58.5136745231103, 90.0353584106097, 88.7890539339681, 82.6531147843648, 85.5315028451663, 107.486065585877, 115.188926317826, 88.5805804644887, 94.0736037456068, 70.9824270336269, 104.904710824485, 56.6492033519975, 49.5998637369, 87.9079770735456, 85.6153352921363, 109.84039662402, 109.70332488243, 109.013245187645, 95.3425794566713, 37.0734004156272, 68.2576366613219, 83.0475074812988, 76.0326431319754, 99.6616980295915, 73.6251299875544, 58.5130416264008, 125.024141331457, 20.7835530562455, 55.6776014119347, 40.421171413744, 109.884427624024, 93.0410150569105, 98.1134492367469, 103.110412789979, 70.3696410672021, 87.9378660295918, 88.828217701758, 145.785916448632, 124.799405281758, 61.6454234189007, 74.1550632448535, 53.7379563429746, 103.133087738044, 64.6916624574593, 21.1435390667273, 77.6876521796421, 75.6302369353579, 71.5323052064531, 160.885271053561, 95.1756368878511, 119.852069907319, 68.0321254425055, 74.3625711976756, 79.8206997832619, 84.3152823622181, 113.17252314479, 54.5297187685205, 82.1517568341444, 87.0940228919095, 69.4555843222865, 81.7696161911127, 73.1591225199505, 93.8908358293613 };

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
            BivariateCopula copula = new FrankCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(7.718761, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new FrankCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            Assert.AreEqual(8.130259, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new FrankCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(8.231434, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(101.4276, ((Normal)copula.MarginalDistributionX).Mu, 1E-2);
            Assert.AreEqual(14.88317, ((Normal)copula.MarginalDistributionX).Sigma, 1E-2);
            // Marginal-Y
            Assert.AreEqual(81.08905, ((Normal)copula.MarginalDistributionY).Mu, 1E-2);
            Assert.AreEqual(26.93537, ((Normal)copula.MarginalDistributionY).Sigma, 1E-2);
        }

        /// <summary>
        /// Test the tail dependence coefficients.
        /// Frank copula: λ_U = λ_L = 0 (no tail dependence).
        /// </summary>
        [TestMethod]
        public void Test_TailDependence()
        {
            var copula = new FrankCopula(5.0);
            Assert.AreEqual(0.0, copula.UpperTailDependence, 1E-10, "Frank copula should have no upper tail dependence.");
            Assert.AreEqual(0.0, copula.LowerTailDependence, 1E-10, "Frank copula should have no lower tail dependence.");

            var copulaNeg = new FrankCopula(-5.0);
            Assert.AreEqual(0.0, copulaNeg.UpperTailDependence, 1E-10, "Frank copula should have no upper tail dependence for negative θ.");
            Assert.AreEqual(0.0, copulaNeg.LowerTailDependence, 1E-10, "Frank copula should have no lower tail dependence for negative θ.");
        }

        /// <summary>
        /// Test that ParametersValid holds across the Frank copula's unbounded θ range.
        /// </summary>
        [TestMethod]
        public void Test_ParametersValid()
        {
            Assert.IsTrue(new FrankCopula(4.2).ParametersValid);
            Assert.IsTrue(new FrankCopula(-7.5).ParametersValid);
        }

        /// <summary>
        /// Test Clone produces an independent copy with deep-copied marginals.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new FrankCopula(4.2, new Normal(100, 10), new Gumbel(50, 5));
            var clone = copula.Clone() as FrankCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(copula.Theta, clone.Theta);
            Assert.IsTrue(clone.ParametersValid);

            // Marginals are deep-copied (distributions memoize lazily, so clones must
            // not share marginal instances), with identical quantiles
            Assert.AreNotSame(copula.MarginalDistributionX, clone.MarginalDistributionX);
            Assert.AreNotSame(copula.MarginalDistributionY, clone.MarginalDistributionY);
            Assert.AreEqual(copula.MarginalDistributionY.InverseCDF(0.9), clone.MarginalDistributionY.InverseCDF(0.9));

            // Mutating the clone should not affect the original
            clone.Theta = 2.0;
            Assert.AreEqual(4.2, copula.Theta);
        }

        /// <summary>
        /// Test the forward conditional CDF (h-function), h(v|u) = ∂C(u,v)/∂u, covering both
        /// the positive- and negative-dependence branches. Reference values were computed with
        /// the Python package pyvinecopulib 0.7.6 (Bicop.hfunc1, family frank), cross-checked
        /// at generation against a 50-digit numerical ∂C/∂u of the Frank CDF computed with
        /// mpmath 1.4.1 (agreement ≤ 1E-15). The analytic conditional is also cross-checked
        /// against a central finite difference of the closed-form CDF (ε = 1E-6; noise floor
        /// ~1E-9, so 1E-7 is asserted) and pinned at the exact upper edge h(1|u) = 1.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDF()
        {
            var copula = new FrankCopula(4d);
            Assert.AreEqual(0.8693978167836881, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.1306021832163123, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.9888149481553941, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(0.0061499217658669445, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            copula = new FrankCopula(-3d);
            Assert.AreEqual(0.5965731714099825, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.4034268285900173, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.755959800956752, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(0.11288571261377746, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            foreach (double theta in new[] { 4d, -3d })
            {
                var c = new FrankCopula(theta);
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
        /// Test the scalar inverse conditional CDF, covering both the positive- and
        /// negative-dependence branches of the reflected closed form. Reference values were
        /// computed with the Python package pyvinecopulib 0.7.6 (Bicop.hinv1, family frank),
        /// whose internal inversion carries ~3E-11 residual, so 1E-8 deltas are asserted for
        /// the pins. The closed-form round trip h(InverseConditionalCDF(u, t) | u) = t is
        /// asserted at 1E-10, and InverseCDF(u, t)[1] must recompose the scalar bit-for-bit.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDF()
        {
            var copula = new FrankCopula(4d);
            Assert.AreEqual(0.5090047384437639, copula.InverseConditionalCDF(0.3, 0.7), 1E-8);
            Assert.AreEqual(0.2597601357556414, copula.InverseConditionalCDF(0.9, 0.05), 1E-8);

            copula = new FrankCopula(-3d);
            Assert.AreEqual(0.7771017148916144, copula.InverseConditionalCDF(0.3, 0.7), 1E-8);
            Assert.AreEqual(0.02170136771746911, copula.InverseConditionalCDF(0.9, 0.05), 1E-8);

            foreach (double theta in new[] { 4d, -3d })
            {
                var c = new FrankCopula(theta);
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
        /// The length of the rank fixtures used by the method of moments tests.
        /// </summary>
        private const int RankFixtureLength = 101;

        /// <summary>
        /// The number of pairs in a rank fixture, which is the denominator of Kendall's tau without ties.
        /// </summary>
        private const double RankFixturePairs = RankFixtureLength * (RankFixtureLength - 1) / 2d;

        /// <summary>
        /// Returns a permutation of 0 to n - 1 carrying exactly the requested number of inversions.
        /// </summary>
        /// <param name="n">The length of the permutation.</param>
        /// <param name="inversions">The number of inverted pairs, between 0 and n * (n - 1) / 2.</param>
        /// <returns>The permutation, as a double array.</returns>
        /// <remarks>
        /// Paired against the ascending sequence 0 to n - 1, a permutation carrying d inversions has exactly d
        /// discordant pairs and no ties, so Kendall's tau of the pair is exactly 1 - 2 d / (n (n - 1) / 2). This
        /// gives a fixture whose tau is a chosen rational value, which is what lets the method of moments fit be
        /// exercised at a prescribed dependence. The permutation is built from its Lehmer code, taking as many
        /// inversions as the remaining positions allow at each step.
        /// </remarks>
        private static double[] PermutationWithInversions(int n, int inversions)
        {
            var available = new List<int>();
            for (int i = 0; i < n; i++) available.Add(i);

            var result = new double[n];
            int remaining = inversions;
            for (int i = 0; i < n; i++)
            {
                int take = Math.Min(remaining, n - 1 - i);
                result[i] = available[take];
                available.RemoveAt(take);
                remaining -= take;
            }
            return result;
        }

        /// <summary>
        /// Returns the ascending sequence 0 to n - 1.
        /// </summary>
        /// <param name="n">The length of the sequence.</param>
        /// <returns>The ascending sequence, as a double array.</returns>
        private static double[] AscendingRanks(int n)
        {
            var result = new double[n];
            for (int i = 0; i < n; i++) result[i] = i;
            return result;
        }

        /// <summary>
        /// Test Kendall's tau as a function of theta against an external oracle.
        /// </summary>
        /// <remarks>
        /// The reference values are from pyvinecopulib 0.7.6, the Python bindings of vinecopulib, which is the
        /// C++ engine behind the R package rvinecopulib. They are an independent implementation of the Frank tau
        /// relation and are not derived from the Debye evaluation used here. The negative rows are the ones that
        /// pin the Debye reflection: statsmodels 0.14.6 omits it and returns 2.7E6 at theta = -30, while it
        /// reproduces the positive rows below to about 3.5E-12.
        /// </remarks>
        [TestMethod]
        public void Test_KendallsTauFromTheta()
        {
            Assert.AreEqual(-0.8739774847415348, FrankCopula.KendallsTauFromTheta(-30.0), 1E-12);
            Assert.AreEqual(-0.6657773862719784, FrankCopula.KendallsTauFromTheta(-10.0), 1E-12);
            Assert.AreEqual(-0.4567009581601168, FrankCopula.KendallsTauFromTheta(-5.0), 1E-12);
            Assert.AreEqual(-0.11001853644899295, FrankCopula.KendallsTauFromTheta(-1.0), 1E-12);
            Assert.AreEqual(0.11001853644899295, FrankCopula.KendallsTauFromTheta(1.0), 1E-12);
            Assert.AreEqual(0.4567009581601168, FrankCopula.KendallsTauFromTheta(5.0), 1E-12);
            Assert.AreEqual(0.6657773862719784, FrankCopula.KendallsTauFromTheta(10.0), 1E-12);
            Assert.AreEqual(0.8739774847415348, FrankCopula.KendallsTauFromTheta(30.0), 1E-12);
        }

        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new FrankCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(7.7385956, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Test that the method of moments fit recovers the theta implied by the sample tau, on both dependency branches.
        /// </summary>
        /// <remarks>
        /// Each fixture has a Kendall's tau that is an exact rational, so the fitted theta must reproduce that
        /// tau when it is pushed back through the tau relation. The tolerance reflects the root-finding
        /// tolerance of Brent's method, not the accuracy of the tau relation.
        /// </remarks>
        [TestMethod]
        public void Test_SetThetaFromTau_RoundTrip()
        {
            var ranks = AscendingRanks(RankFixtureLength);
            foreach (int inversions in new[] { 200, 1200, 2400, 3800, 4800 })
            {
                var permuted = PermutationWithInversions(RankFixtureLength, inversions);
                double expected = 1d - 2d * inversions / RankFixturePairs;
                Assert.AreEqual(expected, Correlation.KendallsTau(ranks, permuted), 1E-12, $"The fixture with {inversions} inversions must carry the expected tau.");

                var copula = new FrankCopula();
                copula.SetThetaFromTau(ranks, permuted);
                Assert.AreEqual(Math.Sign(expected), Math.Sign(copula.Theta), "The fitted theta must take the sign of the sample tau.");
                Assert.IsLessThanOrEqualTo(100d, Math.Abs(copula.Theta));
                Assert.AreEqual(expected, FrankCopula.KendallsTauFromTheta(copula.Theta), 1E-7, $"The fit with {inversions} inversions did not recover the sample tau.");
            }
        }

        /// <summary>
        /// Test that a tau the fitting bracket cannot reach is rejected.
        /// </summary>
        [TestMethod]
        public void Test_SetThetaFromTau_UnattainableTau()
        {
            var ranks = AscendingRanks(RankFixtureLength);
            var copula = new FrankCopula();

            // Perfect concordance and perfect discordance give tau = 1 and tau = -1, both outside the range the
            // bracket theta in [-100, 100] reaches.
            var concordant = PermutationWithInversions(RankFixtureLength, 0);
            var positive = Assert.ThrowsExactly<ArgumentException>(() => copula.SetThetaFromTau(ranks, concordant));
            StringAssert.Contains(positive.Message, "~= [-0.96065, 0.96065]");

            var discordant = PermutationWithInversions(RankFixtureLength, (int)RankFixturePairs);
            var negative = Assert.ThrowsExactly<ArgumentException>(() => copula.SetThetaFromTau(ranks, discordant));
            StringAssert.Contains(negative.Message, "~= [-0.96065, 0.96065]");
        }

        /// <summary>
        /// Test that an independent sample is fitted at the smallest dependency the bracket admits.
        /// </summary>
        /// <remarks>
        /// The Frank independence limit is theta = 0, at which the generator and the distribution functions are
        /// indeterminate, so the fit returns the bracket endpoint nearest independence instead. A tau of exactly
        /// zero takes the negative bracket, matching ParameterConstraints, so the endpoint is -0.001. The
        /// resulting copula must be indistinguishable from independence over the unit square.
        /// </remarks>
        [TestMethod]
        public void Test_SetThetaFromTau_Independence()
        {
            var ranks = AscendingRanks(RankFixtureLength);
            var permuted = PermutationWithInversions(RankFixtureLength, (int)RankFixturePairs / 2);
            Assert.AreEqual(0d, Correlation.KendallsTau(ranks, permuted), 0d, "Half of the pairs discordant gives a tau of exactly zero.");

            var copula = new FrankCopula();
            copula.SetThetaFromTau(ranks, permuted);
            Assert.AreEqual(-0.001d, copula.Theta, 0d);
            Assert.AreEqual(0d, FrankCopula.KendallsTauFromTheta(copula.Theta), 1E-3);
            foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
            {
                foreach (double v in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    Assert.AreEqual(u * v, copula.CDF(u, v), 1E-4, $"The near-independent fit must reproduce the independence copula at ({u}, {v}).");
                }
            }
        }

        /// <summary>
        /// LogPDF and PDF against the closed-form log density evaluated in mpmath 1.3.0 at 120 digits, for θ from -500 to 500
        /// and u, v from 1E-12 to 1 - 1E-12, where the direct density underflows or overflows. The log density is within 1E-12
        /// (relative above 1 in magnitude) and PDF is its exponential.
        /// </summary>
        [TestMethod]
        public void Test_LogPDFMatchesMpmath()
        {
            double[,] cases =
            {
            { -100.0, 1e-12, 1e-12, -95.39482981381191 },
            { -100.0, 1e-12, 0.5, -45.39482981391191 },
            { -100.0, 1e-12, 0.999999999999, 4.605170185788094 },
            { -100.0, 0.01, 0.99, 3.6254099346985913 },
            { -100.0, 0.2, 0.8, 3.2188758269293545 },
            { -100.0, 0.5, 0.5, 3.2188758248682006 },
            { -100.0, 0.999999999999, 0.5, -45.39482981391191 },
            { -100.0, 0.999999999999, 0.999999999999, -95.39482981381191 },
            { -5.0, 1e-12, 1e-12, -3.383801338106411 },
            { -5.0, 1e-12, 0.5, -0.8838013381121697 },
            { -5.0, 1e-12, 0.999999999999, 1.616198661873589 },
            { -5.0, 0.01, 0.99, 1.5209938145165462 },
            { -5.0, 0.2, 0.8, 0.6454877727941052 },
            { -5.0, 0.5, 0.5, 0.38768376934879756 },
            { -5.0, 0.999999999999, 0.5, -0.8838013381121698 },
            { -5.0, 0.999999999999, 0.999999999999, -3.3838013381064114 },
            { 5.0, 1e-12, 1e-12, 1.6161986618735888 },
            { 5.0, 1e-12, 0.5, -0.8838013381121697 },
            { 5.0, 1e-12, 0.999999999999, -3.3838013381064114 },
            { 5.0, 0.01, 0.99, -3.283837002529299 },
            { 5.0, 0.2, 0.8, -1.4234628465060422 },
            { 5.0, 0.5, 0.5, 0.38768376934879756 },
            { 5.0, 0.999999999999, 0.5, -0.8838013381121698 },
            { 5.0, 0.999999999999, 0.999999999999, 1.616198661873589 },
            { 100.0, 1e-12, 1e-12, 4.605170185788092 },
            { 100.0, 1e-12, 0.5, -45.39482981391191 },
            { 100.0, 1e-12, 0.999999999999, -95.39482981381191 },
            { 100.0, 0.01, 0.99, -93.39482981401191 },
            { 100.0, 0.2, 0.8, -55.39482981401191 },
            { 100.0, 0.5, 0.5, 3.2188758248682006 },
            { 100.0, 0.999999999999, 0.5, -45.39482981391191 },
            { 100.0, 0.999999999999, 0.999999999999, 4.605170185788096 },
            { 500.0, 1e-12, 1e-12, 6.214608097422191 },
            { 500.0, 1e-12, 0.5, -243.78539190107782 },
            { 500.0, 1e-12, 0.999999999999, -493.7853919005778 },
            { 500.0, 0.01, 0.99, -483.7853919015778 },
            { 500.0, 0.2, 0.8, -293.78539190157784 },
            { 500.0, 0.5, 0.5, 4.8283137373023015 },
            { 500.0, 0.999999999999, 0.5, -243.78539190107782 },
            { 500.0, 0.999999999999, 0.999999999999, 6.2146080974222135 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new FrankCopula(cases[i, 0]);
                double u = cases[i, 1], v = cases[i, 2], expected = cases[i, 3];
                double actual = copula.LogPDF(u, v);
                if (double.IsNegativeInfinity(expected))
                {
                    Assert.IsTrue(double.IsNegativeInfinity(actual), $"case {i}");
                    Assert.AreEqual(0d, copula.PDF(u, v), 0d, $"case {i}");
                    continue;
                }
                Assert.AreEqual(expected, actual, 1E-12 * Math.Max(1d, Math.Abs(expected)), $"case {i}: theta={cases[i, 0]}, u={u:R}, v={v:R}");
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(Math.Exp(actual)), BitConverter.DoubleToInt64Bits(copula.PDF(u, v)), $"case {i}");
            }
        }

        /// <summary>
        /// θ = 0 is the independence limit: PDF 1, LogPDF 0, CDF uv, ConditionalCDF v and InverseConditionalCDF t, on the
        /// closed unit square.
        /// </summary>
        [TestMethod]
        public void Test_ThetaZeroIsIndependence()
        {
            var copula = new FrankCopula(0d);
            double[] grid = { 0d, 1E-12, 0.01, 0.3, 0.5, 0.7, 0.99, 1d - 1E-12, 1d };
            foreach (double u in grid)
            {
                foreach (double v in grid)
                {
                    Assert.AreEqual(1d, copula.PDF(u, v), 0d);
                    Assert.AreEqual(0d, copula.LogPDF(u, v), 0d);
                    Assert.AreEqual(u * v, copula.CDF(u, v), 0d);
                    Assert.AreEqual(v, copula.ConditionalCDF(u, v), 0d);
                    Assert.AreEqual(v, copula.InverseConditionalCDF(u, v), 0d);
                }
            }
        }

        /// <summary>
        /// The members are continuous through θ = 0: at θ = ±1E-12 and ±1E-9 they agree with mpmath 1.3.0 at 80 digits to
        /// 5E-14 (relative for CDF, ConditionalCDF and InverseConditionalCDF, absolute for LogPDF), and so with the θ = 0
        /// limits to within θ times a bounded factor.
        /// </summary>
        [TestMethod]
        public void Test_ContinuousThroughThetaZero()
        {
            // theta, u, v, then mpmath CDF(u, v), LogPDF(u, v), ConditionalCDF(u, v) and InverseConditionalCDF(u, v).
            double[,] cases =
            {
            { -1e-09, 1e-12, 0.5, 4.99999999875e-13, -4.166666666641667e-20, 0.499999999875, 0.500000000125 },
            { -1e-09, 0.3, 0.7, 0.20999999997795, 8.000000000243333e-11, 0.699999999958, 0.7000000000419999 },
            { -1e-09, 0.99, 0.01, 0.009899999999950995, 4.801999999584313e-10, 0.010000000004851, 0.009999999995149 },
            { -1e-09, 0.5, 1e-06, 4.999999998750001e-07, -4.1666416666916673e-20, 1e-06, 1e-06 },
            { -1e-12, 1e-12, 0.5, 4.99999999999875e-13, -4.1666666666416664e-26, 0.499999999999875, 0.500000000000125 },
            { -1e-12, 0.3, 0.7, 0.20999999999997793, 8.000000000000242e-14, 0.699999999999958, 0.7000000000000419 },
            { -1e-12, 0.99, 0.01, 0.00989999999999995, 4.801999999999584e-13, 0.01000000000000485, 0.00999999999999515 },
            { -1e-12, 0.5, 1e-06, 4.999999999998749e-07, -4.166641666691667e-26, 1e-06, 1e-06 },
            { 1e-12, 1e-12, 0.5, 5.00000000000125e-13, -4.1666666666416664e-26, 0.500000000000125, 0.499999999999875 },
            { 1e-12, 0.3, 0.7, 0.21000000000002203, -7.999999999999755e-14, 0.7000000000000419, 0.699999999999958 },
            { 1e-12, 0.99, 0.01, 0.00990000000000005, -4.802000000000415e-13, 0.00999999999999515, 0.01000000000000485 },
            { 1e-12, 0.5, 1e-06, 5.00000000000125e-07, -4.166641666691667e-26, 1e-06, 1e-06 },
            { 1e-09, 1e-12, 0.5, 5.00000000125e-13, -4.166666666641667e-20, 0.500000000125, 0.499999999875 },
            { 1e-09, 0.3, 0.7, 0.21000000002204997, -7.999999999756666e-11, 0.7000000000419999, 0.699999999958 },
            { 1e-09, 0.99, 0.01, 0.009900000000049005, -4.802000000415687e-10, 0.009999999995149, 0.010000000004851 },
            { 1e-09, 0.5, 1e-06, 5.000000001249998e-07, -4.1666416666916673e-20, 1e-06, 1e-06 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new FrankCopula(cases[i, 0]);
                double u = cases[i, 1], v = cases[i, 2];
                Assert.AreEqual(cases[i, 3], copula.CDF(u, v), 5E-14 * cases[i, 3], $"case {i}: CDF");
                Assert.AreEqual(cases[i, 4], copula.LogPDF(u, v), 5E-14, $"case {i}: LogPDF");
                Assert.AreEqual(cases[i, 5], copula.ConditionalCDF(u, v), 5E-14 * cases[i, 5], $"case {i}: ConditionalCDF");
                Assert.AreEqual(cases[i, 6], copula.InverseConditionalCDF(u, v), 5E-14 * cases[i, 6], $"case {i}: InverseConditionalCDF");
                // Within 1E-12 of independence the values sit within about 30·|θ| of the limits.
                Assert.AreEqual(u * v, copula.CDF(u, v), 1E-7 * u * v);
                Assert.AreEqual(0d, copula.LogPDF(u, v), 1E-7);
            }
        }

        /// <summary>
        /// Every generator member against mpmath 1.3.0 at 150 digits, for φ(t) = −ln((e^(−θt) − 1)/(e^(−θ) − 1)), with the independence limit −ln t at θ = 0: Generator, GeneratorInverse (at
        /// the rounded generator value), GeneratorPrime, GeneratorPrime2 and GeneratorPrimeInverse (at the rounded slope),
        /// on t from 1E-12 to 1 - 1E-12 and θ including the independence region, to 1E-13 relative.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorMembersMatchMpmath()
        {
            // theta, t, phi(t), phi^-1(phi), phi'(t), phi''(t), (phi')^-1(phi').
            double[,] cases =
            {
            { -50.0, 1e-12, 73.7189981104754, 1.0000000000000052e-12, -1000000000025.0, 1e+24, 1e-12 },
            { -50.0, 0.01, 50.43275212956719, 0.009999999999999978, -127.07470412683992, 9794.245222581909, 0.01 },
            { -50.0, 0.5, 25.000000000013888, 0.5, -50.0000000006944, 3.471985966337443e-08, 0.4999999420636232 },
            { -50.0, 0.99, 0.5000000000000004, 0.99, -50.0, 7.949927250494377e-19, double.PositiveInfinity },
            { -50.0, 0.999999999999, 4.9998893913993925e-11, 0.999999999999, -50.0, 4.821874620150883e-19, double.PositiveInfinity },
            { -1e-12, 1e-12, 27.631021115929048, 1e-12, -1000000000000.0, 1e+24, 1e-12 },
            { -1e-12, 0.01, 4.6051701859885865, 0.009999999999999998, -100.0000000000005, 10000.0, 0.01 },
            { -1e-12, 0.5, 0.6931471805601953, 0.5, -2.0000000000005, 4.0, 0.5 },
            { -1e-12, 0.99, 0.01005033585350645, 0.99, -1.0101010101015102, 1.0203040506070808, 0.9899999999999999 },
            { -1e-12, 0.999999999999, 9.999778782808785e-13, 0.999999999999, -1.0000000000015, 1.000000000002, 0.9999999999990001 },
            { 0.0, 1e-12, 27.631021115928547, 1.000000000000001e-12, -1000000000000.0, 1e+24, 1e-12 },
            { 0.0, 0.01, 4.605170185988091, 0.010000000000000004, -100.0, 10000.0, 0.01 },
            { 0.0, 0.5, 0.6931471805599453, 0.5, -2.0, 4.0, 0.5 },
            { 0.0, 0.99, 0.01005033585350145, 0.99, -1.0101010101010102, 1.0203040506070808, 0.99 },
            { 0.0, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.000000000001, 1.000000000002, 0.9999999999989999 },
            { 1e-12, 1e-12, 27.63102111592805, 9.999999999999984e-13, -1000000000000.0, 1e+24, 1e-12 },
            { 1e-12, 0.01, 4.605170185987596, 0.010000000000000002, -99.9999999999995, 10000.0, 0.01 },
            { 1e-12, 0.5, 0.6931471805596953, 0.5, -1.9999999999995, 4.0, 0.5 },
            { 1e-12, 0.99, 0.01005033585349645, 0.99, -1.0101010101005101, 1.0203040506070808, 0.99 },
            { 1e-12, 0.999999999999, 9.999778782798785e-13, 0.999999999999, -1.0000000000005, 1.000000000002, 0.9999999999989999 },
            { 5.0, 1e-12, 26.014822454047458, 1.0000000000000014e-12, -999999999997.5, 1e+24, 1e-12 },
            { 5.0, 0.01, 3.0138673596078887, 0.009999999999999998, -97.52083246532945, 9997.9169270575, 0.01 },
            { 5.0, 0.5, 0.07888973429254963, 0.5, -0.44712744916926006, 2.4355602016469096, 0.5 },
            { 5.0, 0.99, 0.0003478659228176354, 0.99, -0.03566970777178796, 0.17962086691146456, 0.99 },
            { 5.0, 0.999999999999, 3.3917524201030334e-14, 0.999999999999, -0.03391827453169189, 0.17074182200566668, 0.999999999999 },
            { 50.0, 1e-12, 23.7189981105254, 1.0000000000000008e-12, -999999999975.0, 1e+24, 1e-12 },
            { 50.0, 0.01, 0.9327521295671886, 0.01, -77.07470412683992, 9794.245222581909, 0.01 },
            { 50.0, 0.5, 1.3887943864867583e-11, 0.5, -6.943971932578447e-10, 3.471985966337443e-08, 0.5 },
            { 50.0, 0.99, 1.2512210522338331e-22, 0.99, -1.5899854500988755e-20, 7.949927250494377e-19, 0.99 },
            { 50.0, 0.999999999999, 9.643535903739066e-33, 0.999999999999, -9.643749240301766e-21, 4.821874620150883e-19, 0.999999999999 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new FrankCopula(cases[i, 0]);
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
        /// The generator slope where θt is subnormal, φ′(t) = −θ/Expm1(θt) = −1/t + θ/2 to double precision (mpmath
        /// 1.3.0 at 50 digits gives −1E300 for θ = 1E-14 and t = 1E-300), and its θ = 0 limit −1/t.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorPrimeSubnormalProduct()
        {
            Assert.AreEqual(-1E300, new FrankCopula(1E-14).GeneratorPrime(1E-300), 1E-15 * 1E300);
            Assert.AreEqual(-1E300, new FrankCopula(-1E-14).GeneratorPrime(1E-300), 1E-15 * 1E300);
            Assert.AreEqual(-1d / 0.3d, new FrankCopula(0d).GeneratorPrime(0.3d), 0d);
        }
    }
}
