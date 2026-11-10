using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Distributions.Copulas;
using Numerics.Mathematics.RootFinding;

namespace Distributions.BivariateCopulas
{
    /// <summary>
    /// Unit tests for the Gumbel Copula. Reference values come from the R 'copula' package for the original methods and from pyvinecopulib 0.7.6 and mpmath 1.4.1 for the conditional-distribution methods.
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
    public class Test_GumbelCopula
    {
        /// <summary>
        /// Test PDF
        /// </summary>
        [TestMethod]
        public void Test_PDF()
        {
            var copula = new GumbelCopula(10);
            Assert.AreEqual(1.559981e-07, copula.PDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(7.954156e-16, copula.PDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(0.9004251, copula.PDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(0.001609717, copula.PDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new GumbelCopula(10);
            Assert.AreEqual(1.559981e-07, copula.PDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(7.954156e-16, copula.PDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(0.9004251, copula.PDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(0.001609717, copula.PDF(0.8, 0.2), 1E-6);

        }

        /// <summary>
        /// Test CDF
        /// </summary>
        [TestMethod]
        public void Test_CDF()
        {
            var copula = new GumbelCopula(10);
            Assert.AreEqual(0.2, copula.CDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(0.2, copula.CDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(0.1694668, copula.CDF(0.2, 0.8), 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(0.1999967, copula.CDF(0.2, 0.8), 1E-6);

            // Test symmetry
            copula = new GumbelCopula(10);
            Assert.AreEqual(0.2, copula.CDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(20);
            Assert.AreEqual(0.2, copula.CDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(1.1);
            Assert.AreEqual(0.1694668, copula.CDF(0.8, 0.2), 1E-6);

            copula = new GumbelCopula(5);
            Assert.AreEqual(0.1999967, copula.CDF(0.8, 0.2), 1E-6);

        }

        private double[] data1 = new double[] { 112.472889925526, 108.197032461095, 113.64666688938, 113.965212991957, 117.333562148345, 100.559345907851, 102.180053944561, 110.97067894608, 104.077049999534, 102.473973961948, 107.113945888296, 106.082540016705, 87.0978719961949, 75.1150155332284, 112.802882878613, 138.120743213339, 79.1789109173273, 106.623115297302, 110.92063481489, 109.543820203885, 73.138959971343, 98.3431545632918, 102.952511132798, 93.0465814076382, 111.688341287489, 79.2581834425453, 103.249260938183, 85.8316965436919, 84.6132137548954, 110.962644400059, 84.0087454428379, 107.223748328677, 121.506255972552, 123.73645847787, 79.0481183395772, 146.264725342603, 95.3276983471542, 107.603210240998, 99.9974503587899, 90.8379276433018, 104.67679932997, 73.4569986162982, 113.260163491206, 97.0605626966906, 89.0374440090811, 133.729158410007, 78.9190177023254, 108.345190856042, 108.953076977455, 97.0704392084539, 86.2765535028879, 102.219521216094, 63.2511232136125, 76.7529568677982, 74.1684311723453, 84.3685879177063, 94.7059136752836, 92.5916995697557, 107.923148783474, 105.261642868879, 84.7155767681601, 105.132628194701, 102.500956437438, 88.5763924520429, 95.9492723947193, 131.375734529198, 123.260221112319, 73.5681907570129, 93.1887735334748, 106.007277431115, 100.356675515534, 91.4930966213857, 122.640033750305, 107.143857399285, 91.0914154936143, 69.1482456631603, 122.434236175527, 111.906044824347, 111.632367993799, 84.6813208809395, 99.4489764877588, 106.884486176956, 95.3274691248045, 107.839259276592, 90.3967465603013, 76.5485833533036, 98.4847265047609, 96.2342050639354, 91.963420241494, 101.351308748786, 95.4725609428008, 145.010119874994, 88.3665385050844, 109.704283131981, 81.2884314113531, 109.836099370176, 102.386450396465, 93.0092098477109, 83.7158456419459, 89.7486244037052 };
        private double[] data2 = new double[] { 101.222248773938, 81.7313758478523, 123.411631899952, 93.2413907814405, 102.50494017697, 96.5306123268502, 18.1227322903384, 74.307009995493, 82.2705528425937, 89.6693353083137, 87.2256413770307, 81.1744664866471, 81.2165207417765, 81.4927985171313, 100.651709224046, 130.599135003894, 45.4442163124226, 66.6919540534553, 83.5881809451671, 88.4770416926437, 75.1039197272015, 52.1553833221993, 49.528146406319, 66.3900605051234, 97.2816608378485, 26.5831021290638, 88.5107383079526, 78.4085272967254, 70.7569274155408, 109.703373315445, 87.4165170640218, 85.6182893046897, 110.738685438064, 119.768351328997, 21.6161419679033, 154.357497576718, 53.8795395966417, 75.943071862762, 105.35244272571, 86.5868773757059, 76.4752235253919, 48.5832730671902, 84.6605006576596, 77.1067893953142, 46.6141185543174, 126.161273608663, 67.6305173666572, 94.9907772124782, 112.133285294602, 72.9506266499651, 37.128394091496, 63.4492833651505, 59.0759448691672, 56.0401720089048, 74.4223335518519, 43.9584062661215, 64.1266855983013, 103.095921282224, 107.649317301385, 72.9193080139833, 91.8847445373411, 75.1454442126996, 76.9540119546376, 64.9699411254905, 55.4012215444807, 118.83302426895, 106.302264714345, 48.1415318738189, 62.2914940564682, 84.5812336961348, 80.3688979491351, 47.5308700972038, 107.759633266225, 95.6723525891828, 27.2716316191787, 99.1054148896251, 82.2568018516306, 96.0432032275742, 96.3923733319187, 62.287551336599, 47.541023973346, 95.9802695805108, 80.2482799857865, 88.743034379118, 51.9583360652995, 73.0816610190799, 33.729130931116, 55.4771104609277, 93.3789037599054, 90.0234978870232, 62.3635712815712, 146.062306846747, 54.4992333698561, 89.1013264737916, 80.1329566564616, 123.235472480686, 73.4015267834993, 49.6295804684941, 74.6326460427444, 64.7088795224836 };


        /// <summary>
        /// Estimate using the method of moments.
        /// </summary>
        [TestMethod]
        public void Test_MOM_Fit()
        {
            var copula = new GumbelCopula();
            copula.SetThetaFromTau(data1, data2);
            Assert.AreEqual(1.978417, copula.Theta, 1E-4);
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
            BivariateCopula copula = new GumbelCopula();
            BivariateCopulaEstimation.Estimate(ref copula, rank1, rank2, CopulaEstimationMethod.PseudoLikelihood);
            Assert.AreEqual(2.097953, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the inference from margins method.
        /// </summary>
        [TestMethod]
        public void Test_IFM_Fit()
        {
            BivariateCopula copula = new GumbelCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit marginals
            ((IEstimation)copula.MarginalDistributionX).Estimate(data1, ParameterEstimationMethod.MaximumLikelihood);
            ((IEstimation)copula.MarginalDistributionY).Estimate(data2, ParameterEstimationMethod.MaximumLikelihood);
            // Fit copula
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.InferenceFromMargins);
            // IFM with exact Normal MLE margins (sample mean, population SD) maximized in Python 3.14 with scipy 1.17.1
            // (scipy.optimize.minimize_scalar on the Gumbel copula log-likelihood of the margin CDFs): theta = 2.0311985.
            Assert.AreEqual(2.0311985, copula.Theta, 1E-4);
        }

        /// <summary>
        /// Estimate using the method of maximum likelihood estimation.
        /// </summary>
        [TestMethod]
        public void Test_MLE_Fit()
        {
            BivariateCopula copula = new GumbelCopula();
            copula.MarginalDistributionX = new Normal();
            copula.MarginalDistributionY = new Normal();
            // Fit copula and marginals
            BivariateCopulaEstimation.Estimate(ref copula, data1, data2, CopulaEstimationMethod.FullLikelihood);
            // Theta
            Assert.AreEqual(1.999167, copula.Theta, 1E-3);
            // Marginal-X
            Assert.AreEqual(99.90953, ((Normal)copula.MarginalDistributionX).Mu, 1E-1);
            Assert.AreEqual(15.80668, ((Normal)copula.MarginalDistributionX).Sigma, 1E-1);
            // Marginal-Y
            Assert.AreEqual(78.36136, ((Normal)copula.MarginalDistributionY).Mu, 1E-1);
            Assert.AreEqual(24.83673, ((Normal)copula.MarginalDistributionY).Sigma, 1E-1);
        }

        /// <summary>
        /// Test the tail dependence coefficients.
        /// Gumbel copula: λ_U = 2 - 2^(1/θ), λ_L = 0.
        /// </summary>
        [TestMethod]
        public void Test_TailDependence()
        {
            // θ = 2: λ_U = 2 - 2^(1/2) ≈ 0.58579
            var copula = new GumbelCopula(2.0);
            Assert.AreEqual(0.0, copula.LowerTailDependence, 1E-10, "Gumbel copula should have no lower tail dependence.");
            Assert.AreEqual(2.0 - Math.Pow(2.0, 0.5), copula.UpperTailDependence, 1E-6, "λ_U should equal 2 - 2^(1/θ) for θ=2.");

            // θ = 4: λ_U = 2 - 2^(1/4) ≈ 0.81079
            var copula4 = new GumbelCopula(4.0);
            Assert.AreEqual(2.0 - Math.Pow(2.0, 0.25), copula4.UpperTailDependence, 1E-6, "λ_U should equal 2 - 2^(1/θ) for θ=4.");

            // θ = 1 (independence): λ_U = 2 - 2 = 0
            var copulaIndep = new GumbelCopula(1.0);
            Assert.AreEqual(0.0, copulaIndep.UpperTailDependence, 1E-10, "λ_U should be 0 for θ=1 (independence).");

            // θ → ∞: λ_U → 1
            var copulaLarge = new GumbelCopula(100.0);
            Assert.IsGreaterThan(0.99, copulaLarge.UpperTailDependence, "λ_U should approach 1 for large θ.");
        }

        /// <summary>
        /// Test that ParametersValid tracks the dependency parameter's valid range.
        /// ValidateParameter returns null for a valid dependence parameter and a non-null
        /// exception otherwise, so ParametersValid follows the parameter.
        /// </summary>
        [TestMethod]
        public void Test_ParametersValid()
        {
            var copula = new GumbelCopula(2.0);
            Assert.IsTrue(copula.ParametersValid);

            copula.Theta = 0.5; // below the minimum of 1
            Assert.IsFalse(copula.ParametersValid);

            copula.Theta = 3.0;
            Assert.IsTrue(copula.ParametersValid);
        }

        /// <summary>
        /// Test Clone produces an independent copy with deep-copied marginals.
        /// </summary>
        [TestMethod]
        public void Test_Clone()
        {
            var copula = new GumbelCopula(2.0, new Normal(100, 10), new Gumbel(50, 5));
            var clone = copula.Clone() as GumbelCopula;
            Assert.IsNotNull(clone);
            Assert.AreEqual(copula.Theta, clone.Theta);
            Assert.IsTrue(clone.ParametersValid);

            // Marginals are deep-copied (distributions memoize lazily, so clones must
            // not share marginal instances), with identical quantiles
            Assert.AreNotSame(copula.MarginalDistributionX, clone.MarginalDistributionX);
            Assert.AreNotSame(copula.MarginalDistributionY, clone.MarginalDistributionY);
            Assert.AreEqual(copula.MarginalDistributionY.InverseCDF(0.9), clone.MarginalDistributionY.InverseCDF(0.9));

            // Mutating the clone should not affect the original
            clone.Theta = 3.0;
            Assert.AreEqual(2.0, copula.Theta);
        }

        /// <summary>
        /// Test the forward conditional CDF (h-function), h(v|u) = ∂C(u,v)/∂u.
        /// Reference values were computed with the Python package pyvinecopulib 0.7.6
        /// (Bicop.hfunc1, family gumbel), cross-checked at generation against a 50-digit
        /// numerical ∂C/∂u of the Gumbel CDF computed with mpmath 1.4.1 (agreement ≤ 1E-15).
        /// The analytic conditional is also cross-checked against a central finite difference
        /// of the closed-form CDF (ε = 1E-6; noise floor ~1E-9, so 1E-7 is asserted) and
        /// pinned at the exact upper edge h(1|u) = 1.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDF()
        {
            var copula = new GumbelCopula(2d);
            Assert.AreEqual(0.9104803864754554, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.11559784394154599, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.9975327563905108, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(0.001949079483034616, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            copula = new GumbelCopula(3.5);
            Assert.AreEqual(0.9852294418984343, copula.ConditionalCDF(0.3, 0.7), 1E-10);
            Assert.AreEqual(0.0201697450039637, copula.ConditionalCDF(0.7, 0.3), 1E-10);
            Assert.AreEqual(0.9999871895213038, copula.ConditionalCDF(0.05, 0.9), 1E-10);
            Assert.AreEqual(1.2887217392734802E-05, copula.ConditionalCDF(0.9, 0.05), 1E-10);

            foreach (double theta in new[] { 2d, 3.5 })
            {
                var c = new GumbelCopula(theta);
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
        /// Test the scalar inverse conditional CDF. Reference values were computed with the
        /// Python package pyvinecopulib 0.7.6 (Bicop.hinv1, family gumbel). The Gumbel inverse
        /// is a Brent solve with tolerance 1E-8, so 1E-6 absolute deltas are asserted for the
        /// reference pins and for the round trip h(InverseConditionalCDF(u, t) | u) = t.
        /// InverseCDF(u, t)[1] must recompose the scalar bit-for-bit.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDF()
        {
            var copula = new GumbelCopula(2d);
            Assert.AreEqual(0.48403043854946104, copula.InverseConditionalCDF(0.3, 0.7), 1E-6);
            Assert.AreEqual(0.39821464719708766, copula.InverseConditionalCDF(0.9, 0.05), 1E-6);

            copula = new GumbelCopula(3.5);
            Assert.AreEqual(0.39761414717848886, copula.InverseConditionalCDF(0.3, 0.7), 1E-6);
            Assert.AreEqual(0.7272049794715586, copula.InverseConditionalCDF(0.9, 0.05), 1E-6);

            foreach (double theta in new[] { 2d, 3.5 })
            {
                var c = new GumbelCopula(theta);
                foreach (double u in new[] { 0.1, 0.3, 0.5, 0.7, 0.9 })
                {
                    foreach (double t in new[] { 0.1, 0.5, 0.9 })
                    {
                        double v = c.InverseConditionalCDF(u, t);
                        Assert.AreEqual(t, c.ConditionalCDF(u, v), 1E-6);
                        Assert.AreEqual(v, c.InverseCDF(u, t)[1], 0d);
                    }
                }
            }
        }

        /// <summary>
        /// Test the inverse conditional at conditional probabilities within rounding distance
        /// of the boundaries. In exact arithmetic h(1|u) = 1, but the floating-point objective
        /// rounds below 1 by roughly |ln u| ulps, so a requested level of t = 1 − 1E-16 can exceed the
        /// attainable maximum of the objective for u values with |ln u| ≳ 1; the solver must still
        /// complete. The asserted contract is the inverse relationship, not a
        /// particular value: every call returns a probability in [0, 1] whose value is monotone
        /// in the conditional level. That phrasing is deliberate: in the far tail under
        /// strong dependence the conditional CDF is numerically saturated across a band of v,
        /// so several values satisfy h(v|u) = 1 − 1E-16 to full double precision and the
        /// solver may legitimately return any of them. The asserted contract is therefore
        /// completion, range, and monotonicity in the conditional level, deliberately not a
        /// round trip through <see cref="BivariateCopula.ConditionalCDF(double, double)"/>, which the
        /// Archimedean families evaluate through the generator ratio φ′(u)/φ′(C(u,v)) rather
        /// than the closed form this inverse solves, and whose precision degrades at v = 1 for
        /// large θ. Interior levels are pinned against reference values above, which is where
        /// the round trip belongs.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDF_BoundaryConditionals()
        {
            foreach (double theta in new[] { 1.2, 2d, 3.5, 10d, 20d })
            {
                var copula = new GumbelCopula(theta);
                for (int i = 1; i < 2000; i++)
                {
                    double u = i / 2000d;
                    double top = copula.InverseConditionalCDF(u, 1d - 1E-16);
                    Assert.IsTrue(top >= 0d && top <= 1d,
                        $"Top-edge inverse left [0, 1] at θ = {theta}, u = {u}: {top}.");
                    Assert.IsGreaterThanOrEqualTo(copula.InverseConditionalCDF(u, 0.999), top,
                        $"The top-edge inverse must not fall below an interior level at θ = {theta}, u = {u}.");
                    double bottom = copula.InverseConditionalCDF(u, 1E-16);
                    Assert.IsTrue(bottom >= 0d && bottom <= 1d,
                        $"Bottom-edge inverse left [0, 1] at θ = {theta}, u = {u}: {bottom}.");
                    Assert.IsLessThanOrEqualTo(copula.InverseConditionalCDF(u, 0.001), bottom,
                        $"The bottom-edge inverse must not exceed an interior level at θ = {theta}, u = {u}.");
                }
            }
        }

        /// <summary>
        /// The generator derivatives match central differences of φ(t) = (−ln t)^θ and have the signs of a
        /// decreasing convex generator: φ′ &lt; 0 and φ″ &gt; 0 on (0, 1).
        /// </summary>
        [TestMethod]
        public void Test_GeneratorDerivatives()
        {
            foreach (double theta in new[] { 1d, 1.5d, 2d, 5d, 20d })
            {
                var copula = new GumbelCopula(theta);
                foreach (double t in new[] { 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99 })
                {
                    double h = 1E-5 * t;
                    double d1 = (copula.Generator(t + h) - copula.Generator(t - h)) / (2d * h);
                    double d2 = (copula.GeneratorPrime(t + h) - copula.GeneratorPrime(t - h)) / (2d * h);
                    double a = -Math.Log(t);
                    Assert.AreEqual(-theta * Math.Pow(a, theta - 1d) / t, copula.GeneratorPrime(t), 1E-12 * Math.Abs(copula.GeneratorPrime(t)));
                    Assert.AreEqual(d1, copula.GeneratorPrime(t), 1E-4 * Math.Abs(d1));
                    Assert.AreEqual(d2, copula.GeneratorPrime2(t), 1E-4 * Math.Abs(d2));
                    Assert.IsLessThan(0d, copula.GeneratorPrime(t));
                    if (theta > 1d) Assert.IsGreaterThan(0d, copula.GeneratorPrime2(t));
                }
            }
        }

        /// <summary>
        /// GeneratorPrimeInverse inverts GeneratorPrime on (0, 1).
        /// </summary>
        [TestMethod]
        public void Test_GeneratorPrimeInverse_RoundTrip()
        {
            foreach (double theta in new[] { 1.2d, 2d, 5d, 20d })
            {
                var copula = new GumbelCopula(theta);
                foreach (double t in new[] { 0.05, 0.2, 0.5, 0.8, 0.95 })
                    Assert.AreEqual(t, copula.GeneratorPrimeInverse(copula.GeneratorPrime(t)), 1E-6);
            }
        }

        /// <summary>
        /// GeneratorPrimeInverse round-trips GeneratorPrime to 1E-13 relative on t from 1E-12 to 1 - 1E-12, for θ from 1
        /// (independence) to 50, wherever the slope is a normal double.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorPrimeInverse_RoundTripFullPrecision()
        {
            double[] grid = { 1E-12, 1E-6, 0.001, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999, 1d - 1E-6, 1d - 1E-12 };
            foreach (double theta in new[] { 1d, 1.0001d, 1.5d, 2d, 5d, 20d, 50d })
            {
                var copula = new GumbelCopula(theta);
                foreach (double u in grid)
                {
                    double slope = copula.GeneratorPrime(u);
                    // At θ = 50 and u = 1 - 1E-12 the slope, about -5E-587, underflows to zero and has no inverse.
                    if (!(slope < -1E-300)) continue;
                    Assert.AreEqual(u, copula.GeneratorPrimeInverse(slope), 1E-13 * u, $"theta={theta}, u={u}");
                }
            }
        }

        /// <summary>
        /// LogPDF and PDF against the closed-form log density evaluated in mpmath 1.3.0 at 120 digits, for θ from 1 to 100
        /// and u, v from 1E-12 to 1 - 1E-12, where the direct density underflows or overflows. The log density is within 1E-12
        /// (relative above 1 in magnitude) and PDF is its exponential.
        /// </summary>
        [TestMethod]
        public void Test_LogPDFMatchesMpmath()
        {
            double[,] cases =
            {
            { 1.0, 1e-12, 1e-12, 0.0 },
            { 1.0, 1e-12, 0.5, 0.0 },
            { 1.0, 1e-12, 0.999999999999, 0.0 },
            { 1.0, 0.01, 0.99, 0.0 },
            { 1.0, 0.2, 0.8, 0.0 },
            { 1.0, 0.5, 0.5, 0.0 },
            { 1.0, 0.999999999999, 0.5, 0.0 },
            { 1.0, 0.999999999999, 0.999999999999, 0.0 },
            { 2.0, 1e-12, 1e-12, 15.517999325301986 },
            { 2.0, 1e-12, 0.5, -2.966085951287999 },
            { 2.0, 1e-12, 0.999999999999, -30.914430643547444 },
            { 2.0, 0.01, 0.99, -5.920784877891539 },
            { 2.0, 0.2, 0.8, -1.307505257001309 },
            { 2.0, 0.5, 0.5, 0.41605557909055346 },
            { 2.0, 0.999999999999, 0.5, -26.371428362589487 },
            { 2.0, 0.999999999999, 0.999999999999, 26.59132246705494 },
            { 10.0, 1e-12, 1e-12, 24.665548278435484 },
            { 10.0, 1e-12, 0.5, -32.19396460213969 },
            { 10.0, 1e-12, 0.999999999999, -278.26788463849937 },
            { 10.0, 0.01, 0.99, -54.05263908198108 },
            { 10.0, 0.2, 0.8, -15.673422314090208 },
            { 10.0, 0.5, 0.5, 1.9694693663022758 },
            { 10.0, 0.999999999999, 0.5, -242.7428407756116 },
            { 10.0, 0.999999999999, 0.999999999999, 28.51128817216623 },
            { 100.0, 1e-12, 1e-12, 27.58332425102854 },
            { 100.0, 1e-12, 0.5, -362.6442639484146 },
            { 100.0, 1e-12, 0.999999999999, -3062.525912542782 },
            { 100.0, 0.01, 0.99, -603.4820984619677 },
            { 100.0, 0.2, 0.8, -191.24816860109675 },
            { 100.0, 0.5, 0.5, 4.2776212915497664 },
            { 100.0, 0.999999999999, 0.5, -2694.2198915530175 },
            { 100.0, 0.999999999999, 0.999999999999, 30.84680019871416 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new GumbelCopula(cases[i, 0]);
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
        /// At u = 0 and u = 1 the inverse conditional CDF is its limit: for θ > 1 it tends to 0 as u → 0 and to 1 as u → 1
        /// for 0 &lt; t &lt; 1; at θ = 1 (independence) and for t = 0 or t = 1 it is t.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDFEdges()
        {
            var copula = new GumbelCopula(2d);
            foreach (double t in new[] { 1E-9, 0.3, 0.5, 0.9, 1d - 1E-9 })
            {
                Assert.AreEqual(0d, copula.InverseConditionalCDF(0d, t), 0d);
                Assert.AreEqual(1d, copula.InverseConditionalCDF(1d, t), 0d);
                Assert.AreEqual(t, new GumbelCopula(1d).InverseConditionalCDF(0d, t), 0d);
                Assert.AreEqual(t, new GumbelCopula(1d).InverseConditionalCDF(1d, t), 0d);
            }
            foreach (double u in new[] { 0d, 1d })
            {
                Assert.AreEqual(0d, copula.InverseConditionalCDF(u, 0d), 0d);
                Assert.AreEqual(1d, copula.InverseConditionalCDF(u, 1d), 0d);
            }
            // The limits are approached from inside the square.
            Assert.IsLessThan(1E-6, copula.InverseConditionalCDF(1E-300, 0.5));
            Assert.IsGreaterThan(1d - 1E-6, copula.InverseConditionalCDF(1d - 1E-15, 0.5));
        }

        /// <summary>
        /// Every generator member against mpmath 1.3.0 at 150 digits, for φ(t) = (−ln t)^θ, the independence generator −ln t at θ = 1: Generator, GeneratorInverse (at
        /// the rounded generator value), GeneratorPrime, GeneratorPrime2 and GeneratorPrimeInverse (at the rounded slope),
        /// on t from 1E-12 to 1 - 1E-12 and θ including the independence region, to 1E-13 relative.
        /// </summary>
        [TestMethod]
        public void Test_GeneratorMembersMatchMpmath()
        {
            // theta, t, phi(t), phi^-1(phi), phi'(t), phi''(t), (phi')^-1(phi').
            double[,] cases =
            {
            { 1.0, 1e-12, 27.631021115928547, 1.000000000000001e-12, -1000000000000.0, 1e+24, 1e-12 },
            { 1.0, 0.01, 4.605170185988091, 0.010000000000000004, -100.0, 10000.0, 0.01 },
            { 1.0, 0.5, 0.6931471805599453, 0.5, -2.0, 4.0, 0.5 },
            { 1.0, 0.99, 0.01005033585350145, 0.99, -1.0101010101010102, 1.0203040506070808, 0.99 },
            { 1.0, 0.999999999999, 9.999778782803785e-13, 0.999999999999, -1.000000000001, 1.000000000002, 0.9999999999989999 },
            { 1.000000000001, 1e-12, 27.63102111602026, 1.0000000000000014e-12, -1000000000004.3193, 1.0000000000043556e+24, 1e-12 },
            { 1.000000000001, 0.01, 4.605170185995125, 0.009999999999999997, -100.00000000025274, 10000.000000027445, 0.01 },
            { 1.000000000001, 0.5, 0.6931471805596913, 0.5, -2.000000000001267, 4.000000000008305, 0.5 },
            { 1.000000000001, 0.99, 0.010050335853455213, 0.99, -1.0101010100973733, 1.0203040507049357, 0.99 },
            { 1.000000000001, 0.999999999999, 9.999778782527455e-13, 0.999999999999, -0.9999999999743666, 2.0001110247087515, 0.999999999999 },
            { 2.0, 1e-12, 763.4733279088894, 9.999999999999992e-13, -55262042231857.09, 5.72620422318571e+25, 1e-12 },
            { 2.0, 0.01, 21.20759244191359, 0.010000000000000002, -921.0340371976182, 112103.40371976182, 0.01 },
            { 2.0, 0.5, 0.48045301391820144, 0.5, -2.772588722239781, 13.545177444479563, 0.5 },
            { 2.0, 0.99, 0.00010100925076817673, 0.99, -0.020303708794952424, 2.0611168979767402, 0.99 },
            { 2.0, 0.999999999999, 9.999557570501274e-25, 0.999999999999, -1.9999557565627567e-12, 2.000000000006, 0.999999999999 },
            { 20.0, 1e-12, 6.728816655279284e+28, 1e-12, -4.870479905210814e+40, 8.219582271379219e+52, 1e-12 },
            { 20.0, 0.01, 18404061983603.156, 0.01, -7992782564084265.0, 4.0969385509110313e+18, 0.01 },
            { 20.0, 0.5, 0.0006554042011709004, 0.5, -0.037821935632282026, 2.1491335922581287, 0.5 },
            { 20.0, 0.99, 1.1056342854075848e-40, 0.99, -2.2224178868677817e-37, 4.246129361276269e-34, 0.99 },
            { 20.0, 0.999999999999, 9.995576585756203e-241, 0.999999999999, -1.9991595420001065e-227, 3.798487158872116e-214, 0.999999999999 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new GumbelCopula(cases[i, 0]);
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
        /// ConditionalCDF against its closed form evaluated in mpmath 1.3.0 at 150 digits, for θ = 1, 2 and 20 and u, v from
        /// 1E-12 to 1 - 1E-12, to 1E-13 relative; near u, v = 1 the generic ratio φ′(u)/φ′(C(u,v)) loses up to 1E-5 there.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalCDFMatchesMpmath()
        {
            double[,] cases =
            {
            { 1.0, 1e-12, 0.5, 0.5 },
            { 1.0, 0.01, 0.99, 0.99 },
            { 1.0, 0.5, 0.5, 0.5 },
            { 1.0, 0.99, 0.01, 0.01 },
            { 1.0, 0.999999, 0.999999, 0.999999 },
            { 1.0, 0.999999999999, 0.3, 0.3 },
            { 2.0, 1e-12, 0.5, 0.9910331735604143 },
            { 2.0, 0.01, 0.99, 0.9999866517269862 },
            { 2.0, 0.5, 0.5, 0.530633048967315 },
            { 2.0, 0.99, 0.01, 2.2044177473918812e-05 },
            { 2.0, 0.999999, 0.999999, 0.707106488293243 },
            { 2.0, 0.999999999999, 0.3, 2.491695513441184e-13 },
            { 20.0, 1e-12, 0.5, 1.0 },
            { 20.0, 0.01, 0.99, 1.0 },
            { 20.0, 0.5, 0.5, 0.5051329565483142 },
            { 20.0, 0.99, 0.01, 2.780530896529405e-53 },
            { 20.0, 0.999999, 0.999999, 0.5176324436664106 },
            { 20.0, 0.999999999999, 0.3, 8.814990507136217e-231 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new GumbelCopula(cases[i, 0]);
                Assert.AreEqual(cases[i, 3], copula.ConditionalCDF(cases[i, 1], cases[i, 2]), 1E-13 * cases[i, 3], $"case {i}");
            }
        }

        /// <summary>
        /// InverseConditionalCDF against the root of the closed-form conditional found by mpmath 1.3.0 at 60 digits
        /// (bisection in the logit of v to 1E-40), for θ = 1.5, 5 and 20, u = 1E-12, 0.3 and 0.99, and t = 1E-12, 0.5 and
        /// 1 - 1E-6, to 1E-14 relative. Over the wider grid θ ∈ {1, 1.5, 2, 5, 20}, u from 1E-12 to 1 - 1E-6 and t from
        /// 1E-12 to 1 - 1E-12 the largest error is under 7E-15.
        /// </summary>
        [TestMethod]
        public void Test_InverseConditionalCDFMatchesMpmath()
        {
            double[,] cases =
            {
            { 1.5, 1e-12, 1e-12, 1.6961550041269132e-18 },
            { 1.5, 1e-12, 0.5, 0.046013397739818315 },
            { 1.5, 1e-12, 0.999999, 0.9996086344514264 },
            { 1.5, 0.3, 1e-12, 1.6906187490020034e-12 },
            { 1.5, 0.3, 0.5, 0.3849293386077215 },
            { 1.5, 0.3, 0.999999, 0.999889419484749 },
            { 1.5, 0.99, 1e-12, 4.813944486883003e-11 },
            { 1.5, 0.99, 0.5, 0.9661015204764247 },
            { 1.5, 0.99, 0.999999, 0.9999979369986796 },
            { 5.0, 1e-12, 1e-12, 2.0211889901654108e-23 },
            { 5.0, 1e-12, 0.5, 1.6548943148881274e-08 },
            { 5.0, 1e-12, 0.999999, 0.29954125343061105 },
            { 5.0, 0.3, 1e-12, 1.5031792744697492e-08 },
            { 5.0, 0.3, 0.5, 0.30554100079965096 },
            { 5.0, 0.3, 0.999999, 0.9274090972900307 },
            { 5.0, 0.99, 1e-12, 0.021298717949478126 },
            { 5.0, 0.99, 0.5, 0.9893493053257416 },
            { 5.0, 0.99, 0.999999, 0.9993374772379477 },
            { 20.0, 1e-12, 1e-12, 1.419598714742475e-20 },
            { 20.0, 1e-12, 0.5, 4.195089815401857e-12 },
            { 20.0, 1e-12, 0.999999, 1.7179239668158156e-06 },
            { 20.0, 0.3, 1e-12, 0.012719034125384528 },
            { 20.0, 0.3, 0.5, 0.30027942975471444 },
            { 20.0, 0.3, 0.999999, 0.5471076122358633 },
            { 20.0, 0.99, 1e-12, 0.9579550947835443 },
            { 20.0, 0.99, 0.5, 0.9899646558980333 },
            { 20.0, 0.99, 0.999999, 0.9949628275101948 },
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                var copula = new GumbelCopula(cases[i, 0]);
                double v = copula.InverseConditionalCDF(cases[i, 1], cases[i, 2]);
                Assert.AreEqual(cases[i, 3], v, 1E-14 * cases[i, 3], $"case {i}");
                // The root reproduces t through the conditional CDF, to the change a few ulps of v make in h, plus the
                // rounding of h's exponent, about |ln t| ulps of t.
                double slack = 8.9E-16 * (v * copula.PDF(cases[i, 1], v) + cases[i, 2] * (1d - Math.Log(cases[i, 2])));
                Assert.AreEqual(cases[i, 2], copula.ConditionalCDF(cases[i, 1], v), slack, $"case {i}: round trip");
            }
        }
    }
}
