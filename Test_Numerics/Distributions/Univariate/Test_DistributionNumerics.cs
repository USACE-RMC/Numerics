using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>Independent numerical contracts for shared distribution primitives.</summary>
    [TestClass]
    public class Test_DistributionNumerics
    {
        /// <summary>Subnormal derivative scaling must precede the final exponential rounding.</summary>
        [TestMethod]
        public void GammaSubnormalDerivativeAndTinyShapeTailsRetainScale()
        {
            Assert.AreEqual(-5.597735947761608E-21, DistributionNumerics.GammaLogCDF(1E-20, .5), 1E-34);
            Assert.AreEqual(-5.597735947761608E-101, DistributionNumerics.GammaLogCDF(1E-100, .5), 1E-113);
            Assert.AreEqual(4.090863547565521E-321, DistributionNumerics.GammaQuantileShapeDerivative(.9, double.Epsilon), 2 * double.Epsilon);
        }

        /// <summary>R qgamma Richardson derivatives cross-checked independently with pgamma at fixed x.</summary>
        [TestMethod]
        public void GammaShapeDerivativesMatchIndependentQuantileOracles()
        {
            int count = 0;
            foreach (var row in DistributionOracle.Read("extreme-positive.csv"))
            {
                if (row["family"] != "Gamma" || row["quantity"] != "GradientShape" || row["status"] != "finite") continue;
                double a = DistributionOracle.Number(row["shape"]), p = DistributionOracle.Number(row["p"]);
                double q = DistributionNumerics.GammaInverseCDF(a, p);
                double actual = DistributionNumerics.GammaQuantileShapeDerivative(a, q);
                double expected = DistributionOracle.Number(row["value"]);
                Assert.AreEqual(expected, actual, Math.Abs(expected) * 5E-8, $"shape={a},p={p}");
                count++;
            }
            Assert.AreEqual(49, count);
        }

        /// <summary>R 4.4.3 pgamma(log.p=TRUE), including its direct complemented tail.</summary>
        [TestMethod]
        public void GammaLogTailsMatchROracles()
        {
            double[,] cases = {
                { .001, 1E-10, -.022449457331756989, -3.8076925670559727 },
                { .001, 1, -.0002196324750319099, -8.4236647908310722 },
                { .001, 1000, 0, -1013.8090239153952 },
                { .5, 1E-10, -11.392143227368317, -1.1283855333035131E-5 },
                { .5, 1, -.17114331524104021, -1.8496055099332522 },
                { .5, 1000, 0, -1004.0267419589519 },
                { 1, 1E-10, -23.025850929990458, -1E-10 },
                { 1, 1, -.45867514538708193, -1 },
                { 100, 80, -4.0681907911970132, -.017256351106458189 },
                { 100, 100, -.66689715058528931, -.72010489302547409 },
                { 100, 120, -.028259299048920376, -3.5804290809275314 },
                { 10000, 9800, -3.8073232363866576, -.022457843956297584 },
                { 10000, 10000, -.69049109440197243, -.69581034030382005 },
                { 10000, 10200, -.023562756305381644, -3.7598461809322816 },
                { 1E8, 99980000, -3.7834216952685238, -.023007384281740352 },
                { 1E8, 1E8, -.69312058476158833, -.69317377706565741 },
                { 1E8, 100020000, -.023018433854231277, -3.7829470521492174 }
            };
            for (int i = 0; i < cases.GetLength(0); i++)
            {
                Assert.AreEqual(cases[i, 2], DistributionNumerics.GammaLogCDF(cases[i, 0], cases[i, 1]), 3E-12, $"lower case {i}");
                Assert.AreEqual(cases[i, 3], DistributionNumerics.GammaLogSurvival(cases[i, 0], cases[i, 1]), 3E-12, $"upper case {i}");
            }
        }

        /// <summary>Independent closed forms and R references, preserving tiny input probabilities.</summary>
        [TestMethod]
        public void GammaInverseUsesBothTailsDirectly()
        {
            Assert.AreEqual(46.051701859880914, DistributionNumerics.GammaInverseCDF(1, 1E-20, true), 1E-13);
            Assert.AreEqual(1E-20, DistributionNumerics.GammaInverseCDF(1, 1E-20), 1E-35);
            double q = DistributionNumerics.GammaInverseCDF(.001, .5);
            Assert.AreEqual(5.244206408274966E-302, q, 2E-313);
            foreach (double a in new[] { .1, .5, 2, 100, 10000, 1E8 })
            foreach (double p in new[] { 1E-10, .1, .5 })
            foreach (bool upper in new[] { false, true })
            {
                double value = DistributionNumerics.GammaInverseCDF(a, p, upper);
                double actual = upper ? DistributionNumerics.GammaLogSurvival(a, value) : DistributionNumerics.GammaLogCDF(a, value);
                Assert.AreEqual(Math.Log(p), actual, a >= 1E8 ? 2E-11 : 2E-12, $"a={a},p={p},upper={upper}");
            }
        }

        /// <summary>R pnorm log-tail references and exact exponential divided-difference limits.</summary>
        [TestMethod]
        public void LogProbabilityPrimitivesPreserveLimits()
        {
            Assert.AreEqual(-804.6084420137538, DistributionNumerics.NormalLogCDF(-40), 2E-13);
            Assert.AreEqual(-43.62814911333212, DistributionNumerics.NormalLogSurvival(9), 2E-13);
            Assert.AreEqual(double.NegativeInfinity, DistributionNumerics.Log1mExp(0));
            Assert.AreEqual(0, DistributionNumerics.Log1mExp(double.NegativeInfinity));
            Assert.AreEqual(double.NegativeInfinity, DistributionNumerics.LogDifference(-1, -1));
            Assert.AreEqual(1, DistributionNumerics.Exprel(0));
            Assert.AreEqual(.5, DistributionNumerics.ExprelDerivative(0));
            Assert.AreEqual(2, DistributionNumerics.Standardize(1E308, -1E308, 1E308));
        }

        /// <summary>Scaled pivoting preserves sign, exact duplicate-row singularity, and finite logs.</summary>
        [TestMethod]
        public void LogDeterminantAvoidsProductOverflowAndArtificialPivots()
        {
            double log = DistributionNumerics.LogAbsDeterminant(new double[,] { { 1E200, 1E200 }, { 1E200, -1E200 } }, out int sign);
            Assert.AreEqual(400 * Math.Log(10) + Math.Log(2), log, 2E-13);
            Assert.AreEqual(-1, sign);
            Assert.AreEqual(double.NegativeInfinity, DistributionNumerics.LogAbsDeterminant(new double[,] { { 1, 2 }, { 1, 2 } }, out sign));
            Assert.AreEqual(0, sign);
            Assert.AreEqual(0, DistributionNumerics.LogAbsDeterminant(new double[,] { { 1E308, 1E-308 }, { 1E308, 2E-308 } }, out sign), 1E-12);
            Assert.AreEqual(1, sign);
            Assert.AreEqual(double.NegativeInfinity, DistributionNumerics.LogAbsDeterminant(new double[,] { { 1, 2, 3 }, { 4, 5, 6 }, { 5, 7, 9 } }, out sign));
            Assert.AreEqual(0, sign);
            Assert.AreEqual(double.NegativeInfinity, new KappaFour().LogAbsQuantileJacobian(new[] { .5, .5, .5, .5 }));
            Assert.AreEqual(Math.Log(2 * Normal.StandardZ(.9)), new Normal().LogAbsQuantileJacobian(new[] { .1, .9 }), 1E-14);
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new Normal().LogAbsQuantileJacobian(new[] { 0d, .9 }));
        }
    }
}
