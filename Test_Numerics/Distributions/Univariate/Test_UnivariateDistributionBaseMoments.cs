using System;
using System.Linq;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling;

namespace Distributions.Univariate
{
    /// <summary>
    /// Tests for the fixed-step bin-expectation moments inherited from <see cref="UnivariateDistributionBase"/>:
    /// <see cref="UnivariateDistributionBase.CentralMoments(int)"/> and
    /// <see cref="UnivariateDistributionBase.ConditionalMoments(double, double)"/>.
    /// </summary>
    [TestClass]
    public class Test_UnivariateDistributionBaseMoments
    {
        /// <summary>Distributions covering unbounded, bounded, skewed, log-space, and table-backed families.</summary>
        /// <returns>The distributions evaluated in every case.</returns>
        private static UnivariateDistributionBase[] Distributions() =>
        [
            new Normal(3d, 2d),
            new Weibull(5d, 1.4d),
            new GeneralizedExtremeValue(8d, 2d, 0.1d),
            new LogPearsonTypeIII(1d, 0.2d, 0.1d),
            new GammaDistribution(2d, 3d),
            new Uniform(-0d, 4d),
            new EmpiricalDistribution(new[] { 1d, 2d, 4d, 8d }, new[] { 0.05d, 0.4d, 0.8d, 0.99d }),
        ];

        /// <summary>
        /// The stratified bins are contiguous bit for bit, so each interior bin's lower CDF is the previous bin's
        /// upper CDF.
        /// </summary>
        [TestMethod]
        public void StratifiedBinsAreBitwiseContiguous()
        {
            foreach (var (a, b) in new[] { (-3d, 9d), (0.1d, 0.7d), (-1E-8, 1E-8), (1E3, 1E3 + 1d) })
            {
                foreach (int steps in new[] { 2, 7, 300 })
                {
                    var bins = Stratify.XValues(new StratificationOptions(a, b, steps));
                    for (int i = 1; i < bins.Count; i++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(bins[i - 1].UpperBound), BitConverter.DoubleToInt64Bits(bins[i].LowerBound));
                }
            }
        }

        /// <summary>
        /// The fixed-step central moments equal the expression that evaluates both bin edges of every bin, bit for bit,
        /// for several step counts. The standard Normal(3, 2) case also matches its closed-form mean and standard
        /// deviation to the accuracy of the 300-bin expectation.
        /// </summary>
        [TestMethod]
        public void CentralMoments_MatchesTwoEdgeExpression()
        {
            foreach (var distribution in Distributions())
            {
                foreach (int steps in new[] { 2, 3, 7, 300 })
                {
                    var expected = TwoEdgeCentralMoments(distribution, steps);
                    var actual = distribution.CentralMoments(steps);
                    for (int k = 0; k < 4; k++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[k]), BitConverter.DoubleToInt64Bits(actual[k]), $"{distribution.DisplayName} steps={steps} moment {k}");
                }
            }
            var normal = new Normal(3d, 2d).CentralMoments(300);
            Assert.AreEqual(3d, normal[0], 1E-9);
            Assert.AreEqual(2d, normal[1], 1E-3);
        }

        /// <summary>
        /// The conditional moments equal the expression that evaluates both bin edges of every bin, bit for bit,
        /// for the families that inherit the stratified expectation.
        /// </summary>
        [TestMethod]
        public void ConditionalMoments_MatchesTwoEdgeExpression()
        {
            // Normal and Gamma override this member with closed forms.
            foreach (var distribution in Distributions().Where(d => d is not Normal && d is not GammaDistribution).Append(new Logistic(2d, 0.5d)))
            {
                foreach (var (pl, pu) in new[] { (0.2d, 0.7d), (0.01d, 0.99d), (0.5d, 0.5000001d) })
                {
                    double a = distribution.InverseCDF(pl), b = distribution.InverseCDF(pu);
                    var expected = TwoEdgeConditionalMoments(distribution, a, b);
                    var actual = distribution.ConditionalMoments(a, b);
                    for (int k = 0; k < 4; k++)
                        Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[k]), BitConverter.DoubleToInt64Bits(actual[k]), $"{distribution.DisplayName} [{pl},{pu}] moment {k}");
                }
            }
        }

        /// <summary>The fixed-step central moments evaluating the CDF at both edges of every interior bin.</summary>
        /// <param name="d">The distribution.</param>
        /// <param name="steps">The number of bins.</param>
        /// <returns>Mean, standard deviation, skewness, and kurtosis.</returns>
        private static double[] TwoEdgeCentralMoments(UnivariateDistributionBase d, int steps)
        {
            double a = d.InverseCDF(1E-8);
            double b = d.InverseCDF(1 - 1E-8);
            if (a >= b) return [a, double.NaN, double.NaN, double.NaN];
            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));
            var dFx = new double[steps];
            var coordinates = new double[steps];
            double reference = a / 2 + b / 2;
            double scale = Math.Max(Math.Abs(a - reference), Math.Abs(b - reference));
            dFx[0] = d.CDF(bins[0].UpperBound);
            coordinates[0] = DistributionNumerics.Standardize(bins[0].UpperBound, reference, scale);
            for (int i = 1; i < steps - 1; i++)
            {
                dFx[i] = d.CDF(bins[i].UpperBound) - d.CDF(bins[i].LowerBound);
                coordinates[i] = DistributionNumerics.Standardize(bins[i].Midpoint, reference, scale);
            }
            dFx[steps - 1] = 1 - d.CDF(bins.Last().LowerBound);
            coordinates[steps - 1] = DistributionNumerics.Standardize(bins.Last().LowerBound, reference, scale);
            double offset = 0;
            for (int i = 0; i < steps; i++) offset += coordinates[i] * dFx[i];
            double second = 0, third = 0, fourth = 0;
            for (int i = 0; i < steps; i++)
            {
                double centered = coordinates[i] - offset, square = centered * centered;
                second += square * dFx[i];
                third += square * centered * dFx[i];
                fourth += square * square * dFx[i];
            }
            return [reference + scale * offset, scale * Math.Sqrt(second), third / second / Math.Sqrt(second), fourth / second / second];
        }

        /// <summary>The conditional moments evaluating the CDF at both edges of every bin.</summary>
        /// <param name="d">The distribution.</param>
        /// <param name="a">The lower limit.</param>
        /// <param name="b">The upper limit.</param>
        /// <returns>The conditional mean and the second, third, and fourth conditional central moments about the mean.</returns>
        private static double[] TwoEdgeConditionalMoments(UnivariateDistributionBase d, double a, double b)
        {
            if (a >= b) return [double.NaN, double.NaN, double.NaN, double.NaN];
            int steps = 300;
            double l = d.CDF(a), u = d.CDF(b), totalProb = u - l;
            if (totalProb <= 0.0 || double.IsNaN(totalProb)) return [double.NaN, double.NaN, double.NaN, double.NaN];
            double mean = d.Mean;
            var bins = Stratify.XValues(new StratificationOptions(a, b, steps));
            double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;
            for (int i = 0; i < steps; i++)
            {
                double dF = d.CDF(bins[i].UpperBound) - d.CDF(bins[i].LowerBound);
                double xMid = bins[i].Midpoint;
                double delta = xMid - mean;
                double delta2 = delta * delta;
                sum1 += xMid * dF;
                sum2 += delta2 * dF;
                sum3 += delta2 * delta * dF;
                sum4 += delta2 * delta2 * dF;
            }
            return [sum1 / totalProb, sum2 / totalProb, sum3 / totalProb, sum4 / totalProb];
        }
    }
}
