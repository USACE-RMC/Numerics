using Numerics.Distributions;

namespace Test_Numerics.Distributions
{
    /// <summary>Preserves usable constraint envelopes for nonpositive observations from the pre-hardening baseline.</summary>
    [TestClass]
    public class Test_NonpositiveConstraintRegressions
    {
        /// <summary>Checks literal d80bfa8 initialization and rounded bounds without changing estimation data.</summary>
        [TestMethod]
        [DataRow("GammaDistribution", 0d)]
        [DataRow("GammaDistribution", -1d)]
        [DataRow("LnNormal", 0d)]
        [DataRow("LnNormal", -1d)]
        [DataRow("LogNormal", 0d)]
        [DataRow("LogNormal", -1d)]
        [DataRow("LogPearsonTypeIII", 0d)]
        [DataRow("LogPearsonTypeIII", -1d)]
        [DataRow("Weibull", 0d)]
        [DataRow("Weibull", -1d)]
        public void ConstraintsPreserveUsableLegacyNonpositiveSamples(string family, double first)
        {
            var distribution = (IMaximumLikelihoodEstimation)System.Activator.CreateInstance(
                typeof(Normal).Assembly.GetType("Numerics.Distributions." + family))!;
            double[] sample = { first, 1d, 10d, 100d };
            var actual = distribution.GetParameterConstraints(sample);
            double epsilon = Numerics.Tools.DoubleMachineEpsilon;
            double[] initial;
            double[] lower;
            double[] upper;
            switch (family)
            {
                case "GammaDistribution":
                    initial = first == 0d ? new[] { 84.33333333333334, 0.3290513833992095 }
                        : new[] { 85.78181818181818, 0.3205807545570157 };
                    lower = new[] { epsilon, epsilon };
                    upper = new[] { 1000d, 10d };
                    break;
                case "LnNormal":
                    initial = first == 0d ? new[] { 27.75, 48.376130477746976 }
                        : new[] { 27.5, 48.569537778323564 };
                    lower = new[] { epsilon, epsilon };
                    upper = new[] { 1000d, 1000d };
                    break;
                case "LogNormal":
                    initial = new[] { 0.5000000000000002, 1.2909944487358054 };
                    lower = new[] { -10d, epsilon };
                    upper = new[] { 10d, 3d };
                    break;
                case "LogPearsonTypeIII":
                    initial = new[] { 0.25, 1.707825127659933, -0.7528371991317255 };
                    lower = new[] { -10d, epsilon, -6d };
                    upper = new[] { 10d, 3d, 6d };
                    break;
                default:
                    initial = new[] { 12.336441557126482, 0.493577181580963 };
                    lower = new[] { epsilon, epsilon };
                    upper = new[] { 1000d, 10d };
                    break;
            }
            CollectionAssert.AreEqual(initial, actual.Item1, "Legacy initialization");
            CollectionAssert.AreEqual(lower, actual.Item2, "Legacy lower bounds");
            CollectionAssert.AreEqual(upper, actual.Item3, "Legacy upper bounds");
            CollectionAssert.AreEqual(new[] { first, 1d, 10d, 100d }, sample, "Observations are never modified");
        }
    }
}
