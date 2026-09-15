using Numerics.Distributions;

namespace Distributions
{
    /// <summary>Pins approved prior envelopes and initial values for the legacy constraint fixtures.</summary>
    [TestClass]
    public class Test_LegacyParameterConstraints
    {
        // Original fixtures captured from d80bfa8621c48a78cf4ebad7f01326841fac37aa on net8.0.
        // Gamma initializers and LogNormal/LogPearsonTypeIII location bounds repinned to 202095a.
        // Exponential location upper bounds repinned to 2bba500.
        // These fixtures freeze valid initialization and family-specific prior envelopes.
        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_ordinary_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 9d, 4d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_skewed_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { -0.9d, 11.4d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1.9d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1.0000000000000002d, 1000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_nearUnity_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0005d, 0.002000000000000076d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -8.9995d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1.0010000000000001d, 0.1d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_negative_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -17d, 4d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -117d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { -16d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_zeroMean_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { -4d, 4d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -14d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { -3d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Exponential_subUnity_PreservesBaselineConstraints()
        {
            var result = new Exponential().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.00833333333333334d, 0.3666666666666667d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -0.00166666666666666d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 0.10000000000000012d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GammaDistribution_ordinary_PreservesBaselineConstraints()
        {
            var result = new GammaDistribution().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 0.5128205128205128d, 25.349999999999998d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 1000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GammaDistribution_skewed_PreservesBaselineConstraints()
        {
            var result = new GammaDistribution().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 13.4d, 0.7835820895522387d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GammaDistribution_nearUnity_PreservesBaselineConstraints()
        {
            var result = new GammaDistribution().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.6625103906900188e-06d, 603003.7499999721d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 0.0001d, 10000000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GammaDistribution_subUnity_PreservesBaselineConstraints()
        {
            var result = new GammaDistribution().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.25555555555555554d, 1.467391304347826d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_ordinary_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 11.964905241795506d, 2.944098222941484d, 0.2846308210256667d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_skewed_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 3.5251212632019016d, 5.1182032350575986d, -0.44798569084540746d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_nearUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0019824526208967d, 0.0014720491114694771d, 0.28463082102423143d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_negative_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -14.035094758204492d, 2.9440982229414914d, 0.28463082102566944d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_zeroMean_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { -1.035094758204496d, 2.9440982229414825d, 0.284630821025665d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedExtremeValue_subUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedExtremeValue().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.1832069588359371d, 0.18598234830157992d, -0.31866881623469623d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_ordinary_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13.000000000000002d, 1.666666666666666d, 8.881784197001252e-16d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_skewed_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 5.685723293384427d, 4.328789084801283d, -0.49253731343283613d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_nearUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0024999999999988d, 0.0008333333333332416d, -7.993605777301127e-13d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.01d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_negative_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -12.999999999999995d, 1.6666666666666679d, 2.220446049250313e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_zeroMean_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { 1.11022302462516e-16d, 1.6666666666666667d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1e-14d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1e-14d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedLogistic_subUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedLogistic().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.26062346440010725d, 0.1469106309971434d, -0.3913043478260869d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedNormal_ordinary_PreservesBaselineConstraints()
        {
            var result = new GeneralizedNormal().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13d, 3.03207911111111d, 1.817793382485888e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedNormal_skewed_PreservesBaselineConstraints()
        {
            var result = new GeneralizedNormal().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 5.1845392376947785d, 7.317708020900799d, -1.0739272660014685d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedNormal_nearUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedNormal().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0025d, 0.001477372521655275d, -1.6360140442372994e-12d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedNormal_negative_PreservesBaselineConstraints()
        {
            var result = new GeneralizedNormal().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -13d, 3.1009900000000026d, 4.5444834562147206e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedNormal_subUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedNormal().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.24860780884509534d, 0.25425864998736986d, -0.8320397358919761d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_ordinary_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 7.999999999999996d, 10.000000000000027d, 1.0000000000000036d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -2.0000000000000044d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 1000d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_skewed_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { -0.7559999999999949d, 7.654079999999991d, -0.32000000000000056d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1.755999999999995d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1.0000000000000002d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_nearUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0000000000000029d, 0.004999999999986127d, 0.9999999999968026d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -8.999999999999996d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1.0010000000000001d, 0.1d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_negative_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -18.000000000000018d, 10.000000000000082d, 1.0000000000000089d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -118.00000000000001d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { -16d, 1000d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_zeroMean_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { -5d, 10d, 1d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -15d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { -3d, 100d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void GeneralizedPareto_subUnity_PreservesBaselineConstraints()
        {
            var result = new GeneralizedPareto().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.015625d, 0.314453125d, -0.12499999999999992d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -0.084375d, 1.11022302462516e-16d, -10d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 0.10000000000000012d, 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_ordinary_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 11.612089704538555d, 2.4044917348149384d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_skewed_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 4.92060061224499d, 9.666056773956054d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_nearUnity_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0018060448522692d, 0.0012022458674073372d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_negative_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -14.387910295461445d, 2.4044917348149406d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_zeroMean_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { -1.3879102954614455d, 2.4044917348149393d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Gumbel_subUnity_PreservesBaselineConstraints()
        {
            var result = new Gumbel().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.21539031602193381d, 0.27651654950371796d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_ordinary_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 8.000006666483934d, 9.999987917531367d, 0.9999994001012955d, 0.9999991500014682d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_skewed_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { -41.20730147684273d, 41.160111238748186d, 0.2828390688884041d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 1000d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_nearUnity_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0000000033331053d, 0.004999993959029831d, 0.9999994001312362d, 0.9999991500350589d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_negative_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -17.99999333350232d, 9.999987917504608d, 0.9999994000997355d, 0.9999991499997695d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_zeroMean_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { -4.999993333505345d, 9.999987917511882d, 0.9999994001003065d, 0.9999991500001029d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void KappaFour_subUnity_PreservesBaselineConstraints()
        {
            var result = new KappaFour().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -0.739635718854716d, 1.0608365097007928d, 0.36224410221013503d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d, -10d, -2d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 100d, 10d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LnNormal_ordinary_PreservesBaselineConstraints()
        {
            var result = new LnNormal().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13d, 2.581988897471611d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LnNormal_skewed_PreservesBaselineConstraints()
        {
            var result = new LnNormal().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 10.5d, 11.861703081766969d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 1000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LnNormal_nearUnity_PreservesBaselineConstraints()
        {
            var result = new LnNormal().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0025d, 0.0012909944487358356d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LnNormal_subUnity_PreservesBaselineConstraints()
        {
            var result = new LnNormal().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.375d, 0.30956959368344517d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Logistic_ordinary_PreservesBaselineConstraints()
        {
            var result = new Logistic().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13d, 1.4235250868343539d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Logistic_skewed_PreservesBaselineConstraints()
        {
            var result = new Logistic().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 10.5d, 6.539699657891849d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Logistic_nearUnity_PreservesBaselineConstraints()
        {
            var result = new Logistic().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0025d, 0.0007117625434171935d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.01d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Logistic_subUnity_PreservesBaselineConstraints()
        {
            var result = new Logistic().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.375d, 0.17067466214166682d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogNormal_ordinary_PreservesBaselineConstraints()
        {
            var result = new LogNormal().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 1.1073573160954466d, 0.08791220351278327d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 0d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 3d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogNormal_skewed_PreservesBaselineConstraints()
        {
            var result = new LogNormal().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 0.7525749891599528d, 0.5631755534583315d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 2d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogNormal_nearUnity_PreservesBaselineConstraints()
        {
            var result = new LogNormal().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.0010841112099910272d, 0.0005592740172164043d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 2d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogNormal_subUnity_PreservesBaselineConstraints()
        {
            var result = new LogNormal().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -0.5484550065040281d, 0.38862805330516337d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -2d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1d, 2d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_ordinary_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 1.1073573160954466d, 0.08791220351278327d, -0.2899042849970034d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 0d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 3d, 2d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_skewed_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 0.7525749891599528d, 0.5631755534583315d, -1.1187971499007316e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 2d, 2d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_nearUnity_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.0010841112099910272d, 0.0005592740172164043d, -0.0018543970617275146d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 2d, 2d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_subUnity_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -0.5484550065040281d, 0.38862805330516337d, -1.2610041890269048e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -2d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1d, 2d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_ordinary_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13d, 2.581988897471611d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_skewed_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 10.5d, 11.861703081766969d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 1000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_nearUnity_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0025d, 0.0012909944487358356d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_negative_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { -16d, -14d, -12d, -10d });
            AssertArray(new double[] { -13d, 2.581988897471611d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_zeroMean_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            AssertArray(new double[] { 1.11022302462516e-16d, 2.581988897471611d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1e-14d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1e-14d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Normal_subUnity_PreservesBaselineConstraints()
        {
            var result = new Normal().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.375d, 0.30956959368344517d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void PearsonTypeIII_ordinary_PreservesBaselineConstraints()
        {
            var result = new PearsonTypeIII().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
            AssertArray(new double[] { 13d, 2.581988897471611d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void PearsonTypeIII_skewed_PreservesBaselineConstraints()
        {
            var result = new PearsonTypeIII().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 10.5d, 11.861703081766969d, 1.4996929578604363d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -1000d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 1000d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void PearsonTypeIII_nearUnity_PreservesBaselineConstraints()
        {
            var result = new PearsonTypeIII().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0025d, 0.0012909944487358356d, -1.025170155751457e-15d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -100d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 0.1d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void PearsonTypeIII_subUnity_PreservesBaselineConstraints()
        {
            var result = new PearsonTypeIII().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.375d, 0.30956959368344517d, 1.1376243669576893d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -10d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 10d, 6d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Weibull_ordinary_PreservesBaselineConstraints()
        {
            var result = new Weibull().GetParameterConstraints(new double[] { 10d, 12d, 14d, 16d });
#if NETFRAMEWORK
            // Captured independently from d80bfa8 on CLR 4.0.30319.42000; retain exact equality.
            AssertArray(new double[] { 13.949101633065432d, 6.6844918509346618d }, result.Item1, "initial");
#else
            AssertArray(new double[] { 13.949101633065432d, 6.684491850934669d }, result.Item1, "initial");
#endif
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 100d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Weibull_skewed_PreservesBaselineConstraints()
        {
            var result = new Weibull().GetParameterConstraints(new double[] { 1d, 2d, 4d, 8d, 16d, 32d });
            AssertArray(new double[] { 10.19351618226867d, 0.9402756427497648d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 1000d, 10d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Weibull_nearUnity_PreservesBaselineConstraints()
        {
            var result = new Weibull().GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 1.0030553641004294d, 1007.3379396499552d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 100d, 100000d }, result.Item3, "upper");
        }

        /// <summary>Checks initialization and bounds against the literal baseline fixture.</summary>
        [TestMethod]
        public void Weibull_subUnity_PreservesBaselineConstraints()
        {
            var result = new Weibull().GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { 0.4158025497714148d, 1.4492910264321852d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { 1.11022302462516e-16d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 10d, 100d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogNormal_nearUnity_BaseTwo_PreservesBaselineConstraints()
        {
            var result = new LogNormal { Base = 2d }.GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.003601339486451527d, 0.0018578680705316926d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -4d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 7d, 7d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogNormal_subUnity_BaseTwo_PreservesBaselineConstraints()
        {
            var result = new LogNormal { Base = 2d }.GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -1.8219280948873622d, 1.2909944487358056d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -7d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 4d, 7d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogNormal_nearUnity_NaturalBase_PreservesBaselineConstraints()
        {
            var result = new LogNormal { Base = 2.718281828459045d }.GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.002496258311273077d, 0.0012877760149413882d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -3d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 5d, 5d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogNormal_subUnity_NaturalBase_PreservesBaselineConstraints()
        {
            var result = new LogNormal { Base = 2.718281828459045d }.GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -1.2628643221541278d, 0.8948491622597645d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -5d, 1.11022302462516e-16d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 3d, 5d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_nearUnity_BaseTwo_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII { Base = 2d }.GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.003601339486451527d, 0.0018578680705316926d, -0.0018543970617339045d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -4d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 7d, 7d, 6d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_subUnity_BaseTwo_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII { Base = 2d }.GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -1.8219280948873622d, 1.2909944487358056d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -7d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 4d, 7d, 6d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_nearUnity_NaturalBase_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII { Base = 2.718281828459045d }.GetParameterConstraints(new double[] { 1.001d, 1.002d, 1.003d, 1.004d });
            AssertArray(new double[] { 0.002496258311273077d, 0.0012877760149413882d, -0.001854397061727415d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -3d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 5d, 5d, 6d }, result.Item3, "upper");
        }

        /// <summary>Preserves physical-scale decade rounding before conversion to the configured log base.</summary>
        [TestMethod]
        public void LogPearsonTypeIII_subUnity_NaturalBase_PreservesBaselineConstraints()
        {
            var result = new LogPearsonTypeIII { Base = 2.718281828459045d }.GetParameterConstraints(new double[] { 0.1d, 0.2d, 0.4d, 0.8d });
            AssertArray(new double[] { -1.2628643221541278d, 0.8948491622597645d, 0d }, result.Item1, "initial");
            CollectionAssert.AreEqual(new double[] { -5d, 1.11022302462516e-16d, -6d }, result.Item2, "lower");
            CollectionAssert.AreEqual(new double[] { 3d, 5d, 6d }, result.Item3, "upper");
        }

        /// <summary>Retains hardened zero-center bounds when this family's old location envelope was invalid.</summary>
        [TestMethod]
        public void Logistic_ZeroMean_RetainsExceptionalFallback()
        {
            var result = new Logistic().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            Assert.AreEqual(0d, result.Item1[0]);
            Assert.AreEqual(-100d, result.Item2[0]);
            Assert.AreEqual(100d, result.Item3[0]);
        }

        /// <summary>Retains hardened zero-center bounds when this family's old location envelope was invalid.</summary>
        [TestMethod]
        public void PearsonTypeIII_ZeroMean_RetainsExceptionalFallback()
        {
            var result = new PearsonTypeIII().GetParameterConstraints(new double[] { -3d, -1d, 1d, 3d });
            Assert.AreEqual(0d, result.Item1[0]);
            Assert.AreEqual(-100d, result.Item2[0]);
            Assert.AreEqual(100d, result.Item3[0]);
        }

        /// <summary>Preserves a finite ordered legacy envelope whose initializer equals the scale lower bound.</summary>
        [TestMethod]
        public void Normal_ScaleAtLegacyLowerBound_PreservesBoundaryInitializer()
        {
            double epsilon = Numerics.Tools.DoubleMachineEpsilon;
            var result = new Normal().GetParameterConstraints(new double[] { 0d, 0d, 0d, 2d * epsilon });
            Assert.AreEqual(epsilon, result.Item1[1]);
            Assert.AreEqual(epsilon, result.Item2[1]);
            Assert.AreEqual(1E-14, result.Item3[1]);
        }

        /// <summary>Preserves the fallback exception identity and sample parameter expected by fitting callers.</summary>
        [TestMethod]
        public void FailedInitializers_PreserveFallbackSampleException()
        {
            var legacyFailure = new ArithmeticException("Legacy moments overflowed.");
            var fallbackFailure = new ArgumentException("No finite bounds exist.", "sample");
            var actual = Assert.ThrowsExactly<ArgumentException>(() => DistributionNumerics.PreferLegacyConstraints(
                () => throw legacyFailure, () => throw fallbackFailure));
            Assert.AreSame(fallbackFailure, actual);
            Assert.AreEqual("sample", actual.ParamName);
            Assert.AreSame(legacyFailure, actual.Data["LegacyParameterConstraintsFailure"]);
        }

        private static void AssertArray(double[] expected, double[] actual, string label)
        {
            Assert.HasCount(expected.Length, actual, label);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], $"{label}[{i}]");
        }
    }
}
