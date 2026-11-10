using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.MachineLearning;
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MachineLearning
{
    /// <summary>
    /// Unit tests for Gaussian process regression.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    /// </list>
    /// </para>
    /// <para>
    /// Oracle values come from docs/machine-learning/oracles: gpr_sklearn_oracle.py (scikit-learn
    /// GaussianProcessRegressor, zero mean) and gpr_dicekriging_oracle.R (DiceKriging km, constant mean), with the
    /// versions pinned in manifest.json. Both libraries form the same matrices from the same data, so predictions,
    /// variances and likelihoods agree to rounding; the tolerance of 1E-10 leaves several orders of magnitude for
    /// differences in summation order and fused multiply-add.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_GaussianProcessRegression
    {
        private static readonly double[,] X = {
                { 2.17738, 4.468612 },
                { 0.422221, 0.85312 },
                { 2.075278, 1.681632 },
                { 2.589626, 1.854974 },
                { 1.116576, 1.362185 },
                { 2.400869, 3.308961 },
                { 4.865415, 1.700204 },
                { 0.325616, 2.560091 },
                { 1.276959, 4.135343 },
                { 1.341062, 2.710866 },
                { 2.306571, 3.137435 },
                { 4.136313, 4.379834 },
                { 1.833662, 3.562892 },
                { 0.146983, 2.986456 },
                { 4.691899, 4.698327 },
                { 2.858201, 3.11807 },
                { 4.388905, 2.004548 },
                { 4.923774, 0.655181 },
                { 4.544085, 3.590644 },
                { 0.604429, 2.354407 }
            };

        private static readonly double[] Y = {
                0.351994, 0.40821, 0.47164, 0.07751, 0.478611,
                1.263938, -1.495783, 0.700854, 0.722551, 1.36381,
                1.176089, -1.27139, 1.273614, 0.376775, -1.4263,
                0.704797, -1.399868, -0.755925, -0.561013, 0.765501
            };

        private static readonly double[,] XTest = {
                { 4.895609, 4.584463 },
                { 3.311492, 4.835238 },
                { 4.614585, 2.685911 },
                { 3.686731, 4.749443 },
                { 1.109237, 2.813448 }
            };

        /// <summary>
        /// Returns the first feature of a two-column array.
        /// </summary>
        private static double[] FirstColumn(double[,] x)
        {
            var result = new double[x.GetLength(0)];
            for (int i = 0; i < result.Length; i++) result[i] = x[i, 0];
            return result;
        }

        /// <summary>
        /// Asserts that two arrays agree elementwise within an absolute tolerance.
        /// </summary>
        private static void AssertArray(double[] expected, double[] actual, double tolerance, string label)
        {
            Assert.HasCount(expected.Length, actual, label);
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], actual[i], tolerance, $"{label}[{i}]");
        }

        #region scikit-learn, Fixed Hyperparameters

        /// <summary>
        /// Predictive means and latent variances agree with scikit-learn GaussianProcessRegressor
        /// (ConstantKernel(1.3) * kernel, alpha = 0.05, optimizer = None) for each covariance function.
        /// </summary>
        [TestMethod]
        [DataRow(GaussianProcessKernel.SquaredExponential)]
        [DataRow(GaussianProcessKernel.Matern32)]
        [DataRow(GaussianProcessKernel.Matern52)]
        [DataRow(GaussianProcessKernel.Exponential)]
        public void Test_FixedHyperparameters_MatchScikitLearn(GaussianProcessKernel kernel)
        {
            var expected = new Dictionary<GaussianProcessKernel, (double[] Mean, double[] Variance)>
            {
                [GaussianProcessKernel.SquaredExponential] = (new double[] { -1.1127663537671013, -0.4604556753573034, -0.9536015918309001, -0.9766035131606234, 1.1917048089092517 }, new double[] { 0.11480053638022512, 0.685352997352126, 0.07300488298698184, 0.35138274363018357, 0.05880634591730805 }),
                [GaussianProcessKernel.Matern32] = (new double[] { -1.1480876728372387, -0.25263722765765345, -0.9973917506709256, -0.7100386759579167, 1.1570820272031264 }, new double[] { 0.25936234741888115, 0.9942478326403716, 0.34203666010974504, 0.7019147522439588, 0.22461735955131393 }),
                [GaussianProcessKernel.Matern52] = (new double[] { -1.1539773351302574, -0.30133111685036734, -0.9766047557680853, -0.7891942848632747, 1.1717630496275309 }, new double[] { 0.1883266934648995, 0.927063103820427, 0.226011800382395, 0.5868157606259463, 0.14255578199131747 }),
                [GaussianProcessKernel.Exponential] = (new double[] { -1.0327096143735268, -0.16621273974225487, -0.9812034746376436, -0.5324662582286298, 1.037597307547525 }, new double[] { 0.5982533184422584, 1.1111268666271539, 0.7003965830836418, 0.9616314294727016, 0.5738773196012075 }),
            }[kernel];
            var gp = new GaussianProcessRegression(X, Y, kernel, GaussianProcessMean.Zero);
            gp.SetHyperparameters(1.3, new[] { 0.7, 1.4 }, 0.05);
            var (mean, variance) = gp.PredictWithVariance(XTest);
            AssertArray(expected.Mean, mean, 1E-10, "mean");
            AssertArray(expected.Variance, variance, 1E-10, "variance");
            AssertArray(expected.Mean, gp.Predict(XTest), 1E-10, "Predict");
        }

        /// <summary>
        /// The log marginal likelihood and its analytic gradient agree with scikit-learn's
        /// log_marginal_likelihood(theta, eval_gradient=True) for ConstantKernel * kernel + WhiteKernel.
        /// </summary>
        [TestMethod]
        [DataRow(GaussianProcessKernel.SquaredExponential)]
        [DataRow(GaussianProcessKernel.Matern32)]
        [DataRow(GaussianProcessKernel.Matern52)]
        [DataRow(GaussianProcessKernel.Exponential)]
        public void Test_LogMarginalLikelihoodAndGradient_MatchScikitLearn(GaussianProcessKernel kernel)
        {
            var expected = new Dictionary<GaussianProcessKernel, (double[] Theta, double Value, double[] Gradient)>
            {
                [GaussianProcessKernel.SquaredExponential] = (new double[] { 0.26236426446749106, -0.35667494393873245, 0.3364722366212129, -2.995732273553991 }, -17.353567599429812, new double[] { -3.9556504319712587, 8.074717132218396, -1.0716460834247847, -1.3476845680630494 }),
                [GaussianProcessKernel.Matern32] = (new double[] { 0.26236426446749106, -0.35667494393873245, 0.3364722366212129, -2.995732273553991 }, -19.729743150304376, new double[] { -5.592410606350765, 5.1937886718929045, 1.832652000703173, -0.8565991207251593 }),
                [GaussianProcessKernel.Matern52] = (new double[] { 0.26236426446749106, -0.35667494393873245, 0.3364722366212129, -2.995732273553991 }, -18.874995948460196, new double[] { -5.2049826156163475, 6.161658242918673, 1.5715597014179064, -1.034042542484141 }),
                [GaussianProcessKernel.Exponential] = (new double[] { 0.26236426446749106, -0.35667494393873245, 0.3364722366212129, -2.995732273553991 }, -22.068772851111326, new double[] { -6.074820885721332, 2.77041194457008, 1.5381834894073994, -0.5005473080236931 }),
            }[kernel];
            var gp = new GaussianProcessRegression(X, Y, kernel, GaussianProcessMean.Zero);
            var gradient = new double[4];
            double value = gp.EvaluateLogMarginalLikelihood(expected.Theta, gradient);
            Assert.AreEqual(expected.Value, value, 1E-10);
            AssertArray(expected.Gradient, gradient, 1E-10, "gradient");
            Assert.IsFalse(gp.IsTrained);
        }

        /// <summary>
        /// The analytic gradient agrees with central differences for every covariance function and mean, with and
        /// without automatic relevance determination.
        /// </summary>
        [TestMethod]
        public void Test_Gradient_MatchesFiniteDifferences()
        {
            foreach (GaussianProcessKernel kernel in Enum.GetValues(typeof(GaussianProcessKernel)))
            {
                foreach (GaussianProcessMean mean in Enum.GetValues(typeof(GaussianProcessMean)))
                {
                    foreach (bool ard in new[] { true, false })
                    {
                        var gp = new GaussianProcessRegression(X, Y, kernel, mean) { UseARD = ard };
                        var theta = ard ? new[] { 0.2, -0.4, 0.3, -2.5 } : new[] { 0.2, -0.1, -2.5 };
                        var gradient = new double[theta.Length];
                        gp.EvaluateLogMarginalLikelihood(theta, gradient);
                        for (int k = 0; k < theta.Length; k++)
                        {
                            double h = 1E-5;
                            var up = (double[])theta.Clone(); up[k] += h;
                            var down = (double[])theta.Clone(); down[k] -= h;
                            double fd = (gp.EvaluateLogMarginalLikelihood(up) - gp.EvaluateLogMarginalLikelihood(down)) / (2 * h);
                            Assert.AreEqual(fd, gradient[k], 1E-6 * Math.Max(1d, Math.Abs(fd)), $"{kernel} {mean} ARD={ard} k={k}");
                        }
                    }
                }
            }
        }

        #endregion

        #region DiceKriging, Constant Mean

        /// <summary>
        /// With a constant mean and fixed covariance parameters, the generalized least squares mean, log likelihood,
        /// predictions, standard deviations (including the nugget) and leave-one-out values agree with
        /// DiceKriging::km and predict.km(type = "UK"). DiceKriging forms multidimensional Matérn and exponential
        /// kernels as tensor products, so those cases use the first feature only.
        /// </summary>
        [TestMethod]
        [DataRow(GaussianProcessKernel.SquaredExponential)]
        [DataRow(GaussianProcessKernel.Matern52)]
        [DataRow(GaussianProcessKernel.Exponential)]
        public void Test_ConstantMean_MatchesDiceKriging(GaussianProcessKernel kernel)
        {
            var expected = new Dictionary<GaussianProcessKernel, (double Beta, double LogLikelihood, double[] Mean, double[] Sd, double[] LooMean, double[] LooSd)>
            {
                [GaussianProcessKernel.SquaredExponential] = (-0.2533227568495846, -17.207036286389712, new double[] { -1.1422014152327056, -0.5523740019155226, -0.9642646761734249, -1.042565749217851, 1.1890743562663684 }, new double[] { 0.409581468109446, 0.8741755774374709, 0.3512731975770205, 0.645158473740027, 0.3298938559062771 },
                    new double[] {
                    0.8234025631502082, 0.12353915287096223, 0.3537219341535704, 0.40800789700140994, 0.6453737588862574, 1.0671735702824139,
                    -0.9606917798948434, 0.5616720197502254, 0.6864736955206641, 1.0920436391405692, 1.26883683538361, -0.5872929775506722,
                    1.2052106684610255, 0.5241027114460828, -0.5461232422221721, 0.5628768138690396, -0.9617124724879549, -1.3744813865865877,
                    -1.5708138958621498, 0.8658940686821572
                }, new double[] {
                    0.7089292633250356, 0.8527408050038996, 0.5923579420765933, 0.5988085526466984, 0.6583576795945112, 0.3005334408195215,
                    0.5713049053640077, 0.3101177047062815, 0.7427663862634855, 0.5229551068018224, 0.3071767086645191, 0.7423989796914751,
                    0.4085698955945149, 0.48548464872876895, 0.7513924606111904, 0.573973927483676, 0.6920845484163884, 0.7740862512250658,
                    0.6220663488510881, 0.3820101604888891
                }),
                [GaussianProcessKernel.Matern52] = (0.06874590423084427, -26.48166773684785, new double[] { -1.1043774142201899, 0.17746283023219653, -1.1165668475697144, -0.5805645007934437, 0.6072695913584623 }, new double[] { 0.2684408794860683, 0.6876181856393597, 0.26473139574229554, 0.6848645268177032, 0.28485650094887865 },
                    new double[] {
                    0.9121135625835328, 0.7083244604686663, 0.8814667911805174, 1.1615025032291313, 0.8245537633231134, 0.5491465993465364,
                    -0.9021266896786133, 0.4551578676707135, 1.0406669220914162, 0.7958203640935679, 0.6411278455244476, -1.0206535445942537,
                    0.6371858304600188, 0.5544688201326462, -0.9677649328886871, -0.3273620895721353, -0.9410827486724371, -1.4580930832409866,
                    -1.4711323181560059, 0.3783537021266793
                }, new double[] {
                    0.27321176380067663, 0.2838637712192886, 0.2879023749052121, 0.3218198703932474, 0.35314352125455345, 0.27969139716546865,
                    0.28425493972742294, 0.2827589149363555, 0.28184097141343734, 0.3061171008931187, 0.27093970674262013, 0.5285521441558658,
                    0.4061486551242107, 0.4469091101054316, 0.2869968771122185, 0.5452083944126658, 0.30907380111753247, 0.33055036492141954,
                    0.29030849260255565, 0.3713797131696663
                }),
                [GaussianProcessKernel.Exponential] = (0.07911300306020555, -19.627741141988526, new double[] { -1.1037220192685548, 0.05808157217322882, -0.9957927640229468, -0.43535931349581014, 0.5058570211833321 }, new double[] { 0.3578122487026808, 0.9840989508446468, 0.45761938902734434, 0.9827359769723157, 0.3470212962411389 },
                    new double[] {
                    0.796509558738645, 0.6923737425209513, 0.6637206663529045, 0.9642371474162564, 0.7394952173710758, 0.7459217663123714,
                    -0.8979182017868791, 0.4012930098812178, 1.072652616735857, 0.741832865709982, 0.8246178487709575, -0.8173973353459854,
                    0.6428735714721238, 0.5342797309135372, -0.9390757077484818, -0.012614046942782009, -0.8193821185584054, -1.3390454811663268,
                    -1.3647745925728931, 0.41670093206129283
                }, new double[] {
                    0.5317521221551125, 0.5537005596668014, 0.5812299449747842, 0.685957268653342, 0.7113757060045433, 0.552639051401738,
                    0.49036520380840076, 0.5523529359724013, 0.4966834409673532, 0.5397596319774406, 0.5232539789468355, 0.8659341997197199,
                    0.7950629648154715, 0.7839053123416456, 0.6023963347931631, 0.8823881982446433, 0.6491451115405622, 0.5398548420803767,
                    0.5911962627965036, 0.7389540319326197
                }),
            }[kernel];
            bool twoDimensional = kernel == GaussianProcessKernel.SquaredExponential;
            var gp = twoDimensional
                ? new GaussianProcessRegression(X, Y, kernel, GaussianProcessMean.Constant)
                : new GaussianProcessRegression(FirstColumn(X), Y, kernel, GaussianProcessMean.Constant);
            gp.SetHyperparameters(1.3, twoDimensional ? new[] { 0.7, 1.4 } : new[] { 0.7 }, 0.05);
            Assert.AreEqual(expected.Beta, gp.MeanCoefficients[0], 1E-10);
            Assert.AreEqual(expected.LogLikelihood, gp.LogMarginalLikelihood, 1E-10);
            var (mean, variance) = twoDimensional ? gp.PredictWithVariance(XTest, true) : gp.PredictWithVariance(FirstColumn(XTest), true);
            AssertArray(expected.Mean, mean, 1E-10, "mean");
            AssertArray(expected.Sd, variance.Select(Math.Sqrt).ToArray(), 1E-10, "sd");
            var (looMean, looVariance) = gp.LeaveOneOut();
            AssertArray(expected.LooMean, looMean, 1E-10, "loo mean");
            AssertArray(expected.LooSd, looVariance.Select(Math.Sqrt).ToArray(), 1E-10, "loo sd");
        }

        #endregion

        #region Hyperparameter Estimation

        /// <summary>
        /// Training reaches the scikit-learn optimum for a Matérn 5/2 kernel with automatic relevance determination
        /// and an estimated noise level (50 L-BFGS-B restarts in scikit-learn).
        /// </summary>
        /// <remarks>
        /// The optimizer stops when the largest projected gradient of the log likelihood per observation is below
        /// 1E-6, which places the log hyperparameters within about 1E-5 of the optimum on this well-conditioned
        /// problem; 1E-4 allows for scikit-learn's own stopping error. The log likelihood is flat at the optimum, so
        /// its tolerance is much tighter.
        /// </remarks>
        [TestMethod]
        public void Test_Train_MatchesScikitLearnOptimum()
        {
            var expectedTheta = new double[] { -0.21507325548282558, 0.7533513496347233, 0.4717007385151898, -5.149215521352981 };
            var gp = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52, GaussianProcessMean.Zero) { NumberOfRestarts = 10 };
            gp.Train();
            Assert.AreEqual(-8.686455631164588, gp.LogMarginalLikelihood, 1E-7);
            var theta = new[] { Math.Log(gp.SignalVariance), Math.Log(gp.LengthScales[0]), Math.Log(gp.LengthScales[1]), Math.Log(gp.NoiseVariance) };
            AssertArray(expectedTheta, theta, 1E-4, "theta");
            Assert.AreEqual(0d, gp.Jitter);
        }

        /// <summary>
        /// Training with a constant mean reaches the DiceKriging maximum likelihood fit (covtype "gauss",
        /// nugget.estim = TRUE, multistart = 20).
        /// </summary>
        [TestMethod]
        public void Test_Train_MatchesDiceKrigingOptimum()
        {
            var gp = new GaussianProcessRegression(X, Y, GaussianProcessKernel.SquaredExponential, GaussianProcessMean.Constant) { NumberOfRestarts = 10 };
            gp.Train();
            Assert.AreEqual(-6.802974579515503, gp.LogMarginalLikelihood, 1E-7);
        }

        /// <summary>
        /// Training is reproducible for a fixed seed, and a fixed noise variance is held.
        /// </summary>
        [TestMethod]
        public void Test_Train_SeededAndFixedNoise()
        {
            var a = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern32) { NumberOfRestarts = 3, PRNGSeed = 7 };
            var b = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern32) { NumberOfRestarts = 3, PRNGSeed = 7 };
            a.Train();
            b.Train();
            Assert.AreEqual(a.LogMarginalLikelihood, b.LogMarginalLikelihood, 0d);
            CollectionAssert.AreEqual(a.LengthScales, b.LengthScales);

            var c = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern32);
            c.Train(0.01);
            Assert.AreEqual(0.01, c.NoiseVariance, 0d);
            Assert.IsLessThanOrEqualTo(a.LogMarginalLikelihood + 1E-6, c.LogMarginalLikelihood);
        }

        /// <summary>
        /// A shared length scale gives one length scale, and the isotropic fit cannot exceed the ARD fit.
        /// </summary>
        [TestMethod]
        public void Test_Train_SharedLengthScale()
        {
            var iso = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52, GaussianProcessMean.Zero) { UseARD = false, NumberOfRestarts = 5 };
            iso.Train();
            var ard = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52, GaussianProcessMean.Zero) { NumberOfRestarts = 5 };
            ard.Train();
            Assert.HasCount(1, iso.LengthScales);
            Assert.HasCount(2, ard.LengthScales);
            Assert.IsLessThanOrEqualTo(ard.LogMarginalLikelihood + 1E-6, iso.LogMarginalLikelihood);
        }

        #endregion

        #region Behavior and Edge Cases

        /// <summary>
        /// Without noise the model interpolates: the predictive mean reproduces each training value and the latent
        /// variance there is zero.
        /// </summary>
        [TestMethod]
        public void Test_ZeroNoise_Interpolates()
        {
            var gp = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52);
            gp.SetHyperparameters(1.0, new[] { 0.5, 0.5 }, 0.0);
            var (mean, variance) = gp.PredictWithVariance(X);
            AssertArray(Y, mean, 1E-8, "mean");
            AssertArray(new double[Y.Length], variance, 1E-8, "variance");
        }

        /// <summary>
        /// A linear mean recovers an exactly linear response, with the coefficients in the order [β₀, β₁, β₂], and
        /// extrapolates it.
        /// </summary>
        [TestMethod]
        public void Test_LinearMean_RecoversPlane()
        {
            var y = new double[X.GetLength(0)];
            for (int i = 0; i < y.Length; i++) y[i] = 2.0 - 0.5 * X[i, 0] + 3.0 * X[i, 1];
            var gp = new GaussianProcessRegression(X, y, GaussianProcessKernel.SquaredExponential, GaussianProcessMean.Linear);
            gp.SetHyperparameters(1.0, new[] { 1.0, 1.0 }, 1E-6);
            AssertArray(new[] { 2.0, -0.5, 3.0 }, gp.MeanCoefficients, 1E-8, "beta");
            var far = new double[,] { { 20.0, -10.0 } };
            Assert.AreEqual(2.0 - 10.0 - 30.0, gp.Predict(far)[0], 1E-6);
        }

        /// <summary>
        /// Duplicate inputs with zero noise make the covariance matrix singular; the factorization adds jitter, which
        /// is reported, and predictions stay finite.
        /// </summary>
        [TestMethod]
        public void Test_DuplicateInputs_UseJitter()
        {
            var x = new double[] { 0.0, 1.0, 1.0, 2.0, 3.0 };
            var y = new double[] { 0.0, 0.8, 0.8, 0.9, 0.1 };
            var gp = new GaussianProcessRegression(x, y, GaussianProcessKernel.SquaredExponential);
            gp.SetHyperparameters(1.0, new[] { 1.0 }, 0.0);
            Assert.IsGreaterThan(0d, gp.Jitter);
            var (mean, variance) = gp.PredictWithVariance(new[] { 1.0, 1.5 });
            Assert.AreEqual(0.8, mean[0], 1E-3);
            Assert.IsFalse(double.IsNaN(mean[1]) || double.IsInfinity(mean[1]));
            Assert.IsTrue(variance.All(v => v >= 0d));
        }

        /// <summary>
        /// A single observation gives the conditional normal in closed form.
        /// </summary>
        [TestMethod]
        public void Test_SingleObservation()
        {
            var gp = new GaussianProcessRegression(new[] { 0.0 }, new[] { 2.0 }, GaussianProcessKernel.SquaredExponential, GaussianProcessMean.Zero);
            gp.SetHyperparameters(1.0, new[] { 1.0 }, 0.5);
            var (mean, variance) = gp.PredictWithVariance(new[] { 1.0 });
            double k = Math.Exp(-0.5);
            Assert.AreEqual(k * 2.0 / 1.5, mean[0], 1E-14);
            Assert.AreEqual(1.0 - k * k / 1.5, variance[0], 1E-14);
            Assert.AreEqual(-0.5 * 4.0 / 1.5 - 0.5 * Math.Log(1.5) - 0.5 * Math.Log(2 * Math.PI), gp.LogMarginalLikelihood, 1E-14);
        }

        /// <summary>
        /// The prediction interval is centered on the mean with the normal quantile of the predictive standard deviation.
        /// </summary>
        [TestMethod]
        public void Test_PredictionInterval()
        {
            var gp = new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52);
            gp.SetHyperparameters(1.3, new[] { 0.7, 1.4 }, 0.05);
            var interval = gp.Predict(new Matrix(XTest), 0.05, true);
            var (mean, variance) = gp.PredictWithVariance(XTest, true);
            for (int i = 0; i < mean.Length; i++)
            {
                Assert.AreEqual(mean[i], interval[i, 1], 1E-14);
                Assert.AreEqual(mean[i] - 1.959963984540054 * Math.Sqrt(variance[i]), interval[i, 0], 1E-12);
                Assert.AreEqual(mean[i] + 1.959963984540054 * Math.Sqrt(variance[i]), interval[i, 2], 1E-12);
            }
        }

        /// <summary>
        /// The blocked parallel path (n above the parallel threshold) agrees with a direct dense solution.
        /// </summary>
        [TestMethod]
        public void Test_LargeProblem_MatchesDirectSolution()
        {
            int n = 900;
            var x = new double[n, 2];
            var y = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i, 0] = (i * 0.6180339887498949) % 1.0;
                x[i, 1] = (i * 0.7548776662466927) % 1.0;
                y[i] = Math.Sin(6 * x[i, 0]) + x[i, 1];
            }
            var gp = new GaussianProcessRegression(x, y, GaussianProcessKernel.Matern52, GaussianProcessMean.Zero);
            gp.SetHyperparameters(1.0, new[] { 0.3, 0.4 }, 0.01);
            var k = new Matrix(n);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                {
                    double r = Math.Sqrt(Math.Pow((x[a, 0] - x[b, 0]) / 0.3, 2) + Math.Pow((x[a, 1] - x[b, 1]) / 0.4, 2));
                    k[a, b] = (1 + Math.Sqrt(5) * r + 5 * r * r / 3) * Math.Exp(-Math.Sqrt(5) * r) + (a == b ? 0.01 : 0);
                }
            var chol = new CholeskyDecomposition(k);
            var alpha = chol.Solve(new Vector(y));
            double expectedLml = -0.5 * Vector.DotProduct(new Vector(y), alpha) - 0.5 * chol.LogDeterminant() - 0.5 * n * Math.Log(2 * Math.PI);
            Assert.AreEqual(expectedLml, gp.LogMarginalLikelihood, 1E-8 * Math.Abs(expectedLml));
            var point = new double[,] { { 0.25, 0.75 } };
            double rp, expectedMean = 0;
            for (int b = 0; b < n; b++)
            {
                rp = Math.Sqrt(Math.Pow((0.25 - x[b, 0]) / 0.3, 2) + Math.Pow((0.75 - x[b, 1]) / 0.4, 2));
                expectedMean += (1 + Math.Sqrt(5) * rp + 5 * rp * rp / 3) * Math.Exp(-Math.Sqrt(5) * rp) * alpha[b];
            }
            Assert.AreEqual(expectedMean, gp.Predict(point)[0], 1E-8);
        }

        /// <summary>
        /// Invalid arguments and use before training throw.
        /// </summary>
        [TestMethod]
        public void Test_InvalidArguments()
        {
            Assert.ThrowsExactly<ArgumentNullException>(() => new GaussianProcessRegression((double[,])null!, Y));
            Assert.ThrowsExactly<ArgumentException>(() => new GaussianProcessRegression(X, new double[3]));
            Assert.ThrowsExactly<ArgumentException>(() => new GaussianProcessRegression(new double[0, 2], new double[0]));
            Assert.ThrowsExactly<ArgumentException>(() => new GaussianProcessRegression(new[] { 1.0, double.NaN }, new[] { 1.0, 2.0 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new GaussianProcessRegression(X, Y, (GaussianProcessKernel)4));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => new GaussianProcessRegression(X, Y, GaussianProcessKernel.Matern52, (GaussianProcessMean)3));
            Assert.ThrowsExactly<InvalidOperationException>(() => new GaussianProcessRegression(X, Y) { NumberOfRestarts = -1 }.Train());
            var gp = new GaussianProcessRegression(X, Y);
            Assert.ThrowsExactly<InvalidOperationException>(() => gp.Predict(XTest));
            Assert.ThrowsExactly<InvalidOperationException>(() => gp.LeaveOneOut());
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.SetHyperparameters(0.0, new[] { 1.0, 1.0 }, 0.1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.SetHyperparameters(1.0, new[] { 1.0, 1.0, 1.0 }, 0.1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.SetHyperparameters(1.0, new[] { 1.0, -1.0 }, 0.1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.SetHyperparameters(1.0, new[] { 1.0, 1.0 }, -0.1));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.Train(-1.0));
            Assert.ThrowsExactly<ArgumentException>(() => gp.EvaluateLogMarginalLikelihood(new double[3]));
            gp.SetHyperparameters(1.0, new[] { 1.0, 1.0 }, 0.1);
            Assert.ThrowsExactly<ArgumentException>(() => gp.Predict(new double[,] { { 1.0 } }));
            Assert.ThrowsExactly<ArgumentException>(() => gp.Predict(new double[] { 1.0 }));
            Assert.ThrowsExactly<ArgumentOutOfRangeException>(() => gp.Predict(new Matrix(XTest), 1.5));
        }

        #endregion
    }
}
