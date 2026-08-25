using System;
using System.Collections.Generic;
using System.Reflection;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Characterization tests for the NUTS gradient memo.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// Consecutive leapfrog leaves chain: a leaf ends at a position whose gradient it has just
    /// evaluated, and the next leaf opens by evaluating the gradient at that same position. The
    /// sampler memoizes those values, which is only defensible if it changes nothing. These tests are
    /// the evidence for that claim, so their assertions are bitwise on
    /// <see cref="BitConverter.DoubleToInt64Bits(double)"/> rather than on a tolerance. An
    /// epsilon-based assertion would not separate "identical" from "very close", which is the only
    /// distinction being made here.
    /// </para>
    /// <para>
    /// The reference bit patterns were captured from the sampler before the memo was added. They are
    /// the same on .NET Framework 4.8.1 and on .NET 8, 9, and 10: this target's gradient is exact and
    /// the trajectory arithmetic is IEEE-754 addition, multiplication, and division, so no
    /// transcendental last-place difference reaches the recorded draws on this configuration.
    /// </para>
    /// <para>
    /// Both metric settings are covered. With <see cref="NUTS.AdaptMassMatrix"/> enabled, which is the
    /// default, the metric changes at the end of every adaptation window and each change re-runs the
    /// step-size heuristic, so the memo is consulted under several different metrics within one run.
    /// The gradient of the log-density is a function of position alone and never sees the metric, so an
    /// entry stays valid across a metric change; this test is what holds that reasoning to account.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_NUTS_GradientReuse
    {
        /// <summary>The number of leading recorded draws compared against the reference, per chain.</summary>
        private const int GoldenDraws = 20;

        /// <summary>The gradient evaluation count of the adapting-metric configuration.</summary>
        private const int AdaptedMetricGradientCalls = 4611;

        /// <summary>The gradient evaluation count of the fixed-metric configuration.</summary>
        private const int FixedMetricGradientCalls = 11809;

        /// <summary>The per-coordinate standard deviations of the target.</summary>
        /// <remarks>
        /// Powers of two, so that the gradient <c>-x / s^2</c> is exact and the reference values cannot
        /// drift on a division.
        /// </remarks>
        private static readonly double[] Sd = { 0.25d, 1d, 2d, 8d };

        /// <summary>
        /// The log-density of the zero-mean diagonal Gaussian target, up to an additive constant.
        /// </summary>
        /// <param name="x">The parameter vector.</param>
        /// <returns>The log-density.</returns>
        private static double LogLikelihood(double[] x)
        {
            double sum = 0d;
            for (int j = 0; j < Sd.Length; j++)
            {
                double z = x[j] / Sd[j];
                sum += -0.5d * z * z;
            }
            return sum;
        }

        /// <summary>
        /// The exact gradient of <see cref="LogLikelihood(double[])"/>.
        /// </summary>
        /// <param name="x">The parameter vector.</param>
        /// <returns>The gradient.</returns>
        private static Vector Gradient(IList<double> x)
        {
            var g = new Vector(Sd.Length);
            for (int j = 0; j < Sd.Length; j++)
                g[j] = -x[j] / (Sd[j] * Sd[j]);
            return g;
        }

        /// <summary>
        /// Builds the fixed configuration the reference values were captured from.
        /// </summary>
        /// <param name="adapt">Whether to adapt the diagonal mass matrix.</param>
        /// <param name="gradient">The gradient delegate.</param>
        /// <returns>An unsampled sampler.</returns>
        /// <remarks>
        /// Two chains, run serially so that ordering cannot vary, and a seeded PRNG. The analytic
        /// gradient is supplied so that the delegate invocation count is exactly the number of gradient
        /// evaluations the sampler asked for.
        /// </remarks>
        private static NUTS BuildSampler(bool adapt, HMC.Gradient gradient)
        {
            var priors = new List<IUnivariateDistribution>();
            for (int j = 0; j < Sd.Length; j++) priors.Add(new Uniform(-50d, 50d));

            return new NUTS(priors, LogLikelihood, maxTreeDepth: 6, gradientFunction: gradient)
            {
                NumberOfChains = 2,
                ParallelizeChains = false,
                InitialIterations = 8,
                ThinningInterval = 1,
                WarmupIterations = 60,
                Iterations = 150,
                OutputLength = 100,
                PRNGSeed = 12345,
                AdaptMassMatrix = adapt
            };
        }

        /// <summary>
        /// Runs the fixed configuration and asserts every compared draw against its reference bit pattern.
        /// </summary>
        /// <param name="adapt">Whether to adapt the diagonal mass matrix.</param>
        /// <param name="expected">The reference bit patterns, by chain, then draw, then coordinate.</param>
        private static void AssertGoldenDraws(bool adapt, long[] expected)
        {
            var sampler = BuildSampler(adapt, Gradient);
            sampler.Sample();

            int k = 0;
            for (int c = 0; c < sampler.Output.Length; c++)
            {
                Assert.IsGreaterThanOrEqualTo(GoldenDraws, sampler.Output[c].Count,
                    $"Chain {c} recorded only {sampler.Output[c].Count} draws; the reference covers {GoldenDraws}.");

                for (int i = 0; i < GoldenDraws; i++)
                {
                    var values = sampler.Output[c][i].Values;
                    for (int j = 0; j < values.Length; j++, k++)
                    {
                        long actual = BitConverter.DoubleToInt64Bits(values[j]);
                        Assert.AreEqual(expected[k], actual,
                            $"Chain {c} draw {i} coordinate {j} moved: expected {BitConverter.Int64BitsToDouble(expected[k]):R}, got {values[j]:R}. " +
                            "The gradient memo must not move a single draw.");
                    }
                }
            }

            Assert.AreEqual(expected.Length, k, "The reference length and the compared draw count disagree.");
        }

        /// <summary>
        /// Runs the fixed configuration with a counting gradient delegate.
        /// </summary>
        /// <param name="adapt">Whether to adapt the diagonal mass matrix.</param>
        /// <returns>The number of gradient delegate invocations.</returns>
        private static int CountGradientEvaluations(bool adapt)
        {
            int count = 0;
            var sampler = BuildSampler(adapt, (x) => { count++; return Gradient(x); });
            sampler.Sample();
            return count;
        }

        /// <summary>
        /// With the metric adapting, which is the default, the recorded draws must match the reference
        /// bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_GradientReuse_AdaptedMetric_ReproducesReferenceDrawsBitwise()
        {
            AssertGoldenDraws(true, AdaptedMetricDraws);
        }

        /// <summary>
        /// With a fixed metric the recorded draws must match the reference bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_GradientReuse_FixedMetric_ReproducesReferenceDrawsBitwise()
        {
            AssertGoldenDraws(false, FixedMetricDraws);
        }

        /// <summary>
        /// The memo must remove the redundant gradient evaluation at the join between consecutive
        /// leapfrog leaves, and must remove nothing else.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Before the memo this configuration evaluated the gradient 8,568 times with the
        /// metric adapting and 22,280 times with a fixed metric, against 4611 and
        /// 11809 after it: a reduction of 1.86x and 1.89x, close to the 2x ceiling of removing
        /// one of the two gradient evaluations every leapfrog step makes. The counts are asserted
        /// exactly rather than as an inequality: the companion bitwise tests already pin the trajectory,
        /// so with the trajectory fixed the count is decided by the memo policy alone, and an exact
        /// assertion catches a silent change in hit rate that a threshold would not.
        /// </para>
        /// <para>
        /// The counts are the same on all four target frameworks, as are the draws.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_GradientReuse_ReducesGradientEvaluations()
        {
            Assert.AreEqual(AdaptedMetricGradientCalls, CountGradientEvaluations(true),
                "The adapted-metric gradient evaluation count changed; before the memo it was 8,568.");
            Assert.AreEqual(FixedMetricGradientCalls, CountGradientEvaluations(false),
                "The fixed-metric gradient evaluation count changed; before the memo it was 22,280.");
        }

        /// <summary>
        /// The memo must serve a repeat of the identical position and must not serve a position that
        /// differs only in the sign of a zero.
        /// </summary>
        /// <remarks>
        /// <para>
        /// This is a white-box probe of the memo's lookup, driven through reflection because the memo is
        /// private and the situation it guards against cannot be provoked reliably through a real chain.
        /// It exists because the obvious simplification of the lookup, comparing with <c>==</c> instead
        /// of on the bit pattern, would pass every other test in this class: <c>+0 == -0</c> is true, so
        /// an <c>==</c> lookup would answer a query at negative zero with the gradient at positive zero.
        /// For this target those two gradients differ, in the sign of their own zero, and the difference
        /// would then propagate into the momentum.
        /// </para>
        /// <para>
        /// Negative zero is constructed from its bit pattern rather than written as a literal so that
        /// nothing about the constant depends on how the compiler folds a unary minus.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_GradientReuse_TreatsNegativeZeroAsADistinctPosition()
        {
            int count = 0;
            var sampler = BuildSampler(false, (x) => { count++; return Gradient(x); });
            sampler.Sample();

            // Empty chain zero's memo so the probe cannot hit an entry the run left behind.
            var occupiedField = typeof(NUTS).GetField("_gradientCacheOccupied", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(occupiedField, "The private field _gradientCacheOccupied was not found on NUTS.");
            var occupied = occupiedField!.GetValue(sampler) as bool[][];
            Assert.IsNotNull(occupied, "The private field _gradientCacheOccupied was not a bool[][].");
            Array.Clear(occupied![0], 0, occupied[0].Length);

            var evaluate = typeof(NUTS).GetMethod("EvaluateGradient", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(evaluate, "The private method EvaluateGradient was not found on NUTS.");

            double negativeZero = BitConverter.Int64BitsToDouble(long.MinValue);
            var atPositiveZero = new double[Sd.Length];
            var atNegativeZero = new double[Sd.Length];
            for (int j = 0; j < Sd.Length; j++)
            {
                atPositiveZero[j] = 0d;
                atNegativeZero[j] = negativeZero;
            }

            count = 0;
            var first = evaluate!.Invoke(sampler, new object[] { atPositiveZero, 0 }) as double[];
            Assert.IsNotNull(first, "EvaluateGradient did not return a double[].");
            var reference = (double[])first!.Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            var repeated = evaluate.Invoke(sampler, new object[] { (double[])atPositiveZero.Clone(), 0 }) as double[];
            Assert.IsNotNull(repeated, "EvaluateGradient did not return a double[].");
            Assert.AreEqual(1, count, "Repeating the identical position must be served from the memo.");
            for (int j = 0; j < Sd.Length; j++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(reference[j]), BitConverter.DoubleToInt64Bits(repeated![j]),
                    $"Coordinate {j} came back from the memo with different bits than the evaluation that filled it.");
            }

            var flipped = evaluate.Invoke(sampler, new object[] { atNegativeZero, 0 }) as double[];
            Assert.IsNotNull(flipped, "EvaluateGradient did not return a double[].");
            Assert.AreEqual(2, count,
                "Negative zero is a different position from positive zero and must not be answered from the memo.");
            for (int j = 0; j < Sd.Length; j++)
            {
                Assert.AreNotEqual(BitConverter.DoubleToInt64Bits(reference[j]), BitConverter.DoubleToInt64Bits(flipped![j]),
                    $"Coordinate {j} returned the gradient at positive zero for a query at negative zero.");
            }
        }

        /// <summary>
        /// Reference draw bit patterns with the diagonal mass matrix adapting, captured before the
        /// gradient memo was added.
        /// </summary>
        private static readonly long[] AdaptedMetricDraws =
        {
            -4628717957276900768L, -4614344352116551053L, 4603611391276814710L, 4622225145915450648L,
            -4628717957276900768L, -4614344352116551053L, 4603611391276814710L, 4622225145915450648L,
            4595421539249117930L, 4609283367556249551L, -4619553774880583316L, 4610874494819852500L,
            4595421539249117930L, 4609283367556249551L, -4619553774880583316L, 4610874494819852500L,
            4595268059699645008L, -4614515372391026388L, -4615358376842143412L, 4618502018394842424L,
            -4626965350435671000L, -4617098977443399769L, -4613529980504667989L, -4611552404961281768L,
            -4626965350435671000L, -4617098977443399769L, -4613529980504667989L, -4611552404961281768L,
            -4625117031024711957L, -4618975755754580538L, 4610911898721537679L, -4619588748106623692L,
            -4626833903520120148L, -4618788725822154797L, -4617624839434109446L, -4615653340781801740L,
            -4628586846217905589L, -4624617629968939766L, 4607978285003208601L, 4621524210250680202L,
            -4628586846217905589L, -4624617629968939766L, 4607978285003208601L, 4621524210250680202L,
            -4631048037799758788L, 4598606562225593095L, 4592275416877696720L, -4598889720843471700L,
            -4631048037799758788L, 4598606562225593095L, 4592275416877696720L, -4598889720843471700L,
            -4630614677190815564L, 4599829856108105510L, 4602035654696767498L, 4622517052814789598L,
            -4630614677190815564L, 4599829856108105510L, 4602035654696767498L, 4622517052814789598L,
            -4628044555644799119L, -4615544604587728724L, -4605324568804896625L, 4615289631917756594L,
            4586141649926781016L, -4615659434416111236L, -4606765315056840162L, 4619070661431980192L,
            -4623994016138285324L, 4604694158858439959L, 4600862861138731824L, 4619894943440118595L,
            -4624997371801789776L, 4606985563011054683L, 4609022620809182642L, 4622044230572136596L,
            -4624997371801789776L, 4606985563011054683L, 4609022620809182642L, 4622044230572136596L,
            4589218988890646141L, 4601529455578585566L, 4612990821753793139L, 4599039655050837472L,
            4589218988890646141L, 4601529455578585566L, 4612990821753793139L, 4599039655050837472L,
            4591202640630172407L, 4611551731127349298L, -4615612418592602340L, 4618519724767551228L,
            4591202640630172407L, 4611551731127349298L, -4615612418592602340L, 4618519724767551228L,
            4601900831011982165L, -4653844673522834944L, -4609404599743355790L, 4613069599123696405L,
            4594877573097428800L, 4600646389199820926L, -4617865255456648110L, 4623949388794256850L,
            -4629536297211435422L, 4596843827571823539L, -4615522644836689135L, 4617146743729294822L,
            -4629536297211435422L, 4596843827571823539L, -4615522644836689135L, 4617146743729294822L,
            -4664804292599706624L, -4626728832200013664L, 4612284455809246765L, 4621384453701015334L,
            4592195647805789165L, -4626623762681769642L, 4595481563303407552L, -4600832715318815973L,
            4592195647805789165L, -4626623762681769642L, 4595481563303407552L, -4600832715318815973L,
            4592195647805789165L, -4626623762681769642L, 4595481563303407552L, -4600832715318815973L,
            4592195647805789165L, -4626623762681769642L, 4595481563303407552L, -4600832715318815973L,
            4589036769505293890L, 4607652900249041740L, 4612330264484591542L, -4598014163315121310L,
            4589036769505293890L, 4607652900249041740L, 4612330264484591542L, -4598014163315121310L,
            -4630619197172725760L, 4607461456558991839L, 4616368273643204144L, 4624443866127377968L,
            -4625034992050830473L, -4615566671929947463L, -4615817665372380310L, -4614254007529137496L,
            -4633193357778560152L, -4615566018286652080L, 4615300278428133156L, -4602192214965783574L,
            -4633193357778560152L, -4615566018286652080L, 4615300278428133156L, -4602192214965783574L,
            -4634323884716628752L, 4607777917423490830L, -4611052078354232235L, 4617262316064901896L,
        };

        /// <summary>
        /// Reference draw bit patterns with a fixed identity metric, captured before the gradient memo
        /// was added.
        /// </summary>
        private static readonly long[] FixedMetricDraws =
        {
            4594975404019207961L, -4621711899945443789L, -4611952112987838184L, -4611500315287884938L,
            -4628256621123045462L, -4619287892620029858L, -4614232458884030478L, -4612488514663140402L,
            4595761952827806092L, 4600751122092123102L, 4613575793801799735L, 4624992171110232310L,
            -4633345527287176903L, -4624171562507790115L, -4615839794461924276L, 4626190194575195907L,
            -4625189100783473061L, -4618536524957887010L, -4615614065028290685L, 4626263576848944373L,
            -4626259148557215207L, 4607595391508306782L, -4613028951595474552L, 4624957502953452018L,
            4600173594495267507L, -4614617300851356384L, -4624189941312186700L, -4606082020433159272L,
            4597716792141818872L, -4614241438875454975L, -4621289566791541011L, -4606369784396057216L,
            -4624226452947540604L, -4632705290174055120L, 4600388820663058666L, -4603567199260082839L,
            4595030593496324994L, -4622399403773605971L, -4607493426312042384L, -4603197985821305989L,
            -4624299062973987490L, 4596732760167193896L, -4614965948047856994L, -4601755088583506868L,
            4598764634486221846L, -4618308293726385899L, -4611492594737092578L, -4602437470651895100L,
            4598398684193455921L, -4628829198274197006L, 4607299749312890132L, -4602334364034222595L,
            -4638029692586604976L, 4601898400980151357L, 4608354046106682907L, -4602296705490575856L,
            -4640156845448054636L, -4621005638293633998L, 4605491041205431264L, -4601745596801913269L,
            4580794748882690968L, -4630632412424480978L, 4599148763407829511L, -4601027697149221266L,
            4595370424142256248L, 4589932628438169415L, -4622895218238951981L, 4612827535321258879L,
            -4645439079761559424L, -4636650727009762936L, -4619918594546849499L, 4612789608014404110L,
            -4645439079761559424L, -4636650727009762936L, -4619918594546849499L, 4612789608014404110L,
            -4628794829880472607L, 4579388095185454488L, -4618882751248197878L, 4612846310818148172L,
            -4625710289702829098L, 4603182349032855910L, -4621938938232760059L, -4617799282642129463L,
            4583418010629183568L, -4616037850412337266L, 4613423682573598488L, 4606801000288630048L,
            -4633885671984424928L, -4613768735541523568L, 4613915903054395787L, 4597687828877945318L,
            -4633899563109034347L, 4609344896051078928L, -4609273226815769781L, -4604900278076377306L,
            -4628638329354169868L, -4614623471288870178L, -4631850860977198964L, -4602141301970495431L,
            -4628638329354169868L, -4614623471288870178L, -4631850860977198964L, -4602141301970495431L,
            -4629008656552818505L, -4612834134734929207L, -4627303695876047340L, -4602210861703509402L,
            -4641736452535869472L, 4604706186948601715L, -4605196490149347326L, -4603642802112773977L,
            -4624276442547037690L, 4607005947237830839L, -4605441975554914160L, -4603435550651477501L,
            4595696255895012454L, 4602594494893305626L, -4609367765923473060L, -4602244045001584135L,
            -4634704556634067176L, -4617306856281886408L, -4614470594343065330L, -4601884596191411264L,
            -4635705013991812368L, 4585025165220297424L, -4615059973892280391L, -4601753650874562102L,
            -4635641498537308012L, -4632169441432132214L, 4611340934453912532L, -4603328932724678622L,
            -4635680531729715384L, 4571273438108619616L, -4614340023978301181L, -4604687755738060340L,
            -4635680531729715384L, 4571273438108619616L, -4614340023978301181L, -4604687755738060340L,
            -4622031799990434282L, 4604340633299106995L, -4614172774610055267L, -4604855512262537533L,
            4596885573681383132L, 4607201543281707569L, -4614327399435879440L, -4605153905120656034L,
            -4626105143823512956L, -4620259699054322338L, 4597645397151358219L, -4604402949357508067L,
            4600037919184754832L, -4613908372478531545L, 4608276699040816935L, -4604149254759479244L,
            4560105449497328640L, -4613686306857414974L, 4606789166016571720L, -4603618241203553513L,
        };
    }
}
