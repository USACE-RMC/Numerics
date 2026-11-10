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
    /// The reference bit patterns were captured from the sampler without the memo, and the memo
    /// reproduces them bit for bit. They use the step size fixed at the smoothed dual-averaging value
    /// before the first sampling transition. They are the same on .NET Framework 4.8.1 and on .NET 8, 9, and 10: this target's gradient is exact and
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
        private const int AdaptedMetricGradientCalls = 4599;

        /// <summary>The gradient evaluation count of the fixed-metric configuration.</summary>
        private const int FixedMetricGradientCalls = 12041;

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
        /// one of the two gradient evaluations every leapfrog step makes. With the step size fixed
        /// before the first sampling transition the counts are 4599 and 12041. The counts are asserted
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
                "The adapted-metric gradient evaluation count changed; before the memo it was 8,568. " +
                "Check the companion bitwise tests first: if they also moved the trajectory changed and " +
                "this count followed it, which is a different finding from a change in memo hit rate.");
            Assert.AreEqual(FixedMetricGradientCalls, CountGradientEvaluations(false),
                "The fixed-metric gradient evaluation count changed; before the memo it was 22,280. " +
                "Check the companion bitwise tests first: if they also moved the trajectory changed and " +
                "this count followed it, which is a different finding from a change in memo hit rate.");
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
            ClearMemo(sampler);

            double negativeZero = BitConverter.Int64BitsToDouble(long.MinValue);
            var atPositiveZero = new double[Sd.Length];
            var atNegativeZero = new double[Sd.Length];
            for (int j = 0; j < Sd.Length; j++)
            {
                atPositiveZero[j] = 0d;
                atNegativeZero[j] = negativeZero;
            }

            count = 0;
            var reference = (double[])InvokeEvaluateGradient(sampler, atPositiveZero).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            var repeated = InvokeEvaluateGradient(sampler, (double[])atPositiveZero.Clone());
            Assert.AreEqual(1, count, "Repeating the identical position must be served from the memo.");
            AssertSameBits(reference, repeated,
                "came back from the memo with different bits than the evaluation that filled it");

            var flipped = InvokeEvaluateGradient(sampler, atNegativeZero);
            Assert.AreEqual(2, count,
                "Negative zero is a different position from positive zero and must not be answered from the memo.");
            for (int j = 0; j < Sd.Length; j++)
            {
                Assert.AreNotEqual(BitConverter.DoubleToInt64Bits(reference[j]), BitConverter.DoubleToInt64Bits(flipped[j]),
                    $"Coordinate {j} returned the gradient at positive zero for a query at negative zero.");
            }
        }

        /// <summary>
        /// The memo must store a copy of the queried position, not a reference to the caller's array.
        /// </summary>
        /// <remarks>
        /// The requirement is real and not hypothetical: <see cref="NUTS"/>'s in-place leapfrog mutates
        /// the very array it passed in, one statement after the opening half-step. A memo holding that
        /// array by reference would find its key rewritten under it, so the entry would stop matching
        /// the point it was computed at and start matching a point it was not. This probe reproduces
        /// that directly: evaluate at a position, mutate the caller's array, and check that the entry
        /// still answers the original position and does not answer the mutated one.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_GradientReuse_StoresACopyOfTheQueriedPosition()
        {
            int count = 0;
            var sampler = BuildSampler(false, (x) => { count++; return Gradient(x); });
            sampler.Sample();
            ClearMemo(sampler);

            var probe = new double[] { 0.5d, -1.5d, 2.25d, -4d };
            var original = (double[])probe.Clone();

            count = 0;
            var reference = (double[])InvokeEvaluateGradient(sampler, probe).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            // Mutate the caller's array in place, exactly as LeapfrogInPlace does to its position.
            for (int j = 0; j < probe.Length; j++) probe[j] += 1d;

            var atOriginal = InvokeEvaluateGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count,
                "The memo stopped recognising the position it was filled at after the caller's array was mutated; " +
                "it is holding a reference to that array rather than a copy.");
            AssertSameBits(reference, atOriginal, "was not returned unchanged from the memo");

            InvokeEvaluateGradient(sampler, probe);
            Assert.AreEqual(2, count,
                "The memo answered a position it never evaluated at; its stored key followed the caller's mutation.");
        }

        /// <summary>
        /// The memo must record the position before the delegate runs, so a delegate that writes through
        /// its argument cannot leave the memo keyed by the mutated position.
        /// </summary>
        /// <remarks>
        /// <see cref="NUTS.GradientFunction"/> receives the caller's position array by reference, and
        /// nothing prevents an implementation from using it as scratch space during the evaluation. Were
        /// the position recorded after such a delegate returned, the entry would be keyed by the mutated
        /// position while holding the gradient of the original one: a later query at the mutated position
        /// would falsely hit, and a later query at the original position would falsely miss. The delegate
        /// here computes the gradient of its argument and then overwrites the argument, and the assertions
        /// observe which keys the memo answers to afterwards.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_GradientReuse_RecordsThePositionBeforeTheDelegateRuns()
        {
            int count = 0;
            bool hostile = false;
            var sampler = BuildSampler(false, (x) =>
            {
                count++;
                var g = Gradient(x);
                if (hostile)
                {
                    // Overwrite the argument after computing its gradient, as a delegate that uses its
                    // argument as scratch space would.
                    for (int j = 0; j < x.Count; j++) x[j] = 1000d + j;
                }
                return g;
            });
            sampler.Sample();
            ClearMemo(sampler);
            hostile = true;

            var original = new double[] { 0.5d, -1.5d, 2.25d, -4d };
            var mutated = new double[] { 1000d, 1001d, 1002d, 1003d };

            count = 0;
            InvokeEvaluateGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            InvokeEvaluateGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count,
                "A repeat of the original position must be served from the memo; the stored key followed " +
                "the delegate's in-call mutation of its argument.");

            InvokeEvaluateGradient(sampler, (double[])mutated.Clone());
            Assert.AreEqual(2, count,
                "The memo answered the mutated position, which was never evaluated; the position must be " +
                "recorded before the delegate runs.");
        }

        /// <summary>
        /// The memo must store a copy of the gradient the delegate returned, not a reference to the
        /// <see cref="Vector"/> it came back in.
        /// </summary>
        /// <remarks>
        /// <see cref="NUTS.GradientFunction"/> is a caller-supplied delegate and nothing obliges it to
        /// allocate a fresh <see cref="Vector"/> per call; returning one reused buffer is a reasonable
        /// thing for a performance-minded implementation to do. A memo holding that buffer by reference
        /// would have every entry aliased to the same array, so each entry would return the most recent
        /// gradient rather than its own. The delegate here does exactly that, and the probe evaluates at
        /// two distinct positions before re-querying the first.
        /// </remarks>
        [TestMethod]
        public void Test_NUTS_GradientReuse_StoresACopyOfTheReturnedGradient()
        {
            int count = 0;
            var reused = new Vector(Sd.Length);
            var sampler = BuildSampler(false, (x) =>
            {
                count++;
                for (int j = 0; j < Sd.Length; j++) reused[j] = -x[j] / (Sd[j] * Sd[j]);
                return reused;
            });
            sampler.Sample();
            ClearMemo(sampler);

            var first = new double[] { 0.5d, -1.5d, 2.25d, -4d };
            var second = new double[] { 1.5d, -0.5d, 3.25d, -3d };

            count = 0;
            var atFirst = (double[])InvokeEvaluateGradient(sampler, first).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            // This overwrites the delegate's one buffer with the gradient at a different position.
            InvokeEvaluateGradient(sampler, second);
            Assert.AreEqual(2, count, "A different position must reach the gradient delegate.");

            var repeated = InvokeEvaluateGradient(sampler, (double[])first.Clone());
            Assert.AreEqual(2, count, "Repeating the first position must be served from the memo.");
            AssertSameBits(atFirst, repeated,
                "came back as the gradient at a later position; the memo is aliased to the delegate's reused buffer");
        }

        /// <summary>
        /// Empties chain zero's gradient memo so a probe cannot hit an entry a sampling run left behind.
        /// </summary>
        /// <param name="sampler">A sampled sampler, so that its memo is allocated.</param>
        private static void ClearMemo(NUTS sampler)
        {
            var field = typeof(NUTS).GetField("_gradientCacheOccupied", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(field, "The private field _gradientCacheOccupied was not found on NUTS.");
            var occupied = field!.GetValue(sampler) as bool[][];
            Assert.IsNotNull(occupied, "The private field _gradientCacheOccupied was not a bool[][].");
            Array.Clear(occupied![0], 0, occupied[0].Length);
        }

        /// <summary>
        /// Calls the sampler's private gradient memo for chain zero.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="position">The position to evaluate at.</param>
        /// <returns>The gradient the memo returned. This array is owned by the memo.</returns>
        private static double[] InvokeEvaluateGradient(NUTS sampler, double[] position)
        {
            var method = typeof(NUTS).GetMethod("EvaluateGradient", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(method, "The private method EvaluateGradient was not found on NUTS.");
            var result = method!.Invoke(sampler, new object[] { position, 0 }) as double[];
            Assert.IsNotNull(result, "EvaluateGradient did not return a double[].");
            return result!;
        }

        /// <summary>
        /// Asserts two gradients agree on every bit of every coordinate.
        /// </summary>
        /// <param name="expected">The reference gradient.</param>
        /// <param name="actual">The gradient to check.</param>
        /// <param name="what">A description of the failure, completing "Coordinate {j} ...".</param>
        private static void AssertSameBits(double[] expected, double[] actual, string what)
        {
            for (int j = 0; j < expected.Length; j++)
            {
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected[j]), BitConverter.DoubleToInt64Bits(actual[j]),
                    $"Coordinate {j} {what}.");
            }
        }

        /// <summary>
        /// Reference draw bit patterns with the diagonal mass matrix adapting.
        /// </summary>
        private static readonly long[] AdaptedMetricDraws =
        {
            -4628803838224638915L, -4614259607649821298L, 4603662010922677129L, 4622225994844806520L,
            -4628803838224638915L, -4614259607649821298L, 4603662010922677129L, 4622225994844806520L,
            4595458641086650182L, 4609367443672335074L, -4619529541039500751L, 4610870404054544500L,
            4595458641086650182L, 4609367443672335074L, -4619529541039500751L, 4610870404054544500L,
            4595305003787025864L, -4614440700824921740L, -4615354297077804638L, 4618502925082415793L,
            -4626984312180987512L, -4616975988322010404L, -4613527265422439307L, -4611554042944395538L,
            -4626984312180987512L, -4616975988322010404L, -4613527265422439307L, -4611554042944395538L,
            -4625121897131996482L, -4618874470052611135L, 4610910091834939243L, -4619594666327540816L,
            -4626838898629949809L, -4618705313620266169L, -4617622434468878474L, -4615656013683198695L,
            -4628584482840734910L, -4624558171093502998L, 4607978422350908266L, 4621524421351339728L,
            -4628584482840734910L, -4624558171093502998L, 4607978422350908266L, 4621524421351339728L,
            -4631045995772088548L, 4598665552169201780L, 4592274364807965272L, -4598889593688713694L,
            -4631045995772088548L, 4598665552169201780L, 4592274364807965272L, -4598889593688713694L,
            -4630612643845469317L, 4599777464636201594L, 4602035566138858794L, 4622517165546398768L,
            -4630612643845469317L, 4599777464636201594L, 4602035566138858794L, 4622517165546398768L,
            -4628045036671456721L, -4615539936298900246L, -4605324569437090019L, 4615289916824242970L,
            4586142637484533128L, -4615655589922340999L, -4606765315477564720L, 4619070790106958044L,
            -4623993957731875713L, 4604691418380775706L, 4600862860369848520L, 4619895024740155237L,
            -4624997341824395762L, 4606983306134992640L, 4609022620681259727L, 4622044267290431596L,
            -4624997341824395762L, 4606983306134992640L, 4609022620681259727L, 4622044267290431596L,
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
            -4634323884716628752L, 4607777917423490830L, -4611052078354232235L, 4617262316064901896L
        };

        /// <summary>
        /// Reference draw bit patterns with a fixed identity metric.
        /// </summary>
        private static readonly long[] FixedMetricDraws =
        {
            4602855732565887118L, -4626732739343075846L, -4612640245923166201L, 4618160763352123902L,
            4602277981951387782L, -4624170161668334083L, -4611632774089779095L, 4618057012897232669L,
            -4623249125958758910L, -4620820913024834320L, -4614671634791966127L, 4618057360074831607L,
            4594694098991460071L, 4602333501525376794L, -4618848157220165321L, 4614234000573107879L,
            -4628872047729854308L, 4602992273501767377L, -4613807833955270148L, 4614577237255493914L,
            4594371021117727522L, 4599012158553767440L, -4611862696687201453L, 4613840420062162322L,
            4594371021117727522L, 4599012158553767440L, -4611862696687201453L, 4613840420062162322L,
            4585942060083074248L, -4623541182162221402L, -4615167008481005206L, 4615794807005820017L,
            -4622232487202885015L, 4607282377186175520L, -4611367955206881653L, 4618191520223763672L,
            -4629442422201471640L, -4615221373931854389L, 4614711153012019039L, 4618719011867316051L,
            -4629086247559603457L, -4618965650770462370L, 4612378001016413307L, 4620100085384932417L,
            -4629589268209790094L, 4599741382662876812L, 4608735578110874438L, 4620905874111809390L,
            -4634889749221276880L, -4620446802616142140L, 4590359127208608044L, -4605835094595141303L,
            -4635648676814282136L, 4605584867229772981L, -4638878224322713542L, -4603896226126509613L,
            -4635648676814282136L, 4605584867229772981L, -4638878224322713542L, -4603896226126509613L,
            -4636316580429629860L, -4613787283286451675L, 4597605812169304155L, -4604921121440515874L,
            -4626488985499837484L, -4615853935785127962L, 4609806407786666930L, -4605841635179778464L,
            -4626488985499837484L, -4615853935785127962L, 4609806407786666930L, -4605841635179778464L,
            4594882461668625214L, -4619409017741037873L, 4610542978496808683L, -4603778269510759822L,
            -4641072358233889216L, -4613794822582646629L, -4620611617444064458L, -4622071051485566838L,
            4601987273358156671L, 4602908484748259883L, -4610853594350856974L, -4610552338434986915L,
            -4624863264972230947L, 4603488527033387895L, -4612339621013770322L, -4610761186109034125L,
            -4624313688219029240L, 4603333031522645150L, 4610925749764192479L, 4615775574197128862L,
            4598310966234587802L, 4602335174080545758L, 4610624107463044967L, 4615660793982695869L,
            -4634563610865581452L, -4621959175993784810L, -4625767549460945847L, -4605894928342714723L,
            -4635226506452549581L, 4603871466606767542L, -4627142099148941709L, -4605433845155440249L,
            4575442739891999968L, -4618980653469952298L, -4613771600142108599L, 4620504106975429796L,
            4599394657935878806L, -4618669077952662458L, -4614443385872535939L, 4620439789851435806L,
            4599394657935878806L, -4618669077952662458L, -4614443385872535939L, 4620439789851435806L,
            -4635943366103064112L, -4620619672836288196L, -4613825757693661039L, 4620378420858309417L,
            -4651930806027160848L, -4624513524686730203L, 4608848571156319697L, -4630745878194924520L,
            4585607483370533042L, 4602053152345504632L, 4614605257054486957L, -4625491455196346480L,
            -4624280826809067643L, 4594509666975602688L, 4614678843090264315L, 4578262949440779904L,
            -4629293809891282576L, 4597936457069193237L, 4617119945919246943L, -4614605792663109336L,
            4594934226206846383L, -4619229571899503158L, -4606273985741680109L, -4601368005314235483L,
            4599277370276691630L, 4603712062625318153L, -4610711097034143054L, -4600963876379625556L,
            4599277370276691630L, 4603712062625318153L, -4610711097034143054L, -4600963876379625556L,
            -4638304639444041088L, 4597773802553742812L, -4609817400500959216L, -4600830434503728492L,
            4586029131743391262L, 4598266153658111436L, 4616180677486921257L, -4607086704898421080L,
            -4642801670973452168L, -4632970049218469985L, 4617027462772607117L, -4605623437022856945L
        };
    }
}
