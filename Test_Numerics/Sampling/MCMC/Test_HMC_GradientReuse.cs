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
    /// Characterization tests for the HMC start-gradient memo.
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
    /// An HMC trajectory closes with a gradient evaluation at its final position, and the next
    /// transition opens by evaluating the gradient at that same position when the proposal was accepted,
    /// or at the unchanged chain state when it was rejected. The sampler memoizes the closing evaluation
    /// per chain, which is only defensible if it changes nothing. These tests are the evidence for that
    /// claim, so their assertions are bitwise on <see cref="BitConverter.DoubleToInt64Bits(double)"/>
    /// rather than on a tolerance: an epsilon-based assertion would not separate "identical" from "very
    /// close", which is the only distinction being made here.
    /// </para>
    /// <para>
    /// The reference bit patterns were captured from the sampler before the memo was added. They are the
    /// same on .NET Framework 4.8.1 and on .NET 8, 9, and 10: this target's gradient is exact, the
    /// standard deviations are powers of two, and the trajectory arithmetic is IEEE-754 addition,
    /// multiplication, and division, so no transcendental last-place difference reaches the recorded
    /// draws on this configuration.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_HMC_GradientReuse
    {
        /// <summary>The number of leading recorded draws compared against the reference, per chain.</summary>
        private const int GoldenDraws = 20;

        /// <summary>The gradient evaluation count of the fixed configuration.</summary>
        private const int GradientCalls = 5143;

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
        /// <param name="gradient">The gradient delegate.</param>
        /// <returns>An unsampled sampler.</returns>
        /// <remarks>
        /// Two chains, run serially so that ordering cannot vary, and a seeded PRNG. The analytic
        /// gradient is supplied so that the delegate invocation count is exactly the number of gradient
        /// evaluations the sampler asked for. The acceptance rates of this configuration are 89.5% and
        /// 87.5%, so both the accept path, which stores the closing evaluation, and the reject path,
        /// which relies on the entry surviving, are exercised.
        /// </remarks>
        private static HMC BuildSampler(HMC.Gradient gradient)
        {
            var priors = new List<IUnivariateDistribution>();
            for (int j = 0; j < Sd.Length; j++) priors.Add(new Uniform(-50d, 50d));

            return new HMC(priors, LogLikelihood, stepSize: 0.25, steps: 12, gradientFunction: gradient)
            {
                NumberOfChains = 2,
                ParallelizeChains = false,
                InitialIterations = 8,
                ThinningInterval = 1,
                WarmupIterations = 60,
                Iterations = 150,
                OutputLength = 100,
                PRNGSeed = 12345
            };
        }

        /// <summary>
        /// The recorded draws must match the reference bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_HMC_GradientReuse_ReproducesReferenceDrawsBitwise()
        {
            var sampler = BuildSampler(Gradient);
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
                        Assert.AreEqual(ReferenceDraws[k], actual,
                            $"Chain {c} draw {i} coordinate {j} moved: expected {BitConverter.Int64BitsToDouble(ReferenceDraws[k]):R}, got {values[j]:R}. " +
                            "The start-gradient memo must not move a single draw.");
                    }
                }
            }

            Assert.AreEqual(ReferenceDraws.Length, k, "The reference length and the compared draw count disagree.");
        }

        /// <summary>
        /// The memo must remove the opening gradient evaluation of every transition after a chain's
        /// first, and must remove nothing else.
        /// </summary>
        /// <remarks>
        /// <para>
        /// Before the memo this configuration evaluated the gradient 5,541 times; each of the two chains
        /// runs 200 transitions, and the memo removes one evaluation from all but the first of them,
        /// which is 5,541 - 2 x 199 = 5,143. The count is asserted exactly rather than as an inequality:
        /// the companion bitwise test already pins the trajectory, so with the trajectory fixed the count
        /// is decided by the memo policy alone, and an exact assertion catches a silent change in hit
        /// rate that a threshold would not.
        /// </para>
        /// <para>
        /// The count is the same on all four target frameworks, as are the draws.
        /// </para>
        /// </remarks>
        [TestMethod]
        public void Test_HMC_GradientReuse_ReducesGradientEvaluations()
        {
            int count = 0;
            var sampler = BuildSampler((x) => { count++; return Gradient(x); });
            sampler.Sample();
            Assert.AreEqual(GradientCalls, count,
                "The gradient evaluation count changed; before the memo it was 5,541. Check the companion " +
                "bitwise test first: if it also moved the trajectory changed and this count followed it, " +
                "which is a different finding from a change in memo hit rate.");
        }

        /// <summary>
        /// A gradient delegate that writes through its argument must not be able to corrupt the chain
        /// state or the recorded draws.
        /// </summary>
        /// <remarks>
        /// The trajectory works on a private copy of the chain state, so the delegate never receives the
        /// chain's own array. Without that copy, the delegate's write would land in the current state's
        /// array: a rejected transition would then return, and record, the sentinel the delegate wrote
        /// rather than the state the chain was actually at. The sentinel sits inside the prior support on
        /// purpose, so nothing downstream can clamp or reject it away; only the private copy keeps it out
        /// of the record.
        /// </remarks>
        [TestMethod]
        public void Test_HMC_GradientReuse_DelegateWritingThroughItsArgumentCannotCorruptTheChain()
        {
            const double sentinel = 43d;
            var sampler = BuildSampler((x) =>
            {
                var g = Gradient(x);
                // Overwrite the argument after computing its gradient, as a delegate that uses its
                // argument as scratch space would.
                for (int j = 0; j < x.Count; j++) x[j] = sentinel;
                return g;
            });
            sampler.Sample();

            for (int c = 0; c < sampler.MarkovChains.Length; c++)
            {
                for (int i = 0; i < sampler.MarkovChains[c].Count; i++)
                {
                    var values = sampler.MarkovChains[c][i].Values;
                    for (int j = 0; j < values.Length; j++)
                    {
                        Assert.AreNotEqual(sentinel, values[j],
                            $"Chain {c} draw {i} coordinate {j} recorded the delegate's sentinel; the " +
                            "delegate was handed the chain state's own array rather than a private copy.");
                    }
                }
            }
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
        public void Test_HMC_GradientReuse_TreatsNegativeZeroAsADistinctPosition()
        {
            int count = 0;
            var sampler = BuildSampler((x) => { count++; return Gradient(x); });
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
            var reference = (double[])InvokeStartGradient(sampler, atPositiveZero).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            var repeated = InvokeStartGradient(sampler, (double[])atPositiveZero.Clone());
            Assert.AreEqual(1, count, "Repeating the identical position must be served from the memo.");
            AssertSameBits(reference, repeated,
                "came back from the memo with different bits than the evaluation that filled it");

            var flipped = InvokeStartGradient(sampler, atNegativeZero);
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
        /// A gradient delegate that writes through its argument, or a caller that mutates the array it
        /// passed, must not be able to corrupt the memo key: a memo holding the queried array by
        /// reference would find its key rewritten under it, so the entry would stop matching the point
        /// it was computed at and start matching a point it was not. Storing a copy makes the key
        /// immutable. This probe reproduces the hazard directly: evaluate at a position, mutate the
        /// caller's array, and check that the entry still answers the original position and does not
        /// answer the mutated one.
        /// </remarks>
        [TestMethod]
        public void Test_HMC_GradientReuse_StoresACopyOfTheQueriedPosition()
        {
            int count = 0;
            var sampler = BuildSampler((x) => { count++; return Gradient(x); });
            sampler.Sample();
            ClearMemo(sampler);

            var probe = new double[] { 0.5d, -1.5d, 2.25d, -4d };
            var original = (double[])probe.Clone();

            count = 0;
            var reference = (double[])InvokeStartGradient(sampler, probe).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            // Mutate the caller's array in place, exactly as the trajectory does to its working copy.
            for (int j = 0; j < probe.Length; j++) probe[j] += 1d;

            var atOriginal = InvokeStartGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count,
                "The memo stopped recognising the position it was filled at after the caller's array was mutated; " +
                "it is holding a reference to that array rather than a copy.");
            AssertSameBits(reference, atOriginal, "was not returned unchanged from the memo");

            InvokeStartGradient(sampler, probe);
            Assert.AreEqual(2, count,
                "The memo answered a position it never evaluated at; its stored key followed the caller's mutation.");
        }

        /// <summary>
        /// The memo must record the position before the delegate runs, so a delegate that writes through
        /// its argument cannot leave the memo keyed by the mutated position.
        /// </summary>
        /// <remarks>
        /// <see cref="HMC.GradientFunction"/> receives a position array by reference, and nothing
        /// prevents an implementation from using it as scratch space during the evaluation. Were the
        /// position recorded after such a delegate returned, the entry would be keyed by the mutated
        /// position while holding the gradient of the original one: a later query at the mutated position
        /// would falsely hit, and a later query at the original position would falsely miss. The delegate
        /// here computes the gradient of its argument and then overwrites the argument, and the
        /// assertions observe which keys the memo answers to afterwards.
        /// </remarks>
        [TestMethod]
        public void Test_HMC_GradientReuse_RecordsThePositionBeforeTheDelegateRuns()
        {
            int count = 0;
            bool hostile = false;
            var sampler = BuildSampler((x) =>
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
            InvokeStartGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            InvokeStartGradient(sampler, (double[])original.Clone());
            Assert.AreEqual(1, count,
                "A repeat of the original position must be served from the memo; the stored key followed " +
                "the delegate's in-call mutation of its argument.");

            InvokeStartGradient(sampler, (double[])mutated.Clone());
            Assert.AreEqual(2, count,
                "The memo answered the mutated position, which was never evaluated; the position must be " +
                "recorded before the delegate runs.");
        }

        /// <summary>
        /// The memo must store a copy of the gradient the delegate returned, not a reference to the
        /// <see cref="Vector"/> it came back in.
        /// </summary>
        /// <remarks>
        /// <see cref="HMC.GradientFunction"/> is a caller-supplied delegate and nothing obliges it to
        /// allocate a fresh <see cref="Vector"/> per call; returning one reused buffer is a reasonable
        /// thing for a performance-minded implementation to do. A memo holding that buffer by reference
        /// would find its entry rewritten by whatever the buffer holds next. The probe fills the memo,
        /// overwrites the delegate's one buffer directly, and re-queries the position.
        /// </remarks>
        [TestMethod]
        public void Test_HMC_GradientReuse_StoresACopyOfTheReturnedGradient()
        {
            int count = 0;
            var reused = new Vector(Sd.Length);
            var sampler = BuildSampler((x) =>
            {
                count++;
                for (int j = 0; j < Sd.Length; j++) reused[j] = -x[j] / (Sd[j] * Sd[j]);
                return reused;
            });
            sampler.Sample();
            ClearMemo(sampler);

            var position = new double[] { 0.5d, -1.5d, 2.25d, -4d };

            count = 0;
            var reference = (double[])InvokeStartGradient(sampler, position).Clone();
            Assert.AreEqual(1, count, "The first evaluation at a fresh position must reach the gradient delegate.");

            // Overwrite the delegate's one buffer, as its next evaluation anywhere would.
            for (int j = 0; j < Sd.Length; j++) reused[j] = double.NaN;

            var repeated = InvokeStartGradient(sampler, (double[])position.Clone());
            Assert.AreEqual(1, count, "Repeating the identical position must be served from the memo.");
            AssertSameBits(reference, repeated,
                "came back holding the delegate's buffer contents; the memo is aliased to the returned vector");
        }

        /// <summary>
        /// Empties chain zero's start-gradient memo so a probe cannot hit an entry a sampling run left
        /// behind.
        /// </summary>
        /// <param name="sampler">A sampled sampler, so that its memo is allocated.</param>
        private static void ClearMemo(HMC sampler)
        {
            var field = typeof(HMC).GetField("_startGradientOccupied", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(field, "The private field _startGradientOccupied was not found on HMC.");
            var occupied = field!.GetValue(sampler) as bool[];
            Assert.IsNotNull(occupied, "The private field _startGradientOccupied was not a bool[].");
            occupied![0] = false;
        }

        /// <summary>
        /// Calls the sampler's private start-gradient memo for chain zero.
        /// </summary>
        /// <param name="sampler">A sampled sampler.</param>
        /// <param name="position">The position to evaluate at.</param>
        /// <returns>The gradient the memo returned. On a hit the array is owned by the memo.</returns>
        private static double[] InvokeStartGradient(HMC sampler, double[] position)
        {
            var method = typeof(HMC).GetMethod("StartGradient", BindingFlags.NonPublic | BindingFlags.Instance);
            Assert.IsNotNull(method, "The private method StartGradient was not found on HMC.");
            var result = method!.Invoke(sampler, new object[] { 0, position }) as Vector;
            Assert.IsNotNull(result, "StartGradient did not return a Vector.");
            return result!.Array;
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
        /// Reference draw bit patterns, captured before the start-gradient memo was added. Ordered by
        /// chain, then draw, then coordinate.
        /// </summary>
        private static readonly long[] ReferenceDraws =
        {
            4598790761935071958L, -4618981460931078777L, -4633255124910629032L, -4617162068170500676L,
            -4624529944371995155L, -4617749420758822309L, -4626235032242337149L, -4615164204388161845L,
            -4628983265917522267L, -4618908764902722070L, -4625110010327077502L, -4614951484470846344L,
            4598954263896678186L, -4617097697679525926L, -4649749985272669732L, -4618414882460949285L,
            -4624860888538835852L, -4623810609633856776L, 4592026127623997532L, -4615929701390368898L,
            -4621191476304810502L, 4604801501498482217L, -4613202572418810060L, -4605665759072173374L,
            -4629068366436945104L, -4618626112238693940L, 4614750598294948781L, -4605596466180995177L,
            4591244730581707978L, 4596090026709588294L, 4613171933681998900L, -4606128366636146336L,
            4595930585635372142L, -4629060288817124607L, 4612735735303026817L, -4604204730577397544L,
            4594366844203156640L, -4627079147440009976L, -4610561632631668905L, -4600361212597767362L,
            4594366844203156640L, -4627079147440009976L, -4610561632631668905L, -4600361212597767362L,
            -4636155716436921976L, -4624339647361356827L, 4612338231985503694L, -4601005187574736344L,
            -4634973842059191088L, -4620627630749799633L, 4612267695448330660L, -4600836608777488564L,
            -4627659637541364119L, -4622850209188883201L, 4609583962443052758L, -4600612900520531445L,
            -4628528992263354328L, 4600694626607139696L, -4610677816921661479L, -4601882932559591974L,
            -4638620556040634800L, -4619916645768072882L, -4611517246567096082L, -4601317546231322976L,
            -4626808121326507195L, -4620250639947265398L, 4612009310764256968L, -4606366362579769275L,
            -4626808121326507195L, -4620250639947265398L, 4612009310764256968L, -4606366362579769275L,
            4598410045685915280L, -4622157382559180487L, 4611388786152615160L, -4605968032698001500L,
            -4624449941804639026L, 4605092849502181424L, 4609566208919303122L, -4606738776715140986L,
            -4627906937310126020L, -4616720237793338336L, 4611908397114699630L, 4613421123650436916L,
            -4637441883597509066L, 4609533003703851103L, 4611105707947881364L, 4611214443729672069L,
            4592728123375942821L, 4610830707918042185L, 4607185877474830950L, 4613316316417758951L,
            -4628049799448325686L, -4619060982015456801L, -4612859193130813537L, -4608328719210657963L,
            -4631200667548008512L, 4607232683977472203L, -4618021602528479124L, -4603454475928677768L,
            4586358561313321022L, -4627629114948788311L, -4613640643895376738L, -4602369865473133708L,
            4592987823185115278L, -4624290272479213034L, -4612261317450630715L, -4601833935265813524L,
            4599949205595933104L, -4635305139867655008L, 4610930562709882120L, -4599080917220753613L,
            -4623348452407011981L, 4606151590101388335L, -4619642956006844966L, -4599718869619806443L,
            -4623348452407011981L, 4606151590101388335L, -4619642956006844966L, -4599718869619806443L,
            4600515961342593418L, -4615720506001707236L, -4607014327469146078L, -4602706230880192716L,
            4600199425511467748L, -4615856985277434193L, 4616366937767583751L, -4609338645062503337L,
            4593081008405489748L, -4615981226981559728L, 4614939328883772403L, -4609524976034779626L,
            -4629787526051756144L, 4605721730834660398L, -4617157942930840285L, -4600762894074213924L,
            -4640220860144088144L, -4615965449033942918L, -4613160886502937701L, -4600612414082812880L,
            4591101748398142419L, 4600537153890364586L, 4609624430212270331L, -4598857016835340532L,
            4582956471600378512L, 4599167091731948549L, 4608573541217824228L, -4598921277526999901L,
            -4621782872919811902L, 4598058065089355576L, -4617920349762308280L, -4598646684179451416L,
            4600419081759212744L, -4619502348791681535L, -4608041254443229700L, -4597129485098255032L,
            -4621790777160238445L, -4615838884775658076L, -4606515638330080586L, -4597060943261626173L,
        };
    }
}
