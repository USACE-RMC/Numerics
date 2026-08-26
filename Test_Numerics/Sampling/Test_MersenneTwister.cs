using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Sampling;
using System.Linq;
using System.Threading.Tasks;

namespace Sampling
{
    /// <summary>
    /// Unit tests for the Mersenne Twister class. 
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item>Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil</item>
    ///     </list>
    /// </para>
    /// <b> References: </b>
    /// Mersenne Twister with improved initialization. <see href="http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html"/>
    /// </remarks>
    [TestClass]
    public class Test_MersenneTwister
    {
        /// <summary>
        /// Test against output file found here: <see href="http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html"/> 
        /// </summary>
        [TestMethod]
        public void Test_MT19937()
        {
            // Initialize from array
            int[] seeds = { 291, 564, 837, 1110 };
            var prng = new MersenneTwister(seeds);

            // Test Int32
            var trueInt32 = new uint[] { 1067595299, 955945823, 477289528, 4107218783, 4228976476, 3344332714, 3355579695, 227628506, 810200273, 2591290167 };
            for (int i = 0; i < 10; i++)
            {
               Assert.AreEqual(trueInt32[i], prng.GenRandInt32(), 1E-8);
            }
            // run to 1000
            for (int i = 10; i < 1000; i++)
            {
                prng.GenRandInt32();
            }
            // Test NextDouble()
            var trueDouble = new double[] { 0.76275443, 0.99000644, 0.98670464, 0.10143112, 0.27933125, 0.69867227, 0.94218740, 0.03427201, 0.78842173, 0.28180608 };
            for (int i = 0; i < 10; i++)
            {
                Assert.AreEqual(trueDouble[i], prng.NextDouble(), 1E-8);
            }

        }
        /// <summary>
        /// The array-seeded generator reproduces the canonical reference stream, pinning the
        /// generator itself in place.
        /// </summary>
        /// <remarks>
        /// Reference values verified against Python numpy 2.4.2, whose legacy RandomState seeds
        /// with the canonical mt19937ar init_by_array (Matsumoto and Nishimura, 2002):
        /// numpy.random.RandomState([0x123, 0x234, 0x345, 0x456]) reproduces this stream.
        /// </remarks>
        [TestMethod]
        public void Test_ArraySeed_CanonicalReferenceStream()
        {
            var rng = new MersenneTwister(new int[] { 0x123, 0x234, 0x345, 0x456 });
            var expected = new uint[] { 1067595299U, 955945823U, 477289528U, 4107218783U, 4228976476U, 3344332714U };
            for (int i = 0; i < expected.Length; i++)
                Assert.AreEqual(expected[i], rng.GenRandInt32(), $"draw {i}");
        }

        /// <summary>
        /// Unseeded generators constructed back to back or concurrently produce distinct streams,
        /// even inside a single clock tick.
        /// </summary>
        [TestMethod]
        public void Test_UnseededConstruction_DistinctStreams()
        {
            int count = 200;
            var first = new double[count];
            for (int i = 0; i < count; i++)
                first[i] = new MersenneTwister().NextDouble();
            Assert.HasCount(count, first.Distinct().ToList());

            var parallelFirst = new double[count];
            Parallel.For(0, count, i => { parallelFirst[i] = new MersenneTwister().NextDouble(); });
            Assert.HasCount(count, parallelFirst.Distinct().ToList());
        }

    }
}
