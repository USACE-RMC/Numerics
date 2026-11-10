using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Distributions;
using Numerics.Sampling.MCMC;

namespace Sampling.MCMC
{
    /// <summary>
    /// Unit tests for the uniform draw the samplers use in Metropolis tests and momentum draws.
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
    public class Test_MCMCPositiveUniform
    {
        /// <summary>
        /// A generator that returns a fixed sequence of values from <see cref="Random.NextDouble"/>.
        /// </summary>
        private sealed class SequenceRandom : Random
        {
            private readonly double[] _values;
            private int _next;

            /// <summary>
            /// Creates the generator.
            /// </summary>
            /// <param name="values">The values to return, in order.</param>
            public SequenceRandom(params double[] values) => _values = values;

            /// <summary>Gets the number of draws taken.</summary>
            public int Draws => _next;

            /// <inheritdoc/>
            public override double NextDouble() => _values[_next++];
        }

        /// <summary>
        /// An exact zero becomes 2^-32 and every other value passes through unchanged, one draw each.
        /// </summary>
        [TestMethod]
        public void Test_NextPositiveUniform_MapsOnlyExactZero()
        {
            var random = new SequenceRandom(0d, 0.25d, double.Epsilon, 1d - 1d / 4294967296d);
            Assert.AreEqual(Math.Pow(2d, -32d), MCMCSampler.NextPositiveUniform(random));
            Assert.AreEqual(0.25d, MCMCSampler.NextPositiveUniform(random));
            Assert.AreEqual(double.Epsilon, MCMCSampler.NextPositiveUniform(random));
            Assert.AreEqual(1d - 1d / 4294967296d, MCMCSampler.NextPositiveUniform(random));
            Assert.AreEqual(4, random.Draws);
        }

        /// <summary>
        /// With the mapped value a Metropolis test rejects a zero-density proposal, and the momentum quantile
        /// is the standard normal quantile of 2^-32, about -6.23, rather than the -38.5 clamp at zero.
        /// </summary>
        [TestMethod]
        public void Test_NextPositiveUniform_ZeroNoLongerAcceptsZeroDensity()
        {
            double logU = Math.Log(MCMCSampler.NextPositiveUniform(new SequenceRandom(0d)));
            Assert.AreEqual(-32d * Math.Log(2d), logU, 1E-12);
            Assert.IsGreaterThan(double.NegativeInfinity, logU, "A zero-density proposal must be rejected.");

            double z = Normal.StandardZ(MCMCSampler.NextPositiveUniform(new SequenceRandom(0d)));
            Assert.AreEqual(Math.Pow(2d, -32d), new Normal().CDF(z), 1E-6 * Math.Pow(2d, -32d));
            Assert.IsGreaterThan(-6.3, z);
            Assert.IsLessThan(-6.2, z);
        }
    }
}
