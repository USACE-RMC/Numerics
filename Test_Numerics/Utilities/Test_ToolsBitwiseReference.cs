using System;
using System.Collections.Generic;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Bitwise regression tests for <see cref="Tools.Expm1(double)"/>, <see cref="Tools.Log1p(double)"/> and
    /// <see cref="Tools.Log1mExp(double)"/> against frozen copies of their reference formulas.
    /// </summary>
    /// <remarks>
    /// The library methods may be restructured for speed, but every result must keep the exact bit pattern of the
    /// reference formulas below, including the sign of zero and NaN payloads. The grid covers log-spaced magnitudes
    /// from 1e-320 to 1e308 of both signs, subnormals, signed zeros, infinities, NaNs with payloads, dense walks across
    /// every internal branch boundary, and random bit patterns plus random values in the working range of each method.
    /// The unit tests take 10 magnitudes per decade, walks of 256 bit patterns and 20,000 random inputs of each kind; the
    /// full sweep, with 1000 magnitudes per decade, walks of 4096 bit patterns and one million random inputs of each kind,
    /// runs with <c>dotnet run -c Release --project Benchmarks_Numerics -- sweep ToolsBitwise</c>. Both sides run on the same runtime, so the comparison is exact on every target, and
    /// the library methods are called from lambdas so that the comparison also covers their inlined form.
    /// <para>
    /// Expm1 takes the direct difference <c>exp(x) - 1</c> for <c>|x| &gt; 1/2</c>; against mpmath 1.3.0 its maximum
    /// error there is 1.06 units in the last place, against 2.02 for the compensated formula, and its mean error is
    /// lower in every range tested. The reference therefore uses the direct difference outside <c>|x| &lt;= 1/2</c> and
    /// the compensated formula inside it and for NaN. Log1mExp follows Expm1 for <c>-ln 2 &lt;= a &lt; -1/2</c>, where it
    /// calls Expm1 with <c>|a| &gt; 1/2</c>.
    /// </para>
    /// </remarks>
    [TestClass]
    public partial class Test_ToolsBitwiseReference
    {
        /// <summary>Calls the internal <see cref="Tools.Log1mExp(double)"/> directly, so that it can be inlined.</summary>
        private static double LibraryLog1mExp(double a) => Tools.Log1mExp(a);

        /// <summary>Expm1 matches its reference bit for bit on a reduced grid.</summary>
        [TestMethod]
        public void Expm1MatchesReferenceBitwise()
        {
            Expm1MatchesReferenceBitwiseSweep(10, 256, 20000);
        }

        /// <summary>Log1p matches its reference bit for bit on a reduced grid.</summary>
        [TestMethod]
        public void Log1pMatchesReferenceBitwise()
        {
            Log1pMatchesReferenceBitwiseSweep(10, 256, 20000);
        }

        /// <summary>Log1mExp matches its reference bit for bit on a reduced grid.</summary>
        [TestMethod]
        public void Log1mExpMatchesReferenceBitwise()
        {
            Log1mExpMatchesReferenceBitwiseSweep(10, 256, 20000);
        }
    }
}
