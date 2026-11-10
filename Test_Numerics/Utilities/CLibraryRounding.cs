using System;
using System.Runtime.InteropServices;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Utilities
{
    /// <summary>
    /// Comparisons for reference values recorded bit for bit on Windows.
    /// </summary>
    /// <remarks>
    /// <c>Math.Exp</c>, <c>Log</c>, <c>Pow</c>, <c>Sin</c>, <c>Cos</c>, <c>Tan</c> and <c>Atan2</c> call the platform C math
    /// library: the Universal C Runtime on Windows and, for example, glibc on Linux. The two round some arguments to
    /// neighbouring doubles, so a result pinned bit for bit on Windows can differ in its last bits elsewhere. The pinned
    /// references were recorded on Windows and stay exact there; on other platforms a test states the tolerance, in
    /// units in the last place, that covers the difference measured on Linux.
    /// </remarks>
    internal static class CLibraryRounding
    {
        /// <summary>
        /// Gets whether this process runs on the platform the bitwise references were recorded on (Windows).
        /// </summary>
        public static bool IsReferencePlatform { get; } = RuntimeInformation.IsOSPlatform(OSPlatform.Windows);

        /// <summary>
        /// Returns the number of representable doubles between two finite values.
        /// </summary>
        /// <param name="a">The first value.</param>
        /// <param name="b">The second value.</param>
        /// <returns>The distance in units in the last place; zero for equal values, including 0 and -0.</returns>
        public static long Ulps(double a, double b)
        {
            return Math.Abs(Ordered(a) - Ordered(b));
        }

        /// <summary>
        /// Asserts that a value equals its reference bit for bit on Windows and within <paramref name="ulps"/> units in the
        /// last place on other platforms.
        /// </summary>
        /// <param name="expected">The reference value recorded on Windows.</param>
        /// <param name="actual">The computed value.</param>
        /// <param name="ulps">The tolerance on other platforms, in units in the last place.</param>
        /// <param name="message">The message shown on failure.</param>
        public static void AreEqual(double expected, double actual, long ulps, string message = "")
        {
            if (IsReferencePlatform)
                Assert.AreEqual(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual), message);
            else
                Assert.IsLessThanOrEqualTo(ulps, Ulps(expected, actual), $"{message} expected {expected:R}, actual {actual:R}");
        }

        /// <summary>
        /// Maps a double to a signed integer that increases with the value, so that differences count representable doubles.
        /// </summary>
        /// <param name="value">A finite double.</param>
        /// <returns>The ordered integer.</returns>
        private static long Ordered(double value)
        {
            long bits = BitConverter.DoubleToInt64Bits(value);
            return bits < 0 ? long.MinValue - bits : bits;
        }
    }
}
