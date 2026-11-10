using System.Reflection;
using Numerics;

namespace Utilities
{
    /// <summary>
    /// Supplies the library member that the bitwise reference sweep cannot call directly from this project.
    /// </summary>
    public partial class Test_ToolsBitwiseReference
    {
        /// <summary>The internal <c>Tools.Log1mExp</c>, bound once by reflection.</summary>
        private static readonly Func<double, double> Log1mExpDelegate = typeof(Tools)
            .GetMethod("Log1mExp", BindingFlags.NonPublic | BindingFlags.Static, new[] { typeof(double) })!
            .CreateDelegate<Func<double, double>>();

        /// <summary>Calls the internal <c>Tools.Log1mExp</c> through a delegate; unlike the unit test, this call is not inlined.</summary>
        /// <param name="a">The argument, at most zero.</param>
        /// <returns>log(1 - exp(a)).</returns>
        private static double LibraryLog1mExp(double a) => Log1mExpDelegate(a);
    }
}
