using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data;
using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;

namespace Utilities
{
    /// <summary>
    /// Guards the argument order used when constructing <see cref="ArgumentException"/> in the Numerics library.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     <list type="bullet">
    ///     <item> Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil </item>
    ///     </list>
    /// </para>
    /// <para>
    /// <see cref="ArgumentException"/> takes <c>(message, paramName)</c> while
    /// <see cref="ArgumentOutOfRangeException"/> takes <c>(paramName, message)</c>. The two are easy to
    /// confuse, and the inverted form compiles cleanly because both parameters are strings. When inverted,
    /// the identifier lands in <see cref="Exception.Message"/> and the sentence lands in
    /// <see cref="ArgumentException.ParamName"/>, which makes diagnostics misleading and breaks any caller
    /// that keys on <see cref="ArgumentException.ParamName"/>. The defect propagated through the library by
    /// copy-paste, so this class scans the library sources for the inverted shape rather than pinning a
    /// handful of individual call sites.
    /// </para>
    /// </remarks>
    [TestClass]
    public class Test_ArgumentExceptionOrder
    {
        /// <summary>
        /// Scans every Numerics source file and fails if any <c>new ArgumentException(...)</c> passes
        /// <c>nameof(...)</c> as the first argument.
        /// </summary>
        /// <remarks>
        /// The scan masks comments, string literals and character literals before parsing so that text
        /// inside documentation or literals cannot produce a match, and so that parentheses and commas
        /// inside literals cannot confuse the argument split. Only two-argument constructions are
        /// considered; the single-argument overload takes a message and is always correct.
        /// </remarks>
        [TestMethod]
        public void Test_NoInvertedArgumentExceptionArguments()
        {
            string sourceRoot = LocateNumericsSourceRoot();
            if (sourceRoot == null) Assert.Inconclusive("The Numerics source tree was not found next to the compiled test assembly, so the source scan could not run.");

            var offenders = new List<string>();
            foreach (string file in EnumerateSourceFiles(sourceRoot))
            {
                string text = File.ReadAllText(file);
                string masked = MaskCommentsAndLiterals(text);
                foreach (int start in FindConstructions(masked, "ArgumentException"))
                {
                    int open = masked.IndexOf('(', start);
                    if (open < 0) continue;
                    int close = FindMatchingParenthesis(masked, open);
                    if (close < 0) continue;

                    var arguments = SplitTopLevelArguments(masked.Substring(open + 1, close - open - 1));
                    if (arguments.Count < 2) continue;
                    if (arguments[1].Trim().Length == 0) continue;
                    if (!arguments[0].TrimStart().StartsWith("nameof", StringComparison.Ordinal)) continue;

                    int line = masked.Take(start).Count(c => c == '\n') + 1;
                    offenders.Add(file.Substring(sourceRoot.Length).TrimStart(Path.DirectorySeparatorChar) + ":" + line);
                }
            }

            Assert.IsEmpty(offenders,
                "ArgumentException takes (message, paramName). These sites pass nameof(...) first, which puts the identifier in Message and the sentence in ParamName: "
                + string.Join(", ", offenders));
        }

        /// <summary>
        /// Verifies the corrected argument order on a representative site in <see cref="ExtensionMethods"/>.
        /// </summary>
        [TestMethod]
        public void Test_ExtensionMethods_ReportsParameterName()
        {
            var exception = AssertThrows<ArgumentException>(() => new double[] { 1d, 2d }.Add(new double[] { 1d }));
            Assert.AreEqual("array", exception.ParamName);
            StringAssert.Contains(exception.Message, "The arrays must be the same length.");
        }

        /// <summary>
        /// Verifies the corrected argument order on a representative site in <see cref="Vector"/>.
        /// </summary>
        [TestMethod]
        public void Test_Vector_ReportsParameterName()
        {
            var a = new Vector(new[] { 1d, 2d });
            var b = new Vector(new[] { 1d });
            var exception = AssertThrows<ArgumentException>(() => Vector.DotProduct(a, b));
            Assert.AreEqual("Length", exception.ParamName);
            StringAssert.Contains(exception.Message, "The vectors must be the same length.");
        }

        /// <summary>
        /// Verifies the corrected argument order on a representative site in <see cref="Interpolater"/>.
        /// </summary>
        [TestMethod]
        public void Test_Interpolater_ReportsParameterName()
        {
            var exception = AssertThrows<ArgumentException>(() => new Linear(new[] { 1d, 2d, 3d }, new[] { 1d, 2d }));
            Assert.AreEqual("xValues", exception.ParamName);
            StringAssert.Contains(exception.Message, "The x and y lists must be the same length.");
        }

        /// <summary>
        /// Verifies the corrected argument order on a representative site in <see cref="Polynomial"/>.
        /// </summary>
        [TestMethod]
        public void Test_Polynomial_ReportsParameterName()
        {
            var exception = AssertThrows<ArgumentException>(() => new Polynomial(3, new[] { 1d, 2d, 3d }, new[] { 1d, 2d, 3d }));
            Assert.AreEqual("order", exception.ParamName);
            StringAssert.Contains(exception.Message, "The order must be less than the length of the x value list.");
        }

        /// <summary>
        /// Verifies the corrected argument order on a representative site in <see cref="TimeSeries"/>.
        /// </summary>
        [TestMethod]
        public void Test_TimeSeries_ReportsParameterName()
        {
            var series = new TimeSeries(TimeInterval.OneDay, new DateTime(2000, 1, 1), new[] { 1d, 2d, 3d });
            var exception = AssertThrows<ArgumentException>(() => series.MovingAverage(3));
            Assert.AreEqual("period", exception.ParamName);
            StringAssert.Contains(exception.Message, "The period must be less than the length of the time-series.");
        }

        #region Helpers

        /// <summary>
        /// Invokes an action and returns the exception of the expected type, failing the test otherwise.
        /// </summary>
        /// <typeparam name="T">The expected exception type.</typeparam>
        /// <param name="action">The action expected to throw.</param>
        /// <returns>The thrown exception.</returns>
        private static T AssertThrows<T>(Action action) where T : Exception
        {
            try
            {
                action();
            }
            catch (T expected)
            {
                return expected;
            }
            catch (Exception unexpected)
            {
                throw new AssertFailedException("Expected " + typeof(T).Name + " but caught " + unexpected.GetType().Name + ".", unexpected);
            }

            throw new AssertFailedException("Expected " + typeof(T).Name + " but no exception was thrown.");
        }

        /// <summary>
        /// Locates the Numerics library source directory relative to this source file.
        /// </summary>
        /// <param name="callerFilePath">Supplied by the compiler; the full path of this source file.</param>
        /// <returns>The Numerics project directory, or <see langword="null"/> when the source tree is unavailable.</returns>
        private static string LocateNumericsSourceRoot([CallerFilePath] string callerFilePath = "")
        {
            if (string.IsNullOrEmpty(callerFilePath)) return null;

            var directory = new DirectoryInfo(Path.GetDirectoryName(callerFilePath));
            while (directory != null)
            {
                string candidate = Path.Combine(directory.FullName, "Numerics");
                if (File.Exists(Path.Combine(candidate, "Numerics.csproj"))) return candidate;
                directory = directory.Parent;
            }

            return null;
        }

        /// <summary>
        /// Enumerates the C# sources of the library, skipping build output directories.
        /// </summary>
        /// <param name="sourceRoot">The Numerics project directory.</param>
        /// <returns>The full paths of the library sources.</returns>
        private static IEnumerable<string> EnumerateSourceFiles(string sourceRoot)
        {
            string separator = Path.DirectorySeparatorChar.ToString();
            return Directory.GetFiles(sourceRoot, "*.cs", SearchOption.AllDirectories)
                            .Where(f => !f.Contains(separator + "bin" + separator) && !f.Contains(separator + "obj" + separator))
                            .OrderBy(f => f, StringComparer.Ordinal);
        }

        /// <summary>
        /// Replaces the contents of comments, string literals and character literals with spaces.
        /// </summary>
        /// <param name="text">The source text.</param>
        /// <returns>Text of identical length whose comment and literal contents cannot be mistaken for code.</returns>
        /// <remarks>
        /// Delimiters are preserved so that an empty literal remains distinguishable from a missing argument,
        /// and the length is preserved so that offsets and line numbers still refer to the original text.
        /// </remarks>
        private static string MaskCommentsAndLiterals(string text)
        {
            var masked = new StringBuilder(text);
            int i = 0;
            while (i < text.Length)
            {
                char c = text[i];

                if (c == '/' && i + 1 < text.Length && text[i + 1] == '/')
                {
                    while (i < text.Length && text[i] != '\n')
                    {
                        if (text[i] != '\r') masked[i] = ' ';
                        i++;
                    }
                    continue;
                }

                if (c == '/' && i + 1 < text.Length && text[i + 1] == '*')
                {
                    masked[i] = ' ';
                    masked[i + 1] = ' ';
                    i += 2;
                    while (i < text.Length && !(text[i] == '*' && i + 1 < text.Length && text[i + 1] == '/'))
                    {
                        if (text[i] != '\n' && text[i] != '\r') masked[i] = ' ';
                        i++;
                    }
                    if (i < text.Length)
                    {
                        masked[i] = ' ';
                        masked[i + 1] = ' ';
                        i += 2;
                    }
                    continue;
                }

                if (c == '@' && i + 1 < text.Length && text[i + 1] == '"')
                {
                    i += 2;
                    while (i < text.Length)
                    {
                        if (text[i] == '"')
                        {
                            if (i + 1 < text.Length && text[i + 1] == '"')
                            {
                                masked[i] = ' ';
                                masked[i + 1] = ' ';
                                i += 2;
                                continue;
                            }
                            i++;
                            break;
                        }
                        if (text[i] != '\n' && text[i] != '\r') masked[i] = ' ';
                        i++;
                    }
                    continue;
                }

                if (c == '"' || c == '\'')
                {
                    char quote = c;
                    i++;
                    while (i < text.Length && text[i] != quote)
                    {
                        if (text[i] == '\\' && i + 1 < text.Length)
                        {
                            masked[i] = ' ';
                            masked[i + 1] = ' ';
                            i += 2;
                            continue;
                        }
                        if (text[i] != '\n' && text[i] != '\r') masked[i] = ' ';
                        i++;
                    }
                    i++;
                    continue;
                }

                i++;
            }

            return masked.ToString();
        }

        /// <summary>
        /// Finds the start offset of every <c>new {typeName}</c> construction in masked source text.
        /// </summary>
        /// <param name="masked">Masked source text.</param>
        /// <param name="typeName">The exception type name to search for.</param>
        /// <returns>The offsets at which each construction begins.</returns>
        private static IEnumerable<int> FindConstructions(string masked, string typeName)
        {
            string token = "new " + typeName;
            int index = 0;
            while ((index = masked.IndexOf(token, index, StringComparison.Ordinal)) >= 0)
            {
                int after = index + token.Length;
                while (after < masked.Length && char.IsWhiteSpace(masked[after])) after++;
                if (after < masked.Length && masked[after] == '(') yield return index;
                index += token.Length;
            }
        }

        /// <summary>
        /// Finds the offset of the parenthesis matching the one at <paramref name="open"/>.
        /// </summary>
        /// <param name="masked">Masked source text.</param>
        /// <param name="open">The offset of the opening parenthesis.</param>
        /// <returns>The offset of the closing parenthesis, or -1 when unbalanced.</returns>
        private static int FindMatchingParenthesis(string masked, int open)
        {
            int depth = 0;
            for (int i = open; i < masked.Length; i++)
            {
                if (masked[i] == '(') depth++;
                else if (masked[i] == ')')
                {
                    depth--;
                    if (depth == 0) return i;
                }
            }
            return -1;
        }

        /// <summary>
        /// Splits an argument list on its top-level commas.
        /// </summary>
        /// <param name="arguments">The text between the parentheses of a construction.</param>
        /// <returns>The individual arguments, in order.</returns>
        private static List<string> SplitTopLevelArguments(string arguments)
        {
            var parts = new List<string>();
            var current = new StringBuilder();
            int depth = 0;

            foreach (char c in arguments)
            {
                if (c == '(' || c == '[') depth++;
                else if (c == ')' || c == ']') depth--;
                else if (c == ',' && depth == 0)
                {
                    parts.Add(current.ToString());
                    current.Clear();
                    continue;
                }
                current.Append(c);
            }

            parts.Add(current.ToString());
            return parts;
        }

        #endregion
    }
}
