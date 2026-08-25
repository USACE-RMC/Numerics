using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Data;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling.MCMC;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Text.RegularExpressions;

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
        /// Scans every Numerics source file and fails if any <c>new ArgumentException(...)</c> has its
        /// arguments inverted: either <c>nameof(...)</c> passed as the first argument, or a quoted bare
        /// identifier (for example <c>"stepSize"</c>) passed as the first argument.
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
            if (sourceRoot == null)
                Assert.Fail("The Numerics source tree was not found next to the compiled test assembly, so this guard could not scan the library sources. " +
                    "This fails hard rather than reporting Inconclusive: a guard that silently stands down when its own lookup breaks (for example, because the " +
                    "checkout was renamed or moved) could let an inverted ArgumentException ship undetected, which defeats its purpose.");

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

                    // The first split argument always starts at the opening parenthesis, in both the masked
                    // and the original text (masking preserves length and delimiters), so this slice recovers
                    // the un-masked source for the first argument alone.
                    string firstArgumentOriginal = text.Substring(open + 1, arguments[0].Length);

                    bool inverted = IsNameofExpression(arguments[0]) || IsBareIdentifierLiteral(firstArgumentOriginal);
                    if (!inverted) continue;

                    int line = masked.Take(start).Count(c => c == '\n') + 1;
                    offenders.Add(file.Substring(sourceRoot.Length).TrimStart(Path.DirectorySeparatorChar) + ":" + line);
                }
            }

            Assert.IsEmpty(offenders,
                "ArgumentException takes (message, paramName). These sites pass the parameter name first, which puts it in Message and the sentence in ParamName: "
                + string.Join(", ", offenders));
        }

        /// <summary>
        /// Determines whether an argument, taken in its entirety once trimmed, is a single <c>nameof(...)</c>
        /// expression.
        /// </summary>
        /// <param name="argument">A masked, comma-split constructor argument.</param>
        /// <returns><see langword="true"/> when the whole trimmed argument is exactly one <c>nameof(...)</c> call.</returns>
        /// <remarks>
        /// Anchoring on the whole argument, rather than only its prefix, avoids flagging a message built as
        /// <c>nameof(Foo) + " is invalid"</c>: that argument starts with <c>nameof</c> but, once the trailing
        /// concatenation is accounted for, is not itself a <c>nameof</c> expression.
        /// </remarks>
        private static bool IsNameofExpression(string argument)
        {
            string trimmed = argument.Trim();
            if (!trimmed.StartsWith("nameof", StringComparison.Ordinal)) return false;

            int i = "nameof".Length;
            while (i < trimmed.Length && char.IsWhiteSpace(trimmed[i])) i++;
            if (i >= trimmed.Length || trimmed[i] != '(') return false;

            int close = FindMatchingParenthesis(trimmed, i);
            return close == trimmed.Length - 1;
        }

        /// <summary>
        /// Determines whether an argument's original (unmasked) source text is a quoted, bare identifier,
        /// such as <c>"stepSize"</c>.
        /// </summary>
        /// <param name="originalArgument">The corresponding argument text taken from the unmasked source.</param>
        /// <returns><see langword="true"/> when the argument is a string literal containing a single identifier
        /// with no spaces or sentence punctuation.</returns>
        /// <remarks>
        /// A genuine message is a sentence: it has spaces and terminal punctuation. A quoted bare identifier in
        /// the message slot is almost certainly a parameter name that was written as a string literal instead
        /// of passed as <c>nameof(...)</c>, or otherwise placed in the wrong argument position. This mirrors
        /// the shape of <c>new ArgumentException("stepSize", "The leapfrog step size must be positive.")</c>,
        /// which the <c>nameof</c>-prefix check alone cannot see because its first argument is not a
        /// <c>nameof</c> expression at all.
        /// </remarks>
        private static bool IsBareIdentifierLiteral(string originalArgument)
        {
            return Regex.IsMatch(originalArgument.Trim(), "^\"[A-Za-z_][A-Za-z0-9_]*\"$", RegexOptions.None);
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

        /// <summary>
        /// Verifies the corrected argument order on the <c>stepSize</c> validation site in <see cref="NUTS"/>,
        /// which the <c>nameof</c>-prefix scan could not see because its first argument was the string literal
        /// <c>"stepSize"</c> rather than a <c>nameof(...)</c> expression.
        /// </summary>
        [TestMethod]
        public void Test_NUTS_ReportsParameterName()
        {
            var priors = new List<IUnivariateDistribution> { new Uniform(-50d, 50d), new Uniform(-50d, 50d) };
            static double logLH(double[] x) => -0.5d * (x[0] * x[0] + x[1] * x[1]);
            var sampler = new NUTS(priors, logLH, stepSize: 0d);

            var exception = AssertThrows<ArgumentException>(() => sampler.Sample());
            Assert.AreEqual("stepSize", exception.ParamName);
            StringAssert.Contains(exception.Message, "The leapfrog step size must be positive.");
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

                if (c == '"')
                {
                    i++;
                    while (i < text.Length && text[i] != '"')
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

                if (c == '\'')
                {
                    int closeQuote = FindCharLiteralClose(text, i);
                    if (closeQuote < 0)
                    {
                        // Not a valid char literal (for example, an apostrophe inside a preprocessor directive
                        // or other non-literal text). Leave it as ordinary text rather than masking everything
                        // up to the next unrelated quote character.
                        i++;
                        continue;
                    }

                    for (int k = i + 1; k < closeQuote; k++)
                    {
                        if (text[k] != '\n' && text[k] != '\r') masked[k] = ' ';
                    }
                    i = closeQuote + 1;
                    continue;
                }

                i++;
            }

            return masked.ToString();
        }

        /// <summary>
        /// Finds the offset of the closing quote of the char literal that opens at <paramref name="openQuote"/>,
        /// if the text at that position is actually a well-formed char literal.
        /// </summary>
        /// <param name="text">The source text.</param>
        /// <param name="openQuote">The offset of the opening <c>'</c>.</param>
        /// <returns>The offset of the matching closing <c>'</c>, or -1 when the content is not a valid char
        /// literal (most commonly a stray apostrophe that is not part of any literal).</returns>
        /// <remarks>
        /// A well-formed char literal contains exactly one ordinary character, or one backslash escape
        /// (<c>\n</c>, <c>\t</c>, <c>\\</c>, <c>\'</c>, and so on, or a <c>\xH..H</c>, <c>\uHHHH</c>, or
        /// <c>\UHHHHHHHH</c> hexadecimal escape), and never spans a line. Requiring this exact shape, rather
        /// than treating any later <c>'</c> as the terminator, keeps an unrelated apostrophe (for example, one
        /// inside a preprocessor directive, which is not itself masked as a comment or string) from being
        /// mistaken for the start of a char literal and masking everything up to the next unrelated quote.
        /// </remarks>
        private static int FindCharLiteralClose(string text, int openQuote)
        {
            int i = openQuote + 1;
            if (i >= text.Length || text[i] == '\n' || text[i] == '\r') return -1;

            if (text[i] == '\\')
            {
                i++;
                if (i >= text.Length || text[i] == '\n' || text[i] == '\r') return -1;
                char escape = text[i];
                i++;

                int maxDigits = escape == 'U' ? 8 : escape == 'u' || escape == 'x' ? 4 : 0;
                int digits = 0;
                while (digits < maxDigits && i < text.Length && Uri.IsHexDigit(text[i]))
                {
                    i++;
                    digits++;
                }
            }
            else
            {
                i++;
            }

            return i < text.Length && text[i] == '\'' ? i : -1;
        }

        /// <summary>
        /// Finds the start offset of every <c>new {typeName}</c> construction in masked source text.
        /// </summary>
        /// <param name="masked">Masked source text.</param>
        /// <param name="typeName">The exception type name to search for.</param>
        /// <returns>The offsets at which each construction begins.</returns>
        /// <remarks>
        /// Matches on a regular expression rather than the literal token <c>"new " + typeName</c> so that
        /// unusual but legal spacing (extra spaces, a line break between <c>new</c> and the type) and a
        /// namespace-qualified type name (for example <c>System.ArgumentException</c>) are still found.
        /// </remarks>
        private static IEnumerable<int> FindConstructions(string masked, string typeName)
        {
            string pattern = @"\bnew\s+(?:[A-Za-z_][A-Za-z0-9_]*\.)*" + Regex.Escape(typeName) + @"\s*\(";
            foreach (Match match in Regex.Matches(masked, pattern))
                yield return match.Index;
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
