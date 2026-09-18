using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Distributions.Univariate
{
    /// <summary>Reads frozen, runtime-independent distribution oracle CSV resources.</summary>
    internal static class DistributionOracle
    {
        /// <summary>Reads the quoted scalar-column fixture format emitted by the standalone R generators.</summary>
        internal static IEnumerable<Dictionary<string, string>> Read(string suffix)
        {
            var assembly = typeof(DistributionOracle).Assembly;
            string name = assembly.GetManifestResourceNames().Single(x => x.EndsWith(suffix, StringComparison.Ordinal));
            using var reader = new StreamReader(assembly.GetManifestResourceStream(name)!);
            string[] columns = reader.ReadLine()!.Split(',').Select(x => x.Trim('"')).ToArray();
            while (!reader.EndOfStream)
            {
                string[] cells = reader.ReadLine()!.Split(',').Select(x => x.Trim('"')).ToArray();
                if (cells.Length != columns.Length) throw new InvalidDataException("Unexpected distribution oracle CSV columns.");
                yield return columns.Select((key, index) => new { key, value = cells[index] }).ToDictionary(x => x.key, x => x.value);
            }
        }

        /// <summary>Parses invariant R numeric output and explicit nonfinite classifications.</summary>
        internal static double Number(string text) => text == "Inf" ? double.PositiveInfinity : text == "-Inf"
            ? double.NegativeInfinity : text == "NA" || text == "NaN" ? double.NaN : double.Parse(text, CultureInfo.InvariantCulture);
    }
}
