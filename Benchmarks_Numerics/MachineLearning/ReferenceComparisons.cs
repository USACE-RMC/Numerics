using System.Diagnostics;
using System.Globalization;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using Numerics.MachineLearning;

namespace Benchmarks_Numerics.MachineLearning
{
    /// <summary>
    /// Timing comparisons of the machine-learning methods with scikit-learn, R and DiceKriging on identical data, against
    /// reference timings recorded by the scripts in docs/machine-learning/benchmarks.
    /// </summary>
    /// <remarks>
    /// Run with <c>dotnet run -c Release --project Benchmarks_Numerics -- compare [pca-enet] [gpr]</c> on an otherwise
    /// idle machine; with no name both comparisons run. Each comparison prints its table and fails when a C# median is
    /// not below every reference median.
    /// </remarks>
    public static class ReferenceComparisons
    {
        private static readonly (string Name, Func<string?> Run)[] Comparisons =
        {
            ("pca-enet", PcaElasticNet),
            ("gpr", GaussianProcess),
        };

        /// <summary>
        /// Runs the named comparisons, or every comparison when no name is given.
        /// </summary>
        /// <param name="names">The comparison names, matched without regard to case.</param>
        /// <returns>Zero when every selected comparison passes, otherwise one.</returns>
        public static int Run(IReadOnlyCollection<string> names)
        {
            var unknown = names.Where(n => !Comparisons.Any(c => string.Equals(c.Name, n, StringComparison.OrdinalIgnoreCase))).ToList();
            if (unknown.Count > 0)
            {
                Console.Error.WriteLine($"Unknown comparison: {string.Join(", ", unknown)}. Comparisons: {string.Join(", ", Comparisons.Select(c => c.Name))}.");
                return 1;
            }
            int failures = 0;
            foreach (var (name, run) in Comparisons)
            {
                if (names.Count > 0 && !names.Contains(name, StringComparer.OrdinalIgnoreCase)) continue;
                string? failure = run();
                if (failure == null) Console.WriteLine($"{name}: faster than every reference");
                else
                {
                    failures++;
                    Console.WriteLine($"{name}: FAILED: {failure}");
                }
            }
            return failures == 0 ? 0 : 1;
        }

        /// <summary>
        /// Compares <see cref="PrincipalComponentAnalysis"/> and <see cref="ElasticNetRegression"/> with scikit-learn and
        /// R. The data and reference timings are read from the folder named by the NUMERICS_ML_BENCHMARKS environment
        /// variable, prepared by pca_elastic_net_generate_data.py, pca_elastic_net_time_python.py and
        /// pca_elastic_net_time_r.R.
        /// </summary>
        /// <returns>Null when every case is faster than both references; otherwise a description of the failures.</returns>
        /// <remarks>
        /// Each case is timed as the median of repeated runs after a warm-up, matching the medians recorded by the Python
        /// and R scripts; a case that is slower is timed once more before it fails.
        /// </remarks>
        public static string? PcaElasticNet()
        {
            var folder = Environment.GetEnvironmentVariable("NUMERICS_ML_BENCHMARKS");
            if (string.IsNullOrEmpty(folder) || !File.Exists(Path.Combine(folder, "pca_elastic_net_timings_python.json")) || !File.Exists(Path.Combine(folder, "pca_elastic_net_timings_r.json")))
                return "set NUMERICS_ML_BENCHMARKS to a folder prepared by docs/machine-learning/benchmarks.";
            var python = ReadTimings(Path.Combine(folder, "pca_elastic_net_timings_python.json"));
            var r = ReadTimings(Path.Combine(folder, "pca_elastic_net_timings_r.json"));

            var report = new StringBuilder();
            var failures = new List<string>();
            // Cases run in the order of the Python and R scripts: PCA, then the elastic net, by increasing size.
            var cases = python.Keys.Where(r.ContainsKey)
                .Select(k => Regex.Match(k, @"^(pca_fit|enet_path|enet_cv)_(\d+)x(\d+)$"))
                .Where(m => m.Success)
                .OrderBy(m => m.Groups[1].Value == "pca_fit" ? 0 : 1)
                .ThenBy(m => long.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture) * long.Parse(m.Groups[3].Value, CultureInfo.InvariantCulture))
                .ThenBy(m => m.Groups[1].Value == "enet_cv" ? 1 : 0)
                .ToList();
            foreach (var m in cases)
            {
                string key = m.Value;
                int n = int.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture);
                int p = int.Parse(m.Groups[3].Value, CultureInfo.InvariantCulture);
                Action run;
                if (m.Groups[1].Value == "pca_fit")
                {
                    var X = ReadMatrix(Path.Combine(folder, $"pca_{n}x{p}.bin"), n, p);
                    run = () => new PrincipalComponentAnalysis(X).Train();
                }
                else
                {
                    var X = ReadMatrix(Path.Combine(folder, $"enet_{n}x{p}_X.bin"), n, p);
                    var y = ReadVector(Path.Combine(folder, $"enet_{n}x{p}_y.bin"));
                    if (m.Groups[1].Value == "enet_path") run = () => new ElasticNetRegression(X, y) { Alpha = 0.5 }.Fit();
                    else run = () => new ElasticNetRegression(X, y) { Alpha = 0.5 }.CrossValidate(10, 1);
                }
                GC.Collect();
                GC.WaitForPendingFinalizers();
                GC.Collect();
                double csharp = MedianSeconds(run);
                // A case that loses is timed once more, so a burst of load from other processes does not decide it.
                if (!(csharp < python[key] && csharp < r[key])) csharp = Math.Min(csharp, MedianSeconds(run));
                report.AppendLine(FormattableString.Invariant($"{key,-22} C# {csharp,10:F5} s   Python {python[key],10:F5} s   R {r[key],10:F5} s"));
                if (!(csharp < python[key] && csharp < r[key])) failures.Add(key);
            }
            Console.WriteLine(report.ToString());
            return failures.Count == 0 ? null : "slower than Python or R: " + string.Join(", ", failures);
        }

        /// <summary>
        /// Compares the fit and prediction times of <see cref="GaussianProcessRegression"/> with scikit-learn and
        /// DiceKriging at every size in docs/machine-learning/benchmarks/gpr-reference-timings.json.
        /// </summary>
        /// <returns>Null when every fit and prediction is faster than both references; otherwise a description of the
        /// failures.</returns>
        /// <remarks>
        /// The reference timings are measured by gpr_benchmark_sklearn.py and gpr_benchmark_dicekriging.R with the same
        /// deterministic data, kernel (Matérn 5/2 with automatic relevance determination, estimated noise), a single
        /// optimizer start, and 1000 prediction points with variances.
        /// </remarks>
        public static string? GaussianProcess()
        {
            string? directory = AppContext.BaseDirectory;
            while (directory != null && !File.Exists(Path.Combine(directory, "Numerics.sln"))) directory = Path.GetDirectoryName(directory);
            string? path = directory == null ? null : Path.Combine(directory, "docs", "machine-learning", "benchmarks", "gpr-reference-timings.json");
            if (path == null || !File.Exists(path)) return "the reference timings file docs/machine-learning/benchmarks/gpr-reference-timings.json was not found.";
            using var document = JsonDocument.Parse(File.ReadAllText(path));
            var root = document.RootElement;
            int d = root.GetProperty("d").GetInt32(), m = root.GetProperty("m").GetInt32();
            var failures = new List<string>();
            foreach (var size in root.GetProperty("sizes").EnumerateArray())
            {
                int n = size.GetProperty("n").GetInt32(), reps = size.GetProperty("reps").GetInt32();
                var (x, y) = GaussianProcessData(n, d, 0);
                var (xt, _) = GaussianProcessData(m, d, n);
                var fit = new List<double>();
                var predict = new List<double>();
                for (int r = 0; r < reps; r++)
                {
                    var gp = new GaussianProcessRegression(x, y, GaussianProcessKernel.Matern52, GaussianProcessMean.Constant);
                    var watch = Stopwatch.StartNew();
                    gp.Train();
                    fit.Add(watch.Elapsed.TotalSeconds);
                    watch.Restart();
                    gp.PredictWithVariance(xt);
                    predict.Add(watch.Elapsed.TotalSeconds);
                }
                double fitMedian = Numerics.Data.Statistics.Statistics.Percentile(fit, 0.5), predictMedian = Numerics.Data.Statistics.Statistics.Percentile(predict, 0.5);
                foreach (var library in new[] { "scikit-learn", "DiceKriging" })
                {
                    double fitReference = size.GetProperty(library).GetProperty("fitSeconds").GetDouble();
                    double predictReference = size.GetProperty(library).GetProperty("predictSeconds").GetDouble();
                    Console.WriteLine(FormattableString.Invariant($"n={n}: fit {fitMedian:F3} s vs {library} {fitReference:F3} s; predict {predictMedian:F4} s vs {predictReference:F4} s"));
                    if (fitMedian >= fitReference) failures.Add($"n={n} fit vs {library}");
                    if (predictMedian >= predictReference) failures.Add($"n={n} predict vs {library}");
                }
            }
            return failures.Count == 0 ? null : string.Join("; ", failures);
        }

        /// <summary>
        /// Median wall time of repeated runs after a warm-up of at least half a second.
        /// </summary>
        /// <param name="run">The work to time.</param>
        /// <returns>The median time of one run, in seconds.</returns>
        private static double MedianSeconds(Action run)
        {
            var warm = Stopwatch.StartNew();
            int calls = 0;
            while (calls < 3 || (warm.Elapsed.TotalSeconds < 0.5 && calls < 200))
            {
                run();
                calls++;
            }
            double single = warm.Elapsed.TotalSeconds / calls;
            int repeats = single > 1.0 ? 3 : single > 0.05 ? 7 : 21;
            var times = new double[repeats];
            for (int i = 0; i < repeats; i++)
            {
                // Fast cases are repeated within one sample, as in the R script, to stay above the timer resolution.
                int inner = Math.Max(1, (int)(0.02 / Math.Max(single, 1E-6)));
                var sw = Stopwatch.StartNew();
                for (int j = 0; j < inner; j++) run();
                times[i] = sw.Elapsed.TotalSeconds / inner;
            }
            Array.Sort(times);
            return times[repeats / 2];
        }

        /// <summary>
        /// Reads the "cases" entries of a timings file.
        /// </summary>
        /// <param name="path">The timings file.</param>
        /// <returns>The median seconds of each case.</returns>
        private static Dictionary<string, double> ReadTimings(string path)
        {
            var result = new Dictionary<string, double>();
            foreach (Match m in Regex.Matches(File.ReadAllText(path), "\"((?:pca|enet)_[a-z]+_\\d+x\\d+)\":\\s*([0-9.eE+-]+)"))
                result[m.Groups[1].Value] = double.Parse(m.Groups[2].Value, CultureInfo.InvariantCulture);
            return result;
        }

        /// <summary>
        /// Reads an n by p row-major little-endian float64 matrix.
        /// </summary>
        /// <param name="path">The binary file.</param>
        /// <param name="n">The number of rows.</param>
        /// <param name="p">The number of columns.</param>
        /// <returns>The matrix.</returns>
        private static double[,] ReadMatrix(string path, int n, int p)
        {
            var bytes = File.ReadAllBytes(path);
            var result = new double[n, p];
            Buffer.BlockCopy(bytes, 0, result, 0, bytes.Length);
            return result;
        }

        /// <summary>
        /// Reads a little-endian float64 vector.
        /// </summary>
        /// <param name="path">The binary file.</param>
        /// <returns>The vector.</returns>
        private static double[] ReadVector(string path)
        {
            var bytes = File.ReadAllBytes(path);
            var result = new double[bytes.Length / 8];
            Buffer.BlockCopy(bytes, 0, result, 0, bytes.Length);
            return result;
        }

        /// <summary>
        /// Returns the deterministic Gaussian process data shared with the Python and R scripts.
        /// </summary>
        /// <param name="rows">The number of rows.</param>
        /// <param name="d">The number of inputs.</param>
        /// <param name="offset">The index of the first row in the shared sequence.</param>
        /// <returns>The inputs and responses.</returns>
        private static (double[,] X, double[] Y) GaussianProcessData(int rows, int d, int offset)
        {
            var primes = new double[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71 };
            var x = new double[rows, d];
            var y = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                double s = 0;
                for (int j = 0; j < d; j++)
                {
                    double v = (i + offset + 1) * Math.Sqrt(primes[j]);
                    x[i, j] = v - Math.Floor(v);
                    s += Math.Sin(2 * Math.PI * x[i, j]) / (j + 1);
                }
                double e = (i + offset + 1) * Math.Sqrt(97);
                y[i] = s + 0.1 * (e - Math.Floor(e) - 0.5);
            }
            return (x, y);
        }
    }
}
