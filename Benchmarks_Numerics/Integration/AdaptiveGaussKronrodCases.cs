using System.Globalization;
using System.Text;
using Numerics;
using Numerics.Distributions;
using Numerics.Mathematics;
using Numerics.Mathematics.Integration;
using Numerics.Sampling;

namespace Benchmarks_Numerics.Integration
{
    /// <summary>
    /// One integration problem shared by the accuracy comparison and the runtime benchmark: the integrand,
    /// the range or stratification bins, and the integrator settings of the call site it represents.
    /// </summary>
    public sealed class AdaptiveGaussKronrodCase
    {
        /// <summary>Gets the case name, matched by the mpmath oracle script.</summary>
        public string Name { get; init; } = "";
        /// <summary>Gets the call site or test the case represents.</summary>
        public string Source { get; init; } = "";
        /// <summary>Gets the integrand.</summary>
        public Func<double, double> F { get; init; } = x => 0d;
        /// <summary>Gets the lower limit.</summary>
        public double A { get; init; }
        /// <summary>Gets the upper limit.</summary>
        public double B { get; init; }
        /// <summary>Gets the stratification bins, or null to integrate [A, B].</summary>
        public List<StratificationBin>? Bins { get; init; }
        /// <summary>Gets the relative tolerance.</summary>
        public double Rel { get; init; } = 1E-8;
        /// <summary>Gets the absolute tolerance.</summary>
        public double Abs { get; init; } = 1E-8;
        /// <summary>Gets the minimum depth.</summary>
        public int MinDepth { get; init; }
        /// <summary>Gets the maximum depth.</summary>
        public int MaxDepth { get; init; } = 100;
        /// <summary>Gets the evaluation budget.</summary>
        public int MaxEvaluations { get; init; } = 10_000_000;
        /// <summary>Gets a value indicating whether a recorder is attached, as at the TotalRisk sampler sites.</summary>
        public bool Record { get; init; }

        /// <summary>Gets the weighted sum of the recorded nodes from the last run, when <see cref="Record"/> is set.</summary>
        public double Composite { get; private set; } = double.NaN;

        /// <inheritdoc/>
        public override string ToString() => Name;

        /// <summary>
        /// Runs AdaptiveGaussKronrod with the QAGS method on the case.
        /// </summary>
        /// <returns>The integrator after integration.</returns>
        public AdaptiveGaussKronrod RunNew() => Run(AdaptiveGaussKronrodMethod.Qags);

        /// <summary>
        /// Runs AdaptiveGaussKronrod with its default method, local bisection, on the case.
        /// </summary>
        /// <returns>The integrator after integration.</returns>
        public AdaptiveGaussKronrod RunOld() => Run(AdaptiveGaussKronrodMethod.LocalBisection);

        private AdaptiveGaussKronrod Run(AdaptiveGaussKronrodMethod method)
        {
            var gk = new AdaptiveGaussKronrod(F, A, Bins == null ? B : Math.Max(B, A + 1)) { ReportFailure = false, Method = method };
            Configure(gk);
            if (method == AdaptiveGaussKronrodMethod.LocalBisection)
            {
                // Local bisection keeps the 1E-15 tolerance floor, so an epsabs-only request runs at that floor.
                gk.RelativeTolerance = Math.Max(Rel, 1E-15);
                gk.AbsoluteTolerance = Math.Max(Abs, 1E-15);
            }
            gk.MinDepth = MinDepth;
            gk.MaxDepth = MaxDepth;
            Composite = 0;
            if (Record) gk.Recorder = (x, w, f) => Composite += w * f;
            if (Bins == null) gk.Integrate(); else gk.Integrate(Bins);
            return gk;
        }

        private void Configure(Integrator gk)
        {
            gk.RelativeTolerance = Rel;
            gk.AbsoluteTolerance = Abs;
            gk.MaxFunctionEvaluations = MaxEvaluations;
        }
    }

    /// <summary>
    /// The comparison cases: the Test_Numerics integrands, the internal Numerics callers, representative
    /// TotalRisk integrands and classic hard cases.
    /// </summary>
    public static class AdaptiveGaussKronrodCases
    {
        private static double Phi(double z) => Normal.StandardCDF(z);

        // Hosking GEV quantile and the TotalRisk-style hazard used by the representative sampler integrands.
        private const double GevXi = 100, GevAlpha = 30, GevKappa = -0.1;
        private static double GevQuantile(double p) => GevXi + GevAlpha / GevKappa * (1 - Math.Pow(-Math.Log(p), GevKappa));
        private static double Fragility(double x) => Phi(Math.Log(x / 250) / 0.25);
        private static double Damage(double x) => 1E6 * Math.Max(0, x - 150) / (x + 50);

        /// <summary>
        /// Builds the TotalRisk 50-bin hazard stratification in probability space (RiskAnalysis.BuildHazardBins).
        /// </summary>
        /// <returns>The bins.</returns>
        public static List<StratificationBin> HazardBins()
        {
            var gev = new GeneralizedExtremeValue(GevXi, GevAlpha, GevKappa);
            var bins = Stratify.XValues(new StratificationOptions(gev.InverseCDF(1E-16), gev.InverseCDF(1 - 1E-16), 50), true);
            return Stratify.XToProbability(bins, gev.CDF, false);
        }

        /// <summary>
        /// Returns all cases.
        /// </summary>
        /// <returns>The cases.</returns>
        public static List<AdaptiveGaussKronrodCase> All()
        {
            var list = new List<AdaptiveGaussKronrodCase>();
            var gamma = new GammaDistribution(10, 5);
            var ln = new LnNormal(10, 2);
            var gev = new GeneralizedExtremeValue(GevXi, GevAlpha, GevKappa);

            // Test_Numerics Test_AdaptiveGaussKronrod and Test_AdaptiveGaussKronrodRecorder.
            list.Add(new() { Name = "t_x3", Source = "Test_AdaptiveGaussKronrod.Test_FX3", F = x => Math.Pow(x, 3), A = 0, B = 1 });
            list.Add(new() { Name = "t_cos", Source = "Test_AdaptiveGaussKronrod.Test_Cosine", F = Math.Cos, A = -1, B = 1 });
            list.Add(new() { Name = "t_sin", Source = "Test_AdaptiveGaussKronrod.Test_Sine", F = Math.Sin, A = 0, B = 1 });
            list.Add(new() { Name = "t_fxx", Source = "Test_AdaptiveGaussKronrod.Test_FXX", F = x => 0.5 + 24 * x + 3 * x * x, A = 0, B = 2 });
            list.Add(new() { Name = "t_fxxx", Source = "Test_AdaptiveGaussKronrod.Test_FXXX", F = x => 0.5 + 24 * x + 3 * x * x + 8 * x * x * x, A = 0, B = 2 });
            list.Add(new() { Name = "t_gamma_mean", Source = "Test_AdaptiveGaussKronrod.Test_Gamma", F = x => x * gamma.PDF(x), A = gamma.InverseCDF(1E-16), B = gamma.InverseCDF(1 - 1E-16) });
            list.Add(new() { Name = "t_cvar", Source = "Test_AdaptiveGaussKronrod.Test_CVaR", F = x => x * ln.PDF(x), A = ln.InverseCDF(0.99), B = ln.InverseCDF(1 - 1E-16) });
            list.Add(new() { Name = "t_expcos", Source = "Test_AdaptiveGaussKronrod.Test_SmoothIntegrand", F = x => Math.Exp(-x) * Math.Cos(3 * x), A = 0, B = 2 });
            list.Add(new() { Name = "t_step_d3", Source = "Test_AdaptiveGaussKronrod.Test_ForcedAcceptance (step)", F = x => x < 0.3 ? 0.0 : 1.0, A = 0, B = 1, MaxDepth = 3 });
            list.Add(new() { Name = "t_inv_x", Source = "Test_AdaptiveGaussKronrod.Test_ForcedAcceptance (1/x, divergent)", F = x => 1.0 / x, A = 0, B = 1 });
            list.Add(new() { Name = "t_inv_sqrt_abs", Source = "Test_AdaptiveGaussKronrod.Test_NaNIntegrand (1/sqrt|x|)", F = x => 1.0 / Math.Sqrt(Math.Abs(x)), A = -1, B = 1 });
            list.Add(new() { Name = "t_sinc", Source = "Test_AdaptiveGaussKronrod.Test_RemovableNaN", F = x => Math.Sin(x - 0.5) / (x - 0.5), A = 0, B = 1 });
            list.Add(new() { Name = "t_peak_rec", Source = "Test_AdaptiveGaussKronrodRecorder (MinDepth 2, recorder)", F = SharpPeak, A = 0, B = 1, Rel = 1E-10, MinDepth = 2, Record = true });
            list.Add(new() { Name = "t_peak_bins", Source = "Test_AdaptiveGaussKronrodRecorder.Test_Recorder_StratifiedBins", F = SharpPeak, Bins = new() { new(0, 0.25), new(0.25, 0.6), new(0.6, 1) }, A = 0, B = 1, Rel = 1E-10, Record = true });
            list.Add(new() { Name = "t_peak", Source = "Test_AdaptiveGaussKronrodRecorder.Test_Recorder_OffIsByteIdentical", F = SharpPeak, A = 0, B = 1, Rel = 1E-10 });
            list.Add(new() { Name = "t_x2_d2", Source = "Test_AdaptiveGaussKronrodRecorder.Test_Recorder_IsSnapshot", F = x => x * x, A = 0, B = 1, MinDepth = 2, Record = true });

            // Internal Numerics callers.
            double gevMean = 100;
            list.Add(new() { Name = "n_moment_gev_mean", Source = "UnivariateDistributionBase.cs:521", F = x => x * gev.PDF(x), A = gev.InverseCDF(1E-16), B = gev.InverseCDF(1 - 1E-16) });
            list.Add(new() { Name = "n_moment_gev_m2", Source = "UnivariateDistributionBase.cs:525", F = x => Math.Pow(x - gevMean, 2) * gev.PDF(x), A = gev.InverseCDF(1E-16), B = gev.InverseCDF(1 - 1E-16) });
            list.Add(new() { Name = "n_moment_gamma_m4", Source = "UnivariateDistributionBase.cs:534", F = x => Math.Pow((x - 50) / Math.Sqrt(500), 4) * gamma.PDF(x), A = gamma.InverseCDF(1E-16), B = gamma.InverseCDF(1 - 1E-16) });
            foreach (int order in new[] { 0, 2, 4 })
            {
                int k = order;
                list.Add(new() { Name = "n_dmi_normal_k" + k, Source = "DistributionMomentIntegration.cs:50", F = t => MomentTransform(t, k, x => -0.5 * x * x - 0.5 * Math.Log(2 * Math.PI)), A = 0, B = 1, Rel = 1E-8, Abs = 1E-10, MaxEvaluations = 200000 });
            }
            double lgT = Numerics.Mathematics.SpecialFunctions.Gamma.LogGamma(5.5) - Numerics.Mathematics.SpecialFunctions.Gamma.LogGamma(5) - 0.5 * Math.Log(10 * Math.PI);
            list.Add(new() { Name = "n_dmi_t10_k4", Source = "DistributionMomentIntegration.cs:50", F = t => MomentTransform(t, 4, x => lgT - 5.5 * Math.Log(1 + x * x / 10)), A = 0, B = 1, Rel = 1E-8, Abs = 1E-10, MaxEvaluations = 200000 });
            var chi = new GammaDistribution(1, 2.5);
            list.Add(new() { Name = "n_mvt_radial", Source = "MultivariateNormal.cs:2206", F = u => Phi(1.3 * Math.Sqrt(chi.InverseCDF(u) / 2.5)), A = 0, B = 1, Rel = 1E-11, Abs = 1E-13 });
            list.Add(new() { Name = "n_kappa4_mean", Source = "KappaFour.cs:1090", F = u => KappaTransform(u, lp => -Math.Log(-lp)), A = 0, B = 1, Rel = 1E-10, Abs = 1E-12 });
            list.Add(new() { Name = "n_kappa4_m2", Source = "KappaFour.cs:1090", F = u => KappaTransform(u, lp => Math.Pow(Math.Log(-lp), 2)), A = 0, B = 1, Rel = 1E-10, Abs = 1E-12 });
            list.Add(new() { Name = "n_nct", Source = "NoncentralT.cs:610", F = x => x * x / 9 * Math.Exp(-0.5 * (x - 3) * (x - 3)), A = 0, B = 12, Rel = 1E-14, Abs = 1E-15 });
            double betaFunction = Numerics.Mathematics.SpecialFunctions.Beta.Function(2.2, 3.8);
            list.Add(new() { Name = "n_pert_mean", Source = "PertPercentile.cs:533", F = x => x * Math.Pow(x, 1.2) * Math.Pow(10 - x, 2.8) / (betaFunction * 1E5), A = 0, B = 10, Rel = 1E-12, Abs = 1E-15 });
            list.Add(new() { Name = "n_pertz_mean", Source = "PertPercentileZ.cs:522", F = z => Phi(z) * Math.Pow(z + 2, 1.2) * Math.Pow(3 - z, 2.8) / (betaFunction * Math.Pow(5, 5)), A = -2, B = 3, Rel = 1E-12, Abs = 1E-15 });
            double zt = Phi(2) - Phi(-1);
            list.Add(new() { Name = "n_truncated_m2", Source = "TruncatedDistribution.cs:382", F = t => t * t * Normal.StandardPDF(0.5 + 1.5 * t) / zt * 1.5, A = -1, B = 1, Rel = 1E-12, Abs = 1E-15 });
            double w = 8 / Math.Sqrt(10);
            list.Add(new() { Name = "n_vonmises_mid", Source = "VonMises.cs:452", F = t => VonMisesKernel(t), A = 0.5 - w, B = 0.5 });
            list.Add(new() { Name = "n_vonmises_tail", Source = "VonMises.cs:452", F = t => VonMisesKernel(t), A = -Math.PI, B = 0.5 - w });
            var thresholds = new[] { Normal.StandardZ(0.01), Normal.StandardZ(0.02), Normal.StandardZ(0.05) };
            list.Add(new() { Name = "n_union_single_factor", Source = "Probability.cs:810", F = u => Union(u, thresholds, 0.5), A = 1E-16, B = 1 - 1E-16, Rel = 1E-8, Abs = 1E-15, MinDepth = 2 });

            // Representative TotalRisk integrands (read-only scan of C:\GIT\RMC-TotalRisk).
            var bins = HazardBins();
            double lo = bins[0].LowerBound, hi = bins[bins.Count - 1].UpperBound;
            list.Add(new() { Name = "tr_probe_pf", Source = "RiskAnalysis.cs:5554 (failure-probability probe)", F = p => Fragility(GevQuantile(p)), Bins = bins, A = lo, B = hi, Rel = 1E-8, MinDepth = 2, MaxEvaluations = 1_000_000 });
            list.Add(new() { Name = "tr_risk_mean", Source = "RiskAnalysis.cs:5037 (mean pass, recorder)", F = p => { double x = GevQuantile(p); return Fragility(x) * Damage(x); }, Bins = bins, A = lo, B = hi, Rel = 1E-8, MinDepth = 2, MaxEvaluations = 1_000_000, Record = true });
            list.Add(new() { Name = "tr_risk_ensemble", Source = "RiskAnalysis.cs:5037 (ensemble pass, recorder)", F = p => { double x = GevQuantile(p); return Fragility(x) * Damage(x); }, Bins = bins, A = lo, B = hi, Rel = 1E-4, MinDepth = 0, MaxEvaluations = 1_000_000, Record = true });
            double zf = Normal.StandardZ(1E-16), zc = Normal.StandardZ(1 - 1E-16);
            list.Add(new() { Name = "tr_conditional", Source = "SampledComponent.cs:1992 (conditional sweep, recorder)", F = z => Normal.StandardPDF(z) * Phi((z - 1.5) / 0.3), A = zf, B = zc, Rel = 1E-8, Abs = 1E-15, MinDepth = 2, MaxEvaluations = 21 * 50, Record = true });

            // Classic hard cases.
            list.Add(new() { Name = "h_inv_sqrt", Source = "endpoint singularity x^-0.5", F = x => 1 / Math.Sqrt(x), A = 0, B = 1 });
            list.Add(new() { Name = "h_log", Source = "endpoint singularity log x", F = Math.Log, A = 0, B = 1 });
            list.Add(new() { Name = "h_x_m09", Source = "endpoint singularity x^-0.9", F = x => Math.Pow(x, -0.9), A = 0, B = 1 });
            list.Add(new() { Name = "h_log_over_sqrt", Source = "endpoint singularity log(x)/sqrt(x)", F = x => Math.Log(x) / Math.Sqrt(x), A = 0, B = 1 });
            list.Add(new() { Name = "h_lorentz_peak", Source = "peak 1/((x-0.3)^2+1e-4)", F = x => 1 / ((x - 0.3) * (x - 0.3) + 1E-4), A = 0, B = 1 });
            list.Add(new() { Name = "h_cos100", Source = "oscillatory cos(100x)", F = x => Math.Cos(100 * x), A = 0, B = 1 });
            list.Add(new() { Name = "h_xsin50", Source = "oscillatory x sin(50x) on [0, 2pi]", F = x => x * Math.Sin(50 * x), A = 0, B = 2 * Math.PI });
            list.Add(new() { Name = "h_abs_kink", Source = "interior kink |x-1/3|", F = x => Math.Abs(x - 1.0 / 3), A = 0, B = 1 });
            list.Add(new() { Name = "h_exp", Source = "smooth exp(x)", F = Math.Exp, A = 0, B = 1 });
            list.Add(new() { Name = "h_runge", Source = "smooth 1/(1+x^2)", F = x => 1 / (1 + x * x), A = 0, B = 1 });
            list.Add(new() { Name = "h_inv_sqrt_tight", Source = "x^-0.5, epsabs-only 1e-14", F = x => 1 / Math.Sqrt(x), A = 0, B = 1, Rel = 0, Abs = 1E-14 });
            list.Add(new() { Name = "h_gauss_wide", Source = "exp(-x^2) on [-10, 10] (pynumerics callback test)", F = x => Math.Exp(-x * x), A = -10, B = 10 });
            return list;
        }

        /// <summary>The Test_AdaptiveGaussKronrodRecorder sharp peak.</summary>
        public static double SharpPeak(double x)
        {
            double t = (x - 0.3) * 60d;
            return Math.Exp(-t * t) + 0.1 * x;
        }

        // DistributionMomentIntegration: x = t / (1 - t), center 0 and scale 1, offset 0.
        private static double MomentTransform(double t, int order, Func<double, double> logDensity)
        {
            double x = t / (1 - t);
            double log = logDensity(x);
            if (double.IsNegativeInfinity(log)) return 0;
            if (order != 0)
            {
                if (x == 0) return 0;
                log += order * Math.Log(Math.Abs(x));
            }
            return Math.Exp(log - 2 * Tools.Log1p(-t));
        }

        // KappaFour.IntegrateProbability: p = u^8 / 2 on each half.
        private static double KappaTransform(double u, Func<double, double> integrand)
        {
            double logU = Math.Log(u);
            double logTail = 8d * logU - Math.Log(2d);
            double jacobian = Math.Exp(Math.Log(4d) + 7d * logU);
            return integrand(logTail) * jacobian + integrand(Tools.Log1p(-Math.Exp(logTail))) * jacobian;
        }

        // VonMises.CDF kernel with mu = 0.5, kappa = 10 and c = kappa.
        private static double VonMisesKernel(double t)
        {
            double s = Math.Sin(0.5 * (t - 0.5));
            return Math.Exp(10 - 2d * 10 * s * s);
        }

        // Probability.UnionSingleFactor integrand.
        private static double Union(double u, double[] thresholds, double rho)
        {
            double z = Normal.StandardZ(u);
            double logSurvival = 0;
            for (int i = 0; i < thresholds.Length; i++)
                logSurvival += Tools.Log1p(-Phi((thresholds[i] - Math.Sqrt(rho) * z) / Math.Sqrt(1 - rho)));
            return -Tools.Expm1(logSurvival);
        }

        /// <summary>
        /// Writes the comparison table (local bisection and QAGS results, evaluations and status per case) and the ranges the
        /// oracle script integrates over.
        /// </summary>
        /// <param name="path">The output CSV path.</param>
        public static void WriteComparison(string path)
        {
            var sb = new StringBuilder();
            sb.AppendLine("name;source;rel;abs;old_result;old_evals;old_status;old_se;new_result;new_evals;new_status;new_se;new_iterations;ranges;old_composite;new_composite");
            foreach (var c in All())
            {
                var o = c.RunOld();
                double oldComposite = c.Composite;
                var n = c.RunNew();
                double newComposite = c.Composite;
                string ranges = c.Bins == null
                    ? R(c.A) + ":" + R(c.B)
                    : string.Join(",", c.Bins.Select(bn => R(bn.LowerBound) + ":" + R(bn.UpperBound)));
                sb.AppendLine(string.Join(";", c.Name, c.Source, R(c.Rel), R(c.Abs), R(o.Result), o.FunctionEvaluations, o.Status,
                    R(o.StandardError), R(n.Result), n.FunctionEvaluations, n.Status, R(n.StandardError), n.Iterations, ranges, c.Record ? R(oldComposite) : "", c.Record ? R(newComposite) : ""));
            }
            File.WriteAllText(path, sb.ToString());
        }

        private static string R(double v) => v.ToString("R", CultureInfo.InvariantCulture);
    }
}
