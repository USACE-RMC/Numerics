using Numerics.Data.Statistics;
using Numerics.Distributions;
using Numerics.Mathematics.LinearAlgebra;
using Numerics.Mathematics.Optimization;
using Numerics.Sampling;
using System;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
#if NET8_0_OR_GREATER
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
#endif

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Exact Gaussian process regression with stationary covariance functions.
    /// </summary>
    /// <remarks>
    /// <para>
    ///     <b> Authors: </b>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </para>
    /// <para>
    /// <b> Description: </b>
    /// </para>
    /// <para>
    /// The response is modeled as <c>y = m(x) + f(x) + ε</c>, where f is a zero-mean Gaussian process with
    /// covariance <c>σ_f² ρ(r)</c> (see <see cref="GaussianProcessKernel"/>), ε is independent Gaussian noise with
    /// variance <c>σ_n²</c> (the nugget), and m is a zero, constant, or linear mean whose coefficients are
    /// estimated by generalized least squares. The hyperparameters <c>θ = (ln σ_f², ln ℓ_1, ..., ln ℓ_q, ln σ_n²)</c>,
    /// with q = 1 for a shared length scale or q = d with automatic relevance determination, are chosen by
    /// maximizing the log marginal likelihood
    /// </para>
    /// <para>
    /// <c>ln p(y | θ) = -½ (y - Fβ)ᵀ K⁻¹ (y - Fβ) - ½ ln|K| - (n/2) ln 2π,</c>
    /// </para>
    /// <para>
    /// with <c>K = σ_f² R + σ_n² I</c> and β the generalized least squares coefficients, using
    /// <see cref="LBFGSB"/> with the analytic gradient <c>∂/∂θ_k = ½ tr((ααᵀ - K⁻¹) ∂K/∂θ_k)</c>, where
    /// <c>α = K⁻¹(y - Fβ)</c>. Because β maximizes the likelihood for fixed θ, the gradient needs no term for β.
    /// Two data-based starts are always run; optional restarts drawn uniformly from a box of plausible values use a
    /// seeded Mersenne Twister.
    /// </para>
    /// <para>
    /// The covariance matrix is factored with a blocked, vectorized Cholesky factorization. When it is not
    /// numerically positive-definite, the factorization falls back to
    /// <see cref="CholeskyDecomposition.TryCholesky(Matrix, out CholeskyDecomposition, out double)"/>, and the
    /// diagonal jitter it adds is reported in <see cref="Jitter"/> and treated as additional noise. Training costs
    /// O(n³) time and O(n²) memory per likelihood evaluation, which suits data sets up to a few thousand points.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <list type="bullet">
    /// <item><description>
    /// Rasmussen, C. E., and Williams, C. K. I. (2006). Gaussian Processes for Machine Learning. MIT Press.
    /// Algorithm 2.1, equations 2.29-2.30 (mean basis functions), 5.9 (gradient) and 5.12 (leave-one-out).
    /// </description></item>
    /// <item><description>
    /// Dubrule, O. (1983). Cross validation of kriging in a unique neighborhood. Mathematical Geology, 15(6), 687-699.
    /// </description></item>
    /// </list>
    /// </remarks>
    public class GaussianProcessRegression
    {

        #region Construction

        /// <summary>
        /// Creates a Gaussian process regression model for a single feature.
        /// </summary>
        /// <param name="x">The training values of the feature.</param>
        /// <param name="y">The training response values.</param>
        /// <param name="kernel">The covariance function. Default = Matérn 5/2.</param>
        /// <param name="mean">The mean function. Default = constant.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the data are empty, mismatched, or not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the kernel or mean function is not a defined value.</exception>
        public GaussianProcessRegression(double[] x, double[] y, GaussianProcessKernel kernel = GaussianProcessKernel.Matern52, GaussianProcessMean mean = GaussianProcessMean.Constant)
            : this(ToColumn(x), y, kernel, mean) { }

        /// <summary>
        /// Creates a Gaussian process regression model.
        /// </summary>
        /// <param name="x">The training features, with one row per observation.</param>
        /// <param name="y">The training response values.</param>
        /// <param name="kernel">The covariance function. Default = Matérn 5/2.</param>
        /// <param name="mean">The mean function. Default = constant.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the data are empty, mismatched, or not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the kernel or mean function is not a defined value.</exception>
        public GaussianProcessRegression(double[,] x, double[] y, GaussianProcessKernel kernel = GaussianProcessKernel.Matern52, GaussianProcessMean mean = GaussianProcessMean.Constant)
            : this(Flatten(x), x?.GetLength(0) ?? 0, x?.GetLength(1) ?? 0, y, ToCovarianceForm(kernel), mean, null, false) { }

        /// <summary>
        /// Creates a Gaussian process regression model.
        /// </summary>
        /// <param name="x">The training features, with one row per observation.</param>
        /// <param name="y">The training response values.</param>
        /// <param name="kernel">The covariance function. Default = Matérn 5/2.</param>
        /// <param name="mean">The mean function. Default = constant.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the data are empty, mismatched, or not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the kernel or mean function is not a defined value.</exception>
        public GaussianProcessRegression(Matrix x, Vector y, GaussianProcessKernel kernel = GaussianProcessKernel.Matern52, GaussianProcessMean mean = GaussianProcessMean.Constant)
            : this(x?.Array!, y?.Array!, kernel, mean) { }

        /// <summary>
        /// Creates a model from flat row-major features.
        /// </summary>
        /// <param name="x">The n-by-d row-major features, which are copied.</param>
        /// <param name="n">The number of observations.</param>
        /// <param name="d">The number of features.</param>
        /// <param name="y">The response values, which are copied.</param>
        /// <param name="form">The covariance form.</param>
        /// <param name="mean">The mean function.</param>
        /// <param name="distance">Optional. A distance between feature vectors that replaces the scaled Euclidean distance.</param>
        /// <param name="exactNugget">
        /// When true, the nugget also enters the covariance between a prediction point and a training point at zero
        /// distance, so predictions reproduce the data there, as in classical kriging.
        /// </param>
        internal GaussianProcessRegression(double[] x, int n, int d, double[] y, CovarianceForm form, GaussianProcessMean mean,
            Func<double[], double[], double>? distance, bool exactNugget)
        {
            if (y == null) throw new ArgumentNullException(nameof(y));
            if (!Enum.IsDefined(typeof(GaussianProcessMean), mean)) throw new ArgumentOutOfRangeException(nameof(mean), "The mean function is not a defined value.");
            if (n < 1 || d < 1) throw new ArgumentException("There must be at least one observation and one feature.", nameof(x));
            if (y.Length != n) throw new ArgumentException("The y vector must be the same length as the number of rows in x.", nameof(y));
            for (int i = 0; i < x.Length; i++)
                if (!Tools.IsFinite(x[i])) throw new ArgumentException("The features must be finite.", nameof(x));
            for (int i = 0; i < n; i++)
                if (!Tools.IsFinite(y[i])) throw new ArgumentException("The response values must be finite.", nameof(y));

            _n = n;
            _d = d;
            _x = (double[])x.Clone();
            _y = (double[])y.Clone();
            _form = form;
            _distance = distance;
            _exactNugget = exactNugget;
            Kernel = form == CovarianceForm.Spherical ? GaussianProcessKernel.Exponential : (GaussianProcessKernel)form;
            MeanFunction = mean;
            _p = mean == GaussianProcessMean.Zero ? 0 : mean == GaussianProcessMean.Constant ? 1 : d + 1;
            _basis = new double[_p * n];
            for (int i = 0; i < n; i++)
            {
                if (_p > 0) _basis[i] = 1d;
                if (_p > 1)
                    for (int j = 0; j < d; j++) _basis[(j + 1) * n + i] = x[i * d + j];
            }
            if (distance != null)
            {
                _rows = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    _rows[i] = new double[d];
                    Array.Copy(x, i * d, _rows[i], 0, d);
                }
            }
        }

        /// <summary>Converts a public kernel to its covariance form.</summary>
        /// <param name="kernel">The kernel.</param>
        /// <returns>The covariance form of the kernel.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the kernel is not a defined value.</exception>
        private static CovarianceForm ToCovarianceForm(GaussianProcessKernel kernel)
        {
            if (!Enum.IsDefined(typeof(GaussianProcessKernel), kernel))
                throw new ArgumentOutOfRangeException(nameof(kernel), "The kernel is not a defined value.");
            return (CovarianceForm)kernel;
        }

        /// <summary>Copies a vector into a single-column row-major array.</summary>
        /// <param name="x">The values.</param>
        /// <returns>The n-by-1 array.</returns>
        private static double[,] ToColumn(double[] x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var result = new double[x.Length, 1];
            for (int i = 0; i < x.Length; i++) result[i, 0] = x[i];
            return result;
        }

        /// <summary>Copies a two-dimensional array into a flat row-major array.</summary>
        /// <param name="x">The array.</param>
        /// <returns>The flat copy.</returns>
        private static double[] Flatten(double[,] x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            int n = x.GetLength(0), d = x.GetLength(1);
            var result = new double[n * d];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < d; j++) result[i * d + j] = x[i, j];
            return result;
        }

        #endregion

        #region Members

        private readonly int _n, _d, _p;
        private readonly double[] _x, _y, _basis;
        private readonly double[][]? _rows;
        private readonly CovarianceForm _form;
        private readonly Func<double[], double[], double>? _distance;
        private readonly bool _exactNugget;
        private Workspace? _fit;
        private const double Ln10 = 2.302585092994046;

        /// <summary>The relative-reduction stopping tolerance of the optimizer, 10⁷ times the double precision epsilon 2⁻⁵², the reference L-BFGS-B default.</summary>
        private const double LikelihoodRelativeTolerance = 2.220446049250313E-9;
        private double[] _ell = Array.Empty<double>();
        private double _sf2, _sn2;

        /// <summary>
        /// Gets the covariance function.
        /// </summary>
        public GaussianProcessKernel Kernel { get; private set; }

        /// <summary>
        /// Gets the mean function.
        /// </summary>
        public GaussianProcessMean MeanFunction { get; private set; }

        /// <summary>
        /// Gets the number of training observations.
        /// </summary>
        public int SampleSize => _n;

        /// <summary>
        /// Gets the number of features.
        /// </summary>
        public int NumberOfFeatures => _d;

        /// <summary>
        /// Gets or sets whether each feature has its own length scale (automatic relevance determination).
        /// Default = true. When false, one length scale is shared by all features.
        /// </summary>
        public bool UseARD { get; set; } = true;

        /// <summary>
        /// Gets or sets the number of optimizer starts drawn uniformly from a box of plausible values, in addition
        /// to the two data-based starts. Default = 0.
        /// </summary>
        public int NumberOfRestarts { get; set; } = 0;

        /// <summary>
        /// Gets or sets the seed of the pseudo random number generator for the restarts. Default = 12345.
        /// </summary>
        public int PRNGSeed { get; set; } = 12345;

        /// <summary>
        /// Gets or sets the maximum number of L-BFGS-B iterations per start. Default = 200.
        /// </summary>
        public int MaxIterations { get; set; } = 200;

        /// <summary>
        /// Gets or sets the convergence tolerance on the largest projected gradient component of the
        /// log marginal likelihood with respect to the log hyperparameters. Default = 1E-5.
        /// </summary>
        public double Tolerance { get; set; } = 1E-5;

        /// <summary>
        /// Gets the signal variance σ_f².
        /// </summary>
        public double SignalVariance => _sf2;

        /// <summary>
        /// Gets the length scales: one per feature with automatic relevance determination, otherwise one shared value.
        /// </summary>
        public double[] LengthScales => (double[])_ell.Clone();

        /// <summary>
        /// Gets the noise (nugget) variance σ_n².
        /// </summary>
        public double NoiseVariance => _sn2;

        /// <summary>
        /// Gets the diagonal jitter added to make the covariance matrix numerically positive-definite; zero when none
        /// was needed. It acts as additional noise variance.
        /// </summary>
        public double Jitter => _fit?.Jitter ?? 0d;

        /// <summary>
        /// Gets the generalized least squares coefficients of the mean function: empty for a zero mean,
        /// [β₀] for a constant mean, and [β₀, β₁, ..., β_d] for a linear mean.
        /// </summary>
        public double[] MeanCoefficients => _fit == null ? Array.Empty<double>() : (double[])_fit.Beta.Clone();

        /// <summary>
        /// Gets the log marginal likelihood at the current hyperparameters.
        /// </summary>
        public double LogMarginalLikelihood { get; private set; } = double.NaN;

        /// <summary>
        /// Gets the status of the best optimizer run of the last <see cref="Train()"/>.
        /// </summary>
        public OptimizationStatus Status { get; private set; } = OptimizationStatus.None;

        /// <summary>
        /// Gets whether hyperparameters have been set by training or by <see cref="SetHyperparameters"/>.
        /// </summary>
        public bool IsTrained => _fit != null;

        #endregion

        #region Training

        /// <summary>
        /// Estimates all hyperparameters, including the noise variance, by maximizing the log marginal likelihood.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when <see cref="NumberOfRestarts"/> is negative, or no optimizer start produces a finite likelihood.</exception>
        /// <remarks>
        /// <para>
        /// The covariance is written <c>σ²(R + gI)</c> with the noise-to-signal ratio <c>g = σ_n²/σ_f²</c>, and σ² is
        /// profiled out at its maximum likelihood value <c>σ̂² = (y - Fβ)ᵀ(R + gI)⁻¹(y - Fβ) / n</c>, which gives the
        /// concentrated log likelihood <c>-(n/2) ln σ̂² - ½ ln|R + gI| - (n/2)(1 + ln 2π)</c>. Its maximum is the
        /// maximum of the full likelihood, without the poorly scaled direction along which σ_f² and the length scales
        /// trade off.
        /// </para>
        /// <para>
        /// Each start runs <see cref="LBFGSB"/> with the analytic gradient on the negative concentrated log likelihood
        /// in <c>(ln ℓ_1, ..., ln ℓ_q, ln g)</c>, stopping at the projected gradient <see cref="Tolerance"/> or a relative
        /// reduction of 10⁷ machine epsilon, the defaults of the reference code. The likelihood usually has several
        /// local maxima that differ in which features are judged irrelevant, so two data-based starts are always run:
        /// <c>ℓ_j = range_j / 2</c> and <c>ℓ_j = range_j</c>, both with <c>g = 0.01</c>, where range_j is the range of
        /// feature j; <see cref="NumberOfRestarts"/> adds random starts and the best result is kept. The bounds are
        /// <c>[10⁻³ range_j, 10³ range_j]</c> for ℓ_j and <c>[10⁻¹⁰, 10²]</c> for g. Trial points whose covariance
        /// matrix cannot be factored are rejected by the line search.
        /// </para>
        /// </remarks>
        public void Train()
        {
            TrainCore(true, 0d);
        }

        /// <summary>
        /// Estimates the signal variance and length scales by maximizing the log marginal likelihood, holding the
        /// noise variance fixed.
        /// </summary>
        /// <param name="noiseVariance">The fixed noise (nugget) variance. Zero gives an interpolating model.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the noise variance is negative or not finite.</exception>
        /// <exception cref="InvalidOperationException">Thrown when <see cref="NumberOfRestarts"/> is negative, or no optimizer start produces a finite likelihood.</exception>
        /// <remarks>
        /// With zero noise, σ_f² is profiled out as in <see cref="Train()"/>. With a positive noise, the optimizer works
        /// on <c>(ln σ_f², ln ℓ_1, ..., ln ℓ_q)</c>, starting from <c>σ_f² = s²</c>, the sample variance of y (its mean
        /// square for a zero mean), within <c>[10⁻⁶ s², 10⁶ s²]</c>; the length scales start and are bounded as in
        /// <see cref="Train()"/>.
        /// </remarks>
        public void Train(double noiseVariance)
        {
            if (!Tools.IsFinite(noiseVariance) || noiseVariance < 0d)
                throw new ArgumentOutOfRangeException(nameof(noiseVariance), "The noise variance must be non-negative and finite.");
            TrainCore(false, noiseVariance);
        }

        /// <summary>
        /// Sets the hyperparameters without optimization and conditions the model on the training data.
        /// </summary>
        /// <param name="signalVariance">The signal variance σ_f², positive.</param>
        /// <param name="lengthScales">
        /// The length scales: one per feature, which enables <see cref="UseARD"/>, or a single shared value,
        /// which disables it.
        /// </param>
        /// <param name="noiseVariance">The noise (nugget) variance σ_n², non-negative.</param>
        /// <exception cref="ArgumentNullException">Thrown when the length scales are null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a hyperparameter is outside its domain.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the covariance matrix cannot be factored even with jitter, or the mean basis is rank deficient.</exception>
        public void SetHyperparameters(double signalVariance, double[] lengthScales, double noiseVariance)
        {
            if (lengthScales == null) throw new ArgumentNullException(nameof(lengthScales));
            if (!Tools.IsFinite(signalVariance) || signalVariance <= 0d)
                throw new ArgumentOutOfRangeException(nameof(signalVariance), "The signal variance must be positive and finite.");
            if (!Tools.IsFinite(noiseVariance) || noiseVariance < 0d)
                throw new ArgumentOutOfRangeException(nameof(noiseVariance), "The noise variance must be non-negative and finite.");
            if (lengthScales.Length != 1 && lengthScales.Length != _d)
                throw new ArgumentOutOfRangeException(nameof(lengthScales), "There must be one length scale or one per feature.");
            if (_distance != null && lengthScales.Length != 1)
                throw new ArgumentOutOfRangeException(nameof(lengthScales), "A model with a custom distance has one length scale.");
            for (int j = 0; j < lengthScales.Length; j++)
                if (!Tools.IsFinite(lengthScales[j]) || lengthScales[j] <= 0d)
                    throw new ArgumentOutOfRangeException(nameof(lengthScales), "The length scales must be positive and finite.");
            UseARD = lengthScales.Length == _d && _d > 1;
            _sf2 = signalVariance;
            _ell = (double[])lengthScales.Clone();
            if (!UseARD && _ell.Length != 1) _ell = new[] { _ell[0] };
            _sn2 = noiseVariance;
            Condition(_fit ?? new Workspace(_n, _d, _p));
            Status = OptimizationStatus.None;
        }

        /// <summary>
        /// Evaluates the log marginal likelihood, and optionally its gradient, at the given log hyperparameters
        /// without changing the model.
        /// </summary>
        /// <param name="theta">
        /// The log hyperparameters <c>(ln σ_f², ln ℓ_1, ..., ln ℓ_q, ln σ_n²)</c>, with q = d when <see cref="UseARD"/>
        /// is true and q = 1 otherwise.
        /// </param>
        /// <param name="gradient">Optional. An array of the same length that receives the gradient.</param>
        /// <returns>The log marginal likelihood, or NaN when the covariance matrix is not numerically positive-definite.</returns>
        /// <exception cref="ArgumentNullException">Thrown when theta is null.</exception>
        /// <exception cref="ArgumentException">Thrown when an array has the wrong length.</exception>
        /// <exception cref="InvalidOperationException">Thrown for a model with a custom distance.</exception>
        public double EvaluateLogMarginalLikelihood(double[] theta, double[]? gradient = null)
        {
            if (theta == null) throw new ArgumentNullException(nameof(theta));
            if (_distance != null) throw new InvalidOperationException("The likelihood gradient requires the Euclidean distance.");
            int q = UseARD ? _d : 1;
            if (theta.Length != q + 2) throw new ArgumentException("The hyperparameter vector has the wrong length.", nameof(theta));
            if (gradient != null && gradient.Length != q + 2) throw new ArgumentException("The gradient must have the same length as theta.", nameof(gradient));
            var w = new Workspace(_n, _d, _p);
            var hp = Hyperparameters.FromTheta(theta, q);
            if (!Evaluate(w, hp, false)) return double.NaN;
            if (gradient != null) Gradient(w, hp, gradient);
            return w.LogLikelihood;
        }

        /// <summary>
        /// Runs the likelihood maximization.
        /// </summary>
        /// <param name="estimateNoise">Whether the noise variance is estimated.</param>
        /// <param name="fixedNoise">The noise variance when it is held fixed.</param>
        private void TrainCore(bool estimateNoise, double fixedNoise)
        {
            if (_distance != null) throw new InvalidOperationException("Training requires the Euclidean distance.");
            if (NumberOfRestarts < 0) throw new InvalidOperationException("The number of restarts cannot be negative.");
            int q = UseARD ? _d : 1;
            // With the noise estimated or zero, the covariance is σ²(R + gI) and σ² is profiled out in closed form;
            // the optimizer works on (ln ℓ, ln g). A fixed positive noise keeps (ln σ_f², ln ℓ).
            bool profiled = estimateNoise || fixedNoise == 0d;
            int o = profiled ? 0 : 1;
            int m = o + q + (estimateNoise ? 1 : 0);

            // Data-based scales for the starting values and bounds.
            double scale = _p > 0 ? Statistics.PopulationVariance(_y) : Tools.SumProduct(_y, _y) / _n;
            if (!(scale > 0d) || !Tools.IsFinite(scale)) scale = 1d;
            double lnScale = Math.Log(scale);
            var start = new double[m];
            // The second data-based start differs only in the length scales.
            var second = new double[m];
            var lower = new double[m];
            var upper = new double[m];
            // Restarts are drawn from a narrower box of plausible values inside the bounds.
            var restartLower = new double[m];
            var restartUpper = new double[m];
            if (!profiled)
            {
                start[0] = lnScale; lower[0] = lnScale - 6d * Ln10; upper[0] = lnScale + 6d * Ln10;
                restartLower[0] = lnScale - 2d * Ln10; restartUpper[0] = lnScale + 2d * Ln10;
            }
            for (int k = 0; k < q; k++)
            {
                double range = 0d;
                for (int j = UseARD ? k : 0; j < (UseARD ? k + 1 : _d); j++)
                {
                    double lo = double.PositiveInfinity, hi = double.NegativeInfinity;
                    for (int i = 0; i < _n; i++)
                    {
                        lo = Math.Min(lo, _x[i * _d + j]);
                        hi = Math.Max(hi, _x[i * _d + j]);
                    }
                    range = Math.Max(range, hi - lo);
                }
                if (!(range > 0d)) range = 1d;
                double ln = Math.Log(range);
                start[o + k] = ln - Math.Log(2d);
                second[o + k] = ln;
                lower[o + k] = ln - 3d * Ln10;
                upper[o + k] = ln + 3d * Ln10;
                restartLower[o + k] = ln - Ln10;
                restartUpper[o + k] = ln + Ln10;
            }
            if (estimateNoise)
            {
                // The noise-to-signal ratio g.
                start[m - 1] = -2d * Ln10;
                lower[m - 1] = -10d * Ln10;
                upper[m - 1] = 2d * Ln10;
                restartLower[m - 1] = -6d * Ln10;
                restartUpper[m - 1] = 0d;
            }

            // The evaluations overwrite the workspace, so the model is untrained until it is conditioned again.
            var w = _fit ?? new Workspace(_n, _d, _p);
            _fit = null;
            LogMarginalLikelihood = double.NaN;
            var cache = new double[m];
            bool cached = false;
            double profiledVariance = double.NaN;
            Func<double[], Hyperparameters> expand = v =>
            {
                var h = new Hyperparameters { SignalVariance = profiled ? 1d : Math.Exp(v[0]), LengthScales = new double[q] };
                for (int k = 0; k < q; k++) h.LengthScales[k] = Math.Exp(v[o + k]);
                h.NoiseVariance = estimateNoise ? Math.Exp(v[m - 1]) : profiled ? 0d : fixedNoise;
                return h;
            };
            Func<double[], bool> evaluate = v =>
            {
                Array.Copy(v, cache, m);
                cached = Evaluate(w, expand(v), true);
                if (cached && profiled)
                {
                    // ln p = -n/2 ln σ̂² - ½ ln|R + gI| - n/2 (1 + ln 2π) with σ̂² = rᵀ(R + gI)⁻¹r / n.
                    profiledVariance = w.Quadratic / _n;
                    w.LogLikelihood = -0.5 * _n * (Math.Log(profiledVariance) + 1d + Math.Log(2d * Math.PI)) - w.HalfLogDeterminant;
                    cached = profiledVariance > 0d && Tools.IsFinite(w.LogLikelihood);
                }
                return cached;
            };
            Func<double[], double> objective = v => evaluate(v) ? -w.LogLikelihood : double.NaN;
            Func<double[], double[]> gradient = v =>
            {
                bool same = cached;
                for (int i = 0; i < m && same; i++) same = cache[i] == v[i];
                if (!same) evaluate(v);
                var full = new double[q + 2];
                if (!cached)
                {
                    for (int i = 0; i < full.Length; i++) full[i] = double.NaN;
                }
                else
                {
                    // For the profiled likelihood the trace term uses ααᵀ/σ̂² with α = (R + gI)⁻¹r.
                    if (profiled)
                    {
                        double s = 1d / Math.Sqrt(profiledVariance);
                        for (int i = 0; i < _n; i++) w.Alpha[i] *= s;
                    }
                    Gradient(w, expand(v), full);
                    cached = false;
                }
                var g = new double[m];
                if (!profiled) g[0] = -full[0];
                for (int k = 0; k < q; k++) g[o + k] = -full[1 + k];
                if (estimateNoise) g[m - 1] = -full[q + 1];
                return g;
            };

            double[]? best = null;
            double bestValue = double.PositiveInfinity;
            var status = OptimizationStatus.None;
            var prng = new MersenneTwister(PRNGSeed);
            var values = new double[m];
            for (int s = 0; s < NumberOfRestarts + 2; s++)
            {
                if (s == 0) Array.Copy(start, values, m);
                else if (s == 1)
                {
                    Array.Copy(start, values, m);
                    for (int k = 0; k < q; k++) values[o + k] = second[o + k];
                }
                else
                    for (int k = 0; k < m; k++) values[k] = restartLower[k] + (restartUpper[k] - restartLower[k]) * prng.NextDouble();
                var solver = new LBFGSB(objective, m, values, lower, upper, gradient)
                {
                    AbsoluteTolerance = Tolerance,
                    RelativeTolerance = LikelihoodRelativeTolerance,
                    MaxIterations = MaxIterations,
                    RecordTraces = false,
                    ComputeHessian = false,
                    ReportFailure = false
                };
                cached = false;
                solver.Minimize();
                var result = solver.BestParameterSet;
                double runBest = result.Values != null && Tools.IsFinite(result.Fitness) ? result.Fitness : double.PositiveInfinity;
                var point = result.Values;
                var runStatus = solver.Status;
                if (runBest < bestValue)
                {
                    bestValue = runBest;
                    best = point;
                    status = runStatus;
                }
            }
            if (best == null) throw new InvalidOperationException("The log marginal likelihood could not be evaluated at any optimizer start.");

            var final = expand(best);
            double signal = final.SignalVariance;
            if (profiled)
            {
                if (!evaluate(best)) throw new InvalidOperationException("The log marginal likelihood could not be evaluated at the optimum.");
                signal = profiledVariance;
            }
            _sf2 = signal;
            _ell = final.LengthScales;
            _sn2 = profiled ? final.NoiseVariance * signal : final.NoiseVariance;
            Condition(w);
            Status = status;
        }

        #endregion

        #region Likelihood

        /// <summary>
        /// The hyperparameters of one evaluation.
        /// </summary>
        private sealed class Hyperparameters
        {
            /// <summary>The signal variance.</summary>
            public double SignalVariance;

            /// <summary>The length scales (one or one per feature).</summary>
            public double[] LengthScales = Array.Empty<double>();

            /// <summary>The noise variance.</summary>
            public double NoiseVariance;

            /// <summary>Converts log hyperparameters.</summary>
            /// <param name="theta">The log hyperparameters (ln σ_f², ln ℓ, ln σ_n²).</param>
            /// <param name="q">The number of length scales.</param>
            /// <returns>The hyperparameters.</returns>
            public static Hyperparameters FromTheta(double[] theta, int q)
            {
                var h = new Hyperparameters { SignalVariance = Math.Exp(theta[0]), LengthScales = new double[q], NoiseVariance = Math.Exp(theta[q + 1]) };
                for (int k = 0; k < q; k++) h.LengthScales[k] = Math.Exp(theta[1 + k]);
                return h;
            }
        }

        /// <summary>
        /// Buffers for one conditioned model, reused across likelihood evaluations.
        /// </summary>
        private sealed class Workspace
        {
            /// <summary>Creates the buffers.</summary>
            /// <param name="n">The number of observations.</param>
            /// <param name="d">The number of features.</param>
            /// <param name="p">The number of mean basis functions.</param>
            public Workspace(int n, int d, int p)
            {
                L = new double[n * n];
                Diagonal = new double[n];
                Scaled = new double[n * d];
                G = new double[p * n];
                A = new double[p * p];
                Beta = new double[p];
                Alpha = new double[n];
            }

            /// <summary>The Cholesky factor (lower) and, after inversion, the inverse (strict upper).</summary>
            public double[] L;

            /// <summary>The inverse factor transposed, allocated when first needed.</summary>
            public double[]? T;

            /// <summary>The original diagonal, then the diagonal of the inverse.</summary>
            public double[] Diagonal;

            /// <summary>The features divided by the length scales.</summary>
            public double[] Scaled;

            /// <summary>The whitened mean basis, L⁻¹F, one row per basis function.</summary>
            public double[] G;

            /// <summary>The Cholesky factor of GᵀG.</summary>
            public double[] A;

            /// <summary>The generalized least squares coefficients.</summary>
            public double[] Beta;

            /// <summary>The weights α = K⁻¹(y - Fβ).</summary>
            public double[] Alpha;

            /// <summary>The log marginal likelihood.</summary>
            public double LogLikelihood;

            /// <summary>The quadratic form (y - Fβ)ᵀK⁻¹(y - Fβ).</summary>
            public double Quadratic;

            /// <summary>Half the log determinant of K, the sum of the logarithms of the factor's diagonal.</summary>
            public double HalfLogDeterminant;

            /// <summary>The jitter added to the diagonal.</summary>
            public double Jitter;

            /// <summary>Whether L holds the inverse in its strict upper triangle.</summary>
            public bool Inverted;

            /// <summary>Whether T holds the current inverse factor.</summary>
            public bool HasInverseFactor;
        }

        /// <summary>
        /// Conditions the model on the data at the current hyperparameters, with jitter if needed.
        /// </summary>
        /// <param name="w">The workspace to fill.</param>
        /// <exception cref="InvalidOperationException">Thrown when the covariance matrix cannot be factored even with jitter, or the mean basis is rank deficient.</exception>
        private void Condition(Workspace w)
        {
            var hp = new Hyperparameters { SignalVariance = _sf2, LengthScales = _ell, NoiseVariance = _sn2 };
            _fit = null;
            LogMarginalLikelihood = double.NaN;
            if (!Evaluate(w, hp, false))
                throw new InvalidOperationException("The covariance matrix is not positive-definite even with jitter, or the mean basis is rank deficient.");
            _fit = w;
            LogMarginalLikelihood = w.LogLikelihood;
        }

        /// <summary>
        /// Builds and factors the covariance matrix, then computes the mean coefficients, the weights and the log
        /// marginal likelihood.
        /// </summary>
        /// <param name="w">The workspace.</param>
        /// <param name="hp">The hyperparameters.</param>
        /// <param name="strict">When true, a matrix that is not numerically positive-definite fails; otherwise jitter is added.</param>
        /// <returns>True on success.</returns>
        private bool Evaluate(Workspace w, Hyperparameters hp, bool strict)
        {
            int n = _n;
            w.Inverted = false;
            w.HasInverseFactor = false;
            w.Jitter = 0d;
            BuildCovariance(w, hp);
            double tolerance = CholeskyDecomposition.DefaultRelativeTolerance(n);
            if (!CholeskyKernel.TryFactor(w.L, 0, n, n, tolerance, w.Diagonal, out _, out _))
            {
                if (strict) return false;
                BuildCovariance(w, hp);
                var k = new Matrix(n);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j <= i; j++) k[i, j] = k[j, i] = w.L[i * n + j];
                if (!CholeskyDecomposition.TryCholesky(k, out var chol, out double jitter) || chol == null) return false;
                var lower = chol.L;
                for (int i = 0; i < n; i++)
                    for (int j = 0; j <= i; j++) w.L[i * n + j] = lower[i, j];
                w.Jitter = jitter;
            }
            return Solve(w);
        }

        /// <summary>
        /// Computes the generalized least squares mean, the weights and the log marginal likelihood from the factor.
        /// </summary>
        /// <param name="w">The workspace with a valid factor.</param>
        /// <returns>False when the mean basis is rank deficient.</returns>
        private bool Solve(Workspace w)
        {
            int n = _n, p = _p;
            var z = w.Alpha;
            Array.Copy(_y, z, n);
            CholeskyKernel.SolveLowerRows(w.L, 0, n, n, z, 0, 1, null);
            if (p > 0)
            {
                Array.Copy(_basis, w.G, p * n);
                CholeskyKernel.SolveLowerRows(w.L, 0, n, n, w.G, 0, p, null);
                var rhs = new double[p];
                for (int a = 0; a < p; a++)
                {
                    rhs[a] = CrossProductKernel.Dot(w.G, a * n, z, 0, n);
                    for (int b = 0; b <= a; b++)
                        w.A[a * p + b] = CrossProductKernel.Dot(w.G, a * n, w.G, b * n, n);
                }
                if (!CholeskyKernel.TryFactor(w.A, 0, p, p, CholeskyDecomposition.DefaultRelativeTolerance(p), null, out _, out _))
                    return false;
                CholeskyKernel.SolveLowerRows(w.A, 0, p, p, rhs, 0, 1, null);
                CholeskyKernel.SolveUpper(w.A, 0, p, p, rhs, 0);
                Array.Copy(rhs, w.Beta, p);
                for (int a = 0; a < p; a++)
                {
                    double beta = rhs[a];
                    CrossProductKernel.Axpy(-beta, w.G, a * n, z, 0, n);
                }
            }
            double quad = CrossProductKernel.Dot(z, 0, z, 0, n);
            double logDet = 0d;
            for (int i = 0; i < n; i++) logDet += Math.Log(w.L[i * n + i]);
            CholeskyKernel.SolveUpper(w.L, 0, n, n, z, 0);
            w.Quadratic = quad;
            w.HalfLogDeterminant = logDet;
            w.LogLikelihood = -0.5 * quad - logDet - 0.5 * n * Math.Log(2d * Math.PI);
            return Tools.IsFinite(w.LogLikelihood);
        }

        /// <summary>
        /// Writes the covariance matrix of the training data into the lower triangle of the workspace factor.
        /// </summary>
        /// <param name="w">The workspace.</param>
        /// <param name="hp">The hyperparameters.</param>
        private void BuildCovariance(Workspace w, Hyperparameters hp)
        {
            int n = _n;
            if (_distance == null) Scale(_x, n, hp.LengthScales, w.Scaled);
            int chunks = (n + 31) / 32;
            if (n >= CholeskyKernel.ParallelThreshold && _distance == null)
                Parallel.For(0, chunks, c => CovarianceRows(w, hp, c * 32, Math.Min(n, c * 32 + 32)));
            else
                CovarianceRows(w, hp, 0, n);
        }

        /// <summary>
        /// Writes rows [r0, r1) of the lower triangle of the training covariance matrix.
        /// </summary>
        /// <param name="w">The workspace, with the scaled features.</param>
        /// <param name="hp">The hyperparameters.</param>
        /// <param name="r0">The first row.</param>
        /// <param name="r1">One past the last row.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private void CovarianceRows(Workspace w, Hyperparameters hp, int r0, int r1)
        {
            int n = _n, d = _d;
            var l = w.L;
            var xs = w.Scaled;
            double sf2 = hp.SignalVariance;
            for (int a = r0; a < r1; a++)
            {
                int ra = a * n;
                if (_distance == null)
                {
                    int xa = a * d;
                    for (int b = 0; b < a; b++)
                    {
                        int xb = b * d;
                        double r2 = 0d;
                        for (int j = 0; j < d; j++)
                        {
                            double diff = xs[xa + j] - xs[xb + j];
                            r2 += diff * diff;
                        }
                        l[ra + b] = sf2 * Correlation(_form, r2);
                    }
                }
                else
                {
                    double ell = hp.LengthScales[0];
                    for (int b = 0; b < a; b++)
                    {
                        double r = _distance(_rows![a], _rows[b]) / ell;
                        l[ra + b] = sf2 * Correlation(_form, r * r);
                    }
                }
                l[ra + a] = sf2 + hp.NoiseVariance;
            }
        }

        /// <summary>Divides each feature by its length scale.</summary>
        /// <param name="x">The row-major features.</param>
        /// <param name="rows">The number of rows.</param>
        /// <param name="ell">The length scales (one or one per feature).</param>
        /// <param name="scaled">Receives the scaled features.</param>
        private void Scale(double[] x, int rows, double[] ell, double[] scaled)
        {
            int d = _d;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < d; j++)
                    scaled[i * d + j] = x[i * d + j] / ell[ell.Length == 1 ? 0 : j];
        }

        /// <summary>
        /// Returns the correlation at squared scaled distance r².
        /// </summary>
        /// <param name="form">The covariance form.</param>
        /// <param name="r2">The squared scaled distance.</param>
        /// <returns>The correlation ρ(r).</returns>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static double Correlation(CovarianceForm form, double r2)
        {
            double r;
            switch (form)
            {
                case CovarianceForm.SquaredExponential:
                    return Math.Exp(-0.5 * r2);
                case CovarianceForm.Matern32:
                    r = Math.Sqrt(3d * r2);
                    return (1d + r) * Math.Exp(-r);
                case CovarianceForm.Matern52:
                    r = Math.Sqrt(5d * r2);
                    return (1d + r + 5d * r2 / 3d) * Math.Exp(-r);
                case CovarianceForm.Exponential:
                    return Math.Exp(-Math.Sqrt(r2));
                default:
                    r = Math.Sqrt(r2);
                    return r >= 1d ? 0d : 1d - r * (1.5 - 0.5 * r2);
            }
        }

        /// <summary>
        /// Returns h(r) such that <c>∂k/∂ln ℓ_j = k(r) h(r) (Δ_j / ℓ_j)²</c>.
        /// </summary>
        /// <param name="form">The covariance form.</param>
        /// <param name="r2">The squared scaled distance, positive.</param>
        /// <returns>The derivative factor.</returns>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static double DerivativeFactor(CovarianceForm form, double r2)
        {
            double r;
            switch (form)
            {
                case CovarianceForm.SquaredExponential:
                    return 1d;
                case CovarianceForm.Matern32:
                    return 3d / (1d + Math.Sqrt(3d * r2));
                case CovarianceForm.Matern52:
                    r = Math.Sqrt(5d * r2);
                    return 5d / 3d * (1d + r) / (1d + r + 5d * r2 / 3d);
                default:
                    return 1d / Math.Sqrt(r2);
            }
        }

        /// <summary>
        /// Computes the gradient of the log marginal likelihood with respect to the log hyperparameters.
        /// </summary>
        /// <param name="w">The workspace, holding the factor and weights at <paramref name="hp"/>.</param>
        /// <param name="hp">The hyperparameters.</param>
        /// <param name="gradient">Receives (∂/∂ln σ_f², ∂/∂ln ℓ_1..q, ∂/∂ln σ_n²).</param>
        /// <remarks>Overwrites the strict upper triangle of the factor buffer with the inverse covariance.</remarks>
        private void Gradient(Workspace w, Hyperparameters hp, double[] gradient)
        {
            int n = _n, d = _d, q = hp.LengthScales.Length;
            w.T ??= new double[n * n];
            CholeskyKernel.Invert(w.L, n, w.T, w.Diagonal);
            w.Inverted = true;
            w.HasInverseFactor = true;
            var alpha = w.Alpha;
            double sf2 = hp.SignalVariance;
            bool ard = q > 1;
            int width = q + 1;
            const int rowsPerChunk = 32;
            int chunks = (n + rowsPerChunk - 1) / rowsPerChunk;
            var partial = new double[chunks * width];
            if (n >= CholeskyKernel.ParallelThreshold)
                Parallel.For(0, chunks, c => GradientRows(w, sf2, ard, c * rowsPerChunk, Math.Min(n, c * rowsPerChunk + rowsPerChunk), partial, c * width));
            else
                for (int c = 0; c < chunks; c++)
                    GradientRows(w, sf2, ard, c * rowsPerChunk, Math.Min(n, c * rowsPerChunk + rowsPerChunk), partial, c * width);

            Array.Clear(gradient, 0, gradient.Length);
            for (int c = 0; c < chunks; c++)
                for (int k = 0; k < width; k++) gradient[k] += partial[c * width + k];
            double diagonal = 0d;
            for (int a = 0; a < n; a++) diagonal += alpha[a] * alpha[a] - w.Diagonal[a];
            gradient[0] += 0.5 * diagonal * sf2;
            gradient[q + 1] = 0.5 * diagonal * hp.NoiseVariance;
        }

        /// <summary>
        /// Accumulates the off-diagonal trace terms of the gradient for rows [r0, r1).
        /// </summary>
        /// <param name="w">The workspace, with the weights and the inverse in the strict upper triangle.</param>
        /// <param name="sf2">The signal variance.</param>
        /// <param name="ard">Whether each feature has its own length scale.</param>
        /// <param name="r0">The first row.</param>
        /// <param name="r1">One past the last row.</param>
        /// <param name="partial">The partial sums: the signal term, then the length scale terms.</param>
        /// <param name="o">The offset of these rows' partial sums.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private void GradientRows(Workspace w, double sf2, bool ard, int r0, int r1, double[] partial, int o)
        {
            int n = _n, d = _d;
            var l = w.L;
            var alpha = w.Alpha;
            var xs = w.Scaled;
            double gs = 0d;
            for (int a = r0; a < r1; a++)
            {
                int xa = a * d;
                double aa = alpha[a];
                for (int b = 0; b < a; b++)
                {
                    int xb = b * d;
                    double wab = aa * alpha[b] - l[b * n + a];
                    double r2 = 0d;
                    for (int j = 0; j < d; j++)
                    {
                        double diff = xs[xa + j] - xs[xb + j];
                        r2 += diff * diff;
                    }
                    double k = sf2 * Correlation(_form, r2);
                    gs += wab * k;
                    if (r2 <= 0d) continue;
                    double factor = wab * k * DerivativeFactor(_form, r2);
                    if (ard)
                    {
                        for (int j = 0; j < d; j++)
                        {
                            double diff = xs[xa + j] - xs[xb + j];
                            partial[o + 1 + j] += factor * diff * diff;
                        }
                    }
                    else
                    {
                        partial[o + 1] += factor * r2;
                    }
                }
            }
            partial[o] = gs;
        }

        #endregion

        #region Prediction

        /// <summary>
        /// Returns the predictive mean for a single-feature model.
        /// </summary>
        /// <param name="x">The feature values of the prediction points.</param>
        /// <returns>The predictive mean at each point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the model has more than one feature.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        public double[] Predict(double[] x)
        {
            return Predict(ToColumn(x));
        }

        /// <summary>
        /// Returns the predictive mean.
        /// </summary>
        /// <param name="x">The prediction points, with one row per point.</param>
        /// <returns>The predictive mean at each point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when x has the wrong number of columns.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        public double[] Predict(double[,] x)
        {
            return PredictCore(Flatten(x), x.GetLength(0), x.GetLength(1), false, false).Mean;
        }

        /// <summary>
        /// Returns the predictive mean.
        /// </summary>
        /// <param name="x">The prediction points, with one row per point.</param>
        /// <returns>The predictive mean at each point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when x has the wrong number of columns.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        public double[] Predict(Matrix x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            return Predict(x.Array);
        }

        /// <summary>
        /// Returns the predictive mean with a central prediction interval, in a 2D array with columns: lower, mean, upper.
        /// </summary>
        /// <param name="x">The prediction points, with one row per point.</param>
        /// <param name="alpha">One minus the interval coverage. Default = 0.1, a 90% interval.</param>
        /// <param name="includeNoise">Whether the interval covers a new noisy observation rather than the latent function.</param>
        /// <returns>The interval bounds and predictive mean.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when alpha is not in (0, 1).</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        public double[,] Predict(Matrix x, double alpha, bool includeNoise = false)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (!(alpha > 0d && alpha < 1d)) throw new ArgumentOutOfRangeException(nameof(alpha), "Alpha must be between 0 and 1.");
            var (mean, variance) = PredictWithVariance(x.Array, includeNoise);
            double z = Normal.StandardZ(1d - alpha / 2d);
            var result = new double[mean.Length, 3];
            for (int i = 0; i < mean.Length; i++)
            {
                double se = Math.Sqrt(variance[i]);
                result[i, 0] = mean[i] - z * se;
                result[i, 1] = mean[i];
                result[i, 2] = mean[i] + z * se;
            }
            return result;
        }

        /// <summary>
        /// Returns the predictive mean and variance.
        /// </summary>
        /// <param name="x">The prediction points, with one row per point.</param>
        /// <param name="includeNoise">
        /// When false (default), the variance is that of the latent function f; when true, the noise variance and
        /// any jitter are added, giving the variance of a new observation.
        /// </param>
        /// <returns>The predictive mean and variance at each point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when x has the wrong number of columns.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        /// <remarks>
        /// With a constant or linear mean, the variance includes the uncertainty of the estimated mean coefficients,
        /// <c>uᵀ(FᵀK⁻¹F)⁻¹u</c> with <c>u = f(x*) - FᵀK⁻¹k*</c>, as in universal kriging.
        /// </remarks>
        public (double[] Mean, double[] Variance) PredictWithVariance(double[,] x, bool includeNoise = false)
        {
            var (mean, variance) = PredictCore(Flatten(x), x.GetLength(0), x.GetLength(1), true, includeNoise);
            return (mean, variance!);
        }

        /// <summary>
        /// Returns the predictive mean and variance for a single-feature model.
        /// </summary>
        /// <param name="x">The feature values of the prediction points.</param>
        /// <param name="includeNoise">Whether to add the noise variance; see <see cref="PredictWithVariance(double[,], bool)"/>.</param>
        /// <returns>The predictive mean and variance at each point.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the model has more than one feature.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        public (double[] Mean, double[] Variance) PredictWithVariance(double[] x, bool includeNoise = false)
        {
            return PredictWithVariance(ToColumn(x), includeNoise);
        }

        /// <summary>
        /// Computes predictions in blocks of points.
        /// </summary>
        /// <param name="x">The flat row-major prediction points.</param>
        /// <param name="m">The number of points.</param>
        /// <param name="d">The number of columns.</param>
        /// <param name="withVariance">Whether to compute the variance.</param>
        /// <param name="includeNoise">Whether to add the noise variance.</param>
        /// <returns>The mean and, when requested, the variance.</returns>
        internal (double[] Mean, double[]? Variance) PredictCore(double[] x, int m, int d, bool withVariance, bool includeNoise)
        {
            var w = _fit ?? throw new InvalidOperationException("The model must be trained before predicting.");
            if (d != _d) throw new ArgumentException("The prediction points must have the same number of features as the training data.", nameof(x));
            for (int i = 0; i < x.Length; i++)
                if (!Tools.IsFinite(x[i])) throw new ArgumentException("The prediction points must be finite.", nameof(x));
            int n = _n, p = _p;
            var mean = new double[m];
            var variance = withVariance ? new double[m] : null;
            double sf2 = _sf2;
            double noise = _sn2 + w.Jitter;
            double prior = sf2 + (includeNoise || _exactNugget ? noise : 0d);
            const int blockSize = 1024;
            int block = Math.Min(m, blockSize);
            var kstar = new double[block * n];
            var scaled = new double[block * _d];
            var basis = new double[p];
            var u = new double[p];
            for (int start = 0; start < m; start += blockSize)
            {
                int count = Math.Min(blockSize, m - start);
                CrossCovariance(x, start, count, scaled, kstar, noise);
                for (int i = 0; i < count; i++)
                {
                    double mu = CrossProductKernel.Dot(kstar, i * n, w.Alpha, 0, n);
                    for (int a = 0; a < p; a++)
                        mu += w.Beta[a] * (a == 0 ? 1d : x[(start + i) * d + a - 1]);
                    mean[start + i] = mu;
                }
                if (!withVariance) continue;
                CholeskyKernel.SolveLowerRows(w.L, 0, n, n, kstar, 0, count, null);
                for (int i = 0; i < count; i++)
                {
                    int vi = i * n;
                    double v = prior - CrossProductKernel.Dot(kstar, vi, kstar, vi, n);
                    if (p > 0)
                    {
                        for (int a = 0; a < p; a++)
                            u[a] = (a == 0 ? 1d : x[(start + i) * d + a - 1]) - CrossProductKernel.Dot(w.G, a * n, kstar, vi, n);
                        Array.Copy(u, basis, p);
                        CholeskyKernel.SolveLowerRows(w.A, 0, p, p, basis, 0, 1, null);
                        v += CrossProductKernel.Dot(basis, 0, basis, 0, p);
                    }
                    variance![start + i] = Math.Max(v, 0d);
                }
            }
            return (mean, variance);
        }

        /// <summary>
        /// Adds <c>(src[k] - c)²</c> to <c>acc[k]</c> over a contiguous segment, with a separate multiply and add
        /// so that each element is rounded exactly as in the scalar loop.
        /// </summary>
        /// <param name="src">The source array.</param>
        /// <param name="offset">The start of the source segment.</param>
        /// <param name="c">The value subtracted from each source element.</param>
        /// <param name="acc">The accumulator array.</param>
        /// <param name="accOffset">The start of the accumulator segment.</param>
        /// <param name="length">The segment length.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private static void AddSquaredDifferences(double[] src, int offset, double c, double[] acc, int accOffset, int length)
        {
            int k = 0;
#if NET8_0_OR_GREATER
            if (Vector256.IsHardwareAccelerated && length >= 4)
            {
                ref double rs = ref MemoryMarshal.GetArrayDataReference(src);
                ref double ra = ref MemoryMarshal.GetArrayDataReference(acc);
                var vc = Vector256.Create(c);
                for (; k <= length - 4; k += 4)
                {
                    var diff = Vector256.LoadUnsafe(ref rs, (nuint)(offset + k)) - vc;
                    var sum = Vector256.LoadUnsafe(ref ra, (nuint)(accOffset + k)) + diff * diff;
                    sum.StoreUnsafe(ref ra, (nuint)(accOffset + k));
                }
            }
#endif
            for (; k < length; k++)
            {
                double diff = src[offset + k] - c;
                acc[accOffset + k] += diff * diff;
            }
        }

        /// <summary>
        /// Computes the covariance between prediction points and the training points.
        /// </summary>
        /// <param name="x">The flat row-major prediction points.</param>
        /// <param name="start">The first point of the block.</param>
        /// <param name="count">The number of points in the block.</param>
        /// <param name="scaled">A buffer for the scaled points.</param>
        /// <param name="kstar">Receives the count-by-n covariances.</param>
        /// <param name="noise">The nugget added at zero distance when the nugget is exact.</param>
        private void CrossCovariance(double[] x, int start, int count, double[] scaled, double[] kstar, double noise)
        {
            int n = _n, d = _d;
            double[]? transposed = null;
            if (_distance == null)
            {
                for (int i = 0; i < count; i++)
                    for (int j = 0; j < d; j++)
                        scaled[i * d + j] = x[(start + i) * d + j] / _ell[_ell.Length == 1 ? 0 : j];
                // Feature-major training features, so the squared distances to all training points
                // accumulate over contiguous memory, one feature at a time as in BuildCovariance.
                var xs = _fit!.Scaled;
                transposed = new double[d * n];
                for (int b = 0; b < n; b++)
                    for (int j = 0; j < d; j++) transposed[j * n + b] = xs[b * d + j];
            }
            if ((long)count * n >= 1L << 16 && _distance == null)
                Parallel.For(0, count, i => CrossCovarianceRow(x, start, i, scaled, transposed, kstar, noise));
            else
                for (int i = 0; i < count; i++) CrossCovarianceRow(x, start, i, scaled, transposed, kstar, noise);
        }

        /// <summary>
        /// Computes the covariances between one prediction point and the training points.
        /// </summary>
        /// <param name="x">The flat row-major prediction points.</param>
        /// <param name="start">The first point of the block.</param>
        /// <param name="i">The point's index within the block.</param>
        /// <param name="scaled">The scaled points of the block.</param>
        /// <param name="transposed">The feature-major scaled training features, or null with a custom distance.</param>
        /// <param name="kstar">Receives the covariances in row i.</param>
        /// <param name="noise">The nugget added at zero distance when the nugget is exact.</param>
#if NET8_0_OR_GREATER
        [MethodImpl(MethodImplOptions.AggressiveOptimization)]
#endif
        private void CrossCovarianceRow(double[] x, int start, int i, double[] scaled, double[]? transposed, double[] kstar, double noise)
        {
            int n = _n, d = _d;
            double sf2 = _sf2;
            int ki = i * n;
            if (transposed != null)
            {
                Array.Clear(kstar, ki, n);
                for (int j = 0; j < d; j++)
                    AddSquaredDifferences(transposed, j * n, scaled[i * d + j], kstar, ki, n);
            }
            else
            {
                var local = new double[d];
                Array.Copy(x, (start + i) * d, local, 0, d);
                for (int b = 0; b < n; b++)
                {
                    double r = _distance!(local, _rows![b]) / _ell[0];
                    kstar[ki + b] = r * r;
                }
            }
            for (int b = 0; b < n; b++)
            {
                double r2 = kstar[ki + b];
                double k = sf2 * Correlation(_form, r2);
                if (_exactNugget && r2 == 0d) k += noise;
                kstar[ki + b] = k;
            }
        }

        #endregion

        #region Cross-Validation

        /// <summary>
        /// Returns the leave-one-out predictive mean and variance of each training response.
        /// </summary>
        /// <returns>
        /// For each observation i, the mean and variance of <c>y_i</c> predicted from the other observations with the
        /// same hyperparameters, re-estimating the mean coefficients. The residuals are <c>y_i - Mean_i</c>.
        /// </returns>
        /// <exception cref="InvalidOperationException">Thrown when the model has not been trained.</exception>
        /// <remarks>
        /// Uses the closed form of Dubrule (1983): with <c>Q = K⁻¹ - K⁻¹F(FᵀK⁻¹F)⁻¹FᵀK⁻¹</c> and <c>α = Q y</c>,
        /// the residual is <c>α_i / Q_ii</c> and the variance is <c>1 / Q_ii</c>, which includes the noise variance.
        /// For a zero mean this is equation 5.12 of Rasmussen and Williams (2006). Costs one additional O(n³) inversion.
        /// </remarks>
        public (double[] Mean, double[] Variance) LeaveOneOut()
        {
            var w = _fit ?? throw new InvalidOperationException("The model must be trained before cross-validation.");
            int n = _n, p = _p;
            if (!w.HasInverseFactor)
            {
                w.T ??= new double[n * n];
                CholeskyKernel.InvertFactor(w.L, n, w.T);
                w.HasInverseFactor = true;
            }
            var t = w.T!;
            var mean = new double[n];
            var variance = new double[n];
            var h = new double[p];
            for (int i = 0; i < n; i++)
            {
                int ti = i * n + i;
                double qii = CrossProductKernel.Dot(t, ti, t, ti, n - i);
                if (p > 0)
                {
                    for (int a = 0; a < p; a++)
                        h[a] = CrossProductKernel.Dot(t, ti, w.G, a * n + i, n - i);
                    CholeskyKernel.SolveLowerRows(w.A, 0, p, p, h, 0, 1, null);
                    qii -= CrossProductKernel.Dot(h, 0, h, 0, p);
                }
                variance[i] = 1d / qii;
                mean[i] = _y[i] - w.Alpha[i] / qii;
            }
            return (mean, variance);
        }

        #endregion

    }
}
