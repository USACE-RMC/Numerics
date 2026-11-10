using Numerics.Mathematics.LinearAlgebra;
using Numerics.Sampling;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Elastic net regression for a Gaussian response, fitted along a path of penalties by coordinate descent,
    /// with k-fold cross-validation for penalty selection.
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
    /// For each penalty λ on the path the method solves
    /// </para>
    /// <para>
    /// min over (β₀, β) of  1/(2n) Σᵢ (yᵢ − β₀ − xᵢᵀβ)² + λ Σⱼ vⱼ [ (1 − α)/2 βⱼ² + α |βⱼ| ],
    /// </para>
    /// <para>
    /// where α ∈ [0, 1] mixes the ridge (α = 0) and lasso (α = 1) penalties and vⱼ are penalty factors. The
    /// algorithm follows glmnet (Friedman, Hastie and Tibshirani, 2010) in its covariance-update form: the predictors
    /// are standardized to unit population variance and the response to unit population variance, each coordinate
    /// update is a soft-thresholding step on the gradient, the gradient is kept current from cached inner products
    /// of the predictors, full sweeps alternate with sweeps over the active set, and each solution warm-starts the
    /// next. The default penalty sequence is 100 values log-spaced from λ_max, the smallest penalty with all
    /// coefficients zero, down to λ_max·ε with ε = 10⁻⁴ (or 10⁻² when n &lt; p). A sweep has converged when the
    /// largest weighted squared coefficient change, xᵥⱼ·Δβⱼ², is below <see cref="Tolerance"/>. The path stops
    /// early when the fraction of deviance explained exceeds 0.999 or changes by less than 10⁻⁵ relative to its
    /// value, after at least five penalties.
    /// </para>
    /// <para>
    /// Because the response is standardized before fitting, as in glmnet, the reported λ and coefficients reproduce
    /// glmnet for every α. The ridge part of the penalty is then effectively scaled by 1/sd(y), so for α &lt; 1 the
    /// solution equals scikit-learn's <c>ElasticNet</c> fitted to standardized X and standardized y with
    /// <c>alpha</c> = λ/sd(y) and <c>l1_ratio</c> = α (population standard deviations), back-transformed. For α = 1
    /// the response scaling cancels and <c>alpha</c> = λ on standardized X.
    /// </para>
    /// <para>
    /// Inner products are formed only for predictors that enter the model. When one is needed, it is computed together
    /// with the inactive predictors most likely to enter next (largest gradient relative to the penalty factor), in a
    /// batch that doubles with the cache size, by a cache-blocked SIMD kernel. The inner products are shared between
    /// the cross-validation folds through per-fold partial sums, so k-fold cross-validation costs about one extra
    /// pass over the data instead of k refits. Each fold is solved exactly at the full-data penalties, as in
    /// scikit-learn's <c>ElasticNetCV</c>; cv.glmnet instead fits each fold on its own penalty sequence and interpolates
    /// the predictions, so its errors can differ slightly. Cross-validation error is the fold-size weighted mean squared
    /// prediction error, its standard error follows cv.glmnet, λ_min is the largest penalty with the smallest error,
    /// and λ_1se is the largest penalty whose error is within one standard error of that minimum.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths for generalized linear models via
    /// coordinate descent. Journal of Statistical Software, 33(1), 1-22. https://doi.org/10.18637/jss.v033.i01
    /// </para>
    /// <para>
    /// Zou, H., and Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the Royal
    /// Statistical Society: Series B, 67(2), 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
    /// </para>
    /// </remarks>
    public class ElasticNetRegression
    {
        #region Construction

        /// <summary>
        /// Constructs an elastic net regression.
        /// </summary>
        /// <param name="x">The n × p matrix of predictor values.</param>
        /// <param name="y">The response vector of length n.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the dimensions disagree, there are fewer than two
        /// observations or no predictors, or a value is not finite.</exception>
        public ElasticNetRegression(double[,] x, double[] y)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (y == null) throw new ArgumentNullException(nameof(y));
            int n = x.GetLength(0), p = x.GetLength(1);
            if (y.Length != n) throw new ArgumentException("X and Y must have the same number of rows.", nameof(y));
            if (n < 2) throw new ArgumentException("There must be at least two observations.", nameof(y));
            if (p < 1) throw new ArgumentException("There must be at least one predictor.", nameof(x));
            for (int i = 0; i < n; i++)
            {
                if (!Tools.IsFinite(y[i])) throw new ArgumentException("The response must be finite.", nameof(y));
                for (int j = 0; j < p; j++)
                    if (!Tools.IsFinite(x[i, j])) throw new ArgumentException("The predictors must be finite.", nameof(x));
            }
            X = x;
            Y = y;
        }

        /// <summary>
        /// Constructs an elastic net regression.
        /// </summary>
        /// <param name="x">The n × p matrix of predictor values.</param>
        /// <param name="y">The response vector of length n.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the dimensions disagree, there are fewer than two
        /// observations or no predictors, or a value is not finite.</exception>
        public ElasticNetRegression(Matrix x, Vector y)
            : this(x?.Array!, y?.Array!)
        {
        }

        #endregion

        #region Members

        private const double Big = 9.9e35;
        private const double MinimumRatio = 1e-6;
        private const double MinimumDevianceChange = 1e-5;
        private const double MaximumDevianceRatio = 0.999;
        private const int MinimumPathLength = 5;

        /// <summary>
        /// Below this number of data values the data preparation runs on one thread.
        /// </summary>
        private const long ParallelThreshold = 1L << 18;

        /// <summary>
        /// Up to this many multiply-adds, all inner products are computed in a single pass.
        /// </summary>
        private const long FullGramBudget = 1L << 25;

        /// <summary>
        /// The n × p matrix of predictor values.
        /// </summary>
        public double[,] X { get; private set; }

        /// <summary>
        /// The response vector.
        /// </summary>
        public double[] Y { get; private set; }

        /// <summary>
        /// The number of observations, n.
        /// </summary>
        public int SampleSize => Y.Length;

        /// <summary>
        /// The number of predictors, p.
        /// </summary>
        public int NumberOfPredictors => X.GetLength(1);

        /// <summary>
        /// The elastic net mixing parameter α ∈ [0, 1]: 1 is the lasso and 0 is ridge regression. Default = 1.
        /// </summary>
        public double Alpha { get; set; } = 1.0;

        /// <summary>
        /// The number of penalties on the default path. Default = 100. Ignored when <see cref="UserLambdas"/> is set.
        /// </summary>
        public int NumberOfLambdas { get; set; } = 100;

        /// <summary>
        /// The ratio λ_min/λ_max of the default path. NaN (the default) uses 10⁻⁴ when n ≥ p and 10⁻² when n &lt; p.
        /// </summary>
        public double LambdaMinRatio { get; set; } = double.NaN;

        /// <summary>
        /// An optional user-supplied penalty sequence. It is sorted into decreasing order and the path does not stop early.
        /// </summary>
        public double[]? UserLambdas { get; set; }

        /// <summary>
        /// Determines whether the predictors are standardized to unit variance before fitting. Coefficients are always
        /// reported on the original scale. Default = true.
        /// </summary>
        public bool Standardize { get; set; } = true;

        /// <summary>
        /// Determines whether an unpenalized intercept is fitted. Default = true.
        /// </summary>
        public bool HasIntercept { get; set; } = true;

        /// <summary>
        /// The convergence threshold for coordinate descent: the largest weighted squared coefficient change in a sweep,
        /// relative to a unit-variance response. Default = 1E-7.
        /// </summary>
        public double Tolerance { get; set; } = 1E-7;

        /// <summary>
        /// The maximum total number of coordinate-descent sweeps over all penalties. Default = 100,000.
        /// </summary>
        public int MaxIterations { get; set; } = 100000;

        /// <summary>
        /// Optional nonnegative penalty factors vⱼ, one per predictor. A factor of zero leaves the predictor unpenalized.
        /// The factors are rescaled to sum to p.
        /// </summary>
        public double[]? PenaltyFactors { get; set; }

        /// <summary>
        /// Determines whether the path has been fitted.
        /// </summary>
        public bool IsFitted { get; private set; }

        /// <summary>
        /// Determines whether coordinate descent converged at every penalty. When false, the path is truncated at the
        /// last converged penalty.
        /// </summary>
        public bool Converged { get; private set; }

        /// <summary>
        /// The total number of coordinate-descent sweeps used to fit the path.
        /// </summary>
        public int Iterations { get; private set; }

        /// <summary>
        /// The penalties λ of the fitted path, in decreasing order.
        /// </summary>
        public double[] Lambdas { get; private set; } = null!;

        /// <summary>
        /// The intercept at each penalty.
        /// </summary>
        public double[] Intercepts { get; private set; } = null!;

        /// <summary>
        /// The coefficients on the original scale, indexed [penalty, predictor].
        /// </summary>
        public double[,] Coefficients { get; private set; } = null!;

        /// <summary>
        /// The fraction of the null deviance explained at each penalty.
        /// </summary>
        public double[] DevianceRatio { get; private set; } = null!;

        /// <summary>
        /// The number of nonzero coefficients at each penalty.
        /// </summary>
        public int[] NonzeroCount { get; private set; } = null!;

        /// <summary>
        /// Determines whether cross-validation has been run.
        /// </summary>
        public bool IsCrossValidated { get; private set; }

        /// <summary>
        /// The fold label of each observation used in the last cross-validation.
        /// </summary>
        public int[] FoldIds { get; private set; } = null!;

        /// <summary>
        /// The cross-validated mean squared prediction error at each penalty.
        /// </summary>
        public double[] CrossValidationError { get; private set; } = null!;

        /// <summary>
        /// The standard error of the cross-validated error at each penalty.
        /// </summary>
        public double[] CrossValidationStandardError { get; private set; } = null!;

        /// <summary>
        /// The index of λ_min, the largest penalty with the smallest cross-validated error.
        /// </summary>
        public int LambdaMinIndex { get; private set; } = -1;

        /// <summary>
        /// The index of λ_1se, the largest penalty with cross-validated error within one standard error of the minimum.
        /// </summary>
        public int LambdaOneStandardErrorIndex { get; private set; } = -1;

        /// <summary>
        /// The penalty λ_min.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when cross-validation has not been run on the current path.</exception>
        public double LambdaMin => IsCrossValidated ? Lambdas[LambdaMinIndex] : throw new InvalidOperationException("Cross-validation must be run first.");

        /// <summary>
        /// The penalty λ_1se.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when cross-validation has not been run on the current path.</exception>
        public double LambdaOneStandardError => IsCrossValidated ? Lambdas[LambdaOneStandardErrorIndex] : throw new InvalidOperationException("Cross-validation must be run first.");

        #endregion

        #region Methods

        /// <summary>
        /// Fits the coefficient path.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown when a setting is invalid or the response is constant.</exception>
        public void Fit()
        {
            Validate();
            int n = SampleSize;
            var order = new int[n];
            for (int i = 0; i < n; i++) order[i] = i;
            var cache = new ColumnCache(X, Y, order, new[] { 0, n }, HasIntercept);
            var solver = new PathSolver(this, cache, -1);
            solver.Run(UserLambdas == null ? null : SortDescending(UserLambdas));

            int m = solver.PathLength;
            int p = NumberOfPredictors;
            Lambdas = new double[m];
            Intercepts = new double[m];
            Coefficients = new double[m, p];
            DevianceRatio = new double[m];
            NonzeroCount = new int[m];
            for (int k = 0; k < m; k++)
            {
                Lambdas[k] = solver.Lambda[k];
                Intercepts[k] = solver.Intercept[k];
                DevianceRatio[k] = solver.Rsq[k];
                int nz = 0;
                for (int j = 0; j < p; j++)
                {
                    double b = solver.Beta[k * p + j];
                    Coefficients[k, j] = b;
                    if (b != 0d) nz++;
                }
                NonzeroCount[k] = nz;
            }
            if (UserLambdas == null && m >= 3)
                Lambdas[0] = Math.Exp(2d * Math.Log(Lambdas[1]) - Math.Log(Lambdas[2]));
            Iterations = solver.Passes;
            Converged = solver.Converged;
            IsFitted = true;
            // Cross-validation results belong to the previous path.
            IsCrossValidated = false;
            LambdaMinIndex = -1;
            LambdaOneStandardErrorIndex = -1;
            FoldIds = null!;
            CrossValidationError = null!;
            CrossValidationStandardError = null!;
        }

        /// <summary>
        /// Runs k-fold cross-validation over the fitted penalty path, fitting the path first if needed.
        /// </summary>
        /// <param name="folds">The number of folds, at least 3. Default = 10.</param>
        /// <param name="seed">The seed of the random fold assignment. Default = 12345.</param>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the number of folds is less than 3 or more than n.</exception>
        public void CrossValidate(int folds = 10, int seed = 12345)
        {
            int n = SampleSize;
            if (folds < 3 || folds > n) throw new ArgumentOutOfRangeException(nameof(folds), "The number of folds must be between 3 and n.");
            var ids = new int[n];
            for (int i = 0; i < n; i++) ids[i] = i % folds + 1;
            var rng = new MersenneTwister(seed);
            for (int i = n - 1; i > 0; i--)
            {
                int k = rng.Next(i + 1);
                (ids[i], ids[k]) = (ids[k], ids[i]);
            }
            CrossValidate(ids);
        }

        /// <summary>
        /// Runs cross-validation with user-supplied fold labels, fitting the path first if needed.
        /// </summary>
        /// <param name="foldIds">The fold label of each observation; at least three distinct labels.</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="foldIds"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the length is not n or there are fewer than three folds.</exception>
        public void CrossValidate(int[] foldIds)
        {
            if (foldIds == null) throw new ArgumentNullException(nameof(foldIds));
            int n = SampleSize, p = NumberOfPredictors;
            if (foldIds.Length != n) throw new ArgumentException("There must be one fold label per observation.", nameof(foldIds));
            var labels = foldIds.Distinct().OrderBy(v => v).ToArray();
            int K = labels.Length;
            if (K < 3) throw new ArgumentException("There must be at least three folds.", nameof(foldIds));
            if (!IsFitted) Fit();

            // Sort the rows by fold so that each fold is a contiguous segment.
            var index = new Dictionary<int, int>();
            for (int s = 0; s < K; s++) index[labels[s]] = s;
            var counts = new int[K];
            foreach (var f in foldIds) counts[index[f]]++;
            var segments = new int[K + 1];
            for (int s = 0; s < K; s++) segments[s + 1] = segments[s] + counts[s];
            var next = (int[])segments.Clone();
            var order = new int[n];
            for (int i = 0; i < n; i++) order[next[index[foldIds[i]]]++] = i;

            var cache = new ColumnCache(X, Y, order, segments, HasIntercept);
            var lambdas = (double[])Lambdas.Clone();
            int m = lambdas.Length;
            var foldError = new double[K][];
            // Folds run in order so that the batches of cached inner products, and hence every result, are reproducible.
            for (int f = 0; f < K; f++)
            {
                var solver = new PathSolver(this, cache, f);
                solver.Run(lambdas);
                foldError[f] = cache.HoldoutError(f, solver, m);
            }

            // Fold-size weighted mean and standard error, as in cv.glmnet.
            var cvm = new double[m];
            var cvsd = new double[m];
            for (int k = 0; k < m; k++)
            {
                double sw = 0, s = 0;
                for (int f = 0; f < K; f++) { sw += counts[f]; s += counts[f] * foldError[f][k]; }
                cvm[k] = s / sw;
                double v = 0;
                for (int f = 0; f < K; f++) { double d = foldError[f][k] - cvm[k]; v += counts[f] * d * d; }
                cvsd[k] = Math.Sqrt(v / sw / (K - 1));
            }
            int imin = 0;
            for (int k = 1; k < m; k++) if (cvm[k] < cvm[imin]) imin = k;
            double limit = cvm[imin] + cvsd[imin];
            int i1se = 0;
            while (cvm[i1se] > limit) i1se++;

            FoldIds = (int[])foldIds.Clone();
            CrossValidationError = cvm;
            CrossValidationStandardError = cvsd;
            LambdaMinIndex = imin;
            LambdaOneStandardErrorIndex = i1se;
            IsCrossValidated = true;
        }

        /// <summary>
        /// Returns the coefficients (without the intercept) at a penalty on the path.
        /// </summary>
        /// <param name="index">The penalty index.</param>
        /// <returns>The p coefficients on the original scale.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the path has not been fitted.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the index is outside the path.</exception>
        public double[] GetCoefficients(int index)
        {
            CheckIndex(index);
            var b = new double[NumberOfPredictors];
            for (int j = 0; j < b.Length; j++) b[j] = Coefficients[index, j];
            return b;
        }

        /// <summary>
        /// Returns the indexes of the predictors with nonzero coefficients at a penalty on the path.
        /// </summary>
        /// <param name="index">The penalty index.</param>
        /// <returns>The selected predictor indexes, in increasing order.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the path has not been fitted.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the index is outside the path.</exception>
        public int[] SelectedVariables(int index)
        {
            CheckIndex(index);
            var selected = new List<int>();
            for (int j = 0; j < NumberOfPredictors; j++)
                if (Coefficients[index, j] != 0d) selected.Add(j);
            return selected.ToArray();
        }

        /// <summary>
        /// Predicts the response at a penalty on the path.
        /// </summary>
        /// <param name="x">An m × p matrix of predictor values.</param>
        /// <param name="index">The penalty index.</param>
        /// <returns>The m predictions.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the path has not been fitted.</exception>
        /// <exception cref="ArgumentException">Thrown when the number of columns is not p.</exception>
        public double[] Predict(double[,] x, int index)
        {
            CheckIndex(index);
            if (x == null) throw new ArgumentNullException(nameof(x));
            int p = NumberOfPredictors;
            if (x.GetLength(1) != p) throw new ArgumentException("The number of columns must equal the number of predictors.", nameof(x));
            var result = new double[x.GetLength(0)];
            for (int i = 0; i < result.Length; i++)
            {
                double s = Intercepts[index];
                for (int j = 0; j < p; j++) s += Coefficients[index, j] * x[i, j];
                result[i] = s;
            }
            return result;
        }

        /// <summary>
        /// Predicts the response for one observation at a penalty on the path.
        /// </summary>
        /// <param name="x">The p predictor values.</param>
        /// <param name="index">The penalty index.</param>
        /// <returns>The prediction.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the path has not been fitted.</exception>
        /// <exception cref="ArgumentException">Thrown when the length is not p.</exception>
        public double Predict(double[] x, int index)
        {
            CheckIndex(index);
            if (x == null) throw new ArgumentNullException(nameof(x));
            if (x.Length != NumberOfPredictors) throw new ArgumentException("The length must equal the number of predictors.", nameof(x));
            double s = Intercepts[index];
            for (int j = 0; j < x.Length; j++) s += Coefficients[index, j] * x[j];
            return s;
        }

        /// <summary>
        /// Checks that the path is fitted and the index is on it.
        /// </summary>
        private void CheckIndex(int index)
        {
            if (!IsFitted) throw new InvalidOperationException("The path must be fitted first.");
            if (index < 0 || index >= Lambdas.Length) throw new ArgumentOutOfRangeException(nameof(index), "The index must be on the fitted path.");
        }

        /// <summary>
        /// Validates the settings.
        /// </summary>
        private void Validate()
        {
            int p = NumberOfPredictors;
            if (double.IsNaN(Alpha) || Alpha < 0d || Alpha > 1d) throw new ArgumentException("Alpha must be between 0 and 1.", nameof(Alpha));
            if (!(Tolerance > 0d)) throw new ArgumentException("The tolerance must be positive.", nameof(Tolerance));
            if (MaxIterations < 1) throw new ArgumentException("The maximum iterations must be positive.", nameof(MaxIterations));
            if (UserLambdas == null)
            {
                if (NumberOfLambdas < 2) throw new ArgumentException("The path must have at least two penalties.", nameof(NumberOfLambdas));
                if (!double.IsNaN(LambdaMinRatio) && !(LambdaMinRatio > 0d && LambdaMinRatio < 1d))
                    throw new ArgumentException("The penalty ratio must be between 0 and 1.", nameof(LambdaMinRatio));
            }
            else
            {
                if (UserLambdas.Length < 1) throw new ArgumentException("The penalty sequence must not be empty.", nameof(UserLambdas));
                foreach (var l in UserLambdas)
                    if (!(l >= 0d) || double.IsInfinity(l)) throw new ArgumentException("The penalties must be finite and nonnegative.", nameof(UserLambdas));
            }
            if (PenaltyFactors != null)
            {
                if (PenaltyFactors.Length != p) throw new ArgumentException("There must be one penalty factor per predictor.", nameof(PenaltyFactors));
                if (PenaltyFactors.Any(v => !Tools.IsFinite(v))) throw new ArgumentException("The penalty factors must be finite.", nameof(PenaltyFactors));
                if (PenaltyFactors.Max() <= 0d) throw new ArgumentException("At least one penalty factor must be positive.", nameof(PenaltyFactors));
            }
        }

        /// <summary>
        /// Returns a copy sorted in decreasing order.
        /// </summary>
        private static double[] SortDescending(double[] values)
        {
            var copy = (double[])values.Clone();
            Array.Sort(copy);
            Array.Reverse(copy);
            return copy;
        }

        #endregion

        #region Support

        /// <summary>
        /// Column-major, shifted copy of the data with lazily computed per-segment cross products.
        /// </summary>
        /// <remarks>
        /// Columns 0..p − 1 hold xⱼ − μⱼ and column p holds y − μ_y, where μ are the means over all rows when an
        /// intercept is fitted and zero otherwise. Rows are stored in the order given, so that each cross-validation
        /// fold is a contiguous segment. For each requested column k the cache holds, per segment s, the sums
        /// Σ x̃ⱼ x̃ₖ over the rows of s for every j ≤ p; a training set's cross products are sums over its segments.
        /// </remarks>
        private sealed class ColumnCache
        {
            public ColumnCache(double[,] x, double[] y, int[] order, int[] segments, bool shift)
            {
                int n = y.Length, p = x.GetLength(1);
                N = n;
                P = p;
                Segments = segments;
                S = segments.Length - 1;
                Data = new double[n * (p + 1)];
                Shift = new double[p + 1];
                Sum = new double[S, p + 1];
                SumSquares = new double[S, p + 1];
                _minimum = new double[S, p];
                _maximum = new double[S, p];
                Columns = new double[p + 1][][];
                var data = Data;

                const int block = 8;
                Action<int> prepare = b =>
                {
                    int c0 = b * block, c1 = Math.Min(p + 1, c0 + block);
                    for (int c = c0; c < c1; c++)
                    {
                        // Means accumulate in the original row order so the shift does not depend on the folds.
                        double mean = 0d;
                        if (shift)
                        {
                            double s = 0d;
                            if (c < p) for (int i = 0; i < n; i++) s += x[i, c];
                            else for (int i = 0; i < n; i++) s += y[i];
                            mean = s / n;
                        }
                        Shift[c] = mean;
                        int o = c * n;
                        if (c < p) for (int r = 0; r < n; r++) data[o + r] = x[order[r], c] - mean;
                        else for (int r = 0; r < n; r++) data[o + r] = y[order[r]] - mean;
                        for (int seg = 0; seg < S; seg++)
                        {
                            int r0 = segments[seg], r1 = segments[seg + 1];
                            double s1 = 0d, s2 = 0d;
                            double min = double.PositiveInfinity, max = double.NegativeInfinity;
                            for (int r = r0; r < r1; r++)
                            {
                                double v = data[o + r];
                                s1 += v;
                                s2 += v * v;
                                if (v < min) min = v;
                                if (v > max) max = v;
                            }
                            Sum[seg, c] = s1;
                            SumSquares[seg, c] = s2;
                            if (c < p) { _minimum[seg, c] = min; _maximum[seg, c] = max; }
                        }
                    }
                };
                int blocks = (p + block) / block;
                if ((long)n * (p + 1) < ParallelThreshold) for (int b = 0; b < blocks; b++) prepare(b);
                else Parallel.For(0, blocks, prepare);

                // Small problems compute every inner product in one pass; larger ones compute them on demand.
                if ((long)n * (p + 1) * (p + 1) <= FullGramBudget)
                {
                    var all = new int[p + 1];
                    for (int j = 0; j <= p; j++) all[j] = j;
                    Ensure(all);
                }
                else Ensure(new[] { p });
            }

            private readonly object _lock = new object();
            private readonly double[,] _minimum, _maximum;

            public readonly int N, P, S;
            public readonly int[] Segments;
            public readonly double[] Data;
            public readonly double[] Shift;
            public readonly double[,] Sum;
            public readonly double[,] SumSquares;
            public readonly double[][][] Columns;
            public int CachedCount;

            /// <summary>
            /// Determines whether column c is constant over all segments except <paramref name="excluded"/>.
            /// </summary>
            public bool IsConstant(int c, int excluded)
            {
                double min = double.PositiveInfinity, max = double.NegativeInfinity;
                for (int s = 0; s < S; s++)
                {
                    if (s == excluded || Segments[s + 1] == Segments[s]) continue;
                    min = Math.Min(min, _minimum[s, c]);
                    max = Math.Max(max, _maximum[s, c]);
                }
                return !(max > min);
            }

            /// <summary>
            /// Computes and caches the requested columns (in one pass over the data) if not already cached.
            /// </summary>
            public void Ensure(IList<int> request)
            {
                lock (_lock)
                {
                    var right = request.Where(k => Columns[k] == null).Distinct().ToArray();
                    if (right.Length == 0) return;
                    var left = new List<int>();
                    for (int j = 0; j <= P; j++) if (Columns[j] == null) left.Add(j);
                    var leftArray = left.ToArray();
                    int L = leftArray.Length;
                    var results = new double[S][];
                    for (int s = 0; s < S; s++) results[s] = new double[L * right.Length];
                    bool symmetric = right.Length == L && right.SequenceEqual(leftArray);
                    CrossProductKernel.Compute(Data, N, leftArray, right, Segments, results, symmetric);
                    for (int r = 0; r < right.Length; r++)
                    {
                        int k = right[r];
                        var col = new double[S][];
                        for (int s = 0; s < S; s++)
                        {
                            var c = new double[P + 1];
                            for (int l = 0; l < L; l++) c[leftArray[l]] = results[s][r * L + l];
                            for (int j = 0; j <= P; j++)
                                if (Columns[j] != null) c[j] = Columns[j][s][k];
                            col[s] = c;
                        }
                        Columns[k] = col;
                        CachedCount++;
                    }
                }
            }

            /// <summary>
            /// Mean squared prediction error on the rows of segment f at each penalty of a fold's path.
            /// </summary>
            public double[] HoldoutError(int f, PathSolver solver, int m)
            {
                int r0 = Segments[f], len = Segments[f + 1] - r0;
                var error = new double[m];
                int p = P;
                Parallel.For(0, Math.Min(m, solver.PathLength), k =>
                {
                    var r = new double[len];
                    // Residual in shifted coordinates: y - b0 - Σβx = ỹ + (μ_y - b0 - Σβμ) - Σβx̃.
                    double c0 = Shift[p] - solver.Intercept[k];
                    for (int j = 0; j < p; j++) c0 -= solver.Beta[k * p + j] * Shift[j];
                    Array.Copy(Data, p * N + r0, r, 0, len);
                    for (int i = 0; i < len; i++) r[i] += c0;
                    for (int j = 0; j < p; j++)
                    {
                        double b = solver.Beta[k * p + j];
                        if (b != 0d) CrossProductKernel.Axpy(-b, Data, j * N + r0, r, 0, len);
                    }
                    error[k] = CrossProductKernel.Dot(r, 0, r, 0, len) / len;
                });
                for (int k = solver.PathLength; k < m; k++) error[k] = double.NaN;
                return error;
            }
        }

        /// <summary>
        /// Covariance-update coordinate descent along a penalty path for one training set.
        /// </summary>
        private sealed class PathSolver
        {
            public PathSolver(ElasticNetRegression model, ColumnCache cache, int excluded)
            {
                _model = model;
                _cache = cache;
                _excluded = excluded;
                int p = cache.P;
                _p = p;
                _ju = new bool[p];
                _xs = new double[p];
                _xv = new double[p];
                _mx = new double[p];
                _g = new double[p];
                _a = new double[p];
                _vp = new double[p];
                _mm = new int[p];
                _ia = new int[p];
                _da = new double[p];
                _work = new double[p];
                _c = new double[p][];
                for (int j = 0; j < p; j++) _mm[j] = -1;

                int S = cache.S;
                double nt = 0d;
                for (int s = 0; s < S; s++) if (s != excluded) nt += cache.Segments[s + 1] - cache.Segments[s];
                _nt = nt;
                bool intr = model.HasIntercept, isd = model.Standardize;

                var yCol = TrainingColumn(p);
                double sy = TrainingSum(p);
                _my = intr ? sy / nt : 0d;
                double vy = TrainingSumSquares(p) / nt - _my * _my;
                if (!(vy > 0d)) throw new ArgumentException("The response is constant; the Gaussian elastic net cannot be standardized.");
                _ys = Math.Sqrt(vy);

                for (int j = 0; j < p; j++)
                {
                    _ju[j] = !cache.IsConstant(j, excluded);
                    if (!_ju[j]) { _xs[j] = 1d; continue; }
                    double sj = TrainingSum(j);
                    double mean = sj / nt;
                    _mx[j] = intr ? mean : 0d;
                    double ss = TrainingSumSquares(j) / nt;
                    double vc = ss - mean * mean;
                    _xs[j] = isd ? Math.Sqrt(vc) : 1d;
                    _xv[j] = (ss - _mx[j] * _mx[j]) / (_xs[j] * _xs[j]);
                    _g[j] = (yCol[j] / nt - _mx[j] * _my) / (_xs[j] * _ys);
                }

                var vp = model.PenaltyFactors;
                double sum = 0d;
                for (int j = 0; j < p; j++) { _vp[j] = vp == null ? 1d : Math.Max(0d, vp[j]); sum += _vp[j]; }
                for (int j = 0; j < p; j++) _vp[j] *= p / sum;
            }

            private readonly ElasticNetRegression _model;
            private readonly ColumnCache _cache;
            private readonly int _excluded, _p;
            private readonly double _nt, _my, _ys;
            private readonly bool[] _ju;
            private readonly double[] _xs, _xv, _mx, _g, _a, _vp, _da, _work;
            private readonly int[] _mm, _ia;
            private readonly double[][] _c;
            private int _nin;
            private double _rsq, _dlx;
            private bool _warm;

            public int PathLength;
            public int Passes;
            public bool Converged = true;
            public double[] Lambda = null!, Intercept = null!, Rsq = null!, Beta = null!;

            /// <summary>
            /// Training-set sum over segments of the cached cross products of column k.
            /// </summary>
            private double[] TrainingColumn(int k)
            {
                var col = _cache.Columns[k];
                var t = new double[_cache.P + 1];
                for (int s = 0; s < _cache.S; s++)
                {
                    if (s == _excluded) continue;
                    var c = col[s];
                    for (int j = 0; j < t.Length; j++) t[j] += c[j];
                }
                return t;
            }

            private double TrainingSum(int c)
            {
                double t = 0d;
                for (int s = 0; s < _cache.S; s++) if (s != _excluded) t += _cache.Sum[s, c];
                return t;
            }

            private double TrainingSumSquares(int c)
            {
                double t = 0d;
                for (int s = 0; s < _cache.S; s++) if (s != _excluded) t += _cache.SumSquares[s, c];
                return t;
            }

            /// <summary>
            /// Adds predictor k to the active set and stores its standardized inner products with every predictor.
            /// </summary>
            private void Activate(int k)
            {
                if (_cache.Columns[k] == null) Prefetch(k);
                var t = TrainingColumn(k);
                var c = new double[_p];
                for (int j = 0; j < _p; j++)
                {
                    if (!_ju[j]) continue;
                    if (j == k) c[j] = _xv[j];
                    else if (_mm[j] >= 0) c[j] = _c[_mm[j]][k];
                    else c[j] = (t[j] / _nt - _mx[j] * _mx[k]) / (_xs[j] * _xs[k]);
                }
                _c[_nin] = c;
                _mm[k] = _nin;
                _ia[_nin] = k;
                _nin++;
            }

            /// <summary>
            /// Computes column k together with the inactive predictors most likely to enter next, ranked by
            /// |gradient|/penalty factor. The batch doubles with the cache size, so the number of passes over the data
            /// grows logarithmically with the number of cached columns.
            /// </summary>
            private void Prefetch(int k)
            {
                int batch = Math.Max(16, _cache.CachedCount);
                var candidates = new List<int>();
                for (int j = 0; j < _p; j++)
                    if (j != k && _ju[j] && _mm[j] < 0 && _cache.Columns[j] == null) candidates.Add(j);
                var score = new double[_p];
                foreach (var j in candidates) score[j] = _vp[j] > 0d ? Math.Abs(_g[j]) / _vp[j] : double.PositiveInfinity;
                var request = new List<int> { k };
                request.AddRange(candidates.OrderByDescending(j => score[j]).ThenBy(j => j).Take(batch - 1));
                _cache.Ensure(request);
            }

            /// <summary>
            /// One coordinate update. Returns the coefficient change.
            /// </summary>
            private double Update(int k, double ab, double dem)
            {
                double ak = _a[k];
                double u = _g[k] + ak * _xv[k];
                double v = Math.Abs(u) - _vp[k] * ab;
                double anew = 0d;
                if (v > 0d) anew = Tools.Clamp((u >= 0d ? v : -v) / (_xv[k] + _vp[k] * dem), -Big, Big);
                _a[k] = anew;
                return anew - ak;
            }

            /// <summary>
            /// A full sweep over all predictors. Returns true when converged.
            /// </summary>
            private bool FullSweep(double ab, double dem)
            {
                Passes++;
                _dlx = 0d;
                for (int k = 0; k < _p; k++)
                {
                    if (!_ju[k]) continue;
                    double ak = _a[k];
                    double del = Update(k, ab, dem);
                    if (_a[k] == ak) continue;
                    if (_mm[k] < 0) Activate(k);
                    _dlx = Math.Max(_dlx, _xv[k] * del * del);
                    _rsq += del * (2d * _g[k] - del * _xv[k]);
                    // Excluded predictors have zero inner products, so the update can run over all of them.
                    CrossProductKernel.Axpy(-del, _c[_mm[k]], 0, _g, 0, _p);
                }
                return _dlx < _model.Tolerance;
            }

            /// <summary>
            /// Sweeps over the active set until converged, then brings the inactive gradients up to date.
            /// </summary>
            private void PartialFit(double ab, double dem)
            {
                _warm = true;
                for (int l = 0; l < _nin; l++) _da[l] = _a[_ia[l]];
                while (true)
                {
                    Passes++;
                    _dlx = 0d;
                    for (int l = 0; l < _nin; l++)
                    {
                        int k = _ia[l];
                        double ak = _a[k];
                        double del = Update(k, ab, dem);
                        if (_a[k] == ak) continue;
                        _dlx = Math.Max(_dlx, _xv[k] * del * del);
                        _rsq += del * (2d * _g[k] - del * _xv[k]);
                        var c = _c[l];
                        for (int i = 0; i < _nin; i++) _g[_ia[i]] -= c[_ia[i]] * del;
                    }
                    if (_dlx < _model.Tolerance) break;
                    if (Passes > _model.MaxIterations) throw new MaxIterationsException();
                }
                for (int l = 0; l < _nin; l++) _da[l] -= _a[_ia[l]];
                Array.Clear(_work, 0, _p);
                for (int l = 0; l < _nin; l++)
                    if (_da[l] != 0d) CrossProductKernel.Axpy(_da[l], _c[l], 0, _work, 0, _p);
                for (int j = 0; j < _p; j++)
                    if (_mm[j] < 0) _g[j] += _work[j];
            }

            /// <summary>
            /// Solves at one penalty, warm-started from the previous solution.
            /// </summary>
            private void FitPoint(double ab, double dem)
            {
                if (_warm) PartialFit(ab, dem);
                while (true)
                {
                    if (FullSweep(ab, dem)) return;
                    if (Passes > _model.MaxIterations) throw new MaxIterationsException();
                    PartialFit(ab, dem);
                }
            }

            /// <summary>
            /// Runs the path on the default sequence, or on the given decreasing penalties (original scale).
            /// </summary>
            public void Run(double[]? lambdas)
            {
                var model = _model;
                int p = _p;
                double beta = model.Alpha, omb = 1d - beta;
                bool user = lambdas != null;
                int nlam = user ? lambdas!.Length : model.NumberOfLambdas;
                double alf = 1d;
                if (!user)
                {
                    double ratio = double.IsNaN(model.LambdaMinRatio) ? (_cache.N < p ? 1e-2 : 1e-4) : model.LambdaMinRatio;
                    alf = Math.Pow(Math.Max(MinimumRatio, ratio), 1d / (nlam - 1d));
                }
                Lambda = new double[nlam];
                Intercept = new double[nlam];
                Rsq = new double[nlam];
                Beta = new double[nlam * p];
                int mnl = Math.Min(MinimumPathLength, nlam);
                double alm = 0d;
                for (int m = 0; m < nlam; m++)
                {
                    double alm0 = alm;
                    if (user) alm = lambdas![m] / _ys;
                    else if (m > 1) alm *= alf;
                    else if (m == 0) alm = Big;
                    else
                    {
                        alm0 = 0d;
                        for (int j = 0; j < p; j++)
                            if (_ju[j] && _vp[j] > 0d) alm0 = Math.Max(alm0, Math.Abs(_g[j]) / _vp[j]);
                        alm0 /= Math.Max(beta, 1e-3);
                        alm = alm0 * alf;
                    }
                    double rsq0 = _rsq;
                    try
                    {
                        FitPoint(alm * beta, alm * omb);
                    }
                    catch (MaxIterationsException)
                    {
                        Converged = false;
                        break;
                    }

                    // Store the point on the original scale.
                    double b0 = _cache.Shift[p] + _my;
                    for (int l = 0; l < _nin; l++)
                    {
                        int k = _ia[l];
                        double b = _a[k] * _ys / _xs[k];
                        Beta[m * p + k] = b;
                        b0 -= b * (_cache.Shift[k] + _mx[k]);
                    }
                    Intercept[m] = model.HasIntercept ? b0 : 0d;
                    Lambda[m] = alm * _ys;
                    Rsq[m] = _rsq;
                    PathLength = m + 1;

                    if (PathLength < mnl || user) continue;
                    double change = _rsq == 0d ? double.PositiveInfinity : (_rsq - rsq0) / _rsq;
                    if (change < MinimumDevianceChange || _rsq > MaximumDevianceRatio) break;
                }
                if (!user && PathLength > 0 && PathLength < 3)
                {
                    // Report λ_max in place of the internal sentinel for very short paths.
                    double lmax = 0d;
                    for (int j = 0; j < p; j++)
                        if (_ju[j] && _vp[j] > 0d) lmax = Math.Max(lmax, Math.Abs(_g[j]) / _vp[j]);
                    Lambda[0] = PathLength > 1 ? Lambda[1] / alf : lmax / Math.Max(beta, 1e-3) * _ys;
                }
            }
        }

        /// <summary>
        /// Signals that coordinate descent reached the maximum number of sweeps.
        /// </summary>
        private sealed class MaxIterationsException : Exception
        {
        }

        #endregion
    }
}
