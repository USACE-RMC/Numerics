using Numerics.Mathematics.LinearAlgebra;
using System;
using System.Threading.Tasks;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// The basis on which principal components are computed.
    /// </summary>
    public enum PrincipalComponentBasis
    {
        /// <summary>
        /// Components of the sample covariance matrix. The data are centered but not scaled.
        /// </summary>
        Covariance,

        /// <summary>
        /// Components of the sample correlation matrix. The data are centered and scaled to unit sample
        /// standard deviation.
        /// </summary>
        Correlation
    }

    /// <summary>
    /// Principal component analysis (PCA).
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
    /// PCA finds the orthogonal directions of maximum variance in a data set. With the n × p data matrix X
    /// centered (and, on the correlation basis, scaled to unit standard deviation) to Z, the sample
    /// covariance matrix S = ZᵀZ / (n − 1) has the eigen-decomposition S = V Λ Vᵀ. The columns of V are the
    /// principal axes, the eigenvalues λ₁ ≥ λ₂ ≥ ... ≥ λₚ ≥ 0 are the variances explained by each component,
    /// and the scores are T = Z V.
    /// </para>
    /// <para>
    /// S is formed directly from the row-major data, shifted by their first row, by a SIMD rank-8 update of the cross product,
    /// corrected for the column means as S = (Σ x̃x̃ᵀ − n d dᵀ)/(n − 1) with d the mean of the shifted data (the shift
    /// is a data value, so the correction loses little precision), and decomposed by
    /// Householder tridiagonalization and the implicit QL algorithm (Wilkinson and Reinsch, 1971). The cost is
    /// O(n p²) for the cross product plus O(p³) for the decomposition. Eigenvalues much smaller than ε·λ₁ lose
    /// relative accuracy because S squares the condition number of Z; round-off below zero is set to zero.
    /// </para>
    /// <para>
    /// <b> Sign convention: </b> each principal axis is oriented so that its element with the largest absolute
    /// value is positive (the first such element on ties). This is the convention of scikit-learn's
    /// <c>PCA</c>. R's <c>prcomp</c> does not fix the sign, so its axes and scores agree up to a sign per component.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <para>
    /// Jolliffe, I. T. (2002). Principal Component Analysis, 2nd ed. Springer, New York.
    /// </para>
    /// <para>
    /// Wilkinson, J. H., and Reinsch, C. (1971). Handbook for Automatic Computation, Volume II: Linear Algebra.
    /// Springer, Berlin.
    /// </para>
    /// </remarks>
    public class PrincipalComponentAnalysis
    {
        #region Construction

        /// <summary>
        /// Creates a new principal component analysis.
        /// </summary>
        /// <param name="X">The n × p data matrix, with observations in rows and variables in columns.</param>
        /// <param name="basis">The covariance or correlation basis. Default = covariance.</param>
        /// <param name="numberOfComponents">The number of components to retain. Zero or negative retains min(n, p).</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="X"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when X has fewer than two rows or no columns.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="numberOfComponents"/> exceeds min(n, p).</exception>
        public PrincipalComponentAnalysis(double[,] X, PrincipalComponentBasis basis = PrincipalComponentBasis.Covariance, int numberOfComponents = 0)
        {
            if (X == null) throw new ArgumentNullException(nameof(X));
            NumberOfSamples = X.GetLength(0);
            NumberOfFeatures = X.GetLength(1);
            if (NumberOfSamples < 2) throw new ArgumentException("There must be at least two rows of data.", nameof(X));
            if (NumberOfFeatures < 1) throw new ArgumentException("There must be at least one column of data.", nameof(X));
            int maxComponents = Math.Min(NumberOfSamples, NumberOfFeatures);
            if (numberOfComponents > maxComponents)
                throw new ArgumentOutOfRangeException(nameof(numberOfComponents), "The number of components cannot exceed min(n, p).");
            NumberOfComponents = numberOfComponents <= 0 ? maxComponents : numberOfComponents;
            Basis = basis;
            _data = X;
        }

        /// <summary>
        /// Creates a new principal component analysis.
        /// </summary>
        /// <param name="X">The n × p data matrix, with observations in rows and variables in columns.</param>
        /// <param name="basis">The covariance or correlation basis. Default = covariance.</param>
        /// <param name="numberOfComponents">The number of components to retain. Zero or negative retains min(n, p).</param>
        /// <exception cref="ArgumentNullException">Thrown when <paramref name="X"/> is null.</exception>
        /// <exception cref="ArgumentException">Thrown when X has fewer than two rows or no columns.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when <paramref name="numberOfComponents"/> exceeds min(n, p).</exception>
        public PrincipalComponentAnalysis(Matrix X, PrincipalComponentBasis basis = PrincipalComponentBasis.Covariance, int numberOfComponents = 0)
            : this(X?.Array!, basis, numberOfComponents)
        {
        }

        #endregion

        #region Members

        private readonly double[,] _data;
        private double[,]? _scores;

        /// <summary>
        /// The number of observations (rows), n.
        /// </summary>
        public int NumberOfSamples { get; private set; }

        /// <summary>
        /// The number of variables (columns), p.
        /// </summary>
        public int NumberOfFeatures { get; private set; }

        /// <summary>
        /// The number of retained components, k.
        /// </summary>
        public int NumberOfComponents { get; private set; }

        /// <summary>
        /// The covariance or correlation basis.
        /// </summary>
        public PrincipalComponentBasis Basis { get; private set; }

        /// <summary>
        /// Determines whether the analysis has been trained.
        /// </summary>
        public bool IsTrained { get; private set; }

        /// <summary>
        /// The column means of the training data, of length p.
        /// </summary>
        public double[] Means { get; private set; } = null!;

        /// <summary>
        /// The column scales applied after centering, of length p: the sample standard deviations (n − 1 divisor)
        /// on the correlation basis, and ones on the covariance basis.
        /// </summary>
        public double[] Scales { get; private set; } = null!;

        /// <summary>
        /// The k × p matrix of principal axes. Row i is the unit-length direction of component i
        /// (scikit-learn <c>components_</c>; the transpose of R <c>prcomp</c> <c>rotation</c>).
        /// </summary>
        public double[,] Components { get; private set; } = null!;

        /// <summary>
        /// The p × k matrix of loadings, L[j, i] = V[j, i]·√λᵢ. On the correlation basis these are the
        /// correlations between the variables and the component scores.
        /// </summary>
        public double[,] Loadings { get; private set; } = null!;

        /// <summary>
        /// The variance explained by each retained component, λ₁ ≥ ... ≥ λₖ (n − 1 divisor; scikit-learn
        /// <c>explained_variance_</c>, the square of R <c>prcomp</c> <c>sdev</c>).
        /// </summary>
        public double[] ExplainedVariance { get; private set; } = null!;

        /// <summary>
        /// The fraction of the total variance explained by each retained component.
        /// </summary>
        public double[] ExplainedVarianceRatio { get; private set; } = null!;

        /// <summary>
        /// The cumulative fraction of the total variance explained by the first i + 1 components.
        /// </summary>
        public double[] CumulativeExplainedVarianceRatio { get; private set; } = null!;

        /// <summary>
        /// The singular values of the centered (and scaled) data matrix, √((n − 1)·λᵢ).
        /// </summary>
        public double[] SingularValues { get; private set; } = null!;

        /// <summary>
        /// The total variance, the trace of the covariance or correlation matrix.
        /// </summary>
        public double TotalVariance { get; private set; }

        /// <summary>
        /// The n × k matrix of component scores of the training data. It is computed on first access.
        /// </summary>
        /// <exception cref="InvalidOperationException">Thrown when the analysis has not been trained.</exception>
        public double[,] Scores
        {
            get
            {
                if (!IsTrained) throw new InvalidOperationException("The analysis must be trained first.");
                return _scores ??= Transform(_data);
            }
        }

        #endregion

        #region Methods

        /// <summary>
        /// Computes the principal components of the data.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown when the data contain a non-finite value, or when a column is
        /// constant on the correlation basis.</exception>
        public void Train()
        {
            int n = NumberOfSamples, p = NumberOfFeatures, k = NumberOfComponents;
            var means = new double[p];
            var scales = new double[p];
            bool scale = Basis == PrincipalComponentBasis.Correlation;
            var X = _data;

            // Cross product of the data shifted by their first row, with the column sums of the shifted data.
            // The shift keeps the centering correction below well conditioned.
            var shifts = new double[p];
            for (int c = 0; c < p; c++) shifts[c] = X[0, c];
            var sums = new double[p];
            var S = new double[p * p];
            CrossProductKernel.RowMajorGram(X, shifts, sums, S);
            for (int c = 0; c < p; c++)
                if (!Tools.IsFinite(sums[c])) throw new ArgumentException("The data must be finite.");

            // Covariance corrected for the column means: S_ab = (Σ x̃_a x̃_b − n·d_a·d_b)/(n − 1), with d the mean of
            // the shifted column.
            var d = new double[p];
            for (int c = 0; c < p; c++)
            {
                d[c] = sums[c] / n;
                means[c] = shifts[c] + d[c];
            }
            // Each entry is formed once, from the lower triangle, and mirrored, so S is exactly symmetric.
            for (int a = 0; a < p; a++)
            {
                for (int b = 0; b <= a; b++)
                {
                    double v = (S[a * p + b] - n * d[a] * d[b]) / (n - 1);
                    S[a * p + b] = v;
                    S[b * p + a] = v;
                }
            }
            for (int c = 0; c < p; c++) scales[c] = 1d;
            if (scale)
            {
                for (int c = 0; c < p; c++)
                {
                    double v = S[c * p + c];
                    if (!(v > 0d)) throw new ArgumentException("A constant column cannot be scaled to unit variance on the correlation basis.");
                    scales[c] = Math.Sqrt(v);
                }
                for (int a = 0; a < p; a++)
                    for (int b = 0; b < p; b++)
                        S[a * p + b] = a == b ? 1d : S[a * p + b] / (scales[a] * scales[b]);
            }
            double total = 0d;
            for (int a = 0; a < p; a++) total += S[a * p + a];

            EigenValueDecomposition.DecomposeTridiagonal(S, p, out var eigenvalues, out var axes);

            Components = new double[k, p];
            Loadings = new double[p, k];
            ExplainedVariance = new double[k];
            ExplainedVarianceRatio = new double[k];
            CumulativeExplainedVarianceRatio = new double[k];
            SingularValues = new double[k];
            double cumulative = 0d;
            for (int i = 0; i < k; i++)
            {
                int o = i * p, imax = 0;
                double amax = -1d;
                for (int j = 0; j < p; j++)
                {
                    double a = Math.Abs(axes[o + j]);
                    if (a > amax) { amax = a; imax = j; }
                }
                double sign = axes[o + imax] < 0d ? -1d : 1d;
                // Round-off can leave eigenvalues of a singular matrix slightly negative.
                double lambda = Math.Max(0d, eigenvalues[i]);
                double root = Math.Sqrt(lambda);
                for (int j = 0; j < p; j++)
                {
                    double v = sign * axes[o + j];
                    Components[i, j] = v;
                    Loadings[j, i] = v * root;
                }
                ExplainedVariance[i] = lambda;
                ExplainedVarianceRatio[i] = total > 0d ? lambda / total : 0d;
                cumulative += ExplainedVarianceRatio[i];
                CumulativeExplainedVarianceRatio[i] = cumulative;
                SingularValues[i] = Math.Sqrt(lambda * (n - 1));
            }
            Means = means;
            Scales = scales;
            TotalVariance = total;
            _scores = null;
            IsTrained = true;
        }

        /// <summary>
        /// Projects observations onto the retained principal components.
        /// </summary>
        /// <param name="X">An m × p matrix of observations.</param>
        /// <returns>The m × k matrix of component scores.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the analysis has not been trained.</exception>
        /// <exception cref="ArgumentException">Thrown when the number of columns is not p.</exception>
        public double[,] Transform(double[,] X)
        {
            if (!IsTrained) throw new InvalidOperationException("The analysis must be trained first.");
            if (X == null) throw new ArgumentNullException(nameof(X));
            int m = X.GetLength(0), p = NumberOfFeatures, k = NumberOfComponents;
            if (X.GetLength(1) != p) throw new ArgumentException("The number of columns must equal the number of features.", nameof(X));
            // Axes stored column-contiguous per component for the inner product.
            var axes = new double[k * p];
            for (int i = 0; i < k; i++)
                for (int j = 0; j < p; j++)
                    axes[i * p + j] = Components[i, j] / Scales[j];
            var result = new double[m, k];
            Parallel.For(0, (m + 255) / 256, b =>
            {
                var row = new double[p];
                int r1 = Math.Min(m, (b + 1) * 256);
                for (int r = b * 256; r < r1; r++)
                {
                    for (int j = 0; j < p; j++) row[j] = X[r, j] - Means[j];
                    for (int i = 0; i < k; i++)
                        result[r, i] = CrossProductKernel.Dot(row, 0, axes, i * p, p);
                }
            });
            return result;
        }

        /// <summary>
        /// Projects one observation onto the retained principal components.
        /// </summary>
        /// <param name="x">An observation of length p.</param>
        /// <returns>The k component scores.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the analysis has not been trained.</exception>
        /// <exception cref="ArgumentException">Thrown when the length is not p.</exception>
        public double[] Transform(double[] x)
        {
            if (x == null) throw new ArgumentNullException(nameof(x));
            var X = new double[1, x.Length];
            for (int j = 0; j < x.Length; j++) X[0, j] = x[j];
            var t = Transform(X);
            var result = new double[NumberOfComponents];
            for (int i = 0; i < result.Length; i++) result[i] = t[0, i];
            return result;
        }

        /// <summary>
        /// Maps component scores back to the original variable space, X̂ = T Vₖᵀ·diag(scales) + means.
        /// </summary>
        /// <param name="scores">An m × k matrix of component scores.</param>
        /// <returns>The m × p reconstruction. It is exact when k = p.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the analysis has not been trained.</exception>
        /// <exception cref="ArgumentException">Thrown when the number of columns is not k.</exception>
        public double[,] InverseTransform(double[,] scores)
        {
            if (!IsTrained) throw new InvalidOperationException("The analysis must be trained first.");
            if (scores == null) throw new ArgumentNullException(nameof(scores));
            int m = scores.GetLength(0), p = NumberOfFeatures, k = NumberOfComponents;
            if (scores.GetLength(1) != k) throw new ArgumentException("The number of columns must equal the number of components.", nameof(scores));
            var result = new double[m, p];
            for (int r = 0; r < m; r++)
            {
                for (int j = 0; j < p; j++)
                {
                    double s = 0d;
                    for (int i = 0; i < k; i++) s += scores[r, i] * Components[i, j];
                    result[r, j] = s * Scales[j] + Means[j];
                }
            }
            return result;
        }

        /// <summary>
        /// Maps one vector of component scores back to the original variable space.
        /// </summary>
        /// <param name="scores">The k component scores.</param>
        /// <returns>The reconstructed observation of length p.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the analysis has not been trained.</exception>
        /// <exception cref="ArgumentException">Thrown when the length is not k.</exception>
        public double[] InverseTransform(double[] scores)
        {
            if (scores == null) throw new ArgumentNullException(nameof(scores));
            var T = new double[1, scores.Length];
            for (int i = 0; i < scores.Length; i++) T[0, i] = scores[i];
            var x = InverseTransform(T);
            var result = new double[NumberOfFeatures];
            for (int j = 0; j < result.Length; j++) result[j] = x[0, j];
            return result;
        }

        #endregion
    }
}
