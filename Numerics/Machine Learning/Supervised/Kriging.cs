using System;

namespace Numerics.MachineLearning
{
    /// <summary>
    /// Simple, ordinary and universal kriging with a variogram model.
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
    /// Kriging is the best linear unbiased predictor of a spatial variable from observations at known locations.
    /// The variogram <c>γ(h) = c₀ + c (1 - ρ(h))</c> (see <see cref="VariogramModel"/>) defines the covariance
    /// <c>C(h) = c ρ(h)</c> for h &gt; 0 and <c>C(0) = c₀ + c</c>. The predictor is computed by
    /// <see cref="GaussianProcessRegression"/> with fixed hyperparameters: the partial sill is the signal variance and
    /// the nugget is the noise variance. Simple kriging uses a known mean; ordinary and universal kriging estimate a
    /// constant or a linear trend by generalized least squares, which gives the same predictor and kriging variance
    /// as solving the kriging system with Lagrange multipliers.
    /// </para>
    /// <para>
    /// As in classical geostatistics, the nugget is part of the covariance at zero distance, so the predictor
    /// honors the data at an observed location and the kriging variance there is zero. Elsewhere the kriging
    /// variance includes the nugget; it is the variance of the prediction error <c>Z(x₀) - Ẑ(x₀)</c>.
    /// </para>
    /// <para>
    /// <b> References: </b>
    /// </para>
    /// <list type="bullet">
    /// <item><description>
    /// Isaaks, E. H., and Srivastava, R. M. (1989). An Introduction to Applied Geostatistics. Oxford University Press.
    /// </description></item>
    /// <item><description>
    /// Cressie, N. (1993). Statistics for Spatial Data, revised edition. Wiley.
    /// </description></item>
    /// <item><description>
    /// Pebesma, E. J. (2004). Multivariable geostatistics in S: the gstat package. Computers and Geosciences, 30(7), 683-691.
    /// </description></item>
    /// </list>
    /// </remarks>
    public class Kriging
    {

        #region Construction

        /// <summary>
        /// Creates a kriging predictor using Euclidean distance between coordinates.
        /// </summary>
        /// <param name="coordinates">The observation locations, with one row per observation.</param>
        /// <param name="values">The observed values.</param>
        /// <param name="model">The variogram model.</param>
        /// <param name="nugget">The nugget c₀, non-negative.</param>
        /// <param name="partialSill">The partial sill c, positive.</param>
        /// <param name="range">The range parameter a, positive.</param>
        /// <param name="type">The kriging variant. Default = ordinary.</param>
        /// <param name="mean">The known mean used by simple kriging. Ignored otherwise. Default = 0.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the data are empty, mismatched, or not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a variogram parameter is outside its domain, or the model or kriging type is not a defined value.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the covariance matrix cannot be factored even with diagonal jitter, or the trend basis is rank deficient, for example universal kriging with fewer observations than trend coefficients.</exception>
        public Kriging(double[,] coordinates, double[] values, VariogramModel model, double nugget, double partialSill, double range,
            KrigingType type = KrigingType.Ordinary, double mean = 0d)
            : this(coordinates, values, model, nugget, partialSill, range, null, type, mean) { }

        /// <summary>
        /// Creates a kriging predictor using a custom distance between locations, such as a great-circle distance.
        /// </summary>
        /// <param name="coordinates">The observation locations, with one row per observation.</param>
        /// <param name="values">The observed values.</param>
        /// <param name="model">The variogram model.</param>
        /// <param name="nugget">The nugget c₀, non-negative.</param>
        /// <param name="partialSill">The partial sill c, positive.</param>
        /// <param name="range">The range parameter a, positive, in the units of the distance.</param>
        /// <param name="distance">
        /// Optional. The distance between two locations, given as rows of <paramref name="coordinates"/>. It must be a
        /// metric for which the chosen model is a valid covariance. When null, the Euclidean distance is used.
        /// </param>
        /// <param name="type">The kriging variant. Default = ordinary.</param>
        /// <param name="mean">The known mean used by simple kriging. Ignored otherwise. Default = 0.</param>
        /// <exception cref="ArgumentNullException">Thrown when an argument is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the data are empty, mismatched, or not finite.</exception>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when a variogram parameter is outside its domain, or the model or kriging type is not a defined value.</exception>
        /// <exception cref="InvalidOperationException">Thrown when the covariance matrix cannot be factored even with diagonal jitter, or the trend basis is rank deficient, for example universal kriging with fewer observations than trend coefficients.</exception>
        public Kriging(double[,] coordinates, double[] values, VariogramModel model, double nugget, double partialSill, double range,
            Func<double[], double[], double>? distance, KrigingType type = KrigingType.Ordinary, double mean = 0d)
        {
            if (coordinates == null) throw new ArgumentNullException(nameof(coordinates));
            if (values == null) throw new ArgumentNullException(nameof(values));
            if (!Enum.IsDefined(typeof(VariogramModel), model)) throw new ArgumentOutOfRangeException(nameof(model), "The variogram model is not a defined value.");
            if (!Enum.IsDefined(typeof(KrigingType), type)) throw new ArgumentOutOfRangeException(nameof(type), "The kriging type is not a defined value.");
            if (!Tools.IsFinite(nugget) || nugget < 0d) throw new ArgumentOutOfRangeException(nameof(nugget), "The nugget must be non-negative and finite.");
            if (!Tools.IsFinite(partialSill) || partialSill <= 0d) throw new ArgumentOutOfRangeException(nameof(partialSill), "The partial sill must be positive and finite.");
            if (!Tools.IsFinite(range) || range <= 0d) throw new ArgumentOutOfRangeException(nameof(range), "The range must be positive and finite.");
            if (type == KrigingType.Simple && !Tools.IsFinite(mean)) throw new ArgumentOutOfRangeException(nameof(mean), "The known mean must be finite.");
            if (type == KrigingType.Universal && distance != null)
                throw new ArgumentException("Universal kriging uses a trend in the coordinates and requires the Euclidean distance.", nameof(distance));

            int n = coordinates.GetLength(0), d = coordinates.GetLength(1);
            var x = new double[n * d];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < d; j++) x[i * d + j] = coordinates[i, j];
            var y = (double[])values.Clone();
            if (type == KrigingType.Simple)
                for (int i = 0; i < y.Length; i++) y[i] -= mean;

            Model = model;
            Type = type;
            Nugget = nugget;
            PartialSill = partialSill;
            Range = range;
            KnownMean = type == KrigingType.Simple ? mean : double.NaN;
            var form = model == VariogramModel.Spherical ? CovarianceForm.Spherical
                     : model == VariogramModel.Exponential ? CovarianceForm.Exponential
                     : CovarianceForm.SquaredExponential;
            var trend = type == KrigingType.Simple ? GaussianProcessMean.Zero
                      : type == KrigingType.Ordinary ? GaussianProcessMean.Constant
                      : GaussianProcessMean.Linear;
            // gstat's Gaussian model exp(-(h/a)²) is the squared exponential exp(-r²/2) with ℓ = a/√2.
            double lengthScale = model == VariogramModel.Gaussian ? range / Math.Sqrt(2d) : range;
            _process = new GaussianProcessRegression(x, n, d, y, form, trend, distance, true);
            _process.SetHyperparameters(partialSill, new[] { lengthScale }, nugget);
        }

        #endregion

        #region Members

        private readonly GaussianProcessRegression _process;

        /// <summary>
        /// Gets the variogram model.
        /// </summary>
        public VariogramModel Model { get; }

        /// <summary>
        /// Gets the kriging variant.
        /// </summary>
        public KrigingType Type { get; }

        /// <summary>
        /// Gets the nugget c₀.
        /// </summary>
        public double Nugget { get; }

        /// <summary>
        /// Gets the partial sill c.
        /// </summary>
        public double PartialSill { get; }

        /// <summary>
        /// Gets the total sill, c₀ + c.
        /// </summary>
        public double Sill => Nugget + PartialSill;

        /// <summary>
        /// Gets the range parameter a.
        /// </summary>
        public double Range { get; }

        /// <summary>
        /// Gets the known mean of simple kriging, or NaN for the other variants.
        /// </summary>
        public double KnownMean { get; }

        /// <summary>
        /// Gets the number of observations.
        /// </summary>
        public int SampleSize => _process.SampleSize;

        /// <summary>
        /// Gets the number of coordinate dimensions.
        /// </summary>
        public int NumberOfDimensions => _process.NumberOfFeatures;

        /// <summary>
        /// Gets the generalized least squares trend coefficients: empty for simple kriging, [β₀] for ordinary kriging,
        /// and [β₀, β₁, ..., β_d] for universal kriging.
        /// </summary>
        public double[] TrendCoefficients => _process.MeanCoefficients;

        #endregion

        #region Methods

        /// <summary>
        /// Returns the semivariance γ(h) of the variogram model.
        /// </summary>
        /// <param name="h">The separation distance, non-negative.</param>
        /// <returns>The semivariance, which is zero at h = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when h is negative or NaN.</exception>
        public double Semivariance(double h)
        {
            if (!(h >= 0d)) throw new ArgumentOutOfRangeException(nameof(h), "The distance must be non-negative.");
            return h == 0d ? 0d : Sill - Covariance(h);
        }

        /// <summary>
        /// Returns the covariance C(h) of the variogram model.
        /// </summary>
        /// <param name="h">The separation distance, non-negative.</param>
        /// <returns>The covariance, which is the total sill at h = 0.</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when h is negative or NaN.</exception>
        public double Covariance(double h)
        {
            if (!(h >= 0d)) throw new ArgumentOutOfRangeException(nameof(h), "The distance must be non-negative.");
            if (h == 0d) return Sill;
            double r = h / Range;
            switch (Model)
            {
                case VariogramModel.Spherical:
                    return r >= 1d ? 0d : PartialSill * (1d - r * (1.5 - 0.5 * r * r));
                case VariogramModel.Exponential:
                    return PartialSill * Math.Exp(-r);
                default:
                    return PartialSill * Math.Exp(-r * r);
            }
        }

        /// <summary>
        /// Returns the kriging predictions.
        /// </summary>
        /// <param name="x">The prediction locations, with one row per location.</param>
        /// <returns>The predicted value at each location.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when x has the wrong number of columns.</exception>
        public double[] Predict(double[,] x)
        {
            var mean = _process.Predict(x);
            if (Type == KrigingType.Simple)
                for (int i = 0; i < mean.Length; i++) mean[i] += KnownMean;
            return mean;
        }

        /// <summary>
        /// Returns the kriging predictions and kriging variances.
        /// </summary>
        /// <param name="x">The prediction locations, with one row per location.</param>
        /// <returns>The predicted value and the kriging (prediction error) variance at each location.</returns>
        /// <exception cref="ArgumentNullException">Thrown when x is null.</exception>
        /// <exception cref="ArgumentException">Thrown when x has the wrong number of columns.</exception>
        public (double[] Mean, double[] Variance) PredictWithVariance(double[,] x)
        {
            var (mean, variance) = _process.PredictWithVariance(x, true);
            if (Type == KrigingType.Simple)
                for (int i = 0; i < mean.Length; i++) mean[i] += KnownMean;
            return (mean, variance);
        }

        /// <summary>
        /// Returns the kriging prediction and kriging variance at one location.
        /// </summary>
        /// <param name="location">The coordinates of the location.</param>
        /// <returns>The predicted value and the kriging variance.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the location is null.</exception>
        /// <exception cref="ArgumentException">Thrown when the location has the wrong number of coordinates.</exception>
        public (double Mean, double Variance) PredictWithVariance(double[] location)
        {
            if (location == null) throw new ArgumentNullException(nameof(location));
            var x = new double[1, location.Length];
            for (int j = 0; j < location.Length; j++) x[0, j] = location[j];
            var (mean, variance) = PredictWithVariance(x);
            return (mean[0], variance[0]);
        }

        /// <summary>
        /// Returns the leave-one-out cross-validation predictions and kriging variances.
        /// </summary>
        /// <returns>
        /// For each observation, the prediction from all other observations and its kriging variance. The residuals are
        /// the observed values minus the predictions. Equivalent to <c>krige.cv</c> in gstat with all observations in
        /// the neighborhood.
        /// </returns>
        /// <remarks>Uses the closed form of Dubrule (1983), with the trend re-estimated for each left-out observation.</remarks>
        public (double[] Mean, double[] Variance) CrossValidate()
        {
            var (mean, variance) = _process.LeaveOneOut();
            if (Type == KrigingType.Simple)
                for (int i = 0; i < mean.Length; i++) mean[i] += KnownMean;
            return (mean, variance);
        }

        #endregion

    }
}
